#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate 3D coordinates and write one PDBQT per ligand for AutoDock Vina.

    python scripts/04_prepare_pdbqt.py \
        --in data/anticancer_final.sdf --out data/pdbqt --jobs 8

This is the only place geometry is built. Coordinates are embedded with
ETKDGv3 and optimised with MMFF94, then Meeko writes the explicit torsion
tree that AutoDock reads to decide which bonds rotate.

A 2D conformer has a conformer count of 1 but every z equal to zero, so
counting conformers is not enough to tell whether geometry already exists.
Both are checked, otherwise a flat molecule is docked as though it were real.
"""
from __future__ import print_function
__BUILD__ = "2026-11-12d"
import argparse, os, sys, multiprocessing as mp

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
except ImportError:
    sys.exit("ERROR: RDKit not importable by %s" % sys.executable)
RDLogger.DisableLog("rdApp.*")

try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    HAVE_MEEKO, MEEKO_ERR = True, ""
except Exception as e:
    HAVE_MEEKO, MEEKO_ERR = False, "%s: %s" % (type(e).__name__, e)

SAFE = "".join(chr(c) for c in range(48, 58)) + \
       "".join(chr(c) for c in range(65, 91)) + \
       "".join(chr(c) for c in range(97, 123)) + "-_."


def safe_name(s, fallback):
    out = "".join(c if c in SAFE else "_" for c in (s or "").strip())
    out = out.strip("._") or fallback
    return out[:120]


def needs_3d(mol):
    if mol.GetNumConformers() == 0:
        return True
    c = mol.GetConformer()
    if not c.Is3D():
        return True
    zs = set(round(c.GetAtomPosition(i).z, 3) for i in range(mol.GetNumAtoms()))
    return len(zs) <= 1


def one(args):
    molblock, name, seed = args
    try:
        mol = Chem.MolFromMolBlock(molblock, removeHs=False)
        if mol is None:
            return name, None, "unparseable"
        make3d = needs_3d(mol)
        mh = Chem.AddHs(mol, addCoords=not make3d)
        if make3d:
            mh.RemoveAllConformers()
            par = AllChem.ETKDGv3()
            par.randomSeed = seed
            if AllChem.EmbedMolecule(mh, par) != 0:
                par.useRandomCoords = True
                if AllChem.EmbedMolecule(mh, par) != 0:
                    return name, None, "embed_failed"
        try:
            AllChem.MMFFOptimizeMolecule(mh, maxIters=500)
        except Exception:
            pass
        prep = MoleculePreparation()
        setups = prep.prepare(mh)
        if not setups:
            return name, None, "prepare_failed"
        txt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            return name, None, "write_failed"
        return name, txt, ""
    except Exception as e:
        return name, None, type(e).__name__


def main():
    ap = argparse.ArgumentParser(description=__doc__,
          formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--jobs", type=int, default=0, help="0 means all cores")
    ap.add_argument("--seed", type=int, default=0xC0FFEE)
    ap.add_argument("--failures", default=None, help="TSV of ligands that failed")
    a = ap.parse_args()

    if not HAVE_MEEKO:
        print("=" * 64)
        print("Meeko could not be imported:")
        print("    " + MEEKO_ERR)
        print()
        print("Meeko 0.8.0 declares no dependencies, so pip installs none of")
        print("them. It needs all of these just to import:")
        print()
        print("    pip install meeko gemmi scipy pandas")
        print()
        print("Not falling back to Open Babel, because Open Babel writes no")
        print("torsion tree and the ligands would dock with the wrong degrees")
        print("of freedom without anything reporting a problem.")
        print("=" * 64)
        return 2

    if not os.path.isdir(a.out):
        os.makedirs(a.out)
    jobs = a.jobs if a.jobs > 0 else (mp.cpu_count() or 1)

    def feed():
        n = 0
        supp = Chem.ForwardSDMolSupplier(open(a.inp, "rb"), removeHs=False)
        seen = {}
        for mol in supp:
            n += 1
            if mol is None:
                continue
            nm = ""
            for p in ("deck_id", "src_id", "_Name"):
                if mol.HasProp(p):
                    nm = mol.GetProp(p)
                    if nm:
                        break
            nm = safe_name(nm, "mol%07d" % n)
            if nm in seen:
                seen[nm] += 1
                nm = "%s_%d" % (nm, seen[nm])
            else:
                seen[nm] = 0
            try:
                yield Chem.MolToMolBlock(mol), nm, a.seed
            except Exception:
                continue

    print("Preparing ligands with %d worker(s)" % jobs)
    nok = 0
    failures = []
    if jobs == 1:
        results = (one(x) for x in feed())
    else:
        pool = mp.Pool(jobs)
        results = pool.imap_unordered(one, feed(), chunksize=64)

    for i, (name, txt, err) in enumerate(results, 1):
        if txt is None:
            failures.append((name, err))
        else:
            with open(os.path.join(a.out, name + ".pdbqt"), "w") as fh:
                fh.write(txt)
            nok += 1
        if i % 2000 == 0:
            sys.stderr.write("  ...%d processed, %d written\n" % (i, nok))
    if jobs != 1:
        pool.close(); pool.join()

    total = nok + len(failures)
    if a.failures and failures:
        with open(a.failures, "w") as fh:
            fh.write("ligand\treason\n")
            for n_, e_ in failures:
                fh.write("%s\t%s\n" % (n_, e_))

    print("\nProcessed : %d" % total)
    print("Written   : %d  ->  %s" % (nok, a.out))
    print("Failed    : %d" % len(failures))
    if failures:
        reasons = {}
        for _, e_ in failures:
            reasons[e_] = reasons.get(e_, 0) + 1
        for r, c in sorted(reasons.items(), key=lambda x: -x[1]):
            print("    %-18s %d" % (r, c))
        if a.failures:
            print("  listed in %s" % a.failures)
        if total and len(failures) > 0.02 * total:
            print("\n! More than 2 per cent failed. That is higher than normal "
                  "and worth\n! looking into before docking.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
