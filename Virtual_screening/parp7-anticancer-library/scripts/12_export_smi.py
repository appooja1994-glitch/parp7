#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Export an SDF to SMILES (.smi), the format for depositing a compound set in a
repository. One molecule per line, SMILES then a tab then the identifier.

    python scripts/12_export_smi.py --in data/anticancer_final.sdf \
        --out data/anticancer_final.smi

With --props a second tab separated file is written carrying the per compound
annotations, so that the chemistry and its provenance can both be deposited
without either file becoming unreadable.
"""
from __future__ import print_function
__BUILD__ = "2026-11-12d"
import argparse, os, sys

try:
    from rdkit import Chem, RDLogger
except ImportError:
    sys.exit("ERROR: RDKit not importable by %s" % sys.executable)
RDLogger.DisableLog("rdApp.*")

ANNOT = ("inchikey", "sources", "n_sources", "anticancer_annotated",
         "pains_flag", "sugars_removed", "best_activity_um", "cancer_cell_line")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
          formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--props", default=None,
                    help="also write a TSV of the per compound annotations")
    ap.add_argument("--canonical", action="store_true", default=True)
    ap.add_argument("--no-header", action="store_true",
                    help="omit the leading comment line from the .smi")
    a = ap.parse_args()

    if not os.path.exists(a.inp):
        sys.exit("ERROR: no such file: %s" % a.inp)
    d = os.path.dirname(os.path.abspath(a.out))
    if d and not os.path.isdir(d):
        os.makedirs(d)

    n = ok = bad = 0
    seen = {}
    fh = open(a.out, "w")
    if not a.no_header:
        fh.write("# SMILES\tID   generated from %s\n" % os.path.basename(a.inp))
    pf = open(a.props, "w") if a.props else None
    if pf:
        pf.write("id\tsmiles\t" + "\t".join(ANNOT) + "\n")

    for mol in Chem.ForwardSDMolSupplier(open(a.inp, "rb"), removeHs=True):
        n += 1
        if mol is None:
            bad += 1
            continue
        cid = ""
        for p in ("deck_id", "src_id", "_Name"):
            if mol.HasProp(p) and mol.GetProp(p).strip():
                cid = mol.GetProp(p).strip()
                break
        if not cid:
            cid = "mol%07d" % n
        if cid in seen:
            seen[cid] += 1
            cid = "%s_%d" % (cid, seen[cid])
        else:
            seen[cid] = 0
        # whitespace in an id would break the two-column .smi format
        cid = "_".join(cid.split())
        try:
            smi = Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            bad += 1
            continue
        if not smi:
            bad += 1
            continue
        fh.write("%s\t%s\n" % (smi, cid))
        ok += 1
        if pf:
            vals = [mol.GetProp(k) if mol.HasProp(k) else "" for k in ANNOT]
            pf.write("%s\t%s\t%s\n" % (cid, smi, "\t".join(vals)))
    fh.close()
    if pf:
        pf.close()

    print("Read    : %d" % n)
    print("Written : %d  ->  %s" % (ok, a.out))
    if bad:
        print("Skipped : %d (unparseable or no SMILES)" % bad)
    if a.props:
        print("Props   : %s" % a.props)
    return 0


if __name__ == "__main__":
    sys.exit(main())
