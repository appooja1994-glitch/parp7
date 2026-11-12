#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tag anticancer-annotated compounds, then cut the deck to a target size.

    python scripts/03_select.py --in data/merged.sdf --out data/deck_35k.sdf \\
        --coconut-csv data/raw/coconut_csv_lite-11-2026/coconut_complete.csv \\
        --anticancer-keys data/anticancer_inchikeys.txt \\
        --target 35000

Anticancer tagging comes from two places, either of which may be omitted:

  --coconut-csv      the COCONUT CSV carries a dataSource / collection column.
                     Rows tagged NPACT, NPCARE, AfroCancer or InPACdb are the
                     anticancer-curated collections inside COCONUT.
  --anticancer-keys  a plain list of InChIKeys, one per line, that you derived
                     yourself (NPASS cancer cell-line actives, CancerHSP, etc.).
                     This is the defensible route - you set the activity cutoff.

Selection keeps EVERY anticancer-tagged compound, then fills the remaining
slots with a diversity pick over the untagged ones, so the deck is both
enriched for prior evidence and spread across chemical space.

  --method leader    sphere exclusion (default). Fast; scales to ~10^6.
  --method maxmin    MaxMin picker. Better spread, much slower above ~10^5.
  --method random    seeded random. For a null-model control deck.
"""
from __future__ import print_function
__BUILD__ = "2026-11-12d"
import argparse, sys, os, csv, random

try:
    from rdkit import Chem, RDLogger, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.SimDivFilters import rdSimDivPickers
except ImportError:
    sys.exit("ERROR: RDKit not importable by %s" % sys.executable)
RDLogger.DisableLog("rdApp.*")

ANTICANCER_COLLECTIONS = ("npact", "npcare", "afrocancer", "inpacdb")


def load_coconut_anticancer(path, key_col=None, coll_col=None):
    """InChIKeys whose COCONUT record names an anticancer-curated collection."""
    keys = set()
    if not path:
        return keys
    if not os.path.exists(path):
        print("  ! --coconut-csv not found, skipping: %s" % path)
        return keys
    csv.field_size_limit(10 ** 7)
    with open(path, "r") as fh:
        rdr = csv.DictReader(fh)
        cols = rdr.fieldnames or []
        low = dict((c.lower().strip(), c) for c in cols)
        def pick(names, override):
            if override:
                if override in cols: return override
                return low.get(override.lower())
            for n in names:
                if n in low: return low[n]
            return None
        kcol = pick(("standard_inchi_key", "inchikey", "inchi_key",
                     "standard_inchikey"), key_col)
        scol = pick(("collections", "collection", "datasource", "data_source",
                     "source", "sources", "tags"), coll_col)
        if not kcol or not scol:
            missing = []
            if not kcol: missing.append("InChIKey")
            if not scol: missing.append("collection/dataSource")
            print("  ! CSV has no %s column, so no compound can be tagged."
                  % " and no ".join(missing))
            print("    Columns seen: %s" % ", ".join(cols))
            print("    The COCONUT *lite* CSV does not carry collection tags.")
            print("    Use scripts/05_make_anticancer_keys.py to build a key")
            print("    list from a source that does, then pass it with")
            print("    --anticancer-keys. Skipping this file.")
            return keys
        for row in rdr:
            s = (row.get(scol) or "").lower()
            for name in ANTICANCER_COLLECTIONS:
                if name in s:
                    k = (row.get(kcol) or "").strip()
                    if k:
                        keys.add(k)
                    break
    print("  anticancer keys from COCONUT collections : %d" % len(keys))
    return keys


def load_key_list(path):
    keys = set()
    if not path:
        return keys
    if not os.path.exists(path):
        print("  ! --anticancer-keys not found, skipping: %s" % path)
        return keys
    with open(path) as fh:
        for line in fh:
            k = line.strip().split()[0] if line.strip() else ""
            if k and not k.startswith("#"):
                keys.add(k)
    print("  anticancer keys from user list           : %d" % len(keys))
    return keys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--coconut-csv", default=None)
    ap.add_argument("--csv-key-col", default=None,
                    help="override the InChIKey column name in --coconut-csv")
    ap.add_argument("--csv-collection-col", default=None,
                    help="override the collection column name in --coconut-csv")
    ap.add_argument("--require-annotations", action="store_true",
                    help="abort rather than build an untagged library if no "
                         "anticancer annotation could be loaded")
    ap.add_argument("--anticancer-keys", default=None)
    ap.add_argument("--target", type=int, default=35000)
    ap.add_argument("--method", choices=("leader", "maxmin", "random"),
                    default="leader")
    ap.add_argument("--radius", type=float, default=0.65,
                    help="leader: Tanimoto DISTANCE threshold. Lower = more "
                         "compounds kept. Tune to hit --target.")
    ap.add_argument("--seed", type=int, default=0xC0FFEE)
    ap.add_argument("--manifest", default=None, help="TSV of the final deck")
    a = ap.parse_args()

    print("Loading anticancer annotations")
    ac_keys = (load_coconut_anticancer(a.coconut_csv, a.csv_key_col,
                                       a.csv_collection_col)
               | load_key_list(a.anticancer_keys))
    print("  total unique anticancer InChIKeys        : %d" % len(ac_keys))
    if not ac_keys:
        msg = ("No anticancer annotation was loaded. The resulting library "
               "would be a pure diversity sample with no enrichment for prior "
               "evidence, and the annotated/diversity split could not be "
               "reported.")
        if a.require_annotations:
            sys.exit("\nABORTING: " + msg +
                     "\nBuild a key list with scripts/05_make_anticancer_keys.py, "
                     "or rerun without --require-annotations to accept this.")
        print("\n  ! WARNING: " + msg)
        print("  ! Continuing anyway. Pass --require-annotations to make this "
              "an error.\n")

    print("\nReading %s" % a.inp)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    mols, fps, keys, tagged = [], [], [], []
    supp = Chem.ForwardSDMolSupplier(open(a.inp, "rb"), sanitize=True, removeHs=False)
    for mol in supp:
        if mol is None:
            continue
        k = mol.GetProp("inchikey") if mol.HasProp("inchikey") else ""
        if not k:
            try:
                k = Chem.MolToInchiKey(mol)
            except Exception:
                k = ""
        is_ac = 1 if (k and k in ac_keys) else 0
        mol.SetProp("anticancer_annotated", str(is_ac))
        mols.append(mol); keys.append(k); tagged.append(is_ac)
        fps.append(gen.GetFingerprint(mol))
    n = len(mols)
    n_ac = sum(tagged)
    print("  compounds read                           : %d" % n)
    print("  anticancer-annotated in deck             : %d" % n_ac)

    core = [i for i in range(n) if tagged[i]]
    rest = [i for i in range(n) if not tagged[i]]
    need = a.target - len(core)

    if need <= 0:
        print("\n! Anticancer core alone (%d) already meets the target (%d)."
              % (len(core), a.target))
        chosen = core[:a.target]
    elif need >= len(rest):
        print("\n! Deck is smaller than the target; keeping everything.")
        chosen = core + rest
    else:
        print("\nSelecting %d of %d untagged compounds by %s"
              % (need, len(rest), a.method))
        if a.method == "random":
            rnd = random.Random(a.seed)
            chosen = core + rnd.sample(rest, need)
        elif a.method == "maxmin":
            sub = [fps[i] for i in rest]
            picker = rdSimDivPickers.MaxMinPicker()
            idx = picker.LazyBitVectorPick(sub, len(sub), need, seed=a.seed)
            chosen = core + [rest[i] for i in idx]
        else:
            sub = [fps[i] for i in rest]
            picker = rdSimDivPickers.LeaderPicker()
            idx = picker.LazyBitVectorPick(sub, len(sub), a.radius)
            print("  leader picker at distance %.2f returned %d centroids"
                  % (a.radius, len(idx)))
            idx = list(idx)
            if len(idx) > need:
                idx = idx[:need]
            elif len(idx) < need:
                have = set(idx)
                extra = [i for i in range(len(sub)) if i not in have]
                rnd = random.Random(a.seed); rnd.shuffle(extra)
                idx += extra[:need - len(idx)]
                frac = 100.0 * (need - len(have)) / float(need)
                print("  ! topped up with %d additional compounds chosen at "
                      "RANDOM (%.0f%% of the selection)" % (need - len(have), frac))
                if frac > 5:
                    print("  ! This is not a diversity selection. Lower --radius "
                          "(try %.2f) so the picker returns more than %d "
                          "centroids, and rerun." % (max(0.30, a.radius - 0.10), need))
            chosen = core + [rest[i] for i in idx]

    chosen = sorted(set(chosen))
    w = Chem.SDWriter(a.out)
    for i in chosen:
        w.write(mols[i])
    w.close()

    if a.manifest:
        with open(a.manifest, "w") as fh:
            tw = csv.writer(fh, delimiter="\t")
            tw.writerow(["deck_id", "inchikey", "sources", "anticancer_annotated",
                         "pains_flag", "sugars_removed"])
            for i in chosen:
                m = mols[i]
                g = lambda p: m.GetProp(p) if m.HasProp(p) else ""
                tw.writerow([g("deck_id"), keys[i], g("sources"),
                             g("anticancer_annotated"), g("pains_flag"),
                             g("sugars_removed")])

    print("\nFinal deck : %d compounds (%d anticancer-annotated, %d diversity)"
          % (len(chosen), sum(tagged[i] for i in chosen),
             len(chosen) - sum(tagged[i] for i in chosen)))
    print("Written    : %s" % a.out)
    if a.manifest:
        print("Manifest   : %s" % a.manifest)


if __name__ == "__main__":
    main()
