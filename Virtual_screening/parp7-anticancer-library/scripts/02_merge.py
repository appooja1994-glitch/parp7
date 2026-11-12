#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge filtered SDFs into one deck, deduplicating across sources on InChIKey.

    python scripts/02_merge.py --out data/merged.sdf \\
        --add coconut:data/coconut_filtered.sdf \\
        --add imppat:data/imppat_filtered.sdf

Sources are merged in the order given; the first occurrence of an InChIKey wins
and later ones only contribute their source name to the `sources` property.
That keeps COCONUT's standardised coordinates as the canonical geometry while
still recording that IMPPAT also contains the compound.
"""
from __future__ import print_function
__BUILD__ = "2026-11-12d"
import argparse, sys, os

try:
    from rdkit import Chem, RDLogger
except ImportError:
    sys.exit("ERROR: RDKit not importable by %s" % sys.executable)
RDLogger.DisableLog("rdApp.*")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--add", action="append", required=True, metavar="NAME:PATH",
                    help="repeat per source, e.g. --add coconut:data/a.sdf")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    pairs = []
    for spec in a.add:
        if ":" not in spec:
            sys.exit("ERROR: --add needs NAME:PATH, got %r" % spec)
        name, path = spec.split(":", 1)
        if not os.path.exists(path):
            sys.exit("ERROR: no such file: %s" % path)
        pairs.append((name, path))

    index = {}        # inchikey -> (mol, [sources])
    order = []
    per_source = {}

    for name, path in pairs:
        n_new = n_dup = 0
        supp = Chem.ForwardSDMolSupplier(open(path, "rb"), sanitize=True, removeHs=False)
        for mol in supp:
            if mol is None:
                continue
            key = mol.GetProp("inchikey") if mol.HasProp("inchikey") else ""
            if not key:
                try:
                    key = Chem.MolToInchiKey(mol)
                except Exception:
                    continue
            if key in index:
                srcs = index[key][1]
                if name not in srcs:
                    srcs.append(name)
                n_dup += 1
            else:
                index[key] = (mol, [name])
                order.append(key)
                n_new += 1
        per_source[name] = (n_new, n_dup)
        print("%-12s %7d new, %7d already present" % (name, n_new, n_dup))

    w = Chem.SDWriter(a.out)
    for key in order:
        mol, srcs = index[key]
        mol.SetProp("sources", ",".join(srcs))
        w.write(mol)
    w.close()
    print("\nMerged deck: %d unique compounds -> %s" % (len(order), a.out))


if __name__ == "__main__":
    main()
