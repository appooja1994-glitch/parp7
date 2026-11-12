#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the anticancer InChIKey list that 03_select.py needs, from whatever
source you actually have. Schema-agnostic, because every natural-product
database names its columns differently.

Inspect a file before trusting it:

    python scripts/05_make_anticancer_keys.py --inspect data/raw/<file>.csv
    python scripts/05_make_anticancer_keys.py --inspect data/raw/<file>.sdf

Then build the list from one or more sources:

    # (a) COCONUT SDF or CSV that DOES carry a collection column
    python scripts/05_make_anticancer_keys.py \
        --source data/raw/coconut_3d.sdf --collection-col collections \
        --out data/anticancer_inchikeys.txt

    # (b) a whole file that is itself an anticancer dataset (NPCARE, NPACT)
    python scripts/05_make_anticancer_keys.py \
        --source data/raw/npcare.csv --all-rows \
        --out data/anticancer_inchikeys.txt

    # (c) several sources at once, merged and deduplicated
    python scripts/05_make_anticancer_keys.py \
        --source data/raw/npcare.csv --all-rows \
        --source data/raw/npass_cancer.tsv --all-rows \
        --out data/anticancer_inchikeys.txt

Each source may supply InChIKeys directly, or SMILES or InChI from which the
key is computed. Computing is slower but means a database that publishes only
structures is still usable.
"""
from __future__ import print_function
__BUILD__ = "2026-11-12d"
import argparse, csv, gzip, os, sys

try:
    from rdkit import Chem, RDLogger
except ImportError:
    sys.exit("ERROR: RDKit not importable by %s" % sys.executable)
RDLogger.DisableLog("rdApp.*")

ANTICANCER = ("npact", "npcare", "afrocancer", "inpacdb")

KEY_NAMES = ("standard_inchi_key", "inchikey", "inchi_key", "standard_inchikey",
             "inchi key", "key")
SMI_NAMES = ("canonical_smiles", "smiles", "smile", "isomeric_smiles",
             "canonical smiles", "structure")
INCHI_NAMES = ("standard_inchi", "inchi")
COLL_NAMES = ("collections", "collection", "datasource", "data_source",
              "source", "sources", "tags", "database", "db")


def openf(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")


def pick(cols, names, override=None):
    low = {c.lower().strip(): c for c in cols}
    if override:
        if override in cols:
            return override
        if override.lower() in low:
            return low[override.lower()]
        return None
    for n in names:
        if n in low:
            return low[n]
    return None


def key_from(keyval, smi, inchi):
    if keyval:
        k = keyval.strip()
        if len(k) in (14, 27) and "-" in k or len(k) == 14:
            return k
        if k:
            return k
    m = None
    if smi:
        m = Chem.MolFromSmiles(smi)
    if m is None and inchi:
        try:
            m = Chem.MolFromInchi(inchi)
        except Exception:
            m = None
    if m is None:
        return ""
    try:
        return Chem.MolToInchiKey(m)
    except Exception:
        return ""


def sdf_fields(path, limit=200):
    """Property names present in the first `limit` records of an SDF."""
    names, n = set(), 0
    fh = gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")
    for mol in Chem.ForwardSDMolSupplier(fh, sanitize=False):
        n += 1
        if mol is not None:
            names.update(mol.GetPropNames())
        if n >= limit:
            break
    fh.close()
    return sorted(names)


def inspect(path):
    print("\n%s" % path)
    if not os.path.exists(path):
        print("  ! file not found"); return
    if path.endswith((".sdf", ".sdf.gz")):
        f = sdf_fields(path)
        print("  type   : SDF")
        print("  fields : %s" % (", ".join(f) if f else "(none)"))
        print("  --> key column      : %s" % (pick(f, KEY_NAMES) or "NOT FOUND"))
        print("  --> smiles column   : %s" % (pick(f, SMI_NAMES) or "NOT FOUND"))
        print("  --> collection col  : %s" % (pick(f, COLL_NAMES) or "NOT FOUND"))
    else:
        csv.field_size_limit(10 ** 7)
        with openf(path) as fh:
            sample = fh.read(16384); fh.seek(0)
            delim = "\t" if sample.count("\t") > sample.count(",") else ","
            rdr = csv.DictReader(fh, delimiter=delim)
            cols = rdr.fieldnames or []
            row = next(rdr, None)
        print("  type   : delimited (%s)" % ("TAB" if delim == "\t" else "COMMA"))
        print("  columns: %s" % ", ".join(cols))
        print("  --> key column      : %s" % (pick(cols, KEY_NAMES) or "NOT FOUND"))
        print("  --> smiles column   : %s" % (pick(cols, SMI_NAMES) or "NOT FOUND"))
        print("  --> inchi column    : %s" % (pick(cols, INCHI_NAMES) or "NOT FOUND"))
        print("  --> collection col  : %s" % (pick(cols, COLL_NAMES) or "NOT FOUND"))
        if row:
            print("  first row sample:")
            for c in cols[:8]:
                v = (row.get(c) or "")[:60]
                print("      %-24s %s" % (c, v))


def harvest(path, all_rows, key_col, smi_col, coll_col, terms, report):
    keys, n, matched = set(), 0, 0
    is_sdf = path.endswith((".sdf", ".sdf.gz"))

    if is_sdf:
        fields = sdf_fields(path)
        kc = pick(fields, KEY_NAMES, key_col)
        sc = pick(fields, SMI_NAMES, smi_col)
        cc = pick(fields, COLL_NAMES, coll_col)
        if not all_rows and not cc:
            print("  ! no collection field in this SDF and --all-rows not set")
            print("    fields present: %s" % ", ".join(fields))
            return keys
        fh = gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")
        for mol in Chem.ForwardSDMolSupplier(fh, sanitize=False):
            n += 1
            if mol is None:
                continue
            if not all_rows:
                val = (mol.GetProp(cc) if mol.HasProp(cc) else "").lower()
                if not any(t in val for t in terms):
                    continue
            matched += 1
            kv = mol.GetProp(kc) if (kc and mol.HasProp(kc)) else ""
            sv = mol.GetProp(sc) if (sc and mol.HasProp(sc)) else ""
            if not kv and not sv:
                try:
                    Chem.SanitizeMol(mol); sv = Chem.MolToSmiles(mol)
                except Exception:
                    pass
            k = key_from(kv, sv, "")
            if k:
                keys.add(k)
            if report and n % report == 0:
                sys.stderr.write("    ...%d read, %d matched\n" % (n, matched))
        fh.close()
    else:
        csv.field_size_limit(10 ** 7)
        with openf(path) as fh:
            sample = fh.read(16384); fh.seek(0)
            delim = "\t" if sample.count("\t") > sample.count(",") else ","
            rdr = csv.DictReader(fh, delimiter=delim)
            cols = rdr.fieldnames or []
            kc = pick(cols, KEY_NAMES, key_col)
            sc = pick(cols, SMI_NAMES, smi_col)
            ic = pick(cols, INCHI_NAMES)
            cc = pick(cols, COLL_NAMES, coll_col)
            if not all_rows and not cc:
                print("  ! no collection column in this file and --all-rows not set")
                print("    columns present: %s" % ", ".join(cols))
                return keys
            if not (kc or sc or ic):
                print("  ! no InChIKey, SMILES or InChI column found")
                print("    columns present: %s" % ", ".join(cols))
                return keys
            for row in rdr:
                n += 1
                if not all_rows:
                    val = (row.get(cc) or "").lower()
                    if not any(t in val for t in terms):
                        continue
                matched += 1
                k = key_from(row.get(kc) if kc else "",
                             row.get(sc) if sc else "",
                             row.get(ic) if ic else "")
                if k:
                    keys.add(k)
                if report and n % report == 0:
                    sys.stderr.write("    ...%d read, %d matched\n" % (n, matched))

    print("  read %d records, %d matched, %d unique InChIKeys"
          % (n, matched, len(keys)))
    return keys


def main():
    ap = argparse.ArgumentParser(description=__doc__,
          formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inspect", nargs="+", metavar="FILE",
                    help="print the columns/fields of each file and exit")
    ap.add_argument("--source", action="append", default=[], metavar="FILE",
                    help="repeat per source file (.sdf, .csv, .tsv, optionally .gz)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--all-rows", action="store_true",
                    help="treat every record in the source as anticancer "
                         "(use for a file that IS an anticancer dataset)")
    ap.add_argument("--key-col", default=None, help="override the InChIKey column")
    ap.add_argument("--smiles-col", default=None, help="override the SMILES column")
    ap.add_argument("--collection-col", default=None,
                    help="override the collection/source column")
    ap.add_argument("--terms", default=",".join(ANTICANCER),
                    help="comma-separated collection names counted as anticancer "
                         "(default: %(default)s)")
    ap.add_argument("--report-every", type=int, default=100000)
    a = ap.parse_args()

    if a.inspect:
        for p in a.inspect:
            inspect(p)
        return
    if not a.source or not a.out:
        sys.exit("ERROR: need at least one --source and an --out "
                 "(or use --inspect to look first)")

    terms = [t.strip().lower() for t in a.terms.split(",") if t.strip()]
    allk = set()
    for p in a.source:
        print("\n%s" % p)
        if not os.path.exists(p):
            print("  ! not found, skipping"); continue
        allk |= harvest(p, a.all_rows, a.key_col, a.smiles_col,
                        a.collection_col, terms, a.report_every)

    d = os.path.dirname(os.path.abspath(a.out))
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with open(a.out, "w") as fh:
        for k in sorted(allk):
            fh.write(k + "\n")
    print("\nTotal unique anticancer InChIKeys: %d" % len(allk))
    print("Written: %s" % a.out)
    if not allk:
        print("\n! The list is empty. Pass --inspect on your source to see its "
              "columns, then set --collection-col or --all-rows accordingly.")


if __name__ == "__main__":
    main()
