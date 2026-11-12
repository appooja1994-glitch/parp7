#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Remove revoked COCONUT records from an already-filtered SDF.

The COCONUT 3D SDF distribution contains the whole release, revoked entries
included. Filtering does not remove them, so they must be removed explicitly
against a whitelist of active identifiers. This runs on the OUTPUT of
01_filter.py, so the long filtering pass does not have to be repeated.

Check the whitelist size first:

    python scripts/06_drop_revoked.py --count data/raw/<coconut>.csv

Then subset:

    python scripts/06_drop_revoked.py \
        --in  data/coconut_filtered.sdf \
        --out data/coconut_active.sdf \
        --whitelist data/raw/<coconut>.csv \
        --dropped logs/revoked_dropped.txt

A whitelist may be a CSV/TSV with an identifier column, or a plain text file
with one identifier per line. Matching is on the COCONUT identifier carried
through filtering as deck_id, not on InChIKey, because deglycosylation changes
structures and therefore changes their keys.
"""
from __future__ import print_function
__BUILD__ = "2026-11-12d"
import argparse, csv, gzip, os, re, sys

try:
    from rdkit import Chem, RDLogger
except ImportError:
    sys.exit("ERROR: RDKit not importable by %s" % sys.executable)
RDLogger.DisableLog("rdApp.*")

ID_NAMES = ("identifier", "coconut_id", "id", "cnp_id", "name")
CNP = re.compile(r"\bCNP\d{6,}\b")


def openf(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")


def load_whitelist(path, col=None):
    """Identifiers from a delimited file or a plain list."""
    ids = set()
    with openf(path) as fh:
        head = fh.read(16384)
        fh.seek(0)
        looks_delimited = ("," in head or "\t" in head) and "\n" in head
        if looks_delimited:
            delim = "\t" if head.count("\t") > head.count(",") else ","
            csv.field_size_limit(10 ** 7)
            rdr = csv.DictReader(fh, delimiter=delim)
            cols = rdr.fieldnames or []
            low = dict((c.lower().strip(), c) for c in cols)
            icol = None
            if col:
                icol = col if col in cols else low.get(col.lower())
            else:
                for n in ID_NAMES:
                    if n in low:
                        icol = low[n]
                        break
            if not icol:
                sys.exit("ERROR: no identifier column found in %s\n"
                         "       columns: %s\n"
                         "       pass --id-col to name it explicitly"
                         % (path, ", ".join(cols)))
            print("  identifier column: %s" % icol)
            for row in rdr:
                v = (row.get(icol) or "").strip()
                if v:
                    ids.add(v)
        else:
            for line in fh:
                v = line.strip().split()[0] if line.strip() else ""
                if v and not v.startswith("#"):
                    ids.add(v)
    return ids


def main():
    ap = argparse.ArgumentParser(description=__doc__,
          formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--count", metavar="FILE",
                    help="report how many identifiers a whitelist holds, then exit")
    ap.add_argument("--in", dest="inp")
    ap.add_argument("--out", dest="out")
    ap.add_argument("--whitelist")
    ap.add_argument("--id-col", default=None)
    ap.add_argument("--dropped", default=None,
                    help="write the identifiers that were removed")
    ap.add_argument("--report-every", type=int, default=100000)
    a = ap.parse_args()

    if a.count:
        print("\n%s" % a.count)
        ids = load_whitelist(a.count, a.id_col)
        print("  identifiers: %d" % len(ids))
        print("""
  COCONUT published 738,795 molecules, 38,869 of them revoked, leaving
  699,926 active. A whitelist near 699,926 is the active set and can be
  used. One near 738,795 includes revoked entries and is not a whitelist,
  so another source of active identifiers is needed.""")
        return

    for need, name in ((a.inp, "--in"), (a.out, "--out"), (a.whitelist, "--whitelist")):
        if not need:
            sys.exit("ERROR: %s is required (or use --count)" % name)

    print("Loading whitelist")
    keep = load_whitelist(a.whitelist, a.id_col)
    print("  active identifiers: %d" % len(keep))
    if not keep:
        sys.exit("ERROR: whitelist is empty")

    d = os.path.dirname(os.path.abspath(a.out))
    if d and not os.path.isdir(d):
        os.makedirs(d)

    w = Chem.SDWriter(a.out)
    dropped = []
    n = kept = nomatch = 0
    fh = open(a.inp, "rb")
    for mol in Chem.ForwardSDMolSupplier(fh, sanitize=False, removeHs=False):
        n += 1
        if a.report_every and n % a.report_every == 0:
            sys.stderr.write("  ...%d read, %d kept\n" % (n, kept))
        if mol is None:
            continue
        cid = ""
        for p in ("deck_id", "identifier", "coconut_id"):
            if mol.HasProp(p):
                cid = mol.GetProp(p).strip()
                break
        if not cid and mol.HasProp("_Name"):
            cid = mol.GetProp("_Name").strip()
        if not cid:
            nomatch += 1
            continue
        if cid in keep:
            w.write(mol)
            kept += 1
        else:
            dropped.append(cid)
    w.close()
    fh.close()

    if a.dropped:
        with open(a.dropped, "w") as f:
            for c in dropped:
                f.write(c + "\n")

    print("\nRead            : %d" % n)
    print("Kept (active)   : %d" % kept)
    print("Dropped         : %d" % len(dropped))
    if nomatch:
        print("No identifier   : %d" % nomatch)
    print("Written         : %s" % a.out)
    if a.dropped:
        print("Dropped ids     : %s" % a.dropped)

    if n and len(dropped) > 0.25 * n:
        print("\n! More than a quarter of the input was dropped. That is far more"
              "\n! than the 5 per cent of COCONUT that is revoked, so the"
              "\n! identifiers in the SDF and in the whitelist probably do not"
              "\n! use the same form. Check a few values before accepting this.")
        print("! sample dropped: %s" % ", ".join(dropped[:5]))
    elif len(dropped) == 0:
        print("\n! Nothing was dropped. Either the input was already active only,"
              "\n! or the whitelist contains every identifier including revoked"
              "\n! ones. Confirm with --count before concluding the first.")


if __name__ == "__main__":
    main()
