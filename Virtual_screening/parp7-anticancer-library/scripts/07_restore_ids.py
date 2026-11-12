#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Repair placeholder REC######## identifiers in a filtered SDF.

01_filter.py assigns REC00000001, REC00000002 and so on when it cannot find an
identifier property in the source record. The number is the 1-based index of
the record in the source SDF, so the real identifier can be recovered by
reading the source again, without repeating the filtering pass.

Nothing is re-parsed. Both files are streamed as text, so coordinates,
stereochemistry and every other property are untouched.

    # 1. see what the source SDF actually calls its fields
    python scripts/07_restore_ids.py --list-fields data/raw/<coconut>.sdf

    # 2. rewrite the ids
    python scripts/07_restore_ids.py \
        --source   data/raw/<coconut>.sdf \
        --in       data/coconut_filtered.sdf \
        --out      data/coconut_ids.sdf \
        --id-prop  identifier
"""
from __future__ import print_function
__BUILD__ = "2026-11-12d"
import argparse, gzip, os, re, sys

TAG = re.compile(r"^>\s+<([^>]+)>")
RECNO = re.compile(r"^REC(\d{8})$")


def openf(p, mode="rt"):
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)


def records(path):
    """Yield each SDF record as a list of lines, including the $$$$ line."""
    buf = []
    with openf(path) as fh:
        for line in fh:
            buf.append(line)
            if line.startswith("$$$$"):
                yield buf
                buf = []
    if buf and any(x.strip() for x in buf):
        yield buf


def list_fields(path, limit=50):
    print("\n%s" % path)
    if not os.path.exists(path):
        print("  ! not found"); return
    names, titles, n = {}, [], 0
    for rec in records(path):
        n += 1
        title = rec[0].rstrip("\n") if rec else ""
        if len(titles) < 5:
            titles.append(title)
        for i, line in enumerate(rec):
            m = TAG.match(line)
            if m:
                key = m.group(1)
                val = rec[i + 1].rstrip("\n") if i + 1 < len(rec) else ""
                if key not in names:
                    names[key] = val
        if n >= limit:
            break
    print("  records sampled : %d" % n)
    print("  title lines     : %s"
          % ("  |  ".join(repr(t) for t in titles) if titles else "(none)"))
    if not names:
        print("  ! no property tags found in the sampled records")
        return
    print("  property fields and first value seen:")
    w = max(len(k) for k in names)
    for k in sorted(names):
        print("      %-*s  %s" % (w, k, names[k][:70]))
    print()
    for cand in ("identifier", "coconut_id", "inchikey", "standard_inchi_key"):
        for k in names:
            if k.lower() == cand:
                print("  --> likely id field: %s" % k)
                return
    print("  --> no obvious id field; choose one from the list above")


def build_map(path, id_prop, use_title):
    """index (1-based) -> identifier, for every record in the source."""
    out = [None]
    n = miss = 0
    want = "> <%s>" % id_prop
    for rec in records(path):
        n += 1
        val = ""
        if use_title:
            val = rec[0].strip()
        else:
            for i, line in enumerate(rec):
                m = TAG.match(line)
                if m and m.group(1) == id_prop:
                    val = rec[i + 1].strip() if i + 1 < len(rec) else ""
                    break
        if not val:
            miss += 1
        out.append(val)
        if n % 200000 == 0:
            sys.stderr.write("  ...%d source records scanned\n" % n)
    return out, n, miss


def main():
    ap = argparse.ArgumentParser(description=__doc__,
          formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list-fields", metavar="SDF")
    ap.add_argument("--source")
    ap.add_argument("--in", dest="inp")
    ap.add_argument("--out")
    ap.add_argument("--id-prop", default="identifier")
    ap.add_argument("--from-title", action="store_true",
                    help="take the identifier from the record title line "
                         "instead of a property tag")
    ap.add_argument("--map-out", default=None,
                    help="also write a TSV of placeholder to real identifier")
    a = ap.parse_args()

    if a.list_fields:
        list_fields(a.list_fields)
        return
    for need, name in ((a.source, "--source"), (a.inp, "--in"), (a.out, "--out")):
        if not need:
            sys.exit("ERROR: %s is required (or use --list-fields)" % name)

    print("Scanning source for identifiers")
    idmap, nsrc, miss = build_map(a.source, a.id_prop, a.from_title)
    print("  source records     : %d" % nsrc)
    print("  without identifier : %d" % miss)
    if miss == nsrc:
        sys.exit("\nERROR: no record carried %r. Run --list-fields on the source "
                 "to see the real field names." % a.id_prop)

    d = os.path.dirname(os.path.abspath(a.out))
    if d and not os.path.isdir(d):
        os.makedirs(d)

    mapfh = open(a.map_out, "w") if a.map_out else None
    if mapfh:
        mapfh.write("placeholder\tidentifier\n")

    print("\nRewriting %s" % a.inp)
    n = fixed = already = unresolved = 0
    with open(a.out, "w") as out:
        for rec in records(a.inp):
            n += 1
            # locate the deck_id tag and its value line
            vi = None
            for i, line in enumerate(rec):
                m = TAG.match(line)
                if m and m.group(1) == "deck_id":
                    vi = i + 1
                    break
            if vi is None or vi >= len(rec):
                already += 1
                out.writelines(rec)
                continue
            cur = rec[vi].strip()
            mm = RECNO.match(cur)
            if not mm:
                already += 1
                out.writelines(rec)
                continue
            idx = int(mm.group(1))
            real = idmap[idx] if 0 < idx < len(idmap) else ""
            if not real:
                unresolved += 1
                out.writelines(rec)
                continue
            rec[vi] = real + "\n"
            if rec[0].strip() == cur:
                rec[0] = real + "\n"
            out.writelines(rec)
            fixed += 1
            if mapfh:
                mapfh.write("%s\t%s\n" % (cur, real))
            if n % 100000 == 0:
                sys.stderr.write("  ...%d rewritten, %d fixed\n" % (n, fixed))
    if mapfh:
        mapfh.close()

    print("\nRead                 : %d" % n)
    print("Identifiers restored : %d" % fixed)
    print("Already had an id    : %d" % already)
    if unresolved:
        print("Could not resolve    : %d" % unresolved)
    print("Written              : %s" % a.out)
    if a.map_out:
        print("Mapping              : %s" % a.map_out)
    if fixed == 0:
        print("\n! Nothing was changed. Either the ids were already real, or the"
              "\n! placeholders do not match the REC######## pattern.")


if __name__ == "__main__":
    main()
