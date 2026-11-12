#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Turn the NPASS release into one CSV of anticancer compounds.

NPASS splits structures, activity records and target descriptions across
separate files, so the anticancer subset has to be assembled by joining them.
This selects activity records whose target is a cancer cell line, applies a
potency threshold, joins the surviving compound ids to their structures, and
writes a single CSV that 10_build_from_sources.py can read.

    # see what the columns are called in your release
    python scripts/11_npass_join.py --inspect data/raw/npass_raw/*.txt

    # build the anticancer subset
    python scripts/11_npass_join.py \
        --activities data/raw/npass_raw/npass_activities.txt \
        --structures data/raw/npass_raw/npass_structure.txt \
        --targets    data/raw/npass_raw/npass_target.txt \
        --out        data/raw/npass_anticancer.csv

Column names are auto-detected and every one can be overridden. Units are
normalised, so a threshold given in micromolar also matches records stored in
nanomolar or molar.
"""
from __future__ import print_function
__BUILD__ = "2026-11-12d"
import argparse, csv, gzip, os, re, sys

ID_N   = ("np_id", "npaid", "compound_id", "id")
SMI_N  = ("canonical_smiles", "smiles", "isomeric_smiles", "structure")
INCHI_N= ("standard_inchi", "inchi")
KEY_N  = ("standard_inchi_key", "inchikey", "inchi_key")
TID_N  = ("target_id", "tid", "uniprot_id")
TNAME_N= ("target_name", "name", "target", "pref_name", "cell_line")
TTYPE_N= ("target_type", "type", "target_class", "organism")
ATYPE_N= ("activity_type", "act_type", "standard_type", "type")
AVAL_N = ("activity_value", "act_value", "standard_value", "value")
AUNIT_N= ("activity_units", "activity_unit", "act_unit", "act_units",
          "standard_units", "standard_unit", "unit", "units", "uom",
          "activity_value_unit", "value_unit")

# Common human cancer cell lines, plus generic markers. Used only when the
# release has no usable target-type column to filter on.
CANCER_HINTS = re.compile(
    r"\b(mcf-?7|mda-?mb|a549|hela|hct-?116|ht-?29|sw480|sw620|caco-?2|hepg2|"
    r"huh-?7|pc-?3|du-?145|lncap|skov-?3|ovcar|a2780|u-?87|u-?251|sh-?sy5y|"
    r"k562|hl-?60|jurkat|raji|molt-?4|thp-?1|a375|sk-?mel|b16|4t1|ct26|"
    r"panc-?1|miapaca|bxpc-?3|h1299|h460|nci-?h\d|calu-?3|t47d|bt-?474|"
    r"skbr-?3|saos-?2|u2os|hl-?60|kb\b|hep-?2|siha|c33a|ishikawa|"
    r"carcinoma|sarcoma|leukemi|lymphoma|melanoma|glioma|glioblastoma|"
    r"myeloma|hepatoma|neuroblastoma|adenocarcinoma|tumou?r|cancer|"
    r"cell\s*line)\b", re.I)

UNIT_TO_UM = {"m": 1e6, "mm": 1e3, "um": 1.0, "µm": 1.0, "μm": 1.0,
              "nm": 1e-3, "pm": 1e-6, "ug/ml": None, "µg/ml": None,
              "mg/ml": None, "ng/ml": None, "%": None, "": None}


def openf(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")



# --- delimiter detection (inlined on purpose) ------------------------------
# Counting characters over a sample is unreliable for chemical data. An InChI
# string carries commas in its /c, /h and /t layers, dozens per record, so a
# tab separated file can hold far more commas than tabs and a naive count
# picks the wrong delimiter: the whole header then parses as one column and
# every field lookup fails. This judges candidates by the structure they
# produce instead. It is repeated in each script rather than imported, so that
# a script still works when copied out of this directory on its own.
def sniff_delimiter(path, probe_lines=25, default="\t"):
    import gzip as _gz
    op = _gz.open if path.endswith(".gz") else open
    lines = []
    try:
        fh = op(path, "rt", errors="replace")
    except TypeError:
        fh = op(path, "rt")
    try:
        for _ in range(probe_lines):
            ln = fh.readline()
            if not ln:
                break
            ln = ln.rstrip("\r\n")
            if ln.strip():
                lines.append(ln)
    finally:
        fh.close()
    if not lines:
        return default
    header, body = lines[0], lines[1:]
    best, best_score = None, -1.0
    for d in ("\t", ",", ";", "|"):
        nfields = header.count(d) + 1
        if nfields < 2:
            continue
        if body:
            same = sum(1 for ln in body if ln.count(d) + 1 == nfields)
            consistency = float(same) / len(body)
        else:
            consistency = 1.0
        score = consistency * 1000.0 + min(nfields, 50)
        if score > best_score:
            best, best_score = d, score
    return best if best is not None else default
# ---------------------------------------------------------------------------


def delim_of(p):
    return sniff_delimiter(p)


def rows(p):
    csv.field_size_limit(10 ** 7)
    with openf(p) as fh:
        for r in csv.DictReader(fh, delimiter=delim_of(p)):
            yield r


def cols_of(p):
    with openf(p) as fh:
        return next(csv.reader(fh, delimiter=delim_of(p)), [])


def pick(cols, names, override=None):
    low = dict((c.lower().strip(), c) for c in cols)
    if override:
        return override if override in cols else low.get(override.lower())
    for n in names:
        if n in low:
            return low[n]
    return None


def to_um(val, unit):
    """Convert an activity value to micromolar, or None if not possible."""
    try:
        v = float(str(val).strip().lstrip("<>=~ "))
    except (ValueError, TypeError):
        return None
    u = str(unit or "").strip().lower().replace(" ", "")
    if u in UNIT_TO_UM:
        f = UNIT_TO_UM[u]
        return None if f is None else v * f
    return None


def inspect(paths):
    for p in paths:
        print("\n%s" % p)
        if not os.path.exists(p):
            print("  ! not found"); continue
        c = cols_of(p)
        print("  columns: %s" % ", ".join(c))
        first = next(rows(p), None)
        if first:
            for k in c[:10]:
                print("      %-24s %s" % (k, str(first.get(k, ""))[:58]))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
          formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inspect", nargs="+", metavar="FILE")
    ap.add_argument("--activities")
    ap.add_argument("--structures")
    ap.add_argument("--targets", default=None,
                    help="target description file, used to tell which target "
                         "ids are cancer cell lines")
    ap.add_argument("--out")
    ap.add_argument("--max-um", type=float, default=10.0,
                    help="keep activities at or below this, in micromolar "
                         "(default 10). Records whose units cannot be "
                         "converted are kept and flagged.")
    ap.add_argument("--drop-unconverted", action="store_true",
                    help="discard records whose units cannot be converted to "
                         "micromolar, such as ug/mL or per cent inhibition. "
                         "Gives a cleaner potency criterion to state, at the "
                         "cost of losing real actives reported in mass units.")
    ap.add_argument("--activity-types", default="IC50,GI50,EC50,CC50,LC50,ED50",
                    help="comma separated, case insensitive. Empty means any.")
    ap.add_argument("--np-col"); ap.add_argument("--smiles-col")
    ap.add_argument("--target-col"); ap.add_argument("--target-name-col")
    ap.add_argument("--target-type-col")
    ap.add_argument("--unit-col", help="name of the activity unit column")
    ap.add_argument("--value-col", help="name of the activity value column")
    a = ap.parse_args()

    if a.inspect:
        inspect(a.inspect); return 0
    for v, n in ((a.activities, "--activities"), (a.structures, "--structures"),
                 (a.out, "--out")):
        if not v:
            sys.exit("ERROR: %s is required (or use --inspect)" % n)

    # ---------------------------------------------------------- targets
    cancer_tids, tname = set(), {}
    if a.targets and os.path.exists(a.targets):
        tc = cols_of(a.targets)
        tid = pick(tc, TID_N, a.target_col)
        tnm = pick(tc, TNAME_N, a.target_name_col)
        tty = pick(tc, TTYPE_N, a.target_type_col)
        print("Targets   : id=%s name=%s type=%s" % (tid, tnm, tty))
        if not tid:
            print("  ! no target id column, falling back to name matching")
        else:
            n = 0
            for r in rows(a.targets):
                n += 1
                k = (r.get(tid) or "").strip()
                nm = (r.get(tnm) or "") if tnm else ""
                ty = (r.get(tty) or "") if tty else ""
                if k:
                    tname[k] = nm
                    if CANCER_HINTS.search(nm) or CANCER_HINTS.search(ty):
                        cancer_tids.add(k)
            print("  %d targets, %d look like cancer cell lines" % (n, len(cancer_tids)))
    else:
        print("Targets   : no file given, matching on the activity table alone")

    # ---------------------------------------------------------- activities
    ac = cols_of(a.activities)
    aid = pick(ac, ID_N, a.np_col)
    atid = pick(ac, TID_N, a.target_col)
    atnm = pick(ac, TNAME_N, a.target_name_col)
    atty = pick(ac, TTYPE_N, a.target_type_col)
    aty = pick(ac, ATYPE_N)
    ava = pick(ac, AVAL_N, a.value_col)
    aun = pick(ac, AUNIT_N, a.unit_col)
    print("Activities: np=%s target=%s type=%s value=%s unit=%s"
          % (aid, atid or atnm, aty, ava, aun))
    if ava and not aun:
        print("  ! No unit column was found, so no value can be converted to")
        print("  ! micromolar and the --max-um threshold cannot be applied.")
        print("  ! Every matching record will be kept and flagged instead.")
        print("  ! Columns present: %s" % ", ".join(ac))
        print("  ! Name the right one with --unit-col, or accept the unfiltered")
        print("  ! set and say so in the methods section.")
    if not aid:
        sys.exit("ERROR: no compound id column in the activity file.\n"
                 "       columns: %s\n       use --np-col" % ", ".join(ac))

    want_types = set(t.strip().lower() for t in a.activity_types.split(",") if t.strip())
    keep, n, n_cancer, n_pot, n_nounit = {}, 0, 0, 0, 0
    for r in rows(a.activities):
        n += 1
        npid = (r.get(aid) or "").strip()
        if not npid:
            continue
        # is the target a cancer cell line?
        is_c = False
        if cancer_tids and atid:
            is_c = (r.get(atid) or "").strip() in cancer_tids
        if not is_c:
            blob = " ".join(str(r.get(c) or "") for c in (atnm, atty) if c)
            if blob and CANCER_HINTS.search(blob):
                is_c = True
        if not is_c:
            continue
        n_cancer += 1
        if want_types and aty:
            if (r.get(aty) or "").strip().lower() not in want_types:
                continue
        um = to_um(r.get(ava), r.get(aun)) if ava else None
        if um is None:
            if a.drop_unconverted:
                continue
            n_nounit += 1
            flag = "unconverted"
        elif um <= a.max_um:
            n_pot += 1
            flag = "%.4g" % um
        else:
            continue
        prev = keep.get(npid)
        if prev is None or (flag != "unconverted" and
                            (prev[0] == "unconverted" or float(flag) < float(prev[0]))):
            tn = (r.get(atnm) or "") if atnm else tname.get(
                (r.get(atid) or "").strip(), "")
            keep[npid] = (flag, tn)
        if n % 200000 == 0:
            sys.stderr.write("  ...%d activity rows\n" % n)

    print("  %d rows, %d against cancer cell lines, %d within %g uM, "
          "%d kept with unconvertible units"
          % (n, n_cancer, n_pot, a.max_um, n_nounit))
    print("  %d unique compounds" % len(keep))
    if not keep:
        sys.exit("\nERROR: nothing matched. Run --inspect on the activity and "
                 "target files\n       and set the column options explicitly.")

    # ---------------------------------------------------------- structures
    sc = cols_of(a.structures)
    sid = pick(sc, ID_N, a.np_col)
    ssm = pick(sc, SMI_N, a.smiles_col)
    sin = pick(sc, INCHI_N)
    skey = pick(sc, KEY_N)
    print("Structures: np=%s smiles=%s inchi=%s" % (sid, ssm, sin))
    if not sid or not (ssm or sin):
        sys.exit("ERROR: structure file needs a compound id and a SMILES or "
                 "InChI column.\n       columns: %s" % ", ".join(sc))

    d = os.path.dirname(os.path.abspath(a.out))
    if d and not os.path.isdir(d):
        os.makedirs(d)
    written, missing = 0, 0
    with open(a.out, "w") as fh:
        w = csv.writer(fh)
        w.writerow(["np_id", "smiles", "standard_inchi", "standard_inchi_key",
                    "best_activity_um", "cancer_cell_line"])
        seen = set()
        for r in rows(a.structures):
            k = (r.get(sid) or "").strip()
            if k not in keep or k in seen:
                continue
            seen.add(k)
            smi = (r.get(ssm) or "").strip() if ssm else ""
            inc = (r.get(sin) or "").strip() if sin else ""
            key = (r.get(skey) or "").strip() if skey else ""
            if not smi and not inc:
                missing += 1
                continue
            act, cl = keep[k]
            w.writerow([k, smi, inc, key, act, cl])
            written += 1
    missing += len(keep) - written - missing

    print("\nWritten   : %d compounds -> %s" % (written, a.out))
    if len(keep) - written > 0:
        print("  %d compounds had an activity record but no structure"
              % (len(keep) - written))
    print("\nNext:")
    print("  python scripts/10_build_from_sources.py --auto data/raw \\")
    print("      --out data/anticancer_raw.sdf --report logs/source_overlap.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
