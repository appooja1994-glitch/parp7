#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build one SDF from several anticancer natural product source files.

Each source may be a CSV, TSV or SDF, with or without a header, and may carry
structures as SMILES, InChI or InChIKey. Columns are auto-detected and can be
overridden. Compounds are deduplicated across sources on InChIKey and every
compound records which sources it came from.

Always inspect a file before using it:

    python scripts/10_build_from_sources.py --inspect data/raw/npcare.csv

Then build:

    python scripts/10_build_from_sources.py \
        --source npcare:data/raw/npcare.csv \
        --source npact:data/raw/npact.csv \
        --source cancerhsp:data/raw/cancerhsp.csv \
        --out data/anticancer_raw.sdf \
        --report logs/source_overlap.tsv

No three dimensional coordinates are generated here. Nothing in the curation
step needs them and 04_prepare_pdbqt.sh builds the geometry it actually docks,
so generating them twice would only waste time.

A source may be narrowed with a filter, which is useful for an activity table:

    --source npass:data/raw/npass_act.tsv --filter "activity_value<=10"
    --source npass:data/raw/npass_act.tsv --filter "target_type~cell"

A filter applies to the source that precedes it. Repeat for several filters.
"""
from __future__ import print_function
__BUILD__ = "2026-11-12d"
import argparse, csv, gzip, os, re, sys

try:
    from rdkit import Chem, RDLogger
except ImportError:
    sys.exit("ERROR: RDKit not importable by %s" % sys.executable)
RDLogger.DisableLog("rdApp.*")

KEY_NAMES = ("standard_inchi_key", "inchikey", "inchi_key", "standard_inchikey",
             "inchi key")
SMI_NAMES = ("canonical_smiles", "smiles", "smile", "isomeric_smiles",
             "canonical smiles", "structure", "smiles_structure")
INCHI_NAMES = ("standard_inchi", "inchi")
ID_NAMES = ("np_id", "npaid", "identifier", "compound_id", "id", "cid",
            "coconut_id", "molecule_name", "molecule name", "molecule",
            "compound_name", "compound name", "chemical_name", "chemical name",
            "ingredient_name", "ingredient", "mol_name", "compound", "name",
            "title")

FILT = re.compile(r"^\s*([^<>=~!]+?)\s*(<=|>=|!=|~|=|<|>)\s*(.+?)\s*$")


def openf(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")


def pick(cols, names, override=None):
    low = dict((c.lower().strip(), c) for c in cols)
    if override:
        return override if override in cols else low.get(override.lower())
    for n in names:
        if n in low:
            return low[n]
    return None



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


def sniff(path):
    return sniff_delimiter(path)


def read_xlsx(path):
    try:
        import pandas as pd
    except ImportError:
        sys.exit("ERROR: reading .xlsx needs pandas and openpyxl\n"
                 "       pip install pandas openpyxl")
    try:
        df = pd.read_excel(path, dtype=str)
    except ImportError:
        sys.exit("ERROR: reading .xlsx needs openpyxl\n"
                 "       pip install openpyxl")
    df = df.where(df.notna(), "")
    for rec in df.to_dict("records"):
        yield dict((str(k), "" if v is None else str(v)) for k, v in rec.items())


def read_rows(path):
    """Yield dict rows from a delimited file, an Excel sheet, or an SDF."""
    if path.endswith((".xlsx", ".xlsm")):
        for r in read_xlsx(path):
            yield r
        return
    if path.endswith((".sdf", ".sdf.gz")):
        fh = gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")
        for mol in Chem.ForwardSDMolSupplier(fh, sanitize=False):
            if mol is None:
                continue
            d = dict(mol.GetPropsAsDict())
            if mol.HasProp("_Name"):
                d.setdefault("_Name", mol.GetProp("_Name"))
            try:
                d.setdefault("__smiles__", Chem.MolToSmiles(mol))
            except Exception:
                pass
            yield d
        fh.close()
    else:
        csv.field_size_limit(10 ** 7)
        with openf(path) as fh:
            for row in csv.DictReader(fh, delimiter=sniff(path)):
                yield row


def columns_of(path, limit=300):
    if path.endswith((".xlsx", ".xlsm")):
        first = next(read_xlsx(path), {})
        return list(first.keys())
    if path.endswith((".sdf", ".sdf.gz")):
        names, n = set(), 0
        for d in read_rows(path):
            names.update(d.keys()); n += 1
            if n >= limit:
                break
        return sorted(x for x in names if x != "__smiles__")
    with openf(path) as fh:
        return next(csv.reader(fh, delimiter=sniff(path)), [])


def inspect(path):
    print("\n%s" % path)
    if not os.path.exists(path):
        print("  ! not found"); return
    cols = columns_of(path)
    first = next(read_rows(path), {}) or {}
    print("  columns : %s" % (", ".join(cols) if cols else "(none)"))
    print("  --> smiles   : %s" % (pick(cols, SMI_NAMES) or "NOT FOUND"))
    print("  --> inchi    : %s" % (pick(cols, INCHI_NAMES) or "NOT FOUND"))
    print("  --> inchikey : %s" % (pick(cols, KEY_NAMES) or "NOT FOUND"))
    print("  --> id       : %s" % (pick(cols, ID_NAMES) or "NOT FOUND"))
    if first:
        print("  first row:")
        for c in cols[:12]:
            print("      %-26s %s" % (c, str(first.get(c, ""))[:62]))
    if path.endswith((".sdf", ".sdf.gz")) and "__smiles__" in first:
        print("      %-26s %s" % "(structure read from the SDF itself)".ljust(26, " ")
              if False else "      (structures come from the SDF records themselves)")


SOURCE_EXT = (".csv", ".tsv", ".txt", ".sdf", ".sdf.gz", ".csv.gz",
              ".tsv.gz", ".txt.gz", ".xlsx", ".xlsm")


def usable(path):
    """(ok, reason, n_probe_hits). Can structures actually be read from this?"""
    if path.lower().endswith((".zip", ".gz")) and not path.endswith(
            (".sdf.gz", ".csv.gz", ".tsv.gz", ".txt.gz")):
        return False, "archive, unzip it first", 0
    try:
        cols = columns_of(path)
    except SystemExit:
        raise
    except Exception as e:
        return False, "unreadable (%s)" % type(e).__name__, 0
    is_sdf = path.endswith((".sdf", ".sdf.gz"))
    sc = pick(cols, SMI_NAMES)
    ic = pick(cols, INCHI_NAMES)
    if not (sc or ic or is_sdf):
        return False, "no SMILES or InChI column", 0
    hits = 0
    for i, row in enumerate(read_rows(path)):
        if i >= 25:
            break
        if mol_from(row, None, sc, ic) is not None:
            hits += 1
    if hits == 0:
        return False, "structure column found but nothing parsed", 0
    return True, "ok", hits


def scan_dir(d):
    """Candidate files directly in d. Subfolders are deliberately not
    searched, so a file you want to pass explicitly with --filter can simply
    be kept in one."""
    found = []
    for n in sorted(os.listdir(d)):
        if n.startswith("."):
            continue
        p = os.path.join(d, n)
        if not os.path.isfile(p):
            continue
        if p.lower().endswith(SOURCE_EXT) or p.lower().endswith(".zip"):
            found.append(p)
    return found


def parse_filter(spec):
    m = FILT.match(spec)
    if not m:
        sys.exit("ERROR: cannot parse --filter %r. Use col<=value, col=value, "
                 "col~regex and so on." % spec)
    return m.group(1), m.group(2), m.group(3)


def passes(row, filters):
    for col, op, val in filters:
        got = row.get(col)
        if got is None:
            for k in row:
                if k.lower().strip() == col.lower():
                    got = row[k]; break
        got = "" if got is None else str(got).strip()
        if op == "~":
            if not re.search(val, got, re.I):
                return False
        elif op == "=":
            if got.lower() != val.lower():
                return False
        elif op == "!=":
            if got.lower() == val.lower():
                return False
        else:
            try:
                g, v = float(got), float(val)
            except ValueError:
                return False
            if op == "<"  and not g <  v: return False
            if op == "<=" and not g <= v: return False
            if op == ">"  and not g >  v: return False
            if op == ">=" and not g >= v: return False
    return True


def mol_from(row, kc, sc, ic):
    smi = (row.get(sc) or "").strip() if sc else ""
    if not smi:
        smi = (row.get("__smiles__") or "").strip()
    if smi:
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            return m
    inchi = (row.get(ic) or "").strip() if ic else ""
    if inchi:
        try:
            m = Chem.MolFromInchi(inchi)
            if m is not None:
                return m
        except Exception:
            pass
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
          formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inspect", nargs="+", metavar="FILE")
    ap.add_argument("--auto", metavar="DIR",
                    help="use every usable file in this directory, each named "
                         "after its filename. No column names needed.")
    ap.add_argument("--source", action="append", default=[], metavar="NAME:PATH")
    ap.add_argument("--filter", action="append", default=[], metavar="EXPR",
                    help="narrow the preceding --source")
    ap.add_argument("--smiles-col", action="append", default=[])
    ap.add_argument("--id-col", action="append", default=[])
    ap.add_argument("--out")
    ap.add_argument("--report", default=None, help="TSV of per source counts")
    a = ap.parse_args()

    if a.inspect:
        for p in a.inspect:
            inspect(p)
        return
    if a.auto:
        if not os.path.isdir(a.auto):
            sys.exit("ERROR: not a directory: %s" % a.auto)
        print("Scanning %s" % a.auto)
        cand = scan_dir(a.auto)
        if not cand:
            sys.exit("ERROR: no candidate files found in %s\n"
                     "       expected .csv .tsv .txt .sdf .xlsx" % a.auto)
        chosen = []
        for p in cand:
            ok, why, hits = usable(p)
            rel = os.path.relpath(p, a.auto)
            if ok:
                print("  use  %-44s (%d/25 probe rows parsed)" % (rel, hits))
                stem = os.path.splitext(os.path.basename(p))[0]
                stem = re.sub(r"[^A-Za-z0-9]+", "_", stem).strip("_").lower()
                chosen.append("%s:%s" % (stem or "src", p))
            else:
                print("  skip %-44s %s" % (rel, why))
        if not chosen:
            sys.exit("\nERROR: none of the files could be read. Run --inspect "
                     "on them\n       to see their columns.")
        a.source = chosen + list(a.source)
        print("\n%d source file(s) selected\n" % len(chosen))

    if not a.source or not a.out:
        sys.exit("ERROR: need --source, --auto or --inspect, plus --out")

    # attach each --filter to the --source that precedes it on the command line
    order, si, fi = [], 0, 0
    for tok in sys.argv[1:]:
        if tok == "--source":
            order.append(["s", si]); si += 1
        elif tok == "--filter":
            order.append(["f", fi]); fi += 1
    n_auto = len(a.source) - si
    per_src_filters = dict((i, []) for i in range(len(a.source)))
    cur = None
    for kind, i in order:
        if kind == "s":
            cur = i + n_auto          # explicit sources sit after the auto ones
        elif cur is not None:
            per_src_filters[cur].append(parse_filter(a.filter[i]))

    index, order_keys, stats = {}, [], []
    for i, spec in enumerate(a.source):
        if ":" not in spec:
            sys.exit("ERROR: --source needs NAME:PATH, got %r" % spec)
        name, path = spec.split(":", 1)
        print("\n%s  (%s)" % (name, path))
        if not os.path.exists(path):
            print("  ! not found, skipping"); continue
        filters = per_src_filters.get(i, [])
        if filters:
            print("  filters: %s" % "  ".join("%s%s%s" % f for f in filters))

        cols = columns_of(path)
        sc = pick(cols, SMI_NAMES, a.smiles_col[i] if i < len(a.smiles_col) else None)
        ic = pick(cols, INCHI_NAMES)
        kc = pick(cols, KEY_NAMES)
        idc = pick(cols, ID_NAMES, a.id_col[i] if i < len(a.id_col) else None)
        is_sdf = path.endswith((".sdf", ".sdf.gz"))
        if not (sc or ic or is_sdf):
            print("  ! no SMILES or InChI column found. Columns: %s"
                  % ", ".join(cols))
            print("  ! pass --smiles-col after this --source to name it.")
            continue

        n = kept = bad = dup = filtered = 0
        for row in read_rows(path):
            n += 1
            if filters and not passes(row, filters):
                filtered += 1
                continue
            m = mol_from(row, kc, sc, ic)
            if m is None:
                bad += 1
                continue
            try:
                Chem.SanitizeMol(m)
                key = Chem.MolToInchiKey(m)
            except Exception:
                bad += 1
                continue
            if not key:
                bad += 1
                continue
            if key in index:
                if name not in index[key][1]:
                    index[key][1].append(name)
                dup += 1
                continue
            # carry a few informative columns through to the SDF when the
            # source has them, so the deposited annotation table is not empty
            carried = {}
            for want in ("best_activity_um", "cancer_cell_line",
                         "activity_um", "cell_line"):
                for k in row:
                    if k and k.lower().strip() == want:
                        v = row.get(k)
                        if v not in (None, ""):
                            # NB: do not name this `key`, which already holds
                            # the InChIKey in the enclosing scope
                            prop_name = ("best_activity_um" if "activ" in want
                                         else "cancer_cell_line")
                            carried.setdefault(prop_name, str(v).strip())
                        break
            cid = (str(row.get(idc)).strip() if idc and row.get(idc) else "")
            if not cid:
                cid = "%s_%06d" % (name, kept + 1)
            m.SetProp("_Name", cid)
            m.SetProp("src_id", cid)
            for k, v in carried.items():
                m.SetProp(k, v)
            index[key] = [m, [name]]
            order_keys.append(key)
            kept += 1
            if n % 50000 == 0:
                sys.stderr.write("    ...%d rows\n" % n)
        print("  rows %d, new %d, already seen %d, unusable %d%s"
              % (n, kept, dup, bad,
                 (", filtered out %d" % filtered) if filters else ""))
        stats.append((name, n, kept, dup, bad, filtered))

    d = os.path.dirname(os.path.abspath(a.out))
    if d and not os.path.isdir(d):
        os.makedirs(d)
    w = Chem.SDWriter(a.out)
    for key in order_keys:
        m, srcs = index[key]
        m.SetProp("inchikey", key)
        m.SetProp("sources", ",".join(srcs))
        m.SetProp("n_sources", str(len(srcs)))
        m.SetProp("anticancer_annotated", "1")
        w.write(m)
    w.close()

    if a.report:
        with open(a.report, "w") as fh:
            fh.write("source\trows\tnew\tduplicate\tunusable\tfiltered_out\n")
            for s in stats:
                fh.write("\t".join(str(x) for x in s) + "\n")

    multi = sum(1 for k in order_keys if len(index[k][1]) > 1)
    print("\nUnique compounds : %d" % len(order_keys))
    print("In 2 or more sources: %d" % multi)
    print("Written          : %s" % a.out)
    print("\nEvery compound carries anticancer_annotated=1 by construction, "
          "since each\nsource is an anticancer dataset. Curate next with "
          "01_filter.py.")


if __name__ == "__main__":
    main()
