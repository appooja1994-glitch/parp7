#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the environment, then prove the pipeline works by running it end to end
on built-in test molecules. Run this before touching real data.

    python scripts/00_preflight.py

Exits 0 only if every check passes. Anything that fails prints the exact
command that fixes it.
"""
from __future__ import print_function
__BUILD__ = "2026-11-12d"
import os, shutil, subprocess, sys, tempfile

OK, BAD, WARN = "  [ok]   ", "  [FAIL] ", "  [warn] "
fails, warns = [], []


def head(t):
    print("\n" + t)
    print("  " + "-" * (len(t) + 2))


def need(mod, fix, why=""):
    try:
        __import__(mod)
        try:
            v = sys.modules[mod].__version__
        except Exception:
            v = ""
        print(OK + mod + (" " + str(v) if v else ""))
        return True
    except Exception as e:
        print(BAD + "%s  (%s)" % (mod, type(e).__name__))
        if why:
            print("         " + why)
        print("         fix: " + fix)
        fails.append(mod)
        return False


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


def here_dir():
    return os.path.dirname(os.path.abspath(__file__))


def main():
    print("=" * 66)
    print("  Preflight check for the anticancer screening library pipeline")
    print("=" * 66)

    head("Python")
    print(OK + sys.version.split()[0] + "   " + sys.executable)
    if sys.version_info < (3, 8):
        print(BAD + "Python 3.8 or newer is required")
        fails.append("python")

    head("Core packages")
    need("rdkit", "conda install -c conda-forge rdkit")
    need("numpy", "pip install numpy")

    head("Meeko and its dependencies")
    print("  Meeko declares no dependencies of its own, so pip installs none")
    print("  of them. All five below are needed just to import it.")
    for m in ("numpy", "pandas", "scipy", "gemmi"):
        need(m, "pip install " + m)
    need("meeko", "pip install meeko gemmi scipy pandas")

    head("Optional")
    if shutil.which("obabel"):
        print(OK + "obabel  (fallback only, writes no torsion tree)")
    else:
        print(WARN + "obabel not found. Not needed if Meeko works.")
    try:
        import openpyxl  # noqa
        print(OK + "openpyxl  (lets the builder read .xlsx sources)")
    except Exception:
        print(WARN + "openpyxl missing. Needed only for .xlsx source files.")
        print("         fix: pip install openpyxl")
        warns.append("openpyxl")

    head("Disk space")
    try:
        st = os.statvfs(".")
        gb = st.f_bavail * st.f_frsize / 1e9
        if gb < 2:
            print(BAD + "%.1f GB free. Route A needs about 2 GB." % gb)
            fails.append("disk")
        else:
            print(OK + "%.1f GB free" % gb)
    except Exception:
        print(WARN + "could not determine free space")

    if fails:
        print("\n" + "=" * 66)
        print("  Environment is not ready. Fix the items above, then rerun.")
        print("=" * 66)
        return 1

    # ---------------------------------------------------------------- self test
    head("Script versions")
    import glob as _g, re as _re
    mine = __BUILD__
    builds, unstamped = {}, []
    for f in sorted(_g.glob(os.path.join(here_dir(), "*.py"))):
        txt = open(f).read()
        m = _re.search(r'^__BUILD__ = "([^"]+)"', txt, _re.M)
        if m:
            builds.setdefault(m.group(1), []).append(os.path.basename(f))
        else:
            unstamped.append(os.path.basename(f))
    if len(builds) == 1 and not unstamped:
        print(OK + "all %d scripts are from build %s"
              % (sum(len(v) for v in builds.values()), list(builds)[0]))
    else:
        print(BAD + "the scripts are not all from the same build")
        print("         this preflight is build %s" % mine)
        for b in sorted(builds):
            print("         build %-14s %s" % (b, ", ".join(builds[b])))
        if unstamped:
            print("         no build stamp : %s" % ", ".join(unstamped))
        print()
        print("         A partial extraction leaves old scripts in place and")
        print("         they fail in ways the new ones already fix. Re-extract")
        print("         the whole package over this directory:")
        print()
        print("             tar xzvf <package>.tar.gz --strip-components=1 -C .")
        print()
        print("         The v lists each file as it is written, so you can see")
        print("         that every script was replaced.")
        fails.append("build mismatch")

    head("Delimiter detection")
    if True:
        import tempfile as _tf
        inchi = ("InChI=1S/C42H72O14/c1-20(2)10-9-14-42(8,56-39-35(50)32(47)29(44)"
                 "23(17-43)54-39)21-11-15-40(6)26(21)22(45)18-24-41(7)16-13-25(46)"
                 "38(4,5)34(41)27(48)19-36(24,40)37(3)12-10-28(49)33(31)30(37)51-39/"
                 "h10,20-35,43-50H,9,11-19H2,1-8H3/t21-,22+,23+,24-,25-,26+,27+,28-/m0/s1")
        cases = [
            ("tab file whose InChI column is full of commas",
             "np_id\tInChI\tInChIKey\tSMILES\n" +
             "\n".join(["NPC1\t%s\tKEY-A\tCCO" % inchi] * 20), "\t"),
            ("plain comma csv", "a,b,c\n1,2,3\n4,5,6\n", ","),
            ("comma csv with quoted commas",
             'name,smiles\n"a, b",CCO\n"c, d",CCC\n', ","),
        ]
        for label, text, want in cases:
            f = _tf.NamedTemporaryFile("w", suffix=".txt", delete=False)
            f.write(text); f.close()
            got = sniff_delimiter(f.name)
            passed = got == want
            print((OK if passed else BAD) + label)
            if not passed:
                print("         got %r, expected %r" % (got, want))
                fails.append("delimiter")
            os.unlink(f.name)

    head("End to end self test on built-in molecules")
    here = os.path.dirname(os.path.abspath(__file__))
    tmp = tempfile.mkdtemp(prefix="preflight_")
    try:
        src = os.path.join(tmp, "src.csv")
        with open(src, "w") as fh:
            fh.write("compound,smiles,cell_line,ic50_um\n")
            for name, smi in [
                ("quercetin", "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
                ("rutin", "C[C@@H]1O[C@@H](OC[C@H]2O[C@@H](Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@H](O)[C@H]1O"),
                ("glucose", "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"),
                ("betulinic_acid", "CC(=C)[C@@H]1CC[C@]2(CC[C@]3(C)[C@H](CC[C@@H]4[C@@]5(C)CC[C@H](O)C(C)(C)[C@@H]5CC[C@@]34C)[C@@H]12)C(=O)O"),
                ("naproxen_sodium", "COc1ccc2cc([C@@H](C)C(=O)[O-])ccc2c1.[Na+]"),
                ("bad_row", "NOT_A_SMILES"),
            ]:
                fh.write("%s,%s,MCF-7,10.0\n" % (name, smi))

        raw = os.path.join(tmp, "raw.sdf")
        filt = os.path.join(tmp, "filt.sdf")
        pdbqt = os.path.join(tmp, "pdbqt")

        def run(args, label):
            p = subprocess.Popen([sys.executable] + args,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            out = p.communicate()[0].decode("utf8", "replace")
            if p.returncode != 0:
                print(BAD + label)
                print("\n".join("         " + l for l in out.splitlines()[-12:]))
                fails.append(label)
                return None
            return out

        if run([os.path.join(here, "10_build_from_sources.py"),
                "--source", "test:" + src, "--out", raw], "build") is None:
            return 1
        print(OK + "build")

        if run([os.path.join(here, "01_filter.py"), "--in", raw,
                "--out", filt, "--report-every", "0"], "filter") is None:
            return 1
        print(OK + "filter")

        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        got = {}
        for m in Chem.ForwardSDMolSupplier(open(filt, "rb")):
            if m is None:
                continue
            got[m.GetProp("deck_id")] = m

        checks = [
            ("malformed SMILES rejected", "bad_row" not in got),
            ("sugar-only molecule rejected", "glucose" not in got),
            ("quercetin retained", "quercetin" in got),
            ("quercetin flagged as PAINS",
             "quercetin" in got and got["quercetin"].GetProp("pains_flag") == "1"),
            ("betulinic acid retained despite high logP", "betulinic_acid" in got),
            ("rutin deglycosylated and merged onto quercetin",
             "rutin" not in got and "quercetin" in got),
            ("salt stripped from naproxen sodium",
             "naproxen_sodium" in got and
             all(a.GetSymbol() != "Na" for a in got["naproxen_sodium"].GetAtoms())),
            ("provenance preserved",
             "quercetin" in got and got["quercetin"].GetProp("sources") == "test"),
        ]
        for label, passed in checks:
            print((OK if passed else BAD) + label)
            if not passed:
                fails.append(label)

        out = run([os.path.join(here, "04_prepare_pdbqt.py"),
                   "--in", filt, "--out", pdbqt], "pdbqt")
        if out is None:
            return 1
        n = len([f for f in os.listdir(pdbqt) if f.endswith(".pdbqt")]) \
            if os.path.isdir(pdbqt) else 0
        want = len(got)   # one per molecule that survived filtering
        print((OK if n == want else BAD) +
              "pdbqt written (%d of %d surviving molecules)" % (n, want))
        if n != want:
            fails.append("pdbqt")
        else:
            f0 = os.path.join(pdbqt, sorted(os.listdir(pdbqt))[0])
            txt = open(f0).read()
            zs = set()
            for L in txt.splitlines():
                if L.startswith(("ATOM", "HETATM")):
                    try:
                        zs.add(round(float(L[46:54]), 3))
                    except ValueError:
                        pass
            print((OK if len(zs) > 2 else BAD) +
                  "coordinates are genuinely 3D (%d distinct z)" % len(zs))
            if len(zs) <= 2:
                fails.append("3d")
            print((OK if "ROOT" in txt else BAD) + "torsion tree present")
            if "ROOT" not in txt:
                fails.append("torsion")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("\n" + "=" * 66)
    if fails:
        print("  %d check(s) FAILED: %s" % (len(fails), ", ".join(map(str, fails))))
        print("=" * 66)
        return 1
    print("  All checks passed. The pipeline works on this machine.")
    if warns:
        print("  Optional items missing: %s" % ", ".join(warns))
    print("=" * 66)
    return 0


if __name__ == "__main__":
    sys.exit(main())
