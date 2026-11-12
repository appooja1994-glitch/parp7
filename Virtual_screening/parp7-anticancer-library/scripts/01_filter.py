#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filter a natural-product SDF into a docking-ready screening deck.

    python scripts/01_filter.py --in data/raw/<coconut>.sdf --out data/filtered.sdf

Steps, in order (each rejection is logged with its reason):
    1. parse + sanitize
    2. keep largest fragment, strip salts, neutralise
    3. element whitelist (organic only)
    4. deglycosylate: strip terminal circular sugars, keep the aglycone
    5. MW 150-600
    6. Lipinski: at most N violations (default 1)
    7. PAINS A/B/C
    8. deduplicate on InChIKey

Streams the input, so memory stays flat on a multi-GB file.
No f-strings, so it reports a wrong interpreter rather than dying on SyntaxError.
"""
from __future__ import print_function
__BUILD__ = "2026-11-12d"
import argparse
import gzip
import os
import sys
import csv

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors, SaltRemover
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
except ImportError:
    sys.exit("ERROR: RDKit not importable by %s\n"
             "       conda activate parp7 (or your env) first." % sys.executable)

RDLogger.DisableLog("rdApp.*")

ORGANIC = set([1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53])

# A circular sugar: 5-, 6- or 7-membered ring, exactly one ring oxygen,
# all other ring atoms sp3 carbon. Matches pyranose/furanose cores.
SUGAR_RINGS = [
    Chem.MolFromSmarts("[OX2;R1]1[CX4;R1][CX4;R1][CX4;R1][CX4;R1]1"),
    Chem.MolFromSmarts("[OX2;R1]1[CX4;R1][CX4;R1][CX4;R1][CX4;R1][CX4;R1]1"),
    Chem.MolFromSmarts("[OX2;R1]1[CX4;R1][CX4;R1][CX4;R1][CX4;R1][CX4;R1][CX4;R1]1"),
]
# Glycosidic / ester linkage from a ring-anomeric carbon out to the aglycone.
GLYCOSIDIC = Chem.MolFromSmarts("[CX4;R][OX2;!R][#6,#7,#16;!$([CX4;R])]")


def ring_is_sugar(mol, ring_atoms, min_exo_o=2):
    """A candidate ring is a sugar if it carries enough exocyclic oxygens."""
    exo = 0
    for idx in ring_atoms:
        a = mol.GetAtomWithIdx(idx)
        if a.GetAtomicNum() != 6:
            continue
        for nb in a.GetNeighbors():
            if nb.GetAtomicNum() == 8 and nb.GetIdx() not in ring_atoms:
                exo += 1
    return exo >= min_exo_o


def sugar_ring_atom_sets(mol, min_exo_o=2):
    """Atom-index sets of rings that look like sugars."""
    out = []
    seen = set()
    for patt in SUGAR_RINGS:
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            key = frozenset(match)
            if key in seen:
                continue
            seen.add(key)
            if ring_is_sugar(mol, set(match), min_exo_o):
                out.append(set(match))
    return out


def deglycosylate(mol, min_heavy_aglycone=6, max_rounds=6):
    """Iteratively cut terminal sugar rings off and keep the largest remainder.

    Returns (mol, n_sugars_removed). Returns the input unchanged if the molecule
    is essentially all sugar, which the caller then rejects.
    """
    removed = 0
    for _ in range(max_rounds):
        sugars = sugar_ring_atom_sets(mol)
        if not sugars:
            break
        sugar_atoms = set()
        for s in sugars:
            sugar_atoms |= s
        # candidate bonds: glycosidic O joining a sugar ring to a non-sugar part
        cut = []
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtom(), bond.GetEndAtom()
            ia, ib = a.GetIdx(), b.GetIdx()
            in_a, in_b = ia in sugar_atoms, ib in sugar_atoms
            if in_a == in_b:
                continue
            # the exocyclic partner must be an oxygen bridging out of the ring
            outer = b if in_a else a
            if outer.GetAtomicNum() == 8 and outer.GetDegree() == 2:
                cut.append(bond.GetIdx())
        if not cut:
            break
        # Cut with dummies so each fragment keeps explicit valences, then drop
        # the dummies (leaving implicit H) and sanitise BEFORE judging it.
        # Judging an unsanitised fragment silently mis-types its carbons, which
        # makes the sugar test fail and lets the sugar win the size comparison.
        frag_mol = Chem.FragmentOnBonds(mol, cut, addDummies=True)
        try:
            frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=False)
        except Exception:
            break
        if not frags:
            break

        clean = []
        for f in frags:
            try:
                rw = Chem.RWMol(f)
                for at in sorted([x.GetIdx() for x in rw.GetAtoms()
                                  if x.GetAtomicNum() == 0], reverse=True):
                    rw.RemoveAtom(at)
                g = rw.GetMol()
                Chem.SanitizeMol(g)
            except Exception:
                continue
            if g.GetNumHeavyAtoms() == 0:
                continue
            clean.append(g)
        if not clean:
            break

        # Prefer any non-sugar fragment, largest first. Size alone is wrong:
        # a hexose (11 heavy atoms) outweighs plenty of real aglycones.
        aglycones = [g for g in clean if not ring_is_all_sugar(g)]
        if not aglycones:
            return mol, removed          # everything is sugar - caller rejects
        keep = max(aglycones, key=lambda g: g.GetNumHeavyAtoms())
        if keep.GetNumHeavyAtoms() < min_heavy_aglycone:
            return mol, removed
        if keep.GetNumHeavyAtoms() == mol.GetNumHeavyAtoms():
            break                        # nothing actually came off; stop
        mol = keep
        removed += 1
    return mol, removed


def ring_is_all_sugar(frag):
    """True if the fragment is nothing but sugar: one or more sugar rings plus
    their hydroxyl and hydroxymethyl decoration, and no other skeleton.

    Grows a halo out from the sugar ring atoms, accepting an oxygen always and
    a carbon only when every one of its heavy neighbours is already in the halo
    or is an oxygen. That lets the halo run out along a -CH2-OH arm but stops it
    at the first carbon belonging to a real aglycone. The fragment is all sugar
    when the halo reaches every heavy atom.
    """
    try:
        sugars = sugar_ring_atom_sets(frag)
    except Exception:
        return False
    if not sugars:
        return False
    halo = set()
    for ring in sugars:
        halo |= ring

    changed = True
    while changed:
        changed = False
        for idx in list(halo):
            for nb in frag.GetAtomWithIdx(idx).GetNeighbors():
                j = nb.GetIdx()
                if j in halo:
                    continue
                z = nb.GetAtomicNum()
                if z == 8:
                    halo.add(j); changed = True
                elif z == 6:
                    ok = True
                    for nn in nb.GetNeighbors():
                        if nn.GetAtomicNum() == 1:
                            continue
                        if nn.GetIdx() in halo or nn.GetAtomicNum() == 8:
                            continue
                        ok = False; break
                    if ok:
                        halo.add(j); changed = True

    heavy = set(a.GetIdx() for a in frag.GetAtoms() if a.GetAtomicNum() > 1)
    return len(heavy - halo) == 0


def lipinski_violations(mol):
    v = 0
    if Descriptors.MolWt(mol) > 500:
        v += 1
    if Descriptors.MolLogP(mol) > 5:
        v += 1
    if rdMolDescriptors.CalcNumHBD(mol) > 5:
        v += 1
    if rdMolDescriptors.CalcNumHBA(mol) > 10:
        v += 1
    return v


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", required=True, help="input .sdf or .sdf.gz")
    ap.add_argument("--out", dest="out", required=True, help="output .sdf")
    ap.add_argument("--rejects", default=None, help="TSV of rejected ids + reason")
    ap.add_argument("--props", default=None, help="TSV of kept ids + properties")
    ap.add_argument("--mw-min", type=float, default=150.0)
    ap.add_argument("--mw-max", type=float, default=600.0)
    ap.add_argument("--max-lipinski", type=int, default=1)
    ap.add_argument("--keep-sugars", action="store_true",
                    help="skip deglycosylation entirely")
    ap.add_argument("--pains-mode", choices=("flag", "remove", "off"),
                    default="flag",
                    help="flag (default): keep PAINS hits but tag them, so you "
                         "can decide later; remove: drop them; off: skip the check. "
                         "Default is flag because PAINS A/B/C rejects a large "
                         "fraction of legitimate natural products (catechols, "
                         "quinones, many flavonoids).")
    ap.add_argument("--id-prop", default=None,
                    help="SDF property holding the compound id "
                         "(default: molecule title, then identifier/coconut_id)")
    ap.add_argument("--limit", type=int, default=0, help="stop after N records (testing)")
    ap.add_argument("--report-every", type=int, default=25000)
    a = ap.parse_args()

    for p in (a.out, a.rejects, a.props):
        if p:
            d = os.path.dirname(os.path.abspath(p))
            if d:
                try:
                    os.makedirs(d)
                except OSError:
                    pass

    params = FilterCatalogParams()
    for cat in (FilterCatalogParams.FilterCatalogs.PAINS_A,
                FilterCatalogParams.FilterCatalogs.PAINS_B,
                FilterCatalogParams.FilterCatalogs.PAINS_C):
        params.AddCatalog(cat)
    pains = FilterCatalog(params)

    remover = SaltRemover.SaltRemover()
    uncharger = rdMolStandardize.Uncharger()
    chooser = rdMolStandardize.LargestFragmentChooser()

    counts = {}
    def bump(k):
        counts[k] = counts.get(k, 0) + 1

    seen = set()
    writer = Chem.SDWriter(a.out)
    rej_fh = open(a.rejects, "w") if a.rejects else None
    rej = csv.writer(rej_fh, delimiter="\t") if rej_fh else None
    if rej:
        rej.writerow(["id", "reason"])
    prop_fh = open(a.props, "w") if a.props else None
    prp = csv.writer(prop_fh, delimiter="\t") if prop_fh else None
    if prp:
        prp.writerow(["id", "inchikey", "smiles", "mw", "clogp",
                      "hbd", "hba", "tpsa", "rotb", "rings",
                      "lipinski_viol", "sugars_removed", "pains_flag"])

    n = 0
    fh = open_maybe_gz(a.inp)
    supplier = Chem.ForwardSDMolSupplier(fh, sanitize=False, removeHs=False)

    for mol in supplier:
        n += 1
        if a.limit and n > a.limit:
            break
        if a.report_every and n % a.report_every == 0:
            sys.stderr.write("  ...read %d, kept %d\n" % (n, counts.get("kept", 0)))
            sys.stderr.flush()

        if mol is None:
            bump("unparseable"); 
            if rej: rej.writerow(["<unparseable record %d>" % n, "unparseable"])
            continue

        cid = ""
        if a.id_prop and mol.HasProp(a.id_prop):
            cid = mol.GetProp(a.id_prop)
        if not cid:
            cid = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        if not cid:
            for k in ("identifier", "coconut_id", "COCONUT_ID", "ID"):
                if mol.HasProp(k):
                    cid = mol.GetProp(k); break
        if not cid:
            cid = "REC%08d" % n

        try:
            Chem.SanitizeMol(mol)
        except Exception:
            bump("sanitize_fail")
            if rej: rej.writerow([cid, "sanitize_fail"])
            continue

        try:
            mol = remover.StripMol(mol, dontRemoveEverything=True)
            mol = chooser.choose(mol)
            mol = uncharger.uncharge(mol)
            Chem.SanitizeMol(mol)
        except Exception:
            bump("standardise_fail")
            if rej: rej.writerow([cid, "standardise_fail"])
            continue

        if mol.GetNumHeavyAtoms() == 0:
            bump("empty"); 
            if rej: rej.writerow([cid, "empty"])
            continue

        bad_el = False
        for at in mol.GetAtoms():
            if at.GetAtomicNum() not in ORGANIC:
                bad_el = True; break
        if bad_el:
            bump("non_organic_element")
            if rej: rej.writerow([cid, "non_organic_element"])
            continue

        n_sug = 0
        if not a.keep_sugars:
            try:
                mol2, n_sug = deglycosylate(mol)
            except Exception:
                mol2, n_sug = mol, 0
            if ring_is_all_sugar(mol2):
                bump("sugar_only")
                if rej: rej.writerow([cid, "sugar_only"])
                continue
            mol = mol2

        mw = Descriptors.MolWt(mol)
        if mw < a.mw_min or mw > a.mw_max:
            bump("mw_out_of_range")
            if rej: rej.writerow([cid, "mw_out_of_range"])
            continue

        viol = lipinski_violations(mol)
        if viol > a.max_lipinski:
            bump("lipinski")
            if rej: rej.writerow([cid, "lipinski_%d_violations" % viol])
            continue

        pains_hit = 0
        if a.pains_mode != "off":
            try:
                if pains.HasMatch(mol):
                    pains_hit = 1
            except Exception:
                pains_hit = 0
            if pains_hit and a.pains_mode == "remove":
                bump("pains")
                if rej: rej.writerow([cid, "pains"])
                continue
            if pains_hit:
                bump("pains_flagged")

        try:
            key = Chem.MolToInchiKey(mol)
        except Exception:
            key = ""
        if not key:
            bump("inchikey_fail")
            if rej: rej.writerow([cid, "inchikey_fail"])
            continue
        if key in seen:
            bump("duplicate")
            if rej: rej.writerow([cid, "duplicate"])
            continue
        seen.add(key)

        mol.SetProp("_Name", cid)
        mol.SetProp("deck_id", cid)
        mol.SetProp("inchikey", key)
        mol.SetProp("sugars_removed", str(n_sug))
        mol.SetProp("pains_flag", str(pains_hit))
        writer.write(mol)
        bump("kept")

        if prp:
            prp.writerow([cid, key, Chem.MolToSmiles(mol),
                          "%.2f" % mw, "%.2f" % Descriptors.MolLogP(mol),
                          rdMolDescriptors.CalcNumHBD(mol),
                          rdMolDescriptors.CalcNumHBA(mol),
                          "%.1f" % rdMolDescriptors.CalcTPSA(mol),
                          rdMolDescriptors.CalcNumRotatableBonds(mol),
                          rdMolDescriptors.CalcNumRings(mol),
                          viol, n_sug, pains_hit])

    writer.close()
    fh.close()
    if rej_fh: rej_fh.close()
    if prop_fh: prop_fh.close()

    print("\nRead      : %d" % n)
    for k in ("unparseable", "sanitize_fail", "standardise_fail", "empty",
              "non_organic_element", "sugar_only", "mw_out_of_range",
              "lipinski", "pains", "inchikey_fail", "duplicate"):
        if counts.get(k):
            print("  rejected %-22s %d" % (k, counts[k]))
    if counts.get("pains_flagged"):
        print("  (kept but PAINS-flagged: %d)" % counts["pains_flagged"])
    print("Kept      : %d" % counts.get("kept", 0))
    print("Written   : %s" % a.out)


if __name__ == "__main__":
    main()
