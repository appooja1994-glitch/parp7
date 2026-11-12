# Anticancer natural product screening library for PARP7

Reproducible construction of a docking-ready library of natural products with
reported anticancer activity, for structure-based virtual screening against a
comparative model of the PARP7 catalytic domain.

The library was built from **NPASS 3.0**, so every compound in it carries a
measured growth-inhibitory activity against a human cancer cell line, with the
cell line and the literature citation recorded alongside. `docs/` gives the
full procedure and the reasoning behind every filtering decision.

The pipeline can also read other anticancer compound lists and merge them, and
it can optionally build from the whole of COCONUT for breadth. Neither was used
for the deposited library.

## Repository layout

```
scripts/     the pipeline, each script runnable on its own
run_all.sh   runs the whole thing with a check after every stage
data/        downloads (untracked) and the deposited .smi files (tracked)
logs/        provenance, attrition and per-source counts
docs/        the procedure document, the workflow figure and the code that builds both
```

## Deposited compound sets

| File | Content |
|---|---|
| `data/total_molecules.smi` | Every compound assembled from the sources, before curation |
| `data/final_molecules.smi` | The curated screening library |
| `data/final_molecules_annotations.tsv` | Provenance and annotation per compound |

Both `.smi` files are `SMILES<TAB>identifier`, one compound per line.

## Citing the sources

The deposited library comes from NPASS (Zhao et al., *Nucleic Acids Res* 2023),
which should be cited in any work that uses it. If you merge further sources
with this pipeline, cite those too: NPCARE (Choi et al., *J Cheminform* 2017),
NPACT (Mangal et al., *Nucleic Acids Res* 2013), CancerHSP (Tao et al., *Sci
Rep* 2015) and COCONUT (Chandrasekhar et al., *Nucleic Acids Res* 2025).

---

## Quick start

```bash
# 1. environment
conda env create -f environment.yml
conda activate deck

# 2. prove it works before touching real data
python scripts/00_preflight.py

# 3. download the sources
bash scripts/00_fetch_sources.sh

# 4. build everything
bash run_all.sh
```

`run_all.sh` runs all five stages, checks each one, and stops at the first
real failure rather than carrying a broken file forward. Output lands in
`data/anticancer_final.sdf` and `data/pdbqt/`.

---

## Step 1. Environment

```bash
conda env create -f environment.yml
conda activate deck
```

If you already have an env and only want the missing pieces:

```bash
pip install meeko gemmi scipy pandas openpyxl
```

**Why that list.** Meeko 0.8.0 declares no dependencies at all, so
`pip install meeko` installs nothing else, and each missing package only
surfaces as an ImportError when you finally run it. Meeko needs numpy, pandas,
scipy, rdkit and gemmi simply to import. `openpyxl` is needed only if a source
database ships as .xlsx, which some do.

## Step 2. Preflight

```bash
python scripts/00_preflight.py
```

Checks every package, then runs the whole pipeline on built-in test molecules
and asserts the chemistry is right: a malformed SMILES is rejected, glucose is
removed as sugar-only, rutin deglycosylates and merges onto quercetin,
quercetin is PAINS-flagged but retained, betulinic acid survives despite its
logP, the sodium is stripped from naproxen sodium, provenance is preserved,
and the PDBQT output has genuine 3D coordinates and a torsion tree.

It exits 0 only if everything passes, and prints the exact fix command for
anything that does not. Run it before the real build and you will not meet an
environment problem halfway through.

## Step 3. Download the sources

```bash
bash scripts/00_fetch_sources.sh
```

That fetches everything that can be fetched non-interactively, verifies each
file is data rather than an HTML error page, and prints what remains to be
collected by hand. URLs were verified on 2026-11-12.

### What it downloads

**NPASS 3.0**, the backbone, from `https://bidd.group/NPASS/downloadFiles/`.
204,023 compounds and 1,048,756 activity records, split across several files:

| File | Contents |
|---|---|
| `NPASS3.0_naturalproducts_structure.txt` | SMILES, InChI, InChIKey |
| `NPASS3.0_activities.txt` | Targets, activity values, references |
| `NPASS3.0_target.txt` | Target names and types |
| `NPASS3.0_naturalproducts_generalinfo.txt` | Names and cross-references |

These land in `data/raw/npass_raw/`, a subfolder, because they must be joined
before use and because `--auto` deliberately ignores subfolders.

**COCONUT 2.0 complete CSV** is optional and off by default. Turn it on with
`WANT_COCONUT=1 bash scripts/00_fetch_sources.sh`. It comes from Zenodo at
`https://zenodo.org/records/13897048`, 207 MB. Note this is the *complete*
CSV, not the lite one, and only the complete one carries the collection tags
that identify the NPACT, NPCARE, AfroCancer and InPACdb entries inside
COCONUT.

### What has to be collected by hand

These three have no machine-fetchable download. I checked each one.

| Source | Why | What to do |
|---|---|---|
| **NPCARE** | The site blocks automated access | Use the CSV download tab in a browser |
| **CancerHSP** | No download link on the page | Export from the web interface |
| **NPACT** | No bulk download exists, only browse and search | Export what the interface allows, or skip |

**None of them is required.** NPASS alone gives a usable library, and the
COCONUT complete CSV covers NPACT and NPCARE content through its collections.
If you download any of them on your laptop, copy them across with:

```bash
scp ~/Downloads/npcare.csv you@yourserver:~/parp7/deck/data/raw/
```

Plain compound files go directly in `data/raw/`. Accepted formats are .csv,
.tsv, .txt, .sdf and .xlsx, optionally gzipped.

### The NPASS join

NPASS is three tables rather than a compound list, so the anticancer subset
has to be assembled. `run_all.sh` does this automatically when it finds the
files, but the step can also be run alone:

```bash
# confirm the column names in your release first
python scripts/11_npass_join.py --inspect data/raw/npass_raw/*.txt

python scripts/11_npass_join.py \
    --activities data/raw/npass_raw/npass_activities.txt \
    --structures data/raw/npass_raw/npass_structure.txt \
    --targets    data/raw/npass_raw/npass_target.txt \
    --max-um 10 \
    --out data/raw/npass_anticancer.csv
```

It selects activity records whose target is a cancer cell line, applies the
potency threshold, and joins the survivors to their structures. Units are
normalised, so a threshold in micromolar also matches records stored in
nanomolar or molar, and where a compound has several records the most potent
is kept.

Records in mass units such as µg/mL cannot be converted without a molecular
weight, so they are kept and flagged as `unconverted` in the output. Pass
`--drop-unconverted` for a cleaner criterion, at the cost of losing real
actives reported that way.

**The `--max-um` value is your operational definition of anticancer activity.
State it in the methods section.** 10 µM is a defensible default for a
cell-line screen.

## Step 4. Build

```bash
bash run_all.sh
```

Or stage by stage, which is what the driver does:

```bash
python scripts/10_build_from_sources.py --auto data/raw \
    --out data/anticancer_raw.sdf --report logs/source_overlap.tsv

python scripts/01_filter.py --in data/anticancer_raw.sdf \
    --out data/anticancer_final.sdf \
    --props logs/props.tsv --rejects logs/rejects.tsv

python scripts/04_prepare_pdbqt.py --in data/anticancer_final.sdf \
    --out data/pdbqt --jobs 0 --failures logs/pdbqt_failures.tsv
```

`--auto` finds every usable file in `data/raw`, detects which column holds the
structure, names each source after its filename, and reports what it skipped
and why. You never need to know a column name. If a file is skipped and you
think it should not be:

```bash
python scripts/10_build_from_sources.py --inspect data/raw/thatfile.csv
```

That prints its columns and which one is believed to hold the SMILES, InChI,
InChIKey and compound id. Override with `--smiles-col NAME` placed after the
`--source` it belongs to.

## Step 5. What you get

| Path | Contents |
|---|---|
| `data/anticancer_final.sdf` | The curated library |
| `data/pdbqt/` | One PDBQT per ligand, ready for Vina |
| `logs/provenance.txt` | Counts, versions, checksum |
| `logs/attrition_summary.txt` | How many compounds each filter removed |
| `logs/source_overlap.tsv` | Per-source new, duplicate and unusable counts |
| `logs/props.tsv` | Per-compound MW, cLogP, HBD, HBA, TPSA, PAINS flag |
| `logs/rejects.tsv` | Every removed compound with its reason |

Each compound in the SDF carries `deck_id`, `inchikey`, `sources`,
`n_sources`, `anticancer_annotated`, `pains_flag` and `sugars_removed`.

---

## Curation decisions worth knowing

**PAINS are flagged, not removed.** PAINS A/B/C matches a large fraction of
natural products, because catechols, quinones and many flavonoid cores are
PAINS substructures. Quercetin is a match. Deleting them silently removes much
of the plant chemistry the screen exists to interrogate. Matches stay with
`pains_flag=1`, so the decision can be made on the ranked hits and reported.
`--pains-mode remove` is available if you want it, and the count is printed
either way.

**Deglycosylation** cleaves terminal sugars iteratively and keeps the
aglycone, removing molecules that are nothing but sugar. Four reasons: docking
scoring functions handle large flexible solvated sugars badly; the aglycone is
often the species that circulates, since dietary glycosides are hydrolysed at
the intestinal brush border; the PARP7 nicotinamide pocket is too small for a
bulky glycoside; and it collapses glycoside families onto a shared aglycone,
removing redundancy. Glycosides whose sugar is essential to activity, such as
cardiac glycosides, are out of scope by construction. The originals are kept
in the rejection log.

**The molecular weight window is 150 to 600 Da and is applied after
deglycosylation**, so a 900 Da saponin is reduced into range rather than
discarded. The upper bound comes from the target: RBN-2397 occupies the PARP7
site at 523 Da.

**The rule of five is applied with one violation tolerated**, specifically so
that lipophilic terpenoids are retained. A zero-violation filter would
preferentially delete the scaffold classes that characterise plant secondary
metabolism.

**3D geometry is built once**, in `04_prepare_pdbqt.py`, with ETKDGv3 followed
by MMFF94. Nothing earlier needs coordinates.

---

## This library is not a benchmark set

A compound with a reported cancer cell-line IC50 inhibits cancer cell growth
by an unknown mechanism, and almost none of these act through PARP7. Using
them as the positive class in a retrospective enrichment calculation would
produce a number that measures nothing about the model.

Enrichment actives come from ChEMBL and BindingDB for UniProt Q7Z3E1,
screened against property-matched decoys, and are built separately. The
`anticancer_annotated` flag is for prioritising hits after docking, not for
validating the receptor.

---

## Optional: adding COCONUT for breadth (not used for the deposited library)

The library above covers chemistry already examined for anticancer activity.
To add unexplored natural product space, build a second library from COCONUT
and screen both.

```bash
bash scripts/00_download.sh

python scripts/07_restore_ids.py --list-fields data/raw/<coconut>.sdf   # find the id field

python scripts/01_filter.py --in data/raw/<coconut>.sdf \
    --id-prop <the id field> \
    --out data/coconut_filtered.sdf \
    --props logs/coconut_props.tsv --rejects logs/coconut_rejects.tsv

python scripts/03_select.py --in data/coconut_filtered.sdf \
    --out data/coconut_35k.sdf \
    --anticancer-keys data/anticancer_inchikeys.txt --require-annotations \
    --target 35000 --method leader --radius 0.50 \
    --manifest logs/coconut_manifest.tsv
```

Three things to know.

**Pass `--id-prop`.** The COCONUT SDF uses none of the field names
`01_filter.py` checks by default, so without it every compound gets a
placeholder `REC########`. If a run has already produced placeholders,
`07_restore_ids.py` recovers the real identifiers from them without repeating
the filtering pass, because the number is the source record index.

**Watch the radius.** The line `leader picker at distance D returned N
centroids` is the one that matters. If N is below the target, the shortfall is
filled at random and it is no longer a diversity selection. Lower `--radius`
until N exceeds the target, then the script truncates to the first N
centroids, which are the well-separated ones. Overshooting is correct.

**Revoked records.** The November 2026 SDF holds 737,309 records and the CSV
lists 738,793 identifiers, both close to the 738,795 COCONUT reports in total,
and neither distribution carries a revoked flag. The likely reading is that
the published total already excludes revoked entries. Confirm at
https://coconut.naturalproducts.net/stats before claiming either way.

---

## Scripts

| Script | Purpose |
|---|---|
| `00_preflight.py` | Environment check and end-to-end self test |
| `00_fetch_sources.sh` | Downloads NPASS, and optionally COCONUT, with checks |
| `11_npass_join.py` | Joins the NPASS tables and selects cancer cell-line actives |
| `10_build_from_sources.py` | Several heterogeneous source files to one SDF |
| `01_filter.py` | Chemical curation |
| `04_prepare_pdbqt.py` | 3D generation and PDBQT with torsion trees |
| `run_all.sh` | Runs the whole pipeline with checks |
| `00_download.sh` | Optional. Fetches COCONUT |
| `02_merge.py` | Optional. Merges independently filtered sets |
| `03_select.py` | Optional. Diversity selection for a large library |
| `05_make_anticancer_keys.py` | Builds an InChIKey list from any source |
| `06_drop_revoked.py` | Subsets an SDF against an identifier whitelist |
| `07_restore_ids.py` | Recovers real identifiers from REC placeholders |
| `00b_imppat_fetch.sh` | Optional. Polite per-compound IMPPAT fetch |
