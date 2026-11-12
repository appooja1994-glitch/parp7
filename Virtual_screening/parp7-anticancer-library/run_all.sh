#!/usr/bin/env bash
# Build the anticancer screening library end to end.
#
#   bash run_all.sh
#
# Put the downloaded source databases in data/raw/ first. Anything that cannot
# be read is reported and skipped rather than stopping the run. Every stage is
# checked, and the script stops at the first real failure instead of carrying
# a broken file forward.
set -uo pipefail

RAW="${RAW:-data/raw}"
OUT="${OUT:-data}"
LOGS="${LOGS:-logs}"
JOBS="${JOBS:-0}"
PY="${PY:-python}"

mkdir -p "$OUT" "$LOGS"

stage() { printf '\n\033[1m==> %s\033[0m\n' "$1"; }
die()   { printf '\n\033[1mSTOPPED: %s\033[0m\n' "$1" >&2; exit 1; }
count() { grep -c '^\$\$\$\$' "$1" 2>/dev/null || echo 0; }

# ---------------------------------------------------------------- 0 preflight
stage "Stage 0 of 5  Preflight"
$PY scripts/00_preflight.py || die "preflight failed. Fix the items above and rerun."

# ---------------------------------------------------------------- 1 sources
stage "Stage 1 of 5  Checking source files"
[ -d "$RAW" ] || die "no $RAW directory. Create it and put the downloaded databases there."

NPASS_DIR="$RAW/npass_raw"
NPASS_OUT="$RAW/npass_anticancer.csv"
HAVE_NPASS=0
if [ -s "$NPASS_DIR/npass_activities.txt" ] && [ -s "$NPASS_DIR/npass_structure.txt" ]; then
  HAVE_NPASS=1
fi

NFILES=$(find "$RAW" -maxdepth 1 -type f ! -name '.*' | wc -l | tr -d ' ')
echo "  $NFILES file(s) directly in $RAW"
[ "$HAVE_NPASS" = "1" ] && echo "  NPASS release present in $NPASS_DIR"

# Either kind of input is enough on its own. The NPASS files live in a
# subfolder because they have to be joined first, so counting only the top
# level would wrongly call the directory empty.
if [ "$NFILES" -eq 0 ] && [ "$HAVE_NPASS" = "0" ]; then
  die "no usable input found.
         Nothing directly in $RAW, and no NPASS release in $NPASS_DIR.
         Run:  bash scripts/00_fetch_sources.sh"
fi

# ---------------------------------------------------------------- 1b npass
if [ "$HAVE_NPASS" = "1" ]; then
  if [ -s "$NPASS_OUT" ]; then
    echo "  NPASS subset already built, reusing $NPASS_OUT"
  else
    stage "Stage 1b  Joining NPASS and selecting cancer cell-line actives"
    $PY scripts/11_npass_join.py \
        --activities "$NPASS_DIR/npass_activities.txt" \
        --structures "$NPASS_DIR/npass_structure.txt" \
        --targets    "$NPASS_DIR/npass_target.txt" \
        --max-um "${MAX_UM:-10}" \
        --out "$NPASS_OUT" || die "the NPASS join failed"
    [ -s "$NPASS_OUT" ] || die "the NPASS join produced no output"
  fi
fi

# ---------------------------------------------------------------- 2 build
stage "Stage 2 of 5  Building one library from all usable sources"
$PY scripts/10_build_from_sources.py \
    --auto "$RAW" \
    --out "$OUT/anticancer_raw.sdf" \
    --report "$LOGS/source_overlap.tsv" || die "build failed"
RAWN=$(count "$OUT/anticancer_raw.sdf")
[ "$RAWN" -gt 0 ] || die "the build produced an empty file"
echo "  built: $RAWN compounds"

# ---------------------------------------------------------------- 3 curate
stage "Stage 3 of 5  Chemical curation"
$PY scripts/01_filter.py \
    --in  "$OUT/anticancer_raw.sdf" \
    --out "$OUT/anticancer_final.sdf" \
    --props   "$LOGS/props.tsv" \
    --rejects "$LOGS/rejects.tsv" || die "curation failed"
FINN=$(count "$OUT/anticancer_final.sdf")
[ "$FINN" -gt 0 ] || die "curation removed everything. Check $LOGS/rejects.tsv"
echo "  curated: $FINN compounds"

# ---------------------------------------------------------------- 4 pdbqt
stage "Stage 4 of 5  3D coordinates and PDBQT"
$PY scripts/04_prepare_pdbqt.py \
    --in  "$OUT/anticancer_final.sdf" \
    --out "$OUT/pdbqt" \
    --jobs "$JOBS" \
    --failures "$LOGS/pdbqt_failures.tsv" || die "ligand preparation failed"
NPQ=$(ls "$OUT/pdbqt"/*.pdbqt 2>/dev/null | wc -l | tr -d ' ')
[ "$NPQ" -gt 0 ] || die "no PDBQT files were written"

# ---------------------------------------------------------------- 4b smi
stage "Stage 4b  SMILES export for the repository"
$PY scripts/12_export_smi.py --in "$OUT/anticancer_raw.sdf" \
    --out "$OUT/total_molecules.smi" || die "SMILES export of the raw set failed"
$PY scripts/12_export_smi.py --in "$OUT/anticancer_final.sdf" \
    --out "$OUT/final_molecules.smi" \
    --props "$OUT/final_molecules_annotations.tsv" \
    || die "SMILES export of the final set failed"

# ---------------------------------------------------------------- 5 record
stage "Stage 5 of 5  Provenance record"
{
  echo "built                 : $(date -u '+%Y-%m-%d %H:%M UTC')"
  echo "compounds built       : $RAWN"
  echo "compounds after curation: $FINN"
  echo "ligands prepared      : $NPQ"
  $PY -c "import rdkit; print('RDKit                 :', rdkit.__version__)"
  $PY -c "import meeko; print('Meeko                 :', meeko.__version__)"
  echo "sha256                : $(sha256sum "$OUT/anticancer_final.sdf" | cut -d' ' -f1)"
  echo
  echo "source files used:"
  sed 's/^/  /' "$LOGS/source_overlap.tsv" 2>/dev/null
} > "$LOGS/provenance.txt"

tail -n +2 "$LOGS/rejects.tsv" 2>/dev/null | cut -f2 | sort | uniq -c | sort -rn \
    > "$LOGS/attrition_summary.txt"

cat "$LOGS/provenance.txt"
echo
echo "attrition by reason:"
sed 's/^/  /' "$LOGS/attrition_summary.txt"

printf '\n\033[1m==> Done\033[0m\n'
cat <<SUM

  library   : $OUT/anticancer_final.sdf   ($FINN compounds)
  ligands   : $OUT/pdbqt/                 ($NPQ files, ready for Vina)
  total.smi : $OUT/total_molecules.smi    ($RAWN compounds, before curation)
  final.smi : $OUT/final_molecules.smi    ($FINN compounds, for the repository)
  annots    : $OUT/final_molecules_annotations.tsv
  provenance: $LOGS/provenance.txt
  attrition : $LOGS/attrition_summary.txt
  per source: $LOGS/source_overlap.tsv

SUM
