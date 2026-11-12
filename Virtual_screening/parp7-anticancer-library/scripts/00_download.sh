#!/usr/bin/env bash
# Fetch the source datasets for the PARP7 natural-product screening deck.
# Run from the package root:  bash scripts/00_download.sh
set -euo pipefail

RAW="${RAW_DIR:-data/raw}"
mkdir -p "$RAW"
cd "$RAW"

say() { printf '\n== %s ==\n' "$1"; }

# ---------------------------------------------------------------- COCONUT 2.0
# Release used here: November 2026. Licence: CC0 (no attribution required).
# Check https://coconut.naturalproducts.net/download for a newer month and
# edit COCO_MONTH if you want it; keep the month you actually used in your
# methods section, because COCONUT is re-released monthly.
COCO_MONTH="${COCO_MONTH:-2026-11}"
COCO_TAG="${COCO_TAG:-11-2026}"
COCO_BASE="https://coconut.s3.uni-jena.de/prod/downloads/${COCO_MONTH}"

say "COCONUT 3D SDF (~305 MB zipped)"
curl -L -C - --retry 5 --retry-delay 10 --fail-with-body \
     -o "coconut_sdf_3d-${COCO_TAG}.zip" \
     "${COCO_BASE}/coconut_sdf_3d-${COCO_TAG}.zip"

# The *lite* CSV does NOT carry collection tags. The anticancer annotation
# needs a distribution that does, so try the fuller ones first and report
# which one actually downloaded.
say "COCONUT CSV with collection tags"
GOT_CSV=""
for CAND in "coconut_complete-${COCO_TAG}.csv.zip" \
            "coconut_csv-${COCO_TAG}.zip" \
            "coconut_complete-${COCO_TAG}.zip" \
            "coconut_csv_lite-${COCO_TAG}.zip"; do
  echo "  trying ${CAND}"
  if curl -L -C - --retry 3 --retry-delay 5 --fail --silent --show-error \
          -o "${CAND}" "${COCO_BASE}/${CAND}"; then
    echo "  got ${CAND}"
    GOT_CSV="${CAND}"
    [ "${CAND}" = "coconut_csv_lite-${COCO_TAG}.zip" ] && cat <<'W'

  ! Only the LITE csv was available. It has no collection column, so the
  ! anticancer annotation cannot be built from it. After unpacking, run
  !     python scripts/05_make_anticancer_keys.py --inspect <the csv>
  ! to confirm, then build the key list from the 3D SDF or from NPCARE.
W
    break
  fi
  rm -f "${CAND}"
done
[ -n "${GOT_CSV}" ] || echo "  ! no CSV distribution downloaded; check the names at
    https://coconut.naturalproducts.net/download"

say "Unpacking"
for z in coconut_sdf_3d-${COCO_TAG}.zip ${GOT_CSV}; do
  [ -f "$z" ] || continue
  unzip -o -q "$z" -d "${z%.zip}"
done

say "Record counts in the downloaded SDF"
for f in coconut_sdf_3d-${COCO_TAG}/*.sdf; do
  [ -f "$f" ] || continue
  echo "  $f : $(grep -c '^\$\$\$\$' "$f") records"
done
cat <<'N'
  COCONUT published 738,795 molecules in total, of which 38,869 are revoked,
  leaving 699,926 active. If the count above is near 738,795 the distribution
  still contains revoked entries and they must be removed before filtering.
  If it is near 699,926 they are already excluded. Record which, because the
  methods section states that revoked records were removed.
N

say "Checksums (record these in your thesis methods)"
sha256sum coconut_*.zip | tee ../../logs/source_checksums.txt

cat <<'NOTE'

------------------------------------------------------------------------------
IMPPAT 2.0 is NOT bulk-downloadable over plain HTTP; the web interface serves
one compound at a time. Two options:

  (a) Preferred - request the full phytochemical library from the authors via
      https://cb.imsc.res.in/imppat/  (contact page). They supply the
      17,967-compound set on request for academic use.

  (b) Scrape politely from the per-compound endpoints. IDs run IMPHY000001
      upward. Use scripts/00b_imppat_fetch.sh, which rate-limits to 1 req/s
      and resumes. Expect ~6 hours. Check their terms before running it.

Either way, place the result at:  data/raw/imppat/imppat_3d.sdf
If IMPPAT is absent the pipeline still runs; it just skips the merge step.
------------------------------------------------------------------------------
NOTE
