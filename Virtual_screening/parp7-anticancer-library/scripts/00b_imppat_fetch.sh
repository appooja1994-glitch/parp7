#!/usr/bin/env bash
# Polite, resumable fetch of IMPPAT 2.0 3D SDFs, one compound per request.
# Read https://cb.imsc.res.in/imppat/ terms before running. 1 req/s, ~6 h.
set -euo pipefail
OUT="${OUT:-data/raw/imppat/sdf}"
FIRST="${FIRST:-1}"; LAST="${LAST:-18000}"; DELAY="${DELAY:-1}"
mkdir -p "$OUT"
for i in $(seq "$FIRST" "$LAST"); do
  id=$(printf 'IMPHY%06d' "$i")
  f="$OUT/$id.sdf"
  [ -s "$f" ] && continue
  curl -sS -L --max-time 30 --retry 2 \
       -o "$f" "https://cb.imsc.res.in/imppat/images/3D/SDF/${id}_3D.sdf" || true
  # drop HTML error pages and empties
  if [ ! -s "$f" ] || head -c 200 "$f" | grep -qi '<html'; then rm -f "$f"; fi
  sleep "$DELAY"
  [ $((i % 500)) -eq 0 ] && echo "... $id  ($(ls "$OUT" | wc -l) fetched)"
done
echo "Concatenating"
cat "$OUT"/*.sdf > data/raw/imppat/imppat_3d.sdf
echo "Wrote data/raw/imppat/imppat_3d.sdf ($(grep -c '\$\$\$\$' data/raw/imppat/imppat_3d.sdf) records)"
