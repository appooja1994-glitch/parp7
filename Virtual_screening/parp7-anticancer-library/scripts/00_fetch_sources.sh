#!/usr/bin/env bash
# Download the source databases that can be fetched non-interactively.
#
#   bash scripts/00_fetch_sources.sh
#
# URLs below were verified on 2026-11-12. If one 404s the site has moved it;
# the script says so and carries on with the rest.
set -uo pipefail

RAW="${RAW:-data/raw}"
WANT_COCONUT="${WANT_COCONUT:-0}"     # set to 1 to also fetch COCONUT (207 MB)
mkdir -p "$RAW/npass_raw"

ok=0; bad=0

get() {  # get <url> <destination>
  local url="$1" dest="$2" name
  name="$(basename "$dest")"
  if [ -s "$dest" ]; then
    echo "  have    $name"
    ok=$((ok+1)); return 0
  fi
  printf '  fetch   %s ... ' "$name"
  if curl -L -C - --retry 4 --retry-delay 5 --connect-timeout 30 \
          --fail --silent --show-error -o "$dest" "$url"; then
    # a site that answers an error with an HTML page still returns 200
    if head -c 512 "$dest" | grep -qiE '<!doctype html|<html'; then
      echo "got an HTML page, not data. Skipping."
      rm -f "$dest"; bad=$((bad+1)); return 1
    fi
    echo "$(du -h "$dest" | cut -f1)"
    ok=$((ok+1)); return 0
  fi
  echo "failed"
  rm -f "$dest"; bad=$((bad+1)); return 1
}

echo
echo "NPASS 3.0  (the backbone: 204,023 compounds, 1,048,756 activity records)"
NP="https://bidd.group/NPASS/downloadFiles"
get "$NP/NPASS3.0_naturalproducts_structure.txt" "$RAW/npass_raw/npass_structure.txt"
get "$NP/NPASS3.0_activities.txt"                "$RAW/npass_raw/npass_activities.txt"
get "$NP/NPASS3.0_target.txt"                    "$RAW/npass_raw/npass_target.txt"
get "$NP/NPASS3.0_naturalproducts_generalinfo.txt" "$RAW/npass_raw/npass_generalinfo.txt"

if [ "$WANT_COCONUT" = "1" ]; then
  echo
  echo "COCONUT 2.0 complete CSV from Zenodo (207 MB, carries collection tags)"
  Z="https://zenodo.org/records/13897048/files"
  if get "$Z/coconut_complete-10-2024.csv.zip?download=1" \
         "$RAW/coconut_complete.csv.zip"; then
    echo "  unzip   coconut_complete.csv.zip"
    unzip -o -q "$RAW/coconut_complete.csv.zip" -d "$RAW/coconut_complete" \
      && echo "  ->      $RAW/coconut_complete/"
  fi
fi

echo
echo "------------------------------------------------------------------"
echo "  downloaded or already present: $ok    failed: $bad"
echo "------------------------------------------------------------------"
cat <<'MANUAL'

These three have no machine-fetchable download and must be collected by hand,
if you want them at all. NPASS on its own is enough for a usable library.

  NPCARE      http://silver.sejong.ac.kr/npcare
              Use the CSV download tab in a browser. The site blocks automated
              access, so curl and wget will not work on it.

  CancerHSP   https://old.tcmsp-e.com/CancerHSP.php
              Browse and export from the web interface.

  NPACT       https://webs.iiitd.edu.in/raghava/npact/
              No bulk download exists, only browse and search. Export what the
              interface allows, or skip it. Its 1,574 compounds are largely
              covered by NPASS and by COCONUT's collections.

Downloaded on your laptop instead? Copy them over with:

  scp ~/Downloads/npcare.csv  you@thisserver:~/parp7/deck/data/raw/

Put plain compound files directly in data/raw/. The NPASS files stay in
data/raw/npass_raw/ because they need joining first, which 11_npass_join.py
does, and because --auto deliberately ignores subfolders.
MANUAL
