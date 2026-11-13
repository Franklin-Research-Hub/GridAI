#!/usr/bin/env bash
# Simplify the national pipeline's WGS84 GeoJSON into per-state tract TopoJSON plus national counties and states.
# Simplification changes drawn outlines only; every number comes from the CSV/JSON outputs.
#   bash scripts/make_national_topo.sh <NATIONAL_OUT>
set -euo pipefail
OUT="${1:?national output folder}"
MS="${MAPSHAPER:-npx --yes mapshaper@0.7.101}"
mkdir -p "$OUT/topo"
for f in "$OUT"/geo/tracts_*.geojson; do
  s=$(basename "$f" .geojson); s=${s#tracts_}
  $MS "$f" -simplify 3% keep-shapes -o format=topojson quantization=1e5 "$OUT/topo/t$s.json"
done
$MS "$OUT/geo/us_counties.geojson" -simplify 4% keep-shapes -o format=topojson quantization=1e5 "$OUT/topo/us_counties.json"
$MS "$OUT/geo/us_states.geojson"   -simplify 4% keep-shapes -o format=topojson quantization=1e5 "$OUT/topo/us_states.json"
