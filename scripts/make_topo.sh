#!/usr/bin/env bash
# Simplify the pipeline's EPSG:3310 GeoJSON into TopoJSON for the web map.
# Simplification changes drawn outlines only; every number shown comes from the CSV/JSON outputs.
#   bash scripts/make_topo.sh <PIPELINE_OUT>
set -euo pipefail
OUT="${1:?pipeline output folder}"
MS="${MAPSHAPER:-npx --yes mapshaper@0.7.101}"
$MS "$OUT/geo/tracts.geojson"   -simplify 4% keep-shapes -o format=topojson quantization=1e5 "$OUT/tracts.topo.json"
$MS "$OUT/geo/counties.geojson" -simplify 8% keep-shapes -o format=topojson quantization=1e5 "$OUT/counties.topo.json"
