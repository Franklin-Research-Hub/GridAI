#!/usr/bin/env python3
"""
Build the SEGDA SMART Grid site from one pipeline output folder.

    python scripts/build_site.py <PIPELINE_OUT> <SITE_DIR>

Refuses to build if the pipeline's validation gate did not pass, so the site
can never show numbers that no longer reproduce the paper.
Author: Victoria Love Franklin. Copyright (c) 2026. All rights reserved.
"""
import json
import math
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def num(x, n=4):
    return None if x is None or (isinstance(x, float) and math.isnan(x)) else round(float(x), n)


def main(src, dst):
    src, dst = Path(src), Path(dst)
    R = json.loads((src / "segda_results.json").read_text())
    if R["run"]["gate"] != "PASS":
        sys.exit(f"Validation gate is {R['run']['gate']}: {R['run']['n_mismatch']} locked values did not "
                 "reproduce. The site is not rebuilt; the previous version stays published.")
    c = pd.read_csv(src / "segda_counties.csv")
    t = pd.read_csv(src / "segda_tracts.csv", dtype={"GEOID": str})
    comp = pd.read_csv(src / "segda_climate_components.csv")
    counties = []
    for r in c.itertuples():
        counties.append({k: (num(getattr(r, k)) if isinstance(getattr(r, k), float) else getattr(r, k))
                         for k in c.columns})
    comp_map = {int(r["fips"]): {k: num(v, 3) for k, v in r.items() if k != "fips"} for _, r in comp.iterrows()}
    # tracts packed as GEOID -> [ejbi, svi, ces, pb, pop, fips, masked]
    tr = {}
    for r in t.itertuples():
        tr[r.GEOID] = [num(r.ejbi), num(r.svi), num(r.ces), num(r.pb),
                       None if pd.isna(r.Population) else int(r.Population), int(r.fips), int(bool(r.masked))]
    D = {"run": R["run"], "stats": R["stats"], "corr": R["correlations"], "sens": R["sensitivity"],
         "validation": R["validation"], "audit": R["audit"], "sources": R["sources"], "hazards": R["hazards"], "thr": R["thresholds"],
         "counties": counties, "climate": comp_map, "tracts": tr}
    data = json.dumps(D, separators=(",", ":"), allow_nan=False, default=str)
    html = (ROOT / "site_src" / "template.html").read_text(encoding="utf-8")
    html = (html.replace("__DATA__", data.replace("</", "<\\/"))
                .replace("__TRACTS__", (src / "tracts.topo.json").read_text())
                .replace("__COUNTIES__", (src / "counties.topo.json").read_text()))
    dst.mkdir(parents=True, exist_ok=True)
    (dst / "index.html").write_text(html, encoding="utf-8")
    (dst / ".nojekyll").write_text("")
    dl = dst / "data"
    dl.mkdir(exist_ok=True)
    for f in ["segda_counties.csv", "segda_tracts.csv", "segda_validation.csv",
              "segda_climate_components.csv", "segda_results.json", "segda_audit.csv", "segda_sources.csv"]:
        shutil.copyfile(src / f, dl / f)
    img = ROOT / "site_src" / "img"
    if img.exists():
        shutil.copytree(img, dst / "img", dirs_exist_ok=True)
    print(f"Built {dst/'index.html'} ({(dst/'index.html').stat().st_size/1e6:.2f} MB), gate {R['run']['gate']}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
