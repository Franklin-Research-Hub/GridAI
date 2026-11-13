#!/usr/bin/env python3
"""
Build the GridAI (SEGDA) site from one pipeline output folder.

    python scripts/build_site.py <PIPELINE_OUT> <SITE_DIR> [NATIONAL_BUILD]

NATIONAL_BUILD (optional) is the folder written by pipeline/segda_national.py plus
scripts/make_national_topo.sh. Its gate must also be PASS.

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


def national(nb, dst):
    """Copy the national layer into site/data/national and return the small summary embedded in the page."""
    nb = Path(nb)
    R = json.loads((nb / "national_results.json").read_text())
    if R["run"]["gate"] != "PASS":
        sys.exit(f"National gate is {R['run']['gate']}; the site is not rebuilt.")
    out = dst / "data" / "national"
    out.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(nb / "national_tracts.csv", dtype={"GEOID": str, "st": str, "fips": str})
    # per state: GEOID -> [ejbi, svi, eji, ebm, pop, county fips, masked, svi_hidden]
    for st, g in t.groupby("st"):
        d = {r.GEOID: [num(r.ejbi), num(r.svi), num(r.eji), num(r.ebm), None if pd.isna(r.pop) else int(r.pop),
                       r.fips, int(bool(r.masked)), int(bool(r.svi_hidden))] for r in g.itertuples()}
        (out / f"v{st}.json").write_text(json.dumps(d, separators=(",", ":"), allow_nan=False))
    for f in (nb / "topo").glob("*.json"):
        shutil.copyfile(f, out / f.name)
    for f in ["national_counties.csv", "national_states.csv", "national_tracts.csv", "national_results.json",
              "national_audit.csv", "national_sources.csv"]:
        shutil.copyfile(nb / f, out / f)
    c = pd.read_csv(nb / "national_counties.csv", dtype={"fips5": str, "STATEFP": str})
    keep = ["fips5", "STATEFP", "NAME", "pop", "ejbi", "svi_c", "obi", "cs", "tracts", "tracts_valid", "tracts_high",
            "pop_high", "tract_min", "tract_max", "masked_tracts", "masked_pop", "svi_tracts_high", "svi_hidden_tracts",
            "svi_hidden_pop", "lisa_ejbi", "lisa_obi"]
    rows = [[(num(v) if isinstance(v, float) else v) for v in r] for r in c[keep].itertuples(index=False)]
    s = pd.read_csv(nb / "national_states.csv", dtype={"statefp": str})
    states = [{k: (num(v) if isinstance(v, float) else v) for k, v in r.items()} for r in s.to_dict("records")]
    return {"run": R["run"], "nat": R["national"], "rel": R["relations"], "spatial": R["spatial"], "conc": R["concordance"],
            "audit": R["audit"], "sources": R["sources"], "hazards": R["hazards"], "thr": R["thresholds"],
            "ccols": keep, "counties": rows, "states": states}


def main(src, dst, nb=None):
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
    dst.mkdir(parents=True, exist_ok=True)
    N = national(nb, dst) if nb else None
    data = json.dumps(D, separators=(",", ":"), allow_nan=False, default=str)
    ndata = json.dumps(N, separators=(",", ":"), allow_nan=False, default=str)
    html = (ROOT / "site_src" / "template.html").read_text(encoding="utf-8")
    import base64
    imgd = ROOT / "site_src" / "img"
    uri = lambda f: "data:image/png;base64," + base64.b64encode((imgd / f).read_bytes()).decode()
    html = html.replace("__LOGO__", uri("fhsc-logo.png")).replace("__ICON__", uri("fhsc-icon.png"))
    html = html.replace("__HEADSHOT__", "data:image/jpeg;base64," + base64.b64encode((imgd / "vlf-headshot.jpg").read_bytes()).decode())
    html = (html.replace("__NATIONAL__", ndata.replace("</", "<\\/"))
                .replace("__DATA__", data.replace("</", "<\\/"))
                .replace("__TRACTS__", (src / "tracts.topo.json").read_text())
                .replace("__COUNTIES__", (src / "counties.topo.json").read_text()))
    dst.mkdir(parents=True, exist_ok=True)
    (dst / "index.html").write_text(html, encoding="utf-8")
    (dst / ".nojekyll").write_text("")
    (dst / "CNAME").write_text("gridaitm.com\n")   # custom domain for GitHub Pages
    (dst / "robots.txt").write_text("User-agent: *\nAllow: /\n\nSitemap: https://gridaitm.com/sitemap.xml\n")
    import datetime as _dt
    (dst / "sitemap.xml").write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n'
        f'  <url><loc>https://gridaitm.com/</loc><lastmod>{_dt.date.today().isoformat()}</lastmod><changefreq>weekly</changefreq><priority>1.0</priority></url>\n'
        '</urlset>\n')
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
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    main(*sys.argv[1:])
