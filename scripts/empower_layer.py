#!/usr/bin/env python3
"""
GridAI (SEGDA), national layer: electricity-dependent Medicare beneficiaries where outage burden is high.

    python scripts/empower_layer.py <NATIONAL_BUILD_DIR>

Joins the HHS emPOWER county counts (national_build/empower/) to the audited national county table and
reports how many electricity-dependent residents live in high-outage-burden counties. This describes who
is most exposed when the power goes out. It does not measure health outcomes.

Source: HHS ASPR emPOWER public REST Service, layer "Electricity Dependent DME - ALL - CountyLevel",
retrieved 13 Nov 2026 (data last edited 9 Nov 2026). Counts of 1-10 are masked by HHS and shown as 11.

Author: Victoria Love Franklin. Copyright (c) 2026. All rights reserved.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# FIPS renamed with identical boundaries: emPOWER still uses the old code.
RENAMED = {"46113": "46102",   # Shannon County -> Oglala Lakota County, SD
           "02270": "02158"}   # Wade Hampton -> Kusilvak Census Area, AK
# Not joinable without inventing a split: CT old counties vs 2022 planning regions; AK Valdez-Cordova (02261)
# split into Chugach (02063) and Copper River (02066). Left unmatched, as GridAI does for Connecticut outages.
NOT_JOINED_REASON = {"09": "Connecticut: emPOWER uses the 8 former counties, GridAI the 9 planning regions",
                     "02063": "Valdez-Cordova split in 2019; emPOWER reports it combined",
                     "02066": "Valdez-Cordova split in 2019; emPOWER reports it combined"}
STATES_50_DC = set("01 02 04 05 06 08 09 10 11 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56".split())


def compute(nb):
    nb = Path(nb)
    ed = nb / "empower"
    e = pd.read_csv(ed / "empower_county_2026-11.csv", dtype={"fips_recode": str, "fips": str})
    s = pd.read_csv(ed / "empower_state_2026-11.csv", dtype={"state_fips": str})
    c = pd.read_csv(nb / "national_counties.csv", dtype={"fips5": str, "STATEFP": str})
    R = json.loads((nb / "national_results.json").read_text())

    e["st"] = e["fips_recode"].str[:2]
    e["fips5"] = e["fips_recode"].replace(RENAMED)
    m = c.merge(e[["fips5", "medicare_benes", "power_dependent_dme"]], on="fips5", how="left")
    m["dme_rate"] = 1000 * m["power_dependent_dme"] / m["medicare_benes"]

    # High outage burden = top quarter of counties by Primary OBI (same top-quartile convention as EJBI).
    obi_q75 = float(m["obi"].quantile(0.75))
    m["obi_high"] = m["obi"] >= obi_q75
    m["ejbi_high"] = m["ejbi"] >= R["thresholds"]["county"]
    m["has_hidden"] = m["masked_tracts"] > 0
    j = m[m["power_dependent_dme"].notna()]
    jo = j[j["obi"].notna()]

    dme_us = float(e.loc[e["st"].isin(STATES_50_DC), "power_dependent_dme"].sum())
    dme_joined = float(j["power_dependent_dme"].sum())
    dme_obi = float(jo["power_dependent_dme"].sum())
    hi = jo[jo["obi_high"]]
    both_ej = hi[hi["ejbi_high"]]
    both_hidden = hi[hi["has_hidden"]]
    r_obi = stats.pearsonr(jo["dme_rate"], jo["obi"])
    rs_obi = stats.spearmanr(jo["dme_rate"], jo["obi"])
    je = jo[jo["ejbi"].notna()]
    r_ej = stats.pearsonr(je["dme_rate"], je["ejbi"])
    rate_hi = 1000 * hi["power_dependent_dme"].sum() / hi["medicare_benes"].sum()
    lo = jo[~jo["obi_high"]]
    rate_lo = 1000 * lo["power_dependent_dme"].sum() / lo["medicare_benes"].sum()

    top = hi.sort_values("power_dependent_dme", ascending=False).head(10)
    abbr = dict(zip(e["st"], e["state"]))
    top_rows = [{"county": f"{r.NAMELSAD}, {abbr.get(r.STATEFP, r.STATEFP)}", "fips5": r.fips5,
                 "power_dependent_dme": int(r.power_dependent_dme), "medicare_benes": int(r.medicare_benes),
                 "obi": round(float(r.obi), 3), "ejbi": None if pd.isna(r.ejbi) else round(float(r.ejbi), 3),
                 "hidden_tracts": int(r.masked_tracts)} for r in top.itertuples()]

    # ---------------- checks ----------------
    # HHS de-identifies small cells in the county layer, which can only raise county counts; so county sums
    # may exceed the state totals (most in rural states with many small counties) but should never fall below them.
    e50 = e[e["st"].isin(STATES_50_DC)]
    gs = e50.groupby("st")[["power_dependent_dme", "medicare_benes"]].sum()
    sv = s.set_index("state_fips")
    dd = gs["power_dependent_dme"] - sv.loc[gs.index, "power_dependent_dme"]
    bd = (gs["medicare_benes"] - sv.loc[gs.index, "medicare_benes"]).abs() / sv.loc[gs.index, "medicare_benes"]
    state_dme = float(sv.loc[gs.index, "power_dependent_dme"].sum())
    excess_pct = 100 * (float(gs["power_dependent_dme"].sum()) - state_dme) / state_dme
    worst = (100 * dd / sv.loc[gs.index, "power_dependent_dme"]).sort_values(ascending=False)
    unmatched = m[m["power_dependent_dme"].isna()]
    unexplained = [f for f in unmatched["fips5"] if not (f[:2] in NOT_JOINED_REASON or f in NOT_JOINED_REASON)]
    v = e.dropna(subset=["power_dependent_dme", "medicare_benes"])
    checks = [
        ("Every emPOWER county record has a unique FIPS code", e["fips_recode"].is_unique, f"{len(e):,} records"),
        ("Electricity-dependent counts never exceed Medicare beneficiaries", bool((v["power_dependent_dme"] <= v["medicare_benes"]).all()),
         f"{len(v):,} county records with values; {len(e) - len(v)} territory records without values"),
        ("County Medicare enrollment matches HHS state totals within 0.2% in every state", bool((bd <= 0.002).all()),
         f"largest gap {100 * bd.max():.3f}%"),
        ("County electricity-dependent counts are never below HHS state totals (de-identification only raises them)", bool((dd >= 0).all()),
         f"{int((dd >= 0).sum())} of {len(dd)} states at or above"),
        ("County electricity-dependent counts exceed HHS state totals by under 2% nationally", 0 <= excess_pct < 2,
         f"{excess_pct:.2f}% nationally; largest in {', '.join(abbr_s for abbr_s in e.drop_duplicates('st').set_index('st').loc[worst.index[:3], 'state'])} ({', '.join(f'{x:.0f}%' for x in worst.iloc[:3])})"),
        ("National total matches HHS (over 3 million)", 3_000_000 < dme_us < 3_300_000, f"{int(dme_us):,} in the 50 states and DC (state layer {int(state_dme):,})"),
        ("Every GridAI county without an emPOWER match has a stated reason", not unexplained,
         f"{len(unmatched)} unmatched (Connecticut planning regions; Alaska 02063, 02066); unexplained: {', '.join(unexplained) or 'none'}"),
        ("Joined counties hold at least 98% of electricity-dependent residents in the 50 states and DC", dme_joined / dme_us >= 0.98, f"{100 * dme_joined / dme_us:.2f}%"),
    ]
    audit = [{"check": n, "result": "PASS" if ok else "FAIL", "detail": d} for n, ok, d in checks]

    county_rows = j[["fips5", "STATEFP", "NAMELSAD", "medicare_benes", "power_dependent_dme", "dme_rate", "obi", "obi_high",
                     "ejbi", "ejbi_high", "masked_tracts"]].copy()
    return {
        "gate": "PASS" if all(a["result"] == "PASS" for a in audit) else "FAIL",
        "audit": audit,
        "source": {"name": "HHS emPOWER public REST Service, Electricity Dependent DME - ALL - CountyLevel",
                   "publisher": "U.S. Department of Health and Human Services, ASPR, with CMS",
                   "url": "https://services2.arcgis.com/ZQ4jTQn6k7VPXEwO/arcgis/rest/services/HHS_emPOWER_REST_Service_Public/FeatureServer/2",
                   "retrieved": "2026-11-13", "data_last_edited": "2026-11-09",
                   "note": "Medicare Fee-For-Service and Medicare Advantage beneficiaries with a claim for electricity-dependent DME or devices; counts of 1-10 are shown as 11."},
        "summary": {
            "dme_us": int(dme_us), "dme_joined": int(dme_joined), "dme_with_obi": int(dme_obi),
            "counties_joined": int(len(j)), "counties_with_obi": int(len(jo)),
            "obi_high_threshold": round(obi_q75, 4), "counties_obi_high": int(len(hi)),
            "dme_obi_high": int(hi["power_dependent_dme"].sum()),
            "dme_obi_high_pct": round(100 * hi["power_dependent_dme"].sum() / dme_obi, 2),
            "dme_obi_high_and_ejbi_high": int(both_ej["power_dependent_dme"].sum()), "counties_obi_high_and_ejbi_high": int(len(both_ej)),
            "dme_obi_high_and_hidden": int(both_hidden["power_dependent_dme"].sum()), "counties_obi_high_and_hidden": int(len(both_hidden)),
            "rate_obi_high": round(float(rate_hi), 1), "rate_other": round(float(rate_lo), 1),
            "r_rate_obi": round(float(r_obi[0]), 3), "p_rate_obi": float(r_obi[1]), "rho_rate_obi": round(float(rs_obi[0]), 3),
            "county_excess_over_state_pct": round(float(excess_pct), 2), "state_layer_dme_us": int(state_dme),
            "r_rate_ejbi": round(float(r_ej[0]), 3), "p_rate_ejbi": float(r_ej[1]), "n_rate_ejbi": int(len(je)),
        },
        "top": top_rows,
        "counties": county_rows,
    }


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    E = compute(sys.argv[1])
    for a in E["audit"]:
        print(f"{a['result']}  {a['check']}  ({a['detail']})")
    print(json.dumps(E["summary"], indent=1))
    for r in E["top"]:
        print(r)
    if E["gate"] != "PASS":
        sys.exit("emPOWER checks failed.")
