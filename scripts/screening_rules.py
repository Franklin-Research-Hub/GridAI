#!/usr/bin/env python3
"""
GridAI (SEGDA), national layer: how many high-burden residents does each county screening rule reach?

    python scripts/screening_rules.py <NATIONAL_BUILD_DIR>

Reads the audited national tables (national_tracts.csv, national_counties.csv, national_results.json),
compares county screening rules against a direct tract screen, and checks the result against the
published national totals. Any failed check stops the run, so the site is never built on numbers that
disagree with the audited national layer. Nothing is estimated: every figure is a count of tracts,
counties or residents in the published tables.

Author: Victoria Love Franklin. Copyright (c) 2026. All rights reserved.
"""
import json
import sys
from pathlib import Path

import pandas as pd


def compute(nb):
    nb = Path(nb)
    R = json.loads((nb / "national_results.json").read_text())
    thr = R["thresholds"]
    t = pd.read_csv(nb / "national_tracts.csv", dtype={"GEOID": str, "st": str, "fips": str})
    c = pd.read_csv(nb / "national_counties.csv", dtype={"fips5": str, "STATEFP": str})
    t["pop"] = t["pop"].fillna(0)

    t["hi"] = t["ejbi"] >= thr["tract"]
    t["svi_hi"] = t["svi"] >= thr["tract"]
    g = t.groupby("fips")
    cty = pd.DataFrame({
        "pop": g["pop"].sum(),
        "hi_pop": t[t["hi"]].groupby("fips")["pop"].sum(),
        "hi_tracts": g["hi"].sum(),
        "svi_hi_pop": t[t["svi_hi"]].groupby("fips")["pop"].sum(),
        "svi_hi_tracts": g["svi_hi"].sum(),
    }).fillna(0)
    cty = cty.join(c.set_index("fips5")[["ejbi", "svi_c", "NAMELSAD", "STATEFP"]])
    cty["hi_share"] = cty["hi_pop"] / cty["pop"].where(cty["pop"] > 0)

    tot_pop = float(cty["pop"].sum())
    tot_hi = float(cty["hi_pop"].sum())
    tot_svi = float(cty["svi_hi_pop"].sum())
    n_high_c = int((cty["ejbi"] >= thr["county"]).sum())

    def rule(key, name, flagged, target="ejbi", direct=False):
        if direct:
            m = t["hi"] if target == "ejbi" else t["svi_hi"]
            flagged_pop, reached, n_c, n_t = float(t.loc[m, "pop"].sum()), float(t.loc[m, "pop"].sum()), int(t.loc[m, "fips"].nunique()), int(m.sum())
        else:
            col = "hi_pop" if target == "ejbi" else "svi_hi_pop"
            flagged_pop = float(cty.loc[flagged, "pop"].sum())
            reached = float(cty.loc[flagged, col].sum())
            n_c, n_t = int(flagged.sum()), None
        total = tot_hi if target == "ejbi" else tot_svi
        return {"key": key, "rule": name, "index": "EJBI" if target == "ejbi" else "CDC SVI",
                "counties_flagged": n_c, "tracts_flagged": n_t,
                "population_flagged": int(round(flagged_pop)), "high_burden_reached": int(round(reached)),
                "high_burden_missed": int(round(total - reached)), "share_reached_pct": round(100 * reached / total, 2),
                "share_of_flagged_pop_high_burden_pct": round(100 * reached / flagged_pop, 2) if flagged_pop else None}

    rank = cty["hi_pop"].rank(method="first", ascending=False)
    rules = [
        rule("county_avg", f"County average EJBI ≥ {thr['county']:.2f} (current county rule)", cty["ejbi"] >= thr["county"]),
        rule("same_count", f"Same number of counties ({n_high_c}), ranked by high-burden residents", rank <= n_high_c),
        rule("share25", "County with ≥ 25% of residents in high-burden tracts", cty["hi_share"] >= 0.25),
        rule("share10", "County with ≥ 10% of residents in high-burden tracts", cty["hi_share"] >= 0.10),
        rule("any_tract", "County containing any high-burden tract", cty["hi_tracts"] > 0),
        rule("tract", "Tract screen: flag the high-burden tracts themselves", None, direct=True),
        rule("svi_county", f"CDC county SVI ≥ {thr['svi_county']:.2f}", cty["svi_c"] >= thr["svi_county"], target="svi"),
        rule("svi_any_tract", "County containing any top-quartile SVI tract", cty["svi_hi_tracts"] > 0, target="svi"),
    ]

    nat = R["national"]
    A = rules[0]
    checks = [
        ("Population across all tracts equals the national total", int(round(tot_pop)) == int(nat["pop"]), f"{int(round(tot_pop)):,} vs {int(nat['pop']):,}"),
        # The county table stores pop_high after rounding, so a handful of residents can differ; the tract table is exact.
        ("Residents of high-burden tracts agree with the county table (within 100 people of rounding)",
         abs(tot_hi - float(c["pop_high"].sum())) <= 100, f"tracts {int(round(tot_hi)):,} vs county table {int(round(float(c['pop_high'].sum()))):,}"),
        ("High-burden tract count matches the national results", int(t["hi"].sum()) == int(nat["tracts_high"]), f"{int(t['hi'].sum()):,} vs {int(nat['tracts_high']):,}"),
        ("Current county rule flags the published number of high-burden counties", A["counties_flagged"] == int(nat["counties_high"]), f"{A['counties_flagged']} vs {nat['counties_high']}"),
        ("Current county rule misses exactly the published hidden population", A["high_burden_missed"] == int(nat["masked_pop"]), f"{A['high_burden_missed']:,} vs {int(nat['masked_pop']):,}"),
        ("CDC county SVI rule misses exactly the published SVI-hidden population", rules[6]["high_burden_missed"] == int(nat["svi_hidden_pop"]), f"{rules[6]['high_burden_missed']:,} vs {int(nat['svi_hidden_pop']):,}"),
        ("Any-tract rule reaches every high-burden resident", rules[4]["high_burden_missed"] == 0, f"{rules[4]['high_burden_missed']:,} missed"),
        ("Tract screen reaches every high-burden resident", rules[5]["high_burden_missed"] == 0, f"{rules[5]['high_burden_missed']:,} missed"),
    ]
    audit = [{"check": n, "result": "PASS" if ok else "FAIL", "detail": d} for n, ok, d in checks]
    return {"rules": rules, "audit": audit, "gate": "PASS" if all(a["result"] == "PASS" for a in audit) else "FAIL",
            "totals": {"population": int(round(tot_pop)), "high_burden_residents": int(round(tot_hi)),
                       "high_svi_residents": int(round(tot_svi))},
            "thresholds": thr}


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    S = compute(sys.argv[1])
    for a in S["audit"]:
        print(f"{a['result']}  {a['check']}  ({a['detail']})")
    for r in S["rules"]:
        print(f"{r['index']:8} {r['rule'][:62]:62} counties {r['counties_flagged']:>5}  reached {r['share_reached_pct']:6.2f}%  missed {r['high_burden_missed']:>11,}  pop flagged {r['population_flagged']:>12,}  in-target {r['share_of_flagged_pop_high_burden_pct']}%")
    if S["gate"] != "PASS":
        sys.exit("Screening-rule checks failed.")
