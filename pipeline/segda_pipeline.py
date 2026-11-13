#!/usr/bin/env python3
"""
SEGDA pipeline: reproduces "Socio-Environmental Grid Disruption Analysis"
(Franklin, camera-ready) from its source files and writes the site export.

    python pipeline/segda_pipeline.py <DATA_DIR> <OUT_DIR>

Every number on the site comes from this script. Nothing is imputed: a value
that cannot be computed is written as missing and shown as unavailable.
The run ends with a validation gate that compares each reproduced statistic
with the value locked in the paper; scripts/build_site.py refuses to publish
if a locked value is not reproduced.

Author: Victoria Love Franklin. Copyright (c) 2026. All rights reserved.
"""
import json
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

SEED = 12345
PERMUTATIONS = 9999
BOOTSTRAP = 5000
TRACT_HIGH = 0.75      # study-defined tract High-EJBI threshold (paper 3.2)
COUNTY_HIGH = 0.60     # study-defined county High-EJBI threshold (paper 3.2)
SVI_HIGH = 0.75        # CDC convention for the top vulnerability quartile
YEARS = (2018, 2023)
HAZARDS = {"DRGT_RISKS": "Drought", "HWAV_RISKS": "Heat wave", "WFIR_RISKS": "Wildfire",
           "RFLD_RISKS": "Riverine flood", "ERQK_RISKS": "Earthquake",
           "SWND_RISKS": "Strong wind", "TRND_RISKS": "Tornado"}

# Values locked in the paper. tol is the absolute tolerance for a match;
# "sensitive" marks permutation results that vary between runs by design.
LOCKED = [
    ("Valid tract EJBI", "tracts_valid", 7931, 0, "4.1"),
    ("Tract geometries", "tract_geometries", 8057, 0, "4.1"),
    ("High-EJBI tracts (>= 0.75)", "tracts_high", 1690, 0, "4.1"),
    ("High-EJBI counties (>= 0.60)", "counties_high", 11, 0, "4.1"),
    ("Tract EJBI mean", "tract_mean", 0.500, 0.0005, "Table 2"),
    ("Tract EJBI SD", "tract_sd", 0.255, 0.0005, "Table 2"),
    ("County EJBI mean", "county_mean", 0.437, 0.0005, "Table 2"),
    ("County EJBI SD", "county_sd", 0.144, 0.0005, "Table 2"),
    ("Tract EJBI variance", "var_tract", 0.06487, 0.000005, "Table 3"),
    ("County EJBI variance", "var_county", 0.02088, 0.000005, "Table 3"),
    ("Variance suppressed by county aggregation (%)", "suppressed_pct", 67.82, 0.005, "Table 3"),
    ("High-EJBI tracts inside non-high counties", "masked_tracts", 431, 0, "4.2"),
    ("Counties with Primary OBI", "obi_n", 56, 0, "4.3"),
    ("Primary OBI mean", "obi_mean", 0.275, 0.0005, "Table 2"),
    ("Primary OBI SD", "obi_sd", 0.317, 0.0005, "Table 2"),
    ("Climate Stress mean", "cs_mean", 0.652, 0.0005, "Table 2"),
    ("Climate Stress SD", "cs_sd", 0.149, 0.0005, "Table 2"),
    ("Climate Stress min", "cs_min", 0.125, 0.0005, "Table 2"),
    ("Climate Stress max", "cs_max", 0.886, 0.0005, "Table 2"),
    ("Pearson r, EJBI vs Primary OBI", "r_ej_obi", -0.4419, 0.00005, "Table 4"),
    ("Pearson r, EJBI vs Climate Stress", "r_ej_cs", 0.3936, 0.00005, "Table 4"),
    ("Pearson r, Primary OBI vs Climate Stress", "r_obi_cs", -0.5653, 0.00005, "Table 4"),
    ("Spearman rho, EJBI vs Primary OBI", "rho_ej_obi", -0.5195, 0.00005, "Table 4"),
    ("Spearman rho, EJBI vs Climate Stress", "rho_ej_cs", 0.4840, 0.00005, "Table 4"),
    ("Spearman rho, Primary OBI vs Climate Stress", "rho_obi_cs", -0.5193, 0.00005, "Table 4"),
    ("Moran's I, county EJBI", "moran_ej", 0.2017, 0.00005, "Table 5"),
    ("Moran's I, Primary OBI", "moran_obi", 0.4010, 0.00005, "Table 5"),
    ("County SVI vs EJBI, Pearson r", "r_svi_ej", 0.829, 0.0005, "4.6"),
    ("SVI-EJBI threshold agreement (%)", "svi_agree_pct", 86.2, 0.05, "4.6"),
    ("70/30 EAGLE-I/EIA-861: r EJBI vs OBI", "r_ej_blend", -0.551, 0.0005, "4.6"),
    ("70/30 EAGLE-I/EIA-861: r Climate Stress vs OBI", "r_cs_blend", -0.646, 0.0005, "4.6"),
]
SENSITIVE = [
    ("Moran's I permutation p, county EJBI", "moran_ej_p", 0.006, "Table 5"),
    ("Moran's I permutation p, Primary OBI", "moran_obi_p", 0.002, "Table 5"),
    ("LISA High-High counties, EJBI", "lisa_ej_HH", 4, "Table 5"),
    ("LISA Low-High counties, EJBI", "lisa_ej_LH", 1, "Table 5"),
    ("LISA High-High counties, Primary OBI", "lisa_obi_HH", 6, "Table 5"),
]


COUNTY_SHP = "shp/tl_2019_06_county_from_us.shp"
SOURCES = {
    "CalEnviroScreen_4.0.csv": {
        "name": "CalEnviroScreen 4.0", "publisher": "California Office of Environmental Health Hazard Assessment (OEHHA), CalEPA",
        "year": "2021", "url": "https://oehha.ca.gov/calenviroscreen/report/calenviroscreen-40",
        "use": "Tract CES overall percentile, Pollution Burden percentile, population", "ref": 7},
    "SVI_CA_TRACT.csv": {
        "name": "CDC/ATSDR Social Vulnerability Index, California census tracts", "publisher": "CDC/ATSDR Geospatial Research, Analysis, and Services Program",
        "year": "2018 release (2010 tracts)", "url": "https://www.atsdr.cdc.gov/place-health/php/svi/svi-data-documentation-download.html",
        "use": "Tract SVI overall rank (RPL_THEMES)", "ref": 8},
    "SVI_CA_COUNTY_2022.csv": {
        "name": "CDC/ATSDR Social Vulnerability Index 2022, California counties", "publisher": "CDC/ATSDR Geospatial Research, Analysis, and Services Program",
        "year": "2022", "url": "https://www.atsdr.cdc.gov/place-health/php/svi/svi-data-documentation-download.html",
        "use": "County SVI rank (external check), county population", "ref": 8},
    "eaglei_transformed.csv": {
        "name": "EAGLE-I historic county outage data, reconstructed episodes (California)", "publisher": "DOE CESER / Oak Ridge National Laboratory",
        "year": "2014-2023 (2018-2023 used)", "url": "https://doi.org/10.6084/m9.figshare.24237376",
        "use": "Customer-outage-hours for Primary OBI", "ref": 9},
    "MCC.csv": {
        "name": "Modeled County Customers (MCC)", "publisher": "Oak Ridge National Laboratory (EAGLE-I companion)",
        "year": "2022", "url": "https://doi.org/10.6084/m9.figshare.24237376",
        "use": "Customer denominator for Primary OBI", "ref": 10},
    "NRI_Counties_CA.csv": {
        "name": "National Risk Index, county table (California)", "publisher": "Federal Emergency Management Agency",
        "year": "study vintage", "url": "https://hazards.fema.gov/nri/data-resources",
        "use": "Seven hazard risk scores for Climate Stress", "ref": 12},
    "EIA861_CA_COUNTY_RELIABILITY_6YR.csv": {
        "name": "Form EIA-861 reliability, allocated to California counties (SEGDA derivative)", "publisher": "U.S. Energy Information Administration; county allocation by the SEGDA pipeline",
        "year": "2018-2023", "url": "https://www.eia.gov/electricity/data/eia861/",
        "use": "Historical reliability sensitivity only", "ref": 11},
    "shp/tl_2019_06_tract.shp": {
        "name": "TIGER/Line Shapefiles, California census tracts", "publisher": "U.S. Census Bureau",
        "year": "2019 (2010 tracts)", "url": "https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html",
        "use": "Tract geometry and the 8,057-tract universe", "ref": None},
    COUNTY_SHP: {
        "name": "TIGER/Line Shapefiles, U.S. counties, California rows (STATEFP 06)", "publisher": "U.S. Census Bureau",
        "year": "2019", "url": "https://www2.census.gov/geo/tiger/TIGER2019/COUNTY/tl_2019_us_county.zip",
        "use": "County geometry and contiguity for Moran's I and LISA. Subset of the national file (SHA-256 "
               "0b102f2e292be5a81f46ecb2956136e93da03de692980714fabc16e1925e0b51); rows filtered, geometry unchanged "
               "(checked with exact equality)", "ref": None},
}


def sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def audit_inputs(data, tr, co, ces, svi_t, svi_c, eag, mcc, nri, eia):
    """Checks every input before anything is computed. Any FAIL blocks the build."""
    A = []

    def chk(component, ok, detail, level="FAIL"):
        A.append({"component": component, "status": "PASS" if ok else level, "detail": detail})

    need = {
        "CalEnviroScreen_4.0.csv": (ces, ["Tract", "CIscoreP", "Pollution_Pctl", "Population", "County"]),
        "SVI_CA_TRACT.csv": (svi_t, ["FIPS", "RPL_THEMES", "E_TOTPOP"]),
        "SVI_CA_COUNTY_2022.csv": (svi_c, ["FIPS", "RPL_THEMES", "E_TOTPOP"]),
        "eaglei_transformed.csv": (eag, ["fips", "start_time", "duration", "min_customers", "max_customers", "mean_customers"]),
        "MCC.csv": (mcc, ["County_FIPS", "Customers"]),
        "NRI_Counties_CA.csv": (nri, ["STCOFIPS"] + list(HAZARDS)),
        "EIA861_CA_COUNTY_RELIABILITY_6YR.csv": (eia, ["county", "EIA_RELIABILITY_INDEX_2018_2023"]),
    }
    manifest = []
    for f, meta in SOURCES.items():
        p = data / f
        row = {"file": f, **meta, "present": p.exists()}
        if p.exists():
            row["sha256"] = sha256(p)
            row["bytes"] = p.stat().st_size
        if f in need:
            df, cols = need[f]
            miss = [c for c in cols if c not in df.columns]
            row["rows"] = int(len(df))
            chk(f"Schema: {meta['name']}", not miss, f"required columns present" if not miss else f"missing {miss}")
        manifest.append(row)
        chk(f"File present: {f}", p.exists(), "found" if p.exists() else "not found")
    # coverage and integrity
    chk("Tract universe", len(tr) == 8057, f"{len(tr)} California tract geometries (TIGER 2019)")
    chk("County universe", len(co) == 58, f"{len(co)} California counties")
    chk("CalEnviroScreen tracts", ces["Tract"].nunique() == len(ces) and len(ces) >= 8000, f"{len(ces)} rows, {ces['Tract'].nunique()} unique tracts")
    chk("CalEnviroScreen percentiles in 0-100", ces["CIscoreP"].dropna().between(0, 100).all() and ces["Pollution_Pctl"].dropna().between(0, 100).all(),
        f"{ces['CIscoreP'].isna().sum()} tracts without a CES score (left missing, not imputed)")
    rp = svi_t["RPL_THEMES"]
    chk("SVI tract ranks in 0-1 or suppressed", (rp[rp >= 0] <= 1).all(), f"{int((rp < 0).sum())} suppressed tracts (negative codes) treated as missing")
    chk("SVI tract file is California", svi_t["ST_ABBR"].eq("CA").all() if "ST_ABBR" in svi_t else False, f"{len(svi_t)} rows")
    chk("SVI county file covers 58 counties", svi_c["FIPS"].nunique() == 58, f"{svi_c['FIPS'].nunique()} counties")
    chk("NRI covers 58 counties", nri["STCOFIPS"].nunique() == 58, f"{nri['STCOFIPS'].nunique()} counties")
    chk("NRI hazard scores in 0-100", nri[list(HAZARDS)].stack().between(0, 100).all(), "all seven hazard risk scores in range")
    chk("EAGLE-I file is California only", eag["state"].eq("California").all(), f"{len(eag)} rows, {eag['fips'].nunique()} counties")
    chk("EAGLE-I durations positive", (eag["duration"] > 0).all(), f"duration {eag['duration'].min()} to {eag['duration'].max()} hours")
    chk("EAGLE-I customers non-negative", (eag["mean_customers"] >= 0).all(), f"min customers {eag['min_customers'].min()}")
    mc = mcc[mcc["County_FIPS"].astype(str) != "Grand Total"]
    ca = mc[mc["County_FIPS"].astype(int).between(6001, 6115)]
    chk("MCC California counties", len(ca) >= 56, f"{len(ca)} California counties with a customer count (national file, {len(mc)} counties)")
    chk("MCC customers positive", (ca["Customers"] > 0).all(), "all California customer counts > 0")
    miss_eia = sorted(set(co["NAME"]) - set(eia["county"]))
    chk("EIA-861 county table", not miss_eia,
        f"{eia['county'].nunique()} of 58 counties" + (f"; absent: {', '.join(miss_eia)} (left missing). EIA-861 is a sensitivity branch only, so this does not block the build" if miss_eia else ""),
        level="WARNING")
    return A, manifest


def r6(x):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 6)


def load(data):
    tr = gpd.read_file(data / "shp" / "tl_2019_06_tract.shp")
    tr["tract"] = tr["GEOID"].astype("int64")
    cp = data / COUNTY_SHP
    if not cp.exists():  # the full national TIGER file also works; it is filtered to California below
        cp = data / "shp" / "tl_2019_us_county.shp"
    co = gpd.read_file(cp)
    co = co[co["STATEFP"] == "06"].copy()
    co["fips"] = co["GEOID"].astype(int)
    ces = pd.read_csv(data / "CalEnviroScreen_4.0.csv", low_memory=False, encoding="utf-8-sig")
    svi_t = pd.read_csv(data / "SVI_CA_TRACT.csv", low_memory=False, encoding="utf-8-sig")
    svi_c = pd.read_csv(data / "SVI_CA_COUNTY_2022.csv", encoding="utf-8-sig")
    eag = pd.read_csv(data / "eaglei_transformed.csv", low_memory=False)
    mcc = pd.read_csv(data / "MCC.csv")
    nri = pd.read_csv(data / "NRI_Counties_CA.csv", dtype={"STCOFIPS": str})
    eia = pd.read_csv(data / "EIA861_CA_COUNTY_RELIABILITY_6YR.csv")
    return tr, co, ces, svi_t, svi_c, eag, mcc, nri, eia


def tract_ejbi(tr, ces, svi_t):
    s = svi_t.assign(tract=svi_t["FIPS"].astype("int64"))
    s["svi"] = s["RPL_THEMES"].where(s["RPL_THEMES"] >= 0)   # negative = suppressed
    c = ces.rename(columns={"Tract": "tract"})[["tract", "CIscoreP", "Pollution_Pctl", "Population"]]
    d = (tr[["tract", "GEOID", "COUNTYFP"]]
         .merge(s[["tract", "svi"]], on="tract", how="left")
         .merge(c, on="tract", how="left"))
    d["ces"] = d["CIscoreP"] / 100
    d["pb"] = d["Pollution_Pctl"] / 100
    d["ejbi"] = pd.concat([d["svi"], d["ces"], d["pb"]], axis=1).mean(axis=1, skipna=False)
    d["fips"] = 6000 + d["COUNTYFP"].astype(int)
    return d


def county_ejbi(d):
    v = d.dropna(subset=["ejbi"])
    v = v[v["Population"] > 0]
    return v.groupby("fips").apply(lambda g: np.average(g["ejbi"], weights=g["Population"]))


def primary_obi(eag, mcc):
    e = eag.drop_duplicates(["fips", "start_time", "duration", "min_customers",
                             "max_customers", "mean_customers"]).copy()
    n_raw, n_dedup = len(eag), len(e)
    e["start"] = pd.to_datetime(e["start_time"], format="%m/%d/%Y %H:%M", errors="coerce")
    e = e[e["start"].dt.year.between(*YEARS)].copy()
    e["coh"] = e["mean_customers"] * e["duration"]
    m = mcc[mcc["County_FIPS"].astype(str) != "Grand Total"].copy()
    m["fips"] = m["County_FIPS"].astype(int)
    r = (e.groupby("fips").agg(coh=("coh", "sum"), episodes=("coh", "size"))
          .join(m.set_index("fips")["Customers"], how="inner"))
    r["R"] = r["coh"] / r["Customers"]
    q1, q3 = r["R"].quantile([0.25, 0.75])
    rc = r["R"].clip(q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1))
    r["obi"] = (rc - rc.min()) / (rc.max() - rc.min())
    meta = {"rows_raw": n_raw, "rows_after_dedup": n_dedup, "rows_2018_2023": int(len(e)),
            "counties_with_outages": int(e["fips"].nunique())}
    return r, meta


def robust_norm(r):
    q1, q3 = r.quantile([0.25, 0.75])
    rc = r.clip(q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1))
    return (rc - rc.min()) / (rc.max() - rc.min())


def obi_variants(eag, mcc, county_pop):
    """Sensitivity specifications the paper retains (3.3, 4.6): population
    denominator instead of MCC, and excluding episodes at the 34-hour ceiling
    or the 200-customer floor of the upstream reconstruction."""
    e = eag.drop_duplicates(["fips", "start_time", "duration", "min_customers",
                             "max_customers", "mean_customers"]).copy()
    e["start"] = pd.to_datetime(e["start_time"], format="%m/%d/%Y %H:%M", errors="coerce")
    e = e[e["start"].dt.year.between(*YEARS)].copy()
    e["coh"] = e["mean_customers"] * e["duration"]
    m = mcc[mcc["County_FIPS"].astype(str) != "Grand Total"].copy()
    m["fips"] = m["County_FIPS"].astype(int)
    cust = m.set_index("fips")["Customers"]
    out = {}
    base = e.groupby("fips")["coh"].sum()
    out["Population denominator"] = robust_norm((base / county_pop).dropna().loc[lambda x: x.index.isin(cust.index)])
    ceil = e["duration"] >= e["duration"].max()
    out["Exclude 34-hour ceiling episodes"] = robust_norm((e[~ceil].groupby("fips")["coh"].sum() / cust).dropna())
    floor = e["min_customers"] <= e["min_customers"].min()
    out["Exclude 200-customer floor episodes"] = robust_norm((e[~floor].groupby("fips")["coh"].sum() / cust).dropna())
    out["Exclude both"] = robust_norm((e[~ceil & ~floor].groupby("fips")["coh"].sum() / cust).dropna())
    meta = {"ceiling_hours": float(e["duration"].max()), "ceiling_episodes": int(ceil.sum()),
            "floor_customers": int(e["min_customers"].min()), "floor_episodes": int(floor.sum())}
    return out, meta


def climate_stress(nri):
    n = nri.assign(fips=nri["STCOFIPS"].astype(int)).set_index("fips")
    norm = n[list(HAZARDS)].apply(lambda s: (s - s.min()) / (s.max() - s.min()))
    return norm.mean(axis=1), norm


def boot_ci(x, y, n=BOOTSTRAP):
    rng = np.random.default_rng(SEED)
    x, y = np.asarray(x), np.asarray(y)
    out = []
    for _ in range(n):
        i = rng.integers(0, len(x), len(x))
        if np.std(x[i]) > 0 and np.std(y[i]) > 0:
            out.append(stats.pearsonr(x[i], y[i])[0])
    return [r6(np.percentile(out, 2.5)), r6(np.percentile(out, 97.5))]


def spatial(co, col):
    from libpysal.weights import Queen
    from esda.moran import Moran, Moran_Local
    sub = co[co[col].notna()].sort_values("fips").set_index("fips")
    w = Queen.from_dataframe(sub, use_index=True, silence_warnings=True)
    w.transform = "r"
    np.random.seed(SEED)
    m = Moran(sub[col].values, w, permutations=PERMUTATIONS)
    np.random.seed(SEED)
    loc = Moran_Local(sub[col].values, w, permutations=PERMUTATIONS)
    lab = {1: "HH", 2: "LH", 3: "LL", 4: "HL"}
    cl = pd.Series([lab[q] if p < 0.05 else "NS" for q, p in zip(loc.q, loc.p_sim)], index=sub.index)
    return m, cl


def main(data_dir, out_dir):
    data, out = Path(data_dir), Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tr, co, ces, svi_t, svi_c, eag, mcc, nri, eia = load(data)
    audit, manifest = audit_inputs(data, tr, co, ces, svi_t, svi_c, eag, mcc, nri, eia)
    n_audit_fail = sum(a["status"] == "FAIL" for a in audit)
    n_audit_warn = sum(a["status"] == "WARNING" for a in audit)

    # ---------- tracts ----------
    d = tract_ejbi(tr, ces, svi_t)
    v = d.dropna(subset=["ejbi"])
    cej = county_ejbi(d)
    d["county_ejbi"] = d["fips"].map(cej)
    d["high"] = d["ejbi"] >= TRACT_HIGH
    d["masked"] = d["high"] & (d["county_ejbi"] < COUNTY_HIGH)
    d["svi_high"] = d["svi"] >= SVI_HIGH

    # ---------- counties ----------
    obi, obi_meta = primary_obi(eag, mcc)
    cs, cs_parts = climate_stress(nri)
    sc = svi_c.assign(fips=svi_c["FIPS"].astype(int)).set_index("fips")
    co = co.set_index("fips")
    co["ejbi"] = cej
    co["obi"] = obi["obi"]
    co["obi_raw"] = obi["R"]
    co["cs"] = cs
    co["svi22"] = sc["RPL_THEMES"].where(sc["RPL_THEMES"] >= 0)
    name_to_fips = co.reset_index().set_index("NAME")["fips"]
    e6 = eia.assign(fips=eia["county"].map(name_to_fips)).set_index("fips")
    co["eia_index"] = e6["EIA_RELIABILITY_INDEX_2018_2023"]
    co["blend"] = 0.7 * co["obi"] + 0.3 * co["eia_index"]
    co = co.reset_index()

    # Per-county view of what the county figure hides.
    g = d.dropna(subset=["ejbi"]).groupby("fips")
    pop = d.dropna(subset=["ejbi"]).assign(hp=lambda x: x["Population"].where(x["high"], 0))
    hid = pd.DataFrame({
        "tracts_valid": g.size(),
        "tracts_high": g["high"].sum(),
        "tract_min": g["ejbi"].min(), "tract_max": g["ejbi"].max(),
        "tract_p10": g["ejbi"].quantile(0.10), "tract_p90": g["ejbi"].quantile(0.90),
        "pop_high": pop.groupby("fips")["hp"].sum(),
        "pop_total": pop.groupby("fips")["Population"].sum(),
    })
    sv = d.dropna(subset=["svi"]).groupby("fips")
    hid["svi_tracts"] = sv.size()
    hid["svi_tracts_high"] = sv["svi_high"].sum()
    hid["svi_tract_popmean"] = d.dropna(subset=["svi"]).query("Population > 0").groupby("fips").apply(
        lambda x: np.average(x["svi"], weights=x["Population"]))
    co = co.merge(hid.reset_index(), on="fips", how="left")

    # ---------- statistics ----------
    S = {}
    S["tract_geometries"] = int(len(d))
    S["tracts_valid"] = int(len(v))
    S["tracts_high"] = int((v["ejbi"] >= TRACT_HIGH).sum())
    S["tracts_high_pct"] = r6(100 * (v["ejbi"] >= TRACT_HIGH).mean())
    S["counties_high"] = int((cej >= COUNTY_HIGH).sum())
    S["tract_mean"], S["tract_sd"] = r6(v["ejbi"].mean()), r6(v["ejbi"].std())
    S["tract_min"], S["tract_max"] = r6(v["ejbi"].min()), r6(v["ejbi"].max())
    S["county_mean"], S["county_sd"] = r6(cej.mean()), r6(cej.std())
    S["county_min"], S["county_max"] = r6(cej.min()), r6(cej.max())
    S["var_tract"], S["var_county"] = r6(v["ejbi"].var()), r6(cej.var())
    S["suppressed_pct"] = r6(100 * (1 - cej.var() / v["ejbi"].var()))
    S["retained_pct"] = r6(100 * cej.var() / v["ejbi"].var())
    S["masked_tracts"] = int(d["masked"].sum())
    S["masked_pop"] = int(d.loc[d["masked"], "Population"].sum())
    S["masked_counties"] = int(d.loc[d["masked"], "fips"].nunique())
    o = co.dropna(subset=["obi"])
    S["obi_n"] = int(len(o))
    S["obi_mean"], S["obi_sd"] = r6(o["obi"].mean()), r6(o["obi"].std())
    S["obi_unavailable"] = sorted(co.loc[co["obi"].isna(), "NAME"].tolist())
    S["cs_mean"], S["cs_sd"] = r6(co["cs"].mean()), r6(co["cs"].std())
    S["cs_min"], S["cs_max"] = r6(co["cs"].min()), r6(co["cs"].max())
    pairs = {"ej_obi": ("ejbi", "obi"), "ej_cs": ("ejbi", "cs"), "obi_cs": ("obi", "cs")}
    corr = []
    for k, (a, b) in pairs.items():
        x = co.dropna(subset=[a, b])
        pr, pp = stats.pearsonr(x[a], x[b])
        sr, sp = stats.spearmanr(x[a], x[b])
        S[f"r_{k}"], S[f"rho_{k}"] = r6(pr), r6(sr)
        corr.append({"pair": k, "x": a, "y": b, "n": int(len(x)), "pearson_r": r6(pr), "pearson_p": r6(pp),
                     "spearman_rho": r6(sr), "spearman_p": r6(sp), "boot95": boot_ci(x[a], x[b])})
    gco = gpd.GeoDataFrame(co, geometry="geometry", crs=tr.crs)
    m_ej, l_ej = spatial(gco, "ejbi")
    m_ob, l_ob = spatial(gco, "obi")
    S["moran_ej"], S["moran_ej_p"] = r6(m_ej.I), r6(m_ej.p_sim)
    S["moran_obi"], S["moran_obi_p"] = r6(m_ob.I), r6(m_ob.p_sim)
    for nm, lc in (("ej", l_ej), ("obi", l_ob)):
        for k in ("HH", "LH", "LL", "HL"):
            S[f"lisa_{nm}_{k}"] = int((lc == k).sum())
    co["lisa_ej"] = co["fips"].map(l_ej)
    co["lisa_obi"] = co["fips"].map(l_ob)
    # Sensitivity: alternative outage specifications (direction must not reverse).
    cpop = sc["E_TOTPOP"]
    variants, vmeta = obi_variants(eag, mcc, cpop)
    variants = {"Primary OBI (EAGLE-I / MCC)": co.set_index("fips")["obi"].dropna(), **variants,
                "70% EAGLE-I / 30% EIA-861": co.set_index("fips")["blend"].dropna()}
    sens = []
    cix = co.set_index("fips")
    for nm, ser in variants.items():
        j = pd.DataFrame({"o": ser}).join(cix[["ejbi", "cs"]]).dropna()
        r1, p1 = stats.pearsonr(j["ejbi"], j["o"])
        r2, p2 = stats.pearsonr(j["cs"], j["o"])
        sens.append({"spec": nm, "n": int(len(j)), "r_ejbi": r6(r1), "p_ejbi": r6(p1),
                     "ci_ejbi": boot_ci(j["ejbi"], j["o"]), "r_cs": r6(r2), "p_cs": r6(p2)})
    S["sens_direction_holds"] = all(x["r_ejbi"] < 0 for x in sens)
    S["sens_meta"] = vmeta
    # Welch's unequal-variance test: Primary OBI in High-EJBI vs other counties.
    hi = co.loc[(co["ejbi"] >= COUNTY_HIGH) & co["obi"].notna(), "obi"]
    lo = co.loc[(co["ejbi"] < COUNTY_HIGH) & co["obi"].notna(), "obi"]
    wt = stats.ttest_ind(hi, lo, equal_var=False)
    S["welch"] = {"n_high": int(len(hi)), "n_other": int(len(lo)), "mean_high": r6(hi.mean()),
                  "mean_other": r6(lo.mean()), "t": r6(wt.statistic), "p": r6(wt.pvalue)}
    x = co.dropna(subset=["svi22", "ejbi"])
    S["r_svi_ej"] = r6(stats.pearsonr(x["svi22"], x["ejbi"])[0])
    S["svi_agree_pct"] = r6(100 * ((x["ejbi"] >= COUNTY_HIGH) == (x["svi22"] >= SVI_HIGH)).mean())
    b = co.dropna(subset=["blend", "ejbi", "cs"])
    S["blend_n"] = int(len(b))
    S["r_ej_blend"] = r6(stats.pearsonr(b["ejbi"], b["blend"])[0])
    S["r_cs_blend"] = r6(stats.pearsonr(b["cs"], b["blend"])[0])
    # SVI lens: high-vulnerability tracts that the official county rank does not flag.
    tmp = d.dropna(subset=["svi"]).assign(c22=lambda z: z["fips"].map(co.set_index("fips")["svi22"]))
    S["svi_tracts_valid"] = int(len(tmp))
    S["svi_tracts_high"] = int(tmp["svi_high"].sum())
    S["svi_hidden_tracts"] = int((tmp["svi_high"] & (tmp["c22"] < SVI_HIGH)).sum())
    S["svi_hidden_pop"] = int(tmp.loc[tmp["svi_high"] & (tmp["c22"] < SVI_HIGH), "Population"].fillna(0).sum())
    S["svi_counties_high"] = int((co["svi22"] >= SVI_HIGH).sum())

    # ---------- validation gate ----------
    checks = []
    for label, key, paper, tol, where in LOCKED:
        got = S.get(key)
        ok = got is not None and abs(float(got) - float(paper)) <= tol
        checks.append({"check": label, "paper": paper, "reproduced": got, "section": where,
                       "status": "MATCH" if ok else "MISMATCH"})
    for label, key, paper, where in SENSITIVE:
        checks.append({"check": label, "paper": paper, "reproduced": S.get(key), "section": where,
                       "status": "PERMUTATION-SENSITIVE"})
    n_mis = sum(c["status"] == "MISMATCH" for c in checks)

    run = {"run_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
           "seed": SEED, "permutations": PERMUTATIONS, "bootstrap": BOOTSTRAP,
           "gate": "PASS" if n_mis == 0 and n_audit_fail == 0 else "FAIL", "n_locked": len(LOCKED), "n_mismatch": n_mis,
           "n_audit": len(audit), "n_audit_fail": n_audit_fail, "n_audit_warn": n_audit_warn,
           "eaglei": obi_meta,
           "vintages": {"tract_svi": "CDC/ATSDR SVI tract file supplied (2018 schema: E_POV, E_PCI; 2010 tracts)",
                        "county_svi": "CDC/ATSDR SVI 2022, county", "ces": "CalEnviroScreen 4.0 (2021)",
                        "tracts": "Census TIGER/Line 2019 tracts (2010 vintage)", "nri": "FEMA NRI county table supplied",
                        "eaglei": "EAGLE-I reconstructed episodes 2018-2023 (supplied transformed file)",
                        "mcc": "ORNL Modeled County Customers", "eia": "EIA-861 2018-2023, county allocation (supplied)"}}

    # ---------- write ----------
    tcols = ["GEOID", "fips", "Population", "svi", "ces", "pb", "ejbi", "high", "masked", "svi_high"]
    d[tcols].to_csv(out / "segda_tracts.csv", index=False)
    ccols = ["fips", "NAME", "ejbi", "obi", "obi_raw", "cs", "svi22", "eia_index", "blend", "lisa_ej", "lisa_obi",
             "tracts_valid", "tracts_high", "tract_min", "tract_max", "tract_p10", "tract_p90", "pop_high",
             "pop_total", "svi_tracts", "svi_tracts_high", "svi_tract_popmean"]
    co[ccols].to_csv(out / "segda_counties.csv", index=False)
    cs_parts.rename(columns=HAZARDS).reset_index().to_csv(out / "segda_climate_components.csv", index=False)
    pd.DataFrame(checks).to_csv(out / "segda_validation.csv", index=False)
    pd.DataFrame(audit).to_csv(out / "segda_audit.csv", index=False)
    pd.DataFrame(manifest).to_csv(out / "segda_sources.csv", index=False)
    (out / "segda_results.json").write_text(json.dumps(
        {"run": run, "stats": S, "correlations": corr, "sensitivity": sens, "validation": checks,
         "audit": audit, "sources": manifest, "hazards": HAZARDS,
         "thresholds": {"tract": TRACT_HIGH, "county": COUNTY_HIGH, "svi": SVI_HIGH}}, indent=1, default=str))

    # ---------- geometry for the maps ----------
    geo = out / "geo"
    geo.mkdir(exist_ok=True)
    t3310 = tr.to_crs(3310)[["GEOID", "geometry"]]
    t3310 = t3310[~t3310["GEOID"].str.endswith("990000") | t3310["GEOID"].isin(d.dropna(subset=["ejbi"])["GEOID"])]
    t3310.to_file(geo / "tracts.geojson", driver="GeoJSON")
    co.set_geometry("geometry").set_crs(tr.crs, allow_override=True).to_crs(3310)[["fips", "NAME", "geometry"]] \
        .to_file(geo / "counties.geojson", driver="GeoJSON")

    print(f"gate {run['gate']}: {len(LOCKED) - n_mis}/{len(LOCKED)} locked values reproduced; "
          f"input audit {len(audit) - n_audit_fail - n_audit_warn} pass, {n_audit_warn} warning, {n_audit_fail} fail")
    for a in audit:
        print(f"  {a['status']:5s} {a['component']}: {a['detail']}")
    for c in checks:
        print(f"  {c['status']:22s} {c['check']}: paper {c['paper']} | reproduced {c['reproduced']}")
    return 0 if run["gate"] == "PASS" else 2


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1], sys.argv[2]))
