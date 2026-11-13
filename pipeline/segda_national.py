#!/usr/bin/env python3
"""
GridAI (SEGDA), national layer: the paper's MAUP methods applied to all 50 states and DC.

    python pipeline/segda_national.py <NATIONAL_INPUTS> <CA_STATIC_INPUTS> <CA_PIPELINE_OUT> <OUT_DIR>

NATIONAL_INPUTS holds the official files exactly as published (fingerprinted below).
CA_STATIC_INPUTS supplies the national MCC table. CA_PIPELINE_OUT is the validated
California run, used only for the concordance check against the paper.

What changes from the paper, and why (also stated on the site):
  * CalEnviroScreen exists only for California. The national burden index (national EJBI) is
    mean(CDC SVI 2022 RPL_THEMES, CDC/ATSDR EJI 2024 RPL_EJI, EJI 2024 RPL_EBM), the same
    three-part structure: social vulnerability, a cumulative-impact index, and environmental burden.
  * EJI 2024 publishes no ranks for Alaska and Hawaii, so national EJBI is left missing there
    (not imputed). SVI-only views still cover those states.
  * Outages come from the public OEDI release of the EAGLE-I data (PNNL submission 6458, v2),
    processed exactly as in the paper. The paper's California file is not byte-identical to this
    release; the California concordance check below quantifies the difference.
  * Connecticut: tract and SVI/EJI/TIGER 2022 data use the 9 planning regions; EAGLE-I and MCC use
    the 8 former counties. They do not nest, so Connecticut OBI is left missing.

Every value is computed from the inputs. Nothing is simulated or filled in.
Author: Victoria Love Franklin. Copyright (c) 2026. All rights reserved.
"""
import hashlib
import json
import sys
import time
import warnings
import zipfile
from io import BytesIO
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

T_TRACT, T_COUNTY, T_SVI_COUNTY = 0.75, 0.60, 0.75
YEARS = (2018, 2023)
HAZARDS = {"DRGT_RISKS": "Drought", "HWAV_RISKS": "Heat wave", "WFIR_RISKS": "Wildfire",
           "IFLD_RISKS": "Inland flooding (NRI v1.20 successor to Riverine flooding)", "ERQK_RISKS": "Earthquake", "SWND_RISKS": "Strong wind",
           "TRND_RISKS": "Tornado"}
STATES = ("01 02 04 05 06 08 09 10 11 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 "
          "34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56").split()
SEED = 12345

SOURCES = {
    "SVI_2022_US.csv": dict(
        name="CDC/ATSDR Social Vulnerability Index 2022, U.S. census tracts (national ranking)",
        publisher="CDC/ATSDR Geospatial Research, Analysis, and Services Program", year="2022 (ACS 2018-2022; 2020 tracts)",
        url="https://svi.cdc.gov/Documents/Data/2022/csv/states/SVI_2022_US.csv",
        sha256="58df3dca3d4194879ffbd4b5e753eae48676ef691da1fa669f5fc40c3f79bc6a",
        use="Tract SVI rank (RPL_THEMES) and tract population (E_TOTPOP)"),
    "SVI_2022_US_county.csv": dict(
        name="CDC/ATSDR Social Vulnerability Index 2022, U.S. counties (national ranking)",
        publisher="CDC/ATSDR Geospatial Research, Analysis, and Services Program", year="2022",
        url="https://svi.cdc.gov/Documents/Data/2022/csv/states_counties/SVI_2022_US_county.csv",
        sha256="bc47d244153e359d5c09f621a4bf344a1e593159c16fbd3c28a15461b70a6c0f",
        use="County SVI rank: what a county-level screen sees"),
    "EJI_2024_United_States_CSV.zip": dict(
        name="CDC/ATSDR Environmental Justice Index 2024, United States",
        publisher="CDC/ATSDR Geospatial Research, Analysis, and Services Program", year="2024 (2020 tracts)",
        url="https://eji.cdc.gov/Documents/Data/2024/2024_CSV/EJI_2024_United_States_CSV.zip",
        sha256="cd483729f44b28358376240be8bca0ff985916acd833dece180357609aaba297",
        use="Tract EJI rank (RPL_EJI) and Environmental Burden Module rank (RPL_EBM)"),
    "OEDI6458_Outage_Dataset_R1.zip": dict(
        name="Event-correlated Outage Dataset in America (EAGLE-I merged outages, v2)",
        publisher="Pacific Northwest National Laboratory, via DOE Open Energy Data Initiative", year="2014-2023 (2018-2023 used)",
        url="https://data.openei.org/files/6458/Outage_Dataset_R1.zip",
        sha256="44557bae7a5c8e02c611a14f75fe9e123d81d286f700b9fd88ac538edd910fe7",
        use="County customer-outage-hours for Primary OBI"),
    "NRI_Table_Counties_v120.zip": dict(
        name="National Risk Index, county table, v1.20",
        publisher="Federal Emergency Management Agency", year="December 2025 (v1.20)",
        url="https://www.fema.gov/about/reports-and-data/openfema/nri/v120/NRI_Table_Counties.zip",
        sha256="b99dae90a317b16630d7f9b1faa4a2784d5ca458761c1d3ab9646f046df9d072",
        use="Seven hazard risk scores for Climate Stress"),
    "tiger2022/tl_2022_us_county.zip": dict(
        name="TIGER/Line Shapefiles 2022, U.S. counties", publisher="U.S. Census Bureau", year="2022",
        url="https://www2.census.gov/geo/tiger/TIGER2022/COUNTY/tl_2022_us_county.zip", sha256=None,
        use="County geometry and contiguity"),
    "tiger2022/tl_2022_us_state.zip": dict(
        name="TIGER/Line Shapefiles 2022, U.S. states", publisher="U.S. Census Bureau", year="2022",
        url="https://www2.census.gov/geo/tiger/TIGER2022/STATE/tl_2022_us_state.zip", sha256=None,
        use="State outlines"),
}
for s in STATES:
    SOURCES[f"tiger2022/tl_2022_{s}_tract.zip"] = dict(
        name=f"TIGER/Line Shapefiles 2022, census tracts, state {s}", publisher="U.S. Census Bureau", year="2022 (2020 tracts)",
        url=f"https://www2.census.gov/geo/tiger/TIGER2022/TRACT/tl_2022_{s}_tract.zip", sha256=None,
        use="Tract geometry and the tract universe")

AUDIT = []


def chk(name, ok, detail, level="FAIL"):
    AUDIT.append({"check": name, "result": "PASS" if ok else level, "detail": detail})
    return ok


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def r6(x):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), 6)


def boot(x, y, n=5000):
    rng = np.random.default_rng(SEED)
    x, y = np.asarray(x), np.asarray(y)
    rs = []
    for _ in range(n):
        i = rng.integers(0, len(x), len(x))
        if np.std(x[i]) > 0 and np.std(y[i]) > 0:
            rs.append(np.corrcoef(x[i], y[i])[0, 1])
    return float(np.percentile(rs, 2.5)), float(np.percentile(rs, 97.5))


def corr(df, a, b, ci=True):
    d = df[[a, b]].dropna()
    if len(d) < 5:
        return None
    pr, pp = stats.pearsonr(d[a], d[b])
    sr, sp = stats.spearmanr(d[a], d[b])
    out = {"n": int(len(d)), "pearson": r6(pr), "p": float(pp), "spearman": r6(sr), "spearman_p": float(sp)}
    if ci:
        out["ci"] = [r6(v) for v in boot(d[a], d[b])]
    return out


def read_zip_csv(p, member_suffix, **kw):
    with zipfile.ZipFile(p) as z:
        name = [n for n in z.namelist() if n.endswith(member_suffix)][0]
        return pd.read_csv(BytesIO(z.read(name)), **kw)


def load_outages(p):
    frames = []
    with zipfile.ZipFile(p) as z:
        for n in sorted(z.namelist()):
            base = n.split("/")[-1]
            if base.endswith("_merged.csv") or (base.startswith("eaglei_outages_with_events_") and "lag" not in base):
                yr = int(base.replace("eaglei_outages_with_events_", "").replace("eaglei_outages_", "")[:4])
                if YEARS[0] <= yr <= YEARS[1]:
                    frames.append(pd.read_csv(BytesIO(z.read(n)), low_memory=False,
                                              usecols=["fips", "state", "start_time", "duration", "min_customers",
                                                       "max_customers", "mean_customers"]))
    return pd.concat(frames, ignore_index=True)


def primary_obi(eag, mcc, scope_fips=None):
    """Paper method: dedupe, 2018-2023, COH = mean customers x hours, / MCC, 1.5xIQR clip, min-max."""
    e = eag.drop_duplicates(["fips", "start_time", "duration", "min_customers", "max_customers", "mean_customers"]).copy()
    e["start"] = pd.to_datetime(e["start_time"], format="mixed", errors="coerce")
    e = e[e["start"].dt.year.between(*YEARS)]
    e["coh"] = e["mean_customers"] * e["duration"]
    r = e.groupby("fips").agg(coh=("coh", "sum"), episodes=("coh", "size")).join(mcc, how="inner")
    if scope_fips is not None:
        r = r[r.index.isin(scope_fips)]
    r["R"] = r["coh"] / r["Customers"]
    q1, q3 = r["R"].quantile([0.25, 0.75])
    rc = r["R"].clip(q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1))
    r["obi"] = (rc - rc.min()) / (rc.max() - rc.min())
    return r, len(eag), len(e)


def main(nin, ca_static, ca_out, out):
    t0 = time.time()
    nin, ca_static, ca_out, out = map(Path, (nin, ca_static, ca_out, out))
    out.mkdir(parents=True, exist_ok=True)
    (out / "geo").mkdir(exist_ok=True)

    # ---------------- input audit ----------------
    manifest = []
    for f, meta in SOURCES.items():
        p = nin / f
        row = {"file": f, **{k: v for k, v in meta.items() if k != "sha256"}, "present": p.exists()}
        if p.exists():
            row["sha256"] = sha256(p)
            row["bytes"] = p.stat().st_size
            if meta.get("sha256"):
                chk(f"Fingerprint: {f}", row["sha256"] == meta["sha256"],
                    "matches the file downloaded from the publisher" if row["sha256"] == meta["sha256"]
                    else f"SHA-256 {row['sha256'][:16]}… differs from the recorded download")
        chk(f"File present: {f}", p.exists(), "found" if p.exists() else "not found")
        manifest.append(row)
    mcc_p = ca_static / "MCC.csv"
    manifest.append({"file": "MCC.csv", "name": "Modeled County Customers (MCC)", "publisher": "Oak Ridge National Laboratory (EAGLE-I companion)",
                     "year": "2022", "url": "https://doi.ccs.ornl.gov/dataset/ccec86f0-e144-5de8-aee0-fb26028b26e1",
                     "use": "Customer denominator for Primary OBI", "present": mcc_p.exists(), "sha256": sha256(mcc_p), "bytes": mcc_p.stat().st_size})
    if any(r["result"] == "FAIL" for r in AUDIT):
        return finish(out, manifest, {}, gate_only=True)

    # ---------------- load ----------------
    svi = pd.read_csv(nin / "SVI_2022_US.csv", dtype={"FIPS": str, "STCNTY": str}, encoding="utf-8-sig",
                      usecols=["ST", "ST_ABBR", "STATE", "STCNTY", "COUNTY", "FIPS", "RPL_THEMES", "E_TOTPOP"])
    svi["FIPS"] = svi["FIPS"].str.zfill(11)
    svic = pd.read_csv(nin / "SVI_2022_US_county.csv", dtype={"FIPS": str}, encoding="utf-8-sig",
                       usecols=["ST_ABBR", "COUNTY", "FIPS", "RPL_THEMES", "E_TOTPOP"])
    svic["FIPS"] = svic["FIPS"].str.zfill(5)
    eji = read_zip_csv(nin / "EJI_2024_United_States_CSV.zip", "EJI_2024_United_States.csv", dtype={"GEOID": str},
                       usecols=["GEOID", "STATEABBR", "RPL_EJI", "RPL_EBM", "E_TOTPOP"])
    eji["GEOID"] = eji["GEOID"].str.zfill(11)
    nri = read_zip_csv(nin / "NRI_Table_Counties_v120.zip", "NRI_Table_Counties.csv", dtype={"STCOFIPS": str},
                       usecols=["STCOFIPS", "STATEABBRV", "NRI_VER"] + list(HAZARDS), low_memory=False)
    nri["STCOFIPS"] = nri["STCOFIPS"].str.zfill(5)
    mccr = pd.read_csv(mcc_p)
    mccr = mccr[mccr["County_FIPS"].astype(str) != "Grand Total"]
    mcc = mccr.assign(fips=mccr["County_FIPS"].astype(int)).set_index("fips")[["Customers"]]
    eag = load_outages(nin / "OEDI6458_Outage_Dataset_R1.zip")

    tracts = pd.concat([gpd.read_file(f"zip://{nin / f'tiger2022/tl_2022_{s}_tract.zip'}") for s in STATES], ignore_index=True)
    tracts = gpd.GeoDataFrame(tracts[["GEOID", "STATEFP", "COUNTYFP", "ALAND", "geometry"]], crs=tracts.crs)
    n_water = int((tracts["ALAND"] == 0).sum())
    tracts = tracts[tracts["ALAND"] > 0].reset_index(drop=True)   # water-only tracts have no residents and no SVI/EJI row
    counties = gpd.read_file(f"zip://{nin / 'tiger2022/tl_2022_us_county.zip'}")
    counties = counties[counties["STATEFP"].isin(STATES)][["GEOID", "STATEFP", "NAME", "NAMELSAD", "geometry"]].copy()
    states = gpd.read_file(f"zip://{nin / 'tiger2022/tl_2022_us_state.zip'}")
    states = states[states["STATEFP"].isin(STATES)][["STATEFP", "STUSPS", "NAME", "geometry"]].copy()

    # ---------------- data checks ----------------
    chk("Tract universe (TIGER 2022, 50 states + DC)", 84000 < len(tracts) < 86000, f"{len(tracts):,} land tracts ({n_water} water-only tracts excluded)")
    chk("County universe (TIGER 2022, 50 states + DC)", len(counties) == 3144, f"{len(counties):,} counties and county equivalents")
    chk("SVI tract file is the national release", len(svi) > 84000 and svi["ST_ABBR"].nunique() == 51,
        f"{len(svi):,} tracts, {svi['ST_ABBR'].nunique()} states + DC")
    m = tracts["GEOID"].isin(svi["FIPS"]).mean()
    chk("TIGER tracts found in SVI 2022", m > 0.999, f"{100 * m:.2f}% of TIGER 2022 tracts")
    rp = svi["RPL_THEMES"]
    chk("SVI tract ranks in 0-1 or suppressed", (rp[rp >= 0] <= 1).all(), f"{int((rp < 0).sum()):,} tracts with no rank (CDC code -999), left missing")
    chk("SVI county file covers all counties", svic["FIPS"].isin(counties["GEOID"]).sum() == 3144, f"{svic['FIPS'].nunique():,} counties")
    e_us = eji[eji["GEOID"].str[:2].isin(STATES)]
    m = tracts["GEOID"].isin(e_us["GEOID"]).mean()
    chk("TIGER tracts found in EJI 2024", m > 0.99, f"{100 * m:.2f}% of TIGER 2022 tracts")
    akhi = e_us[e_us["STATEABBR"].isin(["AK", "HI"])]
    chk("EJI 2024 ranks for Alaska and Hawaii", (akhi["RPL_EJI"] < 0).all(),
        "EJI 2024 publishes no ranks for AK and HI (-999); national EJBI is left missing there, SVI views still cover them", level="WARNING")
    for c in ("RPL_EJI", "RPL_EBM"):
        v = e_us[c]
        chk(f"EJI {c} in 0-1 or missing", (v[v >= 0] <= 1).all(), f"{int((v < 0).sum()):,} tracts without a rank, left missing")
    chk("NRI covers all counties", nri["STCOFIPS"].isin(counties["GEOID"]).sum() >= 3140,
        f"{nri['STCOFIPS'].isin(counties['GEOID']).sum():,} of 3,144 counties; version {nri['NRI_VER'].iloc[0]}")
    hv = nri[list(HAZARDS)]
    chk("NRI hazard scores in 0-100", bool((hv.min().min() >= 0) and (hv.max().max() <= 100)),
        "all seven hazard risk scores in range; " + ", ".join(f"{HAZARDS[k].split(' (')[0]} missing for {int(v)} counties" for k, v in hv.isna().sum().items() if v)
        + " (Climate Stress left missing there)")
    chk("Outage episodes have positive durations and customers",
        (eag["duration"] > 0).all() and (eag["mean_customers"] > 0).all(), f"{len(eag):,} episode rows, 2018-2023 files")
    yrs = pd.to_datetime(eag["start_time"], format="mixed", errors="coerce").dt.year
    chk("Outage dates parse", yrs.notna().mean() > 0.999, f"{int(yrs.isna().sum())} rows with unreadable dates")

    # ---------------- tract national EJBI ----------------
    t = pd.DataFrame({"GEOID": tracts["GEOID"], "st": tracts["STATEFP"], "fips": tracts["STATEFP"] + tracts["COUNTYFP"]})
    t = t.merge(svi[["FIPS", "RPL_THEMES", "E_TOTPOP"]].rename(columns={"FIPS": "GEOID", "RPL_THEMES": "svi", "E_TOTPOP": "pop"}), on="GEOID", how="left")
    t = t.merge(e_us[["GEOID", "RPL_EJI", "RPL_EBM"]].rename(columns={"RPL_EJI": "eji", "RPL_EBM": "ebm"}), on="GEOID", how="left")
    for c in ("svi", "eji", "ebm"):
        t.loc[t[c] < 0, c] = np.nan
    t["ejbi"] = t[["svi", "eji", "ebm"]].mean(axis=1, skipna=False)

    # ---------------- county measures ----------------
    def wmean(g, col):
        g = g.dropna(subset=[col, "pop"])
        return np.average(g[col], weights=g["pop"]) if len(g) and g["pop"].sum() > 0 else np.nan

    c = counties.rename(columns={"GEOID": "fips5"}).drop(columns="geometry").copy()
    c["ejbi"] = c["fips5"].map(t.groupby("fips").apply(lambda g: wmean(g, "ejbi")))
    c = c.merge(svic[["FIPS", "RPL_THEMES", "E_TOTPOP"]].rename(columns={"FIPS": "fips5", "RPL_THEMES": "svi_c", "E_TOTPOP": "pop"}),
                on="fips5", how="left")
    c.loc[c["svi_c"] < 0, "svi_c"] = np.nan
    c["fips"] = c["fips5"].astype(int)

    # Primary OBI, national scale (Connecticut excluded: former counties do not nest in planning regions)
    eag_us = eag[~eag["fips"].between(9000, 9999)]
    obi, n_raw, n_used = primary_obi(eag_us, mcc, scope_fips=set(c["fips"]))
    c = c.merge(obi[["coh", "episodes", "Customers", "R", "obi"]].rename(columns={"R": "obi_raw", "Customers": "customers"}),
                left_on="fips", right_index=True, how="left")

    # Climate Stress, national min-max per hazard
    hz = nri.set_index("STCOFIPS")[list(HAZARDS)]
    hz = hz[hz.index.isin(c["fips5"])]
    norm = (hz - hz.min()) / (hz.max() - hz.min())
    c["cs"] = c["fips5"].map(norm.mean(axis=1, skipna=False))

    # ---------------- tract-in-county structure ----------------
    t = t.merge(c[["fips5", "ejbi", "svi_c"]].rename(columns={"fips5": "fips", "ejbi": "ejbi_c"}), on="fips", how="left")
    t["masked"] = (t["ejbi"] >= T_TRACT) & (t["ejbi_c"] < T_COUNTY)
    t["svi_hidden"] = (t["svi"] >= T_TRACT) & (t["svi_c"] < T_SVI_COUNTY)
    g = t.groupby("fips")
    c["tracts"] = c["fips5"].map(g.size())
    c["tracts_valid"] = c["fips5"].map(g["ejbi"].count())
    c["tracts_high"] = c["fips5"].map(g["ejbi"].apply(lambda s: int((s >= T_TRACT).sum())))
    c["pop_high"] = c["fips5"].map(t[t["ejbi"] >= T_TRACT].groupby("fips")["pop"].sum()).fillna(0)
    c["tract_min"] = c["fips5"].map(g["ejbi"].min())
    c["tract_max"] = c["fips5"].map(g["ejbi"].max())
    c["masked_tracts"] = c["fips5"].map(g["masked"].sum()).fillna(0).astype(int)
    c["masked_pop"] = c["fips5"].map(t[t["masked"]].groupby("fips")["pop"].sum()).fillna(0)
    c["svi_tracts_high"] = c["fips5"].map(g["svi"].apply(lambda s: int((s >= T_TRACT).sum())))
    c["svi_hidden_tracts"] = c["fips5"].map(g["svi_hidden"].sum()).fillna(0).astype(int)
    c["svi_hidden_pop"] = c["fips5"].map(t[t["svi_hidden"]].groupby("fips")["pop"].sum()).fillna(0)
    c["high_county"] = c["ejbi"] >= T_COUNTY

    # ---------------- MAUP by state and nationally ----------------
    st_names = states.set_index("STATEFP")[["STUSPS", "NAME"]]
    rows = []
    for s in STATES:
        ts, cs_ = t[t["st"] == s], c[c["STATEFP"] == s]
        vt, vc = ts["ejbi"].var(ddof=1), cs_["ejbi"].var(ddof=1)
        rows.append({
            "statefp": s, "abbr": st_names.loc[s, "STUSPS"], "name": st_names.loc[s, "NAME"],
            "tracts": int(len(ts)), "tracts_valid": int(ts["ejbi"].notna().sum()), "counties": int(len(cs_)),
            "tracts_high": int((ts["ejbi"] >= T_TRACT).sum()), "counties_high": int((cs_["ejbi"] >= T_COUNTY).sum()),
            "var_tract": r6(vt), "var_county": r6(vc),
            "suppressed_pct": r6(100 * (1 - vc / vt)) if len(cs_) > 1 and vt and np.isfinite(vt) and np.isfinite(vc) else None,
            "masked_tracts": int(ts["masked"].sum()), "masked_pop": int(ts.loc[ts["masked"], "pop"].sum()),
            "masked_counties": int(cs_["masked_tracts"].gt(0).sum()),
            "svi_hidden_tracts": int(ts["svi_hidden"].sum()), "svi_hidden_pop": int(ts.loc[ts["svi_hidden"], "pop"].sum()),
            "pop": int(ts["pop"].sum()),
            "obi_counties": int(cs_["obi"].notna().sum()),
            "r_ejbi_obi": (corr(cs_, "ejbi", "obi", ci=False) or {}).get("pearson") if cs_["obi"].notna().sum() >= 10 else None,
        })
    st_df = pd.DataFrame(rows)
    vt, vc = t["ejbi"].var(ddof=1), c["ejbi"].var(ddof=1)
    nat = {
        "tracts": int(len(t)), "tracts_valid": int(t["ejbi"].notna().sum()), "counties": int(len(c)),
        "tracts_high": int((t["ejbi"] >= T_TRACT).sum()), "counties_high": int(c["high_county"].sum()),
        "var_tract": r6(vt), "var_county": r6(vc), "suppressed_pct": r6(100 * (1 - vc / vt)),
        "masked_tracts": int(t["masked"].sum()), "masked_pop": int(t.loc[t["masked"], "pop"].sum()),
        "masked_counties": int(c["masked_tracts"].gt(0).sum()),
        "svi_hidden_tracts": int(t["svi_hidden"].sum()), "svi_hidden_pop": int(t.loc[t["svi_hidden"], "pop"].sum()),
        "svi_hidden_counties": int(c["svi_hidden_tracts"].gt(0).sum()),
        "pop": int(t["pop"].sum()),
        "within_state_suppressed_median": r6(st_df["suppressed_pct"].median()),
        "states_with_masked": int((st_df["masked_tracts"] > 0).sum()),
        "obi_counties": int(c["obi"].notna().sum()), "outage_rows": int(n_raw), "outage_rows_used": int(n_used),
        "obi_mean": r6(c["obi"].mean()), "obi_sd": r6(c["obi"].std()),
        "cs_counties": int(c["cs"].notna().sum()), "cs_mean": r6(c["cs"].mean()),
    }

    # ---------------- relations ----------------
    rel = {"ejbi_obi": corr(c, "ejbi", "obi"), "ejbi_cs": corr(c, "ejbi", "cs"), "cs_obi": corr(c, "cs", "obi"),
           "svi_obi": corr(c, "svi_c", "obi")}
    hi, lo = c.loc[c["high_county"], "obi"].dropna(), c.loc[~c["high_county"] & c["ejbi"].notna(), "obi"].dropna()
    w = stats.ttest_ind(hi, lo, equal_var=False)
    rel["welch"] = {"n_high": int(len(hi)), "n_other": int(len(lo)), "mean_high": r6(hi.mean()), "mean_other": r6(lo.mean()),
                    "t": r6(w.statistic), "p": float(w.pvalue)}

    # ---------------- spatial (contiguous states; island counties dropped by the weights) ----------------
    from esda.moran import Moran, Moran_Local
    from libpysal.weights import Queen
    cg = counties.rename(columns={"GEOID": "fips5"}).merge(c[["fips5", "ejbi", "obi"]], on="fips5")
    cg = cg[~cg["STATEFP"].isin(["02", "15"])].to_crs(5070)
    spatial = {}
    for col in ("ejbi", "obi"):
        d = cg.dropna(subset=[col]).reset_index(drop=True)
        wq = Queen.from_dataframe(d, use_index=False, silence_warnings=True)
        keep = [i for i in range(len(d)) if wq.cardinalities[i] > 0]
        d = d.iloc[keep].reset_index(drop=True)
        wq = Queen.from_dataframe(d, use_index=False, silence_warnings=True)
        wq.transform = "r"
        np.random.seed(SEED)
        mi = Moran(d[col].values, wq, permutations=9999)
        np.random.seed(SEED)
        li = Moran_Local(d[col].values, wq, permutations=9999, seed=SEED)
        lab = np.array(["ns"] * len(d), dtype=object)
        sig = li.p_sim < 0.05
        for q, name in ((1, "HH"), (2, "LH"), (3, "LL"), (4, "HL")):
            lab[sig & (li.q == q)] = name
        c[f"lisa_{col}"] = c["fips5"].map(dict(zip(d["fips5"], lab)))
        spatial[col] = {"n": int(len(d)), "I": r6(mi.I), "p": float(mi.p_sim),
                        "counts": {k: int((lab == k).sum()) for k in ("HH", "LL", "LH", "HL")}}

    # ---------------- California concordance with the validated paper run ----------------
    conc = {}
    paper = pd.read_csv(ca_out / "segda_counties.csv")
    ca = c[c["STATEFP"] == "06"].copy()
    ca_obi, _, _ = primary_obi(eag_us, mcc, scope_fips=set(ca["fips"]))   # paper scale: normalized within California
    mm = paper[["fips", "obi", "ejbi", "svi22", "cs"]].merge(ca_obi[["obi"]].rename(columns={"obi": "obi_pub"}),
                                                            left_on="fips", right_index=True, how="inner")
    mm = mm.merge(ca[["fips", "ejbi", "cs"]].rename(columns={"ejbi": "ejbi_nat", "cs": "cs_nat"}), on="fips")
    d = mm.dropna(subset=["obi", "obi_pub"])
    conc["obi_r"] = r6(stats.pearsonr(d["obi"], d["obi_pub"])[0])
    conc["obi_n"] = int(len(d))
    chk("California Primary OBI: public release vs the paper's file", conc["obi_r"] >= 0.95,
        f"Pearson r = {conc['obi_r']:.3f} across {len(d)} counties, same method, California-only scaling (threshold 0.95)")
    d = mm.dropna(subset=["ejbi", "ejbi_nat"])
    conc["ejbi_r"] = r6(stats.pearsonr(d["ejbi"], d["ejbi_nat"])[0])
    conc["ejbi_n"] = int(len(d))
    chk("California county EJBI: national version vs the paper's CalEnviroScreen version", conc["ejbi_r"] >= 0.80,
        f"Pearson r = {conc['ejbi_r']:.3f} across {len(d)} counties (threshold 0.80; different components and tract vintages)")
    d = mm.dropna(subset=["cs", "cs_nat"])
    conc["cs_r"] = r6(stats.pearsonr(d["cs"], d["cs_nat"])[0])
    chk("California Climate Stress: NRI v1.20 national scaling vs the paper", conc["cs_r"] >= 0.80,
        f"Pearson r = {conc['cs_r']:.3f} across {len(d)} counties (threshold 0.80; newer NRI version, national scaling)")
    ca_row = st_df[st_df["statefp"] == "06"].iloc[0]
    conc["ca_suppressed_national_method"] = ca_row["suppressed_pct"]
    conc["ca_suppressed_paper"] = 67.82
    chk("California variance suppression survives the change of index", ca_row["suppressed_pct"] > 50,
        f"{ca_row['suppressed_pct']:.2f}% with the national EJBI vs 67.82% in the paper", level="WARNING")

    # ---------------- outputs ----------------
    t_out = t[["GEOID", "st", "fips", "svi", "eji", "ebm", "ejbi", "pop", "masked", "svi_hidden"]]
    t_out.to_csv(out / "national_tracts.csv", index=False, float_format="%.6g")
    cols = ["fips5", "STATEFP", "NAME", "NAMELSAD", "pop", "ejbi", "svi_c", "obi", "obi_raw", "coh", "episodes", "customers", "cs",
            "tracts", "tracts_valid", "tracts_high", "pop_high", "tract_min", "tract_max", "masked_tracts", "masked_pop",
            "svi_tracts_high", "svi_hidden_tracts", "svi_hidden_pop", "lisa_ejbi", "lisa_obi"]
    c[cols].to_csv(out / "national_counties.csv", index=False, float_format="%.6g")
    st_df.to_csv(out / "national_states.csv", index=False)

    # geometry for the site: WGS84 lon/lat (the browser projects with Albers USA)
    tg = tracts[["GEOID", "STATEFP", "geometry"]].to_crs(4326)
    for s in STATES:
        tg[tg["STATEFP"] == s][["GEOID", "geometry"]].to_file(out / "geo" / f"tracts_{s}.geojson", driver="GeoJSON")
    counties.to_crs(4326)[["GEOID", "STATEFP", "geometry"]].to_file(out / "geo" / "us_counties.geojson", driver="GeoJSON")
    states.to_crs(4326)[["STATEFP", "STUSPS", "NAME", "geometry"]].to_file(out / "geo" / "us_states.geojson", driver="GeoJSON")
    return finish(out, manifest, {"national": nat, "relations": rel, "spatial": spatial, "concordance": conc,
                                  "runtime_s": round(time.time() - t0, 1)})


def finish(out, manifest, res, gate_only=False):
    n_fail = sum(r["result"] == "FAIL" for r in AUDIT)
    n_warn = sum(r["result"] == "WARNING" for r in AUDIT)
    n_pass = sum(r["result"] == "PASS" for r in AUDIT)
    res["run"] = {"gate": "PASS" if n_fail == 0 and not gate_only else "FAIL", "audit_pass": n_pass, "audit_warn": n_warn,
                  "audit_fail": n_fail, "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    res["audit"] = AUDIT
    res["sources"] = manifest
    res["thresholds"] = {"tract": T_TRACT, "county": T_COUNTY, "svi_county": T_SVI_COUNTY}
    res["hazards"] = HAZARDS
    pd.DataFrame(AUDIT).to_csv(out / "national_audit.csv", index=False)
    pd.DataFrame(manifest).to_csv(out / "national_sources.csv", index=False)
    (out / "national_results.json").write_text(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    for r in AUDIT:
        if r["result"] != "PASS":
            print(f"  {r['result']:8s} {r['check']}: {r['detail']}")
    print(f"national gate {res['run']['gate']}: audit {n_pass} pass, {n_warn} warning, {n_fail} fail")
    return 0 if res["run"]["gate"] == "PASS" else 2


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:5]))
