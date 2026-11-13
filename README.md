# GridAI™ · SEGDA

**Live site: [gridaitm.com](https://gridaitm.com/)**

The Socio-Environmental Grid Disruption Analysis (SEGDA) examines environmental justice, county-level power-outage burden, and climate stress across California's census tracts and counties. It tests how geographic aggregation changes measured burden and considers implications for equitable resilience and public health preparedness.

**Developer and Author:** Victoria Love Franklin, Franklin HUB SMART Center.

**Keywords:** GridAI, SEGDA, Socio-Environmental Grid Disruption Analysis, energy justice, environmental justice, power outage burden, grid resilience, climate resilience, public health equity, census tract analysis, modifiable areal unit problem (MAUP), CDC Social Vulnerability Index, Environmental Justice Index, EAGLE-I, FEMA National Risk Index, geospatial analytics, GIS.


**Socio-Environmental Grid Disruption Analysis: how county averages hide census-tract burden, across the United States.**

The site has two parts. The California case study came first, and it led to the national study:

- **National study** (all 50 states and DC, 84,095 census tracts, 3,144 counties). The paper's methods are applied with national federal data. County averaging erases 72.7% of tract-level variation in the national burden index. It hides 7,803 high-burden tracts, home to 28.3 million people, inside 890 counties that do not look high.
- **California case study.** The original SEGDA analysis, reproduced exactly from the paper's own files (31/31 locked values). Its finding that county averages hide high-burden tracts prompted the national study.

An interactive site built from the SEGDA paper (Franklin). It shows the same environmental-justice burden index at two scales, 7,931 California census tracts and 58 counties. Switching between them shows high-burden neighbourhoods disappearing into moderate county colours: the modifiable areal unit problem (MAUP).

The site is rebuilt by GitHub Actions on every push. Before anything is published:

1. **Input audit.** Every input file is checked: present, expected columns, value ranges, coverage. Each file is fingerprinted with SHA-256.
2. **Validation gate.** 31 values reported in the paper are recomputed from the raw files and must match within tolerance.

If either check fails, the build stops and the previous site stays up.

No value on the site is simulated, estimated to fill a gap, or typed in by hand. Missing values stay missing and are shown as unavailable.

## What it shows

| View | Content |
|---|---|
| What counties hide | Tract vs county maps (EJBI or CDC SVI), hidden high-burden tracts outlined, per-county distribution, variance suppression |
| Outages & climate | Primary Outage Burden Index (EAGLE-I / MCC, 2018–2023) and Climate Stress (FEMA NRI) |
| How they relate | Pearson/Spearman with 5,000-resample bootstrap CIs; Welch's t-test |
| Clusters | Global Moran's I and LISA (queen contiguity, 9,999 permutations) |
| Robustness | Alternative outage specifications, CDC SVI concordance, EIA-861 sensitivity |
| Methods & validation | Locked-value table, formulas, data vintages, downloads, limitations |
| Data & audit | Audit results, source manifest with SHA-256, references |

## Key results (reproduced from the inputs)

| Measure | Value |
|---|---|
| Tract-level variance removed by county averaging | 67.82% |
| High-burden tracts (EJBI ≥ 0.75) inside counties below 0.60 | 431 tracts, 2,340,612 people, 21 counties |
| Top-quartile SVI tracts in counties CDC's county rank does not flag | 728 tracts, 3,800,615 people |

## Methods

The **EJBI** is the mean of three components:

- SVI `RPL_THEMES`
- CalEnviroScreen `CIscoreP`/100
- `Pollution_Pctl`/100

It is computed only when all three components exist. County EJBI is population-weighted.

Thresholds: tract ≥ 0.75, county ≥ 0.60.

MAUP suppression = 1 − Var(county)/Var(tract), using sample variance.

**Primary OBI** is built in five steps:

1. Deduplicate EAGLE-I episodes and keep 2018–2023.
2. Sum customer-outage-hours (mean customers × duration) per county.
3. Divide by Modeled County Customers.
4. Clip at 1.5 × IQR.
5. Min-max scale.

**Climate Stress** is the equal-weight mean of seven min-max scaled NRI hazard scores.

Full detail is on the site's Methods page and in `pipeline/segda_pipeline.py`.

## Data (in `static_inputs/`)

| File | Source |
|---|---|
| CalEnviroScreen_4.0.csv | CalEPA / OEHHA, CalEnviroScreen 4.0 |
| SVI_CA_TRACT.csv | CDC/ATSDR SVI, California tracts (**2018 release, 2010 tracts**) |
| SVI_CA_COUNTY_2022.csv | CDC/ATSDR SVI 2022, California counties |
| eaglei_transformed.csv | DOE CESER / ORNL EAGLE-I historic outages, California |
| MCC.csv | ORNL Modeled County Customers |
| NRI_Counties_CA.csv | FEMA National Risk Index, county table |
| EIA861_CA_COUNTY_RELIABILITY_6YR.csv | EIA Form 861 reliability allocated to counties (sensitivity only) |
| shp/tl_2019_06_tract.* | Census TIGER/Line 2019, California tracts |
| shp/tl_2019_06_county_from_us.* | California rows (STATEFP 06) of Census TIGER/Line 2019 `tl_2019_us_county`; geometry unchanged |

The national county shapefile is 127 MB, over GitHub's 100 MB file limit, so only its 58 California rows are committed. The original national file's SHA-256 is recorded in the source manifest. The pipeline also accepts the full national file at `shp/tl_2019_us_county.shp`.

## Run locally

```bash
pip install -r requirements.txt
python pipeline/segda_pipeline.py static_inputs out   # audit + gate; exits 2 on failure
bash scripts/make_topo.sh out                          # needs Node (uses npx mapshaper)
python scripts/build_site.py out site                  # refuses unless gate = PASS
# open site/index.html
```

## National study

See the site's Methods page for how often each source is reported and updated. EAGLE-I records outages every 15 minutes and publishes yearly files; SVI is released every two years; NRI and EJI periodically.

`pipeline/segda_national.py` applies the paper's methods nationally. Every change from the California study is forced by the available data and is stated on the site:

| Element | California study | National layer |
|---|---|---|
| Burden index | mean(SVI, CalEnviroScreen %ile, Pollution Burden %ile) | mean(CDC SVI 2022 `RPL_THEMES`, CDC/ATSDR EJI 2024 `RPL_EJI`, EJI `RPL_EBM`) |
| Tracts | 2010 (TIGER 2019) | 2020 (TIGER 2022), land tracts only |
| Outages | Paper's EAGLE-I file | Public OEDI release of the same EAGLE-I merged episodes (PNNL submission 6458, v2), with the same method |
| Climate Stress | NRI study vintage, scaled within CA | NRI v1.20, scaled nationally; Inland Flooding replaces Riverine Flooding |

The national inputs total about 655 MB, so they are not committed. They are official files downloaded from their publishers, and their SHA-256 fingerprints are recorded in `national_inputs.sha256` and checked by the pipeline. The committed `national_build/` folder holds the pipeline's outputs and audit, and the site build refuses to run unless its gate is PASS. To regenerate:

```bash
python pipeline/segda_national.py national_inputs static_inputs out national_out
bash scripts/make_national_topo.sh national_out
# copy national_out/*.csv, *.json and topo/ into national_build/
```

**Gate checks against the paper (California counties):**

| Check | Result |
|---|---|
| Primary OBI, public release vs paper file | r = 0.987 |
| County EJBI, national vs CalEnviroScreen version | r = 0.807 |
| Climate Stress | r = 0.997 |
| Variance suppression | 72.27% (paper 67.82%) |

**Known gaps, left missing rather than filled:**

- EJI 2024 has no ranks for Alaska and Hawaii.
- Connecticut planning regions do not nest the former counties that EAGLE-I uses, so Connecticut has no OBI.
- 30 counties lack an NRI heat-wave score.

**National finding:** the EJBI–OBI association is near zero nationally (r = 0.054), and its sign differs by state. California's negative association (r = −0.44) is not a national pattern.

## Known notes

- **SVI vintage.** The tract SVI file supplied is the 2018 release (2010 tract geography). The paper's Table 1 lists SVI 2022. All locked results reproduce from the 2018 file, so the file is right and the table label may need an erratum. The site labels the vintage.
- **LISA cluster counts are permutation-sensitive.** Global Moran's I values match the paper exactly (0.2017, 0.4010). Local cluster counts vary between permutation runs. The site lists the paper's counts alongside this run's and does not claim they match.
- **EIA-861** has 57 counties (no Trinity). It enters only the sensitivity branch, so this is an audit warning, not a failure.

## Roadmap

- **Phase 2:** all U.S. states. CDC SVI 2022 tracts, TIGER tracts as vector tiles, and county outage burden from EAGLE-I + MCC, all downloaded inside GitHub Actions.
- **Phase 3:** live county outages from official open feeds only (e.g. ORNL ODIN, state emergency-management APIs). Every value carries its source and timestamp. No scraping of utility websites.

## Citation

Franklin, V. L. SEGDA: Socio-Environmental Grid Disruption Analysis. Franklin HUB SMART Center. https://gridaitm.com/

© 2026 Victoria Love Franklin. All rights reserved. Source datasets remain under their publishers' terms. Most are U.S. federal or State of California public data.
