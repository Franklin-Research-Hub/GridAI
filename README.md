# SEGDA SMART Grid

**Socio-Environmental Grid Disruption Analysis: how county averages hide census-tract burden.**

An interactive site built from the SEGDA paper (Franklin, Meharry Medical College). It shows the same environmental-justice burden index at two scales, 7,931 California census tracts and 58 counties. Switching between them shows high-burden neighbourhoods disappearing into moderate county colours: the modifiable areal unit problem (MAUP).

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

## Known notes

- **SVI vintage.** The tract SVI file supplied is the 2018 release (2010 tract geography). The paper's Table 1 lists SVI 2022. All locked results reproduce from the 2018 file, so the file is right and the table label may need an erratum. The site labels the vintage.
- **LISA cluster counts are permutation-sensitive.** Global Moran's I values match the paper exactly (0.2017, 0.4010). Local cluster counts vary between permutation runs. The site lists the paper's counts alongside this run's and does not claim they match.
- **EIA-861** has 57 counties (no Trinity). It enters only the sensitivity branch, so this is an audit warning, not a failure.

## Roadmap

- **Phase 2:** all U.S. states. CDC SVI 2022 tracts, TIGER tracts as vector tiles, and county outage burden from EAGLE-I + MCC, all downloaded inside GitHub Actions.
- **Phase 3:** live county outages from official open feeds only (e.g. ORNL ODIN, state emergency-management APIs). Every value carries its source and timestamp. No scraping of utility websites.

## Citation

Franklin, V. L. SEGDA: Socio-Environmental Grid Disruption Analysis. Meharry Medical College, Franklin HUB SMART Center.

© 2026 Victoria Love Franklin. All rights reserved. Source datasets remain under their publishers' terms. Most are U.S. federal or State of California public data.
