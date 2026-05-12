# oilftir

**FTIR-based condition monitoring of lubricating oils — ASTM E2412**

`oilftir` is a Python library for processing and analysing FTIR spectra of
used lubricating oils following the
[ASTM E2412](https://www.astm.org/e2412-23.html) standard.  
It provides spectrum loading, baseline correction, net area integration for
all standard degradation and contamination bands, concentration estimation
from calibration curves, spectral library search (HQI), and publication-ready
plots.

---

## Features

| Feature | Description |
|---|---|
| **Spectrum loading** | CSV / TXT and PerkinElmer `.sp` (optional) |
| **Baseline correction** | Modified polynomial (ModPoly) via `pybaselines` |
| **ASTM E2412 areas** | Oxidation, nitration, sulfation, water, glycol, petrol, diesel, ZDDP, soot |
| **Calibration** | Generic polynomial curve — bring your own standards |
| **Reference subtraction** | Least-squares scaling + difference spectrum |
| **Library search** | Hit Quality Index (HQI), parallel brute-force |
| **Plotting** | Dual-panel (absorbance + transmittance), multi-spectrum overlay |

---

## Installation

```bash
pip install oilftir
```

To read PerkinElmer `.sp` files, install the optional dependency:

```bash
pip install oilftir[specio]
```

### From source

```bash
git clone https://github.com/your-username/oilftir.git
cd oilftir
pip install -e .
```

---

## Quick start

```python
from oilftir.astm_utils import load_spectrum, remove_baseline, astm_areas

# Load and correct a spectrum
df, meta = load_spectrum("used_oil.csv", magnitude="A")
df = remove_baseline(df)

# Compute all ASTM E2412 band areas
areas = astm_areas(df)
print(areas)
# {
#   'oxidation':     0.2341,
#   'nitration':     0.0812,
#   'sulfation':     0.1503,
#   'water':         4.7620,
#   'glycol':        0.0000,
#   'fuel_petrol':   0.0034,
#   'fuel_diesel':   0.0198,
#   'antiwear_zddp': 8.3410,
#   'soot':         12.5000,
# }
```

See [`example.py`](example.py) for a full walkthrough including concentration
estimation, reference subtraction, plotting, and library search.

---

## Input file format

Plain-text files must be two-column (wavenumber, amplitude), comma- or
space-delimited, with `#` for comment lines:

```
# wavenumber (cm-1), absorbance
4000.0  0.0031
3999.0  0.0028
...
```

---

## API reference

### Loading

```python
load_spectrum(path, magnitude='A') → (DataFrame, dict)
```
`magnitude`: `'A'` (absorbance) or `'T'` (transmittance, 0–1 or 0–100).

### Pre-processing

```python
remove_baseline(df, units='A', poly_order=5) → DataFrame
scale_reference(df_sample, df_ref, units='A') → DataFrame
subtract_reference(df_sample, df_ref, wavenumber_range, units='A')
    → (y_diff, k, y_ref_scaled)
```

### ASTM E2412 band areas

```python
astm_areas(df)          → dict   # all bands at once
area_oxidation(df)      → float  # 1670–1800 cm⁻¹
area_nitration(df)      → float  # 1600–1650 cm⁻¹
area_sulfation(df)      → float  # 1120–1180 cm⁻¹
area_water(df)          → float  # 3150–3500 cm⁻¹
area_glycol(df)         → float  # 1030–1100 cm⁻¹
area_fuel_petrol(df)    → float  # 745–755 cm⁻¹
area_fuel_diesel(df)    → float  # 805–815 cm⁻¹
area_antiwear_zddp(df)  → float  # 960–1025 cm⁻¹
soot_load(df)           → float  # absorbance at 2000 cm⁻¹ × 100
```

All area functions return a net integrated area in absorbance·cm⁻¹ (≥ 0).

### Calibration

```python
area_to_concentration(area, calibration_areas, calibration_concentrations,
                       poly_order=3) → float
```

Provide your own (area, concentration) calibration pairs — the function fits
a polynomial and evaluates it at `area`.

### Library search

```python
library_search(query_path, library_dir, top_n=5, extension='*.sp')
    → list[dict]   # sorted by HQI descending
```

```python
hqi(x_query, y_query, x_ref, y_ref) → float   # 0–100 %
```

### Plotting

```python
plot_spectra(data_dict, title, xlim, ylim_A, ylim_T, save_path)
    → (fig, axes)
```

`data_dict` format:
```python
{
    "Label": {
        "df":        <DataFrame>,   # required
        "color":     "steelblue",   # optional
        "alpha":     1.0,           # optional
        "linewidth": 1.5,           # optional
    }
}
```

---

## Dependencies

| Package | Purpose |
|---|---|
| `numpy` | Numerical operations |
| `pandas` | Spectrum data frames |
| `scipy` | Spearman correlation (HQI fallback) |
| `matplotlib` | Plotting |
| `pybaselines` | Baseline correction |
| `specio` *(optional)* | PerkinElmer `.sp` file reading |

---

## Contributing

Contributions are welcome. Please open an issue or pull request on GitHub.

---
## Authors
`oilftir` was developed at the Chemistry Laboratory of CENADIF ([Centro Nacional de Desarrollo e Innovación Ferroviaria](https://www.argentina.gob.ar/transporte/fase/cenadif)), Argentina, in the context of used-oil condition monitoring for the local railway industry.

### Lead developer

**Nahuel Mendez** — R&D Engineer/Physical Data Scientist, CENADIF |

[![GitHub](https://img.shields.io/badge/GitHub-NahueMendez-181717?logo=github)](https://github.com/NahueMendez)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-ingnahuelmendez-0A66C2?logo=linkedin)](https://www.linkedin.com/in/ingnahuelmendez/)
[![Email](https://img.shields.io/badge/Email-nahueldanielmendez%40gmail.com-D14836?logo=gmail)](mailto:nahueldanielmendez@gmail.com)


### Collaborators

| Name | Role |
| :--- | :--- |
| **Hernán Gomez Molino**, Eng. | Head of Laboratories, CENADIF |
| **Leandro Asens**, Eng. | Head of the Chemistry Laboratory, CENADIF |
---

## Citing

If you use `oilftir` in your research or academic work, please cite both the ASTM standard and this software library:

### 1. The Software Library
> Mendez, N., Gomez Molino, H., & Asens, L. (2026). *oilftir: A Python library for FTIR-based condition monitoring of lubricating oils following ASTM E2412*. https://doi.org/10.5281/zenodo.20137146

### 2. The Implemented Standard
> ASTM E2412-23, *Standard Practice for Condition Monitoring of Used Lubricants by Trend Analysis Using Fourier Transform Infrared (FT-IR) Spectrometry*, ASTM International, West Conshohocken, PA.

You can also use the following BibTeX entry for your reference manager:

```bibtex
@software{mendez_oilftir_2026,
  author       = {Mendez, Nahuel and Gomez Molino, Hern{\'a}n and Asens, Leandro},
  title        = {oilftir: A Python library for FTIR-based condition monitoring of lubricating oils following ASTM E2412},
  month        = may,
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.20137146},
  url          = {doi.org}
}
```
---

## License

```
MIT License

Copyright (c) 2026 oilftir contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
