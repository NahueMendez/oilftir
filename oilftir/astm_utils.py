# -*- coding: utf-8 -*-
"""
astm_utils.py
==================
FTIR Spectrum Processing Utilities for Lubricating Oil Analysis
Based on ASTM E2412 standard methodology.

Provides:
  - Spectrum loading (CSV/TXT; optional PerkinElmer .sp via specio)
  - Baseline correction and spectral normalisation
  - Net area calculation for ASTM-defined absorption bands
  - Hit Quality Index (HQI) library search
  - General-purpose plotting

Dependencies
------------
  numpy, pandas, scipy, matplotlib, pybaselines
  specio (optional – only needed for .sp files)

Usage example
-------------
  from oilftir.astm_utils import load_spectrum, astm_areas

  df, meta = load_spectrum("sample.csv")
  results = astm_areas(df)
  print(results)

References
----------
  ASTM E2412-23  Standard Practice for Condition Monitoring of Used Lubricants
                 by Trend Analysis Using Fourier Transform Infrared (FT-IR)
                 Spectrometry
"""

from __future__ import annotations
import re
import struct
import time
import concurrent.futures
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pybaselines import Baseline

# ---------------------------------------------------------------------------
# Optional dependency: specio (PerkinElmer .sp files)
# ---------------------------------------------------------------------------
try:
    from specio import specread
    import specio.plugins.sp as _sp_plugin

    def _decode_5104_smart(data):
        """
        Replacement decoder for specio's _decode_5104.
        Searches the binary blob for text segments rather than relying on
        fixed offsets, which are instrument-firmware–dependent.
        """
        text_items = []
        idx = 0
        while idx < len(data) - 2:
            try:
                possible_size = struct.unpack('<h', data[idx:idx + 2])[0]
            except Exception:
                break
            if 0 < possible_size < 500:
                start_str = idx + 2
                end_str = start_str + possible_size
                if end_str <= len(data):
                    chunk = data[start_str:end_str]
                    try:
                        decoded = chunk.decode('latin1')
                        if any(c.isalnum() for c in decoded):
                            text_items.append(decoded)
                            idx = end_str
                            next_tag = data.find(b'#u', idx)
                            if next_tag != -1:
                                idx = next_tag + 2
                                continue
                    except Exception:
                        pass
            idx += 1

        keys = ['description', 'axis_x', 'axis_y', 'date', 'time', 'sample_info']
        return {k: v for k, v in zip(keys, text_items) if v}

    _sp_plugin._decode_5104 = _decode_5104_smart
    _SPECIO_AVAILABLE = True

except ImportError:
    _SPECIO_AVAILABLE = False


# ===========================================================================
# DATA LOADING
# ===========================================================================

def _clean_binary_metadata(meta_dict: dict) -> dict:
    """
    Remove binary artefacts from a metadata dictionary, keeping only
    human-readable strings.

    Parameters
    ----------
    meta_dict : dict
        Raw metadata as returned by specio.

    Returns
    -------
    dict
        Cleaned metadata; non-string values are kept unchanged.
    """
    pattern = r'[a-zA-Z0-9\s:,\.\(\)-]{3,}'
    cleaned = {}
    for key, val in meta_dict.items():
        if not isinstance(val, str):
            cleaned[key] = val
            continue
        matches = re.findall(pattern, val)
        if not matches:
            cleaned[key] = "N/A" if len(val) > 10 else val
        else:
            text = " | ".join(m.strip() for m in matches)
            text = re.sub(r'\s+', ' ', text).strip()
            cleaned[key] = text
    return cleaned


def load_spectrum(path: str, magnitude: str = 'A') -> tuple[pd.DataFrame, dict]:
    """
    Load an FTIR spectrum from a file.

    Supported formats
    -----------------
    * CSV / TXT  – two-column file (wavenumber, amplitude); delimiter is
                   auto-detected (comma or whitespace).
    * .sp        – PerkinElmer binary format (requires the *specio* package).

    Parameters
    ----------
    path : str
        Path to the spectrum file.
    magnitude : {'A', 'T'}
        Indicates whether amplitudes in the file are Absorbance ('A') or
        Transmittance ('T', expressed as a fraction 0–1 or percentage 0–100).

    Returns
    -------
    df : pd.DataFrame
        Columns: ``cm-1`` (wavenumber), ``A`` (absorbance), ``T``
        (transmittance 0–1).
    meta : dict
        Metadata extracted from the file (empty dict for plain text files).

    Raises
    ------
    ValueError
        If *magnitude* is not 'A' or 'T'.
    ImportError
        If a .sp file is requested but *specio* is not installed.
    """
    if magnitude not in ('A', 'T'):
        raise ValueError("'magnitude' must be 'A' (absorbance) or 'T' (transmittance).")

    ext = Path(path).suffix.lower()

    # --- PerkinElmer .sp ---
    if ext == '.sp':
        if not _SPECIO_AVAILABLE:
            raise ImportError(
                "The 'specio' package is required to read .sp files. "
                "Install it with:  pip install specio"
            )
        spectra = specread(path)
        meta = _clean_binary_metadata(spectra.meta)
        df = pd.DataFrame({'cm-1': spectra.wavelength})
        if magnitude == 'A':
            df['A'] = spectra.amplitudes
            df['T'] = 10 ** (-df['A'].values)
        else:
            raw = np.array(spectra.amplitudes)
            # Accept both fractional (0-1) and percentage (0-100) transmittance
            df['T'] = raw / 100.0 if raw.max() > 1.5 else raw
            df['A'] = -np.log10(np.where(df['T'] > 0, df['T'], np.nan))
        return df, meta

    # --- Plain text (CSV / space-delimited) ---
    with open(path) as fh:
        first_line = fh.readline()
    delimiter = ',' if ',' in first_line else None
    data = np.loadtxt(path, delimiter=delimiter, comments='#')
    meta = {}
    df = pd.DataFrame({'cm-1': data[:, 0]})
    if magnitude == 'A':
        df['A'] = data[:, 1]
        df['T'] = 10 ** (-df['A'].values)
    else:
        raw = data[:, 1]
        df['T'] = raw / 100.0 if raw.max() > 1.5 else raw
        df['A'] = -np.log10(np.where(df['T'] > 0, df['T'], np.nan))
    return df, meta


# ===========================================================================
# PRE-PROCESSING
# ===========================================================================

def remove_baseline(df: pd.DataFrame, units: str = 'A',
                    poly_order: int = 5) -> pd.DataFrame:
    """
    Subtract a polynomial baseline from the spectrum using the modified
    polynomial method (ModPoly).

    Parameters
    ----------
    df : pd.DataFrame
        Spectrum DataFrame with columns ``cm-1`` and *units*.
    units : str
        Column to correct ('A' or 'T').
    poly_order : int
        Degree of the fitting polynomial (default 5).

    Returns
    -------
    pd.DataFrame
        Copy of *df* with the baseline subtracted in-place.
    """
    baseline_fitter = Baseline(x_data=df['cm-1'])
    baseline, _ = baseline_fitter.modpoly(
        df[units], poly_order=poly_order, tol=1e-3, max_iter=150
    )
    df = df.copy()
    df[units] = df[units] - baseline
    return df


def scale_reference(df_sample: pd.DataFrame, df_ref: pd.DataFrame,
                    units: str = 'A') -> pd.DataFrame:
    """
    Scale a reference spectrum to best match a sample spectrum using
    least-squares projection (scalar factor *k*).

    The scaling factor is computed as:

        k = (sample · reference) / (reference · reference)

    Parameters
    ----------
    df_sample : pd.DataFrame
        Sample spectrum.
    df_ref : pd.DataFrame
        Reference spectrum to be scaled.
    units : str
        Column used for comparison.

    Returns
    -------
    pd.DataFrame
        Scaled copy of *df_ref*.
    """
    vec_sample = df_sample[units].values
    vec_ref = df_ref[units].values
    norm_sq = np.dot(vec_ref, vec_ref)
    k = np.dot(vec_sample, vec_ref) / norm_sq if norm_sq != 0 else 1.0
    df_scaled = df_ref.copy()
    df_scaled[units] = df_scaled[units] * k
    return df_scaled


def subtract_reference(df_sample: pd.DataFrame, df_ref: pd.DataFrame,
                       wavenumber_range: tuple[float, float],
                       units: str = 'A') -> tuple[np.ndarray, float, np.ndarray]:
    """
    Scale and subtract a reference spectrum from a sample spectrum over a
    specified wavenumber range.

    The reference is first interpolated onto the sample's wavenumber axis,
    then scaled using the least-squares projection within *wavenumber_range*.

    Parameters
    ----------
    df_sample : pd.DataFrame
        Sample spectrum.
    df_ref : pd.DataFrame
        Reference spectrum (may have a different wavenumber axis).
    wavenumber_range : (float, float)
        ``(k_min, k_max)`` – region used to compute the scaling factor.
    units : str
        Column to use ('A' or 'T').

    Returns
    -------
    y_diff : np.ndarray
        Difference spectrum (sample − scaled reference).
    k : float
        Scaling factor applied to the reference.
    y_ref_scaled : np.ndarray
        Scaled reference interpolated onto the sample axis.
    """
    k_min, k_max = wavenumber_range
    x_sample = df_sample['cm-1'].values
    y_sample = df_sample[units].values
    x_ref = df_ref['cm-1'].values
    y_ref = df_ref[units].values

    # Sort reference in ascending wavenumber order
    sort_idx = np.argsort(x_ref)
    x_ref_s, y_ref_s = x_ref[sort_idx], y_ref[sort_idx]

    y_ref_interp = np.interp(x_sample, x_ref_s, y_ref_s, left=0, right=0)

    mask = (x_sample >= k_min) & (x_sample <= k_max)
    if mask.sum() == 0:
        mask = np.ones_like(x_sample, dtype=bool)

    v_s, v_r = y_sample[mask], y_ref_interp[mask]
    norm_sq = np.dot(v_r, v_r)
    k = np.dot(v_s, v_r) / norm_sq if norm_sq != 0 else 1.0

    y_ref_scaled = y_ref_interp * k
    return y_sample - y_ref_scaled, k, y_ref_scaled


# ===========================================================================
# ASTM E2412 AREA CALCULATIONS
# ===========================================================================
# All functions receive a DataFrame with columns 'cm-1' and 'A'.
# Net area is computed by trapezoidal integration after subtracting a
# two-point linear baseline constructed from local absorbance minima.
#
# Band definitions follow ASTM E2412 Table 1.
# ===========================================================================

def _net_area(df: pd.DataFrame,
              baseline_region_1: tuple[float, float],
              baseline_region_2: tuple[float, float],
              measurement_region: tuple[float, float]) -> float:
    """
    Internal helper: compute the net absorbance area of a band.

    A straight-line baseline is drawn between the absorbance minima found
    in *baseline_region_1* and *baseline_region_2*.  The net area is the
    integral of (absorbance − baseline) over *measurement_region*.

    Returns
    -------
    float
        Net area (≥ 0).  Returns 0.0 if any region is absent from *df*.
    """
    def _min_point(region):
        lo, hi = region
        mask = (df['cm-1'] >= lo) & (df['cm-1'] <= hi)
        if mask.sum() == 0:
            return None, None
        idx = df.loc[mask, 'A'].idxmin()
        return df.loc[idx, 'cm-1'], df.loc[idx, 'A']

    x1, y1 = _min_point(baseline_region_1)
    x2, y2 = _min_point(baseline_region_2)
    if x1 is None or x2 is None:
        return 0.0

    m = (y2 - y1) / (x2 - x1) if (x2 - x1) != 0 else 0.0
    b = y1 - m * x1

    lo, hi = measurement_region
    zone = df.loc[(df['cm-1'] >= lo) & (df['cm-1'] <= hi)].copy()
    if zone.empty:
        return 0.0

    net = zone['A'] - (m * zone['cm-1'] + b)
    area = float(-np.trapezoid(net, zone['cm-1']))
    return max(0.0, area)


def area_oxidation(df: pd.DataFrame) -> float:
    """
    Oxidation band area (1800–1670 cm⁻¹) – ASTM E2412.

    Baseline anchored at the minima in 550–650 cm⁻¹ and 1900–2200 cm⁻¹.

    Parameters
    ----------
    df : pd.DataFrame
        Spectrum with columns ``cm-1`` and ``A``.

    Returns
    -------
    float
        Net integrated area (absorbance·cm⁻¹).
    """
    return _net_area(df,
                     baseline_region_1=(550, 650),
                     baseline_region_2=(1900, 2200),
                     measurement_region=(1670, 1800))


def area_nitration(df: pd.DataFrame) -> float:
    """
    Nitration band area (1600–1650 cm⁻¹) – ASTM E2412.

    Baseline anchored at the minima in 550–650 cm⁻¹ and 1900–2200 cm⁻¹.

    Parameters
    ----------
    df : pd.DataFrame
        Spectrum with columns ``cm-1`` and ``A``.

    Returns
    -------
    float
        Net integrated area (absorbance·cm⁻¹).
    """
    return _net_area(df,
                     baseline_region_1=(550, 650),
                     baseline_region_2=(1900, 2200),
                     measurement_region=(1600, 1650))


def area_sulfation(df: pd.DataFrame) -> float:
    """
    Sulfation band area (1120–1180 cm⁻¹) – ASTM E2412.

    Baseline anchored at the minima in 550–650 cm⁻¹ and 1900–2200 cm⁻¹.

    Parameters
    ----------
    df : pd.DataFrame
        Spectrum with columns ``cm-1`` and ``A``.

    Returns
    -------
    float
        Net integrated area (absorbance·cm⁻¹).
    """
    return _net_area(df,
                     baseline_region_1=(550, 650),
                     baseline_region_2=(1900, 2200),
                     measurement_region=(1120, 1180))


def area_water(df: pd.DataFrame) -> float:
    """
    Water contamination band area (3150–3500 cm⁻¹) – ASTM E2412.

    Baseline anchored at the minima in 1900–2200 cm⁻¹ and 3680–4000 cm⁻¹.

    Parameters
    ----------
    df : pd.DataFrame
        Spectrum with columns ``cm-1`` and ``A``.

    Returns
    -------
    float
        Net integrated area (absorbance·cm⁻¹).
    """
    return _net_area(df,
                     baseline_region_1=(1900, 2200),
                     baseline_region_2=(3680, 4000),
                     measurement_region=(3150, 3500))


def area_glycol(df: pd.DataFrame) -> float:
    """
    Ethylene glycol (coolant) band area (1030–1100 cm⁻¹) – ASTM E2412.

    Baseline anchored at the minima in 1010–1030 cm⁻¹ and 1100–1130 cm⁻¹.

    Parameters
    ----------
    df : pd.DataFrame
        Spectrum with columns ``cm-1`` and ``A``.

    Returns
    -------
    float
        Net integrated area (absorbance·cm⁻¹).
    """
    return _net_area(df,
                     baseline_region_1=(1010, 1030),
                     baseline_region_2=(1100, 1130),
                     measurement_region=(1030, 1100))


def area_fuel_petrol(df: pd.DataFrame) -> float:
    """
    Petrol (gasoline) fuel dilution band area (745–755 cm⁻¹) – ASTM E2412.

    Baseline anchored at the minima in 730–750 cm⁻¹ and 760–780 cm⁻¹.

    Parameters
    ----------
    df : pd.DataFrame
        Spectrum with columns ``cm-1`` and ``A``.

    Returns
    -------
    float
        Net integrated area (absorbance·cm⁻¹).
    """
    return _net_area(df,
                     baseline_region_1=(730, 750),
                     baseline_region_2=(760, 780),
                     measurement_region=(745, 755))


def area_fuel_diesel(df: pd.DataFrame) -> float:
    """
    Diesel fuel dilution band area (805–815 cm⁻¹) – ASTM E2412.

    Baseline anchored at the minima in 795–805 cm⁻¹ and 825–835 cm⁻¹.

    Parameters
    ----------
    df : pd.DataFrame
        Spectrum with columns ``cm-1`` and ``A``.

    Returns
    -------
    float
        Net integrated area (absorbance·cm⁻¹).

    Notes
    -----
    The raw integrated area is returned. Convert to a concentration using
    a laboratory-specific calibration curve.
    """
    return _net_area(df,
                     baseline_region_1=(795, 805),
                     baseline_region_2=(825, 835),
                     measurement_region=(805, 815))


def area_antiwear_zddp(df: pd.DataFrame) -> float:
    """
    ZDDP antiwear additive band area (960–1025 cm⁻¹) – ASTM E2412.

    Baseline anchored at the minima in 550–650 cm⁻¹ and 1900–2200 cm⁻¹.

    Parameters
    ----------
    df : pd.DataFrame
        Spectrum with columns ``cm-1`` and ``A``.

    Returns
    -------
    float
        Net integrated area (absorbance·cm⁻¹).
    """
    return _net_area(df,
                     baseline_region_1=(550, 650),
                     baseline_region_2=(1900, 2200),
                     measurement_region=(960, 1025))


def soot_load(df: pd.DataFrame) -> float:
    """
    Soot (carbon black) load index – ASTM E2412.

    The absorbance at the carbon-black scattering reference point (2000 cm⁻¹)
    is multiplied by 100 to yield a dimensionless index.

    Parameters
    ----------
    df : pd.DataFrame
        Spectrum with columns ``cm-1`` and ``A``.

    Returns
    -------
    float
        Soot index (dimensionless, ≥ 0).
    """
    if df.empty:
        return 0.0
    idx_2000 = int(np.argmin(np.abs(df['cm-1'].values - 2000)))
    value = float(df.iloc[idx_2000]['A']) * 100
    return max(0.0, value)


def astm_areas(df: pd.DataFrame) -> dict:
    """
    Compute all standard ASTM E2412 band areas in a single call.

    Parameters
    ----------
    df : pd.DataFrame
        Spectrum with columns ``cm-1`` and ``A``.

    Returns
    -------
    dict
        Keys: ``oxidation``, ``nitration``, ``sulfation``, ``water``,
        ``glycol``, ``fuel_petrol``, ``fuel_diesel``, ``antiwear_zddp``,
        ``soot``.
    """
    return {
        'oxidation':    area_oxidation(df),
        'nitration':    area_nitration(df),
        'sulfation':    area_sulfation(df),
        'water':        area_water(df),
        'glycol':       area_glycol(df),
        'fuel_petrol':  area_fuel_petrol(df),
        'fuel_diesel':  area_fuel_diesel(df),
        'antiwear_zddp': area_antiwear_zddp(df),
        'soot':         soot_load(df),
    }


# ===========================================================================
# CONCENTRATION ESTIMATION
# ===========================================================================

def area_to_concentration(area: float,
                           calibration_areas: np.ndarray,
                           calibration_concentrations: np.ndarray,
                           poly_order: int = 3) -> float:
    """
    Convert a band area to a concentration using a polynomial calibration curve.

    Build and apply your own calibration curve by providing known
    (area, concentration) pairs measured on your instrument.

    Parameters
    ----------
    area : float
        Measured band area to convert.
    calibration_areas : array-like
        Band areas of the calibration standards.
    calibration_concentrations : array-like
        Corresponding concentrations (same units as desired output).
    poly_order : int
        Degree of the fitting polynomial (default 3).

    Returns
    -------
    float
        Estimated concentration (≥ 0).

    Example
    -------
    >>> # Example: water content calibration
    >>> cal_areas = np.array([0.0, 70.0, 87.0, 120.0, 144.0])
    >>> cal_conc  = np.array([0.0,  2.0,  4.0,   5.0,  10.0])  # % v/v
    >>> conc = area_to_concentration(measured_area, cal_areas, cal_conc)
    """
    cal_a = np.asarray(calibration_areas, dtype=float)
    cal_c = np.asarray(calibration_concentrations, dtype=float)
    coeff = np.polyfit(cal_a, cal_c, poly_order)
    poly = np.poly1d(coeff)
    return max(0.0, float(poly(area)))


# ===========================================================================
# HIT QUALITY INDEX (HQI) LIBRARY SEARCH
# ===========================================================================

def hqi(x_query: np.ndarray, y_query: np.ndarray,
        x_ref: np.ndarray, y_ref: np.ndarray) -> float:
    """
    Compute the Hit Quality Index (HQI) between two spectra.

    HQI is defined as the squared cosine similarity (0–100 %):

        HQI = 100 × (A·B)² / (‖A‖² · ‖B‖²)

    The reference is interpolated onto the query's wavenumber axis using
    the overlapping region only.

    Parameters
    ----------
    x_query, y_query : np.ndarray
        Query spectrum (wavenumber, absorbance).
    x_ref, y_ref : np.ndarray
        Reference spectrum.

    Returns
    -------
    float
        HQI score in percent (0–100).
    """
    x_q = x_query[np.argsort(x_query)]
    y_q = y_query[np.argsort(x_query)]
    x_r = x_ref[np.argsort(x_ref)]
    y_r = y_ref[np.argsort(x_ref)]

    lo = max(x_q[0], x_r[0])
    hi = min(x_q[-1], x_r[-1])
    if lo >= hi:
        return 0.0

    mask = (x_q >= lo) & (x_q <= hi)
    y_q_c = y_q[mask]
    y_r_c = np.interp(x_q[mask], x_r, y_r)

    dot = np.dot(y_q_c, y_r_c)
    norm_q = np.dot(y_q_c, y_q_c)
    norm_r = np.dot(y_r_c, y_r_c)

    if norm_q == 0 or norm_r == 0:
        return 0.0
    return float((dot ** 2) / (norm_q * norm_r) * 100)


def _hqi_worker(filepath: Path,
                x_query: np.ndarray,
                y_query: np.ndarray) -> dict:
    """Process a single candidate file for HQI comparison (worker function)."""
    try:
        ext = filepath.suffix.lower()
        if ext == '.sp':
            data, _ = load_spectrum(str(filepath), 'A')
            data = remove_baseline(data, 'A')
            x_c, y_c = data['cm-1'].values, data['A'].values
        else:
            with open(filepath) as fh:
                first = fh.readline()
            delim = ',' if ',' in first else None
            arr = np.loadtxt(filepath, delimiter=delim, comments='#')
            x_c, y_c = arr[:, 0], arr[:, 1]

        score = hqi(x_query, y_query, x_c, y_c)
        return {'file': str(filepath), 'hqi': score}
    except Exception as exc:
        return {'file': str(filepath), 'hqi': -1.0, 'error': str(exc)}


def library_search(query_path: str, library_dir: str,
                   top_n: int = 5,
                   extension: str = '*.sp') -> list[dict]:
    """
    Parallel brute-force HQI library search.

    Loads the query spectrum, then compares it against every file matching
    *extension* in *library_dir* (and subdirectories) using the Hit Quality
    Index.

    Parameters
    ----------
    query_path : str
        Path to the query spectrum file.
    library_dir : str
        Root directory of the reference library.
    top_n : int
        Number of top matches to display.
    extension : str
        Glob pattern for candidate files (e.g. ``'*.sp'``, ``'*.txt'``).

    Returns
    -------
    list[dict]
        All valid results sorted by HQI (descending), each with keys
        ``'file'`` and ``'hqi'``.
    """
    data_q, _ = load_spectrum(query_path)
    data_q = remove_baseline(data_q)
    x_q, y_q = data_q['cm-1'].values, data_q['A'].values

    candidates = list(Path(library_dir).rglob(extension))
    print(f"Searching {len(candidates)} files ({extension}) …")

    t0 = time.time()
    with concurrent.futures.ThreadPoolExecutor() as pool:
        futures = [pool.submit(_hqi_worker, p, x_q, y_q) for p in candidates]
        results = [f.result() for f in concurrent.futures.as_completed(futures)]
    elapsed = time.time() - t0
    print(f"Search completed in {elapsed:.2f} s\n")

    valid = sorted(
        [r for r in results if r['hqi'] >= 0],
        key=lambda r: r['hqi'],
        reverse=True
    )
    print(f"--- TOP {top_n} MATCHES ---")
    for i, match in enumerate(valid[:top_n], 1):
        print(f"  {i}. {Path(match['file']).name}  →  HQI: {match['hqi']:.2f} %")

    return valid


# ===========================================================================
# PLOTTING
# ===========================================================================

def plot_spectra(data_dict: dict,
                 title: str = "FTIR Spectroscopy",
                 xlim: tuple = (4000, 450),
                 ylim_A: tuple = (0, 0.4),
                 ylim_T: tuple = (0, 1),
                 save_path: str | None = None) -> tuple:
    """
    Plot one or more FTIR spectra (absorbance and transmittance panels).

    Parameters
    ----------
    data_dict : dict
        ``{'Label': {'df': df, 'color': 'black', 'alpha': 1.0, 'linewidth': 2}}``
        Only ``'df'`` is required per entry.
    title : str
        Figure title.
    xlim : (float, float)
        Wavenumber axis limits (default: 4000–450 cm⁻¹).
    ylim_A : (float, float)
        Absorbance axis limits.
    ylim_T : (float, float)
        Transmittance axis limits.
    save_path : str or None
        If given, the figure is saved to this path (PNG/PDF/SVG).

    Returns
    -------
    fig, axes : matplotlib Figure and Axes array
    """
    fig, axes = plt.subplots(ncols=2, figsize=(18, 7), dpi=150)

    for label, cfg in data_dict.items():
        df = cfg['df']
        kw = dict(
            linewidth=cfg.get('linewidth', 1.5),
            label=label,
            color=cfg.get('color', None),
            alpha=cfg.get('alpha', 1.0),
        )
        axes[0].plot(df['cm-1'], df['A'], **kw)
        if 'T' in df.columns:
            axes[1].plot(df['cm-1'], df['T'], **kw)

    fig.suptitle(title, fontsize=16, fontweight='bold')

    for ax in axes:
        ax.set_xlim(*xlim)
        ax.set_xlabel(r'Wavenumber [cm$^{-1}$]', fontsize=12)
        ax.legend(fontsize=11)
        ax.grid(True, linestyle='--', alpha=0.4)

    axes[0].set_ylim(*ylim_A)
    axes[0].set_ylabel('Absorbance', fontsize=12)
    axes[1].set_ylim(*ylim_T)
    axes[1].set_ylabel('Transmittance', fontsize=12)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches='tight')
    plt.show()
    return fig, axes
