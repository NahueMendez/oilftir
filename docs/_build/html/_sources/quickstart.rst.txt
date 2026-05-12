.. _quickstart:

Quick Start
===========

Loading a spectrum
------------------

``oilftir`` accepts two-column CSV/TXT files (wavenumber, amplitude) and
PerkinElmer ``.sp`` binary files:

.. code-block:: python

   from oilftir import load_spectrum

   # From a plain-text CSV
   df, meta = load_spectrum("used_oil.csv", magnitude="A")

   # From a PerkinElmer .sp file (requires specio)
   df, meta = load_spectrum("sample.sp", magnitude="A")

   print(df.head())

The returned DataFrame always has columns ``cm-1``, ``A`` (absorbance),
and ``T`` (transmittance 0–1).

Baseline correction
-------------------

.. code-block:: python

   from oilftir import remove_baseline

   df = remove_baseline(df, units="A", poly_order=5)

Computing ASTM E2412 band areas
--------------------------------

All nine standard parameters in a single call:

.. code-block:: python

   from oilftir import astm_areas

   areas = astm_areas(df)
   # {
   #   'oxidation':      0.2341,
   #   'nitration':      0.0812,
   #   'sulfation':      0.1503,
   #   'water':          4.7620,
   #   'glycol':         0.0000,
   #   'fuel_petrol':    0.0034,
   #   'fuel_diesel':    0.0198,
   #   'antiwear_zddp':  8.3410,
   #   'soot':          12.5000,
   # }

Or individually:

.. code-block:: python

   from oilftir import area_oxidation, area_water, soot_load

   ox  = area_oxidation(df)
   h2o = area_water(df)
   sc  = soot_load(df)

Concentration from calibration curve
-------------------------------------

.. code-block:: python

   import numpy as np
   from oilftir import area_to_concentration

   cal_areas = np.array([0.0,  70.0,  87.0, 120.0, 144.0])
   cal_conc  = np.array([0.0,   2.0,   4.0,   5.0,  10.0])  # % v/v

   water_pct = area_to_concentration(
       areas["water"],
       calibration_areas=cal_areas,
       calibration_concentrations=cal_conc,
       poly_order=3,
   )
   print(f"Water content: {water_pct:.2f} % v/v")

Reference subtraction
---------------------

.. code-block:: python

   from oilftir import subtract_reference
   import pandas as pd

   df_ref, _ = load_spectrum("fresh_oil.csv", magnitude="A")
   df_ref = remove_baseline(df_ref)

   y_diff, k, _ = subtract_reference(
       df_sample=df,
       df_ref=df_ref,
       wavenumber_range=(1000, 1800),
   )

   df_diff = pd.DataFrame({"cm-1": df["cm-1"].values, "A": y_diff})

Plotting
--------

.. code-block:: python

   from oilftir import plot_spectra

   plot_spectra(
       data_dict={
           "Used oil":        {"df": df,      "color": "steelblue"},
           "Fresh oil (ref)": {"df": df_ref,  "color": "gray", "alpha": 0.6},
           "Difference":      {"df": df_diff, "color": "firebrick"},
       },
       title="Engine Oil – FTIR Analysis",
       xlim=(4000, 450),
       ylim_A=(0, 0.5),
       save_path="output/comparison.png",
   )

Library search (HQI)
---------------------

.. code-block:: python

   from oilftir import library_search

   results = library_search(
       query_path="unknown_oil.sp",
       library_dir="library/",
       top_n=5,
       extension="*.sp",
   )

   for rank, match in enumerate(results[:5], 1):
       print(f"{rank}. {match['file']}  →  HQI: {match['hqi']:.2f} %")
