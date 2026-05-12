# example.py
"""
oilftir – Quick-start example
==============================
Demonstrates the main workflows of the oilftir library:

  1. Loading a spectrum from a plain-text file
  2. Baseline correction
  3. Computing all ASTM E2412 band areas at once
  4. Converting an area to concentration with a custom calibration curve
  5. Comparing a used-oil spectrum against a fresh-oil reference
  6. Plotting one or several spectra
  7. Running a library search (HQI)

Replace the file paths below with your own spectra.
"""


import os
import collections
import collections.abc
collections.Iterable = collections.abc.Iterable
import specio
import numpy as np


from oilftir.astm_utils import (
    load_spectrum,
    remove_baseline,
    astm_areas,
    area_to_concentration,
    subtract_reference,
    plot_spectra,
    library_search
)

if __name__ == '__main__':
    # ---------------------------------------------------------------------------
    # 1. Load a spectrum
    #    Supported formats: CSV, space-delimited TXT, or PerkinElmer .sp
    # ---------------------------------------------------------------------------
    script_dir = os.path.dirname(os.path.abspath(__file__))
    example_filepath = os.path.join(script_dir, "data", "ftir_example.sp")
    df, meta = load_spectrum(example_filepath, magnitude="A")
    print(df.head())
    
    # ---------------------------------------------------------------------------
    # 2. Baseline correction (optional but recommended)
    # ---------------------------------------------------------------------------
    df = remove_baseline(df, units="A", poly_order=5)
    
    # ---------------------------------------------------------------------------
    # 3. Compute all ASTM E2412 band areas in one call
    # ---------------------------------------------------------------------------
    areas = astm_areas(df)
    
    print("\n--- ASTM E2412 Band Areas ---")
    for param, value in areas.items():
        print(f"  {param:<20s}: {value:.4f}")
    
    # ---------------------------------------------------------------------------
    # 4. Convert water area to concentration using a calibration curve
    #    Supply your own (area, concentration) calibration pairs.
    # ---------------------------------------------------------------------------
    cal_areas = np.array([0.0,  70.0,  87.0, 120.0, 144.0])   # absorbance·cm⁻¹
    cal_conc  = np.array([0.0,   2.0,   4.0,   5.0,  10.0])   # % v/v water
    
    water_pct = area_to_concentration(
        areas["water"],
        calibration_areas=cal_areas,
        calibration_concentrations=cal_conc,
        poly_order=3,
    )
    print(f"\n  Estimated water content : {water_pct:.2f} % v/v")
    
    # ---------------------------------------------------------------------------
    # 5. Reference subtraction (used oil vs. fresh oil baseline)
    # ---------------------------------------------------------------------------
    freshoil_filepath = os.path.join(script_dir, "data", "fresh_oil.txt")
    df_ref, _ = load_spectrum(freshoil_filepath, magnitude="T")
    
    y_diff, k, y_ref_scaled = subtract_reference(
        df_sample=df,
        df_ref=df_ref,
        wavenumber_range=(1000, 1800),  # region used to compute scaling factor k
        units="A",
    )
    print(f"\n  Reference scaling factor k = {k:.4f}")
    
    # Build a DataFrame for the difference spectrum so it can be plotted
    import pandas as pd
    df_diff = pd.DataFrame({"cm-1": df["cm-1"].values, "A": y_diff})
    
    # ---------------------------------------------------------------------------
    # 6. Plot spectra
    # ---------------------------------------------------------------------------
    output_filepath= os.path.join(script_dir, "output", "comparison.png")
    plot_spectra(
        data_dict={
            "Used oil":       {"df": df,      "color": "steelblue"},
            "Fresh oil (ref)":{"df": df_ref,  "color": "gray",      "alpha": 0.6},
            "Difference":     {"df": df_diff, "color": "firebrick",  "linewidth": 1.0},
        },
        title="Engine Oil – FTIR Analysis",
        xlim=(4000, 450),
        ylim_A=(0, 5),
        save_path=output_filepath,   # omit to skip saving
    )
    
    # ---------------------------------------------------------------------------
    # 7. Library search (HQI) – identify an unknown oil type
    # ---------------------------------------------------------------------------
    # top_n    : number of best matches to display
    # extension: glob pattern; change to "*.txt" for text-based libraries
    
    unknownoil_filepath = os.path.join(script_dir, "data", "unknown_oil.sp")
    library_dir=os.path.join(script_dir, "data", "library")
    results = library_search(
        query_path=unknownoil_filepath,
        library_dir=library_dir,
        top_n=5,
        extension="*",
    )
    
    print("\n--- Top HQI matches ---")
    for rank, match in enumerate(results[:5], 1):
        print(f"  {rank}. {match['file']}  →  {match['hqi']:.2f} %")
