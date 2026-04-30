# __init__.py
from .astm_utils import (
    load_spectrum,
    remove_baseline,
    scale_reference,
    subtract_reference,
    astm_areas,
    area_oxidation,
    area_nitration,
    area_sulfation,
    area_water,
    area_glycol,
    area_fuel_petrol,
    area_fuel_diesel,
    area_antiwear_zddp,
    soot_load,
    area_to_concentration,
    hqi,
    library_search,
    plot_spectra,
)

__version__ = "0.1.0"
__author__  = "Nahuel Mendez"