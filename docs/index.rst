.. oilftir documentation master file

oilftir
=======

**FTIR-based condition monitoring of lubricating oils — ASTM E2412 toolkit**

``oilftir`` is a Python library for processing and analysing FTIR spectra of
used lubricating oils following the
`ASTM E2412 <https://www.astm.org/e2412-23.html>`_ standard.

Developed at the **Chemistry Laboratory of CENADIF**
(Centro Nacional de Desarrollo e Innovación Ferroviaria, Argentina)
for used-oil condition monitoring in the local railway industry.

.. code-block:: python

   from oilftir import load_spectrum, remove_baseline, astm_areas

   df, meta = load_spectrum("used_oil.csv", magnitude="A")
   df = remove_baseline(df)

   areas = astm_areas(df)
   print(areas)

----

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   quickstart
   autoapi/index

----

Authors
-------

- **Nahuel Mendez** — R&D Engineer / Physical Data Scientist
  (`GitHub <https://github.com/NahueMendez>`_ ·
  `LinkedIn <https://www.linkedin.com/in/ingnahuelmendez/>`_)
- **Hernán Gomez Molino**, Eng. — Head of Laboratories, CENADIF
- **Leandro Asens**, Eng. — Head of the Chemistry Laboratory, CENADIF

License
-------

MIT License — see ``LICENSE`` for details.

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
