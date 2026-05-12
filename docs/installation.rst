.. _installation:

Installation
============

From PyPI (recommended)
------------------------

.. code-block:: bash

   pip install oilftir

With PerkinElmer ``.sp`` file support
--------------------------------------

.. code-block:: bash

   pip install oilftir[specio]

From source (development / editable install)
---------------------------------------------

.. code-block:: bash

   git clone https://github.com/NahueMendez/oilftir.git
   cd oilftir
   pip install -e .

Requirements
------------

- Python ≥ 3.9
- ``numpy``
- ``pandas``
- ``scipy``
- ``matplotlib``
- ``pybaselines``
- ``specio`` *(optional — only for .sp files)*
