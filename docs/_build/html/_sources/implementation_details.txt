.. _implementation_details:

Implementation Details
======================



Workflow
--------

* Kick off with 

.. code-block:: bash

    cord from_config my_scenario.conf

or 

.. code-block:: bash

    cord interactive --data-dir=path/to/dir --initial-veg-map=path/to/data/veg.asc ...

where the user will be promted to supply any arguments the user does not 
provide.


Conversions
-----------

* vegetation.asc to n.pol

* DFLOW boundary condition calculation

* shear.nc to shear.asc


Data Export via Hydroshare
--------------------------

* zip only a few important files, can be added to if necessary
* push entire modelrun to a single "resource"
* not yet implemented, but we could then have R, MATLAB scripts to pull data, 
  maybe make different visualizations/do other processing
