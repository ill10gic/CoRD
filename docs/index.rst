.. CoRD: Coupled RipCAS-DFLOW modeling documentation master file, created by
   sphinx-quickstart on Wed May  4 16:42:41 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CoRD's documentation!
================================

CoRD is the Coupled RipCAS-DFLOW model developed as part of the Virtual
Watershed project of the Western Consortium for Watershed Analysis,
Visualization, and Exploration (WC-WAVE). This software contains
the RipCAS model itself, which is an open-source, Riparian-zone specific
version of the 
`CASiMiR vegetation model <http://www.casimir-software.de/ENG/veg_eng.html>`_.
There is no standalone wrapper for DFLOW, but we do provide convenience
functions for executing "model runs" in which DFLOW is first run with some
initial vegetation information to generate a map of shear stress, which is then
taken by RipCAS to determine if a vegetation community is destroyed and set to
age zero. There is quite a lot of data transformations happening under the hood,
which is documented more fully in the :ref:`implementation_details` section.


Contents:

.. toctree::
    :maxdepth: 2
    
    install
    cli
    implementation_details
    api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
