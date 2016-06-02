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

Contents:

.. toctree::
    :maxdepth: 2
    
    install
    usage
    api

.. TODO implementation_details left out for now TODO

There is no standalone wrapper for DFLOW, but we do provide convenience
functions for executing "model runs" in which DFLOW is first run with some
initial vegetation information to generate a map of shear stress, which is then
taken by RipCAS to determine if a vegetation community is destroyed and set to
age zero. 

It should be noted, then, that the DFLOW executable is required to be installed
and available on whatever system will run CoRD. You can get the source code and
find installation instructions for DFLOW on the `DFLOW website <https://www.deltares.nl/en/software/module/d-flow-flexible-mesh/>`_
(maybe...currently their open source software portal seems to be down and the
above link says that users must have a license).

.. There is quite a lot of data transformations happening under the hood, which is documented more fully in the :ref:`implementation_details` section.

To get the CoRD package, please visit the section :doc:`/install` 


Quickstart
==========

This guide is for users of the CARC shared system at UNM. We must load the
anaconda module then activate our Anaconda environment. First

.. code-block:: bash

    $ module load anaconda

Since ``conda`` is the preferred way here, to get started 

This guide is for CARC users who do not have access to the ``virtualenv`` command. 
If you prefer ``virtualenv`` then instead of the above command, run

.. code-block:: bash
    
   $ virtualenv cord

then activate your new virtual environment by running

.. code-block:: bash

   $ source cord/bin/activate


Now that you've created and activated an environment you can install ``cord``.
Note that whether you're using Anaconda or a virtual environment, the
environment name will precede your command line as shown in the snippets below.

.. code-block:: bash

    (cord)$ pip install cord

To see that cord is present, just run cord:

.. code-block:: bash

    (cord)$ cord

You should see the following help message. 

.. code-block:: none

    Usage: cord [OPTIONS] COMMAND [ARGS]...

	Options:
	  --debug
	  --logfile TEXT
	  --help          Show this message and exit.

	Commands:
	  from_config  Run CoRD with params from <config_file>
	  interactive  Run CoRD interactively or with options
	  post_hs      Post the model run data to HydroShare

Follow this or see the section :doc:`/usage` for more information. Note to get 
the help for a specific subcommand, you can run ``cord <subcommand>`` with no
arguments. For example, try ``cord from_config``.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
