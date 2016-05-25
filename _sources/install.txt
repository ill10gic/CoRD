.. _install:

Install Instructions
====================

You will need to install 
HDF5 libraries v1.8.16 and the NetCDF libraries v4.3.3.1. On OS X for example,
this can be done using `Homebrew <http://brew.sh>`_. On Linux, use your package 
manager, making sure to install the version that includes required libraries.
On Debian/Ubuntu's apt-get this is done by appending the ``dev`` tag in the
proper way. Refer to the online documentation for your package management tool
to learn how to install a specific software version. To install the ``cord`` 
package on OS X you may need to set some environment variables by running the 
following commands.  You'll know you have to run these if the installation 
steps below fail.

.. code-block:: bash

    export HDF5_INCDIR=/usr/local/include/
    export USE_NCCONFIG=0

Once you have these two dependencies installed, you can install the package
by running 

.. code-block:: bash

    pip install cord

You can verify that the installation completed successfully by running

.. code-block:: bash

    cord

which should print the following help message

.. code-block:: none

    Usage: cord [OPTIONS] COMMAND [ARGS]

	Options:
	  --debug
	  --logfile TEXT
	  --help          Show this message and exit.

	Commands:
	  from_config  Run CoRD with params from <config_file>
	  interactive  Run CoRD interactively or with options
	  post_hs      Post the model run data to HydroShare



Install with Anaconda
---------------------

On CARC we have used Anaconda to allow individual users to install the proper
version of the dependencies which may differ from the version that is installed
for all system users. Please contact your administrator to help you install the
requirements that are listed in the 
`requirements.txt file <https://github.com/VirtualWatershed/CoRD/blob/master/requirements.txt>`_ 
of the repository. After the requirements have been installed and you have
activated your conda environment appropriately, you can run ``pip install cord``
as explained above.
