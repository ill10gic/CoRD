Install Instructions
--------------------

This package can be installed from the Python Packaging Index (PyPI). Since on
CARC users are recommended to use Anaconda for Python dependencies, 
the instructions below create a new Anaconda environment with all the necessary 
requirements. 

If you'd prefer not to use Anaconda, you will need to install 
HDF5 libraries v1.8.16 and the NetCDF libraries v4.3.3.1. On OS X for example,
this can be done using Homebrew. Refer to the online documentation for how to
install a specific software version. To install the ``cord`` package on OS X you
may need to set some environment variables by running the following commands.
You'll know you have to run these if the installation steps below fail.

.. code-block::

    export HDF5_INCDIR=/usr/local/include/
    export USE_NCCONFIG=0

.. TODO need to do PyPI publish before anaconda, need anaconda before I can
.. write these docs, so hold off on this branch and do the release!
