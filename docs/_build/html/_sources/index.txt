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
      generate_config Generate a configuration file
	  interactive  Run CoRD interactively or with options
	  post_hs      Post the model run data to HydroShare

Follow this or see the section :doc:`/usage` for more information. Note to get 
the help for a specific subcommand, you can run ``cord <subcommand>`` with no
arguments. For example, try ``cord from_config``.

Use CoRD to run simulations for a series of streamflows
-------------------------------------------------------

This will kick off a process that takes days to complete. Each individual run of
DFLOW is given a wall time of eight hours. There are thirty years of DFLOW runs
to finish, so we are looking at 8*30 = 240 hours of modeling. 
As such, we don't want to rely on 
our connection to Ulam or whatever other CARC machine to be sustained this whole
time. We'll use the Unix tool ``screen`` to create a session that will continue
after we disconnect our terminal from the server.

Here is a preview of the steps we go through to use screen to run a CoRD model
series:

1. Create a new screen by name that we can go back to after logging out then logging back in
2. Kick off a new series of CoRD modelruns
3. Exit the screen; the CoRD modelruns will continue until finished
4. Check progress by looking at the logfile and the PBS queue
5. To see if the process has stopped, "re-attach" to the screen by using its name and see if any informative output has been printed in the screen

Once we are sure CoRD has performed the modelruns as expected, we can use the
``post_hs`` utility to post our data to HydroShare for offline analysis.

``screen`` is installed on CARC systems so there is no need to install it. We can
get right to using it.

Creating and working with a new screen
``````````````````````````````````````

To create a new screen with a name, run

.. code-block:: bash
    
    screen -S <screen-name>

It's best to name your new screen something memorable and related to what the
CoRD scenario is you're running. Say our peak flows are simulating 100-year 
floods every 3 years. You might create a new screen called ``100yrEvery3`` like so

.. code-block:: bash

    screen -S 100yrEvery3

After this, you'll be in your newly-created screen named ``100yrEvery3``. You can
exit the screen by pressing ``Ctrl-A`` then ``D`` (don't press shift: if shift is
used that will be explicitly stated). The steps you just performed were

1. Creating and 'attaching' to a new screen named ``100yrEvery3``
2. 'Detaching' from that screen

Now let's add one more step, let's re-attach to that screen. To do this,

.. code-block:: bash

    screen -r 100yrEvery3

Now you will see the environment you had previously detached from. If you had
run any commands in the screen you should be able to see their outputs now.

One useful command for working within a ``screen`` is enabling "copy mode" which
incidentally also allows you to scroll. To do this, press ``Ctrl-A`` then ``[``. To
exit copy mode, press ``Esc``.

To list all available screens: 

.. code-block:: bash

    screen -list

To kill the existing screen ``100yrEvery3`` and all processes running within it: 

.. code-block:: bash

    screen -S 100yrEvery3 -X quit

which basically says "using screen with name 100yrEvery3 execute the quit
command". 

Next we'll show how to start a modelrun series, then continue to talk about
how to 


CoRD modelrun series Part 1: generate config file
````````````````````````````````````````````````````````````

Now that we know how to use screen to create a persistent environment in which
we will run our models, let's get ready to run a CoRD modelrun series. To begin,
we need to generate a configuration file that we will edit and use as input
to the cord command ``cord from_config``. There is a default configuration file
included with the ``cord`` distribution you installed with ``pip install cord``,
and can be found in the 
`GitHub repository: cord/default.conf.template <https://github.com/VirtualWatershed/CoRD/blob/master/cord/default.conf.template>`_.
You can copy this to your current directory and give it the name ``my.conf`` by
running 

.. code-block:: bash

    cord generate_config -n my.conf

If you simply run ``cord generate_config`` it will make a copy of
``default.conf.template`` in your current directory, but give it the name
``default.conf``. Currently this will overwrite any previously generate
configuration file with the same name, so be careful.

The next step in generating a configuration file is modifying the template file
that you've just copied to your current directory. If you look at your
freshly-copied config file, you'll see that most fields are blank. The fields
that are filled in are default stream channel and floodplain parameters. These
must be present or else running ``cord from_config`` will fail.

Two other fields are required: ``DATA_DIR`` and ``PEAK_FLOWS_FILE``.
``DATA_DIR`` must be a directory that already exists--this is where all the
RipCAS and DFLOW files will be stored for the modelrun series. The
``PEAK_FLOWS_FILE`` is the series of peak yearly streamflows in cubic meters per
second. You can :download:`download an example peak flows file here 
<_static/peak.txt>`. The peak flows text file must have ``Peak.Flood`` as the
first header line. 

If not specified, all the other config file fields will be filled in
with the latest data we have for the Jemez watershed. Currently we don't use the
HydroShare configuration info at all.


End-to-end example
------------------

CoRD modelrun series Part 2: Make a new screen and run ``cord from_config``
```````````````````````````````````````````````````````````````````````````

Now, let's assume you have a text file of peak flows that simulate 100-year
floods every three years called ``100yrEvery3.txt`` and you'll use this to 
drive the CoRD simulations.  Let's first create a config file that is 
descriptive of the scenario at hand:

.. code-block:: bash
    
    cord generate_config -n 100yrEvery3.conf

OK, now let's say your data will go to a directory called ``100-3`` to simplify
a little bit. Create this directory (``mkdir 100-3``). Now edit
``100yrEvery3.conf`` to look like this

.. code-block:: cfg

    [General]

    DATA_DIR = 100-3 # must create the directory before running 
    INITIAL_VEGETATION_MAP = 
    VEGZONE_MAP = 
    VEG_ROUGHNESS_SHEARRES_LOOKUP = 
    PEAK_FLOWS_FILE = 100yrEvery3.txt
    GEOMETRY_FILE = 

    STREAMBED_ROUGHNESS = 0.035  # in-channel only
    STREAMBED_FLOODPLAIN_ROUGHNESS = 0.04  # including out-of-channel for BC calculation
    STREAMBED_SLOPE = 0.001  # also used in BC calculation

    #
    # python import style, e.g. my_dflow_module.my_dflow_fun 
    # defaults to a function that calls qsub dflow_mpi.pbs
    DFLOW_RUN_FUN =   
    # will default to the data_dir with dashes replacing slashes if blank
    LOG_F = 
    
Now we will create a screen of the same name and start our modelrun series.

.. code-block:: bash

    screen -S 100yrEvery3

Now in your new screen

.. code-block:: bash

    cord --log-file=100-3.log from_config 100yrEvery3.conf

Note here we set a custom log file so that we can easily see the model progress
as it runs.

Exit the screen. As the model runs you can check two things to get an idea of
the progress. First, check the last few lines of the log file like so

.. code-block:: bash

    tail 100-3.log

and investigate the PBS job queue for your username. For me that command is

.. code-block:: bash

    qstat -u maturner

When the process has finished you can export your data to HydroShare by 
running the following command, substituting your username and password where
appropriate. Also change the resource title and keywords as you see fit. This
will help you, your colleagues, and other researchers discover your data later.

.. code-block:: bash

    cord post_hs --modelrun-dir=100-3 \
        --resource-title='100 year floods every three years' \
        --username=my-username \
        --password=my-password \
        --keyword=Jemez --keyword='climate change'



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
