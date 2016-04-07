# vw-Jemez

Functions for running our modified CASiMiR model and D-FLOW coupled
model on a supercomputer.

First we show some use instructions then installation instructions.

# Usage

To get the first .pol of n-values for use in the first D-FLOW run, use the
`jemez/veg2npol.py` script. For example, if the vegetation ESRI .asc is
`data/vegclass_2z.asc` and we want to write our .pol of n-values to
`initial_n.pol`, we would run

```
python jemez/veg2npol.py data/vegclass_2z.asc initial_n.pol
```

To run (not)CASiMiR to use D-FLOW inputs and output a .pol of n-values, use
the `jemez/dflow_casimir.py` script. For example, if the path to the
output netCDF with shear stress from D-FLOW is `data/jemez_r02_map.nc`
and the path to our vegetation map is `data/vegclass_2z.asc`, we would
run (not)casimir by running

```
python jemez/dflow_casimir.py ~/local_data/dflow_outputs/jemez_r02_map.nc ~/local_data/casimir_out/veg-out-1.asc
```

# Installation

### 1. Clone the repo and cd in to the root directory

```bash
git clone https://github.com/VirtualWatershed/vw-jemez && cd vw-jemez
```

### 2. Use a virtual environment and install dependencies

The `virtualenv` command used below can be installed with pip: `pip install virtualenv`.

Then with `virtualenv` installed, run the following

```bash
virtualenv venv
```

```bash
source venv/bin/activate
```

```bash
pip install -r requirements.txt
```

## Check installation by running unit tests

To check that all is well, try running the unit tests:

```bash
nosetests -v
```


The output should be

```
asc2pol should create proper headers and formatted data ... ok
test_casimir (test.test_dflow_casimir.TestDflow) ... ok
Test conversion of vegetation map to Manning's roughness map ... ok
test_vegmap_properly_read (test.test_dflow_casimir.TestDflow) ... ok

----------------------------------------------------------------------
Ran 4 tests in 0.384s

OK
```
