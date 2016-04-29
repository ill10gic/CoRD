# vw-RipCAS-DFLOW

Functions for running our Riparian Community Alteration and Succession model,
RipCAS, and D-FLOW coupled.
One of the major contributions of this work is a process for automatically
submitting new DFLOW jobs to the Portable Batch System automatically, then
running RipCAS when DFLOW has finished for that time step. Roughly the process
proceeds as shown in the diagram below.

![simplified workflow](work_flow_simplified.png)

The DFLOW model is the one that runs for roughly seven hours at a time and uses
eight processors per model. We could use more cores, with an open issue (#9)
addressing this.


## Coupled ModelRun on a Cluster

As it is, we can vary quite a few things and automatically run a variety of
scenarios. One could vary the flooding by year, or change the geometry of the
cross-section of the stream which is used in calculating the boundary
conditions, or the user could create new vegetation type-to-n and/or vegetation
type-to-shear resistance mappings. This is done by creating new files that match
the structure of their counterparts in the
[data
directory](https://github.com/VirtualWatershed/vw-ripcas-dflow/tree/master/data).
For example, to modify the vegetation conversion rules, one must create an Excel
spreadsheet that matches this structure in the screenshot below.

![ripcas required data screenshot](ripcas-required-data-screenshot.png)


Example:

```
python ripcas_dflow/modelrun.py ~/ripcas-dflow-runs/100yr-every5 \
    data/ripcas_inputs/vegclass_2z.asc data/ripcas_inputs/zonemap_2z.asc \
    data/ripcas_inputs/ripcas-data-requirements.xlsx \
    /users/maturner/100yrFlood_every5.txt data/dflow_inputs/DBC_geometry.xyz \
    0.04 0.001
```

More details on the usage:

```
python ripcas_dflow/modelrun.py data_dir initial_vegetation vegzone_map ripcas_required_data\
        peak_flows_file geometry_file streambed_roughness streambed_slope

    data_dir: directory to hold each time step of the model run
    initial_vegetation: .asc file with initial vegetation map
    vegzone_map: .asc file with vegetation zone information
    ripcas_required_data: .xlsx file with veg type-to-n and shear resistance-per-veg type information
    peak_flows_file: file with a column of peak flood flow in cubic meters per second
    geometry_file: xyz file representing geometry of downstream cross section for boundary conditions calculation
    streambed_roughness: floating point number for the roughness value of the streambed
    streambed_slope: floating point number for the slope of the stream geography
```

Currently this is part of a
[raw executable in modelrun.py](https://github.com/VirtualWatershed/vw-ripcas-dflow/blob/master/ripcas_dflow/modelrun.py#L446),
but soon this will be broken out into its own controller as we add new features.


## Using the ModelRun class directly

Here is an example of using the ModelRun class to set up and execute the coupled
DFLOW/RipCAS model:

```python
import matplotlib.pyplot as plt

from ripcas_dflow import ModelRun, veg2n

mr = ModelRun()

# assume we have read these from a file or elsewhere; set them here for ex
peak_flow = 89.55
streambed_roughness = 0.04
reach_slope = 0.001

geometry = Pol.from_river_geometry_file('data/DBC_geometry.xyz')

mr.calculate_bc(peak_flow, geometry, streambed_roughness, reach_slope)

assert mr.bc_converged

mr.run_dflow('data/dflow-test/', 'data/vegclass_2z.asc')

# the output is an ESRIAsc map of vegetation type (coded integer)
out = mr.run_ripcas('data/zonemap_2z.asc', 'data/casimir-data-requirements.xlsx', 'data/ripcas-test')

# translate to Manning's roughness map, which shows communities a little better
n_out = veg2n(out)

plt.matshow(n_out.as_matrix(replace_nodata_val=0.0))
plt.colorbar()
```

![An example map](example_n_map.png)

### Boundary Condition Solver

We provide a solver for boundary conditions which can be used as follows. There
is some data in the `data/` directory that we use.

The boundary condition we have to calculate is actually an inverse problem. We
are given the peak flow for a given year, but we don't know what the elevation
of the water surface (WS elevation) is at the bottom of the reach under consideration. DFLOW
needs this as a boundary condition, as well as the streamflow that is associated
with the WS elevation, or just WS for short. We use
[scipy.minimize_scalar](http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize_scalar.html)
to accomplish this inverse problem. See the
[source code](https://github.com/VirtualWatershed/vw-ripcas-dflow/blob/master/ripcas_dflow/modelrun.py#L357)
for more details.

```python
from dflow_casimir import BoundaryConditionSovler, Pol

# these are variables we would probably read from a file
target_streamflow = 89.55  # the first value in data/peak.txt
streambed_roughness = 0.04  # empirical n value
reach_slope = 0.001

# use our Pol (polygon) object used for xyz/polygon file dialects of interest
geometry = Pol.from_river_geometry_file('data/DBC_geometry.xyz')

bc_solver = BoundaryConditionSover(target_streamflow,
                                   geometry,
                                   streambed_roughness,
                                   reach_slope)

bc_solution = bc_solver.solve()

print bc_solution
```

Will yield the `BoundaryConditionResult` printing

```
BoundaryConditionResult(ws_elev=1777.3393782057494,
streamflow=89.549822073608453, error=0.00017792639154379231, success=True)
```

which shows that for the given streamflow of 89.55 cubic meters per second,
the WS elevation is 1777.3394 meters. The absolute error between calculated
streamflow (89.549822...) and observed streamflow is 0.0002, or a relative error of
0.0002/89.55 = 2.23e-6.



## Command-line scripts

To get the first .pol of n-values for use in the first D-FLOW run, use the
`jemez/veg2npol.py` script. For example, if the vegetation ESRI .asc is
`data/vegclass_2z.asc` and we want to write our .pol of n-values to
`initial_n.pol`, we would run

```
python jemez/veg2npol.py data/vegclass_2z.asc initial_n.pol
```

To run RipCAS to use D-FLOW inputs and output a .pol of n-values, use
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
git clone https://github.com/VirtualWatershed/vw-ripcas-dflow && cd vw-ripcas-dflow
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

## Developer Notes

### Files required for a DFLOW run

There are quite a few files required for a DFLOW run. Here we have a table
of what input files are required and where the requirement is defined.


|                    File | Where Used/Other notes               | Needs to be generated at each timestep? |
|------------------------:|--------------------------------------|:---------------------------------------:|
|           dflow_mpi.pbs | Job script submitted via qsub        |                    No                   |
|               base.mdu  | dflow_mpi.pbs                        |                    No                   |
|            base_net.nc  | base.mdu                             |                    No                   |
|               base.ext  | base.mdu                             |                    No                   |
|                   n.pol | base.ext                             |                   Yes                   |
|      boundriverdown.pli | base.ext                             |                    No                   |
|       boundriver_up.pli | base.ext                             |                    No                   |
| boundriverdown_0001.cmp | ??? dflow automatically detected ??? |                   Yes                   |
|  boundriver_up_0001.cmp | ??? dflow automatically detected???  |                   Yes                   |

Note that this makes for a total of six files that do not change from one run to
another and three that do. The `boundriver*.cmp` are generated by running the
boundary condition calculator of the `ModelRun`. `n.pol` is generated from a
vegetation .asc either given as an initial input or output from the RipCAS step.

There is one final file in `data/dflow_inputs`, `DBC_geometry.xyz`, which is the
polygon file for the geometry of the cross section, required for generating the
boundary condition files. The `boundriver*.pli` do not change, but at this point
I (MT) am not really sure what is the purpose of these.
