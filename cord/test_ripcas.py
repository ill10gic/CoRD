import matplotlib.pyplot as plt
import os

from modelrun import *

mr = ModelRun()
mr.dflow_has_run = True
mr.dflow_run_directory = os.path.join('data', 'dflow-' + str(0))
vegetation_map = 'data/inputs/zonemap_2z.asc'
mr.vegetation_ascii = vegetation_map
ripcas_directory = 'data/ripcas-test'
#os.mkdir(ripcas_directory) # this line is required so that the dataset error will be fixed on line 273 in stitch partitioned output IOError: [Errno 13] Permission denied: 'data/ripcas-0/stitched-shear.nc'

# assume we have read these from a file or elsewhere; set them here for ex
peak_flow = 89.55
streambed_roughness = 0.04
reach_slope = 0.001

#geometry = Pol.from_river_geometry_file('data/DBC_geometry.xyz')

#mr.calculate_bc(peak_flow, geometry, streambed_roughness, reach_slope)

#assert mr.bc_converged

#mr.run_dflow('data/dflow-test/', 'data/vegclass_2z.asc')

# the output is an ESRIAsc map of vegetation type (coded integer)
out = mr.run_ripcas(vegetation_map, 'data/inputs/veg_roughness_shearres.xlsx', ripcas_directory)

# translate to Manning's roughness map, which shows communities a little better
n_out = veg2n(out)

plt.matshow(n_out.as_matrix(replace_nodata_val=0.0))
plt.colorbar()
