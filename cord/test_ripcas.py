from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
import subprocess
import time

from modelrun import modelrun_series, ModelRun
from ripcas_dflow import ESRIAsc, stitch_partitioned_output, shear_mesh_to_asc, ripcas_with_dflow_io, ripcas

mr = ModelRun()
dflow_shear_output = \
            os.path.join('./testing-inputs/', 'stitched-shear.nc')
zone_map_path = './data/inputs/zonemap_2z.asc'
veg_map_path = './data/inputs/vegclass_2z.asc'
streambed_roughness = 0.035
ripcas_testout = './testing-outputs/ripout'
dflow_testout = './testing-outputs/dflow_Out'
required_ripcas_inputs = './data/inputs/veg_roughness_shearres.xlsx'
shear_asc = None
# pol = ripcas_with_dflow_io(ESRIAsc(veg_map_path), ESRIAsc(zone_map_path), streambed_roughness, dflow_shear_output, required_ripcas_inputs)


hdr = ESRIAsc(veg_map_path).header_dict()
mr.vegetation_ascii = ESRIAsc(veg_map_path)
print(hdr)
hdr = mr.vegetation_ascii.header_dict()
print('hdr dict')
print(hdr)
if shear_asc is None:
    shear_asc = shear_mesh_to_asc(dflow_shear_output, hdr)
else:
    assert isinstance(shear_asc, ESRIAsc),\
        'shear_asc must be of type ESRIAsc if provided'

# shear_asc.write(
#     os.path.join(dflow_testout, 'shear_out.asc')
# )

output_veg_ascii = ripcas(
    mr.vegetation_ascii, zone_map_path,
    shear_asc, ripcas_required_data_path
)

output_vegetation_path = os.path.join(
    ripcas_testout, 'vegetation.asc'
)
output_veg_ascii.write(output_vegetation_path)

mr.ripcas_has_run = True
# shear_asc = shear_mesh_to_asc(dflow_shear_output, hdr)
# print('created and asc for shear!')
# shear_asc.write(
#             os.path.join('./testing-outputs/', 'shear_out.asc')
#         )

print('finished test')


# print("hello")