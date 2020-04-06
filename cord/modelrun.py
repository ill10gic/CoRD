"""
Functions to calculate the downstream water surface elevation by minimizing the
difference between flows calculated via the Manning Formula for discharge and
the historical peak flood values.

(https://en.wikipedia.org/wiki/Manning_formula)
(https://en.wikipedia.org/wiki/Volumetric_flow_rate)


Author:
    Matthew A. Turner <maturner01@gmail.com>
Date:
    19 April 2016
"""
from __future__ import print_function

import numpy as np
import os
import shutil
import subprocess
import time


from collections import namedtuple
from glob import glob
from scipy.optimize import minimize_scalar

try:
    from ripcas_dflow import ESRIAsc, Pol, ripcas, shear_mesh_to_asc, veg2n, stitch_partitioned_output
except ImportError:
    from .ripcas_dflow import ESRIAsc, Pol, ripcas, shear_mesh_to_asc, veg2n, stitch_partitioned_output


class ModelRun(object):
    """
    A single coupled run. First DFLOW then RipCAS. CoupledRunSequence will
    encapsulate a series of coupled runs commencing with preparation of the
    initial vegetation map for DFLOW. For now, assume that the vegetation map
    is provided to the run_dflow method.
    """
    def __init__(self):  # , vegetation_ascii_path):

        # have the boundary conditions been found?
        self.bc_converged = False

        # has ripcas been run yet?
        self.vegetation_ascii = None
        self.ripcas_has_run = False
        self.ripcas_directory = None

        # has DFLOW been run yet?
        self.dflow_has_run = False
        self.dflow_run_directory = None
        self.dflow_shear_output = None

        # generate boundry condition objects
        self.upstream_bc = BoundaryCondition()
        self.downstream_bc = BoundaryCondition()
        self.bc_solution_info = BoundarySolutionInfo()

    def calculate_bc(self, target_streamflow,
                     dbc_geometry_file, streambed_roughness, slope):
        """

        Arguments:
            target_streamflow (float): historical or other streamflow that
                will be used to drive DFLOW model; this calculation recovers
                an estimate for the Water Surface elevation (WS) for this given
                streamflow.
            dbc_geometry_file (str): path to the stream's cross-sectional
                geometry xyz file
            streambed_roughness (float): Manning's n-value for the streambed
            slope (float): slope taken for the reach

        Returns:
            (BoundaryCondition, BoundaryCondition): tuple of upstream and
                downstream BoundaryCondition instances
        """
        dbc_geometry = Pol.from_river_geometry_file(dbc_geometry_file)

        bc_solver = BoundaryConditionSolver(
            target_streamflow, dbc_geometry, streambed_roughness, slope
        )

        bc_solution = bc_solver.solve()

        self.bc_solution_info = bc_solution

        self.bc_converged = bc_solution.success

        self.downstream_bc.amplitude = bc_solution.ws_elev

        self.upstream_bc.amplitude = bc_solution.streamflow

        return (self.upstream_bc, self.downstream_bc)

    def run_dflow(self, dflow_run_directory, vegetation_map,
                  veg_roughness_shearres_lookup, streambed_roughness,
                  clobber=True, pbs_script_name='dflow_mpi.pbs',
                  dflow_run_fun=None):
        """
        Both input and output dflow files will go into the dflow_run_directory,
        but in input/ and output/ subdirectories.

        Arguments:
            dflow_run_directory (str): directory where DFLOW files should be
                put and where the dflow_run_fun will be run from
            vegetation_map (str): path to the input vegetation.asc file.
                Converted to n.pol.
            veg_roughness_shearres_lookup (str):
            clobber (bool): whether or not to overwrite dflow_run_directory if
                it exists
            pbs_script_name (str): name of .pbs script w/o directory
            dflow_run_fun (function): argument-free function to run DFLOW.
                Ex. `dflow_run_fun=f` where `f` defined by
                `def f: subprocess.call(['qsub', 'dflow_mpi.pbs'])`

        Returns:
            None
        """
        if not self.bc_converged:
            raise RuntimeError(
                'Boundary conditions must be calculated before ' +
                'DFLOW can be run'
            )

        if self.dflow_has_run:
            raise RuntimeError(
                'DFLOW has already been run for this CoupledRun'
            )

        if os.path.exists(dflow_run_directory):

            if not clobber:
                raise RuntimeError(
                    'DFLOW has already been run for this CoupledRun'
                )

            shutil.rmtree(dflow_run_directory)

        self.dflow_run_directory = dflow_run_directory

        curdir = os.path.dirname(__file__)

        path_to_dflow_inputs = os.path.join(
            curdir, 'data', 'dflow-partition'
        )

        shutil.copytree(path_to_dflow_inputs, dflow_run_directory)

        # self.dflow_shear_output =\
            # os.path.join(dflow_run_directory,
                         # 'DFM_OUTPUT_base',
                         # 'base_map.nc')

        # write boundary conditions to file
        bc_up_path = os.path.join(dflow_run_directory,
                                  'boundriver_up_0001.cmp')

        bc_down_path = os.path.join(dflow_run_directory,
                                    'boundriverdown_0001.cmp')

        self.upstream_bc.write(bc_up_path)
        self.downstream_bc.write(bc_down_path)

        self.vegetation_ascii = ESRIAsc(vegetation_map)

        roughness_path = os.path.join(dflow_run_directory, 'n.pol')

        # convert the vegetation .asc to roughness .pol, write to veg_path
        Pol.from_ascii(
            veg2n(self.vegetation_ascii,
                  veg_roughness_shearres_lookup,
                  streambed_roughness)
        ).write(roughness_path)

        bkdir = os.getcwd()
        os.chdir(dflow_run_directory)

        if dflow_run_fun is None:

            print('\n*****\nDry Run of DFLOW\n*****\n')

            os.chdir(bkdir)
            example_shear_path = 'stitched-shear.nc'
            
            self.dflow_shear_output = \
            os.path.join(dflow_run_directory, 'stitched-shear.nc')
            print(self.dflow_shear_output)
            # if os.path.exists(example_shear_path):
            #     #os.makedirs(os.path.dirname(self.dflow_shear_output))
            #     #shutil.copyfile(example_shear_path, self.dflow_shear_output)

            # else:
            #     print('Get you a copy of a DFLOW output, yo! ' +
            #           'Can\'t run RipCAS without it!')

            #     with open('not_actually_output.nc', 'w') as f:
            #         f.write('A FAKE NETCDF!!!')

            self.dflow_has_run = True

        else:

            # in the case of running a process on CARC, the ret is a Popen inst
            ret = dflow_run_fun()

            os.chdir(bkdir)
            self.dflow_has_run = True

            return ret

    def run_ripcas(self, zone_map_path, ripcas_required_data_path,
                   ripcas_directory, shear_asc=None, clobber=True):

        if not self.dflow_has_run:
            raise RuntimeError(
                'DFLOW must run before ripcas can be run'
            )

        # if os.path.exists(ripcas_directory):

            # if not clobber:
                # raise RuntimeError(
                    # 'DFLOW has already been run for this CoupledRun'
                # )

            # shutil.rmtree(ripcas_directory)

        self.ripcas_directory = ripcas_directory

        #uncommenting this to fix the bug that was occuring
        if os.path.isdir(ripcas_directory) is False:
            os.mkdir(ripcas_directory) #line 273, in stitch_partitioned_output, IOError: [Errno 13] Permission denied: 'data/ripcas-0/stitched-shear.nc'

        partitioned_outputs = \
            glob(
                os.path.join(self.dflow_run_directory, #this line fixed to 'dflow_run_directory'
                             'DFM_OUTPUT_base',
                             'base*map.nc')
            )

        self.dflow_shear_output = \
            os.path.join(ripcas_directory, 'stitched-shear.nc')

        stitch_partitioned_output(partitioned_outputs, self.dflow_shear_output)
        print ('zone_map_path: ' + zone_map_path)
        print ('dflow_run_directory: ' + self.dflow_run_directory)
        print ('dflow_shear_output: ' + self.dflow_shear_output)
        hdr = self.vegetation_ascii.header_dict()
        print('hdr dict')
        print(hdr)
        if shear_asc is None:
            shear_asc = shear_mesh_to_asc(self.dflow_shear_output, hdr)
        else:
            assert isinstance(shear_asc, ESRIAsc),\
                'shear_asc must be of type ESRIAsc if provided'

        shear_asc.write_unflattened_asc(
            os.path.join(self.dflow_run_directory, 'shear_out.asc')
        )

        output_veg_ascii = ripcas(
            self.vegetation_ascii, zone_map_path,
            shear_asc, ripcas_required_data_path
        )

        output_vegetation_path = os.path.join(
            ripcas_directory, 'vegetation.asc'
        )
        output_veg_ascii.write_unflattened_asc(output_vegetation_path)

        self.ripcas_has_run = True

        return output_veg_ascii


BoundarySolutionInfo = namedtuple(
    'BoundarySolutionInfo', ['ws_elev', 'streamflow', 'error', 'success']
)
BoundarySolutionInfo.__new__.__defaults__ = (None, None, None, None)


class BoundaryConditionSolver:

    def __init__(self,
                 historical_streamflow,
                 dbc_geometry,
                 streambed_roughness,
                 slope):

        self.q_hist = historical_streamflow
        self.geom = dbc_geometry
        self.n = streambed_roughness
        self.slope = slope

    def solve(self):

        def _streamflow_error(ws_elev):

            calc =\
                _calculate_streamflow(self.geom, self.n, ws_elev, self.slope)

            return abs(calc - self.q_hist)

        # generate initial guesses with wide-spaced points
        result = minimize_scalar(_streamflow_error,
                                 bounds=(self.geom.z.min(), self.geom.z.max()),
                                 method='bounded',
                                 options={'xatol': 1e-6, 'maxiter': 1000})

        return BoundarySolutionInfo(
            result.x,
            _calculate_streamflow(self.geom, self.n, result.x, self.slope),
            result.fun,
            result.success
        )


StreamflowTuple = namedtuple('StreamflowTuple', ['ws_elev', 'streamflow'])


def _calculate_streamflow(dbc_geometry, streambed_roughness,
                          water_surface_elevation, slope):
    # have N points; get N-1 distances and
    # N-1 Max/Min over those distances
    x = dbc_geometry.x
    y = dbc_geometry.y
    z = dbc_geometry.z

    dx = np.diff(x)
    dy = np.diff(y)

    xydist = np.sqrt(np.square(dx) + np.square(dy))
    # station = np.cumsum(xydist)

    zmax_by_segment = np.array(
        [max(z[i], z[i+1]) for i in range(len(z)-1)]
    )

    zmin_by_segment = np.array(
        [min(z[i], z[i+1]) for i in range(len(z)-1)]
    )

    # get N-1 vector taken from S = ('below', 'triangle', 'trap')
    # for the three possible positions of the water surface and
    # commensurate calculation methods for wetted perimeter
    ws_location = np.array(
        [
            _get_ws_location(water_surface_elevation, _z[0], _z[1])
            for _z in zip(zmax_by_segment, zmin_by_segment)
        ]
    )

    # calculate the cross-sectional area of the stream
    # at the lower bound
    area_vec = np.zeros(len(ws_location))
    # wetted perimeter
    wp_vec = np.zeros(len(ws_location))

    ws_elev = water_surface_elevation
    for idx, loc in enumerate(ws_location):

        if loc == 'triangle':

            zmin = zmin_by_segment[idx]
            zmax = zmax_by_segment[idx]
            xy = xydist[idx]

            # calculate area
            area_vec[idx] = 0.5 * (ws_elev - zmin) * xy

            # calculate wetted perimeter
            _da = ((ws_elev - zmin)/(zmax - zmin)) * xy
            _db = ws_elev - zmin

            wp_vec[idx] = np.sqrt(_da**2.0 + _db**2.0)

        elif loc == 'trapezoid':

            zmin = zmin_by_segment[idx]
            zmax = zmax_by_segment[idx]
            xy = xydist[idx]

            area_vec[idx] = 0.5 * xy * (2*ws_elev - zmax - zmin)

            wp_vec[idx] = np.sqrt(xy**2.0 + (zmax - zmin)**2.0)

    area_sum = sum(area_vec)
    wp_sum = sum(wp_vec)

    n_inv = (1.0/streambed_roughness)

    Q = n_inv * area_sum * (pow((area_sum/wp_sum), 2/3.0)) * np.sqrt(slope)

    return Q


def _get_ws_location(water_surface_elev, zmax, zmin):
    """
    Return one of three values depending on the location of the water surface
    relative to the elevations of the discretized cross-section points.
    Vectorized below.

    Returns:
        (str) one of the following: 'below' if above zmax (and automatically
            zmin), 'triangle' if between zmax and zmin, and 'trapezoid'
            if below zmin. This corresponds to the geometry that will be used
            to calculate the wetted perimeter and area of the induced polygon.
    """

    if water_surface_elev > zmax:
        return 'trapezoid'
    elif water_surface_elev <= zmax and water_surface_elev > zmin:
        return 'triangle'
    else:
        return 'below'


class BoundaryCondition:

    def __init__(self,
                 period=0.0,  # (minutes)
                 amplitude=0.0,  # (ISO)
                 phase=0.0):  # (deg)

        self.period = period
        self.amplitude = amplitude
        self.phase = phase

    def write(self, out_path):

        with open(out_path, 'w') as f:
            f.write(self.__repr__())

    def __repr__(self):
        return '\n'.join([
                '* COLUMN=3',
                '* COLUMN1=Period (min) or Astronomical Componentname',
                '* COLUMN2=Amplitude (ISO)',
                '* COLUMN3=Phase (deg)',
                '{0}  {1}  {2}'.format(self.period, self.amplitude, self.phase)
            ])


def mr_log(log_f, msg):

    ta = time.asctime
    log_f.write('[{0}] '.format(ta()) + msg)
    log_f.flush()
    os.fsync(log_f.fileno())

def mr_progress_add_entry(progressfile,  flow_idx):
    progressfile.write(str(flow_idx) + '\t0\t0\n')
    progressfile.flush()
    os.fsync(progressfile.fileno())


def mr_progress_update_entry(progressfilepath,  flow_idx, dflow_status, ripcas_status):
    lines = open(progressfilepath).read().splitlines()
    lines[flow_idx + 1] = str(flow_idx) + '\t' + str(dflow_status) +'\t' + str(ripcas_status) + '\n'
    open(progressfilepath,'w').write('\n'.join(lines))


def modelrun_series(data_dir, initial_vegetation_map, vegzone_map,
                    veg_roughness_shearres_lookup, peak_flows_file,
                    geometry_file, streambed_roughness,
                    streambed_floodplain_roughness, streambed_slope,
                    dflow_run_fun=None, log_f=None, progressfilepath='cord_progress.log', debug=False):
    '''
    Run a series of flow and succession models with peak flows given in
    peak_flows_file.

    Arguments:
        data_dir (str): write directory for modelrun series. Must exist
        initial_vegetation_map (str): location of year zero veg map
        vegzone_map (str): vegetation zone map location
        veg_roughness_shearres_lookup (str): Excel spreadsheet containing
            conversion from vegetation code to roughness value and vegetation
            code to shear stress resistance
        peak_flow_file (str): location of text file record of peak flows in
            cubic meters per second
        geometry_file (str): location of channel geometry at the downstream
            location for calculating streamflow
        streambed_roughness (float): streambed roughness in channel only; used
            when converting vegetation map to roughness map
        streambed_floodplain_roughness (float): an average roughness of
            stream channel and floodplain used in calculation of downstream
            boundary condition for DFLOW
        streambed_slope (float): rise over run of the channel used in
            calculation of downstream boundary condition for DFLOW
        dflow_run_fun (function): function delegate for the user to provide a
            custom way to run DFLOW. If none is given, defaults to
            submitting a PBS job as is done on CARC systems
        log_f (str): log file. if none is given, defaults to `data_dir`.log
            with dashes replacing slashes
        debug (bool): whether or not to run in debug mode. If running in debug
            mode, each DFLOW run returns fake data and
            each RipCAS run takes cord/data/shear_out.asc as input

        returns:
            None
    '''
    # If standard run on CARC
    if dflow_run_fun is None:

        def dflow_fun():

            import subprocess
            # Send DFLOW run to the queue
            return subprocess.Popen(
                'qsub dflow_mpi.pbs', shell=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )

    # Create log file is none exists
    if log_f is None:
        # get rid of / at begining of file path
        first_char = data_dir[0]
        root_log_f = first_char if first_char != '/' else ''
        root_log_f += data_dir[1:].replace('/', '-')

        log_f = open(root_log_f + '.log', 'w')

    else:
        log_f = open(log_f, 'w')
    
    # Create the progress file if none exists - append the header for the flow_idx were on, dflow, and ripcas columns to indicate if they finished
    progressfile = open(progressfilepath, 'w')
    progressfile.write('flow_idx\tdflow_completed\tripcas_completed\n')

    # create a list that contains the peak flows from input file
    with open(peak_flows_file, 'r') as f:
        l0 = f.readline().strip()
        assert l0 == 'Peak.Flood', '{} not Peak.Flood'.format(l0)
        peak_flows = [float(l.strip()) for l in f.readlines()]

    # create a directory for global inputs
    inputs_dir = os.path.join(data_dir, 'inputs')
    # remove inputs directory if it already existed
    if os.path.isdir(inputs_dir):
        shutil.rmtree(inputs_dir)
    # create inputs directory
    os.mkdir(inputs_dir)
    # bring all input files into the input directory
    shutil.copy(initial_vegetation_map, inputs_dir)
    shutil.copy(vegzone_map, inputs_dir)
    shutil.copy(veg_roughness_shearres_lookup, inputs_dir)
    shutil.copy(peak_flows_file, inputs_dir)
    shutil.copy(geometry_file, inputs_dir)

    roughness_slope_path = os.path.join(inputs_dir, 'roughness_slope.txt')

    # create a text file with info on both streambed roughness and slope
    with open(roughness_slope_path, 'w') as f:
        f.write('roughness\tslope\n')
        f.write('%s\t%s\n' % (streambed_roughness, streambed_slope))

    # Iterate through all annual peak flows
    for flow_idx, flow in enumerate(peak_flows):
        #create an entry in the progress file for this flow_idx
        mr_progress_add_entry(progressfile, flow_idx)
        # create a ModelRun object
        mr = ModelRun()
        # Run the boundary condition calculation method;
        # produces upstream flow file and downstream stage file for DFLOW to use
        mr.calculate_bc(
            flow, geometry_file,
            streambed_floodplain_roughness, streambed_slope
        )

        # Enter information into log file
        mr_log(
            log_f, 'Boundary conditions for flow index {0} finished\n'.format(
                flow_idx
            )
        )

        # Create new directory for this annual flow iteration of DFLOW
        dflow_dir = os.path.join(data_dir, 'dflow-' + str(flow_idx))

        # Get veg map
        if flow_idx == 0:
            veg_file = initial_vegetation_map
        else:
            # Take RipCAS outputs as DFLOW inputs from previous timestep
            veg_file = os.path.join(
                data_dir, 'ripcas-' + str(flow_idx - 1), 'vegetation.asc'
            )

        # Debug is for running on a local machine
        if debug:
            mr.run_dflow(dflow_dir, veg_file,
                         veg_roughness_shearres_lookup, streambed_roughness)
            job_id = 'debug'

        # If running on CARC
        else:
            # Send DFLOW run to CARC
            p_ref = mr.run_dflow(dflow_dir, veg_file,
                                 veg_roughness_shearres_lookup,
                                 streambed_roughness,
                                 dflow_run_fun=dflow_fun)

            job_id = p_ref.communicate()[0].split('.')[0]
        # Enter run start in log file
        mr_log(log_f, 'Job ID {0} submitted for DFLOW run {1}\n'.format(
                job_id, flow_idx
            )
        )

        # check the status of the job by querying qstat; break loop when
        # job no longer exists, giving nonzero poll() value
        job_not_finished = True
        while job_not_finished:

            mr_log(
                log_f,
                'Job ID {0} not yet finished for DFLOW run {1}\n'.format(
                    job_id, flow_idx
                )
            )

            if debug:
                job_not_finished = False

            else:
                p = subprocess.Popen(
                    'qstat ' + job_id, shell=True,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )

                p.communicate()

                poll = p.poll()
                job_not_finished = poll == 0

                time.sleep(600)
        mr_log(
            log_f, 'DFLOW run {0} finished, starting RipCAS\n'.format(
                flow_idx
            )
        )
        
        mr_progress_update_entry(progressfilepath, flow_idx, 1, 0)

        # Creat a directory for this annual iteration of RipCAS
        ripcas_dir = os.path.join(data_dir, 'ripcas-' + str(flow_idx))

        # Debug is for running on a local machine
        if debug:
            p = _join_data_dir('shear_out.asc')
            mr.run_ripcas(vegzone_map, veg_roughness_shearres_lookup,
                          ripcas_dir, shear_asc=ESRIAsc(p))
            # mr.run_ripcas(vegzone_map, veg_roughness_shearres_lookup,
            #               ripcas_dir)

        else:
            # if no explicit shear_asc is given, the method accesses
            # the dflow_shear_output attribute. XXX TODO this method will
            # need to be updated to build a shear_asc by stitching together
            # the partitioned files using stitch_partitioned_output
            # in cord/ripcas_dflow.py
            mr.run_ripcas(vegzone_map, veg_roughness_shearres_lookup,
                          ripcas_dir)
        # Note end of RipCAS in log file
        mr_log(log_f, 'RipCAS run {0} finished\n'.format(flow_idx))
        mr_progress_update_entry(progressfilepath, flow_idx, 1, 1)

    log_f.close()


def _join_data_dir(f):
    '''
    Join the filename, f, to the default data directory
    '''
    data_dir = os.path.join(os.path.dirname(__file__), 'data')

    return os.path.join(data_dir, f)
