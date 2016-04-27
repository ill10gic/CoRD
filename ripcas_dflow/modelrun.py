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
import numpy as np
import os
import shutil

from collections import namedtuple
from scipy.optimize import minimize_scalar

from ripcas_dflow import ESRIAsc, Pol, ripcas, shear_mesh_to_asc, veg2n


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

        self.dflow_has_run = False
        self.dflow_directory = None
        self.dflow_shear_output = None

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

    def run_dflow(self, dflow_directory, vegetation_map, clobber=True,
                  pbs_script_name='dflow_mpi.pbs', dflow_run_fun=None):
        """
        Both input and output dflow files will go into the dflow_directory,
        but in input/ and output/ subdirectories.

        Arguments:
            dflow_directory (str): directory where DFLOW files should be
                put and where the dflow_run_fun will be run from
            vegetation_map (str): path to the input vegetation.pol file. This
                function assumes this has already been generated in the proper
                format b/c this seems like the best separation of
                responsibilities.
            clobber (bool): whether or not to overwrite dflow_directory if
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

        if os.path.exists(dflow_directory):

            if not clobber:
                raise RuntimeError(
                    'DFLOW has already been run for this CoupledRun'
                )

            shutil.rmtree(dflow_directory)

        self.dflow_directory = dflow_directory

        os.mkdir(dflow_directory)

        # inputs_path = os.path.join(dflow_directory, 'inputs')
        # outputs_path = os.path.join(dflow_directory, 'outputs')

        # os.mkdir(inputs_path)
        # os.mkdir(outputs_path)

        # write boundary conditions to file
        bc_up_path = os.path.join(dflow_directory, 'boundriver_up_0001.cmp')
        bc_down_path = os.path.join(dflow_directory, 'boundriverdown_0001.cmp')

        self.upstream_bc.write(bc_up_path)
        self.downstream_bc.write(bc_down_path)

        self.vegetation_ascii = ESRIAsc(vegetation_map)

        veg_path = os.path.join(dflow_directory, 'n.pol')

        Pol.from_ascii(
            veg2n(self.vegetation_ascii,
                  'data/ripcas_inputs/ripcas-data-requirements.xlsx')
        ).write(veg_path)

        oj = os.path.join

        pbs_path = oj(dflow_directory, pbs_script_name)
        mdu_path = oj(dflow_directory, 'base.mdu')
        net_path = oj(dflow_directory, 'base_net.nc')
        ext_path = oj(dflow_directory, 'base.ext')
        brd_path = oj(dflow_directory, 'boundriverdown.pli')
        bru_path = oj(dflow_directory, 'boundriver_up.pli')

        self.dflow_shear_output =\
            os.path.join(dflow_directory, 'jemez_r02_map.nc')

        with open(pbs_path, 'w') as f:
            s = open('data/dflow_inputs/dflow_mpi.pbs', 'r').read()
            f.write(s)

        with open(mdu_path, 'w') as f:
            s = open('data/dflow_inputs/base.mdu', 'r').read()
            f.write(s)

        with open(net_path, 'w') as f:
            s = open('data/dflow_inputs/base_net.nc', 'r').read()
            f.write(s)

        with open(ext_path, 'w') as f:
            s = open('data/dflow_inputs/base.ext', 'r').read()
            f.write(s)

        with open(brd_path, 'w') as f:
            s = open('data/dflow_inputs/boundriverdown.pli', 'r').read()
            f.write(s)

        with open(bru_path, 'w') as f:
            s = open('data/dflow_inputs/boundriver_up.pli', 'r').read()
            f.write(s)

        bkdir = os.getcwd()
        os.chdir(dflow_directory)

        if dflow_run_fun is None:

            print '\n*****\nDry Run of DFLOW\n*****\n'

            os.chdir(bkdir)
            example_shear_path = 'jemez_r02_map.nc'
            if os.path.exists(example_shear_path):
                shutil.copyfile(example_shear_path, self.dflow_shear_output)

            else:
                print 'Get you a copy of a DFLOW output, yo! ' +\
                      'Can\'t run RipCAS without it!'

                with open('not_actually_output.nc', 'w') as f:
                    f.write('A FAKE NETCDF!!!')
        else:

            dflow_run_fun()
            os.chdir(bkdir)

        self.dflow_has_run = True

        return

    def run_ripcas(self, zone_map_path, ripcas_required_data_path,
                   ripcas_directory, shear_asc=None, clobber=True):

        if not self.dflow_has_run:
            raise RuntimeError(
                'DFLOW must run before ripcas can be run'
            )

        if os.path.exists(ripcas_directory):

            if not clobber:
                raise RuntimeError(
                    'DFLOW has already been run for this CoupledRun'
                )

            shutil.rmtree(ripcas_directory)

        self.ripcas_directory = ripcas_directory

        os.mkdir(ripcas_directory)

        hdr = self.vegetation_ascii.header_dict()

        if shear_asc is None:
            shear_asc = shear_mesh_to_asc(self.dflow_shear_output, hdr)
        else:
            assert isinstance(shear_asc, ESRIAsc),\
                'shear_asc must be of type ESRIAsc if provided'

        output_veg_ascii = ripcas(
            self.vegetation_ascii, zone_map_path,
            shear_asc, ripcas_required_data_path
        )

        output_vegetation_path = os.path.join(
            ripcas_directory, 'vegetation.asc'
        )
        output_veg_ascii.write(output_vegetation_path)

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


class ExtFile:

    def __init__(self):
        pass

    def write(self, out_path):
        pass
