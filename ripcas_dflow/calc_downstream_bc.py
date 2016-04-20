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

from collections import namedtuple

from ripcas_dflow import Pol


class CoupledRun(object):

    def __init__(self):

        # have the boundary conditions been found?
        self.bc_converged = False
        # has ripcas been run yet?
        self.ripcas_has_run = False
        # has dflow?
        self.dflow_has_run = False

        self.process_state = 'not started'

        self.upstream_bc = BoundaryCondition()
        self.downstream_bc = BoundaryCondition()

        self.downstream_bc_ws_step = 0.001

        self.streambed_roughness = 0.04

        self.slope = 0.001

    def calculate_bc(self, dbc_geometry_file):

        try:
            dbc_geometry = Pol.from_river_geometry_file(dbc_geometry_file)

            Q_seq = _streamflow_sequence(
                dbc_geometry, self.streambed_roughness, self.slope, self.slope
            )

            # TODO need to read from file somehow...
            historical_streamflow = 20.5

            streamflow_tuple_opt =\
                _find_minimizing_water_surface(Q_seq, historical_streamflow)

            self.bc_converged = True

            self.upstream_bc.amplitude = streamflow_tuple_opt.ws_elev

            # XXX I don't think this is right...
            self.downstream_bc.amplitude = streamflow_tuple_opt.streamflow

            return

        except:

            self.process_state = 'boundary condition calculation failed'

            return

    def run_dflow(self):

        if not self.bc_converged:
            raise RuntimeError(
                'Boundary conditions must be calculated before ' +
                'DFLOW can be run'
            )

        if self.dflow_has_run:
            raise RuntimeError(
                'DFLOW has already been run for this CoupledRun'
            )

        self.dflow_has_run = True

        return

    def run_ripcas(self):

        if not self.dflow_has_run:
            raise RuntimeError(
                'DFLOW must run before ripcas can be run'
            )

        self.ripcas_has_run = True

        return


def _find_minimizing_water_surface(streamflow_sequence, historic_streamflow):

    Q_hist = historic_streamflow
    err_by_elev = (
                      (Q.ws_elev, abs(Q.streamflow - Q_hist))
                      for Q in streamflow_sequence
                  )

    optimal_ws_elev_and_Q = min(err_by_elev, key=lambda a: a[1])

    return optimal_ws_elev_and_Q


def _streamflow_sequence(dbc_geometry, streambed_roughness, slope,
                         water_surface_step):
    """
    Generator to calculate iterate over
    """
    # TODO get these from dbc_geometry
    water_surface_elev_min = 0.5

    water_surface_elev = water_surface_elev_min

    water_surface_elev_max = 20.0

    while water_surface_elev < water_surface_elev_max:

        yield StreamflowTuple(
            water_surface_elev,
            _calculate_streamflow(
                dbc_geometry, streambed_roughness, water_surface_elev
            )
        )

        water_surface_elev += water_surface_step


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


# _get_ws_location = np.vectorize(_get_ws_location_one)


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
            f.writelines([
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
