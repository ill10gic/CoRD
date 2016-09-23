"""
Tests for posting/ingetsting dflow and ripcas data to each respective model
to/from the virtual watershed.
"""
import glob
import numpy
import os
import responses
import shutil
import traceback
import unittest
import numpy as np

from click.testing import CliRunner
from netCDF4 import Dataset

from cord import (
    ESRIAsc, Pol, ripcas, veg2n, ModelRun,
    stitch_partitioned_output
)
from cord.scripts.cord import cli


class TestRipCASAndHelpers(unittest.TestCase):
    """
    Functions for working with DFLOW inputs and outputs
    """
    def setUp(self):
        self.ascii_veg = 'test/data/vegcode.asc'
        self.ripcas_required_data = 'test/data/resist_manning_lookup.xlsx'
        self.expected_ascii_roughness = \
            'test/data/roughness.asc'

        self.expected_ascii_nvals = \
            ESRIAsc(self.expected_ascii_roughness)

    def test_vegmap_properly_read(self):

        vegmap_mat = ESRIAsc(self.ascii_veg).as_matrix()

        vmat_unique = numpy.unique(vegmap_mat)
        vmat_expected = numpy.array([-9999, 0, 100, 101, 102, 106, 210, 215],
                                    dtype='f8')

        assert (vmat_unique == vmat_expected).all()

    def test_ripcas(self):

        # load the expected ESRIAsc output from running ripcas
        expected_output = ESRIAsc(
            'test/data/expected_veg_output.asc'
        )

        # test results when loaded from file
        veg_map_file = self.ascii_veg
        shear_map_file = 'test/data/shear.asc'
        zone_map_file = 'test/data/zonemap.asc'

        generated_output = ripcas(veg_map_file, zone_map_file, shear_map_file,
                                  self.ripcas_required_data)

        assert expected_output == generated_output, \
            "expected: {}\ngenerated: {}".format(
                expected_output.as_matrix(), generated_output.as_matrix()
            )

        # test results when using ESRIAsc instances
        veg_map = ESRIAsc(veg_map_file)
        zone_map = ESRIAsc(zone_map_file)
        shear_map = ESRIAsc(shear_map_file)

        generated_output = ripcas(veg_map, zone_map,
                                  shear_map, self.ripcas_required_data)

        assert expected_output == generated_output, \
            "expected: {}\ngenerated: {}".format(
                expected_output.as_matrix(), generated_output.as_matrix()
            )

    def test_veg2n(self):
        """
        Test conversion of vegetation map to Manning's roughness map
        """
        expected_nmap = ESRIAsc(
            'test/data/expected_nmap.asc'
        )

        veg_map = ESRIAsc(
            'test/data/vegcode.asc'
        )

        nmap = veg2n(veg_map, self.ripcas_required_data, 0.035)

        assert nmap == expected_nmap, \
            "nmap: {}\nexpected_nmap: {}".format(nmap.data, expected_nmap.data)

    def test_asc2pol(self):
        """
        asc2pol should create proper headers and formatted data
        """
        expected_pol = Pol.from_dflow_file('test/data/expected_n.pol')

        nmap = ESRIAsc('test/data/expected_nmap.asc')

        npol = Pol.from_ascii(nmap)

        assert npol == expected_pol


class TestModelRun(unittest.TestCase):
    """

    """
    def setUp(self):

        self.mr = ModelRun()

        self.tmpdir = 'test/data/tmp'
        if os.path.exists(self.tmpdir):
            shutil.rmtree(self.tmpdir)

        os.mkdir(self.tmpdir)

    def tearDown(self):

        shutil.rmtree(self.tmpdir)

    def test_boundary_calculation(self):
        """
        Calculate a series of boundary conditions for the range we'll see
        """
        geometry = 'cord/data/dflow_inputs/DBC_geometry.xyz'
        roughness = 0.04
        reach_slope = 0.001

        peak_flows = range(12, 120, 2)
        for p in peak_flows:
            self.mr.calculate_bc(p, geometry, roughness, reach_slope)
            assert self.mr.bc_converged, 'failed for peak flow {}'.format(p)
            self.mr.bc_converged = False

        self.mr.bc_converged = True

    def test_run_dflow(self):
        """
        DFLOW create the proper directory and populate with required inputs
        """
        self.mr.bc_converged = True
        d = os.path.join(self.tmpdir, 'dflow-test')
        self.mr.run_dflow(d, 'test/data/vegcode.asc',
                          'test/data/resist_manning_lookup.xlsx',
                          0.35)

        assert os.path.exists(d)
        # there are six required files that should have been copied to the
        # dflow directory. vegcode.asc should have been translated
        # to n.pol and written to this directory. The two boundary
        # condition files should also have been written to the
        # directory

        def ex(f):
            return os.path.exists(os.path.join(d, f))

        assert ex('dflow_mpi.pbs')
        assert ex('base.mdu')
        assert ex('base_net.nc')
        assert ex('base.ext')
        assert ex('boundriverdown.pli')
        assert ex('boundriver_up.pli')
        assert ex('boundriverdown_0001.cmp')
        assert ex('boundriver_up_0001.cmp')

        assert ex('n.pol')

        assert self.mr.dflow_has_run

    def test_run_ripcas(self):
        """
        RipCAS output should match expected
        """
        # no explicit private properties/methods in python, allows this hack
        self.mr.dflow_has_run = True
        self.mr.vegetation_ascii = ESRIAsc('test/data/vegcode.asc')

        self.mr.dflow_run_directory = 'test/data'
        out = self.mr.run_ripcas(
            'test/data/zonemap.asc', 'test/data/resist_manning_lookup.xlsx',
            os.path.join(self.tmpdir, 'ripcas-test'),
            shear_asc=ESRIAsc('test/data/shear.asc')
        )

        assert out == ESRIAsc('test/data/expected_veg_output.asc')


class TestCLI(unittest.TestCase):
    """
    Test `cord` CLI
    """

    def setUp(self):
        pass

    def tearDown(self):
        pass

    @responses.activate
    def test_push_hs(self):

        hs_basedir = 'fakedir'
        _set_up_for_hs_test(hs_basedir)

        base_url = 'https://www.hydroshare.org/hsapi/'
        rid = 'X5A67'

        # if all mock rsps are not used w/in context it raises AssertError
        with responses.RequestsMock() as rsps:

            # set up response objects
            # response to create new resource
            rsps.add(responses.POST, base_url + 'resource/',
                     json={'resource_id': rid},
                     status=201)

            # before resource created, client checks types; auth req'd
            rsps.add(responses.GET, base_url + 'resourceTypes/',
                     json=[{'resource_type': 'GenericResource'}],
                     status=200)

            # response to adding vegetation, inputs, and  file
            file_add_url = base_url + 'resource/' + rid + '/files/'

            for i in range(3):
                rsps.add(responses.POST, file_add_url,
                         json={'resource_id': rid},
                         status=201)

            # prepare cli run
            runner = CliRunner()

            result = runner.invoke(
                cli, ['post_hs', '--username', 'fake', '--password', 'fake',
                      '--modelrun-dir', hs_basedir, '--resource-title',
                      'fake resource that will never get to HydroShare']
            )
            assert not result.exception, result.exception
            assert result.exit_code == 0

            assert len(rsps.calls) == 5

        shutil.rmtree(hs_basedir)

    def test_from_config(self):

        runner = CliRunner()

        mr_dir = 'test/data/modelrun'
        if os.path.isdir(mr_dir):
            shutil.rmtree(mr_dir)

        os.mkdir(mr_dir)

        result = runner.invoke(
            cli, ['--debug', 'from_config', 'test/data/test.conf']
        )

        assert not result.exception,\
            str(result.exception) + '\n' + \
            str(traceback.print_tb(result.exc_info[2]))

        assert result.exit_code == 0, result.exit_code

        _test_modelrun_success(mr_dir)

        shutil.rmtree('test/data/modelrun')

    def test_interactive(self):

        runner = CliRunner()

        mr_dir = 'test/data/modelrun'

        if os.path.isdir(mr_dir):
            shutil.rmtree(mr_dir)

        os.mkdir(mr_dir)

        result = runner.invoke(
            cli, ['--debug', 'interactive', '--data-dir', mr_dir,
                  '--initial-veg-map', 'test/data/vegcode.asc',
                  '--vegzone-map', 'test/data/zonemap.asc',
                  '--ripcas-required-data',
                  'test/data/resist_manning_lookup.xlsx',
                  '--peak-flows-file', 'test/data/floods.txt',
                  '--geometry-file', 'cord/data/dflow_inputs/DBC_geometry.xyz',
                  '--streambed-roughness', '0.04',
                  '--streambed-floodplain-roughness', '0.035',
                  '--streambed-slope', '0.001']
        )

        assert not result.exception,\
            str(result.exception) + '\n' + \
            str(traceback.print_tb(result.exc_info[2]))

        assert result.exit_code == 0, result.exit_code

        _test_modelrun_success(mr_dir)

        shutil.rmtree('test/data/modelrun')


class TestMeshToAsc(unittest.TestCase):
    """

    """
    def setUp(self):
        self.nc_names = [
            os.path.join(
                'test', 'data', 'tmp',
                'jempatmdu_000{}_map.nc'.format(i)
            )
            for i in range(4)
        ]

        os.mkdir('test/data/tmp')

        # on inspection found all the relevant vars are not also ds
        for idx, name in enumerate(self.nc_names):

            nc = Dataset(name, 'w')

            nc.createDimension('nFlowElem', 8)
            nc.createVariable('FlowElem_xcc', float, ('nFlowElem',))
            nc.createVariable('FlowElem_ycc', float, ('nFlowElem',))
            nc.createVariable('FlowElemDomain', float, ('nFlowElem',))
            nc.createVariable('taus', float, ('nFlowElem',))

            vars_ = nc.variables
            # TODO fill in variables
            if idx == 0:
                vars_['FlowElemDomain'][:] = [0, 0, 0, 0, 1, 2, 2, 3]
                vars_['FlowElem_xcc'][:] =\
                    [32.3, 34.5, 36.7, 38.1, 48.3, 33.3, 2.1, 1.0]
                vars_['FlowElem_ycc'][:] =\
                    [132.3, 134.5, 136.7, 138.1, 8.3, 133.3, 12.1, 111.0]

            elif idx == 1:
                vars_['FlowElemDomain'][:] = [0, 1, 1, 1, 1, 2, 2, 3]
                vars_['FlowElem_xcc'][:] =\
                    [32.3, 39.255, 40.25, 40.55, 41.1, 123.3, 4.1, 5.0]
                vars_['FlowElem_ycc'][:] =\
                    [232.3, 139.255, 140.25, 140.55, 141.1, 138.1, 8.3, 133.3]

            elif idx == 2:
                vars_['FlowElemDomain'][:] = [0, 1, 2, 2, 2, 2, 3, 3]
                vars_['FlowElem_xcc'][:] =\
                    [32.3, 29.255, 42.2, 43.01, 44.2, 45.5, 14.1, 52.1]
                vars_['FlowElem_ycc'][:] =\
                    [232.3, 139.255, 142.2, 143.01, 144.2, 145.5, 8.3, 133.3]

            elif idx == 3:
                vars_['FlowElemDomain'][:] = [1, 1, 1, 2, 3, 3, 3, 3]
                vars_['FlowElem_xcc'][:] =\
                    [32.3, 29.255, 42.2, 43.01, 49.2, 50.7, 51.6, 52.2]
                vars_['FlowElem_ycc'][:] =\
                    [232.3, 139.255, 142.2, 143.01, 149.2, 150.7, 151.6, 152.2]

            vars_['taus'][:] = np.arange(8) * (idx + 1)

    def tearDown(self):

        os.remove(os.path.join('test', 'data', 'tmp', 'stitched.nc'))

    def test_meshes_to_asc(self):

        mesh_nc_paths = glob.glob(
            os.path.join('test', 'data', 'tmp', '*')
        )

        # should return a netCDF
        stitched = Dataset(
            stitch_partitioned_output(
                mesh_nc_paths, 'test/data/tmp/stitched.nc'
            )
        )

        xcc_stitched = stitched.variables['FlowElem_xcc'][:]
        ycc_stitched = stitched.variables['FlowElem_ycc'][:]

        xcc_expected = np.array(
            [32.3, 34.5, 36.7, 38.1,
             39.255, 40.25, 40.55, 41.1,
             42.2, 43.01, 44.2, 45.5,
             49.2, 50.7, 51.6, 52.2
             ]
        )
        ycc_expected = np.array(
            [132.3, 134.5, 136.7, 138.1,
             139.255, 140.25, 140.55, 141.1,
             142.2, 143.01, 144.2, 145.5,
             149.2, 150.7, 151.6, 152.2
             ]
        )

        assert (xcc_expected == xcc_stitched).all()
        assert (ycc_expected == ycc_stitched).all(), \
            "\nexpected: {}\nstitched: {}".format(ycc_expected, ycc_stitched)

        taus_expected = np.array(
            [0.0, 1.0, 2.0, 3.0,
             2.0, 4.0, 6.0, 8.0,
             6.0, 9.0, 12.0, 15.0,
             16.0, 20.0, 24.0, 28.0
             ]
        )

        taus_stitched = stitched.variables['taus'][:]

        assert (taus_expected == taus_stitched).all()


def _set_up_for_hs_test(basedir):

    if os.path.isdir(basedir):
        shutil.rmtree(basedir)

    os.mkdir(basedir)

    opj = os.path.join

    input_dir = opj(basedir, 'inputs')
    os.mkdir(input_dir)

    inputs = ['geom.txt', 'flows.txt', 'ripcas-required.xlsx',
              'vegzone.asc', 'init_veg.asc', 'roughness_slope.txt']

    for i in inputs:
        with open(opj(input_dir, i), 'w') as f:
            f.write('fake!')

    ripcas_dirs = [opj(basedir, 'ripcas-0'),
                   opj(basedir, 'ripcas-1')]

    dflow_dirs = [opj(basedir, 'dflow-0'),
                  opj(basedir, 'dflow-1')]

    for d in (dflow_dirs + ripcas_dirs):

        os.mkdir(d)

        if 'dflow' in d:
            with open(opj(d, 'shear_out.asc'), 'w') as f:
                f.write('fake shear out!')

        if 'ripcas' in d:
            with open(opj(d, 'vegetation.asc'), 'w') as f:
                f.write('fake vegetation!')


def _test_modelrun_success(modelrun_dir):

    opj = os.path.join
    mrd = modelrun_dir

    g = glob.glob(opj(mrd, '*'))

    assert len(g) == 5

    df_dirs = [opj(mrd, 'dflow-0'), opj(mrd, 'dflow-1')]

    for d in df_dirs:

        assert d in g

        g_df = glob.glob(opj(d, '*'))
        assert opj(d, 'base.ext') in g_df
        assert opj(d, 'base.mdu') in g_df
        assert opj(d, 'base_net.nc') in g_df
        assert opj(d, 'boundriver_up.pli') in g_df
        assert opj(d, 'boundriverdown.pli') in g_df
        assert opj(d, 'boundriver_up_0001.cmp') in g_df
        assert opj(d, 'boundriverdown_0001.cmp') in g_df
        assert opj(d, 'dflow_mpi.pbs') in g_df
        assert opj(d, 'n.pol') in g_df
        assert opj(d, 'shear_out.asc') in g_df

    ripcas_dirs = [opj(mrd, 'ripcas-0'), opj(mrd, 'ripcas-1')]

    for d in ripcas_dirs:
        assert d in g
        assert os.path.exists(opj(d, 'vegetation.asc'))

    inputs_dir = opj(mrd, 'inputs')
    inputs_list = [
        'DBC_geometry.xyz', 'floods.txt', 'resist_manning_lookup.xlsx',
        'roughness_slope.txt', 'vegcode.asc', 'zonemap.asc'
    ]

    for i in inputs_list:
        assert os.path.exists(opj(inputs_dir, i))

    with open(opj(inputs_dir, 'roughness_slope.txt'), 'r') as f:
        assert f.read() == 'roughness\tslope\n0.04\t0.001\n'
