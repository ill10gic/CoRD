"""
Tests for posting/ingetsting dflow and ripcas data to each respective model
to/from the virtual watershed.
"""
import os
import numpy
import re
import responses
import shutil
import unittest

from click.testing import CliRunner

from ripcas_dflow import ESRIAsc, Pol, ripcas, veg2n, ModelRun
from ripcas_dflow.scripts.cord import cli


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
        vmat_expected = numpy.array([-9999, 100, 101, 102, 106, 210, 215],
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

        nmap = veg2n(veg_map, self.ripcas_required_data)

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

    ### TODO
    # def test_mesh_to_asc(self):
        # assert False


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
        # geometry = Pol.from_river_geometry_file('data/DBC_geometry.xyz')
        geometry = 'data/dflow_inputs/DBC_geometry.xyz'
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
        self.mr.run_dflow(d, 'test/data/vegcode.asc')

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

        self.mr.dflow_directory = 'test/data'
        out = self.mr.run_ripcas(
            'test/data/zonemap.asc', 'test/data/resist_manning_lookup.xlsx',
            os.path.join(self.tmpdir, 'ripcas-test'),
            shear_asc=ESRIAsc('test/data/shear.asc')
        )

        assert out == ESRIAsc('test/data/expected_veg_output.asc')

    def test_modelrun_series(self):
        assert False


class TestCLI(unittest.TestCase):
    """
    Test `cord` CLI
    """

    def setUp(self):
        self.r = CliRunner()

        self.hs_basedir = 'fakedir'
        set_up_for_hs_test(self.hs_basedir)

    def tearDown(self):

        shutil.rmtree(self.hs_basedir)

    @responses.activate
    def test_push_hs(self):

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
            rsps.add(responses.POST, base_url + 'resource/' + rid + '/files/',
                     json={'resource_id': rid},
                     status=201)

            runner = CliRunner()

            runner.invoke(
                cli, ['post_hs', '--username', 'fake', '--password', 'fake',
                      '--modelrun-dir', self.hs_basedir, '--resource-title',
                      'fake resource that will never get to HydroShare']
            )

            assert len(rsps.calls) == 5

    def test_from_config(self):
        assert False

    def test_interactive(self):
        assert False


def set_up_for_hs_test(basedir):

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
