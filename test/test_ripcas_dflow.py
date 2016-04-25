"""
Tests for posting/ingetsting dflow and ripcas data to each respective model
to/from the virtual watershed.
"""
import os
import numpy
import shutil
import unittest

from ripcas_dflow import ESRIAsc, Pol, ripcas, veg2n, ModelRun


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

    # def test_ripcas_with_dflow_io(self):
        # assert False

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
        geometry = 'data/DBC_geometry.xyz'
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
        self.mr.run_dflow(os.path.join(self.tmpdir, 'dflow-test'),
                          'test/data/vegcode.asc')

        assert os.path.exists(os.path.join(self.tmpdir, 'dflow-test'))

        # TODO add more checks for files that should be there.. translate
        # dflow checklist from Angela
        assert self.mr.dflow_has_run

    def test_run_ripcas(self):
        """
        RipCAS output should match expected
        """
        # no explicit private properties/methods in python, allows this hack
        self.mr.dflow_has_run = True
        self.mr.vegetation_ascii = ESRIAsc('test/data/vegcode.asc')

        out = self.mr.run_ripcas(
            'test/data/zonemap.asc', 'test/data/resist_manning_lookup.xlsx',
            os.path.join(self.tmpdir, 'ripcas-test'),
            shear_asc=ESRIAsc('test/data/shear.asc')
        )

        assert out == ESRIAsc('test/data/expected_veg_output.asc')
