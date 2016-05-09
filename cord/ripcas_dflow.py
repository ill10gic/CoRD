"""
Utilities for interacting with dflow and ripcas models via the virtual
watershed.
"""
import copy
import six

from netCDF4 import Dataset
from numpy import (array, concatenate, fromstring, meshgrid,
                   reshape, flipud, isnan)
from scipy.interpolate import griddata
from pandas import read_table, read_excel, Series


def ripcas_with_dflow_io(vegetation_map, zone_map,
                         shear_nc_path, ripcas_required_data):
    """
    Wrapper for using DFLOW input/output with ripcas. Note instead of
    shear_map we have shear_nc_path. Use shear_mesh_to_asc to convert the
    shear_mesh that comes from D-FLOW to a shear_map for input to ripcas.
    When ripcas finishes its vegetation updates, convert the updated
    vegetation map to a Manning n-value map for use by DFLOW. See the
    ripcas function below for more details on the model and arguments.

    Arguments:
        vegetation_map (ESRIAsc): location on disk or ESRIAsc
            representation of the vegetation map
        zone_map (ESRIAsc): location on disk or ESRIAsc representation
            of the zone map.
        shear_nc_path (str): location on disk of DFLOW shear output netCDF
        ripcas_required_data (str): Excel spreadsheet of data needed for
            ripcas run. Encompasses the landscape model for the watershed.
            Can have one or two 'Code' columns and must have exactly one
            'shear_resis' column and exactly one 'n_val' column

    Returns:
        (Pol): polygon representation of the map of nvalues
    """
    if not isinstance(vegetation_map, ESRIAsc):
        raise TypeError('vegetation_map must be an ESRIAsc instance')
    if not isinstance(zone_map, ESRIAsc):
        raise TypeError('vegetation_map must be an ESRIAsc instance')

    shear_map = shear_mesh_to_asc(shear_nc_path, vegetation_map.header_dict())

    return Pol.from_ascii(
        veg2n(
            ripcas(
                vegetation_map, zone_map, shear_map, ripcas_required_data
            ), ripcas_required_data
        )
    )


def ripcas(vegetation_map, zone_map, shear_map, ripcas_required_data):
    """
    Simple version of the CASiMiR model for vegetation succession. Before the
    model is run, we check that all the unique values from vegetation_map are
    present in the shear_resistance_dict. Otherwise the process will fail
    wherever the vegetation map value is not present in the dictionary
    on lookup.

    Arguments:
        vegetation_map (str or ESRIAsc): location on disk or ESRIAsc
            representation of the vegetation map
        zone_map (str or ESRIAsc): location on disk or ESRIAsc representation
            of the zone map.
        shear_map (str or ESRIAsc): location on disk or ESRIAsc representation
            of the shear stress map
        ripcas_required_data (str): Excel spreadsheet of data needed for
            ripcas run. Encompasses the landscape model for the watershed.
            Can have one or two 'Code' columns and must have exactly one
            'shear_resis' column and exactly one 'n_val' column

    Returns:
        (ESRIAsc) vegetation map updated with new values corresponding to
            succession rules
    """
    if isinstance(vegetation_map, six.string_types):
        vegetation_map = ESRIAsc(vegetation_map)
    elif not isinstance(vegetation_map, ESRIAsc):
        raise TypeError('vegetation_map must be type str or ESRIAsc')

    if isinstance(shear_map, six.string_types):
        shear_map = ESRIAsc(shear_map)
    elif not isinstance(shear_map, ESRIAsc):
        raise TypeError('shear_map must be type str or ESRIAsc')

    if isinstance(zone_map, six.string_types):
        zone_map = ESRIAsc(zone_map)
    elif not isinstance(zone_map, ESRIAsc):
        raise TypeError('zone_map must be type str of ESRIAsc, not ' +
                        str(type(zone_map)))

    if isinstance(ripcas_required_data, six.string_types):
        cas_df = read_excel(ripcas_required_data)

        # sanity check to make sure our lookup is correct
        assert 'Code.1' in cas_df  # TODO generalize later
        assert 'shear_resis' in cas_df

        shear_resistance_dict = dict(
            zip(cas_df['Code.1'], cas_df['shear_resis'])
        )

    else:
        raise TypeError('shear_resistance_dict must be type str')

    # init the vegetation map that will be returned
    ret_veg_map = copy.deepcopy(vegetation_map)

    for idx in range(len(shear_map.data)):
        # determine whether or not the vegetation should be reset to age zero
        veg_val = int(vegetation_map.data[idx])
        if veg_val != 0:

            is_not_nodata = (
                shear_map.data[idx] != shear_map.NODATA_value and
                vegetation_map.data[idx] != vegetation_map.NODATA_value
            )

            veg_needs_reset = (
                is_not_nodata and
                shear_map.data[idx] > shear_resistance_dict[veg_val]
            )

            if veg_needs_reset:
                # reset vegetation to age zero with veg type appropriate to zone
                ret_veg_map.data[idx] = zone_map.data[idx]

            # whether or not the vegetation was destroyed, age by one
            if is_not_nodata:

                ret_veg_map.data[idx] += 1

    return ret_veg_map


def shear_mesh_to_asc(shear_nc_path, header_dict):
    """
    Extract flow element values and locations from the dflow output netcdf
    and project these onto the grid defined by the corner of the
    grid defined by the lower-left corner of the bounding box and the
    cell size. The results are saved in ESRI .asc format to asc_out_path.

    Arguments:
        shear_nc_path (str): location of the dflow netcdf on disk
        header_dict (dict): dictionary with six fields as required to build
            ESRIAsc; NODATA_value, cellsize, ncols, nrows, yllcorner, xllcorner

    Returns:
        (ESRIAsc) representation of gridded representation of the mesh shear
            stress data output from dflow
    """
    dflow_ds = Dataset(shear_nc_path, 'r')

    # the mesh locations are the x and y centers of the Flow (Finite) Elements
    mesh_x = dflow_ds.variables['FlowElem_xcc'][:]
    mesh_y = dflow_ds.variables['FlowElem_ycc'][:]
    mesh_shear = dflow_ds.variables['taus'][-1]  # take the last timestep

    cellsize = header_dict['cellsize']
    x = array([header_dict['xllcorner'] + (i*cellsize)
               for i in range(header_dict['ncols'])])

    y = array([header_dict['yllcorner'] + (i*cellsize)
               for i in range(header_dict['nrows'])])

    grid_x, grid_y = meshgrid(x, y)

    # use linear interp so we don't have to install natgrid
    asc_mat = griddata((mesh_x, mesh_y), mesh_shear, (grid_x, grid_y))

    # not sure why, but this makes it align with the original vegetation map
    asc_mat = flipud(asc_mat)

    data = reshape(asc_mat, (header_dict['ncols'] * header_dict['nrows']))

    data[isnan(data)] = int(header_dict['NODATA_value'])

    return ESRIAsc(data=Series(data), **header_dict)


def veg2n(veg_map, ripcas_required_data):
    """
    Creat an ESRIAsc representation of an ESRI .asc file that contains roughness
    values substituted for vegetation codes. The translation is found in the
    Excel file found at lookup_path.

    Arguments:
        veg_map (ESRIAsc): path to ESRI .asc file with vegetation codes
        ripcas_required_data (str): path to Excel file with vegetation codes mapped
            to Manning's roughness n-values. The Excel file must have four
            columns with headers
                Code	shear_resis	Code	n_val
            on the first sheet.

    Raises:
        (ValueError) if there is a vegetation code in the .asc that is not
            found in the lookup table

    Returns:
        (ESRIAsc) ESRI .asc map of Manning's n-values in place of veg codes
    """
    assert isinstance(veg_map, ESRIAsc), \
        "veg_map must be an instance of ESRIAsc"

    cas_df = read_excel(ripcas_required_data)
    veg2n_dict = dict(
        zip(cas_df['Code.1'], cas_df['n_val'])
    )

    ret = copy.deepcopy(veg_map)

    ret.data = veg_map.data.replace(veg2n_dict)

    return ret


class ESRIAsc:

    def __init__(self, file_path=None, ncols=None, nrows=None,
                 xllcorner=None, yllcorner=None, cellsize=1,
                 NODATA_value=-9999, data=None):

        self.file_path = file_path
        self.ncols = ncols
        self.nrows = nrows
        self.xllcorner = xllcorner
        self.yllcorner = yllcorner
        self.cellsize = cellsize
        self.NODATA_value = NODATA_value

        if data is not None and not isinstance(data, Series):
            raise RuntimeError('data value must be type pandas.Series')

        self.data = data

        # if a file is provided, the file metadata will overwrite any
        # user-provided kwargs
        if file_path:

            def getnextval(f):
                return f.readline().strip().split()[1]

            f = open(file_path, 'r')

            self.ncols = int(getnextval(f))
            self.nrows = int(getnextval(f))
            self.xllcorner = float(getnextval(f))
            self.yllcorner = float(getnextval(f))
            self.cellsize = int(getnextval(f))
            self.NODATA_value = float(getnextval(f))

            # should not be necessary for well-formed ESRI files, but
            # seems to be for CASiMiR
            data_str = ' '.join([l.strip() for l in f.readlines()])

            self.data = Series(fromstring(data_str, dtype=float, sep=' '))

            colrow_prod = self.nrows*self.ncols

            assert len(self.data) == colrow_prod, \
                "length of .asc data does not equal product of ncols * nrows" \
                "\nncols: {}, nrows: {}, ncols*nrows: {} len(data): {}".format(
                    self.ncols, self.nrows, colrow_prod, len(self.data))

    def header_dict(self):
        return dict(ncols=self.ncols,
                    nrows=self.nrows,
                    xllcorner=self.xllcorner,
                    yllcorner=self.yllcorner,
                    cellsize=self.cellsize,
                    NODATA_value=self.NODATA_value)

    def as_matrix(self, replace_nodata_val=None):
        """
        Convenience method to give 2D numpy.ndarray representation. If
        replace_nodata_val is given, replace all NODATA_value entries with
        it.

        Arguments:
            replace_nodata_val (float): value with which to replace
                NODATA_value entries

        Returns:
            (numpy.ndarray) matrix representation of the data in the .asc
        """
        ret = copy.copy(reshape(self.data, (self.nrows, self.ncols)))
        if replace_nodata_val is not None:
            ret[ret == self.NODATA_value] = replace_nodata_val

        return ret

    def write(self, write_path):
        # replace nan with NODATA_value
        self.data = self.data.fillna(self.NODATA_value)

        with open(write_path, 'w+') as f:
            f.write("ncols {}\n".format(self.ncols))
            f.write("nrows {}\n".format(self.nrows))
            f.write("xllcorner {}\n".format(self.xllcorner))
            f.write("yllcorner {}\n".format(self.yllcorner))
            f.write("cellsize {}\n".format(self.cellsize))
            f.write("NODATA_value {}\n".format(self.NODATA_value))

            # prob not most efficient, but CASiMiR requires
            # ESRI Ascii w/ newlines
            f.write(
                '\n'.join(
                    [
                        ' '.join([str(v) for v in row])
                        for row in self.as_matrix()
                    ]
                )
            )

    def __eq__(self, other):

        if isinstance(other, ESRIAsc):
            ret = self.ncols == other.ncols
            ret = self.nrows == other.nrows and ret
            ret = self.xllcorner == other.xllcorner and ret
            ret = self.yllcorner == other.yllcorner and ret
            ret = self.cellsize == other.cellsize and ret
            ret = self.NODATA_value == other.NODATA_value and ret
            ret = all(self.data == other.data) and ret

            return ret

        return NotImplemented


class Pol:
    """
    Wrapper creating and/or reading the .pol files of n-values for DFLOW
    """
    df = None

    x = None
    y = None
    n = None

    @classmethod
    def from_dflow_file(cls, pol_path=None):

        c = cls()
        if pol_path is not None:
            c.df = read_table(
                pol_path, skiprows=2, skipinitialspace=True, sep='   ',
                header=None, names=['x', 'y', 'z'], engine='python'
            )

            c.x = array(c.df.x)
            c.y = array(c.df.y)
            c.z = array(c.df.z)

        return c

    @classmethod
    def from_ascii(cls, asc):
        """
        Generate a polygon file from ESRIAsc
        """
        assert isinstance(asc, ESRIAsc), "asc input must be ESRIAsc"

        c = cls()
        grid = meshgrid(
            array(
                [asc.xllcorner + asc.cellsize*i for i in range(asc.ncols)]
            ),
            array(
                [asc.yllcorner + asc.cellsize*i for i in range(asc.nrows)]
            )
        )

        c.x = concatenate(grid[0])
        c.y = concatenate(grid[1])
        c.z = array(asc.data)

        return c

    @classmethod
    def from_river_geometry_file(cls, geom_file):

        c = cls()

        df = read_table(geom_file, sep='\s+', engine='python')
        c.df = df

        c.x = array(df.X)
        c.y = array(df.Y)
        c.z = array(df.Z)

        return c

    def write(self, out_path):
        """
        write the Pol to file

        Returns: None
        """
        if self.x is not None and self.y is not None and self.z is not None:

            assert len(self.x) == len(self.y) and len(self.y) == len(self.z), \
                "Lengths of columns is not equal!"

            with open(out_path, 'w') as f:
                f.write('L1\n')

                # weird, but apparently what DFLOW requires
                f.write(' '*4 + str(len(self.x)) + ' '*5 + '3\n')
                lines = (
                    ' '*3 + (' '*3).join(["{:1.7e}".format(val) for val in el])
                    for el in zip(self.x, self.y, self.z)
                )

                for l in lines:
                    f.write(l + '\n')

        else:
            raise Exception("Trying to write an incomplete polygon file")

    def __eq__(self, other):

        if isinstance(other, Pol):
            ret = all(self.x == other.x)
            ret = all(ret and self.y == other.y)
            ret = all(ret and self.z == other.z)

            return ret

        return NotImplemented


if __name__ == '__main__':

    help_msg = '''
dflow_ripcas v1.0  -  Matthew Turner <maturner01@gmail.com>

Usage:

    python jemez/dflow_ripcas.py <path-to-shear-nc> <path-to-save-output-vegmap>

Example:

    $ python jemez/dflow_ripcas.py ~/local_data/dflow_outputs/jemez_r02_map.nc ~/local_data/ripcas_out/veg-out-1.asc
    '''

    import sys
    try:
        shear_nc = sys.argv[1]
        vegout_path = sys.argv[2]

        if shear_nc == '-h':
            print help_msg
            sys.exit(0)
    except:
        print help_msg
        sys.exit(1)

    vegetation_map = ESRIAsc('data/vegclass_2z.asc')
    zone_map = ESRIAsc('data/zonemap_2z.asc')
    ripcas_data = 'data/ripcas-data-requirements.xlsx'

    vegout = ripcas_with_dflow_io(
        vegetation_map, zone_map, shear_nc, ripcas_data
    )

    vegout.write(vegout_path)
