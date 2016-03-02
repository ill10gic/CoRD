"""
Utilities for interacting with dflow and casimir models via the virtual
watershed.
"""
import copy
import json

from numpy import fromstring, reshape
from pandas import Series


def casimir(vegetation_map, shear_map, shear_resistance_dict):
    """
    Simple version of the CASiMiR model for vegetation succession. Before the
    model is run, we check that all the unique values from vegetation_map are
    present in the shear_resistance_dict. Otherwise the process will fail
    wherever the vegetation map value is not present in the dictionary
    on lookup.

    Arguments:
        vegetation_map (str or ESRIAsc): location on disk or ESRIAsc
            representation of the vegetation map
        shear_map (str or ESRIAsc): location on disk or ESRIAsc representation
            of the shear stress map
        shear_resistance_dict (str or dict): location on disk or dictionary
            representation of the resistance dictionary that maps
            vegetation type to shear resistance.

    Returns:
        (ESRIAsc) vegetation map updated with new values corresponding to
            succession rules
    """
    if type(vegetation_map) is str:
        vegetation_map = ESRIAsc(vegetation_map)
    elif not isinstance(vegetation_map, ESRIAsc):
        raise TypeError('vegetation_map must be type str or ESRIAsc')

    if type(shear_map) is str:
        shear_map = ESRIAsc(shear_map)
    elif not isinstance(shear_map, ESRIAsc):
        raise TypeError('shear_map must be type str or ESRIAsc')

    if type(shear_resistance_dict) is str:
        try:
            shear_resistance_dict = json.load(open(shear_resistance_dict))
        except ValueError:
            raise ValueError(
                'The shear_resistance_dict file is not valid JSON!'
            )
    elif not isinstance(shear_resistance_dict, dict):
        raise TypeError('shear_resistance_dict must be type str or dict')

    # init the vegetation map that will be returned
    ret_veg_map = copy.deepcopy(vegetation_map)

    for idx in range(len(shear_map.data)):
        # determine whether or not the vegetation should be reset to age zero
        shear_val_int = str(int(vegetation_map.data[idx]))

        is_not_nodata = shear_map.data[idx] != shear_map.NODATA_value

        veg_needs_reset = (
            is_not_nodata and
            shear_map.data[idx] > shear_resistance_dict[shear_val_int]
        )

        if veg_needs_reset:
            # reset vegetation to age zero while retaining veg type
            ret_veg_map.data[idx] -= vegetation_map.data[idx] % 100

        # whether or not the vegetation was destroyed, age by one
        if is_not_nodata:

            ret_veg_map.data[idx] += 1

    return ret_veg_map


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
