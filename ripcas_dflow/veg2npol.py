"""
Script to convert an ESRI .asc vegetation map to an equivalent polygon (.pol)
"""
import sys

from dflow_casimir import veg2n, Pol, ESRIAsc

if __name__ == '__main__':

    help_message = \
'''
Run this by typing

python jemez/veg2npol.py <veg-map-name> <output-name>

Example:

    python jemez/veg2npol.py data/vegclass_2z.asc data/initial_n.pol
'''

    if len(sys.argv) != 3:

        print help_message

        sys.exit(1)

    if sys.argv[1] == '-h':

        print help_message

        sys.exit(0)

    asc = veg2n(ESRIAsc(sys.argv[1]), 'data/casimir-data-requirements.xlsx')
    pol = Pol.from_ascii(asc)

    pol.write(sys.argv[2])
