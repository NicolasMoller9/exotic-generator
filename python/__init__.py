# $Id$

import sys
import math
from icecube import icetray

try:
    icetray.load('libexotic-generator', False)
except RuntimeError:
    sys.stderr.write("ERROR: Could not load libexotic-generator (%s)." % sys.exc_info()[1])
del sys


def calculate_weight(EXOInfoDict_WeightFilewideNormalized, exogen_DiskRadius_in_meter):
    """Returns weight so you only have to multiply the weight with the expected flux in cm^-2 s^-1 sr^-1
    and divide by the number of files in a speed region (i.e. if you have overlapping simulation, the overlapped region needs to be
    divided by the number of overlapping datasets). If you are not using the MonopoleRenormalizeWeights module, you
    need to adjust your weights by the integrated generated weight you will find in your logs.
    See Anna Pollmanns thesis, page 84, for an explanation of the formula"""
    return 4 * math.pi ** 2 * 10000 * exogen_DiskRadius_in_meter**2 * EXOInfoDict_WeightFilewideNormalized

