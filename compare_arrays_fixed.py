import sys
import numpy as np
import numpy.testing as npt

fname1, fname2 = sys.argv[1], sys.argv[2]
arrays = []

for fname in fname1, fname2:
    f = open(fname, "rb")
    f.seek(8)   # Skip the two first integers.
    arrays.append(np.fromfile(f, dtype = float, count = -1, sep = ""))
    f.close()

npt.assert_approx_equal(actual = arrays[0], desired = arrays[1], significant = 12)

print("Arrays match to within tolerance.")
