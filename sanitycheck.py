"""
I'm not proud I made this.

Reads a matrix file, writes it to ASCII.
"""

import sys
import numpy as np
import numpy.testing as npt

fname_whole_matrix, fname2 = sys.argv[1], sys.argv[2]

# Read large matrix.
f = open(fname_whole_matrix, "rb")
m, n = np.fromfile(f, dtype = np.int32, count = 2, sep = "")
whole_matrix = np.fromfile(f, dtype = np.double, count = m * n, sep = "")
whole_matrix.shape = (m, n)
f.close()

f = open(fname2, "w")
for i in xrange(m):
    for j in xrange(n):
        f.write("%02.03f " %whole_matrix[i, j])
    f.write("\n")

f.close()
