"""
To run:
python checkmatrixdist.py <WHOLE_MATRIX> <SMALL_MATRICES> <PROCS_PER_DIM>

where <SMALL_MATRICES> is the common part of the filenames of the block matrices.
example:
python checkmatrixdist.py small_matrix_a.bin small_matrices_a 2
for a (2 x 2) setup.
"""

import sys
import numpy as np
import numpy.testing as npt

# Assuming a matrix is composed into p^2 blocks, each written to file. Reads the whole matrix and all the little matrices, verify that they match.

try:
    fname_whole_matrix, fname2_common = sys.argv[1], sys.argv[2]
    p = int(sys.argv[3])

except IndexError:
    print "to run %s <WHOLE_MATRIX> <SMALL_MATRICES_COMMON> <PROCS_PER_DIM>"
    exit(0)

# Read large matrix.
f = open(fname_whole_matrix, "rb")
m, n = np.fromfile(f, dtype = np.int32, count = 2, sep = "")
whole_matrix = np.fromfile(f, dtype = np.double, count = m * n, sep = "")
whole_matrix.shape = (m, n)
f.close()

small_matrices = []
small_dims = []

# Read small matrices.
for i in xrange(p):
    small_matrices.append([])
    small_dims.append([])
    for j in xrange(p):
        fname = "%s_%02d_%02d.bin" %(fname2_common, i, j)

        f = open(fname, "rb")

        m_local, n_local = np.fromfile(f, dtype = np.int32, count = 2, sep = "")

        small_dims[i].append((m_local, n_local))
        small_matrices[i].append(np.fromfile(f, dtype = np.double, count = m_local * n_local, sep = ""))
        small_matrices[i][j].shape = (m_local, n_local)

        f.close()

# Finally, compare. The matrices should match exactly.
try:
    for i in xrange(p):

        dispy = i * int(m / p) + min(i, m % p)
        for j in xrange(p):

            dispx = j * int(n / p) + min(j, n % p)
            m_local, n_local = small_dims[i][j]

            if ~np.all((whole_matrix[dispy : dispy + m_local, dispx : dispx + n_local] == small_matrices[i][j])):
                print("Mismatch in matrix with coords (%d, %d)" %(i, j))
                raise ValueError

    print "Matrices were distributed properly!"

except ValueError:
    print "Matrices were not distributed properly after all. SAD"
