#include <stdio.h>
#include <stdlib.h>
#include <string.h>s

void allocate_image (image *u, int m, int n)
{

}

accept two file names at run-time,
• let process 0 read the A and B matrices from the two data files,
• let process 0 distribute the pieces of A and B, after a 2D partitioning,
to all the other processes,
• let each process calculate its piece of C = A ∗ B in parallel,
• let process 0 gather, from all the other processes, the different pieces
of C,
• let process 0 write out the entire C matrix to an output data file.
