#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void read_matrix_binaryformat (char* filename, double*** matrix, int* num_rows, int* num_cols) {
    int i;
    FILE* fp = fopen (filename,"rb");
    fread (num_rows, sizeof(int), 1, fp);
    fread (num_cols, sizeof(int), 1, fp);
    /* storage allocation of the matrix */
    *matrix = (double**)malloc((*num_rows)*sizeof(double*));
    (*matrix)[0] = (double*)malloc((*num_rows)*(*num_cols)*sizeof(double));
    for (i=1; i<(*num_rows); i++)
    (*matrix)[i] = (*matrix)[i-1]+(*num_cols);
    /* read in the entire matrix */
    fread ((*matrix)[0], sizeof(double), (*num_rows)*(*num_cols), fp);
    fclose (fp);
}

void write_matrix_binaryformat (char* filename, double** matrix, int num_rows, int num_cols) {
  FILE *fp = fopen (filename,"wb");
  fwrite (&num_rows, sizeof(int), 1, fp);
  fwrite (&num_cols, sizeof(int), 1, fp);
  fwrite (matrix[0], sizeof(double), num_rows*num_cols, fp);
  fclose (fp);
}

int main(int argc, char *argv[]) {
    if (argc == 2) {
      matrix_a = argv[1];
      matrix_b = argv[2];
    }
    else{
      printf("2 arguments expected");
      exit (EXIT_FAILURE);
    }

accept two file names at run-time,
• let process 0 read the A and B matrices from the two data files,
• let process 0 distribute the pieces of A and B, after a 2D partitioning,
to all the other processes,
• let each process calculate its piece of C = A ∗ B in parallel,
• let process 0 gather, from all the other processes, the different pieces
of C,
• let process 0 write out the entire C matrix to an output data file.
