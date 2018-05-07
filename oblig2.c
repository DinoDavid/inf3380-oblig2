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

void matrix_mult(double** matrix_a, int num_rows_a, int num_cols_a, double** matrix_b, int num_rows_b, int num_cols_b, double **matrix_c) {
  for (int i = 0; i < num_rows_a; i++){
    for (int j = 0; j < num_cols_b; j++) {
      matrix_c[i][j] = 0;
      for (int k = 0; k < num_cols_a; k++) {
        matrix_c[i][j] += matrix_a[j][k] * matrix_b[k][j];
      }
    }
  }
}

int main(int argc, char *argv[]) {
    double **matrix_a, **matrix_b, **matrix_c;
    int *num_rows_a, *num_cols_a, *num_rows_b, *num_cols_b;

    if (argc ==! 2) {
      printf("2 arguments expected");
      exit (EXIT_FAILURE);
    }

    read_matrix_binaryformat((argv[1]), matrix_a, num_rows_a, num_cols_a);
    read_matrix_binaryformat((argv[2]), matrix_b, num_rows_b, num_cols_b);
    *matrix_c = malloc((*num_rows_a)*(*num_cols_b)*sizeof(double*));

    matrix_mult(matrix_a, num_rows_a, num_cols_a, matrix_b, num_rows_a, num_cols_b, matrix_c);
    write_matrix_binaryformat(argv[3], matrix_c, num_rows_a, num_cols_b);
}

/*
• let process 0 read the A and B matrices from the two data files,
• let process 0 distribute the pieces of A and B, after a 2D partitioning,
to all the other processes,
• let each process calculate its piece of C = A ∗ B in parallel,
• let process 0 gather, from all the other processes, the different pieces
of C,
• let process 0 write out the entire C matrix to an output data file.


MPI_Send(buff)

Ikke bruke mpi scatter til å fordele


kartesisk kommunikator


MPI_COMM_WORLD

MPI_CART...
______________________________________________________________

matrixfuncs.c

matrix_dist.c

inf main(){
  int myrnakp
  int my_rank, num_procs, my2drank, mycoords[2];
  MPI_Comm comm_2d, comm_col, comm_row;
  int m, int n, my_m, my_n;
  double **a_global, **a_local;
  double **senddata_collwise, **senddata_rowise;

  int *displs_y, *sendcounts_y, everyones_y;
  MPI_Datatype collumtype,
  MPI_Init(&argc, &argv);
  MPI_Comm_size();
  MPI_Comm_rank();

  num_procs_sqrt = sqrt(num_procs);

  dims[0] = dims[1] = num_procs;
  periods[0] = periods[1] = 1;
  MPI_Cart_Create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_2d);

  MPI_Comm_rank(comm_2d, &my2drank);
  MPI_comm_rank(comm_2d, &my_rank);

  MPI_Cart_sub(comm_2d, (int[]){0,1}, &comm_row);
  MPI_Cart_sub(comm_2d, (int[]){0,1}, &comm_col);

  if(my2drank == 0){
    read_matrix_binaryformat(argc[1], &a, &m, &n);
  }
  MPI_Bcast(&m, 1, MPI_INT, 0, comm_2d);
  MPI_Bcast(&n, 1, MPI_INT, 0, comm_2d);

  my_m = m/num_procs_sqrt + (mycoords[0] < m%num_procs_sqrt);
  my_n = m/num_procs_sqrt + (mycoords[0] < m%num_procs_sqrt);

  if(mycoord[1] == 0) {
    if(mycoords[0] == 0) {
      everyones_m = (int*) calloc(num_procs_sqrt, sizeof(int));
    }
    MPI_Gather(&my_m, 1, MPI_INT, everyones_m, 1, MPI_INT, 0, comm_col);
    if(mycoords[0] == 0) {
      sendcounts_y = (int *) calloc(num_procs_sqrt, sizeof(int));
      displs_y = (int*)
      for (int i = 0; i < num_procs_sqrt; i++){
        sendcounts_y[i] = n*everyones_m[i];
        displs_y[i+1] = displs_y[i] * sendcounts_y[i];
      }
    }
    sendcounts_rowwise = (double *) calloc(my_m*n, sizeof(double));
    if(mycoords[0] == 0) {
      senddata_columnwise = a[0];
    }
    MPI_Scatter(senddata_colummnwise, );
  }

  MPI_type_vector(my_m, 1, n, MPI_DOUBLE, &columntype);
  MPI_Type_commit(&Columntype);
  MPI_Type_create_resized(columntype, 0, sizeof(double), &columntype_scatter);
  MPI_Type_commit(&coulmntype_Scatter);

  if (mycoords[1] == 0) {
    everyones_n = (int*) calloc(num_procs_sqrt, sizeof(int));
    sendcount_x = (int*) calloc(num_procs_sqrt, sizeof(int));
    displs_x = (int*) calloc(num_procs_sqrt + 1, sizeof(int));
  }

  MPI_Gather(&my_n, 1, MPI_INT, everyones_n, 1, MPI_INT, 0, comm_row);
  if(mycoords[1] == 0) {
    displs_x[1] == 0)
    for (int = i 0; i < num_procs_sqrt; i++){
      sendcounts_x[i] = everyones_n[i];
      displs_x[i+1] = displs_x[i] + sendcounts_x[i];
    }
  }
  allocate_matrix(&my_a, my_m, my_m);
  MPI_Scatterv(senddata_rowwise, sendcounts_x, displs_x, my_a, my_m * my_n, MPI_DOUBLE, 0, comm_row);

  if (mycoords[1] = 0){
    free(displs_x);
    free(sendcounts_x);
  }
  MPI_Finalize();
  return 0;
}


side 254 i boka
*/
