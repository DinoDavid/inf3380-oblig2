#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>

//matrix_io
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

void allocate_matrix(double** matrix, int num_rows, int num_cols) {
  matrix_c = malloc(rows_a*sizeof(double*));
  matrix_c[0] = malloc(rows_a * cols_b * sizeof(double));
  for (int i=1; i < rows_a; i++){
    matrix_c[i] = (matrix_c)[i-1] + (cols_a);
  }
}

void matrix_dist(double** matrix_a, double** matrix_b, ) {
  /* Buffers. */
  double *senddata_columnwise, *senddata_rowwise;

  /* Scatterv variables, step 1. */
  int *displs_y, *sendcounts_y, *everyones_m;

  /* Scatterv variables, step 2. */
  int *displs_x, *sendcounts_x, *everyones_n;
  MPI_Datatype columntype, columntype_scatter, columntype_recv, columntype_recv_scatter;


  everyones_m = (int *) calloc(procs_per_dim, sizeof(int));
  sendcounts_y = (int *)calloc(procs_per_dim, sizeof(int));
  displs_y = (int *)calloc(procs_per_dim + 1, sizeof(int));


  displs_y[0] = 0;
  for (int i = 0; i < procs_per_dim; ++i)
  {
      sendcounts_y[i] = n * everyones_m[i];
      displs_y[i + 1] = displs_y[i] + sendcounts_y[i];
  }
}

void gather_matrix() {

}
//matrix_io ends

//Cannons algorithm and matrix matrix multiply

/* This matrix performs a serial matrix-matrix multiplication c = a * b. */
void matrix_mult(double** matrix_a, int rows_a, int cols_a, double** matrix_b, int cols_b, double **matrix_c) {
  for (int i = 0; i < rows_a; i++) {
    for (int j = 0; j < cols_b; j++) {
      matrix_c[i][j] = 0;
      for (int k = 0; k < cols_a; k++) {
        matrix_c[i][j] += matrix_a[j][k] * matrix_b[k][j];
      }
    }
  }
}

void cannonMultiply(int my_rows_a, int my_cols_a, int my_cols_b, double *matrix_a, double *matrix_b, double *matrix_c, MPI_Comm comm) {
    int num_procs, dims[2], periods[2];
    int myrank, my2drank, mycoords[2];
    int uprank, downrank, leftrank, rightrank;
    int shiftsource, shiftdest;
    MPI_Status status;
    MPI_Comm comm_2d;

	/* Get the communicator related information */
    MPI_Comm_size(comm, &num_procs);
    MPI_Comm_rank(comm, &myrank);

	/* Set up the Cartesian topology */
    dims[0] = dims[1] = sqrt(num_procs);

	/* Set the periods for wraparound connections, 1 == true*/
    periods[0] = periods[1] = 1;

	/* Create the Cartesian topology, with rank reordering */
    MPI_Cart_create(comm, 2, dims, periods, 1, &comm_2d);

	/* Get the rank and coordinates with respect to the new topology */
    MPI_Comm_rank(comm_2d, &my2drank);
    MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);

	/* Compute ranks of the up and left shifts */
    MPI_Cart_shift(comm_2d, 1, -1, &rightrank, &leftrank);
    MPI_Cart_shift(comm_2d, 0, -1, &downrank, &uprank);

	/* Perform the initial matrix alignment. First for A and then for B */
    MPI_Cart_shift(comm_2d, 1, -mycoords[0], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(matrix_a, my_rows_a*my_cols_a, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_2d, &status);

    MPI_Cart_shift(comm_2d, 0, -mycoords[1], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(matrix_b, my_cols_a*my_cols_b, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_2d, &status);

	/* Get into the main computation loop */
    for (int i = 0; i < dims[0]; i++) {
        matrix_mult(matrix_a, rows_a, cols_a, matrix_b, cols_b, matrix_c);

	    /* Shift matrix a left by one */
        MPI_Sendrecv_replace(matrix_a, my_rows_a*my_cols_a, MPI_DOUBLE, leftrank, 1, rightrank, 1, comm_2d, &status);

	    /* Shift matrix b up by one */
        MPI_Sendrecv_replace(matrix_b, my_cols_a*my_cols_b, MPI_DOUBLE, uprank, 1, downrank, 1, comm_2d, &status);
    }
  /* Free up communicator */
    MPI_Comm_free(&comm_2d);
}

int main(int argc, char *argv[]) {
  double **matrix_a, **matrix_b, **matrix_c;
  int rows_a, cols_a, rows_b, cols_b;

  //mpi variables
  int my_rank, num_procs, num_procs_sqrt;

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);

  if (argc != 4) {
    if (my_rank == 0) {
      printf("3 arguments expected.\n");
      MPI_Abort(MPI_COMM_WORLD, 0);
      exit (EXIT_FAILURE);
    }
  }

  num_procs_sqrt = sqrt(num_procs);
  if (num_procs_sqrt * num_procs_sqrt != num_procs) {
    if (my_rank == 0) {
      printf("Number of processes must be a square number.\n");
      MPI_Abort(MPI_COMM_WORLD, 0);
    }
  }
  if (my_rank == 0) {
    read_matrix_binaryformat((argv[1]), &matrix_a, &rows_a, &cols_a);
    read_matrix_binaryformat((argv[2]), &matrix_b, &rows_b, &cols_b);
    //matrix_dist();
  }
  if (my_rank == 0) {
    write_matrix_binaryformat(argv[3], matrix_c, rows_a, cols_b);
  }

  free(matrix_c);

  MPI_Finalize ();
  return 0;
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

int main(){
  int my_rank, num_procs, my2drank, mycoords[2];
  MPI_Comm comm_2d, comm_col, comm_row;
  int m, int n, my_rows_a, my_cols_b;
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

  my_rows_a = m/num_procs_sqrt + (mycoords[0] < m%num_procs_sqrt);
  my_cols_b = m/num_procs_sqrt + (mycoords[0] < m%num_procs_sqrt);

  if(mycoord[1] == 0) {
    if(mycoords[0] == 0) {
      everyones_m = (int*) calloc(num_procs_sqrt, sizeof(int));
    }
    MPI_Gather(&my_rows_a, 1, MPI_INT, everyones_m, 1, MPI_INT, 0, comm_col);
    if(mycoords[0] == 0) {
      sendcounts_y = (int *) calloc(num_procs_sqrt, sizeof(int));
      displs_y = (int*)
      for (int i = 0; i < num_procs_sqrt; i++){
        sendcounts_y[i] = n*everyones_m[i];
        displs_y[i+1] = displs_y[i] * sendcounts_y[i];
      }
    }
    sendcounts_rowwise = (double *) calloc(my_rows_a*n, sizeof(double));
    if(mycoords[0] == 0) {
      senddata_columnwise = a[0];
    }
    MPI_Scatter(senddata_colummnwise, );
  }

  MPI_type_vector(my_rows_a, 1, n, MPI_DOUBLE, &columntype);
  MPI_Type_commit(&Columntype);
  MPI_Type_create_resized(columntype, 0, sizeof(double), &columntype_scatter);
  MPI_Type_commit(&coulmntype_Scatter);

  if (mycoords[1] == 0) {
    everyones_n = (int*) calloc(num_procs_sqrt, sizeof(int));
    sendcount_x = (int*) calloc(num_procs_sqrt, sizeof(int));
    displs_x = (int*) calloc(num_procs_sqrt + 1, sizeof(int));
  }

  MPI_Gather(&my_cols_b, 1, MPI_INT, everyones_n, 1, MPI_INT, 0, comm_row);
  if(mycoords[1] == 0) {
    displs_x[1] == 0)
    for (int = i 0; i < num_procs_sqrt; i++){
      sendcounts_x[i] = everyones_n[i];
      displs_x[i+1] = displs_x[i] + sendcounts_x[i];
    }
  }
  allocate_matrix(&my_a, my_rows_a, my_rows_a);
  MPI_Scatterv(senddata_rowwise, sendcounts_x, displs_x, my_a, my_rows_a * my_cols_b, MPI_DOUBLE, 0, comm_row);

  if (mycoords[1] = 0){
    free(displs_x);
    free(sendcounts_x);
  }
  MPI_Finalize();
  return 0;
}


side 254 i boka
*/
