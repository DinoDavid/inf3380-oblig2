#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

void allocate_matrix(double*** matrix, int m, int n) {
  (*matrix) = malloc(m * sizeof(double*));
  (*matrix)[0] = malloc(m * n * sizeof(double));
  for (int i = 1; i < m; i++){
    (*matrix)[i] = (*matrix)[i-1] + n;
  }
}
/*
void deallocate_matrix(double* matrix) {
  free(matrix);
}
*/

/*
void gather_matrix() {

}
*/
//Cannons algorithm and matrix matrix multiply

/* This matrix performs a serial matrix-matrix multiplication c = a * b. */
void matrix_mult(double** matrix_a, int rows_a, int cols_a, double** matrix_b, int cols_b, double** matrix_c) {
  for (int i = 0; i < rows_a; i++) {
    for (int j = 0; j < cols_b; j++) {
      matrix_c[i][j] = 0;
      for (int k = 0; k < cols_a; k++) {
        matrix_c[i][j] += matrix_a[i][k] * matrix_b[k][j];
      }
    }
  }
}

void cannon_mult(int rows_apart, int cols_apart, int cols_bpart, double **matrix_a, double **matrix_b, double **matrix_c, MPI_Comm comm) {
    int num_procs, dims[2], periods[2];
    int myrank, my2drank, mycoords[2];
    int uprank, downrank, leftrank, rightrank;
    int shiftsource, shiftdest;
    MPI_Status status;
    MPI_Comm comm_2d;

	//Get the communicator related information
    MPI_Comm_size(comm, &num_procs);
    MPI_Comm_rank(comm, &myrank);

	//Set up the Cartesian topology
    dims[0] = dims[1] = sqrt(num_procs);

	//Set the periods for wraparound connections, 1 == true
    periods[0] = periods[1] = 1;

	//Create the Cartesian topology, with rank reordering
    MPI_Cart_create(comm, 2, dims, periods, 1, &comm_2d);

	//Get the rank and coordinates with respect to the new topology
    MPI_Comm_rank(comm_2d, &my2drank);
    MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);

	//Compute ranks of the up and left shifts
    MPI_Cart_shift(comm_2d, 1, -1, &rightrank, &leftrank);
    MPI_Cart_shift(comm_2d, 0, -1, &downrank, &uprank);

	//Perform the initial matrix alignment. First for A and then for B
    MPI_Cart_shift(comm_2d, 1, -mycoords[0], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(matrix_a, rows_apart * cols_apart, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_2d, &status);

    MPI_Cart_shift(comm_2d, 0, -mycoords[1], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(matrix_b, cols_apart * cols_bpart, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_2d, &status);

	// Get into the main computation loop
    for (int i = 0; i < dims[0]; i++) {
        matrix_mult(matrix_a, rows_apart, cols_apart, matrix_b, cols_bpart, matrix_c);

	    // Shift matrix a left by one
        MPI_Sendrecv_replace(matrix_a, rows_apart * cols_apart, MPI_DOUBLE, leftrank, 1, rightrank, 1, comm_2d, &status);

	    // Shift matrix b up by one
        MPI_Sendrecv_replace(matrix_b, cols_apart * cols_bpart, MPI_DOUBLE, uprank, 1, downrank, 1, comm_2d, &status);
    }
  // Free up communicator
    MPI_Comm_free(&comm_2d);
}



int main(int argc, char *argv[]) {
  //variables
  double **matrix_a, **matrix_b, **matrix_c;
  int rows_a, cols_a, rows_b, cols_b, rows_c, cols_c;
  double **A_part, **B_part, **C_part;
  int *row_cnt_a = 0, *row_cnt_b = 0, *col_cnt_a = 0, *col_cnt_b = 0;
  int *row_displ_a = 0, *row_displ_b = 0, *col_displ_a = 0, *col_displ_b = 0;

  //mpi variables
  int rows_apart, cols_apart, rows_bpart, cols_bpart, rows_cpart, cols_cpart;
  int my_rank, num_procs, num_procs_sqrt;
  int dims[2], periods[2];
  MPI_Comm comm_2d, comm_rows, comm_cols;
  int rank2d, rankrows, rankcols;
  int mycoords[2];

  //mpi begin
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);

  //arguments check
  if (argc != 4) {
    printf("3 arguments expected.\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
    exit (EXIT_FAILURE);
  }

  //length of process mesh
  num_procs_sqrt = sqrt(num_procs);
  if (num_procs_sqrt * num_procs_sqrt != num_procs) {
    printf("Number of processes must be a square number.\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  //read files and allocate a,b,c
  if (my_rank == 0) {
    read_matrix_binaryformat((argv[1]), &matrix_a, &rows_a, &cols_a);
    read_matrix_binaryformat((argv[2]), &matrix_b, &rows_b, &cols_b);
    allocate_matrix(&matrix_c, rows_a, cols_b);
  }

  //calculate displs
  if (my_rank == 0){
    row_cnt_a = malloc(num_procs_sqrt * sizeof(int));
    row_cnt_b = malloc(num_procs_sqrt * sizeof(int));
    col_cnt_a = malloc(num_procs_sqrt * sizeof(int));
    col_cnt_b = malloc(num_procs_sqrt * sizeof(int));

    row_displ_a = calloc((num_procs_sqrt + 1), sizeof(int));
    row_displ_b = calloc((num_procs_sqrt + 1), sizeof(int));
    col_displ_a = calloc((num_procs_sqrt + 1), sizeof(int));
    col_displ_b = calloc((num_procs_sqrt + 1), sizeof(int));

    for (int i = 0; i < num_procs_sqrt; i++){
      row_cnt_a[i] = rows_a/num_procs_sqrt + (i < rows_a % num_procs_sqrt);
      row_cnt_b[i] = rows_b/num_procs_sqrt + (i < rows_b % num_procs_sqrt);
      col_cnt_a[i] = cols_a/num_procs_sqrt + (i < cols_a % num_procs_sqrt);
      col_cnt_b[i] = cols_b/num_procs_sqrt + (i < cols_b % num_procs_sqrt);

      row_displ_a[i+1] = row_displ_a[i] + row_cnt_a[i];
      row_displ_b[i+1] = row_displ_b[i] + row_cnt_b[i];
      col_displ_a[i+1] = col_displ_a[i] + col_cnt_a[i];
      col_displ_a[i+1] = col_displ_a[i] + col_cnt_a[i];
    }
  }

  //m and n share
  rows_c = rows_a;
  cols_c = cols_b;
  MPI_Bcast (&rows_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&cols_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&rows_b, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&cols_b, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&rows_c, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&cols_c, 1, MPI_INT, 0, MPI_COMM_WORLD);

  //Create Cartesian topology
  dims[0] = dims[1] = num_procs_sqrt;
  periods[0] = periods[1] = 1;

  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_2d);
  MPI_Comm_rank(comm_2d, &rank2d);
  MPI_Cart_coords(comm_2d, rank2d, 2, mycoords);

  MPI_Cart_sub(comm_2d, (int[]){1,0}, &comm_rows);
  MPI_Cart_sub(comm_2d, (int[]){0,1}, &comm_cols);
  MPI_Comm_rank(comm_rows, &rankrows);
  MPI_Comm_rank(comm_cols, &rankcols);

  //find partion sizes
  rows_apart = rows_a / num_procs_sqrt + (rankrows < rows_a % num_procs_sqrt);
  cols_apart = cols_a / num_procs_sqrt + (rankcols < cols_a % num_procs_sqrt);
  rows_bpart = rows_b / num_procs_sqrt + (rankrows < rows_b % num_procs_sqrt);
  cols_bpart = cols_b / num_procs_sqrt + (rankcols < cols_b % num_procs_sqrt);
  rows_cpart = rows_c / num_procs_sqrt + (rankrows < rows_c % num_procs_sqrt);
  cols_cpart = cols_c / num_procs_sqrt + (rankcols < cols_c % num_procs_sqrt);

  //allocate partitions
  allocate_matrix(&A_part, rows_apart+1, cols_apart+1);
  allocate_matrix(&B_part, rows_bpart+1, cols_bpart+1);
  allocate_matrix(&C_part, rows_cpart+1, cols_cpart+1);

  //distribute partitions
  if (my_rank == 0) {
    for (int i = 0; i < rows_apart; i++){
      memcpy(A_part[i], matrix_a[i], cols_apart * sizeof(double));
    }
    for (int i = 0; i < rows_bpart; i++){
      memcpy(B_part[i], matrix_b[i], cols_bpart * sizeof(double));
    }
    for (int i = 1; i < num_procs; i++) {
      int Am, An, Bm, Bn;
      MPI_Recv(&Am, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&An, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&Bm, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&Bn, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      int moff, noff;
      for (int j = 0; j < row_cnt_a[i / num_procs_sqrt]; j++) {
        moff = row_displ_a[i / num_procs_sqrt];
        noff = col_displ_a[i % num_procs_sqrt];
        MPI_Send(&(matrix_a[j + moff][noff]), An, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
      }
      for (int j = 0; j <  row_cnt_b[i / num_procs_sqrt]; j++){
        moff = row_displ_b[i / num_procs_sqrt];
        noff = col_displ_b[i % num_procs_sqrt];
        MPI_Send(&(matrix_b[j + moff][noff]), Bn, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
      }
    }
  }else{
    MPI_Send(&rows_apart, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&cols_apart, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&rows_bpart, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&cols_bpart, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

    for (int j = 0; j < rows_apart; j++){
      MPI_Recv(A_part[j], cols_apart, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    for (int j = 0; j < rows_bpart; j++){
      MPI_Recv(B_part[j], cols_bpart, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }

  //save result
  if (my_rank == 0){
    //cannon_mult(rows_a, cols_a, cols_b, matrix_a, matrix_b, matrix_c, MPI_COMM_WORLD);
    matrix_mult(matrix_a, rows_a, cols_a, matrix_b, cols_b, matrix_c);
    write_matrix_binaryformat(argv[3], matrix_c, rows_c, cols_c);
  }

  //cleanup
  MPI_Finalize();
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
