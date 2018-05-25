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

void deallocate_matrix(double** matrix) {
  free(matrix[0]);
  free(matrix);
}

// This matrix performs a serial matrix-matrix multiplication c = a * b.
void matrix_mult(double** matrix_a, double** matrix_b, double** matrix_c, int rows_a, int cols_a, int cols_b) {
#pragma omp parallel for num_threads(4)
  for (int i = 0; i < rows_a; i++) {
    for (int j = 0; j < cols_b; j++) {
      //matrix_c[i][j] = 0;
      for (int k = 0; k < cols_a; k++) {
        matrix_c[i][j] += matrix_a[i][k] * matrix_b[k][j];
      }
    }
  }
}

int main(int argc, char *argv[]) {
  //variables
  double **matrix_a, **matrix_b, **matrix_c;
  int rows_a = 0, cols_a = 0, rows_b = 0, cols_b = 0, rows_c = 0, cols_c = 0;
  double **A_part, **B_part, **C_part;
  int *row_cnt_a = 0, *row_cnt_b = 0, *row_cnt_c = 0, *col_cnt_a = 0, *col_cnt_b = 0, *col_cnt_c = 0;
  int *row_displ_a = 0, *row_displ_b = 0, *row_displ_c = 0, *col_displ_a = 0, *col_displ_b = 0, *col_displ_c = 0;

  //mpi variables
  int rows_apart, cols_apart, rows_bpart, cols_bpart, rows_cpart, cols_cpart;
  int my_rank, num_procs, num_procs_sqrt;
  int dims[2], periods[2];
  MPI_Comm comm_2d, comm_rows, comm_cols;
  int rank2d, rankrows, rankcols;
  int mycoords[2];
  int uprank, downrank, leftrank, rightrank;
  int shiftsource, shiftdest;
  int A_part_m_max, A_part_n_max, B_part_m_max, B_part_n_max, C_part_m_max, C_part_n_max;

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
    read_matrix_binaryformat(argv[1], &matrix_a, &rows_a, &cols_a);
    read_matrix_binaryformat(argv[2], &matrix_b, &rows_b, &cols_b);
    allocate_matrix(&matrix_c, rows_a, cols_b);
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

  //calculate displs
  if (my_rank == 0){
    row_cnt_a = malloc(num_procs_sqrt * sizeof(int));
    row_cnt_b = malloc(num_procs_sqrt * sizeof(int));
    row_cnt_c = malloc(num_procs_sqrt * sizeof(int));
    col_cnt_a = malloc(num_procs_sqrt * sizeof(int));
    col_cnt_b = malloc(num_procs_sqrt * sizeof(int));
    col_cnt_c = malloc(num_procs_sqrt * sizeof(int));

    row_displ_a = calloc((num_procs_sqrt + 1), sizeof(int));
    row_displ_b = calloc((num_procs_sqrt + 1), sizeof(int));
    row_displ_c = calloc((num_procs_sqrt + 1), sizeof(int));
    col_displ_a = calloc((num_procs_sqrt + 1), sizeof(int));
    col_displ_b = calloc((num_procs_sqrt + 1), sizeof(int));
    col_displ_c = calloc((num_procs_sqrt + 1), sizeof(int));

    for (int i = 0; i < num_procs_sqrt; i++){
      row_cnt_a[i] = rows_a/num_procs_sqrt + (i < rows_a % num_procs_sqrt);
      row_cnt_b[i] = rows_b/num_procs_sqrt + (i < rows_b % num_procs_sqrt);
      row_cnt_c[i] = rows_c/num_procs_sqrt + (i < rows_c % num_procs_sqrt);
      col_cnt_a[i] = cols_a/num_procs_sqrt + (i < cols_a % num_procs_sqrt);
      col_cnt_b[i] = cols_b/num_procs_sqrt + (i < cols_b % num_procs_sqrt);
      col_cnt_c[i] = cols_c/num_procs_sqrt + (i < cols_c % num_procs_sqrt);

      row_displ_a[i+1] = row_displ_a[i] + row_cnt_a[i];
      row_displ_b[i+1] = row_displ_b[i] + row_cnt_b[i];
      row_displ_c[i+1] = row_displ_c[i] + row_cnt_c[i];
      col_displ_a[i+1] = col_displ_a[i] + col_cnt_a[i];
      col_displ_b[i+1] = col_displ_b[i] + col_cnt_b[i];
      col_displ_c[i+1] = col_displ_c[i] + col_cnt_c[i];
    }
  }

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

  //calculate matrix max
  A_part_m_max = rows_a / num_procs_sqrt + !!(cols_c % num_procs_sqrt);
  A_part_n_max = cols_a / num_procs_sqrt + !!(cols_c % num_procs_sqrt);
  B_part_m_max = rows_b / num_procs_sqrt + !!(cols_c % num_procs_sqrt);
  B_part_n_max = cols_b / num_procs_sqrt + !!(cols_c % num_procs_sqrt);
  C_part_m_max = rows_c / num_procs_sqrt + !!(cols_c % num_procs_sqrt);
  C_part_n_max = cols_c / num_procs_sqrt + !!(cols_c % num_procs_sqrt);

  //find partion sizes
  rows_apart = rows_a / num_procs_sqrt + (rankrows < rows_a % num_procs_sqrt);
  cols_apart = cols_a / num_procs_sqrt + (rankcols < cols_a % num_procs_sqrt);
  rows_bpart = rows_b / num_procs_sqrt + (rankrows < rows_b % num_procs_sqrt);
  cols_bpart = cols_b / num_procs_sqrt + (rankcols < cols_b % num_procs_sqrt);
  rows_cpart = rows_c / num_procs_sqrt + (rankrows < rows_c % num_procs_sqrt);
  cols_cpart = cols_c / num_procs_sqrt + (rankcols < cols_c % num_procs_sqrt);

  //allocate partitions
  allocate_matrix(&A_part, A_part_m_max, A_part_n_max);
  allocate_matrix(&B_part, B_part_m_max, B_part_n_max);
  allocate_matrix(&C_part, C_part_m_max, C_part_n_max);

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
      moff = row_displ_a[i / num_procs_sqrt];
      noff = col_displ_a[i % num_procs_sqrt];
      for (int j = 0; j < Am; j++) {
        MPI_Send(&(matrix_a[j + moff][noff]), An, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
      }
      moff = row_displ_b[i / num_procs_sqrt];
      noff = col_displ_b[i % num_procs_sqrt];
      for (int j = 0; j <  Bm; j++){
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

  //Cannon_mult

	// Set up the Cartesian topology
  dims[0] = dims[1] = sqrt(num_procs);

	// Set the periods for wraparound connections, 1 == true
  periods[0] = periods[1] = 1;

	// Get the rank and coordinates with respect to the new topology
  MPI_Comm_rank(comm_2d, &rank2d);
  MPI_Cart_coords(comm_2d, rank2d, 2, mycoords);

  //Compute ranks of the up and left shifts
  MPI_Cart_shift(comm_2d, 1, -1, &rightrank, &leftrank);
  MPI_Cart_shift(comm_2d, 0, -1, &downrank, &uprank);

  //Perform the initial matrix alignment. First for A and then for B
  MPI_Cart_shift(comm_2d, 1, -mycoords[0], &shiftsource, &shiftdest);
  MPI_Sendrecv_replace(&rows_apart, 1, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, MPI_STATUS_IGNORE);
  MPI_Sendrecv_replace(&cols_apart,1, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, MPI_STATUS_IGNORE);
  for (int i = 0; i < A_part_m_max; i++)
    MPI_Sendrecv_replace(A_part[i], A_part_n_max, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_2d, MPI_STATUS_IGNORE);

  MPI_Cart_shift(comm_2d, 0, -mycoords[1], &shiftsource, &shiftdest);
  MPI_Sendrecv_replace(&rows_bpart, 1, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, MPI_STATUS_IGNORE);
  MPI_Sendrecv_replace(&cols_bpart, 1, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, MPI_STATUS_IGNORE);
  for (int i = 0; i < B_part_m_max; i++)
    MPI_Sendrecv_replace(B_part[i], B_part_n_max, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_2d, MPI_STATUS_IGNORE);

  // Get into the main computation loop
  for (int i = 0; i < dims[0]; i++) {
    matrix_mult(A_part, B_part, C_part, rows_apart, cols_apart, cols_bpart);

    // Shift matrix a left by one
    MPI_Sendrecv_replace(&rows_apart, 1, MPI_INT, leftrank, 1, rightrank, 1, comm_2d, MPI_STATUS_IGNORE);
    MPI_Sendrecv_replace(&cols_apart, 1, MPI_INT, leftrank, 1, rightrank, 1, comm_2d, MPI_STATUS_IGNORE);
    for (int j = 0; j < A_part_m_max; j++)
      MPI_Sendrecv_replace(A_part[j], A_part_n_max, MPI_DOUBLE, leftrank, 1, rightrank, 1, comm_2d, MPI_STATUS_IGNORE);

	  // Shift matrix b up by one
    MPI_Sendrecv_replace(&rows_bpart, 1, MPI_INT, uprank, 1, downrank, 1, comm_2d, MPI_STATUS_IGNORE);
    MPI_Sendrecv_replace(&cols_bpart, 1, MPI_INT, uprank, 1, downrank, 1, comm_2d, MPI_STATUS_IGNORE);
    for (int j = 0; j < B_part_m_max; j++)
      MPI_Sendrecv_replace(B_part[j], B_part_n_max, MPI_DOUBLE, uprank, 1, downrank, 1, comm_2d, MPI_STATUS_IGNORE);
  }

  //gather
  if (my_rank == 0) {
    for (int i = 0; i < rows_cpart; i++){
      memcpy(matrix_c[i], C_part[i], cols_cpart * sizeof(double));
    }

    for (int i = 1; i < num_procs; i++) {
      int Cm, Cn;
      MPI_Recv(&Cm, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&Cn, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      int moff = row_displ_c[i / num_procs_sqrt];
      int noff = col_displ_c[i % num_procs_sqrt];
      for (int j = 0; j < row_cnt_c[i / num_procs_sqrt]; j++) {
        MPI_Recv(&(matrix_c[j + moff][noff]), Cn, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
  } else{
    MPI_Send(&rows_cpart, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&cols_cpart, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

    for (int j = 0; j < rows_cpart; j++){
      MPI_Send(C_part[j], cols_cpart, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
  }

  //save result
  if (my_rank == 0){
    write_matrix_binaryformat(argv[3], matrix_c, rows_c, cols_c);
  }

  // Free up communicator
  if (my_rank == 0)
    deallocate_matrix(matrix_c);
  deallocate_matrix(A_part);
  deallocate_matrix(B_part);
  deallocate_matrix(C_part);
  MPI_Comm_free(&comm_2d);
  free(row_cnt_a);
  free(row_cnt_b);
  free(row_cnt_c);
  free(col_cnt_a);
  free(col_cnt_b);
  free(col_cnt_c);
  free(row_displ_a);
  free(row_displ_b);
  free(row_displ_c);
  free(col_displ_a);
  free(col_displ_b);
  free(col_displ_c);

  //cleanup
  MPI_Finalize();
  return 0;
}
