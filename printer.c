#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct matrix {
        double **data;
        int rowsiz;
        int colsiz;
        int m;
        int n;
};

struct partdat {
        int *rowcnt;
        int *colcnt;
        int *rowdspl;
        int *coldspl;
};

void
die(char *msg, ...)
{
        va_list ap;

        va_start(ap, msg);
        vfprintf(stderr, msg, ap);
        va_end(ap);

        exit(EXIT_FAILURE);
}

void
allocmat(int rowsiz, int colsiz, int m, int n, struct matrix *mat)
{
        int i;

        /*
         * TODO check m <= rowsiz etc
         */

        mat->rowsiz = rowsiz;
        mat->colsiz = colsiz;
        mat->m = m;
        mat->n = n;

        mat->data = malloc(rowsiz * sizeof(*mat->data));
        *mat->data = malloc(rowsiz * colsiz * sizeof(**mat->data));
        for (i = 0; i < rowsiz; i++)
                *(mat->data + i) = (*mat->data) + i * colsiz;
}

void
deallocmat(struct matrix *mat)
{
        mat->rowsiz = 0;
        mat->colsiz = 0;
        mat->m = 0;
        mat->n = 0;
        free(*mat->data);
        free(mat->data);
        mat->data = 0;
}

void
readmat(char *file, struct matrix *mat)
{
        FILE *fp;
        int rows;
        int cols;

        fp = fopen(file, "r");
        if (!fp)
                die("failed to open '%s' for reading\n", file);

        fread(&rows, sizeof(rows), 1, fp);
        fread(&cols, sizeof(cols), 1, fp);
        allocmat(rows, cols, rows, cols, mat);
        fread(*mat->data, sizeof(**mat->data), mat->m * mat->n, fp);

        fclose(fp);
}

void
calclens(int **lens, int **displs, int procswid, int fulllen)
{
        int i;
        int minlen;
        int remain;

        *lens = malloc(procswid * sizeof(*lens));
        *displs = malloc((procswid + 1) * sizeof(*displs));

        minlen = fulllen / procswid;
        remain = fulllen % procswid;
        (*displs)[0] = 0;
        for (i = 0; i < procswid; i++) {
                (*lens)[i] = minlen + (i < remain);
                (*displs)[i + 1] = (*displs)[i] + (*lens)[i];
        }
}

void
printmat(struct matrix *A)
{
        int i;
        int j;

        printf("rows: %d\n", A->m);
        printf("cols: %d\n", A->n);

        for (i = 0; i < A->m; i++)
                for (j = 0; j < A->n; j++)
                        printf("[%d][%d] = %f\n", i, j, A->data[i][j]);
}

void
usage(char *argv0)
{
        die("usage: %s matrix full|(part np)\n", argv0);
}

void
printpartition(struct matrix *M, int mp, int np, int moff, int noff)
{
        int i;
        int j;

        for (i = 0; i < mp; i++)
                for (j = 0; j < np; j++)
                        printf("%f\n", M->data[i + moff][j + noff]);
}

void
printparts(struct matrix *M, int meshwid)
{
        int *rowcnt;
        int *colcnt;
        int *rowdspl;
        int *coldspl;
        int pi;
        int pj;

        calclens(&rowcnt, &rowdspl, meshwid, M->m);
        calclens(&colcnt, &coldspl, meshwid, M->n);
        for (pi = 0; pi < meshwid; pi++) {
                for (pj = 0; pj < meshwid; pj++) {
                        printf("partition %d\n", pi * meshwid + pj);

                        printpartition(M, rowcnt[pi], colcnt[pj],
                                       rowdspl[pi], coldspl[pj]);
                }
        }

        printf("%d\n", meshwid);
}

int
main(int argc, char **argv)
{
        struct matrix A;

        if (argc != 3 && argc != 4)
                usage(argv[0]);

        readmat(argv[1], &A);
        if (!strcmp(argv[2], "full"))
                printmat(&A);
        else if (argc == 4 && !strcmp(argv[2], "part"))
                printparts(&A, sqrt(atoi(argv[3])));
        else
                usage(argv[0]);

        deallocmat(&A);
        return EXIT_SUCCESS;
}
