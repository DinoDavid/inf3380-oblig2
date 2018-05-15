#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

struct matrix {
        double **data;
        int rowsiz;
        int colsiz;
        int m;
        int n;
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
assertequal(float a, float b, int i, int j)
{
        float eps = 0.1;
        float diff;

        diff = a - b;
        diff = (diff < 0) ? -diff : diff;
        if (diff > eps || isnan(diff) || isinf(diff))
                die("not equal by %f @%d,%d\n", a - b, i, j);
}

int
main(int argc, char **argv)
{
        struct matrix A;
        struct matrix B;
        int i;
        int j;

        if (argc != 3)
                die("usage: %s <mat1> <mat2>\n", argv[0]);

        readmat(argv[1], &A);
        readmat(argv[2], &B);

        for (i = 0; i < A.m; i++)
                for (j = 0; j < A.n; j++)
                        assertequal(A.data[i][j], B.data[i][j], i, j);

        printf("comparison passed\n");

        return 0;
}
