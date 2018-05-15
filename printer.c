#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lib3380.h"

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
