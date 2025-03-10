#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <mpi.h>
#include "adi.h"

double second();

DATA_TYPE X[N][N];
DATA_TYPE A[N][N];
DATA_TYPE B[N][N];

static void init_array(int n, DATA_TYPE X[N][N], DATA_TYPE A[N][N], DATA_TYPE B[N][N]) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            X[i][j] = ((DATA_TYPE) i*(j+1) + 1) / n;
            A[i][j] = ((DATA_TYPE) i*(j+2) + 2) / n;
            B[i][j] = ((DATA_TYPE) i*(j+3) + 3) / n;
        }
}

static void print_array(int n, DATA_TYPE X[N][N]) {
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            printf(DATA_PRINTF_MODIFIER, X[i][j]);
        }
        printf("\n");
    }
}

static void kernel_adi(int tsteps, int rank, int size) {
    int n = N;
    int chunk_size = (n + size - 1) / size;
    int start_row = rank * chunk_size;
    int end_row = (start_row + chunk_size > n) ? n : (start_row + chunk_size);

    for (int t = 0; t < _PB_TSTEPS; t++) {
        // Этап 1: Обновление X и B по строкам
        for (int i1 = start_row; i1 < end_row; i1++) {
            for (int i2 = 1; i2 < _PB_N; i2++) {
                X[i1][i2] = X[i1][i2] - X[i1][i2-1] * A[i1][i2] / B[i1][i2-1];
                B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1][i2-1];
            }
        }

        // Синхронизация данных между процессами
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, X, chunk_size * N, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, B, chunk_size * N, MPI_DOUBLE, MPI_COMM_WORLD);

        // Этап 2: Нормализация последнего столбца X
        for (int i1 = start_row; i1 < end_row; i1++) {
            X[i1][_PB_N-1] = X[i1][_PB_N-1] / B[i1][_PB_N-1];
        }

        // Синхронизация данных между процессами
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, X, chunk_size * N, MPI_DOUBLE, MPI_COMM_WORLD);

        // Этап 3: Обратное обновление X
        for (int i1 = start_row; i1 < end_row; i1++) {
            for (int i2 = 0; i2 < _PB_N-2; i2++) {
                X[i1][_PB_N-i2-2] = (X[i1][_PB_N-2-i2] - X[i1][_PB_N-2-i2-1] * A[i1][_PB_N-i2-3]) / B[i1][_PB_N-3-i2];
            }
        }

        // Синхронизация данных между процессами
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, X, chunk_size * N, MPI_DOUBLE, MPI_COMM_WORLD);

        // Этап 4: Обновление X и B по столбцам
        for (int i1 = 1; i1 < _PB_N; i1++) {
            for (int i2 = 0; i2 < _PB_N; i2++) {
                X[i1][i2] = X[i1][i2] - X[i1-1][i2] * A[i1][i2] / B[i1-1][i2];
                B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1-1][i2];
            }
        }

        // Синхронизация данных между процессами
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, X, chunk_size * N, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, B, chunk_size * N, MPI_DOUBLE, MPI_COMM_WORLD);

        // Этап 5: Нормализация последней строки X
        if (rank == size - 1) { // Только последний процесс обрабатывает последнюю строку
            for (int i2 = 0; i2 < _PB_N; i2++) {
                X[_PB_N-1][i2] = X[_PB_N-1][i2] / B[_PB_N-1][i2];
            }
        }

        // Синхронизация данных между процессами
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, X, chunk_size * N, MPI_DOUBLE, MPI_COMM_WORLD);

        // Этап 6: Обратное обновление X по строкам
        for (int i1 = 0; i1 < _PB_N-2; i1++) {
            for (int i2 = 0; i2 < _PB_N; i2++) {
                X[_PB_N-2-i1][i2] = (X[_PB_N-2-i1][i2] - X[_PB_N-i1-3][i2] * A[_PB_N-3-i1][i2]) / B[_PB_N-2-i1][i2];
            }
        }

        // Синхронизация данных между процессами
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, X, chunk_size * N, MPI_DOUBLE, MPI_COMM_WORLD);
    }
}

int main(int argc, char** argv) {
    int n = N;
    int tsteps = TSTEPS;
    int rank, size;
    double time0, time1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        init_array(n, X, A, B);
    }
    MPI_Bcast(X, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(A, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(B, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    time0 = MPI_Wtime();
    
    kernel_adi(tsteps, rank, size);
    
    MPI_Barrier(MPI_COMM_WORLD);
    time1 = MPI_Wtime();

    if (rank == 0) {
        printf("\nn=%d\n", n);
        printf("\ntime=%f\n", time1 - time0);
        print_array(n, X);
    }
    
    MPI_Finalize();
    return 0;
}
