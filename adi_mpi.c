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
    int i, j;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            X[i][j] = ((DATA_TYPE) i*(j+1) + 1) / n;
            A[i][j] = ((DATA_TYPE) i*(j+2) + 2) / n;
            B[i][j] = ((DATA_TYPE) i*(j+3) + 3) / n;
        }
}

static void print_array(int n, DATA_TYPE X[N][N]) {
    int i, j;

    for (i = 1; i < 2; i++) {
        for (j = 2020; j < 2048; j++) {
            printf("%f\n ", X[i][j]);
        }
        printf("\n");
    }
}

static void kernel_adi(int tsteps, int n, DATA_TYPE X[N][N], DATA_TYPE A[N][N], DATA_TYPE B[N][N], int rank, int size) {
    int t, i1, i2;

    // Определение диапазона строк для текущего процесса
    int chunk = n / size;// Количество строк на процесс
    int remainder = n % size; // Остаток строк, если n не делится нацело на size
    int start_i1 = rank * chunk + (rank < remainder ? rank : remainder); // Начальная строка для текущего процесса
    int end_i1 = start_i1 + chunk + (rank < remainder ? 1 : 0);// Конечная строка для текущего процесса
    int local_rows = end_i1 - start_i1;// Количество строк, обрабатываемых текущим процессом

    // Буферы для передачи граничных строк
    DATA_TYPE *send_X = (DATA_TYPE *)malloc(n * sizeof(DATA_TYPE));
    DATA_TYPE *send_B = (DATA_TYPE *)malloc(n * sizeof(DATA_TYPE));
    DATA_TYPE *recv_X = (DATA_TYPE *)malloc(n * sizeof(DATA_TYPE));
    DATA_TYPE *recv_B = (DATA_TYPE *)malloc(n * sizeof(DATA_TYPE));

    // Подготовка параметров для Allgatherv
    int *sendcounts = (int *)malloc(size * sizeof(int));// Количество элементов для отправки каждым процессом
    int *displs = (int *)malloc(size * sizeof(int));// Смещения для каждого процесса в общем массиве
    for (int i = 0; i < size; ++i) {
        int s = i * chunk + (i < remainder ? i : remainder);// Начальная строка для процесса i
        int e = s + chunk + (i < remainder ? 1 : 0);// Конечная строка для процесса i
        sendcounts[i] = (e - s) * n; // Количество элементов для процесса i
        displs[i] = s * n; // Смещение для процесса i
    }

    for (t = 0; t < tsteps; t++) {
        // Горизонтальные обновления
        for (i1 = start_i1; i1 < end_i1; i1++) {
            for (i2 = 1; i2 < n; i2++) {
                X[i1][i2] = X[i1][i2] - X[i1][i2-1] * A[i1][i2] / B[i1][i2-1];
                B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1][i2-1];
            }
        }

        // Обновление X[i1][n-1] для локальных строк
        for (i1 = start_i1; i1 < end_i1; i1++) {
            X[i1][n-1] = X[i1][n-1] / B[i1][n-1];
        }

        // Обратная подстановка для локальных строк
        for (i1 = start_i1; i1 < end_i1; i1++) {
            for (i2 = 0; i2 < n-2; i2++) {
                int idx = n - i2 - 2;
                X[i1][idx] = (X[i1][idx] - X[i1][idx - 1] * A[i1][idx - 1]) / B[i1][idx -1];
            }
        }

        // Обмен граничными строками для вертикальных обновлений
        if (size > 1) {
            MPI_Request reqs[4];
            if (rank != 0) {
                MPI_Irecv(recv_X, n, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &reqs[0]);
                MPI_Irecv(recv_B, n, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &reqs[1]);
            }
            if (rank != size-1) {
                memcpy(send_X, X[end_i1-1], n*sizeof(DATA_TYPE));
                memcpy(send_B, B[end_i1-1], n*sizeof(DATA_TYPE));
                MPI_Isend(send_X, n, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &reqs[2]);
                MPI_Isend(send_B, n, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &reqs[3]);
            }

            // Ожидание завершения операций приема и отправки
            if (rank != 0) MPI_Waitall(2, reqs, MPI_STATUS_IGNORE);
            if (rank != size-1) MPI_Waitall(2, &reqs[2], MPI_STATUS_IGNORE);

            if (rank != 0) {
                // Копирование полученных данных в локальные массивы
                memcpy(X[start_i1-1], recv_X, n*sizeof(DATA_TYPE));
                memcpy(B[start_i1-1], recv_B, n*sizeof(DATA_TYPE));
            }
        }

        // Вертикальные обновления
        for (i1 = start_i1; i1 < end_i1; i1++) {
            if (i1 == 0) continue; // Пропустить первую строку
            for (i2 = 0; i2 < n; i2++) {
                X[i1][i2] = X[i1][i2] - X[i1-1][i2] * A[i1][i2] / B[i1-1][i2];
                B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1-1][i2];
            }
        }

        for (i2 = 0; i2 < n; i2++) {
            X[n-1][i2] = X[n-1][i2] / B[n-1][i2];
        }

        // Обратная подстановка для вертикальных обновлений
        for (i1 = 0; i1 < n-2; i1++) {
            for (i2 = 0; i2 < n; i2++) {
                int idx = n - i1 - 2;
                X[idx][i2] = (X[idx][i2] - X[idx - 1][i2] * A[idx - 1][i2]) / B[idx][i2];
            }
        }

        // Синхронизация данных между всеми процессами
        MPI_Allgatherv(MPI_IN_PLACE, local_rows*n, MPI_DOUBLE, X, sendcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(MPI_IN_PLACE, local_rows*n, MPI_DOUBLE, B, sendcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
    }

    free(send_X);
    free(send_B);
    free(recv_X);
    free(recv_B);
    free(sendcounts);
    free(displs);
}

int main(int argc, char** argv) {
    int n = N;
    int tsteps = TSTEPS;
    int rank, size;
    double time, time0, time1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        init_array(n, X, A, B);
    }

    MPI_Bcast(X, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(A, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(B, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    time0 = MPI_Wtime();
    kernel_adi(tsteps, n, X, A, B, rank, size);
    time1 = MPI_Wtime();

    time = time1 - time0;

    if (rank == 0) {
        printf("\nn=%d\n", n);
        printf("\n\n\ntime=%f\n", time);
        print_array(n, X);
        printf("Total execution time: %f seconds\n", time);
    }

    MPI_Finalize();

    return 0;
}

double second() {
    struct timeval tm;
    double t;

    static int base_sec = 0, base_usec = 0;

    gettimeofday(&tm, NULL);

    if (base_sec == 0 && base_usec == 0) {
        base_sec = tm.tv_sec;
        base_usec = tm.tv_usec;
        t = 0.0;
    } else {
        t = (double)(tm.tv_sec - base_sec) + (double)(tm.tv_usec - base_usec) / 1.0e6;
    }

    return t;
}
