#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <mpi.h>

#include "adi.h"

double second();

// Удаляем статические массивы и заменяем их на динамические
DATA_TYPE **X;
DATA_TYPE **A;
DATA_TYPE **B;
int local_rows;

// Функция инициализации массивов X, A и B
static void init_array(int n, DATA_TYPE **X, DATA_TYPE **A, DATA_TYPE **B, int rank, int size)
{
    int i, j;
    // int chunk_rows = n / size;     // Количество строк, обрабатываемых каждым процессом
    // int chunk_cols = n / size;     // Количество столбцов, обрабатываемых каждым процессом
    // int remainder_rows = n % size; // Остаток строк
    // int remainder_cols = n % size; // Остаток столбцов

    // int start_row = rank * chunk_rows + (rank < remainder_rows ? rank : remainder_rows);
    // int start_col = rank * chunk_cols + (rank < remainder_cols ? rank : remainder_cols);
    // int local_rows = chunk_rows + (rank < remainder_rows ? 1 : 0);
    // int local_cols = chunk_cols + (rank < remainder_cols ? 1 : 0);

    // int end_row = start_row + chunk_rows + (rank < remainder_rows ? 1 : 0);
    // Инициализируем динамические массивы
    // Инициализация локальной части массивов
    for (i = 0; i < local_rows; i++)
    {
        for (j = 0; j < n; j++)
        {
            X[i][j] = ((DATA_TYPE)(rank * local_rows + i) * (j + 1) + 1) / n;
            A[i][j] = ((DATA_TYPE)(rank * local_rows + i) * (j + 2) + 2) / n;
            B[i][j] = ((DATA_TYPE)(rank * local_rows + i) * (j + 3) + 3) / n;
        }
    }
}

static void print_row(int row)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    row--;

    int valid_rank = row / local_rows;
    row = row % local_rows;

    if (rank == valid_rank)
    {
        printf("Rank %d, row %d:\n", valid_rank, row);

        for (int i = 0; i < N; i++)
            printf("%f \n", X[row][i]);
    }
}

// Функция для вывода массива X
static void print_array(int n, DATA_TYPE **X, int rank, int size)
{
    int chunk = n / size;
    int remainder = n % size;
    int start_i1 = rank * chunk + (rank < remainder ? rank : remainder);
    int end_i1 = start_i1 + chunk + (rank < remainder ? 1 : 0);
    int local_rows = end_i1 - start_i1;
    int row = 16;

    // int i, j;
    // int chunk_rows = n / size; // Количество строк, обрабатываемых каждым процессом
    // int chunk_cols = n / size;
    // int remainder_rows = n % size;// Остаток строк
    // int remainder_cols = n % size;

    // int start_row = rank * chunk_rows + (rank < remainder_rows ? rank : remainder_rows);
    // int start_col = rank * chunk_cols + (rank < remainder_cols ? rank : remainder_cols);
    // int local_rows = chunk_rows + (rank < remainder_rows ? 1 : 0);
    // int local_cols = chunk_cols + (rank < remainder_cols ? 1 : 0);

    // for (i = start_row; i < local_rows; i++) {
    //     for (j = start_col; j < local_cols; j++) {
    //         printf("%f ", X[i][j]);
    //     }
    //     printf("\n");
    // }
    // Печатаем только из процесса 0 после сбора всех данных

    // Для простоты оставим вывод как было, но нужно будет собрать данные на процессе 0
    //    for (int i = 10; i < 11; i++) {
    //        for (int j = 0; j < n; j++) {
    //            printf("%f\n ", X[i][j]);
    //          //  printf("Rank %d: start_i1 = %d, end_i1 = %d, local_rows = %d\n", rank, start_i1, end_i1, local_rows);
    //        }
    //        printf("\n");
    //    }
    if (row >= start_i1 && row < end_i1)
    {
        int local_i = row - start_i1; // Преобразуем глобальный индекс в локальный
        printf("Process %d, row %d:\n", rank, row);
        for (int j = 0; j < n; j++)
        {
            printf("%f ", X[local_i][j]);
            printf("\n");
        }
        printf("\n");
    }
}
static void kernel_adi(int tsteps, int n, DATA_TYPE **X, DATA_TYPE **A, DATA_TYPE **B, int rank, int size, int local_rows)
{
    int t, i1, i2;

    for (t = 0; t < tsteps; t++)
    {
        // Горизонтальные обновления
        for (i1 = 0; i1 < local_rows; i1++)
        {
            for (i2 = 1; i2 < n; i2++)
            {
                X[i1][i2] = X[i1][i2] - X[i1][i2 - 1] * A[i1][i2] / B[i1][i2 - 1];
                B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1][i2 - 1];
                // printf("Rank %d: start_i1 = %d, end_i1 = %d, local_rows = %d, current_row = %d\n", rank, start_i1, end_i1, local_rows, i1);
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }
    

    // Обновление X[i1][n-1] для локальных строк
    for ( i1 = 0; i1 < local_rows; i1++) 
       X[i1][n - 1] = X[i1][n - 1] / B[i1][n - 1];
    
    
    // Обратная подстановка для локальных строк
    for (int i1 = 0; i1 < local_rows; i1++) {
        for (i2 = 0; i2 < n - 2; i2++) {
            int idx = n - i2 - 2;
            X[i1][idx] = (X[i1][idx] - X[i1][idx - 1] * A[i1][idx - 1]) / B[i1][idx - 1];
        }
    }

    // // Вертикальные обновления с обменом данными между процессами
    // for (int i1 = 0; local_i1 < local_rows; local_i1++) {
    //     int global_i1 = start_i1 + local_i1;  // Преобразование локального индекса в глобальный

    //     if (global_i1 == 0) continue; // Пропустить первую строку

    //     // Если это не первый процесс, получаем данные от предыдущего процесса
    //     if (rank != 0 && local_i1 == 0) {
    //         printf("Rank %d: Receiving data from rank %d\n", rank, rank - 1);
    //         MPI_Recv(X[local_i1 - 1], n, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         MPI_Recv(B[local_i1 - 1], n, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     }

    //     for (i2 = 0; i2 < n; i2++) {
    //         X[local_i1][i2] = X[local_i1][i2] - X[local_i1 - 1][i2] * A[local_i1][i2] / B[local_i1 - 1][i2];
    //         B[local_i1][i2] = B[local_i1][i2] - A[local_i1][i2] * A[local_i1][i2] / B[local_i1 - 1][i2];
    //     }

    //     // Если это не последний процесс, отправляем данные следующему процессу
    //     if (rank != size - 1 && local_i1 == local_rows - 1) {
    //         printf("Rank %d: Sending data to rank %d\n", rank, rank + 1);
    //         MPI_Send(X[local_i1], n, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    //         MPI_Send(B[local_i1], n, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
    //     }
    // }

    // // Синхронизация процессов перед обновлением последней строки
    // MPI_Barrier(MPI_COMM_WORLD);

    // // Обновление последней строки только последним процессом
    // if (rank == size - 1) {
    //     for (i2 = 0; i2 < n; i2++) {
    //         X[local_rows - 1][i2] = X[local_rows - 1][i2] / B[local_rows - 1][i2];
    //     }
    // }

    // // Обратная подстановка для вертикальных обновлений с обменом данными
    // for (i1 = 0; i1 < n - 2; i1++) {
    //     int idx = n - i1 - 2;

    //     // Если это не первый процесс, получаем данные от предыдущего процесса
    //     if (rank != 0 && idx == start_i1) {
    //         printf("Rank %d: Receiving data from rank %d\n", rank, rank - 1);
    //         MPI_Recv(X[idx - start_i1 - 1], n, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         MPI_Recv(B[idx - start_i1 - 1], n, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     }

    //     for (i2 = 0; i2 < n; i2++) {
    //         X[idx - start_i1][i2] = (X[idx - start_i1][i2] - X[idx - start_i1 - 1][i2] * A[idx - start_i1 - 1][i2]) / B[idx - start_i1][i2];
    //     }

    //     // Если это не последний процесс, отправляем данные следующему процессу
    //     if (rank != size - 1 && idx == end_i1 - 1) {
    //         printf("Rank %d: Sending data to rank %d\n", rank, rank + 1);
    //         MPI_Send(X[idx - start_i1], n, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    //         MPI_Send(B[idx - start_i1], n, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
    //     }
    //}
    // }
}
}

int main(int argc, char **argv)
{
    int n = N;                 // Размер массива
    int tsteps = TSTEPS;       // Количество временных шагов
    int rank, size;            // Ранг и количество процессов
    double time, time0, time1; // Время выполнения

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Выделение памяти для локальной части массивов
    int chunk_rows = n / size;
    int chunk_cols = n / size;
    int remainder_rows = n % size;
    int remainder_cols = n % size;

    local_rows = chunk_rows + (rank < remainder_rows ? 1 : 0);
    int local_cols = chunk_cols + (rank < remainder_cols ? 1 : 0);

    X = (DATA_TYPE **)malloc(local_rows * sizeof(DATA_TYPE *));
    A = (DATA_TYPE **)malloc(local_rows * sizeof(DATA_TYPE *));
    B = (DATA_TYPE **)malloc(local_rows * sizeof(DATA_TYPE *));

    for (int i = 0; i < local_rows; i++)
    {
        X[i] = (DATA_TYPE *)malloc(n * sizeof(DATA_TYPE));
        A[i] = (DATA_TYPE *)malloc(n * sizeof(DATA_TYPE));
        B[i] = (DATA_TYPE *)malloc(n * sizeof(DATA_TYPE));
    }

    init_array(n, X, A, B, rank, size);

    // Замер времени выполнения
    time0 = MPI_Wtime();
    kernel_adi(tsteps, n, X, A, B, rank, size, local_rows);
    time1 = MPI_Wtime();

    time = time1 - time0;

    if (rank == 0)
    {
        printf("\nn=%d\n", n);
        printf("\n\n\ntime=%f\n", time);
        printf("Total execution time: %f seconds\n", time);
    }

    print_row(30);
    
    // print_array(n, X, rank, size);

    for (int i = 0; i < local_rows; i++)
    {
        free(X[i]);
        free(A[i]);
        free(B[i]);
    }
    free(X);
    free(A);
    free(B);

    MPI_Finalize();

    return 0;
}

double second()
{
    struct timeval tm;
    double t;

    static int base_sec = 0, base_usec = 0;

    gettimeofday(&tm, NULL);

    if (base_sec == 0 && base_usec == 0)
    {
        base_sec = tm.tv_sec;
        base_usec = tm.tv_usec;
        t = 0.0;
    }
    else
    {
        t = (double)(tm.tv_sec - base_sec) + (double)(tm.tv_usec - base_usec) / 1.0e6;
    }

    return t;
}
