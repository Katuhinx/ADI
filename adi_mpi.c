#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <mpi.h>

#include "adi.h"

double second();

DATA_TYPE **X;
DATA_TYPE **A;
DATA_TYPE **B;
int local_rows;
int num_quanta = 1;
int quantum_size;

// Функция инициализации массивов X, A и B
static void init_arrays()
{
    int i, j;
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    X = (DATA_TYPE **)malloc(local_rows * sizeof(DATA_TYPE *));
    A = (DATA_TYPE **)malloc(local_rows * sizeof(DATA_TYPE *));
    B = (DATA_TYPE **)malloc(local_rows * sizeof(DATA_TYPE *));

    for (int i = 0; i < local_rows; i++)
    {
        X[i] = (DATA_TYPE *)malloc(N * sizeof(DATA_TYPE));
        A[i] = (DATA_TYPE *)malloc(N * sizeof(DATA_TYPE));
        B[i] = (DATA_TYPE *)malloc(N * sizeof(DATA_TYPE));
    }

    for (i = 0; i < local_rows; i++)
    {
        for (j = 0; j < N; j++)
        {
            X[i][j] = ((DATA_TYPE)(rank * local_rows + i) * (j + 1) + 1) / N;
            A[i][j] = ((DATA_TYPE)(rank * local_rows + i) * (j + 2) + 2) / N;
            B[i][j] = ((DATA_TYPE)(rank * local_rows + i) * (j + 3) + 3) / N;
        }
    }
}

// Функция освобождения памяти массивов X, A и B
static void free_arrays()
{
    for (int i = 0; i < local_rows; i++)
    {
        free(X[i]);
        free(A[i]);
        free(B[i]);
    }

    free(X);
    free(A);
    free(B);
}

// Функция вывода результирующих значений из строки row массива X
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
    int t, i1, i2, q;
    DATA_TYPE *quantum_X = (DATA_TYPE *)malloc(quantum_size * sizeof(DATA_TYPE));
    DATA_TYPE *quantum_B = (DATA_TYPE *)malloc(quantum_size * sizeof(DATA_TYPE));
    DATA_TYPE *temp_X = (DATA_TYPE *)malloc(N * sizeof(DATA_TYPE));
    DATA_TYPE *temp_A = (DATA_TYPE *)malloc(N * sizeof(DATA_TYPE));
    MPI_Status status;

    for (t = 0; t < tsteps; t++)
    {
        // Горизонтальные обновления
        for (i1 = 0; i1 < local_rows; i1++)
        {
            for (i2 = 1; i2 < n; i2++)
            {
                X[i1][i2] = X[i1][i2] - X[i1][i2 - 1] * A[i1][i2] / B[i1][i2 - 1];
                B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1][i2 - 1];
            }

            X[i1][n - 1] = X[i1][n - 1] / B[i1][n - 1];
        }

        // Обратная подстановка для локальных строк
        for (int i1 = 0; i1 < local_rows; i1++)
        {
            for (i2 = 0; i2 < n - 2; i2++)
            {
                int idx = n - i2 - 2;
                X[i1][idx] = (X[i1][idx] - X[i1][idx - 1] * A[i1][idx - 1]) / B[i1][idx - 1];
            }
        }

        // Вертикальные обновления с конвейерной обработкой
        for (q = 0; q < num_quanta; q++)
        {
            int start = q * quantum_size;
            int end = start + quantum_size;

            end = (end > N) ? N : end;

            // Получение данных от предыдущего процесса
            if (rank > 0)
            {
                MPI_Recv(quantum_X, quantum_size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
                MPI_Recv(quantum_B, quantum_size, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
            }

            // Обработка текущего кванта
            for (i1 = 0; i1 < local_rows; i1++)
            {
                if (i1 == 0)
                {
                    if (rank > 0)
                    {
                        for (i2 = start; i2 < end; i2++)
                        {
                            X[i1][i2] = X[i1][i2] - quantum_X[i2 - start] * A[i1][i2] / quantum_B[i2 - start];
                            B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / quantum_B[i2 - start];
                        }
                    }
                }
                else
                {
                    for (i2 = start; i2 < end; i2++)
                    {
                        X[i1][i2] = X[i1][i2] - X[i1 - 1][i2] * A[i1][i2] / B[i1 - 1][i2];
                        B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1 - 1][i2];
                    }
                }
            }

            // Отправка данных следующему процессу
            if (rank < size - 1)
            {
                for (i2 = start; i2 < end; i2++)
                {
                    quantum_X[i2 - start] = X[local_rows - 1][i2];
                    quantum_B[i2 - start] = B[local_rows - 1][i2];
                }
                MPI_Send(quantum_X, quantum_size, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
                MPI_Send(quantum_B, quantum_size, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
            }
        }

        // Обновление последней строки
        if (rank == size - 1)
        {
            for (i2 = 0; i2 < n; i2++)
            {
                X[local_rows - 1][i2] = X[local_rows - 1][i2] / B[local_rows - 1][i2];
            }
        }

        // Обратная подстановка для вертикальных обновлений с обменом данными
        if (size > 1 && rank > 0)
        {
            MPI_Recv(temp_X, n, MPI_DOUBLE, rank - 1, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(temp_A, n, MPI_DOUBLE, rank - 1, 3, MPI_COMM_WORLD, &status);
        }

        for (i1 = local_rows - 1; i1 >= 0; i1--)
        {
            if (i1 == (local_rows - 1) && size > 1 && rank < size - 1)
            {
                MPI_Send(X[local_rows - 1], n, MPI_DOUBLE, rank + 1, 2, MPI_COMM_WORLD);
                MPI_Send(A[local_rows - 1], n, MPI_DOUBLE, rank + 1, 3, MPI_COMM_WORLD);
            }

            if (i1 == 0)
            {
                if (rank > 0)
                {
                    for (i2 = 0; i2 < n; i2++)
                    {
                        X[i1][i2] = (X[i1][i2] - temp_X[i2] * temp_A[i2]) / B[i1][i2];
                    }
                }
            }
            else if (rank < size - 1 || ((rank == size - 1) && i1 <= local_rows - 2))
            {
                for (i2 = 0; i2 < n; i2++)
                {
                    X[i1][i2] = (X[i1][i2] - X[i1 - 1][i2] * A[i1 - 1][i2]) / B[i1][i2];
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    free(quantum_X);
    free(quantum_B);
    free(temp_X);
    free(temp_A);
}

int main(int argc, char **argv)
{
    int rank, size;      // Ранг и количество процессов
    double time0, time1; // Время выполнения

    if (argc > 1)
    {
        if (atoi(argv[1]) > 0 || atoi(argv[1]) < N)
        {
            num_quanta = atoi(argv[1]);
        }
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    local_rows = (N / size) + (rank < (N % size) ? 1 : 0);
    quantum_size = (N / num_quanta) + ((N % num_quanta > 0) ? 1 : 0);

    init_arrays();

    // Замер времени выполнения
    time0 = MPI_Wtime();
    kernel_adi(TSTEPS, N, X, A, B, rank, size, local_rows);
    time1 = MPI_Wtime();

    if (rank == 0)
    {
        printf("\nN: %d", N);
        printf("\nQuantum size: %d", quantum_size);
        printf("\nNumber of quanta: %d", N / quantum_size);
        printf("\nTotal execution time: %f seconds\n", time1 - time0);
    }

    print_row(31);

    free_arrays();

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
