#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "adi.h"
#include <sys/time.h>


double 
second();

DATA_TYPE X[N][N];
  DATA_TYPE A[N][N];
  DATA_TYPE B[N][N];

static
void init_array (int n,//инициализация массивов
   DATA_TYPE X[N][N],
   DATA_TYPE A[N][N],
   DATA_TYPE B[N][N])
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) 
    {
      X[i][j] = ((DATA_TYPE) i*(j+1) + 1) / n;
      A[i][j] = ((DATA_TYPE) i*(j+2) + 2) / n;
      B[i][j] = ((DATA_TYPE) i*(j+3) + 3) / n;
    }
}

static
void print_array(int n,
   DATA_TYPE X[N][N])
{
  int i, j;

  // for (i = 0; i < n; i++)
  //   for (j = 0; j < n; j++) {
  //     fprintf(stderr, DATA_PRINTF_MODIFIER, X[i][j]);
  //     if ((i * N + j) % 20 == 0) fprintf(stderr, "\n");
  //   }
  // fprintf(stderr, "\n");
  
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      printf("%f ", X[i][j]);
    }
}

static
void kernel_adi//метод ADI для каждого шага времени
(int tsteps,// Количество шагов времени, которое необходимо выполнить для решения уравнения
  int n,// Размер массива, который используется для хранения значений функции
  DATA_TYPE X[N][N],//Массив, который хранит значения функции в текущий момент времени
  DATA_TYPE A[N][N],// Массив, который хранит коэффициенты уравнения по направлению x
  DATA_TYPE B[N][N])//Массив, который хранит коэффициенты уравнения по направлению y
{
  int t, i1, i2;

  for (t = 0; t < _PB_TSTEPS; t++)
    {
      for (i1 = 0; i1 < _PB_N; i1++)
 for (i2 = 1; i2 < _PB_N; i2++)//обновление элементов по направлению x
   {
     X[i1][i2] = X[i1][i2] - X[i1][i2-1] * A[i1][i2] / B[i1][i2-1];//·	Обновление значений массива X по направлению x 
     B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1][i2-1];//·	Обновление значений массива B по направлению x 
   }

      for (i1 = 0; i1 < _PB_N; i1++)//обновление значений на границе по направлению x
 X[i1][_PB_N-1] = X[i1][_PB_N-1] / B[i1][_PB_N-1];

      for (i1 = 0; i1 < _PB_N; i1++)
        for (i2 = 0; i2 < _PB_N-2; i2++)
         X[i1][_PB_N-i2-2] = (X[i1][_PB_N-2-i2] - X[i1][_PB_N-2-i2-1] * A[i1][_PB_N-i2-3]) / B[i1][_PB_N-3-i2];
/*
      for (i1 = 1; i1 < _PB_N; i1++)
 for (i2 = 0; i2 < _PB_N; i2++)//обновление значений по направлению y
  {
   X[i1][i2] = X[i1][i2] - X[i1-1][i2] * A[i1][i2] / B[i1-1][i2];
   B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1-1][i2];
 }

      for (i2 = 0; i2 < _PB_N; i2++)//обновление значений на границе по направлению 
 X[_PB_N-1][i2] = X[_PB_N-1][i2] / B[_PB_N-1][i2];

      for (i1 = 0; i1 < _PB_N-2; i1++)
 for (i2 = 0; i2 < _PB_N; i2++)
   X[_PB_N-2-i1][i2] = (X[_PB_N-2-i1][i2] - X[_PB_N-i1-3][i2] * A[_PB_N-3-i1][i2]) / B[_PB_N-2-i1][i2];*/
    }
}


/*
static void save_array_to_file(int n, DATA_TYPE X[N][N]) {
    FILE *file = fopen("result.txt", "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file for writing.\n");
        return;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(file, "%lf ", X[i][j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}



bool compare_array_with_file(int n, DATA_TYPE X[N][N]) {
    FILE *file = fopen("result.txt", "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening file for reading.\n");
        return false;
    }

    DATA_TYPE value;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fscanf(file, "%lf", &value) != 1) {
                fprintf(stderr, "Error reading value from file.\n");
                fclose(file);
                return false;
            }
            if (X[i][j] != value) {
                fprintf(stderr, "Arrays differ at [%d][%d]: %f (in array) != %f (in file)\n", i, j, X[i][j], value);
                fclose(file);
                return false;
            }
        }
    }

    fclose(file);
    return true;
}

*/

int main(int argc, char** argv)
{

  int n = N;
  int tsteps = TSTEPS;

  DATA_TYPE X[N][N];
  DATA_TYPE A[N][N];
  DATA_TYPE B[N][N];
  
  double time, time0, time1;

  init_array (n, X, A, B);

  time0 = second();
  kernel_adi (tsteps, n, X, A, B);
  time1 = second();



  time = time1 - time0;
  printf("\ntime=%f\n\n\n", time);

 // save_array_to_file(n, X);

 /* if (compare_array_with_file(n, X)) {
    printf("The arrays are identical.\n");
  } else {
    printf("The arrays differ.\n");
  }
 */
  print_array(n, X);
  
  printf("\n\n\ntime=%f\n", time);

  return 0;
}

double
second()
{

  struct timeval tm;
  double t;

  static int base_sec = 0,base_usec = 0;

  gettimeofday(&tm, NULL);
  
  if(base_sec == 0 && base_usec == 0)
    {
      base_sec = tm.tv_sec;
      base_usec = tm.tv_usec;
      t = 0.0;
  } else {
    t = (double) (tm.tv_sec-base_sec) + 
      ((double) (tm.tv_usec-base_usec))/1.0e6 ;
  }

  return t ;
}