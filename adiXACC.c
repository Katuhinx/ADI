#include <stdio.h>
#include "adi.h"
#include <sys/time.h>


double 
second();

 static DATA_TYPE X[N][N];
 static DATA_TYPE A[N][N];
 static DATA_TYPE B[N][N];

#pragma xmp template t[N][N]//шаблон для матрицы
#pragma xmp nodes n [1][1]

#pragma xmp distribute t[block][block] onto n
#pragma xmp align X [i][j] with t[i][j]
#pragma xmp align A [i][j] with t[i][j]
#pragma xmp align B [i][j] with t[i][j]
#pragma xmp shadow X[1:0][1:0]
#pragma xmp shadow B[1:0][1:0]

static
void init_array (int n,//инициализация массивов
   DATA_TYPE X[N][N],
   DATA_TYPE A[N][N],
   DATA_TYPE B[N][N])
{
//#pragma xmp gmove(X, A, B) in//массивы X, A и B должны быть перемещены из глобальной памяти в локальную память
  int i, j;
  #pragma xmp loop on t[i][j] //параллельное выполнение оператора цикла
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
  //  #pragma xmp gmove(X, A, B) out//массивы X, A и B должны быть перемещены из локальной памяти в глобальную память
  int t, i1, i2;

  for (t = 0; t < _PB_TSTEPS; t++)
    {
     #pragma xmp reflect (B)   
     #pragma xmp loop  on t[i1][i2]
      for (i1 = 0; i1 < _PB_N; i1++)
        for (i2 = 1; i2 < _PB_N; i2++)//обновление элементов по направлению x
             {
           //  #pragma xmp depend (in: X[i1][i2-1], out: X[i1][i2]) 
              X[i1][i2] = X[i1][i2] - X[i1][i2-1] * A[i1][i2] / B[i1][i2-1];//·	Обновление значений массива X по направлению x 
              
         //    #pragma xmp depend (in: B[i1][i2-1], out: B[i1][i2]) 
              B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1][i2-1];//·	Обновление значений массива B по направлению x 
             }

    #pragma xmp loop [i1] on t[i1][i2]
      for (i1 = 0; i1 < _PB_N; i1++)//обновление значений на границе по направлению x
         X[i1][_PB_N-1] = X[i1][_PB_N-1] / B[i1][_PB_N-1];

    #pragma xmp loop on t[i1][i2]
      for (i1 = 0; i1 < _PB_N; i1++)
        for (i2 = 0; i2 < _PB_N-2; i2++)
             X[i1][_PB_N-i2-2] = (X[i1][_PB_N-2-i2] - X[i1][_PB_N-2-i2-1] * A[i1][_PB_N-i2-3]) / B[i1][_PB_N-3-i2];
    
    #pragma xmp loop on t[i1][i2]
    for (i1 = 1; i1 < _PB_N; i1++)
        for (i2 = 0; i2 < _PB_N; i2++)//обновление значений по направлению y
            #pragma xmp gmove 
            {
                X[i1][i2] = X[i1][i2] - X[i1-1][i2] * A[i1][i2] / B[i1-1][i2];
                B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1-1][i2];
            }

    #pragma xmp loop [i2] on t[i1][i2]
      for (i2 = 0; i2 < _PB_N; i2++)//обновление значений на границе по направлению 
            X[_PB_N-1][i2] = X[_PB_N-1][i2] / B[_PB_N-1][i2];

    #pragma xmp loop on t[i1][i2]
      for (i1 = 0; i1 < _PB_N-2; i1++)
        for (i2 = 0; i2 < _PB_N; i2++)
         X[_PB_N-2-i1][i2] = (X[_PB_N-2-i1][i2] - X[_PB_N-i1-3][i2] * A[_PB_N-3-i1][i2]) / B[_PB_N-2-i1][i2];
 //    #pragma xmp reflect(X,B)
     }

}


int main(int argc, char** argv)
{

  int n = N;
  int tsteps = TSTEPS;
  double time, time0, time1;

  init_array (n, X, A, B);

  time0 = second();
  kernel_adi (tsteps, n, X, A, B);
  time1 = second();

  time = time1 - time0;
  #pragma xmp reduction(max:time)//вычисляет максимальное значение переменной time, хранящееся в памяти каждого ускорителя в каждом узле.
  #pragma xmp task on t[0][0]//узел p[0][0] выполняет print и выводит указанный текст на экран
 {
  printf("\ntime=%f\n\n\n", time);
  print_array(n, X);
  printf("\n\n\ntime=%f\n", time);
 }
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


