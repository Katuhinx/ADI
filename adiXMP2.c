#include <stdio.h>
#include <stdbool.h>
#include <sys/time.h>
#include "adi.h"
#include "adi_utils/adi_utils.h"

double second();

#pragma xmp template t [N][N]
#pragma xmp nodes p [1][1]
#pragma xmp distribute t[block][block] onto p
DATA_TYPE X[N][N];
#pragma xmp align X [i][j] with t [i][j]  
#pragma xmp shadow X[1][1]
DATA_TYPE B[N][N];
#pragma xmp align B [i][j] with t [i][j] 
#pragma xmp shadow B[1:0][1:0]
DATA_TYPE A[N][N];
#pragma xmp align A [i][j] with t [i][j] 
#pragma xmp shadow A[1:0][1:0]

static void init_array
    ( int n )
{
    int i, j;

    #pragma xmp loop on t [i][j]
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) 
        {
            X[i][j] = ((DATA_TYPE) i*(j+1) + 1) / n;
            A[i][j] = ((DATA_TYPE) i*(j+2) + 2) / n;
            B[i][j] = ((DATA_TYPE) i*(j+3) + 3) / n;
        }
}

static void print_array
    ( int n )
{
    int i, j;
    double A[n][n];

    #pragma xmp loop on t [i][j]
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            A[i][j] = (double)X[i][j];
        }

    if (compare_array_with_file(n, A)) {
        printf("The arrays are identical.\n");
    } else {
        printf("The arrays differ.\n");
    }
}

static void kernel_adi//метод ADI для каждого шага времени
  (
    int tsteps,// Количество шагов времени, которое необходимо выполнить для решения уравнения
    int n// Размер массива, который используется для хранения значений функции
  )
{
    int t, i1, i2;

    for (t = 0; t < tsteps; t++)
    {
       // #pragma xmp reflect(X,B)
        #pragma xmp loop on t [i1][i2] //работает
        for (i1 = 0; i1 < _PB_N; i1++)
            for (i2 = 1; i2 < _PB_N; i2++)//обновление элементов по направлению x
            {
                
                X[i1][i2] = X[i1][i2] - X[i1][i2-1] * A[i1][i2] / B[i1][i2-1];//·	Обновление значений массива X по направлению x 
                B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1][i2-1];//·	Обновление значений массива B по направлению x 
            }
#pragma xmp reflect(X,B)
     //   #pragma xmp loop on t [i1][i2]//не нравиться измерение шаблона
        for (i1 = 0; i1 < _PB_N; i1++)//обновление значений на границе по направлению x
    //          for( i2=_PB_N-1; i2<=_PB_N-1;i2++) 
           X[i1][_PB_N-1] = X[i1][_PB_N-1] / B[i1][_PB_N-1];

//int i3=_PB_N-i2-2;
     //   #pragma xmp reflect(B,X,A)
      //   #pragma xmp loop on t [i1][i3] //не может определить условие цикла
        for (i1 = 0; i1 < _PB_N; i1++)
            for (i2 = 0; i2 < _PB_N-2; i2++)
                X[i1][_PB_N-i2-2] = (X[i1][_PB_N-2-i2] - X[i1][_PB_N-2-i2-1] * A[i1][_PB_N-i2-3]) / B[i1][_PB_N-3-i2];

      //  #pragma xmp reflect(B,X)
        #pragma xmp loop on t [i1][i2]   //работает
        for (i1 = 1; i1 < _PB_N; i1++)
            for (i2 = 0; i2 < _PB_N; i2++)//обновление значений по направлению y
            {
                X[i1][i2] = X[i1][i2] - X[i1-1][i2] * A[i1][i2] / B[i1-1][i2];
                B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1-1][i2];
            }

    //     #pragma xmp loop on t [i1][i2]
    //   for(i1=_PB_N-1;i1<=_PB_N-1;i1++)//не нравиться измерение шаблона
        for (i2 = 0; i2 < _PB_N; i2++)//обновление значений на границе по направлению 
            X[_PB_N-1][i2] = X[_PB_N-1][i2] / B[_PB_N-1][i2];

        // #pragma xmp loop on t [_PB_N-i1-2][i2] reflect(X,A)//не может определить измерение цикла
        for (i1 = 0; i1 < _PB_N-2; i1++)
            for (i2 = 0; i2 < _PB_N; i2++)
                X[_PB_N-2-i1][i2] = (X[_PB_N-2-i1][i2] - X[_PB_N-i1-3][i2] * A[_PB_N-3-i1][i2]) / B[_PB_N-2-i1][i2];  
    }
}

int main(int argc, char** argv)
{
    int n = N;
    double time, time0, time1;

    init_array (n);

    time0 = second();
    kernel_adi(TSTEPS, n);
    time1 = second();

    time = time1 - time0;

    #pragma xmp task on p[0][0]
    printf("\nn=%d\n", n);
    print_array(n);
    printf("\ntime=%f\n", time);

    return 0;
}

double second()
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
        t = (double) (tm.tv_sec-base_sec) + ((double) (tm.tv_usec-base_usec))/1.0e6 ;
    }

    return t ;
}