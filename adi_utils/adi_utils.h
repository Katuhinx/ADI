#ifndef ADI_UTILS_H
#define ADI_UTILS_H

#include <stdio.h>
#include <stdbool.h>
#include "../adi.h"

void save_result(int n, double X[N][N]);

bool compare_result(int n, double X[N][N]);

#endif // ADI_UTILS_H