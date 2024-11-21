#ifndef ADI_UTILS_H
#define ADI_UTILS_H

#include <stdio.h>
#include <stdbool.h>
#include "../adi.h"

void save_array_to_file(int n, double X[N][N]);

bool compare_array_with_file(int n, double X[N][N]);

#endif // ADI_UTILS_H