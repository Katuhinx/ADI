#include "adi_utils.h"
#include <math.h>

void save_array_to_file(int n, double X[n][n]) {
    FILE *file = fopen("result.txt", "w");

    fprintf(file, "\n");

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

bool compare_array_with_file(int n, double X[n][n]) {
    FILE *file = fopen("result.txt", "r");

    fprintf(file, "\n");

    if (file == NULL) {
        fprintf(stderr, "Error opening file for reading.\n");
        return false;
    }

    const float EPSILON = 1e-5;
    double value;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fscanf(file, "%lf", &value) != 1) {
                fprintf(stderr, "Error reading value from file.\n");
                fclose(file);
                return false;
            }
            if (fabs((double)X[i][j] - (double)value) > EPSILON) {
                fprintf(stderr, "Arrays differ at [%d][%d]: %f (in array) != %f (in file)\n", i, j, X[i][j], value);
                fclose(file);
                return false;
            }
        }
    }

    fclose(file);
    return true;
}