#include "adi_utils.h"
#include <math.h>

void save_result(int n, double X[n][n]) {
    FILE *file = fopen("output.txt", "wb");

    if (file == NULL) {
        perror("Error opening file for writing");
        return;
    }

    size_t written = fwrite(X, sizeof(double), n * n, file);
    if (written != n * n) {
        fprintf(stderr, "Error writing to file: only %zu out of %d elements written.\n", written, n * n);
    }

    fclose(file);
}

bool compare_result(int n, double X[n][n]) {
    const double EPSILON = 1e-5;
    double value_from_result, value_from_output;

    FILE *result_file = fopen("result.txt", "wb");
    if (result_file == NULL) {
        fprintf(stderr, "Error opening result.txt for writing.\n");
        return false;
    }

    size_t written = fwrite(X, sizeof(double), n * n, result_file);
    if (written != n * n) {
        fprintf(stderr, "Error writing to result.txt: only %zu out of %d elements written.\n", written, n * n);
        fclose(result_file);
        return false;
    }
    fclose(result_file);

    FILE *output_file = fopen("output.txt", "rb");
    if (output_file == NULL) {
        fprintf(stderr, "Error opening output.txt for reading.\n");
        return false;
    }

    result_file = fopen("result.txt", "rb");
    if (result_file == NULL) {
        fprintf(stderr, "Error opening result.txt for reading.\n");
        fclose(output_file);
        return false;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fread(&value_from_result, sizeof(double), 1, result_file) != 1) {
                fprintf(stderr, "Error reading value from result.txt.\n");
                fclose(result_file);
                fclose(output_file);
                return false;
            }
            if (fread(&value_from_output, sizeof(double), 1, output_file) != 1) {
                fprintf(stderr, "Error reading value from output.txt.\n");
                fclose(result_file);
                fclose(output_file);
                return false;
            }
            if (fabs(value_from_result - value_from_output) > EPSILON) {
                fprintf(stderr, "Arrays differ at [%d][%d]: %f (in result.txt) != %f (in output.txt)\n", i, j, value_from_result, value_from_output);
                fclose(result_file);
                fclose(output_file);
                return false;
            }
        }
    }

    fclose(result_file);
    fclose(output_file);
    return true;
}