#ifndef MAPPRAISER_WEIGHT_H
#define MAPPRAISER_WEIGHT_H

#include <midapack.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum weight_stgy_t {
    BASIC = 0, // consider that there are no gaps, do not iterate
    ITER = 1,  // consider that there are gaps, and iterate to mitigate that
    ITER_IGNORE = 2, // iterate, but ignore the gaps (call non-gappy routines)
} WeightStgy;

typedef struct {
    Tpltz *Nm1; // Approximate Toeplitz inverse noise covariance
    Tpltz *N;   // Toeplitz noise covariance
    Gap *G;     // Timestream gaps
    WeightStgy stgy;
} WeightMatrix;

// Simple constructor
WeightMatrix createWeightMatrix(Tpltz *Nm1, Tpltz *N, Gap *G, WeightStgy stgy);

int applyWeightMatrix(const WeightMatrix *W, double *tod);

void set_tpltz_struct(Tpltz *single_block_struct, const Tpltz *full_struct,
                      Block *block);

#ifdef __cplusplus
}
#endif

#endif // MAPPRAISER_WEIGHT_H
