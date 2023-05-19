#ifndef MAPPRAISER_APPLY_WEIGHTS_H
#define MAPPRAISER_APPLY_WEIGHTS_H

#include <midapack.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    BASIC       = 0, // consider that there are no gaps, do not iterate
    ITER        = 1, // consider that there are gaps, and iterate to mitigate that
    ITER_IGNORE = 2, // iterate, but ignore the gaps (call non-gappy routines)
} weight_stgy_t;

void apply_weights(Tpltz *Nm1, Tpltz *N, Gap *Gaps, double *tod, weight_stgy_t stgy, bool verbose);

void set_tpltz_struct(Tpltz *single_block_struct, Tpltz *full_struct, Block *block);

#ifdef __cplusplus
}
#endif

#endif // MAPPRAISER_APPLY_WEIGHTS_H
