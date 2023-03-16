#ifndef MAPPRAISER_APPLY_WEIGHTS_H
#define MAPPRAISER_APPLY_WEIGHTS_H

#include <midapack.h>

#ifdef __cplusplus
extern "C" {
#endif

void apply_weights(Tpltz *Nm1, Tpltz *N, Gap *Gaps, double *tod);

void set_tpltz_struct(Tpltz *single_block_struct, Tpltz *full_struct, Block *block);

#ifdef __cplusplus
}
#endif

#endif //MAPPRAISER_APPLY_WEIGHTS_H
