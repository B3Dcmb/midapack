#ifndef MAPPRAISER_MAPPING_H
#define MAPPRAISER_MAPPING_H

#include <midapack.h>

int build_pixel_to_time_domain_mapping(Mat *A);

void build_gap_struct(int64_t gif, Gap *Gaps, Mat *A);

#endif //MAPPRAISER_MAPPING_H
