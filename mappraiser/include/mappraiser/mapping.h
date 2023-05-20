#ifndef MAPPRAISER_MAPPING_H
#define MAPPRAISER_MAPPING_H

#include <midapack.h>

#ifdef __cplusplus
extern "C" {
#endif

int build_pixel_to_time_domain_mapping(Mat *A);

void build_gap_struct(int64_t gif, Gap *Gaps, Mat *A);

int compute_global_gap_count(MPI_Comm comm, Gap *gaps);

void print_gap_info(Gap *gaps);

#ifdef __cplusplus
}
#endif

#endif // MAPPRAISER_MAPPING_H
