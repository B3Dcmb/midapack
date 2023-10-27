#ifndef MAPPRAISER_MAPPING_H
#define MAPPRAISER_MAPPING_H

#include <midapack.h>

#ifdef __cplusplus
extern "C" {
#endif

/// @brief Possible choices for taking gaps into account in the map
typedef enum gap_strategy_t {
    // condition on signal being zero inside gaps
    // (original approach)
    COND = 0,

    // marginalize on extra pixels
    // introduce 1 extra pixel per scan
    // the extra pixels are not shared among processes
    MARG_LOCAL_SCAN,

    // inverse the noise covariance matrix anyway
    // despite some rows and columns missing
    NESTED_PCG,

    // inverse the noise covariance matrix iteratively
    // BUT ignore the gaps
    // this will only correct the Toeplitz approximation
    // near the timestream edges
    NESTED_PCG_NO_GAPS,
} GapStrategy;

void print_gap_stgy(GapStrategy gs);

int get_actual_map_size(const Mat *A);

int get_valid_map_size(const Mat *A);

int create_extra_pix(int *indices, double *weights, int nnz, int nb_blocks_loc,
                     const int *local_blocks_sizes, GapStrategy gs);

int build_pixel_to_time_domain_mapping(Mat *A);

void build_gap_struct(int64_t gif, Gap *gaps, Mat *A);

void compute_gaps_per_block(Gap *gaps, int nb_blocks, Block *blocks);

void copy_gap_info(int nb_blocks, Block *src, Block *dest);

int compute_global_gap_count(MPI_Comm comm, Gap *gaps);

void reset_relevant_gaps(double *tod, Tpltz *tmat, Gap *gaps);

void condition_extra_pix_zero(Mat *A);

void point_pixel_to_trash(Mat *A, int ipix);

__attribute__((unused)) void print_gap_info(Gap *gaps);

#ifdef __cplusplus
}
#endif

#endif // MAPPRAISER_MAPPING_H
