/**
 * @file mapping.c
 * @brief Implementation of routines for pixel to time-domain mapping and gap
 * measurement
 * @author Simon Biquard
 * @date Nov 2022
 */

#include "mappraiser/mapping.h"

#ifndef NDEBUG
#include <assert.h>
#endif

#include <stdbool.h>
#include <stdlib.h>

void print_gap_stgy(GapStrategy gs) {
    switch (gs) {
    case COND:
        puts("conditioning");
        break;
    case MARG_LOCAL_SCAN:
        puts("marginalization (1 extra pixel/scan)");
        break;
    case NESTED_PCG:
        puts("nested PCG");
        break;
    case NESTED_PCG_NO_GAPS:
        puts("nested PCG ignoring gaps");
    }
}

int get_actual_map_size(const Mat *A) {
    if (A->flag_ignore_extra)
        return A->lcount - A->trash_pix * A->nnz;
    else
        return A->lcount;
}

int get_valid_map_size(const Mat *A) {
    return A->lcount - A->trash_pix * A->nnz;
}

int create_extra_pix(int *indices, int nnz, int nb_blocks_loc,
                     const int *local_blocks_sizes, GapStrategy gs) {
    switch (gs) {
    case MARG_LOCAL_SCAN: {
        // position in the data
        int offset = 0;

        // loop through local blocks
        for (int i = 0; i < nb_blocks_loc; ++i) {
            int local_size = local_blocks_sizes[i];
            for (int j = 0; j < local_size * nnz; j += nnz) {
                if (indices[offset + j] < 0) {
                    // set a negative index corresponding to local scan number
                    // don't forget the nnz multiplicity of the indices
                    for (int k = 0; k < nnz; ++k) {
                        indices[offset + j + k] = -(i * nnz + k + 1);
                    }
                }
            }
            offset += local_size;
        }
        break;
    }

    default:
        /* nothing to do */
        break;
    }
    return 0;
}

/**
 * @brief Build the pixel-to-time-domain mapping, i.e.
 * i) A->id_last_pix which contains the indexes of the last samples pointing to
 * each pixel ii) A->ll which is a linked list of time sample indexes
 *
 * @param A the pointing matrix structure
 * @return int the number of timestream gaps found
 */
int build_pixel_to_time_domain_mapping(Mat *A) {
    int i, j;
    int ipix;
    int ngap, lengap;

    // index of last sample pointing to each pixel
    A->id_last_pix = malloc((sizeof A->id_last_pix) * A->lcount / A->nnz);

    // linked list of time samples indexes
    A->ll = malloc((sizeof A->ll) * A->m);

    if (A->id_last_pix == NULL || A->ll == NULL) {
        fputs("memory allocation of id_last_pix or ll failed", stderr);
        exit(EXIT_FAILURE);
    }

    // initialize the mapping arrays to -1
    for (i = 0; i < A->m; i++) {
        A->ll[i] = -1;
    }
    for (j = 0; j < A->lcount / A->nnz; j++) {
        A->id_last_pix[j] = -1;
    }

    // build the linked list chain of time samples corresponding to each pixel
    // and compute number of timestream gaps
    ngap = 0;
    lengap = 0;
    for (i = 0; i < A->m; i++) {
        ipix = A->indices[i * A->nnz] / A->nnz;
        if (A->id_last_pix[ipix] == -1) {
            A->id_last_pix[ipix] = i;
        } else {
            A->ll[i] = A->id_last_pix[ipix];
            A->id_last_pix[ipix] = i;
        }

        // compute the number of gaps in the timestream
        if (A->trash_pix > 0) {
            if (A->indices[i * A->nnz] >= A->trash_pix * A->nnz) {
                // valid sample: reset gap length
                lengap = 0;
            } else {
                // flagged sample -> gap
                if (lengap == 0) {
                    // new gap: increment gap count
                    ++ngap;
                }

                // increment current gap size
                ++lengap;
            }
        }
    }
    return ngap;
}

/**
 * @brief Build the gap structure for the local samples.
 * @param gif global row index offset of the local data
 * @param gaps Gap structure (Gaps->ngap must already be computed!)
 * @param A pointing matrix structure
 */
void build_gap_struct(int64_t gif, Gap *gaps, Mat *A) {
    // allocate the arrays

    // only test correct allocation if ngap > 0 because
    // behaviour of malloc(0) is implementation-defined
    // free(NULL) produces no error

    gaps->id0gap = malloc((sizeof gaps->id0gap) * gaps->ngap);
    gaps->lgap = malloc((sizeof gaps->lgap) * gaps->ngap);

    if (gaps->ngap > 0) {
        int i = gaps->ngap - 1;    // index of the gap being computed
        int lengap = 1;            // length of the current gap
        int j = A->id_last_pix[0]; // index to go through linked time samples
        int gap_start = j;         // index of the first sample of the gap

        if (gaps->id0gap == NULL || gaps->lgap == NULL) {
            fputs("malloc of id0gap or lgap failed", stderr);
            exit(EXIT_FAILURE);
        }

        // go through the time samples
        while (j != -1) {
            // go to previous flagged sample
            j = A->ll[j];

            if (j != -1 && gap_start - j == 1) {
                // same gap, and there are flagged samples left
                ++lengap;
            } else {
                // different gap, or no flagged samples remaining
                gaps->id0gap[i] = gif + gap_start; // global row index
                gaps->lgap[i] = lengap;
                lengap = 1;
                --i;
            }
            gap_start = j;
        }
    }
}

bool gap_overlaps_with_block(Gap *gaps, int i_gap, Block *block) {
    if (i_gap < 0 || i_gap > gaps->ngap - 1)
        return false;
    int64_t id0g = gaps->id0gap[i_gap];
    int64_t idv = block->idv;
    int lg = gaps->lgap[i_gap];
    int n = block->n;
    return (idv < id0g + lg) && (id0g < idv + n);
}

void compute_gaps_per_block(Gap *gaps, int nb_blocks, Block *blocks) {
    int i_gap;       // index to go through the gaps
    int first, last; // indexes of first and last relevant gaps
    Block *b;        // pointer to the current block

    if (gaps->ngap > 0) {
        i_gap = 0;
        for (int i = 0; i < nb_blocks; ++i) {
            b = &(blocks[i]);

            // find the first relevant gap
            while (!gap_overlaps_with_block(gaps, i_gap, b)) {
                ++i_gap;
#ifndef NDEBUG
                assert(i_gap < gaps->ngap);
#endif
            }

            // store its index
            first = i_gap;
            // printf("first = %d\n", first);

            if (first == -1) {
                // no relevant gaps found for this block
                last = -1;
            } else {
                // go through relevant gaps
                while (gap_overlaps_with_block(gaps, i_gap, b)) {
                    ++i_gap;
                }

                // store the index of the last relevant gap
                last = i_gap - 1;
            }
            // printf("last = %d\n", last);

            // store the information for this block
            b->first_gap = first;
            b->last_gap = last;
        }
    } else {
        // no local gaps: set everything to -1
        for (int i = 0; i < nb_blocks; ++i) {
            blocks[i].first_gap = -1;
            blocks[i].last_gap = -1;
        }
    }
}

void copy_gap_info(int nb_blocks, Block *src, Block *dest) {
    for (int i = 0; i < nb_blocks; ++i) {
        dest[i].first_gap = src[i].first_gap;
        dest[i].last_gap = src[i].last_gap;
    }
}

int compute_global_gap_count(MPI_Comm comm, Gap *gaps) {
    int gap_count = gaps->ngap;
    MPI_Allreduce(MPI_IN_PLACE, &gap_count, 1, MPI_INT, MPI_SUM, comm);
    return gap_count;
}

void fill_gap_with_zero(double *tod, int n, int64_t idv, int64_t id0g, int lg) {
#ifndef NDEBUG
    // assert that gap is relevant for the given data block
    assert(idv < id0g + lg && id0g < idv + n);
#endif
    // set intersection of tod and gap to zero
    for (int64_t j = id0g; j < id0g + lg; ++j) {
        if (idv <= j && j < n + idv)
            tod[j - idv] = 0;
    }
}

/**
 * Fill all timestream gaps of a vector with zeros. Warning! This routine
 * assumes that relevant gaps have been determined for each data block (e.g.
 * through a call to compute_gaps_per_block).
 * @param tod pointer to the data vector
 * @param tmat pointer to a Tpltz matrix containing information about the data
 * blocks
 * @param gaps pointer to the gaps structure
 */
void reset_relevant_gaps(double *tod, Tpltz *tmat, Gap *gaps) {
    // loop over data blocks
    Block *b;
    double *tod_block;
    int pos = 0;
    for (int i = 0; i < tmat->nb_blocks_loc; ++i) {
        b = &(tmat->tpltzblocks[i]);
        tod_block = (tod + pos);
        // loop over the relevant gaps for this block
        for (int j = b->first_gap; j <= b->last_gap; ++j) {
            fill_gap_with_zero(tod_block, b->n, b->idv, gaps->id0gap[j],
                               gaps->lgap[j]);
        }
        pos += b->n;
    }
}

void condition_extra_pix_zero(Mat *A) {
    // number of extra pixels in the pointing matrix
    int extra = A->trash_pix * A->nnz;
    int nnz = A->nnz;

    // if no extra pixels, there is nothing to do
    if (extra == 0)
        return;

    // set pointing weights for extra pixels to zero
    for (int i = 0; i < A->m * nnz; i += nnz) {
        if (A->indices[i] < extra) {
            for (int j = 0; j < nnz; ++j) {
                A->values[i + j] = 0;
            }
        }
    }
}

void point_pixel_to_trash(Mat *A, int ipix) {
    // last index of time sample pointing to pixel
    int j = A->id_last_pix[ipix];
    int nnz = A->nnz;

    while (j != -1) {
        // point sample to trash pixel
        for (int k = 0; k < nnz; k++) {
            A->indices[j * nnz + k] = k - nnz;
        }
        j = A->ll[j];
    }
}

__attribute__((unused)) void print_gap_info(Gap *gaps) {
    printf("Local Gap structure\n");
    printf("  { ngap: %d\n", gaps->ngap);
    printf("    lgap: ");
    for (int i = 0; i < gaps->ngap; ++i) {
        printf("%d ", gaps->lgap[i]);
    }
    printf("\n");
    printf("    id0gap: ");
    for (int i = 0; i < gaps->ngap; ++i) {
        printf("%ld ", gaps->id0gap[i]);
    }
    printf("}\n");
}
