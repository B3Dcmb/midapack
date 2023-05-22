/**
 * @file mapping.c
 * @brief Implementation of routines for pixel to time-domain mapping and gap measurement
 * @author Simon Biquard
 * @date Nov 2022
 */

#include "mappraiser/mapping.h"

/**
 * @brief Build the pixel-to-time-domain mapping, i.e.
 * i) A->id_last_pix which contains the indexes of the last samples pointing to each pixel
 * ii) A->ll which is a linked list of time sample indexes
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
        printf("memory allocation of id_last_pix or ll failed\n");
        exit(EXIT_FAILURE);
    }

    // initialize the mapping arrays to -1
    for (i = 0; i < A->m; i++) { A->ll[i] = -1; }
    for (j = 0; j < A->lcount / A->nnz; j++) { A->id_last_pix[j] = -1; }

    // build the linked list chain of time samples corresponding to each pixel
    // and compute number of timestream gaps
    ngap   = 0;
    lengap = 0;
    for (i = 0; i < A->m; i++) {
        ipix = A->indices[i * A->nnz] / A->nnz;
        if (A->id_last_pix[ipix] == -1) {
            A->id_last_pix[ipix] = i;
        } else {
            A->ll[i]             = A->id_last_pix[ipix];
            A->id_last_pix[ipix] = i;
        }

        // compute the number of gaps in the timestream
        if (A->trash_pix) {
            if (A->indices[i * A->nnz] != 0) {
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
 *
 * @param gif global row index offset of the local data
 * @param Gaps Gap structure (Gaps->ngap must already be computed!)
 * @param A pointing matrix structure
 */
void build_gap_struct(int64_t gif, Gap *Gaps, Mat *A) {
    // allocate the arrays

    // only test correct allocation if ngap > 0 because
    // behaviour of malloc(0) is implementation-defined
    // free(NULL) produces no error

    Gaps->id0gap = malloc((sizeof Gaps->id0gap) * Gaps->ngap);
    Gaps->lgap   = malloc((sizeof Gaps->lgap) * Gaps->ngap);

    if (Gaps->ngap > 0) {
        int i         = Gaps->ngap - 1;    // index of the gap being computed
        int lengap    = 1;                 // length of the current gap
        int j         = A->id_last_pix[0]; // index to go through linked time samples
        int gap_start = j;                 // index of the first sample of the gap

        if (Gaps->id0gap == NULL || Gaps->lgap == NULL) {
            printf("malloc of id0gap or lgap failed\n");
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
                Gaps->id0gap[i] = gif + gap_start; // global row index
                Gaps->lgap[i]   = lengap;
                lengap          = 1;
                --i;
            }
            gap_start = j;
        }
    }
}

int compute_global_gap_count(MPI_Comm comm, Gap *gaps) {
    int gap_count = gaps->ngap;
    MPI_Allreduce(MPI_IN_PLACE, &gap_count, 1, MPI_INT, MPI_SUM, comm);
    return gap_count;
}

void print_gap_info(Gap *gaps) {
    printf("Local Gap structure\n");
    printf("  { ngap: %d\n", gaps->ngap);
    printf("    lgap: ");
    for (int i = 0; i < gaps->ngap; ++i) { printf("%d ", gaps->lgap[i]); }
    printf("\n");
    printf("    id0gap: ");
    for (int i = 0; i < gaps->ngap; ++i) { printf("%ld ", gaps->id0gap[i]); }
    printf("}\n");
}