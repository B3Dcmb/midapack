#ifndef MAPPRAISER_CREATE_TOEPLITZ_H
#define MAPPRAISER_CREATE_TOEPLITZ_H

#include "toeplitz/toeplitz.h"

int defineTpltz_avg(Tpltz *Nm1, int64_t nrow, int m_cw, int m_rw, Block *tpltzblocks, int nb_blocks_loc,
                    int nb_blocks_tot, int64_t idp, int local_V_size, Flag flag_stgy, MPI_Comm comm);

int defineBlocks_avg(Block *tpltzblocks, double *T, int nb_blocks_loc, void *local_blocks_sizes, int lambda_block_avg,
                     int64_t id0);

#endif // MAPPRAISER_CREATE_TOEPLITZ_H
