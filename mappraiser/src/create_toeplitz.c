/**
 * @file   createToeplitz.c
 * @brief Utility routines for creating the inverse time-time noise
 * correlation matrix
 * @last_update May 2019 by Hamza El Bouhargani
 */

#include <mappraiser/create_toeplitz.h>

int defineTpltz_avg(Tpltz *Nm1, int64_t nrow, int m_cw, int m_rw,
                    Block *tpltzblocks, int nb_blocks_loc, int nb_blocks_tot,
                    int64_t idp, int local_V_size, Flag flag_stgy,
                    MPI_Comm comm) {

    Nm1->nrow = nrow;
    Nm1->m_cw = m_cw;
    Nm1->m_rw = m_rw;
    Nm1->tpltzblocks = tpltzblocks;
    Nm1->nb_blocks_loc = nb_blocks_loc;
    Nm1->nb_blocks_tot = nb_blocks_tot;
    Nm1->idp = idp;
    Nm1->local_V_size = local_V_size;
    Nm1->flag_stgy = flag_stgy;
    Nm1->comm = comm;

    return 0;
}

int defineBlocks_avg(Block *tpltzblocks, double *T, int nb_blocks_loc,
                     void *local_blocks_sizes, int lambda_block_avg,
                     int64_t id0) {

    int i, index0;

    for (i = 0; i < nb_blocks_loc; i++)
        tpltzblocks[i].n = ((int *)local_blocks_sizes)[i];

    for (i = 0; i < nb_blocks_loc; i++)
        tpltzblocks[i].lambda = lambda_block_avg;

    // tpltzblocks[0].idv = (int64_t) (id0/n_block_avg) * n_block_avg ;
    tpltzblocks[0].idv = id0;
    for (i = 1; i < nb_blocks_loc; i++)
        tpltzblocks[i].idv =
            (int64_t)tpltzblocks[i - 1].idv + tpltzblocks[i - 1].n;

    index0 = 0;
    for (i = 0; i < nb_blocks_loc; i++) {
        tpltzblocks[i].T_block = (T + index0);
        index0 += tpltzblocks[i].lambda;
    }

    return 0;
}
