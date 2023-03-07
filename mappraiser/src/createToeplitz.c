// MAPPRAISER utils vdev
// Utilitary routines for creating the inverse time-time noise correlation matrix.

/** @file   createToeplitz.c
    @last_update May 2019 by Hamza El Bouhargani */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
// #include "s2hat.h"
// #include "s2hat_tools.h"
// #include "domain_generalization.h"
#include "mappraiser.h"

char CHAR_RW = '\0'; // global variable for write mode

int defineTpltz_avg(Tpltz *Nm1, int64_t nrow, int m_cw, int m_rw, Block *tpltzblocks, int nb_blocks_loc, int nb_blocks_tot, int64_t idp, int local_V_size, Flag flag_stgy, MPI_Comm comm)
{

    // faire les allocs ici avec la structure Tpltz

    Nm1->nrow = nrow;                   // glob //recup du fichier params apres (en variables globales)
    Nm1->m_cw = m_cw;                   // glob
    Nm1->m_rw = m_rw;                   // glob
    Nm1->tpltzblocks = tpltzblocks;     // toep
    Nm1->nb_blocks_loc = nb_blocks_loc; // toep
    Nm1->nb_blocks_tot = nb_blocks_tot; // toep
    Nm1->idp = idp;                     // comput
    Nm1->local_V_size = local_V_size;   // comput
    Nm1->flag_stgy = flag_stgy;         // param
    Nm1->comm = comm;                   // param

    return 0;
}

int defineBlocks_avg(Block *tpltzblocks, double *T, int nb_blocks_loc, void *local_blocks_sizes, int lambda_block_avg, int64_t id0)
{

    int i, index0;

    for (i = 0; i < nb_blocks_loc; i++)
        tpltzblocks[i].n = ((int *)local_blocks_sizes)[i];

    for (i = 0; i < nb_blocks_loc; i++)
        tpltzblocks[i].lambda = lambda_block_avg;

    // tpltzblocks[0].idv = (int64_t) (id0/n_block_avg) * n_block_avg ;
    tpltzblocks[0].idv = id0;
    for (i = 1; i < nb_blocks_loc; i++)
        tpltzblocks[i].idv = (int64_t)tpltzblocks[i - 1].idv + tpltzblocks[i - 1].n;

    index0 = 0;
    for (i = 0; i < nb_blocks_loc; i++)
    {
        tpltzblocks[i].T_block = (T + index0);
        index0 += tpltzblocks[i].lambda;
    }

    return 0;
}
