/*
@file toeplitz_gappy.c version 1.1b, July 2012
@brief Gappy routines used to compute the Toeplitz product when gaps are defined
@author  Frederic Dauvergne
**
** Project:  Midapack library, ANR MIDAS'09 - Toeplitz Algebra module
** Purpose:  Provide Toeplitz algebra tools suitable for Cosmic Microwave
Background (CMB)
**           data analysis.
**
***************************************************************************
@note Copyright (c) 2010-2012 APC CNRS Université Paris Diderot
@note
@note This program is free software; you can redistribute it and/or modify it
under the terms
@note of the GNU Lesser General Public License as published by the Free Software
Foundation;
@note either version 3 of the License, or (at your option) any later version.
This program is
@note distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even
@note the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU
@note Lesser General Public License for more details.
@note
@note You should have received a copy of the GNU Lesser General Public License
along with this
@note program; if not, see http://www.gnu.org/licenses/lgpl.html
@note
@note For more information about ANR MIDAS'09 project see :
@note http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/index.html
@note
@note ACKNOWLEDGMENT: This work has been supported in part by the French
National Research
@note Agency (ANR) through COSINUS program (project MIDAS no. ANR-09-COSI-009).
***************************************************************************
** Log: toeplitz*.c
**
** Revision 1.0b  2012/05/07  Frederic Dauvergne (APC)
** Official release 1.0beta. The first installement of the library is the
Toeplitz algebra
** module.
**
** Revision 1.1b  2012/07/-  Frederic Dauvergne (APC)
** - mpi_stbmm allows now rowi-wise order per process datas and no-blocking
communications.
** - OMP improvment for optimal cpu time.
** - bug fixed for OMP in the stmm_basic routine.
** - distcorrmin is used to communicate only lambda-1 datas when it is needed.
** - new reshaping routines using transformation functions in stmm. Thus, only
one copy
**   at most is needed.
** - tpltz_init improvement using define_nfft and define_blocksize routines.
** - add Block struture to define each Toeplitz block.
** - add Flag structure and preprocessing parameters to define the computational
strategy.
**
**
***************************************************************************
**
*/

#include "toeplitz.h"

#define max(a, b)                                                              \
    ({                                                                         \
        __typeof__(a) _a = (a);                                                \
        __typeof__(b) _b = (b);                                                \
        _a > _b ? _a : _b;                                                     \
    })

#define min(a, b)                                                              \
    ({                                                                         \
        __typeof__(a) _a = (a);                                                \
        __typeof__(b) _b = (b);                                                \
        _a < _b ? _a : _b;                                                     \
    })

// r1.1 - Frederic Dauvergne (APC)
// this is the gappy routines used when gaps are defined

//====================================================================
#ifdef W_MPI
/// Performs the multiplication of a symmetric, Toeplitz block-diagonal matrix,
/// T, by an arbitrary matrix, V, distributed over processes in the generalized
/// column-wise way. This matrix V contains defined gaps which represents the
/// useless data for the comutation. The gaps indexes are defined in the global
/// time space as the generized toeplitz matrix, meaning the row dimension. Each
/// of his diagonal blocks is a symmetric, band-diagonal Toeplitz matrix, which
/// can be different for each block.
/** @ingroup group12
    We first rebuild the Toeplitz block matrix structure to reduce the
   computation cost and skip the computations of the values on the defined gaps.
   then, each process performs the multiplication sequentially for each of the
   gappy block and based on the sliding window algorithm. Prior to that MPI
   calls are used to exchange data between neighboring process. The parameters
   are : \param V \b [input] distributed data matrix (with the convention
   V(i,j)=V[i+j*n]) ; \b [out] result of the product TV \param nrow number of
   rows of the global data matrix V \param m number of columns for the data
   matrix V in the global rowwise order \param m_rowwise number of columns for
   the data matrix V in the rowwise order per processor \param tpltzblocks list
   of the toeplitz blocks struture with its own parameters (idv, n, T_block,
   lambda) :
    - idv is the global row index defining for each Toeplitz block as stored in
   the vector T ;
    - n size of each Toeplitz block
    - T_block location of each Toeplitz matrix data composed of the non-zero
   entries of the first row ;
    - lambda size of each Toeplitz block data T_block. The bandwith size is then
   equal to lambda*2-1 \param nb_blocks_all number of all Toeplitz block on the
   diagonal of the full Toeplitz matrix \param nb_blocks_local number of
   Toeplitz blocks as stored in T \param idp global index of the first element
   of the local part of V \param local_V_size a number of all elements in local
   V \param id0gap index of the first element of each defined gap \param lgap
   length of each defined gaps \param ngap number of defined gaps \param
   flag_stgy flag strategy for the product computation \param comm MPI
   communicator
*/
int mpi_gstbmm(double **V, int nrow, int m, int m_rowwise, Block *tpltzblocks,
               int nb_blocks_local, int nb_blocks_all, int id0p,
               int local_V_size, int64_t *id0gap, int *lgap, int ngap,
               Flag flag_stgy, MPI_Comm comm) {

    // MPI parameters
    int rank; // process rank
    int size; // process number

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int i, j, k; // some indexes

    int flag_skip_build_gappy_blocks = flag_stgy.flag_skip_build_gappy_blocks;

    FILE *file;
    file       = stdout;
    PRINT_RANK = rank;

    // put zeros at the gaps location
    reset_gaps(V, id0p, local_V_size, m, nrow, m_rowwise, id0gap, lgap, ngap);

    // allocation for the gappy structure of the diagonal block Toeplitz matrix
    int nb_blocks_gappy;

    int nb_blockgappy_max;
    int Tgappysize_max;

    Block *tpltzblocks_gappy;

    // some computation usefull to determine the max size possible for the gappy
    // variables
    int Tsize     = 0;
    int lambdamax = 0;

    if (VERBOSE)
        fprintf(file, "[%d] flag_skip_build_gappy_blocks=%d\n", rank,
                flag_skip_build_gappy_blocks);

    if (flag_skip_build_gappy_blocks == 1) { // no build gappy blocks strategy,
                                             // just put zeros at gaps location

        // compute the product using only the input Toeplitz blocks structure
        // with zeros at the gaps location
        mpi_stbmm(V, nrow, m, m_rowwise, tpltzblocks, nb_blocks_local,
                  nb_blocks_all, id0p, local_V_size, flag_stgy, MPI_COMM_WORLD);
    } else { // build gappy blocks strategy

        for (Tsize = i = 0; i < nb_blocks_local; i++)
            Tsize += tpltzblocks[i].lambda;

        for (i = 0; i < nb_blocks_local; i++) {
            if (tpltzblocks[i].lambda > lambdamax)
                lambdamax = tpltzblocks[i].lambda;
        }

        // compute max size possible for the gappy variables
        nb_blockgappy_max = nb_blocks_local + ngap;
        Tgappysize_max    = Tsize + lambdamax * ngap;

        // allocation of the gappy variables with max size possible
        tpltzblocks_gappy = (Block *) calloc(nb_blockgappy_max, sizeof(Block));

        // build gappy Toeplitz block structure considering significant gaps
        // locations, meaning we skip the gaps lower than the minimum
        // correlation distance. You can also use the flag_param_distmin_fixed
        // parameter which allows you to skip the gap lower than these value.
        // Indeed, sometimes it's better to just put somes zeros than to
        // consider two separates blocks. ps: This criteria could be dependant
        // of the local lambda in futur impovements.
        int flag_param_distmin_fixed = flag_stgy.flag_param_distmin_fixed;
        build_gappy_blocks(nrow, m, tpltzblocks, nb_blocks_local, nb_blocks_all,
                           id0gap, lgap, ngap, tpltzblocks_gappy,
                           &nb_blocks_gappy, flag_param_distmin_fixed);

        if (VERBOSE) {
            fprintf(file, "[%d] nb_blocks_gappy=%d\n", rank, nb_blocks_gappy);
            for (i = 0; i < nb_blocks_gappy; i++)
                fprintf(file, "[%d] idvgappy[%d]=%ld ; ngappy[%d]=%d\n", rank,
                        i, tpltzblocks_gappy[i].idv, i, tpltzblocks_gappy[i].n);
        }
        // ps: we could reallocate the gappy variables to their real size. Not
        // sure it's worth it.

        // compute the product using the freshly created gappy Toeplitz blocks
        // structure
        mpi_stbmm(V, nrow, m, m_rowwise, tpltzblocks_gappy, nb_blocks_local,
                  nb_blocks_all, id0p, local_V_size, flag_stgy, MPI_COMM_WORLD);

    } // end flag_skip_build_gappy_blocks==1

    // put zeros on V for the gaps location again. Indeed, some gaps are just
    // replaced by zeros in input, it's created some fakes results we need to
    // clear after the computation.
    reset_gaps(V, id0p, local_V_size, m, nrow, m_rowwise, id0gap, lgap, ngap);

    return 0;
}

//====================================================================
/// Set the data to zeros at the gaps location.
/** @ingroup group11
    The datas located on a gap are set to zeros. The gaps are defined in the
   time space, meaning their indexes are defined in the row dimension.
*/

// put zeros on V for the gaps location
int reset_gaps(double **V, int64_t id0, int local_V_size, int m, int64_t nrow, int m_rowwise, const int64_t *id0gap,
               const int *lgap, int ngap) {
    int i, j, k, l;

    for (j = 0; j < m; j++) {

#pragma omp parallel for private(i) schedule(dynamic, 1)
        for (k = 0; k < ngap; k++)
            for (i = 0; i < lgap[k]; i++)
                if (id0gap[k] + i + j * nrow >= id0 && id0gap[k] + i + j * nrow < id0 + local_V_size) {
                    for (l = 0; l < m_rowwise; l++)
                        (*V)[id0gap[k] + i + j * nrow - id0
                             + l * local_V_size] = 0.;
                }
    }

    return 0;
}

#endif

//====================================================================
/// Build the gappy Toeplitz block structure to optimise the product computation
/// at gaps location.
/** @ingroup group21
    Considering the significant gaps, the blocks to which they belong are cut
   and split between the gap's edges to reduce the total row size of the
   flotting blocks. It take into consideration the minimum correlation length
   and a parameter that allows us to control the minimum gap size allowed to
   split the blocks. In some cases, the gap can be partially reduce to fit the
   minimum block size needed for computation or just for performance criteria.
   This is based on the fact that the gaps are previously set to zeros before
   calling this routine. \param nrow number of rows of the global data matrix V
    \param m number of columns for the data matrix V in the global rowwise order
    \param tpltzblocks list of the toeplitz blocks struture with its own
   parameters (idv, n, T_block, lambda). \param nb_blocks_local number of
   Toeplitz blocks as stored in T \param nb_blocks_all number of all Toeplitz
   block on the diagonal of the full Toeplitz matrix \param id0gap index of the
   first element of each defined gap \param lgap length of each defined gaps
    \param ngap number of defined gaps
    \param tpltzblocks_gappy list of the gappy toeplitz blocks struture with its
   own parameters \param nb_blocks_gappy_final real number of obtained gappy
   Toeplitz blocks \param flag_param_distmin_fixed flag to defined the minimum
   gap value allowed to split a Toeplitz block
*/
int build_gappy_blocks(int nrow, int m, Block *tpltzblocks, int nb_blocks_local,
                       int nb_blocks_all, const int64_t *id0gap,
                       const int *lgap, int ngap, Block *tpltzblocks_gappy,
                       int *nb_blocks_gappy_final,
                       int  flag_param_distmin_fixed) {

    int i, j, k;
    int id, ib;
    int idtmp;
    int igapfirstblock, igaplastblock;

    int param_distmin = 0;
    if (flag_param_distmin_fixed != 0) param_distmin = flag_param_distmin_fixed;

    int lambdaShft;

    int igaplastblock_prev = -1;
    int lambdaShftgappy    = 0;
    int offset_id          = 0;
    // int offset_id_gappy=0;

    int flag_igapfirstinside, flag_igaplastinside;
    int nbloc         = nb_blocks_local;
    int nblocks_gappy = 0;

    int idvtmp_firstblock;

    int nb_blockgappy_max = nb_blocks_local + ngap;
    int Tgappysize_max;

    int ngappy_tmp;
    int lgap_tmp;

    int flag_gapok = 0;

    int distcorr_min;

    int Tgappysize = 0;
    int k_prev     = -1;

    for (k = 0; k < ngap; k++) {

        // find block for the gap begining
        for (igapfirstblock = -1; igapfirstblock == -1;) {
            idtmp = id0gap[k];

            for (ib = 0; ib < nbloc; ib++) {
                if (tpltzblocks[ib].n != 0
                    && idtmp % nrow < tpltzblocks[ib].idv + tpltzblocks[ib].n)
                    break;
            }

            if (ib < nbloc && tpltzblocks[ib].idv <= idtmp) {
                igapfirstblock       = ib; // the block contained the id0gap
                flag_igapfirstinside = 1;
            } else if (ib < nbloc && tpltzblocks[ib].idv > idtmp) {
                igapfirstblock       = ib; // first block after the id0gap
                flag_igapfirstinside = 0;
            } else {                       // ib=nbloc
                igapfirstblock       = -2; // no block defined
                flag_igapfirstinside = 0;
            }
        }

        // find block for the end of the gap - reverse way
        for (igaplastblock = -1; igaplastblock == -1;) {
            idtmp = id0gap[k] + lgap[k] - 1;

            for (ib = nbloc - 1; ib >= 0; ib--) {
                if (tpltzblocks[ib].n != 0 && tpltzblocks[ib].idv <= idtmp)
                    break;
            }

            if (ib >= 0 && idtmp < tpltzblocks[ib].idv + tpltzblocks[ib].n) {
                igaplastblock       = ib;
                flag_igaplastinside = 1;
            } else if (ib >= 0
                       && tpltzblocks[ib].idv + tpltzblocks[ib].n <= idtmp) {
                igaplastblock       = ib;
                flag_igaplastinside = 0;
            } else {                      // ib=-1
                igaplastblock       = -2; // no block defined.
                flag_igaplastinside = 0;
            }
        }

        if (igapfirstblock == igaplastblock)
            distcorr_min = tpltzblocks[igapfirstblock].lambda
                         - 1; // update for lambda-1
        else
            distcorr_min = 0;

        // igapfirstblock != -2 && igaplastblock != -2 not really need but it's
        // a shortcut
        if (lgap[k] > max(distcorr_min, param_distmin) && igapfirstblock != -2
            && igaplastblock != -2) {

            idvtmp_firstblock = max(tpltzblocks[igapfirstblock].idv,
                                    id0gap[k_prev] + lgap[k_prev]);

            // test if the gap is ok for block reduce/split
            if (igapfirstblock != igaplastblock) {

                flag_gapok = 1; // reduce the gap in each block. no pb if we add
                                // max() inside the ifs.
            } else if (id0gap[k] - idvtmp_firstblock
                               >= tpltzblocks[igapfirstblock].lambda
                       && tpltzblocks[igaplastblock].idv
                                          + tpltzblocks[igaplastblock].n
                                          - (id0gap[k] + lgap[k])
                                  >= tpltzblocks[igaplastblock].lambda) {

                flag_gapok = 1;
            } else if (igapfirstblock == igaplastblock) {

                int ngappyleft_tmp  = id0gap[k] - idvtmp_firstblock;
                int leftadd         = max(0, tpltzblocks[igapfirstblock].lambda
                                                     - ngappyleft_tmp);
                int ngappyright_tmp = tpltzblocks[igaplastblock].idv
                                    + tpltzblocks[igaplastblock].n
                                    - (id0gap[k] + lgap[k]);
                int rightadd = max(0, tpltzblocks[igapfirstblock].lambda
                                              - ngappyright_tmp);
                int restgap  = lgap[k] - (leftadd + rightadd);

                //  flag_gapok = (restgap>=0);
                flag_gapok = (restgap >= max(0, param_distmin));
            } else {
                flag_gapok = 0;
            }

            // create gappy blocks if criteria is fullfill
            if (flag_gapok == 1) {

                // copy the begining blocks
                for (id = igaplastblock_prev + 1; id < igapfirstblock; id++) {

                    tpltzblocks_gappy[nblocks_gappy].T_block =
                            tpltzblocks[id].T_block;
                    tpltzblocks_gappy[nblocks_gappy].lambda =
                            tpltzblocks[id].lambda;
                    tpltzblocks_gappy[nblocks_gappy].n   = tpltzblocks[id].n;
                    tpltzblocks_gappy[nblocks_gappy].idv = tpltzblocks[id].idv;

                    nblocks_gappy = nblocks_gappy + 1;
                }

                // clear last blockgappy if same block again - outside the "if"
                // for border cases with n[]==0
                if (igaplastblock_prev == igapfirstblock && k != 0) {
                    nblocks_gappy = nblocks_gappy - 1;
                    //      idvtmp_firstblock = id0gap[k-1]+lgap[k-1]; //always
                    //      exist because igaplastblock_prev!=-1
                    // so not first turn - it's replace "idv[igapfirstblock]"
                }

                // reduce first block if defined
                if (flag_igapfirstinside == 1
                    && (id0gap[k] - idvtmp_firstblock)
                               > 0) { // check if inside and not on the border -
                                      // meaning n[] not zero

                    tpltzblocks_gappy[nblocks_gappy].T_block =
                            tpltzblocks[igapfirstblock].T_block;
                    tpltzblocks_gappy[nblocks_gappy].lambda =
                            tpltzblocks[id].lambda;
                    tpltzblocks_gappy[nblocks_gappy].n =
                            id0gap[k] - idvtmp_firstblock;
                    tpltzblocks_gappy[nblocks_gappy].n =
                            max(tpltzblocks_gappy[nblocks_gappy].n,
                                tpltzblocks[igapfirstblock].lambda);

                    tpltzblocks_gappy[nblocks_gappy].idv = idvtmp_firstblock;
                    nblocks_gappy                        = nblocks_gappy + 1;
                }

                // reduce last block if defined
                if (flag_igaplastinside == 1
                    && (tpltzblocks[igaplastblock].idv
                        + tpltzblocks[igaplastblock].n - (id0gap[k] + lgap[k]))
                               > 0) { // check if inside and not on the border -
                                      // meaning n[] not zero

                    tpltzblocks_gappy[nblocks_gappy].T_block =
                            tpltzblocks[igaplastblock].T_block;
                    tpltzblocks_gappy[nblocks_gappy].lambda =
                            tpltzblocks[id].lambda;
                    tpltzblocks_gappy[nblocks_gappy].n =
                            tpltzblocks[igaplastblock].idv
                            + tpltzblocks[igaplastblock].n
                            - (id0gap[k] + lgap[k]);
                    int rightadd0 = max(
                            0, tpltzblocks[igapfirstblock].lambda
                                       - tpltzblocks_gappy[nblocks_gappy].n);

                    tpltzblocks_gappy[nblocks_gappy].n =
                            max(tpltzblocks_gappy[nblocks_gappy].n,
                                tpltzblocks[igaplastblock].lambda);

                    tpltzblocks_gappy[nblocks_gappy].idv =
                            id0gap[k] + lgap[k] - rightadd0;

                    nblocks_gappy = nblocks_gappy + 1;
                    lambdaShftgappy =
                            lambdaShftgappy + tpltzblocks[igaplastblock].lambda;
                }

                igaplastblock_prev = igaplastblock;
                k_prev             = k;

            } // end if (flag_gapok)
        }     // end if (lgap[k]>param_distmin)
    }         // end gap loop

    // now continu to copy the rest of the block left
    for (id = igaplastblock_prev + 1; id < nb_blocks_local; id++) {

        tpltzblocks_gappy[nblocks_gappy].T_block = tpltzblocks[id].T_block;
        tpltzblocks_gappy[nblocks_gappy].lambda  = tpltzblocks[id].lambda;
        tpltzblocks_gappy[nblocks_gappy].n       = tpltzblocks[id].n;
        tpltzblocks_gappy[nblocks_gappy].idv     = tpltzblocks[id].idv;
        nblocks_gappy                            = nblocks_gappy + 1;
        lambdaShftgappy = lambdaShftgappy + tpltzblocks[id].lambda;
    }

    *nb_blocks_gappy_final = nblocks_gappy; // just for output

    return 0;
}
