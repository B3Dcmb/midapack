/**
@file toeplitz_rshp.c version 1.1b, July 2012  
@brief Contains reshaping routines to build the optimal data structure when needed for Toeplitz algebra
@author  Frederic Dauvergne
**  
** Project:  Midapack library, ANR MIDAS'09 - Toeplitz Algebra module
** Purpose:  Provide Toeplitz algebra tools suitable for Cosmic Microwave Background (CMB)
**           data analysis.
**
***************************************************************************
@note Copyright (c) 2010-2012 APC CNRS Université Paris Diderot
@note 
@note This program is free software; you can redistribute it and/or modify it under the terms
@note of the GNU Lesser General Public License as published by the Free Software Foundation; 
@note either version 3 of the License, or (at your option) any later version. This program is
@note distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
@note the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
@note Lesser General Public License for more details.
@note 
@note You should have received a copy of the GNU Lesser General Public License along with this
@note program; if not, see http://www.gnu.org/licenses/lgpl.html
@note For more information about ANR MIDAS'09 project see :
@note http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/index.html
@note ACKNOWLEDGMENT: This work has been supported in part by the French National Research 
@note Agency (ANR) through COSINUS program (project MIDAS no. ANR-09-COSI-009).
***************************************************************************
** Log: toeplitz*.c
**
** Revision 1.0b  2012/05/07  Frederic Dauvergne (APC)
** Official release 1.0beta. The first installement of the library is the Toeplitz algebra
** module.
**
** Revision 1.1b  2012/07/-  Frederic Dauvergne (APC)
** - mpi_stbmm allows now rowi-wise order per process datas and no-blocking communications.
** - OMP improvment for optimal cpu time.
** - bug fixed for OMP in the stmm_basic routine.
** - distcorrmin is used to communicate only lambda-1 datas when it is needed.
** - new reshaping routines using transformation functions in stmm. Thus, only one copy
**   at most is needed.
** - tpltz_init improvement using define_nfft and define_blocksize routines.
** - add Block struture to define each Toeplitz block.
** - add Flag structure and preprocessing parameters to define the computational strategy.
**   All the flag parameters are then available directly from the API.
**
** Revision 1.2b  2012/11/30  Frederic Dauvergne (APC)
** - extend the mpi product routine to rowwise order data distribution. This is now allowing
** tree kinds of distribution.
** - add int64 for some variables to extend the global volume of data you can use.
** - Openmp improvments.
** - Add toeplitz_wizard.c, which contains a set of easy to use routines with defined structures.
**
***************************************************************************
**
*/


#include "toeplitz.h"

//r1.1 - Frederic Dauvergne (APC)
//This is the reshaping routines to build the optimal data structure when needed. 
//The index functions find the right index number of the data location for a choosen
//transformation.


int fctid_mat2vect(int i, int id0, int n, int lambda) {
    int I, J, i_out;
    int distcorrmin = lambda - 1;
    int rfirst = id0 % n;

    if (i == -1)
        return (-1);


    I = (i + rfirst) % (n + distcorrmin);
    J = (i + rfirst) / (n + distcorrmin);

    if (I < n)
        i_out = I - rfirst + J * n;
    else
        i_out = -1; //not defined. value is zero.


    return i_out;
}


int fctid_mat2vect_inv(int i, int id0, int n, int lambda) {
    int I, J, i_out;
    int distcorrmin = lambda - 1;
    int rfirst = id0 % n;

    if (i == -1)
        i_out = -1; //not defined. value is zero.

    I = (i + rfirst) % (n);
    J = (i + rfirst) / (n);

    i_out = I - rfirst + J * (n + distcorrmin);

    return i_out;
}


int fctid_concatcol(int i, int id0, int n, int m, int l, int lconc, int lambda, int *nocol, int nbcol) {
    int I, J, i_out;
    int distcorrmin = lambda - 1;
    int rfirst = id0 % n;

    if (i == -1)
        return (-1);

    if (i >= lconc)
        return (-2);//this indice not define. It shouldn't be used

    I = (i + rfirst) % (n);
    J = (i + rfirst) / (n);

    i_out = I - rfirst + nocol[J] * (n);


    return i_out;
}


int fctid_concatcol_inv(int i, int id0, int n, int m, int l, int lconc, int lambda, int *nocol_inv, int nbcol) {
    int I, J, i_out;
    int distcorrmin = lambda - 1;
    int rfirst = id0 % n;

    if (i == -1)
        return (-1);

    if (i >= l)
        return (-2);//this indice not define. It shouldn't be used

    I = (i + rfirst) % (n);
    J = (i + rfirst) / (n);

    if (nocol_inv[J] == (-1))
        i_out = -1;
    else
        i_out = I - rfirst + nocol_inv[J] * (n);


    return i_out;
}


int fctid_vect2nfftblock(int i, int v1_size, int fft_size, int nfft, int lambda) {

    int I, J, i_out;
    int distcorrmin = lambda - 1;

    if (i == -1)
        return (-1);

    I = (i) % (fft_size);
    J = (i) / (fft_size);

    i_out = (I - distcorrmin) + J * (fft_size - 2 * distcorrmin);

    if (i_out < 0 || i_out >= v1_size)
        i_out = -1;


    return i_out;
}


int is_needconcat(int *nocol, int nbcol) {
    int i;
    int ip = nocol[0];
    for (i = 1; i < nbcol; i++) {
        if (nocol[i] != (ip + i))
            return 1;
    }


    return 0;
}


int fctid_vect2nfftblock_inv(int i, int v1_size, int fft_size, int nfft, int lambda) {

    int I, J, i_out;
    int distcorrmin = lambda - 1;

    if (i < 0 || i >= v1_size)
        return (-2);

    I = (i) % (fft_size - 2 * distcorrmin);
    J = (i) / (fft_size - 2 * distcorrmin);

    i_out = (I + distcorrmin) + J * (fft_size);

    return i_out;
}


int define_rshp_size(int flag_format_rshp, int fft_size, int nfft, int v1_size, int vedge_size, int *nrshp, int *mrshp,
                     int *lrshp) {

    if (flag_format_rshp == 2) {
        *nrshp = fft_size;
        *mrshp = nfft;
        *lrshp = (*nrshp) * (*mrshp);
    } else if (flag_format_rshp == 1) {
        *nrshp = v1_size;
        *mrshp = 1;
        *lrshp = (*nrshp) * (*mrshp);
    } else if (flag_format_rshp == 0) { //this case appear only if flag_shortcut_nbcol_eq_1==0
        *nrshp = vedge_size;
        *mrshp = 1;
        *lrshp = vedge_size;
    } else {//error not a good flag_format_rshp
    }

    return 0;
}


int build_nocol_inv(int *nocol, int nbcol, int m)  //ncol_inv to define as parameters
{
    int i;
    int *nocol_inv;
    nocol_inv = (int *) calloc(m, sizeof(double));

    for (i = 0; i < m; i++)
        nocol_inv[i] = -1;
    for (i = 0; i < nbcol; i++)
        nocol_inv[nocol[i]] = i;


    return 0;
}


int build_reshape(double *Vin, int *nocol, int nbcol, int lconc, int n, int m, int id0, int l, int lambda, int nfft,
                  double *Vrshp, int nrshp, int mrshp, int lrshp, int flag_format_rshp) {

    int i;
    int rfirst = id0 % n;
    int i_out1, i_out2, i_out3;
    int distcorrmin = lambda - 1;

    int v1_size;
    int fft_size;

    int idf = id0 + l - 1;
    int lconc0;

    FILE *file;
    file = stdout;

    v1_size = lconc + (distcorrmin) * (nbcol - 1);
    fft_size = ceil(1.0 * v1_size / nfft) + 2 * distcorrmin;

    //used transformation
    if (VERBOSE) {
        fprintf(file, "fctid_concatcol: \t %d\n", (is_needconcat(nocol, nbcol) == 1));
        fprintf(file, "fctid_mat2vect: \t %d\n", (nbcol > 1));
        fprintf(file, "fctid_vect2nfftblock \t %d\n", (nfft > 1));
    }


    for (i = 0; i < lrshp; i++) {

        if (nfft > 1)
            i_out1 = fctid_vect2nfftblock(i, v1_size, fft_size, nfft, lambda);
        else
            i_out1 = i;

        if (nbcol > 1)
            i_out2 = fctid_mat2vect(i_out1, rfirst, n, lambda);
        else
            i_out2 = i_out1;

        if (is_needconcat(nocol, nbcol) == 1)
            i_out3 = fctid_concatcol(i_out2, id0, n, m, l, lconc, lambda, nocol, nbcol);
        else
            i_out3 = i_out2;


        if (i_out3 == -1)
            Vrshp[i] = 0;
        else
            Vrshp[i] = Vin[i_out3];

    }//end for


    return 0;
}


int extract_result(double *Vout, int *nocol, int nbcol, int lconc, int n, int m, int id0, int l, int lambda, int nfft,
                   double *Vrshp, int nrshp, int mrshp, int lrshp, int flag_format_rshp) {

    int i;
    int rfirst = id0 % n;
    int i_out1, i_out2, i_out3;
    int i_in1;
    int distcorrmin = lambda - 1;

    int v1_size;
    int fft_size;

    FILE *file;
    file = stdout;

    v1_size = lconc + (distcorrmin) * (nbcol - 1);
    fft_size = ceil(1.0 * v1_size / nfft) + 2 * distcorrmin;

    //used transformation
    if (VERBOSE) {
        fprintf(file, "fctid_concatcol: \t %d\n", (is_needconcat(nocol, nbcol) == 1));
        fprintf(file, "fctid_mat2vect: \t %d\n", (nbcol > 1));
        fprintf(file, "fctid_vect2nfftblock \t %d\n", (nfft > 1));
    }

    int lcol;
    int j, k;

    for (i = 0; i < lconc; i++) {

        if (is_needconcat(nocol, nbcol) == 1)
            i_in1 = fctid_concatcol(i, id0, n, m, l, lconc, lambda, nocol, nbcol);
        else
            i_in1 = i;

        if (nbcol > 1)
            i_out2 = fctid_mat2vect_inv(i, rfirst, n, lambda);
        else
            i_out2 = i_out1;

        if (nfft > 1)
            i_out3 = fctid_vect2nfftblock_inv(i_out2, v1_size, fft_size, nfft, lambda);
        else
            i_out3 = i_out2;

        if (i_out3 == -1)
            Vout[i] = -1;
        else if (i_out3 == -2)
            Vout[i] = -2;
        else
            Vout[i_in1] = Vrshp[i_out3];
    }


    return 0;
}

