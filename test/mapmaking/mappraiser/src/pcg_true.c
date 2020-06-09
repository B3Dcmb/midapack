// Midapack library
// mapmaking code example using the Midapack library - release 1.2b, Nov 2012
// PCG routine applied to the mapmaking equation
// This can use the diagonal or the block-diagonal jacobi preconditionners

/** @file   pcg_true.c
 @author Frederic Dauvergne
 @date   November 2012
 @Last_update October 2018 by Hamza El Bouhargani*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include "midapack.h"
#include "mappraiser.h"
#include <mkl.h>
#include <assert.h>

#define MAX_LINE      (12*128*128)
#define eps 1.0e-17

#define FILE_I      "input_I.txt"
#define FILE_Q      "input_Q.txt"
#define FILE_U      "input_U.txt"


int get_pixshare_pond(Mat *A, double *pixpond);

//int precondblockjacobilike(Mat *A, Tpltz Nm1, Mat *BJ_inv, Mat *BJ, double *b, double *cond, int *lhits);


double * iqu_to_x(int lcount, const int *indices)
{
    double *xI;
    double *xQ;
    double *xU;
    double *x_real;
    
    int count;
    int i = 0;
    int ind = 0;
    
    FILE *file;
    
    xI = malloc(MAX_LINE * sizeof(double));
    xQ = malloc(MAX_LINE * sizeof(double));
    xU = malloc(MAX_LINE * sizeof(double));
    
    for (i = 0; i < MAX_LINE; i++)
    {
        xI[i] = 0.0;
        xQ[i] = 0.0;
        xU[i] = 0.0;
    }
    
    // Load I
    
    file = fopen(FILE_I, "r");
    if (file == NULL) {
        perror("Error opening file I");
        return NULL;
    }
    
    count = 0;
    while (!feof(file) && (count < MAX_LINE))
    {
        fscanf(file, "%le\n", &(xI[count]));
        count++;
    }
    fclose(file);
    
    // Load Q
    
    file = fopen(FILE_Q, "r");
    if (file == NULL) {
        perror("Error opening file Q");
        return NULL;
    }
    
    count = 0;
    while (!feof(file) && (count < MAX_LINE))
    {
        fscanf(file, "%le\n", &(xQ[count]));
        count++;
    }
    fclose(file);
    
    // Load U
    
    file = fopen(FILE_U, "r");
    if (file == NULL) {
        perror("Error opening file U");
        return NULL;
    }
    
    count = 0;
    while (!feof(file) && (count < MAX_LINE))
    {
        fscanf(file, "%le\n", &(xU[count]));
        count++;
    }
    fclose(file);
    
    // Compute x_real
    
    x_real = malloc(lcount * sizeof(double));
    for (i = 0; i < lcount; i++)
    {
        x_real[i] = 0.0;
    }
    
    // Construct the local exact solution method 2, what is the relation between local and global indices ?
    for (i = 0; i < lcount; i++)
    {
        ind = indices[i]/3;
        if (i%3 == 0)
            x_real[i] = xI[ind];
        else if (i%3 == 1)
            x_real[i] = xQ[ind];
        else
            x_real[i] = xU[ind];
    }
    
    free(xI);
    free(xQ);
    free(xU);
    
    return x_real;
}



void Anorm(double *x_real, double *x, int rank, int k, double *pixpond, Mat *A, Tpltz Nm1)
{
    
    int i = 0;
    int n = 0;
    int m = 0;
    //Declaration of variables of anorm
    double *d;
    double *anorm, *anormt;
    double anormc = 0.0;
    double result = 0.0;
    FILE* file2;
    char filenametxt2 [1024];
    
    n = A->lcount;
    m = A->m;
    d =  malloc(n*sizeof(double));
    anorm = (double *) malloc(m*sizeof(double));
    anormt = (double *) malloc(n*sizeof(double));
    
    //anorm added
    for (i=0; i<n; i++)
        d[i] = fabs(x_real[i]- x[i]);
    
    MatVecProd(A, d, anorm, 0);
    stbmmProd(Nm1, anorm);
    TrMatVecProd(A, anorm, anormt, 0);
    
    
    //Multiplication of pixpond because the vector is overlapped
    for (i=0; i<n; i++)
        anormc += d[i] * anormt[i] *  pixpond[i];
    
    //Use Allreduce, because all processors compute a part of result(Anorm)
    MPI_Allreduce(&anormc, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    //Print the anormc which corresponds to the given iteration, than compare in 2 different precondtioners
    if (rank == 0)
    {
        
        sprintf(filenametxt2, "BD_Anorm%01d.txt", Nm1.tpltzblocks->lambda);
        
        file2 = fopen(filenametxt2, "a");
        
        fprintf(file2, "%d, %25.18e\n", k, result);
        fclose(file2);
        
        //printf("************************\nIteration k = %d, the Anorm is = %le\n", k, result);
    }
    
    free(anorm);
    free(anormt);
    
}



void transpose_nn(double *A, int n)
{
    int i, j;
    double temp;
    
    for(i = 0; i < n-1 ; i++)
        for (j = i+1; j < n; j++)
        {
            temp = A[i*n+j];
            A[i*n+j] = A[j*n+i];
            A[j*n+i] = temp;
        }
}



void inverse_svd(int m, int n, int lda,  double *a)
{
    // A = UDVt
    // A = USVt
    int i, j, k;
    
    // Setup a buffer to hold the singular values:
    int nsv = m < n ? m : n;
    double *s = malloc(nsv * sizeof(double)); // = D
    
    // Setup buffers to hold the matrices U and Vt:
    double *u = malloc(m*m * sizeof(double)); // = U
    double *vt = malloc(n*n * sizeof(double));
    
    // Workspace and status variables:
    /*
     double workSize;//
     double *work = &workSize;//
     int lwork = -1;//
     int *iwork = malloc(8 * nsv * sizeof(int));//
     
     double *superb = malloc((nsv-1) * sizeof(double));
     */
    
    int info = 0;
    
    transpose_nn(a, m);
    mkl_set_num_threads(1);
    
    info = LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'A', m, n, a, lda, s, u, m, vt, n);
    //info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'A', 'A', m, n, a, lda, s, u, m, vt, n, superb);
    /*
     
     // Call dgesdd_ with lwork = -1 to query optimal workspace size:
     //dgesvd_("A", "A", &m, &n, a, &lda, s, u, &m, vt, &n, work, &lwork, &info);
     dgesdd_("A", &m, &n, a, &lda, s, u, &m, vt, &n, work, &lwork, iwork, &info);
     // info = LAPACKE_dgesdd_work(LAPACK_COL_MAJOR, 'A', m, n, a, lda, s, u, m, vt, n, work, lwork, iwork);
     
     //info = LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'A', m, n, a, lda, s, u, m, vt, n);
     if (info > 0) {
     printf( "The algorithm computing SVD failed to converge (1).\n" );
     exit(1);
     }
     if (info < 0)
     { printf( "General error .\n" );
     exit(1);
     }
     
     // Optimal workspace size is returned in work[0].
     lwork = workSize;
     work = malloc(lwork * sizeof(double));
     
     // Call dgesdd_ to do the actual computation:
     //dgesvd_("A", "A", &m, &n, a, &lda, s, u, &m, vt, &n, work, &lwork, &info);
     dgesdd_("A", &m, &n, a, &lda, s, u, &m, vt, &n, work, &lwork, iwork, &info);
     
     */
    
    if (info > 0) {
        printf( "The algorithm computing SVD failed to converge (1).\n" );
        exit(1);
    }
    if (info < 0)
    { printf( "General error .\n" );
        exit(1);
    }
    
    // Cleanup workspace:
    // free(work);
    // free(iwork);
    
    // Computing S-1
    for (k = 0; k < nsv; ++k) {
        if (fabs(s[k]) < 1.0e-7) s[k] = 0.0;
        //if (fabs(s[k]) < 1.0e-12) s[k] = 0.0;
        else s[k] = 1.0 / s[k];
    }
    
    // do something useful with U, S, Vt ...
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            a[i * m + j] = 0.0;
            for (k = 0; k < n; k++) {
                a[i * m + j] += vt[i * m + k] * s[k] * u[k * m + j];
            }
        }
    }
    
    // and then clean them up too:
    free(s);
    free(u);
    free(vt);
    
}

//in vector x overlapped,
//in P distributed among processors
//Z should be contructed in place
// Calculation of E, then E^{-1}, we have the algorithm
//Pt N^{-1} P Z (Zt Pt N^{-1} P Z)^{-1} Zt x = A Z E^{-1} Zt x = A Q x
//v1 = P * Zi column by column
//v2 = N^{-1} * v1
//v3 = Pt * v2
//Ei = Zt * v3, Zt is entire matrix and Ei is computed column by column, Zi is a column




double ** build_Z_random(Mat *A, int Zn, int n, int pflags)
{
    double **Z; // Zn * A->lcount, pointers to columns (ie column-major)
    int j,k;
    
    Z = calloc(Zn, sizeof(double *));
    for (j = 0; j < Zn; j++)
    {
        Z[j] = calloc(n, sizeof(double));
        for (k = 0; k < n; k++) //A.lindices[j] is the Z is dense and not overlapped
            Z[j][k] = 1.0; // TO COMPUTE
    }
    
    return Z;
}



/** Function m2m_sum for "sum map to map"
 Extract values from one map (A1, vA1), and for each pixel shared with an other map (A2, vA2),
 sum pixel value in vA1 to pixel value in vA2.
 @return a number of elements shared between A1 and A2
 @sa m2m
 @ingroup matmap_group22*/
int m2m_sum_i(int *vA1, int *A1, int n1, int *vA2, int *A2, int n2){
    int i=0, j=0, k= 0;
    while( i<n1 && j<n2){
        if(A1[i] < A2[j]){
            i++;
        }
        else if(A1[i] > A2[j]){
            j++;
        }
        else{
            vA2[j]+=vA1[i];
            k++;
            i++;
            j++;
        }
    }
    return k;
}




int m2m_sum_d(double *vA1, int *A1, int n1, double *vA2, int *A2, int n2){
    int i=0, j=0, k= 0;
    while( i<n1 && j<n2){
        if(A1[i] < A2[j]){
            i++;
        }
        else if(A1[i] > A2[j]){
            j++;
        }
        else{
            vA2[j] += vA1[i];
            k++;
            i++;
            j++;
        }
    }
    return k;
}

double **build_Z(Mat *A, int Zn, int pflags)
{
    int i, j, e, k, rank, size;
    int p, rp, sp, tag = 0;
    MPI_Request s_request, r_request;
    MPI_Status status;
    
    int lcount_max = 0, rlcount;
    int *count, *tcount;
    int *rcount, *rindices;
    double *rZ;
    double **Z; // Zn * A->lcount, pointers to columns (ie column-major)
    
    
    MPI_Comm_rank(A->comm, &rank);                //get rank and size of the communicator
    MPI_Comm_size(A->comm, &size);
    
    
    // Calculer le lcount_max (tous proc)
    MPI_Allreduce(&(A->lcount), &(lcount_max), 1, MPI_INT, MPI_MAX, A->comm);
    
    // Allouer les buffers
    count = (int *)calloc(A->lcount, sizeof(int)); //the number of appereance in local processor
    tcount = (int *)calloc(A->lcount, sizeof(int));
    rcount = (int *)malloc(lcount_max * sizeof(int)); // the number of appereance in neighbor processors
    rindices  = (int *)malloc(lcount_max * sizeof(int)); //real indices in neighbor processors
    
    
    // Calculer count local pour un pixel donnÃ©
    for (i = 0; i < A->m * A->nnz; i++)
        count[A->indices[i]]++;
    
    // Recopier count local dans total
    for (i = 0; i < A->m * A->nnz; i++)
        tcount[A->indices[i]] = count[A->indices[i]];
    
    // Compute total counts
    
    for (p=1; p < size; p++) {    //loop : collective global reduce in ring-like fashion
        rp = (size + rank - p)%size;
        sp = (rank + p)%size;
        MPI_Send(&(A->lcount), 1, MPI_INT, sp, 0, A->comm);                //exchange sizes
        MPI_Recv(&rlcount, 1, MPI_INT, rp, 0, A->comm, &status);
        tag++;
        MPI_Irecv(rindices, rlcount, MPI_INT, rp, tag, A->comm, &r_request);//exchange global indices
        MPI_Isend(A->lindices, A->lcount, MPI_INT, sp, tag, A->comm, &s_request);
        MPI_Wait(&r_request, &status);
        MPI_Wait(&s_request, &status);
        tag++;
        MPI_Irecv(rcount, rlcount, MPI_INT, rp, tag, A->comm, &r_request);    //exchange local count/rcount values
        MPI_Isend(count, A->lcount, MPI_INT, sp, tag, A->comm, &s_request);
        tag++;
        MPI_Wait(&r_request, &status);
        m2m_sum_i(rcount, rindices, rlcount, tcount, A->lindices, A->lcount);        //sum in the result
        MPI_Wait(&s_request, &status);
    }
    
    // LibÃ©rer les buffers non utilisÃ©s
    free(rcount);
    
    
    // Allouer Z
    Z = calloc(size, sizeof(double *));
    
    // Calculer le Z correspondant au processus courant
    Z[rank] = calloc(A->lcount, sizeof(double));
    for (i = 0; i < A->lcount; i += A->nnz)
        Z[rank][i] = (double)((double)count[i] / (double)tcount[i]);
    
    
    // LibÃ©rer les buffers non utilisÃ©s
    free(count);
    free(tcount);
    
    // Allouer le buffer pour Ã©changer Z
    rZ = (double *)malloc(lcount_max * sizeof(double));
    
    // Echanger Z
    
    for (p=1; p < size; p++) {    //loop : collective global reduce in ring-like fashion
        rp = (size + rank - p)%size;
        sp = (rank + p)%size;
        MPI_Send(&(A->lcount), 1, MPI_INT, sp, 0, A->comm);                //exchange sizes
        MPI_Recv(&rlcount, 1, MPI_INT, rp, 0, A->comm, &status);
        tag++;
        MPI_Irecv(rindices, rlcount, MPI_INT, rp, tag, A->comm, &r_request);    //exchange global indices
        MPI_Isend(A->lindices, A->lcount, MPI_INT, sp, tag, A->comm, &s_request);
        MPI_Wait(&r_request, &status);
        MPI_Wait(&s_request, &status);
        tag++;
        MPI_Irecv(rZ, rlcount, MPI_DOUBLE, rp, tag, A->comm, &r_request);    //exchange local values
        MPI_Isend(Z[rank], A->lcount, MPI_DOUBLE, sp, tag, A->comm, &s_request);
        tag++;
        MPI_Wait(&r_request, &status);
        Z[rp] = calloc(A->lcount, sizeof(double));
        m2m(rZ, rindices, rlcount, Z[rp], A->lindices, A->lcount);        //copy the interesting value the corresponding Z[rp] value, thanks!!!!!!!!
        MPI_Wait(&s_request, &status);
    }
    
    
    int ratio = size / Zn;
    assert(size >= Zn);
    assert(Zn * ratio == size);
    
    double **ZZ = calloc(Zn, sizeof(double *));
    for (j=0; j < Zn; j++)
    {
        ZZ[j] = Z[j * ratio];
        for (k = 1; k < ratio; k++)
        {
            for (i = 0; i < A->lcount; i++)
            {
                ZZ[j][i] += Z[j * ratio + k][i];
            }
            free(Z[j * ratio + k]);
        }
    }
    free(Z);
    
    for(j=0; j < Zn; j++)
    {
        for (i = 0; i < A->lcount - A->nnz * A->trash_pix; i++)
        {
            ZZ[j][i] = ZZ[j][i + (A->nnz)*(A->trash_pix)];
        }
    }
    
    
    // LibÃ©rer les buffers non utilisÃ©s
    free(rindices);
    free(rZ);
    
    return ZZ;
}



void build_Em1(Mat *A, Tpltz Nm1, double **Z, double *pixpond, int Zn, int n, int precond, double **out_E, double **AZ, int pflags)
{
    double *E, *EO;
    int i,j,k;
    
    int rank;
    MPI_Comm_rank(A->comm, &rank);
    
    E = calloc(Zn * Zn, sizeof(double));
    
    // Ei = Zt * Pt * N^{-1} * P * Zi
    // Em1 = E^{-1}
    
    for (i = 0; i < Zn; i++)
    {
        // E = Zt * v3
        for (j = 0; j < Zn; j++)
        {
            for (k = 0; k < n; k++)
            {
                E[i * Zn + j] += (double)(Z[j][k] * AZ[i][k] * pixpond[k]);
            }
        }
    }
    
    MPI_Allreduce(MPI_IN_PLACE, E, Zn * Zn, MPI_DOUBLE, MPI_SUM, A->comm);
    
    int info;
    double anorm, rcond, *w;
    int *iw;
    iw = calloc(Zn, sizeof(int));
    w = calloc(Zn*Zn*2, sizeof(double));
    
    inverse_svd(Zn, Zn, Zn, E);
    
    EO = calloc(Zn * Zn, sizeof(double));
    
    memcpy(EO, E, sizeof(double) * Zn * Zn);
    
    /* Computes the norm of x */
    anorm = dlange_("1", &Zn, &Zn, EO, &Zn, w);
    
    /* Modifies x in place with a LU decomposition */
    dgetrf_(&Zn, &Zn, EO, &Zn, iw, &info);
    // if (info != 0) fprintf(stderr, "failure with error %d\n", info);
    
    /* Computes the reciprocal norm */
    dgecon_("1", &Zn, EO, &Zn, &anorm, &rcond, w, iw, &info);
    // if (info != 0) fprintf(stderr, "failure with error %d\n", info);
    
    printf("condition number of Einv = %25.18e\n", rcond);
    
    free(EO);
    
    *out_E = E;
}


void build_AZ(Mat *A, Tpltz Nm1, double **Z, int Zn, int n, double ***out_AZ)
{
    double **AZ;
    double *v;
    int i,j,k;
    
    v = calloc(A->m, sizeof(double)); // P * Zi
    
    AZ = calloc(Zn, sizeof(double *));
    for (k = 0; k < Zn; k++)
    {
        AZ[k] = calloc(n, sizeof(double));
    }
    
    for (i = 0; i < Zn; i++)
    {
        // AZ[i] = Pt N^{-1} P * Z[i]
        MatVecProd(A, Z[i], v, 0);
        stbmmProd(Nm1, v); // In-place
        TrMatVecProd(A, v, AZ[i], 0);
    }
    
    free(v);
    
    *out_AZ = AZ;
}



typedef struct _buffers {
    double *w1;
    double *w2;
    double *w3;
} buffers;

buffers * init_buffers(Mat *A, int Zn, int n)
{
    buffers *buf = malloc(sizeof(buffers));
    buf->w1 = calloc(Zn, sizeof(double)); // Zt * x
    buf->w2 = calloc(Zn, sizeof(double)); // Em1 * Zt * x
    //buf->w3 = calloc(n, sizeof(double)); // Z * Em1 * Zt * x
    return buf;
}

void free_buffers(buffers *buf)
{
    free(buf->w1);
    free(buf->w2);
    //free(buf->w3);
    free(buf);
}

double* build_Qx(Mat *A, Tpltz Nm1, double **Z, double *Em1, double *x, double *pixpond, buffers *buf, int Zn, int n, int pflags)
{
    int i,j,k;
    
    int rank;
    MPI_Comm_rank(A->comm, &rank);
    
    
    // Pt N^{-1} P Z (Zt Pt N^{-1} P Z)^{-1} Zt x
    // Pt N^{-1} P Z (E)^{-1} Zt x
    // Pt N^{-1} P Q x
    
    // w1 = Zt x
    for (j = 0; j < Zn; j++)
    {
        buf->w1[j] = 0.0;
        for (k = 0; k < n; k++)
        {
            buf->w1[j] += Z[j][k] * x[k] * pixpond[k];
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, buf->w1, Zn, MPI_DOUBLE, MPI_SUM, A->comm);
    
    
    // w2 = Em1 * w1 (dense matrix E (Zn*Zn) times a vector w1 (Zn*1)
    for (i = 0; i < Zn; i++)
    {
        buf->w2[i] = 0.0;
        for (j = 0; j < Zn; j++)
        {
            buf->w2[i] += Em1[i * Zn + j] * buf->w1[j];
        }
    }
    
    // w3 = Z * w2 (overlapped times dense);
    for (k = 0; k < n; k++)
    {
        buf->w3[k] = 0.0;
        for (j = 0; j < Zn; j++)
        {
            buf->w3[k] += Z[j][k] * buf->w2[j];
        }
    }
    
    // Vector W3 is similar to vector x, so we can continue the multiplication
    return buf->w3;
}




void build_Qtx(Mat *A, Tpltz Nm1, double **Z, double *Em1, double *x, double *pixpond, buffers *buf, int Zn, int n, int pflags)
{
    int i,j,k;
    
    int rank;
    MPI_Comm_rank(A->comm, &rank);
    
    
    // Pt N^{-1} P Z (Zt Pt N^{-1} P Z)^{-1} Zt x
    // Pt N^{-1} P Z (E)^{-1} Zt x
    // Pt N^{-1} P Q x
    
    // w1 = Zt x
    for (j = 0; j < Zn; j++)
    {
        buf->w1[j] = 0.0;
        for (k = 0; k < n; k++)
        {
            buf->w1[j] += Z[j][k] * x[k] * pixpond[k];
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, buf->w1, Zn, MPI_DOUBLE, MPI_SUM, A->comm);
    
    
    // w2 = Em1 * w1 (dense matrix E (Zn*Zn) times a vector w1 (Zn*1)
    for (i = 0; i < Zn; i++)
    {
        buf->w2[i] = 0.0;
        for (j = 0; j < Zn; j++)
        {
            buf->w2[i] += Em1[i * Zn + j] * buf->w1[j];
        }
    }
    
    // return buf->w2;
}

void mul_ZQtx(double **Z, buffers *buf, double *vec, int Zn, int n, int pflags)
{
    
    int k, j;
    
    // w3 = Z * w2 (overlapped times dense);
    for (k = 0; k < n; k++)
    {
        vec[k] = 0.0;
        for (j = 0; j < Zn; j++)
        {
            vec[k] += Z[j][k] * buf->w2[j];
        }
    }
    
}



void PCG_Lanczos_eig(Mat *A, Tpltz Nm1, Mat *BJ, Mat *BJ_true, double *x, int *lhits, double *cond, double *b, double *noise, double tol, double *pixpond, int K, double ***out_Ritz_vectors, double ***out_Ritz_vectors_AZ)
{
    
    int i, j, k ;            // some indexes
    int m, n, rank, size;
    double st, t;               //timers
    double solve_time = 0.0;
    double beta, alpha, result, dot;
    int info = 0, lwork = -1;
    
    double *Av = NULL, *_g = NULL;
    double *Tt = NULL, *T = NULL;
    double *w = NULL, *v = NULL, *vold = NULL, *tmp = NULL;
    double *V = NULL, *BJw = NULL, *AmulV = NULL;
    
    double *Ritz_values = NULL;
    double *Ritz_vectors_out = NULL;
    double *Ritz_vectors_out_r = NULL;
    double **Ritz_vectors = NULL;
    double **Ritz_vectors_AZ = NULL;
    
    double *work = NULL;
    double wkopt = 0.0;
    
    
    
    m=A->m;                    // Number of local time samples
    //n=A->lcount;               // Number of local pixels
    
    
    n=(A->lcount)-(A->nnz)*(A->trash_pix);
    
    MPI_Comm_rank(A->comm, &rank);
    MPI_Comm_size(A->comm, &size);
    
    st=MPI_Wtime();
    
    tmp = (double *)calloc(n, sizeof(double));
    
    // Map domain
    Av = (double *)malloc(n*sizeof(double));
    _g = (double *)malloc(m*sizeof(double));
    
    T = (double *)calloc((K+1)*(K+1), sizeof(double));
    Tt = (double *)calloc(K*K, sizeof(double));
    
    
    
    Ritz_vectors_out = (double *)calloc(n*K,  sizeof(double));
    // Ritz_vectors_out_r = (double *)calloc(n*K,  sizeof(double));
    w = (double *)malloc(n*sizeof(double));
    v = (double *)calloc(n, sizeof(double));
    vold = (double *)calloc(n, sizeof(double));
    V = (double *)calloc(n * (K+1), sizeof(double));
    
    AmulV = (double *)calloc(n * (K+1), sizeof(double));
    
    
    BJw = (double *)malloc(n*sizeof(double));
    
    Ritz_values = (double *)calloc(K, sizeof(double));
    
    Ritz_vectors = calloc(K, sizeof(double *));
    for (i = 0; i < K; i++)
        Ritz_vectors[i] = calloc(n, sizeof(double));
    
    
    Ritz_vectors_AZ = calloc(K, sizeof(double *));
    for (i = 0; i < K; i++)
        Ritz_vectors_AZ[i] = calloc(n, sizeof(double));

    
    
    //srand (time(NULL));
    for (i = 0; i < n; ++i)
        w[i] = 0.0;
    //w[i] = (double)rand() / (RAND_MAX);
    
    MatVecProd(A, x, _g, 0);
    
    for(i=0; i<m; i++)
        
        _g[i] = b[i] + noise[i] - _g[i];
    
    stbmmProd(Nm1, _g);
    
    TrMatVecProd(A, _g, w, 0);
    
    //MatVecProd(BJ_true, w, BJw, 0);
    memcpy(BJw, w, n * sizeof(double));
    
    
    
    //beta = sqrt(dot(w, BJw))
    dot = 0.0;
    for (i = 0; i < n; i++)
        dot += w[i] * BJw[i] * pixpond[i];
    MPI_Allreduce(&dot, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    beta = sqrt(result);
    
    
    if (rank == 0){
        printf("Before Iteration 0\n");
        printf("beta = %e\n", beta);
    }
    
    if (beta > eps){
        for (i = 0; i < n; i++){
            v[i] = w[i]/beta;
            V[i*(K+1)] = v[i];
        }
    }
    
    
    //Av = A * v = Pt N P * v
    MatVecProd(A, v, _g, 0);
    stbmmProd(Nm1, _g);
    TrMatVecProd(A, _g, Av, 0);
    
    
    // w = inv(BD) * Av
    //MatVecProd(BJ, Av, w, 0);
    memcpy(w, Av, n * sizeof(double));
    
    //alpha = dot(v, Av)
    dot = 0.0;
    for (i = 0; i < n; i++)
        dot += v[i] * Av[i] * pixpond[i];
    MPI_Allreduce(&dot, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    if (rank == 0)
        printf("alpha = %e\n", alpha);
    
    for (i = 0; i < n; i++)
    {
        AmulV[i*(K+1)] = w[i];
        w[i] = w[i] - (alpha * v[i]);
    }
    
    
    t=MPI_Wtime();
    if (rank==0) {
        printf("[rank %d] Lanczos init time=%lf \n", rank, t-st);
        fflush(stdout);
    }
    
    
    st=MPI_Wtime();
    
    for (i = 0; i < K; i++)
    {
        if (rank == 0) printf("Iteration %d\n", i);
        
        //MatVecProd(BJ_true, w, BJw, 0);
        memcpy(BJw, w, n * sizeof(double));
        
        dot = 0.0;
        for (j = 0; j < n; j++)
            dot += w[j] * BJw[j] * pixpond[j];
        MPI_Allreduce(&dot, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        beta = sqrt(result);
        
        if (rank == 0)
            printf("beta = %e\n", beta);
        
        if (beta > eps)
        {
            for (j = 0; j < n; j++)
            {
                v[j] = w[j] / beta;
                V[j * (K+1) + i + 1] = v[j];
            }
        }
        
        else if (rank == 0) printf("division by zero in iteration %d\n", i);
        
        // What should we do to construct this special triangular matrix
        //T(i,i) = alpha
        //T(i,i+1) = beta
        //T(i+1,i) = beta
        
        T[(i * (K+1)) + i] = alpha;
        T[(i * (K+1)) + i + 1] = beta;
        T[((i+1) * (K+1)) + i] = beta;
        
        //Av = A * v = Pt N P * v
        MatVecProd(A, v, _g, 0);
        stbmmProd(Nm1, _g);
        TrMatVecProd(A, _g, Av, 0);
        
        // w = inv(BD) * Av
        //MatVecProd(BJ, Av, w, 0);
        memcpy(w, Av, n * sizeof(double));
        
        
        //alpha = dot(v, Av)
        dot = 0.0;
        for (j = 0; j < n; j++)
            dot += v[j] * Av[j] * pixpond[j];
        MPI_Allreduce(&dot, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        if (rank == 0)
            printf("alpha = %e\n", alpha);
        
        for (j = 0; j < n; j++)
        {
            AmulV[j * (K+1) + i + 1] = w[j];
            w[j] = w[j] - (alpha * v[j]) - (beta * vold[j]);
            vold[j] = v[j];
            
        }
        
        t = MPI_Wtime();
        if (rank==0) {
            printf("Iteration = %d, [rank %d] Lanczos constructing time=%lf \n", i, rank, t-st);
            fflush(stdout); }
        
        st=MPI_Wtime();
    }
    
    
    // Here we reduce the dimention of T from (K+1 * K+1) to K * K;
    // Here we reduce the dimention of V from (N * K+1) to (N * K)
    
    for (i = 0; i < K; i++){
        for (j = 0; j < K; j++){
            T[i*K + j] = T[i*(K+1) + j];
        }
    }
    // memcpy(Tt, T, K * K * sizeof(double));
    
    
    
    for (i = 0; i < n; i++){
        for (j = 0; j < K; j++){
            AmulV[i*K + j] = AmulV[i*(K+1) + j];
        }
    }
    
    for (i = 0; i < n; i++){
        for (j = 0; j < K; j++){
            V[i*K + j] = V[i*(K+1) + j];
        }
    }
    
    
    
    
    //[Ritz_vectors, Ritz_values] = eig(T)
    //Ritz_values contains the eigenvalues of the matrix A in ascending order.
    
    //lapack_int LAPACKE_dsyev(int matrix_layout, char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w);
    
    st = MPI_Wtime();
    
    transpose_nn(T, K);
    
    //int dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info)
    
    
    
    dsyev_("Vectors", "Upper", &K, T, &K , Ritz_values, &wkopt, &lwork, &info);
    
    
    
    
    lwork = (int) wkopt;
    work = (double*)malloc( lwork*sizeof(double));
    
    dsyev_("Vectors", "Upper", &K, T, &K , Ritz_values, work, &lwork, &info);
    
    
    //free(work);
    
    
    
    transpose_nn(T, K);
    
    
    t=MPI_Wtime();
    if (rank==0) {
        printf("[rank %d] Lanczos dsyev time=%lf \n", rank, t-st);
        fflush(stdout);
    }
    
    
    
    st =MPI_Wtime();
    
    memset(Ritz_vectors_out, 0, n*K*sizeof(double));
    
    
    
    for (i = 0; i < n; i++){
        for (k = 0; k < K; k++) {
            for (j = 0; j < K; j++){
                Ritz_vectors_out[i*K + j] += V[i*K + k] * T[k*K + j];
            }
        }
    }
    
    for (i = 0; i < n; i++)
        for (j = 0; j < K; j++)
            Ritz_vectors[j][i] =  Ritz_vectors_out[i*K + j];
    
    
    
    memset(Ritz_vectors_out, 0, n*K*sizeof(double));
    
    for (i = 0; i < n; i++){
        for (k = 0; k < K; k++) {
            for (j = 0; j < K; j++){
                Ritz_vectors_out[i*K + j] += AmulV[i*K + k] * T[k*K + j];
            }
        }
    }
    
    for (i = 0; i < n; i++)
        for (j = 0; j < K; j++)
            Ritz_vectors_AZ[j][i] =  Ritz_vectors_out[i*K + j];
    
    
    
    
    t=MPI_Wtime();
    if (rank==0) {
        printf("[rank %d] Lanczos V*T multiplication time=%lf \n", rank, t-st);
        fflush(stdout);
    }
    
    free(Av); free(_g); free(Tt); free(T); free(w); free(v); free(vold);
    free(V); free(BJw); free(Ritz_values); free(Ritz_vectors_out);
    
    *out_Ritz_vectors = Ritz_vectors;
    *out_Ritz_vectors_AZ = Ritz_vectors_AZ;
}



extern double stbmm_time;
extern double trmv_time;
extern double mv_time;
extern double greedy_time;


int PCG_GLS_true(char *outpath, char *ref, Mat *A, Tpltz Nm1, double *x, double *b, double *noise, double *cond, int *lhits, double tol, int K, int precond)
{
    
    int         i, j, k ;            // some indexes
    int        m, n, rank, size;
    double     localreduce;            //reduce buffer
    double    st, t;                //timers
    double    solve_time = 0.0;
    double    res, res0,  *res_rel;
    double *tmp;
    FILE *fp;
    
    res_rel = (double *) malloc(1*sizeof(double));
    
    m=A->m;                    //number of local time samples
    //n=(A->lcount)-(A->nnz)*(A->trash_pix);           //number of local pixels
    MPI_Comm_rank(A->comm, &rank);            //
    MPI_Comm_size(A->comm, &size);            //
    
    
    double *_g, *ACg, *Ah, *Nm1Ah ;        // time domain vectors
    double *g, *gp, *gt, *Cg, *h;//                     // map domain vectors
    double *AtNm1Ah ;                            // map domain
    double ro, gamma, coeff ;            // scalars
    double g2pix, g2pixp;
    double norm2b;
    
    double *x_real;
    double g2pix_polak;
    
    Mat BJ, BJ_true;
    
    st=MPI_Wtime();
    
    precondblockjacobilike(A, Nm1, &BJ, &BJ_true, b, cond, lhits);
    
    
    n=(A->lcount)-(A->nnz)*(A->trash_pix);
    // Reallocate memory for well-conditioned map
    tmp = realloc(x, n * sizeof(double));
    if(tmp !=NULL){
        x = tmp;
    }
    
    //map domain
    h = (double *) malloc(n*sizeof(double));      //descent direction
    g = (double *) malloc(n*sizeof(double));
    gp = (double *) malloc(n*sizeof(double));
    AtNm1Ah = (double *) malloc(n*sizeof(double));
    
    //time domain
    Ah = (double *) malloc(m*sizeof(double));
    
    
    _g = Ah;
    Cg = AtNm1Ah;
    Nm1Ah = Ah;
    
    printf("n=%d, m=%d, A.nnz=%d \n", n, m, A->nnz );
    
    // For M2lvl
    double **Z; // Zn * A->lcount, pointers to columns (ie column-major)
    double *Em1;
    double **AZ;
    double *Qg;
    double *AQg;
    buffers *buf;
    int Zn = size; // Set the size of Z
    
    buf = init_buffers(A, Zn, n);
    
    double *pixpond;
    pixpond = (double *) malloc(n*sizeof(double));
    
    //compute pixel share ponderation
    get_pixshare_pond( A, pixpond);
    
    //if we want to use the true norm to compute the residual
    int TRUE_NORM=1;  //0: No ; 1: Yes
    
    t=MPI_Wtime();
    if (rank==0) {
        printf("[rank %d] 2lvl compute BD time=%lf \n", rank, t-st);
        fflush(stdout);
    }

    // 2lvl-specific code
    if (precond != 0) {
        Qg = calloc(n, sizeof(double));
        AQg = calloc(n, sizeof(double));
    
        st=MPI_Wtime();
        
        // For M2lvl
        //Z = build_Z_random(A, Zn, n, 0);
        if (precond == 1) {
            Z = build_Z(A, Zn, 0); // 2lvl a priori preconditioner
            build_AZ(A, Nm1, Z, Zn, n, &AZ);
        }
        else if (precond == 2)
            PCG_Lanczos_eig(A, Nm1, &BJ, &BJ_true, x, lhits, cond, b, noise, tol, pixpond, Zn, &Z, &AZ); // 2lvl a posteriori preconditioner
        else {
            printf("Unknown preconditioner: %d\n", precond);
            exit(1);
        }

        t=MPI_Wtime();
        if (rank==0) {
            printf("[rank %d] 2lvl compute Z time=%lf \n", rank, t-st);
            fflush(stdout);
        }
        
        st=MPI_Wtime();
        
        //Em1 = build_Em1(A, Nm1, Z, pixpond, Zn, n, precond, 0);
        if (precond != 0) build_Em1(A, Nm1, Z, pixpond, Zn, n, precond, &Em1, AZ, 0);
        
        t=MPI_Wtime();
        if(rank==0){
            printf("2lvl compute Em1 time=%lf \n", t-st);
            fflush(stdout);
        }

        //memcpy(x, Z[Zn-1], n * sizeof(double));
        memcpy(x, Z[Zn-1], n * sizeof(double));
    }
    
    st=MPI_Wtime();
    
    MatVecProd(A, x, _g, 0);
    
    for(i=0; i<m; i++)
        _g[i] = b[i] + noise[i] - _g[i];
    
    stbmmProd(Nm1, _g);        // _g = Nm1 (Ax-b)
    
    TrMatVecProd(A, _g, g, 0);        //  g = At _g
    
    
    if (precond == 0) {
        MatVecProd(&BJ, g, Cg, 0);
    }
    // 2lvl-specific code
    else {
        build_Qtx(A, Nm1, Z, Em1, g, pixpond, buf, Zn, n, 0);
        mul_ZQtx(Z, buf, Qg, Zn, n, 0);
        mul_ZQtx(AZ, buf, AQg, Zn, n, 0);
        
        for (i = 0; i < n; ++i) {
            AQg[i] = g[i] - AQg[i];
        }
        
        MatVecProd(&BJ, AQg, Cg, 0);
        
        for (i = 0; i < n; ++i) {
            Cg[i] += Qg[i];
        }
    }
    
    for(j=0; j<n; j++)                   //  h = -Cg
        h[j]=Cg[j];
    
    g2pix=0.0;                               //g2 = "res"
    localreduce=0.0;
    for(i=0; i<n; i++)                    //  g2 = (Cg , g)
        localreduce += Cg[i] * g[i] * pixpond[i];
    
    MPI_Allreduce(&localreduce, &g2pix, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    t=MPI_Wtime();
    solve_time += (t-st);
    
    
    //Just to check with the true norm:
    if (TRUE_NORM==1) {
        res=0.0;                               //g2 = "res"
        localreduce=0.0;
        for(i=0; i<n; i++)                    //  g2 = (Cg , g)
            localreduce += g[i] * g[i] * pixpond[i];
        
        MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    else {
        
        res = g2pix;
        
    }//end if
    
    double g2pixB = g2pix;
    double tol2rel = tol*tol*res;//tol*tol*g2pixB;//*g2pixB;//tol;//*tol*g2;
    res0=res;
    //Test if already converged
    if(rank==0) {
        printf("res = %e, res0 = %e, g2pix = %e, g2pixB = %e\n", res, res0, g2pix, g2pixB);
        *res_rel = sqrt(res)/sqrt(res0);
        printf("k=%d res_g2pix=%e res_g2pix_rel=%e res_rel=%e time=%lf\n", 0, g2pix , sqrt(g2pix)/sqrt(g2pixB), sqrt(res)/sqrt(res0), t-st);
        char filename[256];
        sprintf(filename,"%s/pcg_residuals_%s.dat",outpath, ref);
        fp=fopen(filename, "wb");
        fwrite(res_rel, sizeof(double), 1, fp);
    }
    
    if(res <= tol){                         //
        if(rank==0)                      //
            printf("--> converged (%e < %e)\n", res, tol);
        k=K;//to not enter inside the loop
    }
    
    st=MPI_Wtime();
    fflush(stdout);
    
    stbmm_time = 0.0;
    mv_time = 0.0;
    trmv_time = 0.0;
    greedy_time = 0.0;
    
    // PCG Descent Loop *********************************************
    for(k=1; k<K ; k++){
        
        // Swap g backup pointers (Ribière-Polak needs g from previous iteration)
        gt = gp;
        gp = g;
        g = gt;
        
        MatVecProd(A, h, Ah, 0);        // Ah = A h
        stbmmProd(Nm1, Nm1Ah);        // Nm1Ah = Nm1 Ah   (Nm1Ah == Ah)
        TrMatVecProd(A, Nm1Ah, AtNm1Ah, 0); //  AtNm1Ah = At Nm1Ah
        
        coeff=0.0;
        localreduce=0.0;
        for(i=0; i<n; i++)
            localreduce+= h[i]*AtNm1Ah[i] * pixpond[i];
        MPI_Allreduce(&localreduce, &coeff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        
        ro = g2pix/coeff;
        
        
        for(j=0; j<n; j++)            //  x = x + ro*h
            x[j] = x[j] + ro*h[j];        //
        
        
        for(j=0; j<n; j++)                  //   g = g + ro * (At Nm1 A) h
            g[j] = gp[j] - ro * AtNm1Ah[j] ; // Use Ribière-Polak formula
        
        
        if (precond == 0) {
            MatVecProd(&BJ, g, Cg, 0);
        }
        // 2lvl-specific code
        else {
            build_Qtx(A, Nm1, Z, Em1, g, pixpond, buf, Zn, n, 0);
            mul_ZQtx(Z, buf, Qg, Zn, n, 0);
            mul_ZQtx(AZ, buf, AQg, Zn, n, 0);
                        
            for (i = 0; i < n; ++i) {
                AQg[i] = g[i] - AQg[i];
            }
            
            MatVecProd(&BJ, AQg, Cg, 0);
            
            for (i = 0; i < n; ++i) {
                Cg[i] += Qg[i];
            }
        }
        
        
        g2pixp=g2pix;                               // g2p = "res"
        localreduce=0.0;
        for(i=0; i<n; i++)                    // g2 = (Cg , g)
            localreduce += Cg[i] * g[i] *pixpond[i];
        
        MPI_Allreduce(&localreduce, &g2pix, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        
        
        localreduce = 0.0;
        for(i=0; i<n; i++)                    // g2 = (Cg , g)
            localreduce += Cg[i] * (g[i]-gp[i]) *pixpond[i];
        
        MPI_Allreduce(&localreduce, &g2pix_polak, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        
        t=MPI_Wtime();                       //stop timer
        solve_time += (t-st);
        
        //Just to check with the true norm:
        if (TRUE_NORM==1) {
            localreduce=0.0;
            for(i=0; i<n; i++)                    // g2 = (Cg , g)
                localreduce += g[i] * g[i] *pixpond[i];
            
            MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
        else {
            
            res = g2pix_polak;
            
        }//end if
        
        
        if(rank==0){//print iterate info
            printf("res = %e, res0 = %e, g2pix = %e, g2pixB = %e\n", res, res0, g2pix_polak, g2pixB);
            *res_rel = sqrt(res)/sqrt(res0);
            printf("k=%d res_g2pix=%e res_g2pix_rel=%e res_rel=%e time=%lf \n", k, g2pix_polak, sqrt(g2pix_polak)/sqrt(g2pixB), sqrt(res)/sqrt(res0), t-st);
            fwrite(res_rel, sizeof(double), 1, fp);
        }
        //   if(g2pix<tol2rel){                         //
        fflush(stdout);
                
        if(res <= tol2rel){
            if(rank ==0) {                     //
                printf("--> converged (%e < %e) \n", res, tol2rel);
                printf("--> i.e. \t (%e < %e) \n", sqrt(res/res0), tol);
                printf("--> solve time = %lf \n", solve_time);
                fclose(fp);
            }
            break;
        }
        if(g2pix_polak>g2pixp){                         //
            if(rank ==0)                      //
                printf("--> g2pix>g2pixp pb (%e > %e) \n", g2pix, g2pixp);
            //    break;                        //
        }
        
        st=MPI_Wtime();
        
        //gamma = g2pix/g2pixp ;
        gamma = g2pix_polak/g2pixp;
        
        for(j=0; j<n; j++)            // h = h * gamma + Cg
            h[j] = h[j] * gamma + Cg[j] ;
        
        
    }  //End loop
    
    
    if(k==K){                //check unconverged
        if(rank==0){
            printf("--> unconverged, max iterate reached (%le > %le)\n", g2pix, tol2rel);
            fclose(fp);
        }
    }
    
    if(rank ==0)
        printf("--> res_g2pix=%e  \n", g2pix);
    
    if(rank == 0) {
        printf("MatVec time: %f\n", mv_time);
        printf("stbmmProd time: %f\n", stbmm_time);
        printf("TrMatVec time (without comms): %f\n", trmv_time);
        printf("greedy reduce time (comms only): %f\n", greedy_time);
    }
    
    free(h);
    free(Ah);
    free(g);
    free(AtNm1Ah);
    free(res_rel);
    
    return 0;
}

