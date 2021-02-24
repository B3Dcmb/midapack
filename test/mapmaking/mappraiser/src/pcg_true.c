// Midapack library
// mapmaking code example using the Midapack library - release 1.2b, Nov 2012
// PCG routine applied to the map-making equation
// This can use the block-diagonal jacobi or Two-level preconditionners

/** @file   pcg_true.c
 @author Hamza El Bouhargani
 @date   May 2019
 @credit  Adapted from work by Frederic Dauvergne
 @Last_update June 2020 by Aygul Jamal */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <mkl.h>
#include <unistd.h>
#include "fitsio.h"
#include "midapack.h"
#include "mappraiser.h"

int xmap( double *mapI, double *mapQ, double *mapU, int npix, double *x, int *lstid, int xsize);

int PCG_GLS_true(char *outpath, char *ref, Mat *A, Tpltz Nm1, double *x, double *b, double *noise, double *cond, int *lhits, double tol, int K, int precond, int Z_2lvl)
{
    int    i, j, k;     // some indexes
    int    m, n;        // number of local time samples, number of local pixels
    int    rank, size;
    double localreduce; // reduce buffer
    double st, t;       // timers
    double solve_time = 0.0;
    double res, res0, res_rel;

    double *_g, *ACg, *Ah, *Nm1Ah; // time domain vectors
    double *g, *gp, *gt, *Cg, *h;  // map domain vectors
    double *AtNm1Ah;               // map domain
    double ro, gamma, coeff;       // scalars
    double g2pix, g2pixp, g2pix_polak;

    struct Precond *p = NULL;
    double *pixpond;
    MPI_Status status;
    // if we want to use the true norm to compute the residual
    int TRUE_NORM = 1;  //0: No ; 1: Yes

    FILE *fp;

    MPI_Comm_rank(A->comm, &rank);
    MPI_Comm_size(A->comm, &size);
    m = A->m;

    st = MPI_Wtime();

    if (Z_2lvl == 0) Z_2lvl = size;
    build_precond(&p, &pixpond, &n, A, &Nm1, &x, b, noise, cond, lhits, tol, Z_2lvl, precond);

    t = MPI_Wtime();
    if (rank == 0) {
        printf("[rank %d] Preconditioner computation time = %lf \n", rank, t - st);
        fflush(stdout);
    }

    // map domain objects memory allocation
    h = (double *) malloc(n * sizeof(double)); // descent direction
    g = (double *) malloc(n * sizeof(double)); // residual
    gp = (double *) malloc(n * sizeof(double)); // residual of previous iteration
    AtNm1Ah = (double *) malloc(n * sizeof(double));

    // time domain objects memory allocation
    Ah = (double *) malloc(m * sizeof(double));

    _g = Ah;
    Cg = AtNm1Ah;
    Nm1Ah = Ah;

    st = MPI_Wtime();

    // Compute RHS
    MatVecProd(A, x, _g, 0);

    for (i = 0; i < m; i++)
        _g[i] = b[i] + noise[i] - _g[i];

    stbmmProd(Nm1, _g); // _g = Nm1 (Ax-b)

    TrMatVecProd(A, _g, g, 0); // g = At _g

    apply_precond(p, A, &Nm1, g, Cg);

    for (j = 0; j < n; j++) // h = -Cg
        h[j] = Cg[j];

    g2pix = 0.0; // g2 = "res"
    localreduce = 0.0;
    for (i = 0; i < n; i++) // g2 = (Cg, g)
        localreduce += Cg[i] * g[i] * pixpond[i];

    MPI_Allreduce(&localreduce, &g2pix, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    t = MPI_Wtime();
    solve_time += (t - st);

    // Just to check with the true norm:
    if (TRUE_NORM == 1) {
        res = 0.0; // g2 = "res"
        localreduce = 0.0;
        for (i = 0; i < n; i++) // g2 = (Cg, g)
            localreduce += g[i] * g[i] * pixpond[i];

        MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    else {
        res = g2pix;
    }

    double g2pixB = g2pix;
    double tol2rel = tol * tol * res; //tol*tol*g2pixB; //*g2pixB; //tol; //*tol*g2;
    res0 = res;
    //Declarations for the sake of saving iterations maps
    int mapsize, map_id, Nnz, nside, npix, oldsize;
    int *lstid;
    double *xbis, *mapI, *mapQ, *mapU;
    char Imap_name[256];
    char Qmap_name[256];
    char Umap_name[256];
    char nest = 1;
    char *cordsys = "C";
    int ret,w=1;
    // Test if already converged
    if (rank == 0) {

        res_rel = sqrt(res) / sqrt(res0);
	      printf("k = %d, res = %e, g2pix = %e, res_rel = %e, time = %lf\n", 0, res, g2pix, res_rel, t - st);
        char filename[256];
        sprintf(filename,"%s/pcg_residuals_%s.dat", outpath, ref);
        fp = fopen(filename, "wb");
        fwrite(&res_rel, sizeof(double), 1, fp);
        fflush(stdout);
    }

    //write output to fits files:

    st=MPI_Wtime();
    mapsize = n;
    map_id = rank;
    Nnz = A->nnz;

    if (rank == 0)
      xbis = (double *) calloc(mapsize, sizeof(double));

    lstid = (int *) calloc(mapsize, sizeof(int));
    for(i=0; i< mapsize; i++){
      lstid[i] = A->lindices[i+(A->nnz)*(A->trash_pix)];
    }


    if (rank!=0){
      MPI_Send(&mapsize, 1, MPI_INT, 0, 0, A->comm);
      MPI_Send(lstid, mapsize, MPI_INT, 0, 1, A->comm);
      MPI_Send(x, mapsize, MPI_DOUBLE, 0, 2, A->comm);
    }

    if (rank==0){
      nside = 512;
      npix = 12*pow(nside,2);
      oldsize;

      mapI    = (double *) calloc(npix, sizeof(double));
      mapQ    = (double *) calloc(npix, sizeof(double));
      mapU    = (double *) calloc(npix, sizeof(double));

      for (i=0;i<size;i++){
        if (i!=0){
          oldsize = mapsize;
          MPI_Recv(&mapsize, 1, MPI_INT, i, 0, A->comm, &status);
          if (oldsize!=mapsize){
            lstid = (int *) realloc(lstid, mapsize*sizeof(int));
            xbis = (double *) realloc(xbis, mapsize*sizeof(double));
          }
          MPI_Recv(lstid, mapsize, MPI_INT, i, 1, A->comm, &status);
          MPI_Recv(xbis, mapsize, MPI_DOUBLE, i, 2, A->comm, &status);
        }
        else
          memcpy(xbis, x, mapsize*sizeof(double));
        xmap(mapI, mapQ, mapU, npix, xbis, lstid, mapsize);
      }
      //printf("Checking output directory ... old files will be overwritten\n");

      sprintf(Imap_name,"%s/mapI_%s_iter%d.fits", outpath, ref, 0);
      sprintf(Qmap_name,"%s/mapQ_%s_iter%d.fits", outpath, ref, 0);
      sprintf(Umap_name,"%s/mapU_%s_iter%d.fits", outpath, ref, 0);

      if( access( Imap_name, F_OK ) != -1 ) {
        ret = remove(Imap_name);
        if(ret != 0){
          printf("Error: unable to delete the file %s\n",Imap_name);
          w = 0;
        }
      }

      if( access( Qmap_name, F_OK ) != -1 ) {
        ret = remove(Qmap_name);
        if(ret != 0){
          printf("Error: unable to delete the file %s\n",Qmap_name);
          w = 0;
        }
      }

      if( access( Umap_name, F_OK ) != -1 ) {
        ret = remove(Umap_name);
        if(ret != 0){
          printf("Error: unable to delete the file %s\n",Umap_name);
          w = 0;
        }
      }


      if(w==1){
        //printf("Writing HEALPix maps FITS files ...\n");
        write_map(mapI, TDOUBLE, nside, Imap_name, nest, cordsys);
        write_map(mapQ, TDOUBLE, nside, Qmap_name, nest, cordsys);
        write_map(mapU, TDOUBLE, nside, Umap_name, nest, cordsys);
      }
      else{
        printf("IO Error: Could not overwrite old files, map results will not be stored ;(\n");
      }

    }

    t=MPI_Wtime();
    if (rank==0) {
      printf("[rank %d] Write output files time=%lf \n", rank, t-st);
      fflush(stdout);
    }

    if (res <= tol) {
        if (rank == 0)
            printf("--> converged (%e < %e)\n", res, tol);
        k = K; // to not enter inside the loop
    }

    st = MPI_Wtime();
    fflush(stdout);


    // PCG Descent Loop *********************************************
    for (k = 1; k < K ; k++){

        // Swap g backup pointers (Ribière-Polak needs g from previous iteration)
        gt = gp;
        gp = g;
        g = gt;

        MatVecProd(A, h, Ah, 0);            // Ah = A h
        stbmmProd(Nm1, Nm1Ah);              // Nm1Ah = Nm1 Ah   (Nm1Ah == Ah)
        TrMatVecProd(A, Nm1Ah, AtNm1Ah, 0); // AtNm1Ah = At Nm1Ah

        coeff = 0.0;
        localreduce = 0.0;
        for (i = 0; i < n; i++)
            localreduce += h[i] * AtNm1Ah[i] * pixpond[i];
        MPI_Allreduce(&localreduce, &coeff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        ro = g2pix / coeff;


        for (j = 0; j < n; j++) // x = x + ro * h
            x[j] = x[j] + ro * h[j];

        for (j = 0; j < n; j++)             // g = g + ro * (At Nm1 A) h
            g[j] = gp[j] - ro * AtNm1Ah[j]; // Use Ribière-Polak formula

	      apply_precond(p, A, &Nm1, g, Cg);

        g2pixp = g2pix; // g2p = "res"
        localreduce = 0.0;
        for (i = 0; i < n; i++) // g2 = (Cg, g)
            localreduce += Cg[i] * g[i] * pixpond[i];

        MPI_Allreduce(&localreduce, &g2pix, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        localreduce = 0.0;
        for (i = 0; i < n; i++) // g2 = (Cg, g)
            localreduce += Cg[i] * (g[i] - gp[i]) * pixpond[i];

        MPI_Allreduce(&localreduce, &g2pix_polak, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        t = MPI_Wtime();
        solve_time += (t - st);

        // Just to check with the true norm:
        if (TRUE_NORM == 1) {
            localreduce = 0.0;
            for (i = 0; i < n; i++) // g2 = (Cg, g)
                localreduce += g[i] * g[i] * pixpond[i];

            MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
        else {
            res = g2pix_polak;
        }

        if (rank == 0){ //print iterate info
	          res_rel = sqrt(res) / sqrt(res0);
            printf("k = %d, res = %e, g2pix = %e, res_rel = %e, time = %lf\n", k, res, g2pix_polak, res_rel, t - st);
            fwrite(&res_rel, sizeof(double), 1, fp);
        }
        //write output to fits files:
        st=MPI_Wtime();

        if (rank!=0){
          MPI_Send(&mapsize, 1, MPI_INT, 0, 0, A->comm);
          MPI_Send(lstid, mapsize, MPI_INT, 0, 1, A->comm);
          MPI_Send(x, mapsize, MPI_DOUBLE, 0, 2, A->comm);
        }

        if (rank==0){

          for (i=0;i<size;i++){
            if (i!=0){
              oldsize = mapsize;
              MPI_Recv(&mapsize, 1, MPI_INT, i, 0, A->comm, &status);
              if (oldsize!=mapsize){
                lstid = (int *) realloc(lstid, mapsize*sizeof(int));
                xbis = (double *) realloc(xbis, mapsize*sizeof(double));
              }
              MPI_Recv(lstid, mapsize, MPI_INT, i, 1, A->comm, &status);
              MPI_Recv(xbis, mapsize, MPI_DOUBLE, i, 2, A->comm, &status);
            }
            else
              memcpy(xbis, x, mapsize*sizeof(double));
            xmap(mapI, mapQ, mapU, npix, xbis, lstid, mapsize);
          }
          //printf("Checking output directory ... old files will be overwritten\n");

          sprintf(Imap_name,"%s/mapI_%s_iter%d.fits", outpath, ref, k);
          sprintf(Qmap_name,"%s/mapQ_%s_iter%d.fits", outpath, ref, k);
          sprintf(Umap_name,"%s/mapU_%s_iter%d.fits", outpath, ref, k);

          if( access( Imap_name, F_OK ) != -1 ) {
            ret = remove(Imap_name);
            if(ret != 0){
              printf("Error: unable to delete the file %s\n",Imap_name);
              w = 0;
            }
          }

          if( access( Qmap_name, F_OK ) != -1 ) {
            ret = remove(Qmap_name);
            if(ret != 0){
              printf("Error: unable to delete the file %s\n",Qmap_name);
              w = 0;
            }
          }

          if( access( Umap_name, F_OK ) != -1 ) {
            ret = remove(Umap_name);
            if(ret != 0){
              printf("Error: unable to delete the file %s\n",Umap_name);
              w = 0;
            }
          }


          if(w==1){
            //printf("Writing HEALPix maps FITS files ...\n");
            write_map(mapI, TDOUBLE, nside, Imap_name, nest, cordsys);
            write_map(mapQ, TDOUBLE, nside, Qmap_name, nest, cordsys);
            write_map(mapU, TDOUBLE, nside, Umap_name, nest, cordsys);
          }
          else{
            printf("IO Error: Could not overwrite old files, map results will not be stored ;(\n");
          }

        }

        t=MPI_Wtime();
        if (rank==0) {
          printf("[rank %d] Write output files time=%lf \n", rank, t-st);
          fflush(stdout);
        }
        fflush(stdout);

        if (res <= tol2rel) {
            if (rank == 0) {
                printf("--> converged (%e < %e) \n", res, tol2rel);
                printf("--> i.e. \t (%e < %e) \n", res_rel, tol);
                printf("--> solve time = %lf \n", solve_time);
                fclose(fp);
                //to remove when done
                free(lstid);
                free(xbis);
                free(mapI);
                free(mapQ);
                free(mapU);

            }
            break;
        }

        if (g2pix_polak > g2pixp) {
            if (rank == 0)
                printf("--> g2pix > g2pixp pb (%e > %e) \n", g2pix, g2pixp);
        }

        st = MPI_Wtime();

        //gamma = g2pix / g2pixp;
        gamma = g2pix_polak / g2pixp;

        for (j = 0; j < n; j++) // h = h * gamma + Cg
            h[j] = h[j] * gamma + Cg[j];

    } // End loop

    if (k == K) { // check unconverged
        if (rank == 0) {
            printf("--> unconverged, max iterate reached (%le > %le)\n", res, tol2rel);
            fclose(fp);
        }
    }

    if (rank == 0)
        printf("--> g2pix = %e\n", g2pix);

    free(h);
    free(g);
    free(gp);
    free(AtNm1Ah);
    free(Ah);
    free_precond(&p);

    return 0;
}

int xmap( double *mapI, double *mapQ, double *mapU, int npix, double *x, int *lstid, int xsize)
{

  int i;

  for(i=0; i<xsize; i++){
    if(i%3 == 0){
      mapI[(int)(lstid[i]/3)]= x[i];
    }
    else if (i%3 == 1)
      mapQ[(int)(lstid[i]/3)]= x[i];
    else
      mapU[(int)(lstid[i]/3)]= x[i];
  }

  return 0;
}
