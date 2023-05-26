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
#include <unistd.h>

#include "midapack.h"
// #include "s2hat.h"
// #include "s2hat_tools.h"
// #include "domain_generalization.h"
#include "mappraiser.h"

int apply_weights(Tpltz *Nm1, double *tod);

int apply_PNmP(Mat *A, double *h, Tpltz *Nm1, double *AtNm1Ah);
// Apply P^t(N^-1)P to the local pixel map h, with output AtNm1Ah

int apply_sys_matrix(Mat *A, Tpltz *Nm1, struct Precond *p, Harmonic_superstruct *Harmonic_sup, PCG_var *input_variable, PCG_var *output_variable);

int compute_norm(double *norm_to_compute, double *vector_left, double *vector_right, double *weighting_pond, int size);
// Compute norm of vector_left with vector_right

int compute_norm_v2(double *norm_to_compute, double *vector_left, double *vector_right, double *vector_right_2, double *weighting_pond, int size);
// Compute norm of vector_left with (vector_right-vector_right_2)

int swap_pointers(PCG_var *PCG_variable, PCG_var *PCG_variable_2);

/** Perform PCG routine **/
int PCG_GLS_true(char *outpath, char *ref, Mat *A, Tpltz *Nm1, PCG_var *PCG_variable, double *b, double *noise, double *cond, int *lhits, double tol, int K, int precond, int Z_2lvl, int is_pixel_scheme_ring, int nside, Harmonic_superstruct *Harmonic_sup, int lmax, char *c_ell_path, int number_correlations)
{
    int i, j, k; // some indexes
    int m, n;    // number of local time samples, number of local pixels
    int rank, size;
    double localreduce; // reduce buffer
    double st, t;       // timers
    double solve_time = 0.0;
    double res, res0, res_rel;

    double *_g, *ACg, *Ah, *Nm1Ah; // time domain vectors
    double *g, *gp, *gt, *Cg, *h;  // map domain vectors
    double *AtNm1Ah;               // map domain
    double ro, gamma, coeff;       // scalars
    double g2pix, g2pixp, g2pix_polak;

    Precond *p = NULL;
    double *pixpond;
    


    // if we want to use the true norm to compute the residual
    int TRUE_NORM = 1; // 0: No ; 1: Yes

    // Size of the alm tables if no harmonic transform or Wiener filtert is applied, for initialisation
    int size_alm = 0;
    // Changed if harmonic transforms used or Wiener filter applied

    FILE *fp;

    MPI_Comm_rank(A->comm, &rank);
    MPI_Comm_size(A->comm, &size);
    m = A->m;

    st = MPI_Wtime();

    if (Z_2lvl == 0)
        Z_2lvl = size;

    build_precond(&p, &pixpond, &n, A, Nm1, PCG_variable, b, noise, cond, lhits, tol, Z_2lvl, Harmonic_sup, precond);
    // building preconditonner depending on the value of precond (0 for classic, 1 and 2 for two-level precond, 3 for Wiener-filtering precond)

    t = MPI_Wtime();
    if (rank == 0)
    {
        printf("[rank %d] Preconditioner computation time = %lf \n", rank, t - st);
        // printf("[rank %d] trash_pix flag = %d \n", rank, A->trash_pix);
        // printf("[rank %d] nbr sky pixels = %d \n", rank, n);
        fflush(stdout);
    }

    // Wiener filter extension initialization ///////////////    
    if ( (PCG_variable->bool_apply_filter == 1) || (PCG_variable->domain_PCG_computation == 1)){ // Initialize both with WF or for harmonic transforms in general
        // if (rank == 0)
        // {
        //     printf("### Initializing Wiener Filter extension for harmonic operations \n");

        //     // For later: get mask with trash_pix for spherical harmonic transforms !!!! --- TO DO LATER 
        // }
        // int root=0; // Choice to put root rank to 0
        // init_s2hat_parameters_superstruct(S2HAT_params->Files_WF_struct, PCG_variable->S2HAT_parameters, A->comm);
        // Initialization of S2HAT_parameters structure

        int *mask_binary = (int *)calloc(12*nside*nside, sizeof(int));

        get_mask_from_indices(A, mask_binary, nside, 0);

        init_harmonic_superstruct(A, Harmonic_sup, mask_binary, nside, lmax, c_ell_path, number_correlations);
        free(mask_binary);
        // Initialization of S2HAT_parameters structure

        // S2HAT_GLOBAL_parameters *Global_param_s2hat = PCG_variable->S2HAT_parameters->Global_param_s2hat;
        // The S2HAT_parameters structure has been initialized, definition of the varariable corresponding to the global S2HAT parameters which will be known by all mpi-tasks

        // Prepare to allocate non-empty alm
        size_alm = Harmonic_sup->S2HAT_params.size_alm;
    }
    // End of Wiener filter initialization


    // map domain objects memory allocation
    h = (double *)malloc(n * sizeof(double));  // descent direction
    g = (double *)malloc(n * sizeof(double));  // residual
    gp = (double *)malloc(n * sizeof(double)); // residual of previous iteration
    AtNm1Ah = (double *)malloc(n * sizeof(double));

    // time domain objects memory allocation
    Ah = (double *)malloc(m * sizeof(double));

    _g = Ah;
    Cg = AtNm1Ah;
    Nm1Ah = Ah;

    // Create PCG_var struct for h, the descent direction
    // PCG_var *Descent_dir_var = (PCG_var *)malloc(1*sizeof(PCG_var)); // corresponds to g
    PCG_var Descent_dir_var_; // corresponds to g
    initialize_PCG_var_struct(&(Descent_dir_var_), h);
    PCG_var *Descent_dir_var = &(Descent_dir_var_);
    // Define Descent_dir_var = h (in pixel space), will allocate alm only if relevant, ie if PCG_variable->bool_apply_filter==1 or PCG_variable->domain_PCG_computation==1

    // Create PCG_var struct for g, the residual
    // PCG_var *Residual_var = (PCG_var *)malloc(1*sizeof(PCG_var)); // corresponds to g
    PCG_var Residual_var_; // corresponds to g
    initialize_PCG_var_struct(&(Residual_var_), g);
    PCG_var *Residual_var = &(Residual_var_);
    // Define Residual_var = g (in pixel space), will allocate alm only if relevant, ie if PCG_variable->bool_apply_filter==1 or PCG_variable->domain_PCG_computation==1

    // Create PCG_var struct for gp, the residual of previous iteration
    // PCG_var *Last_Iter_Res_var = (PCG_var *)malloc(1*sizeof(PCG_var)); // corresponds to Cg
    PCG_var Last_Iter_Res_var_; // corresponds to Cg
    initialize_PCG_var_struct(&(Last_Iter_Res_var_), gp);
    PCG_var *Last_Iter_Res_var = &(Last_Iter_Res_var_);
    // Define Last_Iter_Res_var = gp (in pixel space), will allocate alm only if relevant, ie if PCG_variable->bool_apply_filter==1 or PCG_variable->domain_PCG_computation==1


    // Create PCG_var struct for Cg, the preconditionned residual
    // PCG_var *PrecRes_var = (PCG_var *)malloc(1*sizeof(PCG_var)); // corresponds to Cg
    PCG_var PrecRes_var_; // corresponds to Cg
    initialize_PCG_var_struct(&(PrecRes_var_), Cg);
    PCG_var *PrecRes_var = &(PrecRes_var_);
    // Define PrecRes_var = Cg (in pixel space), will allocate alm only if relevant, ie if PCG_variable->bool_apply_filter==1 or PCG_variable->domain_PCG_computation==1

    st = MPI_Wtime();

    // RHS computed in pixel for now -- TO DO LATER -> To compute in harmonic as well if asked
    // Compute RHS - initial guess
    MatVecProd(A, PCG_variable->local_map_pix, _g, 0);
    
    for (i = 0; i < m; i++)
        _g[i] = b[i] + noise[i] - _g[i];
    
    apply_weights(Nm1, _g); // _g = Nm1 (d-Ax0)  (d = signal + noise)

    // TrMatVecProd(A, _g, g, 0); // g = At _g
    TrMatVecProd(A, _g, Residual_var->local_map_pix, 0); // Residual_var = g = At _g
    // Residual_var->does_local_alm_need_update = 1; // Prepare to update local_alm
    // update_PCG_var(Residual_var, A); // Update Residual_var

    // apply_precond(p, A, &Nm1, g, Cg);
    apply_precond(p, A, Nm1, Harmonic_sup, Residual_var, PrecRes_var);
    // update_PCG_var(PrecRes_var, A); // Update PrecRes_var

    // for (j = 0; j < n; j++) // h = Cg
    //     h[j] = Cg[j];
    for (j = 0; j < n; j++) // h = Cg
        Descent_dir_var->local_map_pix[j] = PrecRes_var->local_map_pix[j];

    // g2pix = 0.0; // g2pix = "Cg res"
    // localreduce = 0.0;
    // for (i = 0; i < n; i++) // g2pix = (Cg, g)
    //     localreduce += Cg[i] * g[i] * pixpond[i];
    // // Previous version above
    g2pix = 0.0; // g2pix = "Cg res"
    localreduce = 0.0;
    compute_norm(&localreduce, PrecRes_var->local_map_pix, Residual_var->local_map_pix, pixpond, n); 
    //Compute the norm in local reduce with Cg=PrecRes_var->local_map_pix ; g=Residual_var->local_map_pix 


    MPI_Allreduce(&localreduce, &g2pix, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    t = MPI_Wtime();
    solve_time += (t - st);

    if (TRUE_NORM == 1)
    {
        res = 0.0; // g2 = "res"
        // localreduce = 0.0;
        // for (i = 0; i < n; i++) // g2 = (g, g)
        //     localreduce += Residual_var->local_map_pix[i] * Residual_var->local_map_pix[i] * pixpond[i];
        
        localreduce = 0.0;
        compute_norm(&localreduce, Residual_var->local_map_pix, Residual_var->local_map_pix, pixpond, n); 

        MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    else
    {
        res = g2pix;
    }

    double g2pixB = g2pix;
    double tol2rel = tol * tol * res; // tol*tol*g2pixB; //*g2pixB; //tol; //*tol*g2;
    res0 = res;
    // Test if already converged
    if (rank == 0)
    {

        res_rel = sqrt(res) / sqrt(res0);
        printf("k = %d, res = %e, g2pix = %e, res_rel = %e, time = %lf\n", 0, res, g2pix, res_rel, t - st);
        char filename[256];
        sprintf(filename, "%s/pcg_residuals_%s.dat", outpath, ref);
        fp = fopen(filename, "wb");
        fwrite(&res_rel, sizeof(double), 1, fp);
        fflush(stdout);
    }

    if (res <= tol)
    {
        if (rank == 0)
            printf("--> converged (%e < %e)\n", res, tol);
        k = K; // to not enter inside the loop
    }

    st = MPI_Wtime();
    fflush(stdout);

    // PCG Descent Loop *********************************************
    for (k = 1; k < K; k++)
    {

        // Swap g backup pointers (Ribière-Polak needs g from previous iteration)
        // gt = gp;
        // gp = g;
        // g = gt;
        swap_pointers(Residual_var, Last_Iter_Res_var); 
        // Swap maps and alms tables of Residual_var and Last_Iter_Res_var


        MatVecProd(A, h, Ah, 0); // Ah = A h

        apply_weights(Nm1, Nm1Ah); // Nm1Ah = Nm1 Ah   (Nm1Ah == Ah)

        // All 3 steps above should be contained in the following line
        // apply_PNmP(A, h, Nm1, AtNm1Ah);

        // All the steps commented above are executed by the following line
        apply_sys_matrix(A, Nm1, p, Harmonic_sup, Descent_dir_var, PrecRes_var);

        coeff = 0.0;
        localreduce = 0.0;
        // for (i = 0; i < n; i++)
        //     localreduce += h[i] * AtNm1Ah[i] * pixpond[i];
        compute_norm(&localreduce, Descent_dir_var->local_map_pix, PrecRes_var->local_map_pix, pixpond, n);

        MPI_Allreduce(&localreduce, &coeff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        ro = g2pix / coeff;

        for (j = 0; j < n; j++) // x = x + ro * h
            PCG_variable->local_map_pix[j] = PCG_variable->local_map_pix[j] + ro * Descent_dir_var->local_map_pix[j];
        // PCG_variable->does_local_alm_need_update = 1; // Prepare to update local_alm
        // update_PCG_var(PCG_variable, A); // Update Residual_var

        for (j = 0; j < n; j++)             // g = g + ro * (At Nm1 A) h
            Residual_var->local_map_pix[j] = Last_Iter_Res_var->local_map_pix[j] - ro * PrecRes_var->local_map_pix[j]; // Use Ribière-Polak formula
        // Residual_var->does_local_alm_need_update = 1; // Prepare to update local_alm
        // update_PCG_var(Residual_var, A); // Update Residual_var

        // apply_precond(p, A, &Nm1, g, Cg);
        apply_precond(p, A, Nm1, Harmonic_sup, Residual_var, PrecRes_var);



        g2pixp = g2pix; // g2p = "res"
        localreduce = 0.0;
        // for (i = 0; i < n; i++) // g2 = (Cg, g)
        //     localreduce += PrecRes_var->local_map_pix[i] * Residual_var->local_map_pix[i] * pixpond[i];
        compute_norm(&localreduce, PrecRes_var->local_map_pix, Residual_var->local_map_pix, pixpond, n);

        MPI_Allreduce(&localreduce, &g2pix, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        localreduce = 0.0;
        // for (i = 0; i < n; i++) // g2 = (Cg, g)
        //     localreduce += PrecRes_var->local_map_pix[i] * (Residual_var->local_map_pix[i] - Last_Iter_Res_var->local_map_pix[i]) * pixpond[i];
        compute_norm_v2(&localreduce, PrecRes_var->local_map_pix, Residual_var->local_map_pix, Last_Iter_Res_var->local_map_pix, pixpond, n);

        MPI_Allreduce(&localreduce, &g2pix_polak, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        t = MPI_Wtime();
        solve_time += (t - st);

        // Just to check with the true norm:
        if (TRUE_NORM == 1)
        {
            localreduce = 0.0;
            // for (i = 0; i < n; i++) // g2 = (Cg, g)
            //     localreduce += g[i] * g[i] * pixpond[i];
            compute_norm(&localreduce, Residual_var->local_map_pix, Residual_var->local_map_pix, pixpond, n);

            MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
        else
        {
            res = g2pix_polak;
        }

        if (rank == 0)
        { // print iterate info
            res_rel = sqrt(res) / sqrt(res0);
            printf("k = %d, res = %e, g2pix = %e, res_rel = %e, time = %lf\n", k, res, g2pix_polak, res_rel, t - st);
            fwrite(&res_rel, sizeof(double), 1, fp);
        }

        fflush(stdout);

        if (res <= tol2rel)
        {
            if (rank == 0)
            {
                printf("--> converged (%e < %e) \n", res, tol2rel);
                printf("--> i.e. \t (%e < %e) \n", res_rel, tol);
                printf("--> solve time = %lf \n", solve_time);
                fclose(fp);
            }
            break;
        }

        if (g2pix_polak > g2pixp)
        {
            if (rank == 0)
                printf("--> g2pix > g2pixp pb (%e > %e) \n", g2pix, g2pixp);
        }

        st = MPI_Wtime();

        // gamma = g2pix / g2pixp;
        gamma = g2pix_polak / g2pixp;

        for (j = 0; j < n; j++) // h = h * gamma + Cg
            Descent_dir_var->local_map_pix[j] = Descent_dir_var->local_map_pix[j] * gamma + PrecRes_var->local_map_pix[j];

    } // End loop

    if (k == K)
    { // check unconverged
        if (rank == 0)
        {
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


/* Weights TOD data according to the adopted noise model*/
int apply_weights(Tpltz *Nm1, double *tod)
{
    int t_id; // time sample index in local data
    int i, j;

    // Use straightforward loop for white noise model
    if (Nm1->tpltzblocks[0].lambda == 1)
    {
        // Here it is assumed that we use a single bandwidth for all TOD intervals, i.e. lambda is the same for all Toeplitz blocks
        t_id = 0;
        for (i = 0; i < Nm1->nb_blocks_loc; i++)
        {
            for (j = 0; j < Nm1->tpltzblocks[i].n; j++)
            {
                tod[t_id + j] = Nm1->tpltzblocks[i].T_block[0] * tod[t_id + j];
            }
            t_id += Nm1->tpltzblocks[i].n;
        }
    }
    // Use stbmmProd routine for correlated noise model (No det-det correlations for now)
    else
        stbmmProd(Nm1, tod);

    return 0;
}


int apply_PNmP(Mat *A, double *h, Tpltz *Nm1, double *AtNm1Ah){
    /*  Apply P^T N^{-1} P to vector h ; 
        return the result in AtNm1Ah
        
        Note that the application is done in pixel_MAPPRAISER 
    */

    double *Ah, *Nm1Ah; // time domain vectors
    int m = A->m;
    Ah = (double *)malloc(m * sizeof(double));
    Nm1Ah = Ah;

    MatVecProd(A, h, Ah, 0); // Ah = A h

    apply_weights(Nm1, Nm1Ah); // Nm1Ah = Nm1 Ah   (Nm1Ah == Ah)
    
    TrMatVecProd(A, Nm1Ah, AtNm1Ah, 0); // AtNm1Ah = At Nm1Ah

    free(Ah);
    return 0;
}

int apply_sys_matrix(Mat *A, Tpltz *Nm1, struct Precond *p, Harmonic_superstruct *Harmonic_sup, PCG_var *input_variable, PCG_var *output_variable)
{
    // Apply system of the PCG, depending of the value of output_variable->bool_apply_filter 
    // -> bool_apply_filter == 0, classic PCG is done with only application of  P^T N^{-1} P to input_variable in pixel space
    // -> bool_apply_filter == 1, application of  P^T N^{-1} P.input_variable (in pixel space) + C^{-1}.input_variable (in harmonic space)

    int i;
    double *new_local_variable_pix;
    s2hat_dcomplex *local_alm_in, *local_alm_out;

    

    apply_PNmP(A, input_variable->local_map_pix, Nm1, output_variable->local_map_pix);

    if(output_variable->bool_apply_filter==1)
    {
        S2HAT_parameters *S2HAT_params = &(Harmonic_sup->S2HAT_params);

        local_alm_in = (s2hat_dcomplex *) malloc( S2HAT_params->size_alm * sizeof(s2hat_dcomplex));
        local_alm_out = (s2hat_dcomplex *) malloc( S2HAT_params->size_alm * sizeof(s2hat_dcomplex));

        global_map_2_harmonic(input_variable->local_map_pix,local_alm_in, A, Harmonic_sup);

        apply_inv_covariance_matrix_to_alm(local_alm_in, local_alm_out, p->inverse_covariance_matrix, S2HAT_params);
        free(local_alm_in);

        // Do the addition of (C^{-1} + P^T N^{-1} P) in pixel or harmonic space
        switch (output_variable->domain_PCG_computation)
        {
            case 0: // Addition done in pixel domain
            new_local_variable_pix = (double *)malloc(p->n*sizeof(double));
            global_harmonic_2_map(new_local_variable_pix, local_alm_out, A, Harmonic_sup);
            // Transformation of a_lm back into pixel domain

            for (i = 0; i < p->n; ++i)
            {
                output_variable->local_map_pix[i] += new_local_variable_pix[i]; // Adding C^(-1).input_variable and P N^{-1} P.input_variable
            }

            // output_variable->does_local_alm_need_update = 1; // Change done on pixel domain, harmonic domain need update
            
            free(new_local_variable_pix);


            // case 1:
            // local_alm_out = (s2hat_dcomplex *) malloc( A->nnz * S2HAT_params->size_alm * sizeof(s2hat_dcomplex));
    
            // global_map_2_harmonic(output_variable->local_map_pix, local_alm_out, A, S2HAT_params);
            // // Result in pixel transformed in harmonic

            // for (i = 0; i < A->nnz * S2HAT_params->size_alm ; ++i)
            // {
            //     output_variable->local_alm[i].re += local_alm_out[i].re; // Adding C^(-1).input_variable and P N^{-1} P.input_variable
            //     output_variable->local_alm[i].im += local_alm_out[i].im; // Adding C^(-1).input_variable and P N^{-1} P.input_variable
            // }

            // output_variable->does_map_pixel_need_update = 1; // Change done on harmonic domain, pixel domain need update

            // free(local_alm_out);
            break;
        }
    }
    return 0;
}


int compute_norm(double *norm_to_compute, double *vector_left, double *vector_right, double *weighting_pond, int size)
{    
    // Compute norm in pixel domain of vector_left with vector_right
    // BEWARE : the norm_to_compute is reinitialized to 0 !
    int index;
    *norm_to_compute = 0.0;
    for (index = 0; index < size; index++)
        *norm_to_compute += vector_left[index] * vector_right[index] * weighting_pond[index];

    return 0;
}

int compute_norm_v2(double *norm_to_compute, double *vector_left, double *vector_right, double *vector_right_2, double *weighting_pond, int size)
{    
    // Compute norm of vector_left with (vector_right-vector_right_2)
    // BEWARE : the norm_to_compute is reinitialized to 0 !
    int index;
    *norm_to_compute = 0.0;
    for (index = 0; index < size; index++)
        *norm_to_compute += vector_left[index] * (vector_right[index] - vector_right_2[index]) * weighting_pond[index];
    
    return 0;
}

int swap_pointers(PCG_var *PCG_variable, PCG_var *PCG_variable_2)
{
    // Swap tables of PCG_variable and PCG_variable_2
    double *var_exchange_pix;
    s2hat_dcomplex *var_exchange_alm;
    
    // Initial code
    // gt = gp;
    // gp = g;
    // g = gt;
    
    var_exchange_pix = PCG_variable->local_map_pix;
    // var_exchange_alm = PCG_variable->local_alm;
    
    PCG_variable->local_map_pix = PCG_variable_2->local_map_pix ;
    // PCG_variable->local_alm = PCG_variable_2->local_alm;
    
    PCG_variable_2->local_map_pix = var_exchange_pix;
    // PCG_variable_2->local_alm = var_exchange_alm;

    return 0;
}
