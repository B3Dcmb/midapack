void test_cls(int nside, int lmax, int nstokes, double *CMB_map, double *CMB_map_output, double *c_ell_output, double *red_matrix, int iter_alm, float error_alm, double *mask_binary, MPI_Comm worldcomm)
{
    // MPI_Init( &argc, &argv);
    int rank, nprocs;
    MPI_Comm_rank( worldcomm, &rank);
    MPI_Comm_size( worldcomm, &nprocs);
    // printf("Test -1 !!! \n"); fflush(stdout);
    int npix = 12*nside*nside;
    int index, index_2;

    char *c_ell_path = "";
    int number_correlations = ceil((nstokes*nstokes)/2) + floor(nstokes/2) + nstokes%2;
    if (nstokes == 1)
        number_correlations = 1;
    // printf("Test -1a !!! nb_correl %d \n", number_correlations); fflush(stdout);
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);
    // printf("Test -1ab !!! \n"); fflush(stdout);
    // printf("Initializing S2HAT_param \n"); fflush(stdout);
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary, nstokes, iter_alm, error_alm, &S2HAT_params, worldcomm);
    // printf("Finish initializing S2HAT_param !! \n"); fflush(stdout);
    // printf("Test -1bb !!! \n"); fflush(stdout);

    S2HAT_LOCAL_parameters Local_param_s2hat = S2HAT_params.Local_param_s2hat;
    // printf("Test -1b !!! \n"); fflush(stdout);
    if (Local_param_s2hat.gangrank >= 0){
        // printf("Entering Wiener-filter with rank %d \n", rank); fflush(stdout);
        
        double *CMB_map_temp  = (double *)malloc(3*npix*sizeof(double));
        // double CMB_map_temp[3*npix];
        memcpy(CMB_map_temp, CMB_map, 3*npix*sizeof(double));
        // printf("Test 0 !!! \n"); fflush(stdout);

        // double *local_map_pix = (double *) calloc( 3*Local_param_s2hat.map_size, sizeof(double));
        // MPI_Barrier(Local_param_s2hat.gangcomm);
        double local_map_pix[3*Local_param_s2hat.map_size];
        distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_temp, local_map_pix, &S2HAT_params);
        // printf("Test 1a !!! \n"); fflush(stdout);
        free(CMB_map_temp);
        // printf("Test 1b !!! \n"); fflush(stdout);
        
        // s2hat_dcomplex *local_alm, *local_alm_out;
        // local_alm = (s2hat_dcomplex *) calloc( (3*S2HAT_params.size_alm),sizeof(s2hat_dcomplex));
        // local_alm_out = (s2hat_dcomplex *) malloc( (3*S2HAT_params.size_alm)*sizeof(s2hat_dcomplex));
        // s2hat_dcomplex *local_alm_inverted = (s2hat_dcomplex *) calloc( (3*S2HAT_params.size_alm),sizeof(s2hat_dcomplex));
        s2hat_dcomplex local_alm[3*S2HAT_params.size_alm];
        s2hat_dcomplex local_alm_modified[3*S2HAT_params.size_alm];
        
        for (index=0; index<3*S2HAT_params.size_alm; index++){
            local_alm[index].re = 0;
            local_alm[index].im = 0;
            local_alm_modified[index].re = 0;
            local_alm_modified[index].im = 0;
        }
        
        
        // Doing Pix2Alm
        if (iter_alm <= 0)
            apply_pix2alm(local_map_pix, local_alm, &S2HAT_params);
        else{
            apply_pix2alm_iter(local_map_pix, local_alm, &S2HAT_params);
        }
        // free(local_map_pix);
        // printf("Test 1d !!! \n"); fflush(stdout);

        // Getting covariance matrix
        int ell_value;
        double **red_matrix_copy;
        red_matrix_copy = malloc((lmax+1)*sizeof(double *));

        // printf("Test 1e !!! \n"); fflush(stdout);
        for(ell_value=0; ell_value<lmax+1; ell_value++){
            red_matrix_copy[ell_value] = calloc(nstokes*nstokes,sizeof(double));
            // printf("- %d -", ell_value); fflush(stdout);
            memcpy(red_matrix_copy[ell_value], red_matrix+ell_value*nstokes*nstokes, nstokes*nstokes*sizeof(double));   
        }
        
        // printf("Red matrix copy \n"); fflush(stdout);
        // int number_ells_to_probe = 4;
        // int first_ell_to_probe = 0;
        // int nb_stokes_1, nb_stokes_2;
        // for(ell_value=first_ell_to_probe; ell_value<number_ells_to_probe+first_ell_to_probe; ell_value++){
        //     printf("ell : %d \n", ell_value); fflush(stdout);
        //     for(nb_stokes_1=0; nb_stokes_1<nstokes; nb_stokes_1++){
        //         for(nb_stokes_2=0; nb_stokes_2<nstokes; nb_stokes_2++){
        //             printf("- %f -", red_matrix_copy[ell_value][nb_stokes_1*nstokes+nb_stokes_2]);
        //     }
        //     printf("\n");
        //     }
        // }
        // printf("\n");
        // printf("Start applying red matrix \n"); fflush(stdout);
        apply_inv_block_diag_covariance_matrix_to_alm(local_alm, local_alm_modified, red_matrix_copy, &S2HAT_params);
        // printf("End applying red matrix \n"); fflush(stdout);

        free_covariance_matrix(red_matrix_copy, lmax);
        // printf("End freeing red matrix \n"); fflush(stdout);

        int nspec = nstokes*(nstokes+1)/2;
        double *c_ell_array = (double *)calloc((lmax+1)*nspec,sizeof(double));
        alm2cls(local_alm_modified, c_ell_array, nspec, &S2HAT_params);
        if (rank==0)
            memcpy(c_ell_output, c_ell_array, (lmax+1)*nspec*sizeof(double));

        // double *local_map_pix_out = (double *)calloc(3*Local_param_s2hat.map_size, sizeof(double));
        double local_map_pix_out[3*Local_param_s2hat.map_size];
        
        // apply_alm2pix(local_map_pix_out, local_alm, &S2HAT_params);
        apply_alm2pix(local_map_pix_out, local_alm_modified, &S2HAT_params);
        // free(local_alm);
        // free(local_alm_out);

        double *full_sky_map_2;
        full_sky_map_2 = (double *) malloc(nstokes*npix*sizeof(double));
        gather_map(local_map_pix_out, full_sky_map_2, nstokes, &S2HAT_params);
        // free(local_map_pix_out);
        // MPI_Barrier(Local_param_s2hat.gangcomm);

        if (rank==0)
            memcpy(CMB_map_output, full_sky_map_2, nstokes*npix*sizeof(double));


        free(full_sky_map_2);
        // free(local_alm_modified);

        free_s2hat_parameters_struct(&S2HAT_params);

    // MPI_Finalize();
    }
    // printf("Wiener-filter pixel done for rank %d ! \n", rank);
    fflush(stdout);
}

void test_cls(int nside, int lmax, int nstokes, double *CMB_map, double *CMB_map_output, double *c_ell_output, double *red_matrix, int iter_alm, float error_alm, double *mask_binary, MPI_Comm worldcomm)
{
    // MPI_Init( &argc, &argv);
    int rank, nprocs;
    MPI_Comm_rank( worldcomm, &rank);
    MPI_Comm_size( worldcomm, &nprocs);
    // printf("Test -1 !!! \n"); fflush(stdout);
    int npix = 12*nside*nside;
    int index, index_2;

    char *c_ell_path = "";
    int number_correlations = ceil((nstokes*nstokes)/2) + floor(nstokes/2) + nstokes%2;
    if (nstokes == 1)
        number_correlations = 1;
    // printf("Test -1a !!! nb_correl %d \n", number_correlations); fflush(stdout);
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);
    // printf("Test -1ab !!! \n"); fflush(stdout);
    // printf("Initializing S2HAT_param \n"); fflush(stdout);
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary, nstokes, iter_alm, error_alm, &S2HAT_params, worldcomm);
    // printf("Finish initializing S2HAT_param !! \n"); fflush(stdout);
    // printf("Test -1bb !!! \n"); fflush(stdout);

    S2HAT_LOCAL_parameters Local_param_s2hat = S2HAT_params.Local_param_s2hat;
    // printf("Test -1b !!! \n"); fflush(stdout);
    if (Local_param_s2hat.gangrank >= 0){
        // printf("Entering Wiener-filter with rank %d \n", rank); fflush(stdout);
        
        double *CMB_map_temp  = (double *)malloc(3*npix*sizeof(double));
        // double CMB_map_temp[3*npix];
        memcpy(CMB_map_temp, CMB_map, 3*npix*sizeof(double));
        // printf("Test 0 !!! \n"); fflush(stdout);

        // double *local_map_pix = (double *) calloc( 3*Local_param_s2hat.map_size, sizeof(double));
        // MPI_Barrier(Local_param_s2hat.gangcomm);
        double local_map_pix[3*Local_param_s2hat.map_size];
        distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_temp, local_map_pix, &S2HAT_params);
        // printf("Test 1a !!! \n"); fflush(stdout);
        free(CMB_map_temp);
        // printf("Test 1b !!! \n"); fflush(stdout);
        
        // s2hat_dcomplex *local_alm, *local_alm_out;
        // local_alm = (s2hat_dcomplex *) calloc( (3*S2HAT_params.size_alm),sizeof(s2hat_dcomplex));
        // local_alm_out = (s2hat_dcomplex *) malloc( (3*S2HAT_params.size_alm)*sizeof(s2hat_dcomplex));
        // s2hat_dcomplex *local_alm_inverted = (s2hat_dcomplex *) calloc( (3*S2HAT_params.size_alm),sizeof(s2hat_dcomplex));
        s2hat_dcomplex local_alm[3*S2HAT_params.size_alm];
        s2hat_dcomplex local_alm_modified[3*S2HAT_params.size_alm];
        
        for (index=0; index<3*S2HAT_params.size_alm; index++){
            local_alm[index].re = 0;
            local_alm[index].im = 0;
            local_alm_modified[index].re = 0;
            local_alm_modified[index].im = 0;
        }
        
        
        // Doing Pix2Alm
        if (iter_alm <= 0)
            apply_pix2alm(local_map_pix, local_alm, &S2HAT_params);
        else{
            apply_pix2alm_iter(local_map_pix, local_alm, &S2HAT_params);
        }
        // free(local_map_pix);
        // printf("Test 1d !!! \n"); fflush(stdout);

        // Getting covariance matrix
        int ell_value;
        double **red_matrix_copy;
        red_matrix_copy = malloc((lmax+1)*sizeof(double *));

        // printf("Test 1e !!! \n"); fflush(stdout);
        for(ell_value=0; ell_value<lmax+1; ell_value++){
            red_matrix_copy[ell_value] = calloc(nstokes*nstokes,sizeof(double));
            // printf("- %d -", ell_value); fflush(stdout);
            memcpy(red_matrix_copy[ell_value], red_matrix+ell_value*nstokes*nstokes, nstokes*nstokes*sizeof(double));   
        }
        
        // printf("Red matrix copy \n"); fflush(stdout);
        // int number_ells_to_probe = 4;
        // int first_ell_to_probe = 0;
        // int nb_stokes_1, nb_stokes_2;
        // for(ell_value=first_ell_to_probe; ell_value<number_ells_to_probe+first_ell_to_probe; ell_value++){
        //     printf("ell : %d \n", ell_value); fflush(stdout);
        //     for(nb_stokes_1=0; nb_stokes_1<nstokes; nb_stokes_1++){
        //         for(nb_stokes_2=0; nb_stokes_2<nstokes; nb_stokes_2++){
        //             printf("- %f -", red_matrix_copy[ell_value][nb_stokes_1*nstokes+nb_stokes_2]);
        //     }
        //     printf("\n");
        //     }
        // }
        // printf("\n");
        // printf("Start applying red matrix \n"); fflush(stdout);
        apply_inv_block_diag_covariance_matrix_to_alm(local_alm, local_alm_modified, red_matrix_copy, &S2HAT_params);
        // printf("End applying red matrix \n"); fflush(stdout);

        free_covariance_matrix(red_matrix_copy, lmax);
        // printf("End freeing red matrix \n"); fflush(stdout);

        int nspec = nstokes*(nstokes+1)/2;
        double *c_ell_array = (double *)calloc((lmax+1)*nspec,sizeof(double));
        alm2cls(local_alm_modified, c_ell_array, nspec, &S2HAT_params);
        if (rank==0)
            memcpy(c_ell_output, c_ell_array, (lmax+1)*nspec*sizeof(double));

        // double *local_map_pix_out = (double *)calloc(3*Local_param_s2hat.map_size, sizeof(double));
        double local_map_pix_out[3*Local_param_s2hat.map_size];
        
        // apply_alm2pix(local_map_pix_out, local_alm, &S2HAT_params);
        apply_alm2pix(local_map_pix_out, local_alm_modified, &S2HAT_params);
        // free(local_alm);
        // free(local_alm_out);

        double *full_sky_map_2;
        full_sky_map_2 = (double *) malloc(nstokes*npix*sizeof(double));
        gather_map(local_map_pix_out, full_sky_map_2, nstokes, &S2HAT_params);
        // free(local_map_pix_out);
        // MPI_Barrier(Local_param_s2hat.gangcomm);

        if (rank==0)
            memcpy(CMB_map_output, full_sky_map_2, nstokes*npix*sizeof(double));


        free(full_sky_map_2);
        // free(local_alm_modified);

        free_s2hat_parameters_struct(&S2HAT_params);

    // MPI_Finalize();
    }
    // printf("Wiener-filter pixel done for rank %d ! \n", rank);
    fflush(stdout);
}
