import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pymaster as nmt


def generate_full_covariance_matrix(nside, mask_apodized, c_ells_array, delta_ell, lmax=0, new=True):
    
    cl_tt, cl_ee, cl_bb, cl_te = c_ells_array[:4, :]

    if c_ells_array.shape[0] <= 4:
        cl_tb = 0*cl_tt
        cl_eb = 0*cl_tt
    else :
        cl_eb, cl_tb = c_ells_array[4:]

    map_TQU = hp.synfast(c_ells_array, nside, new=new)

    spin_0_field = nmt.NmtField(mask_apodized, [map_TQU[0]])
    spin_2_field = nmt.NmtField(mask_apodized, [map_TQU[1], map_TQU[2]])

    if lmax==0:
        lmax = 3*nside-1
        bin_scheme = nmt.NmtBin.from_nside_linear(nside, delta_ell)
    else:
        bin_scheme = nmt.NmtBin.from_lmax_linear(lmax, delta_ell)

    n_ell = lmax-1

    workspace_NaMS_00 = nmt.NmtWorkspace()
    workspace_NaMS_00.compute_coupling_matrix(spin_0_field, spin_0_field, bin_scheme)
    workspace_NaMS_02 = nmt.NmtWorkspace()
    workspace_NaMS_02.compute_coupling_matrix(spin_0_field, spin_2_field, bin_scheme)
    workspace_NaMS_22 = nmt.NmtWorkspace()
    workspace_NaMS_22.compute_coupling_matrix(spin_2_field, spin_2_field, bin_scheme)


    covar_workspace = nmt.NmtCovarianceWorkspace()
    # This is the time-consuming operation
    # Note that you only need to do this once,
    # regardless of spin
    covar_workspace.compute_coupling_coefficients(spin_0_field, spin_0_field, spin_0_field, spin_0_field)

    covar_dict = dict() # Allo covariances

    # The next few lines show how to extract the covariance matrices
    # for different spin combinations.
    covar_00_00 = nmt.gaussian_covariance(covar_workspace,
                                        0, 0, 0, 0,  # Spins of the 4 fields
                                        [cl_tt],  # TT
                                        [cl_tt],  # TT
                                        [cl_tt],  # TT
                                        [cl_tt],  # TT
                                        workspace_NaMS_00, wb=workspace_NaMS_00).reshape([n_ell, 1,
                                                                n_ell, 1])
    covar_dict['TT_TT'] = covar_00_00[:, 0, :, 0]
    # covar_TT_TT = covar_00_00[:, 0, :, 0]

    covar_02_02 = nmt.gaussian_covariance(covar_workspace, 0, 2, 0, 2,  # Spins of the 4 fields
                                        [cl_tt],  # TT
                                        [cl_te, cl_tb],  # TE, TB
                                        [cl_te, cl_tb],  # ET, BT
                                        [cl_ee, cl_eb,
                                        cl_eb, cl_bb],  # EE, EB, BE, BB
                                        workspace_NaMS_02, wb=workspace_NaMS_02).reshape([n_ell, 2,
                                                                n_ell, 2])
    list_correl = ['TE', 'TB']
    for i in range(2):
        for j in range(2):
            covar_dict[list_correl[i] + '_' + list_correl[j]] = covar_02_02[:, i, :, j]

    # covar_TE_TE = covar_02_02[:, 0, :, 0]
    # covar_TE_TB = covar_02_02[:, 0, :, 1]
    # covar_TB_TE = covar_02_02[:, 1, :, 0]
    # covar_TB_TB = covar_02_02[:, 1, :, 1]


    covar_00_22 = nmt.gaussian_covariance(covar_workspace, 0, 0, 2, 2,  # Spins of the 4 fields
                                        [cl_te, cl_tb],  # TE, TB
                                        [cl_te, cl_tb],  # TE, TB
                                        [cl_te, cl_tb],  # TE, TB
                                        [cl_te, cl_tb],  # TE, TB
                                        workspace_NaMS_00, wb=workspace_NaMS_22).reshape([n_ell, 1,
                                                                n_ell, 4])
    list_correl = ['EE', 'EB', 'BE', 'BB']
    for j in range(4):
        covar_dict['TT' + '_' + list_correl[j]] = covar_00_22[:, 0, :, j]
    
    # covar_TT_EE = covar_00_22[:, 0, :, 0]
    # covar_TT_EB = covar_00_22[:, 0, :, 1]
    # covar_TT_BE = covar_00_22[:, 0, :, 2]
    # covar_TT_BB = covar_00_22[:, 0, :, 3]

    covar_02_22 = nmt.gaussian_covariance(covar_workspace, 0, 2, 2, 2,  # Spins of the 4 fields
                                        [cl_te, cl_tb],  # TE, TB
                                        [cl_te, cl_tb],  # TE, TB
                                        [cl_ee, cl_eb,
                                        cl_eb, cl_bb],  # EE, EB, BE, BB
                                        [cl_ee, cl_eb,
                                        cl_eb, cl_bb],  # EE, EB, BE, BB
                                        workspace_NaMS_02, wb=workspace_NaMS_22).reshape([n_ell, 2,
                                                                n_ell, 4])
    list_correl = ['EE', 'EB', 'BE', 'BB']
    for j in range(4):
        covar_dict['TE' + '_' + list_correl[j]] = covar_02_22[:, 0, :, j]
        covar_dict['TB' + '_' + list_correl[j]] = covar_02_22[:, 1, :, j]
    
    # covar_TE_EE = covar_02_22[:, 0, :, 0]
    # covar_TE_EB = covar_02_22[:, 0, :, 1]
    # covar_TE_BE = covar_02_22[:, 0, :, 2]
    # covar_TE_BB = covar_02_22[:, 0, :, 3]
    # covar_TB_EE = covar_02_22[:, 1, :, 0]
    # covar_TB_EB = covar_02_22[:, 1, :, 1]
    # covar_TB_BE = covar_02_22[:, 1, :, 2]
    # covar_TB_BB = covar_02_22[:, 1, :, 3]


    covar_22_22 = nmt.gaussian_covariance(covar_workspace, 2, 2, 2, 2,  # Spins of the 4 fields
                                        [cl_ee, cl_eb,
                                        cl_eb, cl_bb],  # EE, EB, BE, BB
                                        [cl_ee, cl_eb,
                                        cl_eb, cl_bb],  # EE, EB, BE, BB
                                        [cl_ee, cl_eb,
                                        cl_eb, cl_bb],  # EE, EB, BE, BB
                                        [cl_ee, cl_eb,
                                        cl_eb, cl_bb],  # EE, EB, BE, BB
                                        workspace_NaMS_22, wb=workspace_NaMS_22).reshape([n_ell, 4,
                                                                n_ell, 4])
    list_correl = ['EE', 'EB', 'BE', 'BB']
    for i in range(4):
        for j in range(4):
            covar_dict[list_correl[i] + '_' + list_correl[j]] = covar_22_22[:, i, :, j]

    # covar_EE_EE = covar_22_22[:, 0, :, 0]
    # covar_EE_EB = covar_22_22[:, 0, :, 1]
    # covar_EE_BE = covar_22_22[:, 0, :, 2]
    # covar_EE_BB = covar_22_22[:, 0, :, 3]
    # covar_EB_EE = covar_22_22[:, 1, :, 0]
    # covar_EB_EB = covar_22_22[:, 1, :, 1]
    # covar_EB_BE = covar_22_22[:, 1, :, 2]
    # covar_EB_BB = covar_22_22[:, 1, :, 3]
    # covar_BE_EE = covar_22_22[:, 2, :, 0]
    # covar_BE_EB = covar_22_22[:, 2, :, 1]
    # covar_BE_BE = covar_22_22[:, 2, :, 2]
    # covar_BE_BB = covar_22_22[:, 2, :, 3]
    # covar_BB_EE = covar_22_22[:, 3, :, 0]
    # covar_BB_EB = covar_22_22[:, 3, :, 1]
    # covar_BB_BE = covar_22_22[:, 3, :, 2]
    # covar_BB_BB = covar_22_22[:, 3, :, 3]

    return covar_dict

