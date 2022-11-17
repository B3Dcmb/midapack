# Spherical Harmonics Tools for Wiener Filtering

## Purpose

This extension to MAPPRAISER add Wiener filtering in the PCG routine. 

This specific directory contain all the files necessary to wrap the routine [S2HAT](https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html) for MAPPRAISER.

## Organisation of the directory

In particular, we create two main structures, `S2HAT_GLOBAL_parameters` and `S2HAT_LOCAL_parameters` containing all the relevant information to call `S2HAT` routines.
- `S2HAT_GLOBAL_parameters` contains the information all the processors need : choice of pixelization scheme, information about sky coverage
- `S2HAT_LOCAL_parameters` contains the information needed by each processors to efficiently decompose the operations done in harmonic domain, mostly for pixel sky and $a_{\ell m}$ decomposition, including the $m$ values given to each processor, the first and last ring the processor will cover, etc.

All those details are given in the file [s2hat_init_parameters.c](s2hat_init_parameters.c) and in the documentation, starting [here](https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/Cmanual/Cprelims.html).

The file [file_tools.c](file_tools.c) provides the tools to read the files required by this extension.

The file [covariance_matrix_tools.c](covariance_matrix_tools.c) provides the tools to create and use the CMB covariance matrix which will be used.

Finally, the file [alm_pix_tools.c](alm_pix_tools.c) wraps [S2HAT](https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html) routines centered around transformation from pixel sky maps to $a_{\ell m}$.

## Formulation of the proposed solution

In order to solve the initial mapmaking problem with Wiener-filtering, given by 

$$ d = Ps + n$$

We add the covariance matrix $C$ to the PCG approach of the mapmaking problem, given here by :

$$ (C^{-1} + P^t N^{-1} P) a^{WF}_{\ell m} = P^t N^{-1} d$$

We also apply this solution mostly in harmonic domain, with with $a^{WF}_{\ell m}$ as the variable updated by each iteration.

For now, we consider the preconditioner given by :

$$ M = (P^t diag(N^{-1}) P)^{-1} + C^{-1} $$
