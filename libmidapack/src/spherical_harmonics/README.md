# Spherical Harmonics Tools for Wiener Filtering

## Purpose

This extension to MAPPRAISER add Wiener filtering in the PCG routine. 

This specific directory contain all the files necessary to wrap the routine [S2HAT](https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html) for MAPPRAISER.

## Organisation of the directory

In particular, we create two main structures, `S2HAT_GLOBAL_parameters` and `S2HAT_LOCAL_parameters` containing all the relevant information to call `S2HAT` routines.
- `S2HAT_GLOBAL_parameters` contains the information all the processors need : choice of pixelization scheme, information about sky coverage
- `S2HAT_LOCAL_parameters` contains the information needed by each processors to efficiently decompose the operations done in harmonic domain, mostly for pixel sky and $a_{\ell m}$ decomposition, including the $m$ values given to each processor, the first and last ring the processor will cover, etc.

All those details are given in the file [s2hat_init_parameters.c](s2hat_init_parameters.c) and in the documentation, starting [here](https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/Cmanual/Cprelims.html).

The file [files_io.c](files_io.c) provides the tools to read the files required by this extension.

The file [covariance_matrix_tools.c](covariance_matrix_tools.c) provides the tools to create and use the CMB covariance matrix which will be used.

Finally, the file [alm_pix_tools.c](alm_pix_tools.c) wraps [S2HAT](https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html) routines centered around transformation from pixel sky maps to $a_{\ell m}$.

## Formulation of the proposed solution

In order to solve the initial mapmaking problem with Wiener-filtering, given by 

$$ d = Ps + n$$

We add the covariance matrix $C$ to the PCG approach of the mapmaking problem, given here by :

$$ (C^{-1} + P^t N^{-1} P) a^{WF}_{\ell m} = P^t N^{-1} d$$

We also apply this solution mostly in harmonic domain, with $a^{WF}_{\ell m}$ as the variable updated by each iteration.

For now, we consider the preconditioner given by :

$$ M = (P^t diag(N^{-1}) P)^{-1} + C^{-1} $$

## Use of the harmonic operations

In order to perform the harmonic operations, we use S2HAT which an algorithm to efficiently compute spherical harmonics operations. It proceeds by using the pixel scheme as ring scheme, and dividing them on each MPI process i.e. each MPI process will have at least one ring. Then, it optimizes spherical 

As a reminder, the ring scheme proceeds by dividing the pixels in the sky into $4*NSIDE -1$, the rings being distributed symmetrically to the equatorial ring.


## Notes

- You may need to add a link to S2HAT in your `.bashrc`, in the form `S2HATROOT=path/to/s2hat`
- When applying `apply_pix2alm` and its reverse `apply_alm2pix` to a map you generated, make sure the map was generated with a definite $\ell_{max}$, and in particular, if a pixel white noise is applied to the map, make sure it is band-limited, for instance with a beam with $\ell_{max}$ you defined.
- The two arguments which will be needed by the user are : a mask of the map (if not defined, it will just be a full map of the sky full of 1); a list of $c_\ell$ which will be used as the CMB covariance assumption for the Wiener filtering. It is important to note that **the $c_\ell$ file must be constructed with the same $\ell_{max} which the user is supposed to give as a parameter, otherwise the file won't be read correctly**. It is also worth noting the user can either provide 4 correlations (TT, EE, BB, TE), or 6 (TT, EE, BB, TE, TB, EB), in this order.