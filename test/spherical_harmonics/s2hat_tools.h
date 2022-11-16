 /* Full documentation here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html */ 
    

// typedef struct S2HAT_GLOBAL_parameters S2HAT_GLOBAL_parameters;
typedef struct S2HAT_GLOBAL_parameters{
    /* Global parameters of S2HAT, to give to all processors */
    s2hat_pixeltype pixelization_scheme;
    s2hat_scandef scan_sky_structure_pixel;
    s2hat_pixparameters pixpar;

    int nside;
    int nlmax;
    int nmmax;
} S2HAT_GLOBAL_parameters;

// typedef struct S2HAT_LOCAL_parameters S2HAT_LOCAL_parameters;
typedef struct S2HAT_LOCAL_parameters{
    /* Local parameters of S2HAT, dependent on each processor */

    int gangrank;
    int gangsize;
    int gangroot;


    int nmvals;
    int first_ring;
    int last_ring;
    int map_size;
    int* mvals; // size given by nmvals


    // Tools to precompute Legendre functions, but not taken into account by s2hat if plms=0, which is the default behaviour we choose
    int plms;
    long int nplm;
} S2HAT_LOCAL_parameters;

