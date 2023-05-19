/**
 * @file solver_info.c
 * @author Simon Biquard
 * @brief Implementation of routines for the SolverInfo structure.
 * @version 0.1
 * @date Jan 2023
 */

#include <math.h>
#include <stdio.h>

#include "mappraiser/solver_info.h"

/// @brief Print current solver parameters and information
/// @param si SolverInfo struct
void solverinfo_print(SolverInfo *si) {
    puts("Solver parameters:");
    printf("  max_steps          = %d\n", si->max_steps);
    printf("  abs_res_reduct     = %e\n", si->abs_res_reduct);
    printf("  rel_res_reduct     = %e\n", si->rel_res_reduct);
    printf("  rel_res_growth     = %e\n", si->rel_res_growth);
    printf("  use_exact_residual = %s\n", si->use_exact_residual ? "true" : "false");
    printf("  store_hist         = %s\n", si->store_hist ? "true" : "false");
    printf("  print              = %s\n", si->print ? "true" : "false");

    puts("Solver information:");
    // printf("  start_time    = %lf\n", si->start_time);
    printf("  r0            = %e\n", si->r0);
    printf("  has_converged = %s\n", si->has_converged ? "true" : "false");
    printf("  has_diverged  = %s\n", si->has_diverged ? "true" : "false");
    printf("  has_failed    = %s\n", si->has_failed ? "true" : "false");
    printf("  n_iter        = %d\n", si->n_iter);
    printf("  solve_time    = %lf s\n", si->solve_time);
    printf("  conv_rate     = %lf\n", si->conv_rate);
    printf("  res_norm      = %e\n", si->res_norm);
    if (si->store_hist) {
        fputs("  res_hist      = { ", stdout);
        for (int j = 0; j < (si->n_iter) + 1; ++j) { printf("%lf ", si->res_hist[j]); }
        puts("}");
    }
}

/// @brief Set default values for the solver parameters
/// @param si SolverInfo struct
void solverinfo_set_defaults(SolverInfo *si) {
    si->max_steps          = 100;
    si->abs_res_reduct     = 1e-14;
    si->rel_res_reduct     = 1e-06;
    si->rel_res_growth     = 1e+06;
    si->use_exact_residual = false;

    si->store_hist = false;
    si->print      = false;
}

/// @brief Initialise/allocate output variables of the solver, print initial message if needed
/// @param si SolverInfo struct
void solverinfo_init(SolverInfo *si) {
    si->has_converged = false;
    si->has_diverged  = false;
    si->has_failed    = false;
    si->n_iter        = 0;
    si->conv_rate     = 0.0;
    si->res_norm      = 0.0;
    si->res_hist      = NULL;

    if (si->store_hist) si->res_hist = malloc((sizeof si->res_hist) * si->max_steps);

    if (si->print) {
        printf("PCG solver running with kmax = %d\n", si->max_steps);
        fflush(stdout);
    }
}

/// @brief Update the SolverInfo structure depending on the iteration state
/// and determine if the iteration is to be continued or stopped.
/// @param si [in/out] SolverInfo struct
/// @param stop [out] continue or stop the iteration
/// @param step_nbr [in] current iteration number
/// @param res [in] current residual norm
/// @param wtime [in] current walltime
void solverinfo_update(SolverInfo *si, bool *stop, int step_nbr, double res, double wtime) {
    // if k=0, store the residual for later
    if (step_nbr == 0) si->r0 = res;

    // print information on screen if needed
    if (si->print) {
        printf("  k = %d, |r|^2 = %e, |r|/|r0| = %e\n", step_nbr, res, sqrt(res / si->r0));
        fflush(stdout);
    }

    // store residual history if needed
    if (si->store_hist) si->res_hist[step_nbr] = res;

    // stop the iteration if the maximal number of steps is reached
    if (step_nbr >= si->max_steps) *stop = true;

    // check other stop conditions (convergence/divergence)
    if (res <= si->rel_res_reduct * si->rel_res_reduct * si->r0) {
        si->has_converged = true;
        *stop             = true;
    }

    if (res <= si->abs_res_reduct * si->abs_res_reduct) {
        si->has_converged = true;
        *stop             = true;
    }

    if (res >= si->rel_res_growth * si->rel_res_growth * si->r0) {
        si->has_diverged = true;
        *stop            = true;
    }

    // handle the case where iteration stops
    if (*stop) {
        si->res_norm   = res;
        si->n_iter     = step_nbr;
        si->solve_time = wtime - si->start_time;
    }
}

/// @brief Reallocate history buffer and print results on screen if needed
/// @param si SolverInfo struct
void solverinfo_finalize(SolverInfo *si) {
    if ((si->n_iter < si->max_steps) && (si->store_hist)) {
        // res_history has size n_iter + 1 !!!
        double *tmp;
        tmp = realloc(si->res_hist, (sizeof si->res_hist) * (si->n_iter + 1));
        if (tmp != NULL) { si->res_hist = tmp; }
    }

    if (si->print) {
        if (si->has_converged) {
            if (si->res_norm <= si->abs_res_reduct * si->abs_res_reduct)
                printf("  -> converged (|r| <= %e)\n", si->abs_res_reduct);

            if (si->res_norm <= si->rel_res_reduct * si->rel_res_reduct * si->r0)
                printf("  -> converged (|r|/|r0| <= %e)\n", si->rel_res_reduct);
        } else if (si->has_diverged) {
            printf("  -> diverged (|r|/|r0| >= %e)\n", si->rel_res_growth);
        } else if (si->has_failed) {
            puts("  -> failed (possible cause: memory issue)");
        } else {
            // no convergence nor divergence
            printf("  -> maximal iteration number reached (%d)\n", si->max_steps);
        }
        printf("  -> n_iter = %d, time = %lf s\n", si->n_iter, si->solve_time);
        fflush(stdout);
    }
}

/// @brief Free any allocated buffer inside the structure
/// @param si SolverInfo struct
void solverinfo_free(SolverInfo *si) {
    if (si->res_hist) { free(si->res_hist); }
    si->res_hist = NULL;
}
