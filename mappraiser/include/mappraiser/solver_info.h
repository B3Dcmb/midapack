#ifndef MAPPRAISER_SOLVER_INFO_H
#define MAPPRAISER_SOLVER_INFO_H

#include <stdlib.h>
#include <stdbool.h>

struct solverinfo_t
{
    int max_steps;         // maximal number of iteration steps
    double abs_res_reduct; // absolute residual norm (tolerance) (stop if |r_i| <= epsilon)
    double rel_res_reduct; // relative residual norm (stop if |r_i| <= epsilon * |r_0|)
    double rel_res_growth; // relative residual growth (divergence) (stop if |r_i| >= epsilon * |r_0|)

    bool use_exact_residual; // compute exact residual during iteration
    bool store_hist;         // store complete iteration history
    bool print;              // print this data on screen during iteration

    double start_time;  // start time of the solver
    double r0;          // initial residual of the solver
    bool has_converged; // did the iteration converge
    bool has_diverged;  // did the iteration diverge (based on a stopping criterion)
    bool has_failed;    // did the iteration fail (possible reasons?)
    int n_iter;         // number of iteration steps
    double solve_time;  // time elapsed during iteration
    double conv_rate;   // convergence rate
    double res_norm;    // norm of residual after iteration stopped
    double *res_hist;   // complete iteration history (only used if store_hist is true)
};

typedef struct solverinfo_t SolverInfo;

void solverinfo_print(SolverInfo *si);
void solverinfo_set_defaults(SolverInfo *si);
void solverinfo_init(SolverInfo *si);
void solverinfo_update(SolverInfo *si, bool *stop, int step_nbr, double res, double wtime);
void solverinfo_finalize(SolverInfo *si);
void solverinfo_free(SolverInfo *si);

#endif //MAPPRAISER_SOLVER_INFO_H
