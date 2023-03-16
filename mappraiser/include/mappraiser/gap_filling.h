#ifndef MAPPRAISER_GAP_FILLING_H
#define MAPPRAISER_GAP_FILLING_H

#include <midapack.h>

#ifdef __cplusplus
namespace mappraiser {
    void psd_from_tt(int fftlen, int lambda, int psdlen, const double *tt, double *psd, double rate = 200.0);

    double compute_mean(int samples, double *buf, bool subtract = false);

    double compute_variance(int samples, const double& mean, double *buf);

    void sim_noise_tod(int samples, int lambda, const double *tt, double *buf, double var_goal = 1.0);

    void sim_constrained_noise_block(Tpltz *N_block, Tpltz *Nm1_block, const double *noise,
                                     Gap *gaps, double *out_constrained);

    void sim_constrained_noise(Tpltz *N, Tpltz *Nm1, const double *noise,
                               Gap *gaps, double *out_constrained);
}
#endif

// Here the routine destined to be called by C code
#ifdef __cplusplus
extern "C" {
#endif

void sim_constrained_noise(Tpltz *N, Tpltz *Nm1, const double *noise, Gap *gaps, double *out_constrained);

#ifdef __cplusplus
}
#endif

#endif //MAPPRAISER_GAP_FILLING_H
