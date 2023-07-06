#ifndef MAPPRAISER_GAP_FILLING_H
#define MAPPRAISER_GAP_FILLING_H

#include <midapack.h>

#ifdef __cplusplus

#include <iostream>
#include <vector>

namespace mappraiser {

    class GapFillInfo {
    public:
        GapFillInfo(int, int, int);

        double get_mean_iterations() const;
        double get_mean_seconds() const;

        void print_curr_block() const;

        void set_current_block(int i) { current_block = i; }
        void set_current_size(int n) { current_size = n; }

        void store_nb_iterations(int n) { nb_iterations[current_block] = n; }
        void store_block_time(double milliseconds) { block_times[current_block] = milliseconds * 0.001; }
        void store_pcg_time(double milliseconds) { pcg_times[current_block] = milliseconds * 0.001; }
        void store_valid_frac(double f) { valid_fracs[current_block] = f; }

        int n_gaps;        // number of local timestream gaps
        int n_blocks;      // number of data blocks to treat
        int current_block; // current block index
        int current_size;  // size of current block

        std::vector<int>    nb_iterations; // number of PCG iterations for each block
        std::vector<double> block_times;   // total time for each block
        std::vector<double> pcg_times;     // PCG time for each block
        std::vector<double> valid_fracs;   // proportion of valid samples in each block

        // id for printing
        int id;
    };

    class GapFillRecap {
    public:
        GapFillRecap() = default;
        explicit GapFillRecap(mappraiser::GapFillInfo const &info);

        void send(MPI_Comm comm, int dest) const;

        void receive(MPI_Comm comm, int src);

    private:
        void print(std::ostream &out) const;

        int    n_blocks;
        double mean_iterations;
        int    min_iter;
        int    min_iter_idx;
        int    max_iter;
        int    max_iter_idx;
        double mean_time;
        double min_time;
        double min_time_idx;
        double max_time;
        double max_time_idx;

        friend std::ostream &operator<<(std::ostream &out, GapFillRecap const &recap);
    };

    int find_valid_samples(Gap *gaps, size_t id0, std::vector<uint8_t> &valid);

    void sim_constrained_noise_block(GapFillInfo &gfi, Tpltz *N_block, Tpltz *Nm1_block, double *noise, Gap *gaps,
                                     uint64_t realization, uint64_t detindx, uint64_t obsindx, uint64_t telescope,
                                     double sample_rate);

    void sim_constrained_noise(GapFillInfo &gfi, Tpltz *N, Tpltz *Nm1, double *noise, Gap *gaps, uint64_t realization,
                               const uint64_t *detindxs, const uint64_t *obsindxs, const uint64_t *telescopes,
                               double sample_rate);
} // namespace mappraiser
#endif

// Here the prototypes of the C routines
#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

void psd_from_tt(int fftlen, int lambda, int psdlen, const double *tt, double *psd);

void remove_baseline(int samples, double *buf, double *baseline, const uint8_t *valid, int bandwidth, bool rm);

void sim_noise_tod(int samples, int lambda, const double *tt, double *buf, uint64_t realization, uint64_t detindx,
                   uint64_t obsindx, uint64_t telescope, double sample_rate);

void perform_gap_filling(MPI_Comm comm, Tpltz *N, Tpltz *Nm1, double *noise, Gap *gaps, uint64_t realization,
                         const uint64_t *detindxs, const uint64_t *obsindxs, const uint64_t *telescopes,
                         double sample_rate, bool verbose);

// test routine to call from python
void gap_filling(MPI_Comm comm, const int *data_size_proc, int nb_blocks_loc, int *local_blocks_sizes, int nnz,
                 double *tt, double *inv_tt, int lambda, double *noise, int *indices, uint64_t realization,
                 const uint64_t *detindxs, const uint64_t *obsindxs, const uint64_t *telescopes, double sample_rate);

#ifdef __cplusplus
}
#endif

#endif // MAPPRAISER_GAP_FILLING_H
