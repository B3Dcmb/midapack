#ifndef MAPPRAISER_RNG_H
#define MAPPRAISER_RNG_H

#ifdef __cplusplus

#include <cstddef>
#include <cstdint>

namespace mappraiser {

typedef struct {
    uint64_t x, y;
} uint64_2;

void rng_dist_uint64(size_t n, uint64_t key1, uint64_t key2, uint64_t counter1,
                     uint64_t counter2, uint64_t *data);

void rng_dist_uint64_2(size_t n, uint64_t key1, uint64_t key2,
                       uint64_t counter1, uint64_t counter2, uint64_2 *data);

void rng_dist_uniform_01(size_t n, uint64_t key1, uint64_t key2,
                         uint64_t counter1, uint64_t counter2, double *data);

void rng_dist_uniform_11(size_t n, uint64_t key1, uint64_t key2,
                         uint64_t counter1, uint64_t counter2, double *data);

void rng_dist_normal(size_t n, uint64_t key1, uint64_t key2, uint64_t counter1,
                     uint64_t counter2, double *data);
} // namespace mappraiser

#endif

#endif // MAPPRAISER_RNG_H
