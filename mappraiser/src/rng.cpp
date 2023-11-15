/**
 * @brief Implementation of routines for generating random streams
 * @author Simon Biquard
 * @date March 2023
 * @credit This code was adapted from the TOAST library
 * whose licence can be found under this banner.
 */

// Time Ordered Astrophysics Scalable Tools (TOAST)
//
// Copyright (c) 2015-2018, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "Random123/boxmuller.hpp"
#include "Random123/threefry.h"
#include "Random123/uniform.hpp"
#include <vector>

#include "mappraiser/rng.h"

typedef r123::Threefry2x64 RNG;

// Unsigned 64bit random integers
void mappraiser::rng_dist_uint64(size_t n, uint64_t key1, uint64_t key2,
                                 uint64_t counter1, uint64_t counter2,
                                 uint64_t *data) {
    RNG rng;
    RNG::ukey_type uk = {{key1, key2}};

    for (size_t i = 0; i < n; ++i) {
        data[i] = rng(RNG::ctr_type({{counter1, counter2 + i}}),
                      RNG::key_type(uk))[0];
    }
}

// Unsigned 64bit random integers x 2
void mappraiser::rng_dist_uint64_2(size_t n, uint64_t key1, uint64_t key2,
                                   uint64_t counter1, uint64_t counter2,
                                   mappraiser::uint64_2 *data) {
    RNG rng;
    RNG::ukey_type uk = {{key1, key2}};
    RNG::ctr_type r;

    for (size_t i = 0; i < n; ++i) {
        r = rng(RNG::ctr_type({{counter1, counter2 + i}}), RNG::key_type(uk));
        data[i].x = r[0];
        data[i].y = r[1];
    }
}

// Uniform double precision values on [0.0, 1.0]
void mappraiser::rng_dist_uniform_01(size_t n, uint64_t key1, uint64_t key2,
                                     uint64_t counter1, uint64_t counter2,
                                     double *data) {
    RNG rng;
    RNG::ukey_type uk = {{key1, key2}};

    for (size_t i = 0; i < n; ++i) {
        data[i] = r123::u01<double, uint64_t>(rng(
            RNG::ctr_type({{counter1, counter2 + i}}), RNG::key_type(uk))[0]);
    }
}

// Uniform double precision values on [-1.0, 1.0]
void mappraiser::rng_dist_uniform_11(size_t n, uint64_t key1, uint64_t key2,
                                     uint64_t counter1, uint64_t counter2,
                                     double *data) {
    RNG rng;
    RNG::ukey_type uk = {{key1, key2}};

    for (size_t i = 0; i < n; ++i) {
        data[i] = r123::uneg11<double, uint64_t>(rng(
            RNG::ctr_type({{counter1, counter2 + i}}), RNG::key_type(uk))[0]);
    }
}

// Normal distribution.
void mappraiser::rng_dist_normal(size_t n, uint64_t key1, uint64_t key2,
                                 uint64_t counter1, uint64_t counter2,
                                 double *data) {
    // First compute 64bit random integers
    auto tmp = std::vector<uint64_2>(n);
    mappraiser::rng_dist_uint64_2(n, key1, key2, counter1, counter2,
                                  tmp.data());

    // Now convert pairs of uniform randoms with Box-Muller transformation
    r123::double2 pair;
    for (size_t i = 0; i < n; i++) {
        pair = r123::boxmuller(tmp[i].x, tmp[i].y);
        data[i] = pair.x;
    }
}
