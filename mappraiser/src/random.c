/**
 * @file random.c
 * @author Simon Biquard, following Numerical Recipes The Art of Scientific Computing (Ch. 7)
 * @brief Implementation of routines for the Random and NormalDev structs.
 * @version 0.1
 * @date Feb 2023
 */

#include <stdlib.h>
#include <math.h>
#include "random.h"

//------------------------------------------------------------*
// struct Random routines
//------------------------------------------------------------*

/// @brief Initialise a struct to generate uniform deviates
/// @param r Random struct
/// @param seed any integer except 4101842887655102017LL
void random_init(Random *r, unsigned long long seed)
{
    r->v = 4101842887655102017LL;
    r->w = 1;
    r->u = seed ^ r->v;
    r->v = r->u;
    r->w = r->v;
}

/// @brief Draw a 64-bit random integer
/// @param r Random struct (must be initialised with random_init)
/// @return the generated number
unsigned long long random_getint64(Random *r)
{
    r->u = r->u * 2862933555777941757LL + 7046029254386353087LL;
    r->v ^= r->v >> 17;
    r->v ^= r->v << 31;
    r->v ^= r->v >> 8;
    r->w = 4294957665U * (r->w & 0xffffffff) + (r->w >> 32);
    unsigned long long x = r->u ^ (r->u << 21);
    x ^= x >> 35;
    x ^= x << 4;
    return (x + r->v) ^ r->w;
}

/// @brief Draw a random double-precision floating value in the range 0. to 1.
/// @param r Random struct (must be initialised with random_init)
/// @return the generated number
double random_getdouble(Random *r)
{
    return 5.42101086242752217E-20 * random_getint64(r);
}

//------------------------------------------------------------*
// struct NormalDev routines
//------------------------------------------------------------*

/// @brief Initialise a struct to generate normal (Gaussian) deviates
/// @param nd NormalDev struct
/// @param mu mean of the distribution
/// @param sigma standard deviation
/// @param seed any integer except 4101842887655102017LL
void normaldev_init(NormalDev *nd, double mu, double sigma, unsigned long long seed)
{
    nd->mu = mu;
    nd->sig = sigma;

    Random myran;
    random_init(&myran, seed);
    nd->ran = myran;
}

/// @brief Draw a double-precision floating normal deviate
/// @param nd [in] NormalDev struct
/// @return the generated number
double normaldev_get(NormalDev *nd)
{
    double u, v, x, y, q;
    do
    {
        u = random_getdouble(&nd->ran);
        v = 1.7156 * (random_getdouble(&nd->ran) - 0.5);
        x = u - 0.449871;
        y = fabs(v) + 0.386595;
        q = SQR(x) + y * (0.19600 * y - 0.25472 * x);
    } while (q > 0.27597 && (q > 0.27846 || SQR(v) > -4. * log(u) * SQR(u)));
    return nd->mu + nd->sig * v / u;
}

/// @brief Draw N normal deviates
/// @param nd [in] NormalDev struct
/// @param nb [in] number of values to draw
/// @param buf [out] buffer of size at least N for the generated numbers
void normaldev_getN(NormalDev *nd, int nb, double *buf)
{
    for (int i = 0; i < nb; i++)
    {
        buf[i] = normaldev_get(nd);
    }
}
