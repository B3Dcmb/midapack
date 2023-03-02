/**
 * @file random.h
 * @author Simon Biquard, following Numerical Recipes The Art of Scientific Computing (Ch. 7)
 * @brief Declaration of the Random and NormalDev structs, routine prototypes.
 * @version 0.1
 * @date Feb 2023
 */

#include <stdlib.h>
#include <stdbool.h>

#ifndef RANDOM_H
#define RANDOM_H

#define SQR(x) ((x) * (x))

//------------------------------------------------------------*
// struct Random
//------------------------------------------------------------*

/// @brief Struct for generating uniform deviates
struct random_t
{
    // random generator internal state needs 3 values
    unsigned long long u, v, w;
};
typedef struct random_t Random;

void random_init(Random *r, unsigned long long seed);
unsigned long long random_getint64(Random *r);
double random_getdouble(Random *r);

//------------------------------------------------------------*
// struct NormalDev
//------------------------------------------------------------*

/// @brief Struct for generating normal deviates
struct normaldev_t
{
    // parameters of the normal distribution
    double mu, sig;

    // internal random generator (uniform)
    Random ran;
};
typedef struct normaldev_t NormalDev;

void normaldev_init(NormalDev *nd, double mu, double sigma, unsigned long long seed);
double normaldev_get(NormalDev *nd);
void normaldev_getN(NormalDev *nd, int nb, double *buf);

#endif /* RANDOM_H */
