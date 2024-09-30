/**
 * @file mathroutines.h
 * @author Xaver Robrecht (xaver.robrecht@gmail.com)
 * @brief implements integration routines, multilog functions and bisection
 * @version 1.0
 */

#ifndef MATHROUTINES_H
#define MATHROUTINES_H

// function prototype used for integration (integration.c)
typedef double (*real_function)(double x, void *params);

// signum function
double sign(double x);

/**
 * @brief integrates func(x,prms) on the interval [a,b]
 * using romberg integration based on the midpoint rule(@see numerical recipes in c).
 *
 * @param func points to function of the type @see realfctn
 * @param prms used to pass additional parameters to func
 * @param a lower integration interval bound
 * @param b upper integration interval bound
 * @return integral of func in  [a,b]
 */
double integrate(real_function func, void *prms, double a, double b);


/**
 * @brief bisection search of zero of f(x) in the interval [a,b].
 *
 * Returns x such that f(x) = 0.
 * If the actual zero of f(x) lies at x0, the algorithm assures |x-x0| < eps.
 * @param f function that is evaluated
 * @param a lower interval boundary
 * @param b upper interval boundary
 * @param params additional parameters passed to f
 * @return double zero of f(x) in [a,b]
 */
double bisect(real_function f, double a, double b, void *params);

/**
 * @brief calculates the real valued dilogarithm
 *
 * @param x x in [-infty, 1]
 * @return double \f$ Li_2(x) \f$
 */
double li2(double x);

/**
 * @brief calculates the real valued trilogarithm
 *
 * @param x x in [-infty, +infty]
 * @return double \f$ Li_3(x) \f$
 */
double li3(double x);

#endif// MATHROUTINES_H