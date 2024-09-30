/**
 * @file mathroutines.c
 * @author Xaver Robrecht (xaver.robrecht@gmail.com)
 * @brief 
 * @version 1.0
 */
#include "mathroutines.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>

const static double EPS = 1e-6;
const static double BISECT_EPS = 1e-6;


//wrapper for later use
double integrate(real_function func, void* prms, double a, double b) {
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(500);
    double result, error;

    gsl_function F;
    F.function = func;
    F.params = prms;

    gsl_integration_qags(&F, a, b, 0, EPS, 500, workspace, &result, &error);
    gsl_integration_workspace_free(workspace);
    return result;
}

double sign(double x) {
    if (x < 0) { return -1.0; }
    return 1.0;
}

double bisect(real_function f, double a, double b, void* params) {
    // accuracy and max iterations
    const int maxIterations = 1000;
    const double eps = BISECT_EPS;

    double limits[3] = {a,NAN, b};
    double values[3] = {f(a, params),NAN, f(b, params)};


    if (a >= b) {
        fprintf(stderr, "error in bisect:\n");
        fprintf(stderr, "\t -illegal interval limits in bisect (a >= b)!\n");
        fprintf(stderr, "f(a)=%e, f(b) =%e\n", values[0], values[2]);
        fprintf(stderr, "a=%e, b =%e\n", limits[0], limits[2]);
        exit(1);
    }
    else if (sign(values[0]) == sign(values[2])) {
        fprintf(stderr, "error in bisect:\n");
        fprintf(stderr, "\t-bisecion interval does not contain sign change!\n");
        fprintf(stderr, "f(a)=%e, f(b) =%e\n", values[0], values[2]);
        fprintf(stderr, "a=%e, b =%e\n", limits[0], limits[2]);
        exit(1);
    }

    for (int i = 0; i < maxIterations; i++) {
        limits[1] = (limits[0] + limits[2]) / 2.0;
        values[1] = f(limits[1], params);

        //check accuracy
        if ((limits[2] - limits[0]) / 2.0 < eps) { return limits[1]; }

        //bisection step
        if (sign(values[0]) != sign(values[1])) {
            limits[2] = limits[1];
            values[2] = values[1];
        }
        else {
            limits[0] = limits[1];
            values[0] = values[1];
        }
    }
    perror("error in bisect:\n\t-Desired accuracy in bisect couldn't be met!\n");
    return NAN;
}
