/**
* @file main.c
 * @author Xaver Robrecht (xaver.robrecht@gmail.com)
 * @brief implements routines to generate spin wave dispersion relations
 * @version 1.0
 */

#include "spinwaves.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char* argv[]) {
    double rs, qstep, lambda, qmax;
    rs = atof(argv[1]);
    qstep = atof(argv[2]);
    qmax = atof(argv[3]);
    lambda = atof(argv[4]);

    print_reference_spectrum(rs, qstep, qmax);
    print_spectrum(rs, qstep, qmax, lambda);
    return 0;
}
