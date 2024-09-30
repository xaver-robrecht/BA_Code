/**
 * @file dispersion.c
 * @author Xaver Robrecht (xaver.robrecht@gmail.com)
 * @brief
 * @version 1.0
 */

#include "mathroutines.h"
#include "spinwaves.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define BISECT_OFFSET 1e-11
#define QSTART_LOGLOG 0.0011
typedef struct node_struct {
    double q, w;
    struct node_struct *prev, *next;
} list_node;

typedef struct {
    double q, rs;
    double (*chi)(double q, double w, double rs, void* params);
    double (*f)(double q, double w, double rs, void* params);
    void* params;
} disp_prms;


double upper_stoner_bound(double q, double delta, double kf) {
    return delta + q * q / 2 + kf * fabs(q);
}


double lower_stoner_bound(double q, double delta, double kf) {
    return delta + q * q / 2 - kf * fabs(q);
}


static double bisection_wrapper(double w, void* params) {
    //workaround to pass non double parameters to a realfcnt
    disp_prms* this_prms = ((disp_prms*)params);
    double chi, f;
    chi = this_prms->chi(this_prms->q, w, this_prms->rs, this_prms->params);
    f = this_prms->f(this_prms->q, w, this_prms->rs, this_prms->params);
    double res = 1.0 - chi * f;
    printf("q = %.3e, bisection est for w = %e: %e, chi = %e, f = %e\r", this_prms->q, w, res, chi, f);
    fflush(stdout);
    return res;
}


spinwave_dispersion* dispersion(double q_step, double qmax, double rs, double delta,
                                double (*chi)(double q, double w, double rs, void* params),
                                double (*f)(double q, double w, double rs, void* params), void* params) {
    double kf = // (9/2*PI)^(1/3) * 1/rs evaluated with wolfram alpha
        2.4179879310247044610754502699589274084169211084857923335384381130 / rs;
    //struct to pass additional paramteres into bisection
    disp_prms p = {
        .chi = chi,
        .f = f,
        .rs = rs,
        .q = 0,
        .params = params
    };

    // set up linked lists for dispersion values
    list_node* vals;
    vals = malloc(sizeof(list_node));
    vals->next = NULL;
    vals->prev = NULL;
    vals->w = 0;
    vals->q = 0;


    int n_points = 1; //equal to length of vals
    int cnt=0; // counts "empty steps" where w(q) is too small
    while (1) {
        // append next list element and make it the current list element
        vals->next = malloc(sizeof(list_node));
        vals->next->next = NULL;
        vals->next->prev = vals;
        vals = vals->next;

        //calculate next q and boundary values for bisection search of root
        vals->q = p.q + q_step;
        p.q = vals->q;


        if (lower_stoner_bound(vals->q, delta, kf) <= 0) {
            perror("stoner minimum is lower than 0!");
            exit(1);
        }
        double w_upper = lower_stoner_bound(vals->q, delta, kf) - BISECT_OFFSET;
        //if(p.q < 0.1){w_upper *= 0.5;}
        double upper_sign = sign(bisection_wrapper(w_upper, &p));
        double lower_sign = sign(bisection_wrapper(BISECT_OFFSET, &p));

        //bisection search of root (calculate omega)
        if (upper_sign != lower_sign) {
            vals->w = bisect(&bisection_wrapper, BISECT_OFFSET, w_upper, &p);
            printf("\r");
            fflush(stdout);
            n_points++;
            cnt = 100;
            if (p.q > qmax && (qmax > 0)) { break; }


        }else if(n_points<=10 && cnt < 100 ){// dispersion starts very small-> keep steppin'
            vals = vals->prev;
            free(vals->next);}
        else {// finish dispersion calculation, remove empty list node
            vals = vals->prev;
            free(vals->next);
            break;
        }
    }
    printf("\r");

    /* write values into arrays and return structure that holds all values necessary to plot the dispersion */
    spinwave_dispersion* disp = malloc(sizeof(spinwave_dispersion));
    disp->n = n_points;
    disp->omega = calloc(n_points, sizeof(double));
    disp->q = calloc(n_points, sizeof(double));
    for (n_points--; n_points >= 0; n_points--) {
        disp->q[n_points] = vals->q;
        disp->omega[n_points] = vals->w;

        if (vals->prev != NULL) {
            vals = vals->prev;
            free(vals->next);
        }
        else {
            free(vals);
        }
    }
    return disp;
}

void free_dispersion(spinwave_dispersion* disp) {
    free(disp->omega);
    free(disp->q);
    free(disp);
}


void print_reference_spectrum(double rs, double qstep, double qmax) {
    double delta = 0.7696694631182531658018655526735302207774377924859591819415734686 / rs;
    double kf = // (9/2*PI)^(1/3) * 1/rs evaluated with wolfram alpha
        2.4179879310247044610754502699589274084169211084857923335384381130 / rs;
    //open file for ALDA spectrum
    char fname[20];
    sprintf(fname, "./data/rs_%.3f_alda_ref.txt", rs);
    FILE* file = fopen(fname, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }


    //calculate ALDA spectrum
    spinwave_dispersion* d1 = dispersion(qstep, qmax, rs, delta, &chi_0_ref, &f_alda_ref,NULL);
    fprintf(file, "#Reference spinwave dispersion for rs = %.3f in ALDA approximation\n", rs);
    fprintf(file, "#printed values are q, omega in atomic units\n");
    for (int i = 0; i < d1->n; i++) {
        fprintf(file, "%e,%e\n", d1->q[i], d1->omega[i]);
    }
    fclose(file);
    free_dispersion(d1);


    //calculate & write perturbation dispersion relation
    sprintf(fname, "./data/rs_%.3f_pert_ref.txt", rs);
    file = fopen(fname, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    d1 = dispersion(qstep, qmax, rs, delta, &chi_0_ref, &f_pert_ref,NULL);
    fprintf(file, "#Reference spinwave dispersion for rs = %.3f, using perturbation theory exchange kernel\n", rs);
    fprintf(file, "#printed values are q, omega in atomic units\n");
    for (int i = 0; i < d1->n; i++) {
        fprintf(file, "%e,%e\n", d1->q[i], d1->omega[i]);
    }
    fclose(file);
    free_dispersion(d1);

    //write stoner continnum
    sprintf(fname, "./data/rs_%.3e_stoner_ref.txt", rs);
    file = fopen(fname, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    fprintf(file, "#stoner continuum for rs = %.3f and lambda = %.4e using perturbation theory\n", rs);
    fprintf(file, "#printed values are q,stoner_min, stoner_max\n");
    for (double q = 0; q <= 1; q += qstep) {
        fprintf(file, "%e,%e,%e\n", q, lower_stoner_bound(q, delta, kf), upper_stoner_bound(q, delta, kf));
    }
    fclose(file);
}

void print_spectrum(double rs, double qstep, double qmax, double lambda) {
    double kf = // (9/2*PI)^(1/3) * 1/rs evaluated with wolfram alpha
        2.4179879310247044610754502699589274084169211084857923335384381130 / rs;
    const double PI = 3.14159265358979323846264338327950288419716939937510582097;
    double delta = kf / PI
        + lambda * lambda / (4 * PI * kf) * log((4 * kf * kf) / (lambda * lambda) + 1)
        - lambda / PI * atan(2 * kf / lambda);

    //calculate & write ALDA dispersion relation
    char fname[20];
    sprintf(fname, "./data/rs_%.3f_lambda_%.3e_alda.txt", rs, lambda);
    FILE* file = fopen(fname, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }


    spinwave_dispersion* d1 = dispersion(qstep, qmax, rs, delta, &chi_0, &f_alda, (double[]){delta, lambda});
    fprintf(file, "#spinwave dispersion for rs = %.3f and lambda = %.4e in ALDA approximation\n", rs, lambda);
    fprintf(file, "#printed values are q, omega\n");
    for (int i = 0; i < d1->n; i++) {
        fprintf(file, "%e,%e\n", d1->q[i], d1->omega[i]);
    }
    fclose(file);
    free_dispersion(d1);


    //calculate & write perturbation dispersion relation
    sprintf(fname, "./data/rs_%.3f_lambda_%.3e_pert.txt", rs, lambda);
    file = fopen(fname, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    d1 = dispersion(qstep, qmax, rs, delta, &chi_0, &f_pert, (double[]){delta, lambda});
    fprintf(file, "#spinwave dispersion for rs = %.3f and lambda = %.4e using perturbation theory\n", rs, lambda);
    fprintf(file, "#printed values are q, omega\n");
    for (int i = 0; i < d1->n; i++) {
        fprintf(file, "%e,%e\n", d1->q[i], d1->omega[i]);
    }
    fclose(file);
    free_dispersion(d1);

    //write stoner continnum
    sprintf(fname, "./data/rs_%.3f_lambda_%.3e_stoner.txt", rs, lambda);
    file = fopen(fname, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    fprintf(file, "#stoner continuum for rs = %.3f and lambda = %.4e using perturbation theory\n", rs, lambda);
    fprintf(file, "#printed values are q,stoner_min, stoner_max\n");
    for (double q = 0; q <= 1; q += qstep) {
        fprintf(file, "%e,%e,%e\n", q, lower_stoner_bound(q, delta, kf), upper_stoner_bound(q, delta, kf));
    }
    fclose(file);

    //calculate & write perturbation dispersion relation, with q dependance only
    sprintf(fname, "./data/rs_%.3f_lambda_%.3e_pert_q.txt", rs, lambda);
    file = fopen(fname, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    d1 = dispersion(qstep, qmax, rs, delta, &chi_0, &f_pert_q, (double[]){delta, lambda});
    fprintf(file, "#spinwave dispersion for rs = %.3f and lambda = %.4e using perturbation theory, q dependance only\n", rs, lambda);
    fprintf(file, "#printed values are q, omega\n");
    for (int i = 0; i < d1->n; i++) {
        fprintf(file, "%e,%e\n", d1->q[i], d1->omega[i]);
    }
    fclose(file);
    free_dispersion(d1);

    //calculate & write perturbation dispersion relation, with w dependance only
    sprintf(fname, "./data/rs_%.3f_lambda_%.3e_pert_w.txt", rs, lambda);
    file = fopen(fname, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    d1 = dispersion(qstep, qmax, rs, delta, &chi_0, &f_pert_w, (double[]){delta, lambda});
    fprintf(file, "#spinwave dispersion for rs = %.3f and lambda = %.4e using perturbation theory, w dependance only\n", rs, lambda);
    fprintf(file, "#printed values are q, omega\n");
    for (int i = 0; i < d1->n; i++) {
        fprintf(file, "%e,%e\n", d1->q[i], d1->omega[i]);
    }
    fclose(file);
    free_dispersion(d1);
}