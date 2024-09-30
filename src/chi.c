/**
* @file chi.c
 * @author Xaver Robrecht (xaver.robrecht@gmail.com)
 * @brief implements functions to calculate susceptibilities
 * @version 1.0
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mathroutines.h"
#include "spinwaves.h"


// Riemann-Zeta function values
const static double zeta2 = //=pi^2/6
    1.6449340668482264364724151666460251892189499012067984377355582293;

const static double zeta3 = //zeta(3)
    1.2020569031595942853997381615114499907649862923404988817922715553;

// Other Constants
const static double pi =
    3.1415926535897932384626433832795028841971693993751058209749445923;

double chi_1a(double q, double w, double rs, void* params);
double chi_1a_z(double z, void* params);
double chi_1a_zp(double zp, void* params);

double chi_1b(double q, double w, double rs, void* params);
double chi_1b_k(double k, void* params);
double chi_1b_kd(double kd, void* params);

double chi_1c(double q, double w, double rs, void* params);


double chi_0(double q, double w, double rs, void* params) {
    double kf = // (9/2*pi)^(1/3) * 1/rs evaluated with wolfram alpha
        2.4179879310247044610754502699589274084169211084857923335384381130 / rs;
    double delta = ((double*)params)[0];
    if(q==0){return kf*kf*kf/(6*delta*pi*pi-6*pi*pi*w);}
    double a = (delta + q * q * 0.5 - w) / (q * kf);
    double result = 0.5 * (a * a - 1) * log((1.0 + a) / (a - 1.0)) - a;
    result *= kf * kf / (4.0 * pi * pi * q);
    return result;
}

double chi_1(double q, double w, double rs, void* params) {
    double chi1a = chi_1a(q, w, rs, params);
    double chi1b = chi_1b(q, w, rs, params);
    double chi1c = chi_1c(q, w, rs, params);
    return (chi1a + chi1b) - chi1c;
}


double chi_1a_z(double z, void* params) {
    double I1, w, delta, q, zp, lam, kf;
    w = ((double*)params)[0];
    q = ((double*)params)[1];
    delta = ((double*)params)[2];
    lam = ((double*)params)[3];
    kf = ((double*)params)[4];
    zp = ((double*)params)[5];
    I1 = 0.0;
    double a, c, b;
    c = sqrt(fmax(0.0, kf * kf - z * z));
    b = sqrt(fmax(0.0, kf * kf - zp * zp));
    a = lam * lam + (z - zp) * (z - zp);
    if (c == 0.0 || b == 0.0) {
        return 0.0;
    }

    double I_11 = sqrt((b - c) * (b - c) + a) * sqrt((b + c) * (b + c) + a) / 4;
    double I_12 = ((-2 * b * b - 2 * c * c) * log(a)) / 4;
    double I_13 = ((-2 * b * b - 2 * c * c) * log(2)) / 4;
    double I_14 = ((2 * log(sqrt((b - c) * (b - c) + a) * sqrt((b + c) * (b + c) + a) - b * b + c * c + a) - 1) * b * b)
        / 4;
    double I_15 = ((2 * log(sqrt((b - c) * (b - c) + a) * sqrt((b + c) * (b + c) + a) + b * b - c * c + a) - 1) * c * c)
        / 4;

    // integral of I__16 gets calculated analytically
    I1 = I_11 + I_12 + I_13 + I_14 + I_15;
    return I1 / (w - delta - q * q * 0.5 - q * z);
}

double chi_1a_zp(double zp, void* params) {
    double w, q, lam, delta, kf;
    w = ((double*)params)[0];
    q = ((double*)params)[1];
    delta = ((double*)params)[2];
    lam = ((double*)params)[3];
    kf = ((double*)params)[4];
    double sing_pos = ((w - delta - q * q * 0.5) / q);
    double result = integrate(&chi_1a_z, (double[]){w, q, delta, lam, kf, zp}, -kf, kf);
    return result / (w - delta - q * q * 0.5 - q * zp);
}

double chi_1a(double q, double w, double rs, void* params) {
    double delta = ((double*)params)[0];
    double lam = ((double*)params)[1];
    double kf = 2.4179879310247044610754502699589274084169211084857923335384381130 / rs;

    //sum of I__11 ... I__15
    double result = integrate(&chi_1a_zp, (double[]){w, q, delta, lam, kf}, -kf, kf);

    //I__16
    double d = w - delta - q * q * 0.5;
    double I__16 = (d - kf * q) / (d + kf * q); //argument of logs in I_16
    I__16 = -lam * lam / (4 * q * q) * (log(I__16) * log(I__16))
            + d * kf / (q * q * q) * log(I__16) + 2 * kf * kf / (q * q);
    if(q==0){result += (-2*kf*kf*kf*kf-3*lam*lam*kf*kf)/(3*d*d);}
    else{result += I__16;}
    double factor = 1.0 / (2.0 * pi);
    factor = factor * factor * factor;
    return -factor * result;
}


double chi_1b_k(double k, void* params) {
    double kf, lam, w, q, delta;
    w = ((double*)params)[0];
    q = ((double*)params)[1];
    delta = ((double*)params)[2];
    lam = ((double*)params)[3];
    kf = ((double*)params)[4];

    //g1+g2 = sigma*k
    double g1 = pi * (k * k - kf * kf - lam * lam) * log(
        ((k - kf) * (k - kf) + lam * lam) / ((k + kf) * (k + kf) + lam * lam)
        );
    double g2 = 4 * pi * k * (lam * (atan((k - kf) / lam) - atan((k + kf) / lam)) + kf);
    double d = (w - q * q / 2 - delta);
    double num = (-1 / (4 * pi * pi) * (g1 + g2) * k);
    double denum = (d * d - k * q * k * q);
    return num / denum;
}

double chi_1b(double q, double w, double rs, void* params) {
    double kf, lam, delta;
    delta = ((double*)params)[0];
    lam = ((double*)params)[1];
    kf = // (9/2*pi)^(1/3) * 1/rs evaluated with wolfram alpha
        2.4179879310247044610754502699589274084169211084857923335384381130 / rs;
    double result;
    result = integrate(&chi_1b_k, (double[]){w, q, delta, lam, kf}, 0, kf);
    return -1.0 / (2 * pi * pi) * result;
}

double chi_1c(double q, double w, double rs, void* params) {
    double kf = // (9/2*pi)^(1/3) * 1/rs evaluated with wolfram alpha
        2.4179879310247044610754502699589274084169211084857923335384381130 / rs;
    double delta = ((double*)params)[0];
    double a = (delta + q * q * 0.5 - w) / (q * kf);
    double result = a * 0.5 * log((a + 1.0) / (a - 1.0)) - 1.0;
    result *= kf * delta / (2 * pi * pi * q * q);
    if(q==0){return kf*kf*kf*delta/(6*delta*delta*pi*pi-12*delta*pi*pi*w+6*pi*pi*w*w);}
    return result;
}


double chi_0_ref(double q, double w, double rs, void* params) {
    double delta = // (9/2*pi)^(1/3) * 1/(rs*pi) evaluated with wolfram alpha
        2.4179879310247044610754502699589274084169211084857923335384381130 / (rs * pi);
    double result = chi_0(q, w, rs, &delta);
    return result;
}

double chi_1_ref(double q, double w, double rs, void* params) {
    double kf_up, a, delta;
    kf_up = // (9/2*pi)^(1/3) * 1/rs evaluated with wolfram alpha
        2.4179879310247044610754502699589274084169211084857923335384381130 / rs;
    delta = kf_up / pi;
    a = (delta + q * q * 0.5 - w) / (q * kf_up);
    double result =
    (
        -2
        + (a + 1) * (a + 1) * li2(2.0 / (1.0 + a))
        + (a - 1) * (a - 1) * li2(2.0 / (1.0 - a))
        - 2 * (a * a - 1) *log(1+2/(a-1))
        * (
            zeta2 + li2((a - 1) / (a + 1))
            - 1.0 / 6.0 * li2(2.0 / (1.0 + a))
            - 1.0 / 6.0 * li2(2.0 / (1.0 - a))
        )
        + 4 * (a * a - 1.0) * (zeta3 - li3((a - 1.0) / (a + 1.0)))
    );
    result *= -kf_up * kf_up / (8 * pi * pi * pi * q * q);
    return result;
}

double f_alda_ref(double q, double w, double rs, void* params) {
    double kf_up = // (9/2*pi)^(1/3) * 1/rs evaluated with wolfram alpha
        2.4179879310247044610754502699589274084169211084857923335384381130 / rs;
    double delta = kf_up / pi;
    return -6 * pi * pi * delta / (kf_up * kf_up * kf_up);
}

double f_pert_ref(double q, double w, double rs, void* params) {
    double a, b;
    a = chi_0_ref(q, w, rs,NULL);
    b = chi_1_ref(q, w, rs,NULL);
    return b / (a * a);
}

double f_pert(double q, double w, double rs, void* params) {
    double a, b;
    a = chi_0(q, w, rs, params);
    b = chi_1(q, w, rs, params);
    return b / (a * a);
}


double f_alda(double q, double w, double rs, void* params) {
    double kf_up = // (9/2*pi)^(1/3) * 1/rs evaluated with wolfram alpha
        2.4179879310247044610754502699589274084169211084857923335384381130 / rs;
    double delta = ((double*)params)[0];
    return -6 * pi * pi * delta / (kf_up * kf_up * kf_up);
}

double f_pert_w(double q, double w, double rs,void *params){
    return  f_pert(0,w,rs,params);
}
double f_pert_q(double q, double w, double rs,void *params){
    static double result, lastq=-1;
    if(q != lastq){result =f_pert(q,0,rs,params);}
    return result;
}

