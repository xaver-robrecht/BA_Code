/**
 * @file spinwaves.h
 * @author Xaver Robrecht (xaver.robrecht@gmail.com)
 * @brief implements routines to generate spin wave dispersion relations
 * @version 1.0
 */

#ifndef SPINWAVES_H
#define SPINWAVES_H


//structure to store spin wave dispersion relations
typedef struct{
    double *q;
    double *omega;
    int n;
}spinwave_dispersion;

//functions to calculate spin wave dispersion(dispersion.c)

/**
 * @brief
 *
 * @param disp
 */
void free_dispersion(spinwave_dispersion *disp);

/**
 * @brief
 *
 * @param q_step stepsize in generated dispersion relation
 * @param q_max calculation gets aborted if q>q_max. if q_max<=0  the complete spectrum is calculated.
 * @param rs wiegner seitz radius
 * @param delta energy difference at fermi wave vector
 * @param chi used to calculate spin susceptibility
 * @param f exchange correlation kernel
 * @param params additional parameters passed to chi and f.
 * @return spinwave_dispersion*
 */
spinwave_dispersion *dispersion(double q_step,double qmax,double rs, double delta,
                    double (*chi)(double q,double w, double rs,void *params),
                    double (*f)(double q,double w, double rs,void *params),void *params);

/**
 * @brief returns upper bound of stoner continuum at given q for given kf, delta.
 *
 * @param q magnon wave vector
 * @param delta difference between spin channels at k=kf
 * @param kf Fermi wave vector
 * @return lower bound of stoner continuum at given q
 */
double upper_stoner_bound(double q, double delta,double kf);

/**
 * @brief returns lower bound of stoner continuum at given q for given kf, delta.
 *
 * @param q magnon wave vector
 * @param delta difference between spin channels at k=kf
 * @param kf Fermi wave vector
 * @return lower bound of stoner continuum at given q
 */
double lower_stoner_bound(double q, double delta,double kf);

// susceptibilities & exchange kernels(chi.c)
//reference functions(chi.c)
/**
 * @brief Kohn-Sham spin suceptibility for Jellium @see BA Borrameo
 *
 * @param q magnon wave vector.
 * @param w magnon angular frequency.
 * @param rs Wigner-Seitz radius in atomic units
 * @param params pointer to pass additional parameters
 * @return Kohn-Sham suceptibility.
 */
double chi_0_ref(double q, double w, double rs,void *params);

/**
 * @brief Kohn-Sham spin suceptibility for Jellium @see BA Borrameo
 *
 * @param q magnon wave vector.
 * @param w magnon angular frequency.
 * @param rs Wigner-Seitz radius in atomic units
 * @param params expects pointer to delta value(double)
 * @return Kohn-Sham suceptibility.
 */
double chi_0(double q, double w, double rs,void *params);/**

 * @brief Kohn-Sham spin susceptibility for Jellium
 *
 * @param q magnon wave vector.
 * @param w magnon angular frequency.
 * @param rs Wigner-Seitz radius in atomic units
 * @param params expects pointer to (double[]){delta,lambda}
 * @return perturbation-theory susceptibility.
 */
double chi_1(double q, double w, double rs,void *params);

/**
 * @brief exchange correlation kernel in ALDA. Does not actually depend
 * on q, w or additional parameters.
 *
 * @param q magnon wave vector.
 * @param w magnon angular frequency.
 * @param rs Wigner-Seitz radius in atomic units
 * @param params pointer to pass additional params
 * @return exchange correlation kernel in ALDA approximation .
 */
double f_alda_ref(double q, double w, double rs,void *params);

/**
 * @brief perturbation theory exchange correlation kernel.
 *
 * @param q magnon wave vector.
 * @param w magnon angular frequency.
 * @param rs Wigner-Seitz radius in atomic units
 * @param params pointer to pass additional params
 * @return first order perturbation theory exchange correlation kernel.
 */
double f_pert_ref(double q, double w, double rs,void *params);
double f_pert(double q, double w, double rs,void *params);
double f_alda(double q, double w, double rs,void *params);
double f_pert_w(double q, double w, double rs,void *params);
double f_pert_q(double q, double w, double rs,void *params);

void print_spectrum(double rs, double qstep,double qmax,double lambda);
void print_reference_spectrum(double rs,double qstep,double qmax);
#endif //SPINWAVES_H
