#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>

#ifndef PAV_ANALYSIS_H
#define PAV_ANALYSIS_H


float compute_power(const float *x, unsigned int N);

float compute_am(const float *x, unsigned int N);

float compute_zcr(const float *x, unsigned int N, float fm);

float compute_LCP(float* x, unsigned int N, float* lpc_coeffs);

float compute_normalized_autocorrelation(float* x, unsigned int N);

float compute_normalized_prediction_error(float* x, unsigned int N, float* lpcCoeffs, float energy);

float compute_log_energy(float* x, unsigned int N,float umbral1);

void compute_windowed_frame(const float* x, unsigned int N, float* frame);

void covariance_method(float *input, unsigned int N, float *lpc_coeffs);

float compute_decision(float** covariance, float* median, float* measurements, float det);

float compute_gaussian(const gsl_matrix* covariance, const gsl_vector* median,const gsl_vector* measurements, float det);



#endif	/* PAV_ANALYSIS_H */
