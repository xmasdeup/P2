#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>

#include "pav_analysis.h"
#define ORDER 12

float compute_power(const float *x, unsigned int N) {
    float res = 0;
    //Realitza el sumatori de la potencia 
    for(unsigned int i=0;i<N;i++){
        res += x[i]*x[i];
    }
    return (10*log10(res/N));
}

float compute_am(const float *x, unsigned int N) {
    float res = 0;
    //Realitza el sumatori del valor absolut de la amplitud
    for(unsigned int i=0;i<N;i++){ 
        res+=fabs(x[i]);
    }
    return res/N;
}

float compute_zcr(const float *x, unsigned int N, float fm) {
    float res = 0;
    //Realitza el sumatori dels cops que passa per zero la senyal
    for (long unsigned int i = 1; i<N ; i++) {
        if(x[i]*x[i-1] <0){
            res++;
        }
    }
    //float zcr = (fm / 2)/(N - 1) * res;
    return res;
}


float compute_LCP(float* x, unsigned int N, float* lpc_coeffs)
{
    float a1 = 0;

    covariance_method(x, N, lpc_coeffs);

    a1 = lpc_coeffs[1];
    return a1;
}   

void covariance_method(float *input, unsigned int N, float *lpc_coeffs)
{
    float autocorr[ORDER+1];
    float lpc_temp[ORDER+1];

    for(int lag = 0; lag<ORDER+1; lag++)
    {
        autocorr[lag] = 0.0;
        for(int n =0; n< N-lag; n++)
        {
            autocorr[lag] += input[n]*input[n+lag];
        }
    }

    float E = autocorr[0];
    float k = 0;
    lpc_coeffs[0] = 1.0;
    lpc_temp[0] =1.0;
    float sum = 0;

    for(int i =1; i< ORDER+1; i++)
    {
        sum = autocorr[i];
        for(int j = 1; j<i ; j++)
        {
            sum += lpc_coeffs[j] * autocorr[i-j];
        }

        lpc_coeffs[i] = k = -sum/E;

        E = E * (1-k*k);

        for(int z = 1; z<i;z++)
        {
            lpc_temp[z] = lpc_coeffs[z] + k* lpc_coeffs[i-z];
        }

        for(int q = 1; q<i; q++)
        {
            lpc_coeffs[q] = lpc_temp[q];
        }

    }



}

float compute_normalized_autocorrelation(float* x, unsigned int N)
{
    float numerator = 0.0;
    float denominator1 = 0.0;
    float denominator2 = 0.0;

    for(int n = 0; n<N; n++)
    {
        numerator += x[n]*x[n+1];
       denominator1 += x[n]*x[n];
    }
    for(int n = 0; n<N-1; n++)
    {
        denominator2 += x[n+1]*x[n+1];
    }
    float value = numerator/ sqrt(denominator1*denominator2);
    return value;
}

float compute_normalized_prediction_error(float* x, unsigned int N, float* lpcCoeffs, float energy)
{

    float prediction_error = 0.0;
    float sum_of_covariance = 0.0;
    float covariance[ORDER+1];
    float sum =0;
    for(int k=1; k<= ORDER; k++){
        covariance[k] = 0.0;
        for(int n=k; n<N; n++)
        {
            covariance[k] += x[n]*x[n-k];

        }
        covariance[k] /= N;
    }
    covariance[0] = 0.0;
    for(int q = 0; q<N; q++)
    {
        covariance[0] += x[q]*x[q];
    }
    covariance[0]/= N;
    for(int j = 1; j<=ORDER;j++)
    {
        sum_of_covariance += lpcCoeffs[j] *covariance[j];
        
    }
    sum = sum_of_covariance + covariance[0];
    sum = fabs(sum);

    prediction_error = energy - 10* log10(0.00000000001+sum);
    return prediction_error;
}

float compute_log_energy(float* x, unsigned int N,float umbral1)
{   
    float log_energy = 0.0;

    for(int n=0; n<N; n++)
    {
        log_energy += x[n]*x[n];
    }

    log_energy = 10*log10(log_energy/N + 0.0000000000001); 

    return log_energy;
}

void compute_windowed_frame(const float* x, unsigned int N, float* frame)
{

    for (int n = 0; n < N; n++) {
        frame[n] = x[n]*(0.53836 - 0.46164 * cos(2 * M_PI * n / (N - 1)));
    }

}

float compute_decision(float **covariance, float *median, float *measurements, float det)
{
    int size = 5;

    gsl_vector* mean = gsl_vector_alloc(size);
    gsl_matrix* cov = gsl_matrix_alloc(size,size);
    gsl_vector* meas = gsl_vector_alloc(size);

    for(int i = 0; i< size; i++)
    {
        gsl_vector_set(mean, i, median[i]);
        gsl_vector_set(meas, i, measurements[i]);
        
        for(int j = 0; j<size; j++)
        {
            gsl_matrix_set(cov, i, j, covariance[i][j]);
        }
    }

    float value = compute_gaussian(cov,mean,meas,det);

    gsl_vector_free(mean);
    gsl_matrix_free(cov);
    gsl_vector_free(meas);

    return value;
}

float compute_gaussian(const gsl_matrix *covariance, const gsl_vector *median, const gsl_vector *measurements,float det)
{
    int size = median->size;
    gsl_vector* result = gsl_vector_alloc(size);
    gsl_vector* difference1 = gsl_vector_alloc(size);
    gsl_vector* difference2 = gsl_vector_alloc(size);

    gsl_vector_memcpy(difference1, measurements);
    gsl_vector_memcpy(difference2, measurements);
    gsl_vector_sub(difference1,median);
    gsl_vector_sub(difference2,median);

    double gaussian_exponent = 0;
    
    gsl_blas_dtrmv(CblasUpper,CblasNoTrans,CblasNonUnit,covariance,difference1);

    gsl_blas_ddot(difference2, difference1, &gaussian_exponent);

    float decision = gaussian_exponent;

    gsl_vector_free(result);
    gsl_vector_free(difference1);
    gsl_vector_free(difference2);

    return decision;
}

