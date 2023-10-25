#include <math.h>
#include "pav_analysis.h"

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
    float zcr = (fm / 2)/(N - 1) * res;
    return zcr;
}

float compute_powerham(const float* x, const float* w, unsigned int N)
{
    //Realitza el calcul de potencia en el cas de tenir una finestra de 
    //Hamming
    float den = 0;
    float num = 0;
    float res = 0;
    for (long unsigned int n = 0; n < N; n++)
    {
        num += (x[n] * w[n] * x[n] * w[n]);
        den += (w[n] * w[n]);
    }
    res = 10 * log10(num/den);
    return res;
}
