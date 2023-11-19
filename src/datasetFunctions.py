import numpy as np
import librosa
from scipy.fft import fft, ifft
ORDER = 12

def compute_windowed_frame(x,hamm):
    frame = np.multiply(x,hamm)
    return frame



def compute_zcr(x , N, fm):

    res = 0

    for ii in range(1,N):
        if(x[ii]*x[ii-1]<0):
            res += 1

    #zcr = (fm/2)/(N-1) * res

    return res

def covariance_method(input, N):

    E = np.zeros(ORDER+1)
    correlation = np.zeros(ORDER+1)
    alpha = np.zeros(ORDER+1)
    alpha_t = np.zeros(ORDER+1)

    for lag in range(0,ORDER+1):
        for n in range(0,N-lag):
            correlation[lag] += input[n]*input[n+lag]

    k = 0
    E = correlation[0]
    alpha[0] = 1.0
    alpha_t[0] = 1.0

    for i in range(1,ORDER+1):
        
        sum = correlation[i]
        for j in range(1,i):

            sum += alpha[j] * correlation[i-j]
        
        alpha[i] = k = -sum/E
        
        E = E*(1-k**2)

        for z in range(1,i):
            alpha_t[z] = alpha[z] + k*alpha[i-z]
        
        for y in range(1,i):
            alpha[y] = alpha_t[y]



    lpc_coeffs = alpha

    return lpc_coeffs

        

def compute_LPC(x, N):
    a1 = 0

    lpc_coeffs = covariance_method(x,N)
    
    a1 = lpc_coeffs[1]
    return a1, lpc_coeffs


def compute_normalize_autocorrelation(x,N):
    numerator = 0
    denominator1 = 0
    denominator2 = 0

    for ii in range(0, N-1):
        numerator += x[ii]*x[ii+1]
        denominator2 += x[ii+1]*x[ii+1]

    for jj in range(0,N):
        denominator1 += x[jj]*x[jj]

    value = numerator / np.sqrt(denominator1*denominator2)
    return value

def compute_normalized_prediction_error(x, N, lpc_coeffs, energy):
    prediction_error = 0.0
    sum_of_covariance = 0.0
    covariance = np.empty(ORDER+1)

    for kk in range(1,ORDER+1):
        covariance[kk] = 0.0
        for nn in range(kk,N-kk):
            covariance[kk] += x[nn]*x[nn-kk]
        covariance[kk]/= N-kk

    for qq in range(0,N):
        covariance[0] += x[qq]*x[qq]
    covariance[0] /= N

    for jj in range(1,ORDER+1):
        sum_of_covariance += lpc_coeffs[jj] * covariance[jj]

    sum_of_covariance = np.fabs(sum_of_covariance + covariance[0])
    prediction_error = energy - 10*np.log10(0.000001 +sum_of_covariance)

    return prediction_error


def compute_log_energy(x,N):
    log_energy = 0

    for ii in range(0,N):
        log_energy += x[ii] * x[ii]

    log_energy = 10*np.log10(log_energy/N + 0.00000000001)
    
    return log_energy


