import numpy as np
from scipy.io import wavfile
import glob
import os
import datasetFunctions as dt

hamm_size = 3
frame = 100
N = frame * hamm_size

voiced_path ='./training/silence'
file_type = ".wav"
#files = glob.glob(os.path.join(voiced_path,file_type))
files = os.listdir(voiced_path)
print(len(files))

hamm = []

for jj in range(0,N):
    num = 0.53836 - 0.46164*np.cos(2*np.pi*jj/(N-1))
    hamm.append(num)
hamm = np.array(hamm)

i = 0
files_data =[]

for name in files:
    name = "./training/silence/"+name
    samplerate, data = wavfile.read(name)
    files_data.append(data)

count = 0
total_energy = 0
total_first_LPC = 0
total_error = 0
total_norm_correlation = 0
total_zcr = 0

nsamples = 0
total_cov = np.zeros((5,5))
samples_list = []

for samples in files_data:
    index = 0 
    hamm_window = np.zeros(N)
    while len(samples)-index>N :
        if index == 0:
            for ii in range(0,N):
                hamm_window[ii] = samples[ii]*hamm[ii]
            index = 100
        else:
            for jj in range(0,N):
                hamm_window[jj] = samples[jj+index]*hamm[jj]
            index +=100
        count +=1
        
        energy = dt.compute_log_energy(hamm_window,N)
        [first_LPC,lpc_coeffs] = dt.compute_LPC(hamm_window,N)
        zcr = dt.compute_zcr(hamm_window,N,samplerate)
        norm_correlation = dt.compute_normalize_autocorrelation(hamm_window,N)
        error = dt.compute_normalized_prediction_error(hamm_window,N,lpc_coeffs,energy)

        total_energy += energy
        total_first_LPC += first_LPC
        total_norm_correlation += norm_correlation
        total_error += error
        total_zcr += zcr

        samples_list.append([zcr,energy,norm_correlation,first_LPC,error])

        #print(count)
        #print(index)
    nsamples +=1
    print(nsamples)


print(len(samples_list))
cov = np.array(samples_list)
deviation = np.std(cov,axis = 0)
mean = np.mean(samples_list,axis = 0)

cov = cov#-mean
#trans_cov = np.matrix.getT(cov)
#total_cov += np.multiply(cov,trans_cov)
total_cov = np.cov(cov,rowvar=False,bias=True)

# print(total_cov)
# print(count)

mean_energy = total_energy/count
mean_LPC_first = total_first_LPC/count
mean_error = total_error/count
mean_zcr = total_zcr/count
mean_norm_correlation = total_norm_correlation/count

print("deviation")
print(deviation)
print("mean")
print(mean)

total_inv_cov = np.matrix.getI(total_cov)
#total_inv_cov = np.corrcoef(total_inv_cov)

det_cov = np.linalg.det(total_cov)
print("Silence")
print(det_cov)
print(total_inv_cov)
aa= 0
        