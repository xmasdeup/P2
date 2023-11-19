#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pav_analysis.h"
#include "vad.h"

#define ORDER 12
#define PARAMETERS 5

FILE *fp;

const float FRAME_TIME = 10.0F; 


const char *state_str[] = {
  "S", "S", "V", "S", "S","V"
};

const char *state2str(VAD_STATE st) {
  return state_str[st];
}

typedef struct {
  float zcr;
  float p;
  float am;
  float umbral1;
  float log_energy;
  float norm_correlation;
  float norm_prediction_error;
  float ***covariance;
  float **median;
  float *hamming_frame;
  float first_coeff;
  float *lpc_coeffs;
  float *measurements;
  float voiced;
  float unvoiced;
  float silence;
  float frame_number;
  FILE *fp;

} Features;


Features compute_features(const float *x, int N, float fm, float ***cov, float **med, float *det_cov, float umbral1, FILE *fp, int count) {

  Features feat;

  feat.hamming_frame = (float *) malloc(N * sizeof(float));
  feat.lpc_coeffs = (float *) malloc((ORDER+1)*sizeof(float));
  feat.measurements = (float*) malloc((PARAMETERS)*sizeof(float));

  compute_windowed_frame(x,N,feat.hamming_frame);
  feat.zcr = compute_zcr(feat.hamming_frame,N,fm);
  feat.log_energy = compute_log_energy(feat.hamming_frame,N,umbral1);
  feat.norm_correlation = compute_normalized_autocorrelation(feat.hamming_frame,N);
  feat.first_coeff = compute_LCP(feat.hamming_frame,N,feat.lpc_coeffs);
  feat.norm_prediction_error = compute_normalized_prediction_error(feat.hamming_frame,N,feat.lpc_coeffs,feat.log_energy);
  
  feat.am = compute_am(feat.hamming_frame,N);
  feat.measurements[0] = feat.zcr;
  feat.measurements[1] = feat.log_energy;
  feat.measurements[2] = feat.norm_correlation;
  feat.measurements[3] = feat.first_coeff;
  feat.measurements[4] = feat.norm_prediction_error;


  feat.voiced = abs(compute_decision(cov[0],med[0],feat.measurements,det_cov[0]));
  feat.unvoiced = abs(compute_decision(cov[1],med[1],feat.measurements,det_cov[1]));
  feat.silence = abs(compute_decision(cov[2],med[2],feat.measurements,det_cov[2]));

    fprintf(fp,"Voiced: %.3f\tUnvoiced: %.3f\tSilence: %.3f\t Count: %d\t\n ",feat.voiced,feat.unvoiced,feat.silence,count);
  fprintf(fp,"Energy: %.3f\tZCR: %.3f\t Correl: %.3f\t First:%.3f\t Error:%.3f\t",feat.log_energy,feat.zcr,feat.norm_correlation,feat.first_coeff,feat.norm_prediction_error);
  if((feat.voiced>feat.silence)&&(feat.unvoiced>feat.silence))   fprintf(fp,"silence\n");
  if(feat.silence>feat.voiced) fprintf(fp,"Voice\n");
  else if(feat.silence>feat.unvoiced) fprintf(fp,"Unvoiced\n");
  
  free(feat.hamming_frame);
  free(feat.lpc_coeffs);
  free(feat.measurements);
  return feat;
}


VAD_DATA * vad_open(float rate, float umbral1, float sensitivity) {
  VAD_DATA *vad_data = malloc(sizeof(VAD_DATA));
  vad_data->state = ST_INIT;
  vad_data->sampling_rate = rate;
  vad_data->frame_length = rate * FRAME_TIME * 1e-3;
  vad_data->umbral1 = umbral1;
  vad_data->sensitivity = sensitivity;
  vad_data->count = 0;
  return vad_data;
}

VAD_STATE vad_close(VAD_DATA *vad_data) {

  VAD_STATE state = vad_data->state;

  free(vad_data);
  return state;
}

unsigned int vad_frame_size(VAD_DATA *vad_data) {
  return vad_data->frame_length;
}

VAD_STATE vad(VAD_DATA *vad_data, float *x, float *hamm, unsigned int hamm_size, unsigned int frame_number, float ***inv_cov, float **med, float *det_cov, float umbral1) {


  if(frame_number <= hamm_size)
  {
    for(int i = 0; i< vad_data->frame_length; i++)
    {
      hamm[i+(vad_data->frame_length*(frame_number-1))] = x[i];
    }
  }
  else
  {
    for(int j = 0; j< (vad_data->frame_length)*(hamm_size-1); j++)
    {
      hamm[j] = hamm[j + (vad_data->frame_length)];
    }
    for(int k = 0; k< (vad_data->frame_length); k++)
    {
      hamm[k + (vad_data->frame_length)*(hamm_size-1)] = x[k];
    }
  }
  Features f;

  if(frame_number >= hamm_size)
  { 

    f = compute_features(hamm, (vad_data->frame_length) * hamm_size, vad_data->sampling_rate, inv_cov, med, det_cov, umbral1,vad_data->fp,vad_data->count);
    vad_data->last_feature = f.log_energy; /* save feature, in case you want to show */
  }
  switch (vad_data->state) {
  case ST_INIT:
    
    if(frame_number == hamm_size)
    { 
      vad_data->umbral1 = umbral1 + f.log_energy;
      vad_data->error = f.norm_prediction_error;
      vad_data->state = ST_SILENCE;
      vad_data->noise = f.zcr;
      vad_data->a0 = f.am;

    }
    break;

  case ST_SILENCE:
    if (((f.voiced < f.silence)&& (f.am > vad_data->a0 * 3))||(f.unvoiced < f.silence))
    {
      vad_data->count =0;
      vad_data->state = ST_MYBVOICE;
    }
    break;

  case ST_VOICE:
    if(((f.zcr > 29) && (f.norm_correlation >0.69) && (f.am< 2.3*vad_data->a0))||(f.log_energy<(vad_data->umbral1 + vad_data->sensitivity)))
    {
      vad_data->count = 0;
      vad_data->state = ST_MYBSILENCE;
    }
    else if ((f.silence < f.voiced) && (f.silence < f.unvoiced))
    {
      vad_data->count = 0;
      vad_data->state = ST_MYBSILENCE;
    }
    break;

  case ST_UNDEF:
    break;
  
  case ST_MYBVOICE:

      if(((f.zcr >29) && (f.am < 2.75*vad_data->a0))||(f.log_energy<(vad_data->umbral1 + vad_data->sensitivity))) vad_data->state = ST_SILENCE;

      else if((f.silence<f.voiced) && (f.silence<f.unvoiced)) vad_data->state=ST_SILENCE;

      else if((vad_data->count <4)) vad_data->count++;

      else vad_data->state = ST_VOICE; 
    break;

  case ST_MYBSILENCE:
    if(((f.zcr > 29) && (f.norm_correlation >0.69) && (f.am< 2.3*vad_data->a0))||(f.log_energy<(vad_data->umbral1 + vad_data->sensitivity)))
    {
      if((vad_data->count <13)) vad_data->count++;

      else vad_data->state = ST_SILENCE;  
    }
    
    else if (f.am > 3*vad_data->a0){vad_data->state = ST_VOICE;}

    else if((f.silence>f.voiced)||(f.silence>f.unvoiced)) vad_data->state=ST_VOICE;
    
    else if((vad_data->count < 14)) vad_data->count++;
    
    else vad_data->state = ST_SILENCE;    
    
    break;
  }

    return vad_data->state;
}

void vad_show_state(const VAD_DATA *vad_data, FILE *out) {
  fprintf(out, "%d\t%f\n", vad_data->state, vad_data->last_feature);
}

void assign_matrix_values(const VAD_DATA *vad_data, float ***inv_cov, float **mean, float *cov_det)
{
  //Determinants
  cov_det[0] = 0.1283;
  cov_det[1] = 0.0183;
  cov_det[2] = 349.4;
  
  float energy_distance = 32.84 -vad_data->umbral1;

  //Voiced means
  mean[0][0] = 20;
  mean[0][1] = 62.054 - energy_distance;
  mean[0][2] = 0.9733;
  mean[0][3] = -2.1792 ;
  mean[0][4] = 24.37;
  //Unvoiced means
  mean[1][0] = 110;
  mean[1][1] = 45.67 -energy_distance ;
  mean[1][2] = 0.3397;
  mean[1][3] = -0.8136;
  mean[1][4] = 11.80;

  //Silence means
  mean[2][0] = 45;
  mean[2][1] = 32.99 - energy_distance;
  mean[2][2] = 0.8784;
  mean[2][3] =  -1.34;
  mean[2][4] =  10.98; 


  //Voiced inverse covariance
    inv_cov[0][0][0] = 0.0485;
    inv_cov[0][0][1] = -0.00244;
    inv_cov[0][0][2] = 7.0;
    inv_cov[0][0][3] = 0.142;
    inv_cov[0][0][4] = 0.0229; 

    inv_cov[0][1][0] = 0;
    inv_cov[0][1][1] = 0.0097;
    inv_cov[0][1][2] = -0.23;
    inv_cov[0][1][3] = 0.0322;
    inv_cov[0][1][4] = -0.00214;

    inv_cov[0][2][0] = 0;
    inv_cov[0][2][1] = 0;
    inv_cov[0][2][2] = 1726;
    inv_cov[0][2][3] = 12.76;
    inv_cov[0][2][4] = 0.708;

    inv_cov[0][3][0] = 0;
    inv_cov[0][3][1] = 0;
    inv_cov[0][3][2] = 0;
    inv_cov[0][3][3] = 5.49;
    inv_cov[0][3][4] = 0.346;
    
    inv_cov[0][4][0] = 0;
    inv_cov[0][4][1] = 0;
    inv_cov[0][4][2] = 0;
    inv_cov[0][4][3] = 0;
    inv_cov[0][4][4] = 0.0557;

  //Unvoiced inverse covariance
    inv_cov[1][0][0] = 0.00613;
    inv_cov[1][0][1] = -0.00103;
    inv_cov[1][0][2] = 0.696;
    inv_cov[1][0][3] = 0.0295;
    inv_cov[1][0][4] = 0.00587;

    inv_cov[1][1][0] = 0;
    inv_cov[1][1][1] = 0.00867;
    inv_cov[1][1][2] = -0.0349;
    inv_cov[1][1][3] = 0.00615;
    inv_cov[1][1][4] = 0.103;

    inv_cov[1][2][0] = 0;
    inv_cov[1][2][1] = 0;
    inv_cov[1][2][2] = 88.6;
    inv_cov[1][2][3] = 5.96;
    inv_cov[1][2][4] = 0.649;

    inv_cov[1][3][0] = 0;
    inv_cov[1][3][1] = 0;
    inv_cov[1][3][2] = 0;
    inv_cov[1][3][3] = 2.57;
    inv_cov[1][3][4] = 0.305;

    inv_cov[1][4][0] = 0;
    inv_cov[1][4][1] = 0;
    inv_cov[1][4][2] = 0;
    inv_cov[1][4][3] = 0;
    inv_cov[1][4][4] = 0.068;

    //Silence inverse covariance
    inv_cov[2][0][0] = 0.0148;
    inv_cov[2][0][1] = -0.00373;
    inv_cov[2][0][2] = 1.81;
    inv_cov[2][0][3] = 0.1459;
    inv_cov[2][0][4] = 0.0244;

    inv_cov[2][1][0] = 0;
    inv_cov[2][1][1] = 0.00861;
    inv_cov[2][1][2] = -0.429;
    inv_cov[2][1][3] = -0.0455;
    inv_cov[2][1][4] = -1.53;

    inv_cov[2][2][0] = 0;
    inv_cov[2][2][1] = 0;
    inv_cov[2][2][2] = 331.9;
    inv_cov[2][2][3] = 11.77;
    inv_cov[2][2][4] = 1.13;

    inv_cov[2][3][0] = 0;
    inv_cov[2][3][1] = 0;
    inv_cov[2][3][2] = 0;
    inv_cov[2][3][3] = 10.8;
    inv_cov[2][3][4] = 1.27;

    inv_cov[2][4][0] = 0;
    inv_cov[2][4][1] = 0;
    inv_cov[2][4][2] = 0;
    inv_cov[2][4][3] = 0;
    inv_cov[2][4][4] = 0.199;
 
}


