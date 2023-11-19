#ifndef _VAD_H
#define _VAD_H
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>

typedef enum {ST_UNDEF=0, ST_SILENCE, ST_VOICE, ST_INIT, ST_MYBVOICE, ST_MYBSILENCE} VAD_STATE;

const char *state2str(VAD_STATE st);


typedef struct {
  VAD_STATE state;
  float umbral1;
  float sensitivity;
  float ***covariance;
  float **median;
  float sampling_rate;
  unsigned int frame_length;
  unsigned int count;
  FILE *fp;
  float a0;
  float last_feature; /* for debuggin purposes */
} VAD_DATA;


VAD_DATA *vad_open(float sampling_rate, float umbral1, float len);


unsigned int vad_frame_size(VAD_DATA *);


VAD_STATE vad(VAD_DATA *vad_data, float *x, float *hamm ,unsigned int hamm_size, unsigned int frame_number, float ***inv_cov, float **med, float *det_cov,float umbral1);


VAD_STATE vad_close(VAD_DATA *vad_data);


void vad_show_state(const VAD_DATA *, FILE *);

void assign_matrix_values(const VAD_DATA *vad_data,float ***inv_cov, float **mean, float *cov_det);

#endif
