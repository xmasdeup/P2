#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sndfile.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>

#include "vad.h"
#include "vad_docopt.h"
#include "pav_analysis.h"

#define DEBUG_VAD 0x1
#define PARAMETERS 5
#define DECISIONS 3

int main(int argc, char *argv[]) {
  int verbose = 0; /* To show internal state of vad: verbose = DEBUG_VAD; */

  SNDFILE *sndfile_in, *sndfile_out = 0;
  SF_INFO sf_info;
  FILE *vadfile;
  int n_read = 0, i;

  VAD_DATA *vad_data;
  VAD_STATE state, last_state;

  float *buffer, *buffer_zeros, *hamming_frame;
  unsigned int hamming_size = 3;
  int frame_size;         /* in samples */
  float frame_duration;   /* in seconds */
  unsigned int t, last_t; /* in frames */
  unsigned int temp_t;

  char	*input_wav, *output_vad, *output_wav;

  DocoptArgs args = docopt(argc, argv, /* help */ 1, /* version */ "2.0");

  verbose    = args.verbose ? DEBUG_VAD : 0;
  input_wav  = args.input_wav;
  output_vad = args.output_vad;
  output_wav = args.output_wav;
  float sensitivity = atof(args.sensitivity);
  float umbral1 = 0;
  float ***inverse_covariance;
  float **median;
  float *det_covariance;
  float *zeros;
  int newPosition;
  int number;

  if (input_wav == 0 || output_vad == 0) {
    fprintf(stderr, "%s\n", args.usage_pattern);
    return -1;
  }

  /* Open input sound file */
  if ((sndfile_in = sf_open(input_wav, SFM_READ, &sf_info)) == 0) {
    fprintf(stderr, "Error opening input file %s (%s)\n", input_wav, strerror(errno));
    return -1;
  }

  if (sf_info.channels != 1) {
    fprintf(stderr, "Error: the input file has to be mono: %s\n", input_wav);
    return -2;
  }

  /* Open vad file */
  if ((vadfile = fopen(output_vad, "wt")) == 0) {
    fprintf(stderr, "Error opening output vad file %s (%s)\n", output_vad, strerror(errno));
    return -1;
  }

  /* Open output sound file, with same format, channels, etc. than input */
  if (output_wav) {
    if ((sndfile_out = sf_open(output_wav, SFM_WRITE, &sf_info)) == 0) {
      fprintf(stderr, "Error opening output wav file %s (%s)\n", output_wav, strerror(errno));
      return -1;
    }
  }
  
  vad_data = vad_open(sf_info.samplerate,umbral1,sensitivity);
  
  /* Allocate memory for buffers */
  frame_size   = vad_frame_size(vad_data);
  zeros = (float *) malloc(frame_size * sizeof(float));
  buffer       = (float *) malloc(frame_size * sizeof(float));
  buffer_zeros = (float *) malloc(frame_size * sizeof(float));
  hamming_frame = (float *) malloc(frame_size * hamming_size * sizeof(float));

  median = (float **) malloc(DECISIONS * sizeof(float *));
  inverse_covariance = (float ***) malloc(DECISIONS*sizeof(float *));

  for(int k = 0;k<DECISIONS;k++)
  {
    median[k] = (float *) malloc(PARAMETERS * sizeof(float));
    inverse_covariance[k] = (float **) malloc(PARAMETERS*sizeof(float **));
    for(int z=0; z<PARAMETERS;z++) inverse_covariance[k][z] = (float *) malloc(PARAMETERS * sizeof(float));
  }
  det_covariance = (float *) malloc(DECISIONS * sizeof(float));
  for (i=0; i< frame_size; ++i) buffer_zeros[i] = zeros[i] = 0.0F;


  frame_duration = ((float) frame_size/ (float) sf_info.samplerate);
  last_state = ST_UNDEF;

  double frame_number = 1;

  assign_matrix_values(vad_data,inverse_covariance,median,det_covariance);
  vad_data->fp = fopen("oufile.txt","w");
  
  for (t = last_t = 0; ; t++) { /* For each frame ... */
    /* End loop when file has finished (or there is an error) */
    if  ((n_read = sf_read_float(sndfile_in, buffer, frame_size)) != frame_size) break;
  

    state = vad(vad_data, buffer, hamming_frame, hamming_size, frame_number, inverse_covariance, median, det_covariance, umbral1);
    
    if (verbose & DEBUG_VAD) vad_show_state(vad_data, stdout);

    if(frame_number == hamming_size) 
    {
      assign_matrix_values(vad_data,inverse_covariance,median,det_covariance);
    }
    
    if (sndfile_out != 0) {
        number = sf_write_float(sndfile_out, buffer, frame_size); 
    }


    if ((state != last_state)&& (last_state != ST_INIT)) {
      if ((t != last_t))
      { 
        if((((state == ST_MYBVOICE)&&(last_state ==ST_SILENCE))||((state==ST_MYBSILENCE)&&(last_state==ST_VOICE))))
        {
          if(state ==ST_MYBVOICE && last_state == ST_SILENCE) temp_t = t-1;
          if(state ==ST_MYBSILENCE && last_state == ST_VOICE) temp_t = t-1;

          }
        else if((((last_state == ST_MYBVOICE)&&(state ==ST_SILENCE))||((last_state==ST_MYBSILENCE)&&(state==ST_VOICE))))
        {
        }
        else 
        {
          if(temp_t!=0) 
          {
            fprintf(vadfile, "%.5f\t%.5f\t%s\n", last_t * frame_duration, temp_t * frame_duration, state2str(last_state));
            
          if((last_state==ST_MYBVOICE) && (state==ST_VOICE) && (sndfile_out != 0)){}
          else {
            last_t = temp_t;
            temp_t = 0;
          }          
          }
          else 
          {
            fprintf(vadfile, "%.5f\t%.5f\t%s\n", last_t * frame_duration, t * frame_duration, state2str(last_state));
            last_t = t;
          }

        }
      }

      if (sndfile_out != 0) 
      {
      if((last_state==ST_MYBVOICE) && (state==ST_VOICE)){

        newPosition = sf_seek(sndfile_out, last_t*frame_size, SEEK_SET); 
        
        for(int i=0; i<(temp_t-last_t);i++){
          number = sf_write_float(sndfile_out, zeros, frame_size);
        }
        newPosition = sf_seek(sndfile_out, (t+1)*frame_size, SEEK_SET);
        
        last_t = temp_t;
        temp_t = 0;      
        }
      }


      last_state = state;

    }
    if((state!= last_state)&&(last_state == ST_INIT)) last_state = state;


    frame_number++;
  }

  state = vad_close(vad_data);
  /* TODO: what do you want to print, for last frames? */
  if ((t != last_t))
  {
    fprintf(vadfile, "%.5f\t%.5f\t%s\n", last_t * frame_duration, t * frame_duration + n_read / (float) sf_info.samplerate, state2str(state));
  }

  if(state==ST_SILENCE){
    newPosition = sf_seek(sndfile_out, last_t*frame_size, SEEK_SET); 
     for(int i=0; i<(t-last_t);i++){
          number = sf_write_float(sndfile_out, zeros, frame_size);
        }
  }

  /* clean up: free memory, close open files */
  free(buffer);
  free(buffer_zeros);
  free(inverse_covariance);
  free(median);
  free(det_covariance);
  free(hamming_frame);
  sf_close(sndfile_in);
  fclose(vadfile);
  if (sndfile_out) sf_close(sndfile_out);
  return 0;
}
