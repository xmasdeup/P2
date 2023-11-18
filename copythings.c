#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sndfile.h>

#include "vad.h"
#include "vad_docopt.h"

#define DEBUG_VAD 0x1

int main(int argc, char *argv[]) {
  int verbose = 0; /* To show internal state of vad: verbose = DEBUG_VAD; */

  SNDFILE *sndfile_in, *sndfile_out = 0;
  SF_INFO sf_info;
  FILE *vadfile;
  int n_read = 0, i;

  VAD_DATA *vad_data;
  VAD_STATE state, last_state;

  float *buffer, *buffer_zeros;
  int frame_size;         /* in samples */
  float frame_duration;   /* in seconds */
  unsigned int t, last_t; /* in frames */

  char	*input_wav, *output_vad, *output_wav;

  DocoptArgs args = docopt(argc, argv, /* help */ 1, /* version */ "2.0");

  verbose    = args.verbose ? DEBUG_VAD : 0;
  input_wav  = args.input_wav;
  output_vad = args.output_vad;
  output_wav = args.output_wav;
  float umbral1 = atof(args.umbral1);

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
  vad_data = vad_open(sf_info.samplerate,umbral1);
  /* Allocate memory for buffers */
  frame_size   = vad_frame_size(vad_data);
  buffer       = (float *) malloc(frame_size * sizeof(float));
  buffer_zeros = (float *) malloc(frame_size * sizeof(float));
  for (i=0; i< frame_size; ++i) buffer_zeros[i] = 0.0F;

  frame_duration = (float) frame_size/ (float) sf_info.samplerate;
  last_state = ST_UNDEF;

  for (t = last_t = 0; ; t++) { /* For each frame ... */
    /* End loop when file has finished (or there is an error) */
    if  ((n_read = sf_read_float(sndfile_in, buffer, frame_size)) != frame_size) break;

    if (sndfile_out != 0) {
      /* TODO: copy all the samples into sndfile_out */
    int nframes = sf_write_float(sndfile_out, buffer, frame_size);
    }

    state = vad(vad_data, buffer);
    if (verbose & DEBUG_VAD) vad_show_state(vad_data, stdout);

    /* TODO: print only SILENCE and VOICE labels */
    /* As it is, it prints UNDEF segments but is should be merge to the proper value */
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
            
          if((last_state==ST_MYBVOICE) && (state==ST_VOICE)){}
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
      last_state = state;

    }
    if((state!= last_state)&&(last_state == ST_INIT)) last_state = state;

    if (sndfile_out != 0) {
      /* TODO: go back and write zeros in silence segments */
      if((last_state==ST_MYBVOICE) && (state==ST_VOICE)){
        sf_count_t newPosition = sf_seek(sndfile_out, last_t*frame_size, SEEK_SET); 
        for(int i=0; i<(temp_t-last_t);i++){
          int num = sf_write_float(sndfile_out, zeros, frame_size);
        }
        sf_count_t newPosition = sf_seek(sndfile_out, (t+1)*frame_size, SEEK_END);
        last_t = temp_t;
        temp_t = 0;      
        }
    }
  }

  state = vad_close(vad_data);
  /* TODO: what do you want to print, for last frames? */
  if (t != last_t)
    fprintf(vadfile, "%.5f\t%.5f\t%s\n", last_t * frame_duration, t * frame_duration + n_read / (float) sf_info.samplerate, state2str(state));

  if(state==ST_SILENCE){
    sf_count_t newPosition = sf_seek(sndfile_out, last_t*frame_size, SEEK_SET); 
     for(int i=0; i<(t-last_t);i++){
          int num = sf_write_float(sndfile_out, zeros, frame_size);
        }
  }

  /* clean up: free memory, close open files */
  free(buffer);
  free(buffer_zeros);
  sf_close(sndfile_in);
  fclose(vadfile);
  if (sndfile_out) sf_close(sndfile_out);
  return 0;
}
