#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Use the command
//  &(AUXZONE[ACTIVEOFFSET=1]:time_in_mins)
// to access the time stamp

int main() {
  FILE *fid_in;
  FILE *fid_out;
//  const char *src = "C:\\Users\\ijoynes\\Documents\\run_023\\particle_paths_2013-06-21.dat";
  const char *src = "C:\\Users\\ijoynes\\Documents\\field_of_view_study\\fov_2012.11.24_no_filter.dat";
  const char *dst = "C:\\Users\\ijoynes\\Documents\\field_of_view_study\\fov_2013-13-12_no_filter.dat";
  char temp[999];
  int i = 0;
  int mins;
  int secs;

  fid_in = fopen(src, "r");
  if ( fid_in == NULL ) {
    perror("Error opening file");
    return(EXIT_FAILURE);
  }

  fid_out = fopen(dst, "w");

  i = 900;
  while ( !feof(fid_in) ) {
  //while ( i <  ) {


    fgets(temp, 999, fid_in);
    fprintf(fid_out, "%s", temp);

    if ( strcmp(temp, "ZONE\n") == 0 && i >= 0 ) {
      printf("i = %d\n",i);
      mins = (int) floor(i/60.0);
      secs = i % 60;
    
      if ( mins < 10 & secs < 10 ) {
        fprintf(fid_out, "AUXDATA time_in_mins=\"0%d:0%d\"\n", mins, secs);
      }
      else if ( mins < 10 & secs >= 10 ) {
        fprintf(fid_out, "AUXDATA time_in_mins=\"0%d:%d\"\n", mins, secs);
      }
      else if ( mins >= 10 & secs < 10 ) {
        fprintf(fid_out, "AUXDATA time_in_mins=\"%d:0%d\"\n", mins, secs);
      }
      else if ( mins >= 10 & secs >= 10 ) {
        fprintf(fid_out, "AUXDATA time_in_mins=\"%d:%d\"\n", mins, secs);
      }
      
      i--;

    }

  }

  fclose(fid_in);
  fclose(fid_out);
  return(EXIT_SUCCESS);
}

