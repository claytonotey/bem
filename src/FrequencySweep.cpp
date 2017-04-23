#include "FrequencySweep.h"
#include <stdio.h>

FrequencySweep :: FrequencySweep(char *filename)
{
  FILE *file = fopen(filename, "r");
  
  if( file != NULL ) {
    char line [128];
    while( fgets(line, sizeof line, file ) != NULL ) {
      f.push_back(atof(line));
    }
    fclose(file);
  } else {
    perror(filename);
  }
}
