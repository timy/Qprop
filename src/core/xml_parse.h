#ifndef xml_parse_h
#define xml_parse_h

#include<string.h>
#include<stdio.h>
#include<complex>
#define complex std::complex<double>

int xml_parse_fgetcd(complex *cdres,  FILE *file, 
  const char*c_str_stag, const char*c_str_ftag);

int xml_parse_fgetd(double *dres,  FILE *file, 
  const char*c_str_stag, const char*c_str_ftag);

int xml_parse_fgeti(int *ires,  FILE *file, 
  const char*c_str_stag, const char*c_str_ftag);

int xml_parse_fgetl(long *lres,  FILE *file, 
  const char*c_str_stag, const char*c_str_ftag);

int xml_parse_fgets(char *c_str_res, long isize,  FILE *file, 
   const char*c_str_stag, const char*c_str_ftag);

int xml_parse_fseektag(FILE *file, const char*c_str_tag);

#endif // xml_parse_h
