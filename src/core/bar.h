#ifndef bar_h
#define bar_h bar_h
#include<stdio.h>
#include<math.h>
#include<time.h>

bool is_time(double time_interval_s);
int put_progress_bar(FILE*,char *,long step,long no_of_steps,long len_of_bar);
int put_progress_counter(FILE*, long st, long no_st, long len, int ind);
int put_progress_counter(FILE* file, const char *c_str_format,
                         double acc, double time_interval_s, int ind);


#endif // bar_h
