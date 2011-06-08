#include<bar.h>

//
// is_time or not time to do something ?
//
/*! \fn int 
	
*/
bool is_time(double time_interval_s)
{
   double dt;
   static time_t time_cur;
   static time_t time_prev;
      
   time_cur = time(NULL);   
   dt = difftime(time_cur, time_prev);
   
   if (dt>time_interval_s) {
 
   time_prev = time_cur;
   return true;
   }

   return false;   
}


//
// Progress bar
//
int put_progress_bar(FILE* file, char * str_symbol, 
                     long step, long no_of_steps, long len_of_bar)
{
   long d;
   static long counter;
   
   if (step==0) counter=0;
   counter = counter + 1;
      
   d = no_of_steps / len_of_bar;
   
   if (counter==d)
   {
     fprintf(file, "%s", str_symbol); fflush(file);
     counter = 0;
   };
   
   return(0);
}

//
// Progress counter
//
/*! \fn int put_progress_counter(FILE* file, long step, 
                         long no_of_steps, long len_of_bar, int index)
   
   The procedure puts the string with a progress counter

      1%  2% ...
   
   to the file \a file. Suits to use in a for--cycle
   
   for(long step=0; step<no_of_steps; step++)
   {
   
   }
   
	\param step current step in the for--cycle
	
	\param no_of_steps how much steps in the for--cycle
	
	\param len_of_bar how much output to produce during \a step
	= 0 .. \a no_of_steps.
	
	\param index index of the counter if you want to put several 
	counters (for instance to several files)
	
	\return 0.
	
	\remarks There can be up to 1024 counters which are distinguished 
	by index \a index=0..1023.
*/
int put_progress_counter(FILE* file, long step, 
                         long no_of_steps, long len_of_bar, int index)
{
   long d;
   static long counter[1024], counter_newline[1024];
   
   if (step==0) 
   {
     counter[index]=0; 
     counter_newline[index]=0;
   }
   
   counter[index] = counter[index] + 1;
      
   d = no_of_steps / len_of_bar;
   
   if (counter[index]==d)
   {
     fprintf(file, "%3ld%% ",(long)floor((100.0 * step)/no_of_steps));
     fflush(file); 
     counter[index] = 0;
     counter_newline[index] = counter_newline[index] + 1;
     
     if (counter_newline[index]==14) 
     {
       fprintf(file, "\n"); 
       counter_newline[index]=0;
     }
   };
   
   return(0);
}

//
// Progress counter
//
/*! \fn int put_progress_counter(FILE* file, const char* c_str_format, 
                                 double acc, int time_interval_s)
   
   The procedure puts double with the format string \a c_str_format
   to the streal \a file if the \a time_interval_s is over.
   This siuts for the self-consistent loops with a self-consistency
   criterion \a acc.
  
	\param file output stream
	
	\param c_str_format format string for the \a acc.
	
	\param acc double which is to be put to the \a file if 
	\a time_interval_s is over.
	
	\param time_interval_s time interval in seconds.
	
	\param ind the number of the counter.
	
	\return 0.
	
*/
int put_progress_counter(FILE* file, const char *c_str_format,
                         double acc, double time_interval_s, int ind)
{
   double dt;
   static long counter_newline[1024];
   static time_t time_cur;
   static time_t time_prev;
      
   time_cur = time(NULL);   
   dt = difftime(time_cur, time_prev);
   
   if (dt>time_interval_s) 
   {   
     time_prev = time_cur;   
     fprintf(file, c_str_format, acc);
     counter_newline[ind] = counter_newline[ind] + 1;
     if (counter_newline[ind]==6) 
     {
       fprintf(file, "\n"); 
       counter_newline[ind]=0;
     }
   };
   
   return(0);
}
