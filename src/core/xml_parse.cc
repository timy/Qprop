#include<xml_parse.h>
#define i_buf_len 1024
//
//
//
/*! \fn int xml_parse_fseektag(FILE *file, const char*c_str_tag)
   read characters until c_str_tag is reached. After the 
   file will be set to the next character.

   \param file input stream
   \param c_str_tag tag to be read e.g. <html>
   
   \return number of characters read until \a c_str_tag is met.
           negative integer if error occurs.
   
*/
int xml_parse_fseektag(FILE *file, const char*c_str_tag)
{
   int ires, ic, ilen, inum;   
   ic = 0; inum = 0; ilen = strlen(c_str_tag);
   do 
   {
     inum = inum + 1;
     ires = fgetc(file);
     if (c_str_tag[ic]==ires) {ic=ic+1;} else {ic=0;}
     if (ic==ilen) break;
   }
   while (ires != EOF);
   
   if (ic!=ilen) 
   {
     fprintf(stderr,"err: xml_parse_fseektag: can't find tag! %s\n",c_str_tag);
     inum=-1;
   }
   
   return inum;
};

//
//
//
/*! \fn int xml_parse_fgetcd(complex* cdres, FILE *file, 
   const char*c_str_stag, const char*c_str_ftag)
   
   read the double between tags \a c_str_stag, \a c_str_ftag.
   
   \param cdres result (complex<double> result).
   \param file input stream.
   \param c_str_stag start tag e.g. <duration>.
   \param c_str_ftag finish tag e.g. </duration>.
      
   \return negative integer if a error occurs, zero otherwise.
*/

int xml_parse_fgetcd(complex *cdres,  FILE *file, 
  const char*c_str_stag, const char*c_str_ftag)
{
   int ierr;
   char c_str_buf[i_buf_len];
   double re, im;

   ierr = xml_parse_fgets(c_str_buf, i_buf_len,  file, c_str_stag, c_str_ftag);
   if (ierr<0)
   {
     fprintf(stderr, "err: xml_parse_fgetcd(): wrong field %s\n", c_str_buf);
     *cdres = complex(0.0, 0.0);
     return ierr;
   }
   
   ierr = sscanf(c_str_buf, "(%lf,%lf)", &re, &im);
    
   if (ierr<0) 
   {
     fprintf(stderr, 
       "err: fgetcd(): smth is wrong with complex %s\n", c_str_buf);
     *cdres = complex(0.0, 0.0);
     return ierr;
   }
   
   *cdres = complex(re, im);
   return ierr;
};

//
//
//
/*! \fn int xml_parse_fgets(char *c_str_res, long isize,  FILE *file, 
   const char*c_str_stag, const char*c_str_ftag)
   
   read the string between tags \a c_str_stag, \a c_str_ftag.
   
   \param c_str_res result string.
   \parse isize size of the \a c_str_res.
   \param file input stream.
   \param c_str_stag start tag e.g. <html>.
   \param c_str_ftag finish tag e.g. </html>.
      
   \return negative integer if a error occurs, zero otherwise.
*/

int xml_parse_fgets(char *c_str_res, long isize,  FILE *file, 
   const char*c_str_stag, const char*c_str_ftag)
{
   long ires, ic, ic1, ilen;
   char * c_str_buf;
   
   if (strlen(c_str_stag)<1) 
   {
     fprintf(stderr,"err: xml_parse_fgets(): strlen(c_str_stag)<1!\n");
     return -1;
   };

   if (strlen(c_str_ftag)<1) 
   {
     fprintf(stderr,"err: xml_parse_fgets(): strlen(c_str_ftag)<1!\n");
     return -1;
   };

   if (isize<1) 
   {
     fprintf(stderr,"err: xml_parse_fgets(): isize<1!\n");
     return -1;
   };
   
   c_str_buf = new (char [isize + strlen(c_str_ftag) + 1]);
   
   ires = xml_parse_fseektag(file, c_str_stag); 
   if (ires<0) return (ires);

   ic = 0; ic1 = 0; ilen = strlen(c_str_ftag);
   do 
   {
      ires = fgetc(file);
      if (ic1>(isize-2)) 
      { 
        fprintf(stderr, "err: xml_parse_fgets(): isize is too small.\n");
        return -1;
      }
      c_str_buf[ic1]=ires; ic1 = ic1 + 1;
      if (c_str_ftag[ic]==ires) {ic=ic+1;} else {ic=0;}
      if (ic==ilen) break;
   } 
   while (ires != EOF);
   
   c_str_buf[ic1] = '\0';
   strncpy(c_str_res, c_str_buf, ic1-strlen(c_str_ftag));
   c_str_res[ic1-strlen(c_str_ftag)] = '\0';
   delete c_str_buf;
   
   if (ic==ilen) {ires=0;}
   
   return ires;
};

//
//
//
/*! \fn int xml_parse_fgeti(int* ires, FILE *file, 
   const char*c_str_stag, const char*c_str_ftag)
   
   read the integer between tags \a c_str_stag, \a c_str_ftag.
   
   \param ires result.
   \param file input stream.
   \param c_str_stag start tag e.g. <html>.
   \param c_str_ftag finish tag e.g. </html>.
      
   \return negative integer if a error occurs, zero otherwise.
   
   \remarks Uses \ref xml_parse_fgetl().
*/

int xml_parse_fgeti(int *ires,  FILE *file, 
  const char*c_str_stag, const char*c_str_ftag)
{
   int ierr;
   long lres;
   
   ierr = xml_parse_fgetl(&lres, file, c_str_stag, c_str_ftag);
   if (ierr<0) return ierr;
   *ires = lres;
   return ierr;
};


//
//
//
/*! \fn int xml_parse_fgetl(long* lres, FILE *file, 
   const char*c_str_stag, const char*c_str_ftag)
   
   read the long integer between tags \a c_str_stag, \a c_str_ftag.
   
   \param lres result.
   \param file input stream.
   \param c_str_stag start tag e.g. <html>.
   \param c_str_ftag finish tag e.g. </html>.
      
   \return negative integer if a error occurs, zero otherwise.
*/

int xml_parse_fgetl(long *lres,  FILE *file, 
  const char*c_str_stag, const char*c_str_ftag)
{
   int ierr;
   long lr;
   char c_str_buf[i_buf_len];

   ierr = xml_parse_fgets(c_str_buf, i_buf_len,  file, c_str_stag, c_str_ftag);
   if (ierr<0)
   {
     fprintf(stderr, 
       "err: fgetl(): something is wrong with the long %s\n", c_str_buf);
     *lres = 0;
     return ierr;
   }
   
   ierr = sscanf(c_str_buf, "%ld", &lr);
   
   if (ierr<0) 
   {
     fprintf(stderr, 
       "err: fgetl(): something is wrong with the long %s\n", c_str_buf);
     *lres = 0;
     return ierr;
   }

   *lres = lr;
   return ierr;
};

//
//
//
/*! \fn int xml_parse_fgetd(long* lres, FILE *file, 
   const char*c_str_stag, const char*c_str_ftag)
   
   read the double between tags \a c_str_stag, \a c_str_ftag.
   
   \param dres result.
   \param file input stream.
   \param c_str_stag start tag e.g. <html>.
   \param c_str_ftag finish tag e.g. </html>.
      
   \return negative integer if a error occurs, zero otherwise.
*/

int xml_parse_fgetd(double *dres,  FILE *file, 
  const char*c_str_stag, const char*c_str_ftag)
{
   int ierr;
   char c_str_buf[i_buf_len];

   ierr = xml_parse_fgets(c_str_buf, i_buf_len,  file, c_str_stag, c_str_ftag);
   if (ierr<0)
   {
     fprintf(stderr, 
       "err: xml_parse_fgetd(): wrong double %s\n", c_str_buf);
     *dres = 0.0;
     return ierr;
   }
   
   ierr = sscanf(c_str_buf, "%lf", dres);
   
   if (ierr<0) 
   {
     fprintf(stderr, 
       "err: fgetd(): smth is wrong with double %s\n", c_str_buf);
     *dres = 0.0;
     return ierr;
   }

   return ierr;
};
