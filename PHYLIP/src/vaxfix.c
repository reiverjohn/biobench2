#include<stdio.h>
#include<stdarg.h>
#define printf vax_printf_is_broken

void vax_printf_is_broken(const char *fmt,...);
void vax_fprintf_is_broken(FILE *fp,const char *fmt,...);
void vax_tweak_fmt(char *);


void vax_tweak_fmt(char *format)
{
int inperc,rp,wp;
rp=0;
wp=0;
while (rp <= strlen(format))
  {
   format[rp] = format[wp];
   if (format[rp] == '%')
        inperc = 1;
   if (inperc){
   if ((format[rp] < '0' || format[rp] > '9') && format[rp] != 'h'
                                              && format[rp] != '%')
        inperc=0;
   }
 wp++;
 if (!inperc || (inperc && format[rp] != 'h'))
   rp++;
 }
}

void vax_printf_is_broken(const char *fmt,...)
{
va_list args;
char output[256];
vax_tweak_fmt(fmt);
va_start(args,fmt);
vsprintf(output,fmt,args);
va_end(args);
fputs(output,stdout);
}

void vax_fprintf_is_broken(FILE *fp,const char *fmt,...)
{
va_list args;
char output[256];
vax_tweak_fmt(fmt);
va_start(args,fmt);
vfprintf(fp,fmt,args);
va_end(args);
}
