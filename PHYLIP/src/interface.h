/* interface.h:  access to interface.c, a 2 window text/graphics
    environment, with a scrolling text window and c-like I/O functions.
    This also sets up some defines for the standard c stuff. */

#define printf  macprintf
#define gets    macgets
#define scanf   macscanf
#define puts    macputs
#undef putchar
#define putchar macputchar
#undef getchar
#define getchar() '\n'
#define fflush(file) macflush()
#define exit(x) eventloop()
void macprintf(char *,...);
int macscanf(char *,...);
void macsetup(char *,char *);
void macgets(char *);
void macputs(char *);
void eventloop();
void queryevent();
void fixmacfile(char *);

#ifdef THINK_C
#define MAC
#else
#ifdef GENERATINGPOWERPC
#undef METRO
#undef MAC
#define METRO
#define MAC
#else
#ifdef GENERATING68K
#undef METRO
#undef MAC
#define METRO
#define MAC
#endif
#endif
#endif


