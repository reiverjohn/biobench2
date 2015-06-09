#define VERSION "3.573c"

/* machine-specific stuff:
   based on a number of factors in the library stdlib.h, we will try
   to determine what kind of machine/compiler this program is being
   built on.  However, it doesn't always succeed.  However, if you have
   ANSI conforming C, it will probably work.

 we will try to figure out machine type
 based on defines in stdio, and compiler-defined things as well.: */

#ifdef  GNUDOS
#define DJGPP
#define DOS
#endif

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

#include <stdio.h>
#include <stdlib.h>

#ifdef __CMS_OPEN
#define CMS
#define EBCDIC true
#define INFILE "infile data"
#define OUTFILE "outfile data"
#define TREEFILE "treefile data"
#define FONTFILE "fontfile data"
#define PLOTFILE "plotfile data"
#define INTREE "intree data"
#define OUTTREE "outtree data"
#else
#define EBCDIC false
#define INFILE "infile"
#define OUTFILE "outfile"
#define TREEFILE "treefile"
#define FONTFILE "fontfile" /* on unix this might be /usr/local/lib/fontfile */
#define PLOTFILE "plotfile"
#define INTREE "intree"
#define OUTTREE "outtree"
#endif

#ifdef L_ctermid            /* try and detect for sysV or V7. */
#define SYSTEM_FIVE
#endif

#ifdef sequent
#define SYSTEM_FIVE
#endif

#ifndef SYSTEM_FIVE
#include <stdlib.h>
#if defined(_STDLIB_H_) || defined(_H_STDLIB) || defined(H_SCCSID) || defined(unix)
#define UNIX
#define MACHINE_TYPE "BSD Unix C"
#endif
#endif


#ifdef __STDIO_LOADED
#define VMS
#define MACHINE_TYPE "VAX/VMS C"
#define printf vax_printf_is_broken
#define fprintf vax_fprintf_is_broken
void vax_printf_is_broken(const char *fmt,...);
void vax_fprintf_is_broken(FILE *fp,const char *fmt,...);
void vax_tweak_fmt(char *);
#endif

#ifdef __WATCOMC__
#define QUICKC
#define WATCOM
#define DOS
#include "graph.h"
#endif
/* watcom-c has graphics library calls that are almost identical to    *
 * quick-c, so the "QUICKC" symbol name stays.                         */


#ifdef _QC
#define MACHINE_TYPE "MS-DOS / Quick C"
#define QUICKC
#include "graph.h"
#define DOS
#endif

#ifdef _DOS_MODE
#define MACHINE_TYPE "MS-DOS /Microsoft C "
#define DOS           /* DOS is  always defined if  on a dos machine */
#define MSC           /* MSC is defined for microsoft C              */
#endif

#ifdef __MSDOS__      /* TURBO c compiler, ONLY (no other DOS C compilers) */
#define DOS
#define TURBOC
#include<stdlib.h>
#include<graphics.h>
#endif

#ifdef DJGPP          /* DJ's gnu  C/C++ port */
#include<graphics.h>
#endif

#ifndef MACHINE_TYPE
#define MACHINE_TYPE "ANSI C"
#endif

#ifdef DOS
#define MALLOCRETURN void 
#else
#define MALLOCRETURN void
#endif

#ifdef VMS
#define signed /* signed doesn't exist in VMS */
#endif

#ifdef DJGPP
#undef MALLOCRETURN
#define MALLOCRETURN void
#endif


/* includes: */
#ifdef UNIX
#include<strings.h>
#else
#include<string.h>
#endif

#include <math.h>
#include <ctype.h>

#ifdef MAC
#include "interface.h"
#endif

#define FClose(file) if (file) fclose(file) ; file=NULL
#define Malloc(x) mymalloc((long)x)

typedef void *Anyptr;
#define Signed     signed
#define Void       void      /* Void f() = procedure */
#define Const     const
#define Volatile  volatile
#define Char        char      /* Characters (not bytes) */
#define Static     static     /* Private global funcs and vars */
#define Local      static     /* Nested functions */

typedef unsigned char boolean;

#define true    1
#define false   0
#define SETBITS 32

#ifdef MAC
MALLOCRETURN    *mymalloc(long);
#else
MALLOCRETURN    *mymalloc();
#endif
int      eof();
int      eoln();
void     memerror();

