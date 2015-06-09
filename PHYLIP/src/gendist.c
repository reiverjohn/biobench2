#include "phylip.h"

/* version 3.56c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define namelength      10   /* number of characters max. in species name */
#define epsilon         0.02 /* a small number                            */

#define ibmpc0          false
#define ansi0           true
#define vt520           false


typedef double *phenotype;
typedef Char naym[namelength];

Static FILE *infile, *outfile;
Static short numsp, loci, totalleles, df, datasets, ith;
Static short *alleles;
Static phenotype *x;
Static double **d;
Static naym *nayms;
Static boolean all, cavalli, lower, nei, reynolds,  mulsets, ibmpc,
	       vt52, ansi, firstset, progress;


openfile(fp,filename,mode,application,perm)
FILE **fp;
char *filename;
char *mode;
char *application;
char *perm;
{
  FILE *of;
  char file[100];
  strcpy(file,filename);
  while (1){
    of = fopen(file,mode);
    if (of)
      break;
    else {
      switch (*mode){
      case 'r':
        printf("%s:  can't read %s\n",application,file);
	file[0] = '\0';
        while (file[0] =='\0'){
          printf("Please enter a new filename>");
          gets(file);}
        break;
      case 'w':
        printf("%s: can't write %s\n",application,file);
	file[0] = '\0';
        while (file[0] =='\0'){
          printf("Please enter a new filename>");
          gets(file);}
        break;
      }
    }
  }
  *fp=of;
  if (perm != NULL)
    strcpy(perm,file);
}


Static Void uppercase(ch)
Char *ch;
{  /* convert a character to upper case -- either ASCII or EBCDIC */
   *ch = (islower(*ch) ?  toupper(*ch) : (*ch));
}  /* uppercase */


void getnums()
{
  /* read number of species and loci for first data set */
  short i;
  fscanf(infile, "%hd%hd", &numsp, &loci);
}  /* getnums */

void getoptions()
{
  /* interactively set options */
  Char ch;
  boolean  done1;

  all = false;
  cavalli = false;
  lower = false;
  nei = true;
  reynolds = false;
  lower = false;
  progress = true;
  for (;;) {
    printf(ansi ? "\033[2J\033[H" :
	   vt52 ? "\033E\033H"    : "\n");
    printf("\nGenetic Distance Matrix program, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  A   Input file contains all alleles at each locus?  %s\n",
	   all ? "Yes" : "One omitted at each locus");
    printf("  N                        Use Nei genetic distance?  %s\n",
	   nei ? "Yes" : "No");
    printf("  C                Use Cavalli-Sforza chord measure?  %s\n",
	   cavalli ? "Yes" : "No");
    printf("  R                   Use Reynolds genetic distance?  %s\n",
	   reynolds ? "Yes" : "No");
    printf("  L                         Form of distance matrix?  %s\n",
	   lower ? "Lower-triangular" : "Square");
    printf("  M                      Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2hd sets\n", datasets);
    else
      printf("  No\n");
    printf("  0              Terminal type (IBM PC, VT52, ANSI)?  %s\n",
	   ibmpc ? "IBM PC" :
	   ansi  ? "ANSI"   :
	   vt52  ? "VT52"   : "(none)");
    printf("  1            Print indications of progress of run?  %s\n",
	   progress ? "Yes" : "No");
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (ch == 'A' || ch == 'C' || ch == 'N' || ch == 'M' || ch == 'R' ||
	ch == 'L' || ch == '0' || ch == '1') {
      switch (ch) {
	
      case 'A':
	all = !all;
	break;
	
      case 'C':
	cavalli = true;
	nei = false;
	reynolds = false;
	break;
	
      case 'N':
	cavalli = false;
	nei = true;
	reynolds = false;
	break;
	
      case 'R':
	reynolds = true;
	cavalli = false;
	nei = false;
	break;
	
      case 'L':
	lower = !lower;
	break;
	
      case 'M':
	mulsets = !mulsets;
	if (mulsets) {
	  done1 = false;
	  do {
	    printf("How many data sets?\n");
	    scanf("%hd%*[^\n]", &datasets);
	    getchar();
	    done1 = (datasets >= 1);
	    if (!done1)
	      printf("BAD DATA SETS NUMBER:  it must be greater than 1\n");
	  } while (done1 != true);
	}
	break;
	
      case '0':
	if (ibmpc) {
	  ibmpc = false;
	  vt52 = true;
	} else {
	  if (vt52) {
	    vt52 = false;
	    ansi = true;
	  } else if (ansi)
	    ansi = false;
	  else
	    ibmpc = true;
	}
	break;
	
      case '1':
	progress = !progress;
	break;
      }
    } else
      printf("Not a possible option!\n");
  }
  putchar('\n');
}  /* getoptions */


void doinit()
{
  /* initializes variables */
  short i;

  getnums();
  x = (phenotype *)Malloc(numsp*sizeof(phenotype));
  d = (double **)Malloc(numsp*sizeof(double *));
  for (i = 0; i < (numsp); i++)
    d[i] = (double *)Malloc(numsp*sizeof(double));
  alleles = (short *)Malloc(loci*sizeof(short));
  getoptions();
}  /* doinit */


void getalleles()
{
  short i, cursp, curloc;

  if (!firstset) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%hd%hd", &cursp, &curloc);
    if (cursp != numsp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4hd\n", ith);
      exit(-1);
    }
    loci = curloc;
  }
  totalleles = 0;
  fscanf(infile, "%*[^\n]");
  getc(infile);
  for (i = 0; i < (loci); i++) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%hd", &alleles[i]);
    totalleles += alleles[i];
  }
  df = totalleles - loci;
}  /* getalleles */

void getdata()
{
  /* read allele frequencies */
  short i, j, k, m, n, p;
  double sum;

  for (i = 0; i < numsp; i++)
    x[i] = (phenotype)Malloc(totalleles*sizeof(double));
  for (i = 1; i <= (numsp); i++) {
    fscanf(infile, "%*[^\n]");
    getc(infile);
    for (j = 0; j < namelength; j++)
      nayms[i - 1][j] = getc(infile);
    m = 1;
    p = 1;
    for (j = 1; j <= (loci); j++) {
      sum = 0.0;
      if (all)
	n = alleles[j - 1];
      else
	n = alleles[j - 1] - 1;
      for (k = 1; k <= n; k++) {
	if (eoln(infile)) {
	  fscanf(infile, "%*[^\n]");
	  getc(infile);
	}
	fscanf(infile, "%lf", &x[i - 1][m - 1]);
	sum += x[i - 1][m - 1];
	if (x[i - 1][m - 1] < 0.0) {
	  printf("\nLOCUS%3hd IN SPECIES%3hd: AN ALLELE", j, i);
	  printf(" FREQUENCY IS NEGATIVE\n");
	  exit(-1);
	}
	p++;
	m++;
      }
      if (all && fabs(sum - 1.0) > epsilon) {
	printf("\nLOCUS%3hd IN SPECIES%3hd: FREQUENCIES DO NOT ADD UP TO 1\n",
	       j, i);
	exit(-1);
      }
      if (!all) {
	x[i - 1][m - 1] = 1.0 - sum;
	if (x[i - 1][m - 1] < -epsilon) {
	  printf("\nLOCUS%3hd IN SPECIES%3hd: ",j,i);
	  printf("FREQUENCIES ADD UP TO MORE THAN 1\n");
	  exit(-1);
	}
	m++;
      }
    }
  }
}  /* getdata */


void getinput()
{
  /* read the input data */
  getalleles();
  getdata();
}  /* getinput */


void makedists()
{
  short i, j, k;
  double s, s1, s2, s3, f;
  double TEMP;

  for (i = 0; i < (numsp); i++)
    d[i][i] = 0.0;
  for (i = 1; i <= (numsp); i++) {
    for (j = 0; j <= i - 2; j++) {
      if (cavalli) {
	s = 0.0;
	for (k = 0; k < (totalleles); k++) {
	  f = x[i - 1][k] * x[j][k];
	  if (f > 0.0)
	    s += sqrt(f);
	}
	d[i - 1][j] = 4 * (loci - s) / df;
      }
      if (nei) {
	s1 = 0.0;
	s2 = 0.0;
	s3 = 0.0;
	for (k = 0; k < (totalleles); k++) {
	  s1 += x[i - 1][k] * x[j][k];
	  TEMP = x[i - 1][k];
	  s2 += TEMP * TEMP;
	  TEMP = x[j][k];
	  s3 += TEMP * TEMP;
	}
	if (s1 <= 1.0e-20)
	  d[i - 1][j] = -1.0;
	else
	  d[i - 1][j] = -log(s1 / sqrt(s2 * s3));
      }
      if (reynolds) {
	s1 = 0.0;
	s2 = 0.0;
	for (k = 0; k < (totalleles); k++) {
	  TEMP = x[i - 1][k] - x[j][k];
	  s1 += TEMP * TEMP;
	  s2 += x[i - 1][k] * x[j][k];
	}
	d[i - 1][j] = s1 / (loci * 2 - 2 * s2);
      }
      d[j][i - 1] = d[i - 1][j];
    }
  }
}  /* makedists */

void writedists()
{
  short i, j, k;

  fprintf(outfile, "%5hd\n", numsp);
  for (i = 0; i < (numsp); i++) {
    for (j = 0; j < namelength; j++)
      putc(nayms[i][j], outfile);
    if (lower)
      k = i;
    else
      k = numsp;
    for (j = 1; j <= k; j++) {
      fprintf(outfile, "%8.4f", d[i][j - 1]);
      if ((j + 1) % 9 == 0 && j < k)
	putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  if (progress)
    printf("Distances written to file\n\n");
}  /* writedists */


main(argc, argv)
int argc;
Char *argv[];
{  /* main program */
char infilename[100],outfilename[100];
#ifdef MAC
  macsetup("Gendist","");
  argv[0] = "Gendist";
#endif
  openfile(&infile,INFILE,"r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);

  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  mulsets = false;
  firstset = true;
  datasets = 1;
  doinit();
  nayms = (naym *)Malloc(numsp*sizeof(naym));
  for (ith = 1; ith <= (datasets); ith++) {
    getinput();
    firstset = false;
    if ((datasets > 1) && progress)
      printf("\nData set # %hd:\n\n",ith);
    makedists();
    writedists();
  }
  FClose(infile);
  FClose(outfile);
#ifdef MAC
  fixmacfile(outfilename);
#endif
  exit(0);
}


int eof(f)
FILE *f;
{
    register int ch;

    if (feof(f))
	return 1;
    if (f == stdin)
	return 0;
    ch = getc(f);
    if (ch == EOF)
	return 1;
    ungetc(ch, f);
    return 0;
}


int eoln(f)
FILE *f;
{
    register int ch;

    ch = getc(f);
    if (ch == EOF)
        return 1;
    ungetc(ch, f);
    return (ch == '\n');
}

void memerror()
{
printf("Error allocating memory\n");
exit(-1);
}

MALLOCRETURN *mymalloc(x)
long x;
{
MALLOCRETURN *mem;
mem = (MALLOCRETURN *)malloc(x);
if (!mem)
  memerror();
else
  return (MALLOCRETURN *)mem;
}

