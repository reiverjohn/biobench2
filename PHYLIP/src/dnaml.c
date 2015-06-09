#include "phylip.h"

/* version 3.56c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define maxcategs       9
#define maxtrees        10   /* maximum number of user trees for test        */
#define smoothings      7    /* number of passes through smoothing algorithm */
#define iterations      9    /* number of iterates for each branch           */
#define nmlngth         10   /* max. number of characters in species name    */

#define epsilon         0.0001   /* used in getthree, makenewv */

#define ibmpc0          false
#define ansi0           true
#define vt520           false
#define down            2
#define over            60
#define point           "."


typedef enum {  A, C, G, T} base;
typedef double sitelike[(short)T - (short)A + 1];
typedef sitelike *ratelike;
typedef ratelike *phenotype;
typedef Char **sequence;
typedef double contribarr[maxcategs];
typedef Char naym[nmlngth];
typedef short longer[6];

typedef struct node {
  struct node *next, *back;
  boolean tip, iter;
  short number;
  phenotype x;
  naym nayme;
  double v;
  short xcoord, ycoord, ymin, ymax;
} node;

typedef struct tree {
  node **nodep;
  double likelihood;
  node *start;
} tree;

typedef double *lf[maxtrees];

typedef struct valrec {
  double rat, ratxi, ratxv, zz, z1, y1, ww1, zz1, ww2, zz2, z1zz, z1yy, xiz1,
	 xiy1xv, ww1zz1, vv1zz1, ww2zz2, vv2zz2;
} valrec;

extern short categs,endsite,sites,numsp,numsp2,jumb,lengths,njumble,weightsum,
            outgrno;
extern double *probcat;
extern double *rate;
extern contribarr *contribution;
extern short   *category,*weight,*alias,*ally,*location,*aliasweight;
extern double xi,xv,ttratio,freqa,freqc,freqg,freqt,freqr,freqy,freqar,
              freqcy,freqgr,freqty,lambda,fracchange;
extern short   *enterorder;
extern boolean auto_,usertree,global,progress,treeprint,outgropt,trout,
               ctgry,jumble,lngths;
extern lf l0gf;
extern tree curtree,bestree,bestree2;
extern FILE *infile, *outfile, *treefile;
extern longer seed;

boolean lngths;
double *rate,*probcat;
FILE *infile, *outfile, *treefile;
short numsp, numsp1, numsp2, sites, endsite, weightsum, categs, inseed,
	    outgrno, datasets, ith, i, j, l, jumb, njumble=0;
boolean  freqsfrom, global, jumble,  outgropt, weights,
	       trout, usertree, ctgry, auto_, ttr, printdata, progress,
	       treeprint, mulsets, firstset, interleaved, ibmpc, vt52, ansi;
tree curtree, bestree, bestree2;
double xi, xv, ttratio, ttratio0, freqa, freqc, freqg, freqt, freqr,
	      freqy, freqar, freqcy, freqgr, freqty, fracchange, sumrates,
	      lambda;
longer seed;
short *enterorder;
sequence y;
short *category, *weight, *alias, *ally, *location, *aliasweight;

contribarr *contribution;
lf l0gf;

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
          gets(file);
          }
        break;
      case 'w':
        printf("%s: can't write %s\n",application,file);
	file[0] = '\0';
        while (file[0] =='\0'){
          printf("Please enter a new filename>");
          gets(file);
          }
        break;
      }
    }
  }
  *fp=of;
  if (perm != NULL)
    strcpy(perm,file);
}




void uppercase(ch)
Char *ch;
{
  /* convert ch to upper case -- either ASCII or EBCDIC */
   *ch = isupper(*ch) ? *ch : toupper(*ch);
}  /* uppercase */

Local Void getnums()
{
  /* input number of species, number of sites */
  fprintf(outfile, "\n\n");
  fscanf(infile, "%hd%hd", &numsp, &sites);
  if (printdata)
    fprintf(outfile, "%4hd Species, %4hd Sites\n", numsp, sites);
  numsp1 = numsp + 1;
  numsp2 = numsp * 2 - 2;
}  /* getnums */

void getoptions()
{
  /* interactively set options */
  short i, j, inseed0;
  Char ch;
  char line[256];
  char rest[256];
  int scanned;
  boolean done1, done2, didchangecat;
  double probsum;

  fprintf(outfile, "\nNucleic acid sequence Maximum Likelihood");
  fprintf(outfile, " method, version %s\n\n",VERSION);
  putchar('\n');
  auto_ = false;
  ctgry = false;
  didchangecat = false;
  categs = 1;
  freqsfrom = true;
  global = false;
  jumble = false;
  njumble = 1;
  lngths = false;
  lambda = 1.0;
  outgrno = 1;
  outgropt = false;
  trout = true;
  ttratio = 2.0;
  ttr = false;
  usertree = false;
  weights = false;
  printdata = false;
  progress = true;
  treeprint = true;
  interleaved = true;
  for (;;){
    printf((ansi) ? "\033[2J\033[H" :
	   (vt52) ? "\033E\033H"    : "\n");
    printf("\nNucleic acid sequence Maximum Likelihood");
    printf(" method, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
	   (usertree ? "No, use user trees in input file" : "Yes"));
    if (usertree) {
      printf("  L          Use lengths from user trees?  %s\n",
	     (lngths ? "Yes" : "No"));
    }
    printf("  T        Transition/transversion ratio:%8.4f\n",
	   (ttr ? ttratio : 2.0));
    printf("  F       Use empirical base frequencies?  %s\n",
	   (freqsfrom ? "Yes" : "No"));
    printf("  C   One category of substitution rates?");
    if (!ctgry || categs == 1)
      printf("  Yes\n");
    else {
      printf("  %hd categories\n", categs);
      printf("  R   Rates at adjacent sites correlated?");
      if (!auto_)
	printf("  No, they are independent\n");
      else
	printf("  Yes, mean block length =%6.1f\n", 1.0 / lambda);
    }
    if (!usertree) {
      printf("  G                Global rearrangements?  %s\n",
	     (global ? "Yes" : "No"));
      printf("  J   Randomize input order of sequences?");
      if (jumble)
        printf("  Yes (seed =%8hd,%3hd times)\n", inseed0, njumble);
      else
	printf("  No. Use input order\n");
    }
    printf("  O                        Outgroup root?  %s%3hd\n",
	   (outgropt ? "Yes, at sequence number" :
                       "No, use as outgroup species"),outgrno);
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2hd sets\n", datasets);
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n",
	   (interleaved ? "Yes" : "No, sequential"));
    printf("  0   Terminal type (IBM PC, VT52, ANSI)?  %s\n",
	   (ibmpc ? "IBM PC" :
	    ansi  ? "ANSI"   :
            vt52  ? "VT52"   : "(none)"));
    printf("  1    Print out the data at start of run  %s\n",
	   (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
	   (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
	   (treeprint ? "Yes" : "No"));
    printf("  4       Write out trees onto tree file?  %s\n",
	   (trout ? "Yes" : "No"));
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    uppercase(&ch);
    if (strchr("JOUCRFGLTMI01234",ch) != NULL){
      switch (ch) {
      case 'C':
	ctgry = !ctgry;
	if (!ctgry)
	  auto_ = false;
	if (ctgry) {
	  do {
	    printf("Number of categories (1 - %ld)?\n",(long)maxcategs);
	    gets(line);
	    categs = (short)atoi(line);
	  } while (categs < 1 || categs > maxcategs);
          if (probcat){
            free(probcat);
	    free(rate);
          }
          probcat = (double *)Malloc(categs * sizeof(double));
          rate    = (double *)Malloc(categs * sizeof(double));
          didchangecat = true;
	  for (;;){
            printf("Rate for each category? (use a space to separate)\n");
            gets(line);
            done1 = true;
            for (i = 0; i < categs; i++){
              scanned = sscanf(line,"%lf %[^\n]", &rate[i],rest);
              if ((scanned != 2 && i < (categs - 1)) ||
                 (scanned != 1 && i == (categs - 1))) {
                printf("Please enter exactly %hd values.\n",categs);
                done1 = false;
	        break;
              }
              strcpy(line,rest);
            }
	    if (done1)
	      break;
	  }
        probsum = 0.0;
        for (;;){
          printf("Probability for each category? (use a space to separate)\n");
          gets(line);
          done1 = true;
          probsum = 0.0;
          for (i = 0; i < categs; i++){
            scanned = sscanf(line,"%lf %[^\n]", &probcat[i],rest);
            if ((scanned != 2 && i < (categs - 1)) ||
                (scanned != 1 && i == (categs - 1))){
              printf("Please enter exactly %hd values.\n",categs);
              done1 = false;
              break;
            }
            strcpy(line,rest);
            probsum += probcat[i];
          }
	  if (!done1)
	    continue;
          if (fabs(1.0 - probsum) > 0.001) {
            done1  = false;
	    printf("Probabilities must add up to 1.0, plus or minus 0.001.\n");
	  }
  	  else
	    done1 = true;
  	  if (done1)
	    break;
          }
        }
	break;
	
      case 'R':
	auto_ = !auto_;
	if (auto_) {
	  do {
	    printf(
		"Mean block length of sites having the same rate (greater than 1)?\n");
	    scanf("%lf%*[^\n]", &lambda);
	    getchar();
	    } while (lambda <= 1.0);
	  lambda = 1.0 / lambda;
	}
	break;
	
      case 'F':
	freqsfrom = !freqsfrom;
	if (!freqsfrom) {
	  printf("Base frequencies for A, C, G, T/U (use blanks to separate)?\n");
	  scanf("%lf%lf%lf%lf%*[^\n]", &freqa, &freqc, &freqg, &freqt);
	  getchar();
	}
	break;
	
      case 'G':
	global = !global;
	break;
	
      case 'J':
	jumble = !jumble;
	if (jumble) {
	  printf("Random number seed (must be odd)?\n");
	  scanf("%hd%*[^\n]", &inseed);
	  getchar();
	  inseed0 = inseed;
	  for (i = 0; i <= 5; i++)
	    seed[i] = 0;
	  i = 0;
	  do {
	    seed[i] = inseed & 63;
	    inseed /= 64;
	    i++;
	  } while (inseed != 0);
	  printf("Number of times to jumble?\n");
	  scanf("%hd%*[^\n]", &njumble);
	  getchar();
	}
	else njumble = 1;
	break;
	
      case 'L':
	lngths = !lngths;
	break;
	
      case 'O':
	outgropt = !outgropt;
	if (outgropt) {
	  for (;;) {
	    printf("Type number of the outgroup:\n");
	    scanf("%hd%*[^\n]", &outgrno);
	    getchar();
	    if (outgrno >= 1 || outgrno <= numsp)
	      break;
	    else {
	      printf("BAD OUTGROUP NUMBER: %4hd\n", outgrno);
	      printf("  Must be in range 1 -%2hd\n", numsp);
	    }
	  }
	}
	break;
	
      case 'T':
	ttr = !ttr;
	if (ttr) {
	  do {
	    printf("Transition/transversion ratio?\n");
	    scanf("%lf%*[^\n]", &ttratio);
	    getchar();
	  } while (ttratio < 0.0);
	}
	break;
	
      case 'U':
	usertree = !usertree;
	break;
	
      case 'M':
	mulsets = !mulsets;
	if (mulsets) {
	  for(;;) {
	    printf("How many data sets?\n");
	    scanf("%hd%*[^\n]", &datasets);
	    getchar();
	    if (datasets >= 1)
	      break;
	    else
	      printf("BAD DATA SETS NUMBER:  it must be greater than 1\n");
	  }
	}
	break;
	
      case 'I':
	interleaved = !interleaved;
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
	printdata = !printdata;
	break;
	
      case '2':
	progress = !progress;
	break;
	
      case '3':
	treeprint = !treeprint;
	break;
	
      case '4':
	trout = !trout;
	break;
      }
    } else
      printf("Not a possible option!\n");
  }
if (!didchangecat){
     rate       = Malloc(categs*sizeof(double));
     probcat    = Malloc(categs*sizeof(double));
     rate[0]    = 1.0;
     probcat[0] = 1.0;
   }
}  /* getoptions */


void doinit()
{
  /* initializes variables */
  short i, j;
  node *p, *q;

  getnums();
  getoptions();
  y = (Char **)Malloc(numsp*sizeof(Char *));
  for (i = 0; i < numsp; i++)
    y[i] = (Char *)Malloc(sites*sizeof(Char));
  curtree.nodep = (node **)Malloc(numsp2*sizeof(node *));
  for (i = 0; i < numsp; i++)
    curtree.nodep[i] = (node *)Malloc(sizeof(node));
  for (i = numsp1 - 1; i < numsp2; i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc(sizeof(node));
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    curtree.nodep[i] = p;
  }
  if (usertree)
    return;
  bestree.nodep = (node **)Malloc(numsp2*sizeof(node *));
  for (i = 0; i < numsp; i++)
    bestree.nodep[i] = (node *)Malloc(sizeof(node));
  for (i = numsp1 - 1; i < numsp2; i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc(sizeof(node));
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    bestree.nodep[i] = p;
  }
  if (njumble <= 1)
    return;
  bestree2.nodep = (node **)Malloc(numsp2*sizeof(node *));
  for (i = 0; i < numsp; i++)
    bestree2.nodep[i] = (node *)Malloc(sizeof(node));
  for (i = numsp1 - 1; i < numsp2; i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc(sizeof(node));
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    bestree2.nodep[i] = p;
  }
}  /* doinit */



void inputweights()
{
  /* input the character weights, 0 or 1 */
  Char ch;
  short i;

  for (i = 1; i < nmlngth; i++) {
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
  }
  weightsum = 0;
  for (i = 0; i < sites; i++) {
    do {
      if (eoln(infile)) {
	fscanf(infile, "%*[^\n]");
	getc(infile);
      }
      ch = getc(infile);
      if (ch == '\n')
	ch = ' ';
    } while (ch == ' ');
    weight[i] = 1;
    if (ch == '0' || ch == '1')
      weight[i] = ch - '0';
    else {
      printf("BAD WEIGHT CHARACTER: %c -- WEIGHTS IN DNAML MUST BE 0 OR 1\n",
	     ch);
      exit(-1);
    }
    weightsum += weight[i];
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  weights = true;
}  /* inputweights */

void printweights()
{
  /* print out the weights of sites */
  short i, j;

  fprintf(outfile, "\n   Sites are weighted as follows:\n");
  for (i = 1; i <= sites; i++) {
    if ((i - 1) % 60 == 0) {
      putc('\n', outfile);
      for (j = 1; j <= nmlngth + 3; j++)
	putc(' ', outfile);
    }
    fprintf(outfile, "%hd", weight[i - 1]);
    if (i % 10 == 0 && i % 60 != 0)
      putc(' ', outfile);
  }
  fprintf(outfile, "\n\n");
}  /* printweights */

void inputoptions()
{
  Char ch;
  short i, extranum, cursp, curst;

  if (!firstset) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%hd%hd", &cursp, &curst);
    if (cursp != numsp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4hd\n",
	     ith);
      exit(-1);
    }
    sites = curst;
  }
  for (i = 0; i < sites; i++)
    category[i] = 1;
  for (i = 0; i < sites; i++)
    weight[i] = 1;
  weightsum = sites;
  extranum = 0;
  while (!(eoln(infile))) {
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'C' || ch == 'W')
      extranum++;
    else if (ch != ' ') {
      printf("BAD OPTION CHARACTER: %c\n", ch);
      exit(-1);
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  for (i = 1; i <= extranum; i++) {
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch != 'W'){
      printf("ERROR: INCORRECT AUXILIARY OPTIONS LINE WHICH STARTS WITH %c\n",
	     ch);
      exit(-1);
      }
    else
      inputweights();
  }
  if (categs > 1) {
    fprintf(outfile, "\nSite category   Rate of change    Probability\n\n");
    for (i = 1; i <= categs; i++)
      fprintf(outfile, "%12hd%13.3f%17.3f\n", i, rate[i - 1], probcat[i-1]);
    putc('\n', outfile);
    if (auto_)
      fprintf(outfile,
     "\nExpected length of a patch of sites having the same rate = %8.3f\n",
             1/lambda);
  }
  if (weights && printdata)
    printweights();
}  /* inputoptions */

void getbasefreqs()
{
  double aa, bb;

  putc('\n', outfile);
  if (freqsfrom)
    fprintf(outfile, "Empirical ");
  fprintf(outfile, "Base Frequencies:\n\n");
  fprintf(outfile, "   A    %10.5f\n", freqa);
  fprintf(outfile, "   C    %10.5f\n", freqc);
  fprintf(outfile, "   G    %10.5f\n", freqg);
  fprintf(outfile, "  T(U)  %10.5f\n", freqt);
  freqr = freqa + freqg;
  freqy = freqc + freqt;
  freqar = freqa / freqr;
  freqcy = freqc / freqy;
  freqgr = freqg / freqr;
  freqty = freqt / freqy;
  fprintf(outfile, "\nTransition/transversion ratio = %10.6f\n\n", ttratio);
  aa = ttratio * freqr * freqy - freqa * freqg - freqc * freqt;
  bb = freqa * freqgr + freqc * freqty;
  xi = aa / (aa + bb);
  xv = 1.0 - xi;
  ttratio = xi / xv;
  if (xi <= 0.0) {
    printf("WARNING: This transition/transversion ratio\n");
    printf("is impossible with these base frequencies!\n");
    xi = 3.0 / 5;
    xv = 2.0 / 5;
    fprintf(outfile, " Transition/transversion parameter reset\n\n");
  }
  fprintf(outfile, "(Transition/transversion parameter = %10.6f)\n\n",
	  xi / xv);
  fracchange = xi * (2 * freqa * freqgr + 2 * freqc * freqty) +
      xv * (1.0 - freqa * freqa - freqc * freqc - freqg * freqg - freqt * freqt);
}  /* getbasefreqs */

void setuptree(a)
tree *a;
{
  short i, j;
  node *p;

  for (i = 1; i <= numsp; i++) {
    a->nodep[i - 1]->tip = true;
    a->nodep[i - 1]->iter = true;
    a->nodep[i - 1]->number = i;
  }
  for (i = numsp1; i <= numsp2; i++) {
    p = a->nodep[i - 1];
    for (j = 1; j <= 3; j++) {
      p->tip = false;
      p->iter = true;
      p->number = i;
      p = p->next;
    }
  }
  a->likelihood = -999999.0;
  a->start = a->nodep[0];
}  /* setuptree */

void getdata()
{
  /* read sequences */
  short i, j, k, l, basesread, basesnew;
  Char ch;
  boolean allread, done;

  putc('\n', outfile);
  j = nmlngth + (sites + (sites - 1) / 10) / 2 - 5;
  if (j < nmlngth - 1)
    j = nmlngth - 1;
  if (j > 37)
    j = 37;
  if (printdata) {
    fprintf(outfile, "Name");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "Sequences\n");
    fprintf(outfile, "----");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "---------\n\n");
  }
  setuptree(&curtree);
  basesread = 0;
  allread = false;
  while (!(allread)) {
    allread = true;
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    i = 1;
    while (i <= numsp) {
      if ((interleaved && basesread == 0) || !interleaved) {
	for (j = 0; j < nmlngth; j++) {
	  if (eof(infile) || eoln(infile)){
	    printf("ERROR: END-OF-LINE OR END-OF-FILE");
            printf(" IN THE MIDDLE OF A SPECIES NAME\n");
	    exit(-1);
            }
	  curtree.nodep[i - 1]->nayme[j] = getc(infile);
	}
      }
      if (interleaved)
	j = basesread;
      else
	j = 0;
      done = false;
      while (!done & !eof(infile)) {
	if (interleaved)
	  done = true;
	while (((j < sites) & (!(eoln(infile) | eof(infile))))) {
	  ch = getc(infile);
	  if (ch == '\n')
	    ch = ' ';
	  if (ch == ' ' || (ch >= '0' && ch <= '9'))
	    continue;
	  uppercase(&ch);
	  if (strchr("ABCDGHKMNRSTUVWXY?O-.",ch) == NULL){
	    printf("ERROR: BAD BASE:%c AT POSITION%5hd OF SPECIES %3ld\n",
		   ch, j, i);
	    exit(-1);
	  }
	  j++;
	  if (ch == '.')
	    ch = y[0][j - 1];
	  y[i - 1][j - 1] = ch;
	}
	if (interleaved)
	  continue;
	if (j < sites) {
	  fscanf(infile, "%*[^\n]");
	  getc(infile);
	} else if (j == sites)
	  done = true;
      }
      if (interleaved && i == 1)
	basesnew = j;
      fscanf(infile, "%*[^\n]");
      getc(infile);
      if ((interleaved && j != basesnew) || ((!interleaved) && j != sites)){
	printf("ERROR: SEQUENCES OUT OF ALIGNMENT\n");
	exit(-1);
        }
     i++;
    }
    if (interleaved) {
      basesread = basesnew;
      allread = (basesread == sites);
    } else
      allread = (i > numsp);
  }
  if (printdata) {
    for (i = 1; i <= ((sites - 1) / 60 + 1); i++) {
      for (j = 1; j <= numsp; j++) {
	for (k = 0; k < nmlngth; k++)
	  putc(curtree.nodep[j - 1]->nayme[k], outfile);
	fprintf(outfile, "   ");
	l = i * 60;
	if (l > sites)
	  l = sites;
	for (k = (i - 1) * 60 + 1; k <= l; k++) {
	  if (j > 1 && y[j - 1][k - 1] == y[0][k - 1])
	    ch = '.';
	  else
	    ch = y[j - 1][k - 1];
	  putc(ch, outfile);
	  if (k % 10 == 0 && k % 60 != 0)
	    putc(' ', outfile);
	}
	putc('\n', outfile);
      }
      putc('\n', outfile);
    }
  }
  if (!usertree) {
    setuptree(&bestree);
    if (njumble > 1)
      setuptree(&bestree2);
  }
}  /* getdata */

void sitesort()
{
  /* Shell sort keeping sites, weights in same order */
  short gap, i, j, jj, jg, k, itemp;
  boolean flip, tied;

  gap = sites / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= sites; i++) {
      j = i - gap;
      flip = true;
      while (j > 0 && flip) {
	jj = alias[j - 1];
	jg = alias[j + gap - 1];
        tied = (weight[jj - 1] == weight[jg - 1]);
        flip = (weight[jj - 1] < weight[jg - 1]);
	k = 1;
	while (k <= numsp && tied) {
	  flip = (y[k - 1][jj - 1] > y[k - 1][jg - 1]);
	  tied = (tied && y[k - 1][jj - 1] == y[k - 1][jg - 1]);
	  k++;
	}
	if (!flip)
	  break;
	itemp = alias[j - 1];
	alias[j - 1] = alias[j + gap - 1];
	alias[j + gap - 1] = itemp;
	itemp = aliasweight[j - 1];
	aliasweight[j - 1] = aliasweight[j + gap - 1];
	aliasweight[j + gap - 1] = itemp;
	j -= gap;
      }
    }
    gap /= 2;
  }
}  /* sitesort */

void sitecombine()
{
  /* combine sites that have identical patterns */
  short i, j, k;
  boolean tied;

  i = 1;
  while (i < sites) {
    j = i + 1;
    tied = true;
    while (j <= sites && tied) {
      tied = (weight[alias[i - 1] - 1] == weight[alias[j - 1] - 1]);
      k = 1;
      while (k <= numsp && tied) {
	tied = (tied &&
	    y[k - 1][alias[i - 1] - 1] == y[k - 1][alias[j - 1] - 1]);
	k++;
      }
      if (!tied)
	break;
      aliasweight[i - 1] += aliasweight[j - 1];
      aliasweight[j - 1] = 0;
      ally[alias[j - 1] - 1] = alias[i - 1];
      j++;
    }
    i = j;
  }
}  /* sitecombine */

void sitescrunch()
{
  /* move so positively weighted sites come first */
  short i, j, itemp;
  boolean done, found;

  done = false;
  i = 1;
  j = 2;
  while (!done) {
    found = false;
    if (aliasweight[i - 1] > 0)
      i++;
    else {
      if (j <= i)
	j = i + 1;
      if (j <= sites) {
	found = false;
	do {
	  found = (aliasweight[j - 1] > 0);
	  j++;
	} while (!(found || j > sites));
	if (found) {
	  j--;
	  itemp = alias[i - 1];
	  alias[i - 1] = alias[j - 1];
	  alias[j - 1] = itemp;
	  itemp = aliasweight[i - 1];
	  aliasweight[i - 1] = aliasweight[j - 1];
	  aliasweight[j - 1] = itemp;
	} else
	  done = true;
      } else
	done = true;
    }
    done = (done || i >= sites);
  }
}  /* sitescrunch */

void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  int i,j,k;
  node *p;
  double temp;

  for (i = 1; i <= sites; i++) {
    alias[i - 1] = i;
    ally[i - 1] = 0;
    aliasweight[i - 1] = weight[i - 1];
    location[i - 1] = 0;
  }
  sitesort();
  sitecombine();
  sitescrunch();
  for (i = 1; i <= sites; i++) {
    weight[i - 1] = aliasweight[i - 1];
    if (weight[i - 1] > 0)
      endsite = i;
  }
  for (i = 1; i <= endsite; i++) {
    location[alias[i - 1] - 1] = i;
    ally[alias[i - 1] - 1] = alias[i - 1];
  }
  sumrates = 0.0;
  for (i = 0; i < categs; i++)
    sumrates += probcat[i] * rate[i];
  for (i = 0; i < categs; i++)
    rate[i] /= sumrates;
  if (usertree)
    for (i = 0; i < maxtrees; i++)
      l0gf[i] = (double *)Malloc(endsite*sizeof(double));
  contribution = (contribarr *)Malloc(endsite*sizeof(contribarr));
  for (i = 0; i < numsp; i++){
    curtree.nodep[i]->x = (phenotype)Malloc(endsite*sizeof(ratelike));
    for (j=0;j<endsite;j++)
      curtree.nodep[i]->x[j]  = (ratelike)Malloc(categs*sizeof(sitelike));
  }
  for (i = numsp1 - 1; i < numsp2; i++) {
    p = curtree.nodep[i];
    for (j = 1; j <= 3; j++) {
      p->x = (phenotype)Malloc(endsite*sizeof(ratelike));
      for (k=0;k<endsite;++k)
           p->x[k] = (ratelike)Malloc(categs*sizeof(sitelike));
      p = p->next;
    }
  }
  if (usertree)
    return;
  for (i = 0; i < numsp; i++){
    bestree.nodep[i]->x = (phenotype)Malloc(endsite*sizeof(ratelike));
    for (j=0;j<endsite;++j)
        bestree.nodep[i]->x[j] = (ratelike)Malloc(categs*sizeof(sitelike));
      }

  for (i = numsp1 - 1; i < numsp2; i++) {
    p = bestree.nodep[i];
    for (j = 1; j <= 3; j++) {
      p->x = (phenotype)Malloc(endsite*sizeof(ratelike));
      for (k=0;k<endsite;++k)
        p->x[k] = (ratelike)Malloc(categs*sizeof(sitelike));
      p = p->next;
    }
  }
  if (njumble <= 1)
    return;
  for (i = 0; i < numsp; i++) {
    bestree2.nodep[i]->x = (phenotype)Malloc(endsite*sizeof(ratelike));
    for (j=0;j<endsite;++j)
        bestree2.nodep[i]->x[j] = (ratelike)Malloc(categs*sizeof(sitelike));
  }
  for (i = numsp1 - 1; i < numsp2; i++) {
    p = bestree2.nodep[i];
    for (j = 1; j <= 3; j++) {
      p->x = (phenotype)Malloc(endsite*sizeof(ratelike));
      for (k=0;k<endsite;++k)
        p->x[k] = (ratelike)Malloc(categs*sizeof(sitelike));
      p = p->next;
    }
  }
}  /* makeweights */

void makevalues()
{
  /* set up fractional likelihoods at tips */
  short i, j, k, l;
  base b;

  for (k = 0; k < endsite; k++) {
    j = alias[k];
    for (i = 0; i < numsp; i++) {
      for (l = 0; l < categs; l++) {
	for (b = A; (short)b <= (short)T; b = (base)((short)b + 1))
	  curtree.nodep[i]->x[k][l][(short)b - (short)A] = 0.0;
	switch (y[i][j - 1]) {

	case 'A':
	  curtree.nodep[i]->x[k][l][0] = 1.0;
	  break;

	case 'C':
	  curtree.nodep[i]->x[k][l][(short)C - (short)A] = 1.0;
	  break;

	case 'G':
	  curtree.nodep[i]->x[k][l][(short)G - (short)A] = 1.0;
	  break;

	case 'T':
	  curtree.nodep[i]->x[k][l][(short)T - (short)A] = 1.0;
	  break;

	case 'U':
	  curtree.nodep[i]->x[k][l][(short)T - (short)A] = 1.0;
	  break;

	case 'M':
	  curtree.nodep[i]->x[k][l][0] = 1.0;
	  curtree.nodep[i]->x[k][l][(short)C - (short)A] = 1.0;
	  break;

	case 'R':
	  curtree.nodep[i]->x[k][l][0] = 1.0;
	  curtree.nodep[i]->x[k][l][(short)G - (short)A] = 1.0;
	  break;

	case 'W':
	  curtree.nodep[i]->x[k][l][0] = 1.0;
	  curtree.nodep[i]->x[k][l][(short)T - (short)A] = 1.0;
	  break;

	case 'S':
	  curtree.nodep[i]->x[k][l][(short)C - (short)A] = 1.0;
	  curtree.nodep[i]->x[k][l][(short)G - (short)A] = 1.0;
	  break;

	case 'Y':
	  curtree.nodep[i]->x[k][l][(short)C - (short)A] = 1.0;
	  curtree.nodep[i]->x[k][l][(short)T - (short)A] = 1.0;
	  break;

	case 'K':
	  curtree.nodep[i]->x[k][l][(short)G - (short)A] = 1.0;
	  curtree.nodep[i]->x[k][l][(short)T - (short)A] = 1.0;
	  break;

	case 'B':
	  curtree.nodep[i]->x[k][l][(short)C - (short)A] = 1.0;
	  curtree.nodep[i]->x[k][l][(short)G - (short)A] = 1.0;
	  curtree.nodep[i]->x[k][l][(short)T - (short)A] = 1.0;
	  break;

	case 'D':
	  curtree.nodep[i]->x[k][l][0] = 1.0;
	  curtree.nodep[i]->x[k][l][(short)G - (short)A] = 1.0;
	  curtree.nodep[i]->x[k][l][(short)T - (short)A] = 1.0;
	  break;

	case 'H':
	  curtree.nodep[i]->x[k][l][0] = 1.0;
	  curtree.nodep[i]->x[k][l][(short)C - (short)A] = 1.0;
	  curtree.nodep[i]->x[k][l][(short)T - (short)A] = 1.0;
	  break;

	case 'V':
	  curtree.nodep[i]->x[k][l][0] = 1.0;
	  curtree.nodep[i]->x[k][l][(short)C - (short)A] = 1.0;
	  curtree.nodep[i]->x[k][l][(short)G - (short)A] = 1.0;
	  break;

	case 'N':
	  for (b = A; (short)b <= (short)T; b = (base)((short)b + 1))
	    curtree.nodep[i]->x[k][l][(short)b - (short)A] = 1.0;
	  break;

	case 'X':
	  for (b = A; (short)b <= (short)T; b = (base)((short)b + 1))
	    curtree.nodep[i]->x[k][l][(short)b - (short)A] = 1.0;
	  break;

	case '?':
	  for (b = A; (short)b <= (short)T; b = (base)((short)b + 1))
	    curtree.nodep[i]->x[k][l][(short)b - (short)A] = 1.0;
	  break;

	case 'O':
	  for (b = A; (short)b <= (short)T; b = (base)((short)b + 1))
	    curtree.nodep[i]->x[k][l][(short)b - (short)A] = 1.0;
	  break;

	case '-':
	  for (b = A; (short)b <= (short)T; b = (base)((short)b + 1))
	    curtree.nodep[i]->x[k][l][(short)b - (short)A] = 1.0;
	  break;
	}
      }
    }
  }
}  /* makevalues */

Local Void empiricalfreqs()
{
  /* Get empirical base frequencies from the data */
  short i, j, k;
  double sum, suma, sumc, sumg, sumt, w;

  freqa = 0.25;
  freqc = 0.25;
  freqg = 0.25;
  freqt = 0.25;
  for (k = 1; k <= 8; k++) {
    suma = 0.0;
    sumc = 0.0;
    sumg = 0.0;
    sumt = 0.0;
    for (i = 0; i < numsp; i++) {
      for (j = 0; j < endsite; j++) {
	w = weight[j];
	sum = freqa * curtree.nodep[i]->x[j][0][0];
	sum += freqc * curtree.nodep[i]->x[j][0][(short)C - (short)A];
	sum += freqg * curtree.nodep[i]->x[j][0][(short)G - (short)A];
	sum += freqt * curtree.nodep[i]->x[j][0][(short)T - (short)A];
	suma += w * freqa * curtree.nodep[i]->x[j][0][0] / sum;
	sumc += w * freqc * curtree.nodep[i]->x[j][0][(short)C - (short)A] / sum;
	sumg += w * freqg * curtree.nodep[i]->x[j][0][(short)G - (short)A] / sum;
	sumt += w * freqt * curtree.nodep[i]->x[j][0][(short)T - (short)A] / sum;
      }
    }
    sum = suma + sumc + sumg + sumt;
    freqa = suma / sum;
    freqc = sumc / sum;
    freqg = sumg / sum;
    freqt = sumt / sum;
  }
}  /* empiricalfreqs */


void getinput()
{
  /* reads the input data */
  inputoptions();
  if (!freqsfrom)
    getbasefreqs();
  getdata();
  makeweights();
  makevalues();
  if (freqsfrom) {
    empiricalfreqs();
    getbasefreqs();
  }
}  /* getinput */


main(argc, argv)
int argc;
Char *argv[];
{  /* DNA Maximum Likelihood */
char infilename[100],outfilename[100],trfilename[100];
#ifdef MAC
  macsetup("Dnaml","");
  argv[0] = "Dnaml";
#endif
  openfile(&infile,"infile","r",argv[0],infilename);
  openfile(&outfile,"outfile","w",argv[0],outfilename);
  mulsets = false;
  datasets = 1;
  firstset = true;
  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  doinit();
  ttratio0    = ttratio;
  enterorder  = (short *)Malloc(numsp*sizeof(short));
  category    = (short *)Malloc(sites*sizeof(short));
  weight      = (short *)Malloc(sites*sizeof(short));
  alias       = (short *)Malloc(sites*sizeof(short));
  ally        = (short *)Malloc(sites*sizeof(short));
  location    = (short *)Malloc(sites*sizeof(short));
  aliasweight = (short *)Malloc(sites*sizeof(short));
  if (trout)
    openfile(&treefile,"treefile","w",argv[0],trfilename);
  for (ith = 1; ith <= datasets; ith++) {
    if (datasets > 1) {
      fprintf(outfile, "Data set # %hd:\n", ith);
      printf("\nData set # %hd:\n", ith);
    }
    ttratio = ttratio0;
    getinput();
    if (ith == 1)
      firstset = false;
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
  }
  FClose(infile);
  FClose(outfile);
  FClose(treefile);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(trfilename);
#endif
  exit(0);
}  /* DNA Maximum Likelihood */

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
long  x;
{
MALLOCRETURN *mem;
mem = (MALLOCRETURN *)calloc(1,x);
if (!mem)
  memerror();
else
  return (MALLOCRETURN *)mem;
}


