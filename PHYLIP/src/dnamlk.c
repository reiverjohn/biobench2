#include "phylip.h"

/* version 3.56c. (c) Copyright 1986-1993 by the University of Washington
  and by Joseph Felsenstein.  Written by Joseph Felsenstein.  Permission is
  granted to copy and use this program provided no fee is charged for it
  and provided that this copyright notice is not removed. */

#define maxcategs       9    /* maximum number of site types                 */
#define smoothings      2    /* number of passes through smoothing algorithm */
#define iterations      10   /* number of iterates for each branch           */
#define nmlngth         10   /* number of characters max. in species name    */

#define epsilon         0.0001   /* used in makenewv, getthree, update */

#define point           '.'

#define ibmpc0          false
#define ansi0           true
#define vt520           false


typedef enum {
  A, C, G, T
} base;
typedef double sitelike[(short)T - (short)A + 1];
typedef sitelike **phenotype;
typedef double *contribarr;
typedef char **sequence;


typedef short longer[6];

typedef struct node {
  struct node *next, *back;
  phenotype x;
  double tyme, v;
  short number, xcoord, ycoord, ymin, ymax;
  boolean tip, haslength, initialized;
} node;


typedef struct tree {
  node **nodep;
  double likelihood;
  node *root;
} tree;


extern short      categ,spp,endsite,sites,nonodes,weightsum,datasets,
                  njumble,jumb;
extern double     rate[maxcategs];
extern double     xi,xv,freqa,freqc,freqg,freqt,freqr,freqy,freqar,freqcy,
                  freqgr,freqty,lambda,lambda1,fracchange;
extern boolean    usertree,globle,jumble,lengths,trout,weights,ctgry,auto_,
                  printdata,progress,treeprint;
extern tree       curtree,bestree,bestree2;
extern contribarr probcat;
extern            double **contribution;
extern short       *alias,*ally,*location,*aliasweight;
extern FILE       *infile,*outfile,*treefile;
extern            Char **naym;
extern            short *enterorder;
extern            longer seed;



/* end inclusion of file dnamlk.h                                         */


FILE *infile, *outfile, *treefile;
short spp, nonodes, sites, weightsum, endsite, categs, inseed,
	    datasets, ith, i, j, l, jumb, njumble;
/* spp = number of species
  nonodes = number of nodes in tree
  sites = number of sites in actual sequences
  endsite = number of patterns in data
  numtrees = number of user-defined trees */
boolean freqsfrom, globle, jumble, lengths, trout, usertree,
	       weights, ctgry, ttr, auto_, printdata, progress, treeprint,
	       mulsets, firstset, interleaved, ibmpc, vt52, ansi;
tree curtree, bestree, bestree2;
Char **naym;   /* names of species */
double xi, xv, ttratio, ttratio0, freqa, freqc, freqg, freqt, freqr,
	      freqy, freqar, freqcy, freqgr, freqty, fracchange, sumrates,
	      lambda,lambda1;
longer seed;
short *enterorder;
sequence y;
short *weight,*alias,*aliasweight,*ally,*location;
double rate[maxcategs];
contribarr probcat;
double   **contribution;

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
        file[0]='\0';
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


void setuptree(a)
tree *a;
{
  short i;
  node *p;

  for (i = 0; i < spp; i++) {
    a->nodep[i]->number = i + 1;
    a->nodep[i]->tip = true;
    a->nodep[i]->initialized = false;
    a->nodep[i]->tyme = 0.0;
    a->nodep[i]->next = NULL;
    a->nodep[i]->back = NULL;
  }
  for (i = spp + 1; i <= nonodes; i++) {
    p = a->nodep[i - 1];
    for (j = 1; j <= 3; j++) {
      p->tip = false;
      p->initialized = false;
      p->number = i;
      p->tyme = 0.0;
      p->back = NULL;
      p = p->next;
    }
  }
  a->likelihood = -999999.0;
  a->root = NULL;
}  /* setuptree */


boolean letter(ch)
Char ch;
{
  /* tests whether ch is a letter
   -- works for both ASCII and EBCDIC codes */
  return ((ch >= 'A' && ch <= 'I') || (ch >= 'J' && ch <= 'R') ||
	  (ch >= 'S' && ch <= 'Z') || (ch >= 'a' && ch <= 'i') ||
	  (ch >= 'j' && ch <= 'r') || (ch >= 's' && ch <= 'z'));
}  /* letter */


void uppercase(ch)
Char *ch;
{
  (*ch) = isupper(*ch) ? (*ch) : toupper(*ch);
}  /* uppercase */


void getoptions()
{
  /* interactively set options */
  short i, j, inseed0,scanned;
  Char ch;
  char line[128];
   char rest[128];
  boolean done, done1, done2;
  double probsum;

  fprintf(outfile, "\nNucleic acid sequence\n");
  fprintf(outfile, "   Maximum Likelihood method with molecular ");
  fprintf(outfile, "clock, version %s\n\n",VERSION);
  putchar('\n');
  auto_ = false;
  ctgry = false;
  categs = 1;
  rate[0] = 1.0;
  freqsfrom = true;
  globle = false;
  jumble = false;
  njumble = 1;
  lambda = 1.0;
  lambda1 = 0.0;
  lengths = false;
  trout = true;
  ttratio = 2.0;
  ttr = false;
  usertree = false;
  weights = false;
  printdata = false;
  progress = true;
  treeprint = true;
  interleaved = true;
  probcat = (double *)Malloc(sizeof(double));
  probcat[0] = 1.0;
  do {
    if (ansi)
      printf("\033[2J\033[H");
    else if (vt52)
      printf("\033E\033H");
    else
      putchar('\n');
    printf("\nNucleic acid sequence\n");
    printf("   Maximum Likelihood method with molecular clock, version %s\n\n",
	   VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?");
    if (usertree)
      printf("  No, use user trees in input file\n");
    else
      printf("  Yes\n");
    if (usertree) {
      printf("  L           Use lengths from user tree?");
      if (lengths)
	printf("  Yes\n");
      else
	printf("  No\n");
    }
    printf("  T        Transition/transversion ratio:");
    if (!ttr)
      printf("  2.0\n");
    else
      printf("  %8.4f\n", ttratio);
    printf("  F       Use empirical base frequencies?");
    if (freqsfrom)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  C   One category of substitution rates?");
    if (!ctgry)
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
      printf("  G                Global rearrangements?");
      if (globle)
	printf("  Yes\n");
      else
	printf("  No\n");
      printf("  J   Randomize input order of sequences?");
      if (jumble)
        printf("  Yes (seed =%8hd,%3hd times)\n", inseed0, njumble);
      else
	printf("  No. Use input order\n");
    }
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2hd sets\n", datasets);
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?");
    if (interleaved)
      printf("  Yes\n");
    else
      printf("  No, sequential\n");
    printf("  0   Terminal type (IBM PC, VT52, ANSI)?");
    if (ibmpc)
      printf("  IBM PC\n");
    if (ansi)
      printf("  ANSI\n");
    if (vt52)
      printf("  VT52\n");
    if (!(ibmpc || vt52 || ansi))
      printf("  (none)\n");
    printf("  1    Print out the data at start of run");
    if (printdata)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  2  Print indications of progress of run");
    if (progress)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  3                        Print out tree");
    if (treeprint)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  4       Write out trees onto tree file?");
    if (trout)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      uppercase(&ch);
      if (ch == 'J' || ch == 'U' || ch == 'C' || ch == 'R' || ch == 'M' ||
	  ch == 'F' || ch == 'G' || ch == 'L' || ch == 'T' || ch == 'I' ||
	  ch == '1' || ch == '2' || ch == '3' || ch == '4' || ch == '0') {
	switch (ch) {

	case 'C':
	  ctgry = !ctgry;
	  if (ctgry) {
	    do {
	      printf("Number of categories (1-%ld)?\n",(long)maxcategs);
	      gets(line);
	      categs = (short)atoi(line);
	    } while (categs < 1 || categs > maxcategs);
	    if (probcat){
	      free(probcat);
	      free(rate);
	    }
	    probcat = (double *)Malloc(categs * sizeof(double));
	    for (;;){
	      printf("Rate for each category? (use a space to separate)\n");
	      gets(line);
	      done1 = true;
	      for (i = 0; i < categs; i++){
		scanned = sscanf(line,"%lf %[^\n]", &rate[i],rest);
		if ((scanned != 2 && i < (categs - 1)) ||
		    (scanned != 1 && i == (categs - 1))){
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
          if (!ctgry)
	    auto_ = false;
	  break;
	  
	case 'R':
	  auto_ = !auto_;
	  if (auto_) {
	    do {
	      printf(
		"Mean block length of sites having the same rate (greater than 1)?\n");
	      scanf("%lf%*[^\n]", &lambda);
	      getchar();
	    } while (lambda < 1.0);
	    lambda = 1.0 / lambda;
	    lambda1 = 1.0 - lambda;
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
	  globle = !globle;
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
	  break;

	case 'L':
	  lengths = !lengths;
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
  } while (!done);
}  /* getoptions */

void getnums()
{
  /* input number of species, number of sites */
  fprintf(outfile, "\n\n");
  fscanf(infile, "%hd%hd", &spp, &sites);
  fprintf(outfile, "%4hd Species, %4hd Sites\n", spp, sites);
  nonodes = spp * 2 - 1;
}  /* getnums */


void doinit()
{
  /* initializes variables */
  short i;
  node *p, *q;

  getnums();
  getoptions();
  y     = (Char **)Malloc(spp*sizeof(Char *));
  naym  = (char **)Malloc(spp*sizeof(char *));
  curtree.nodep = (node **)Malloc(nonodes * sizeof(node *));
  for (i=0;i<spp;i++)
    naym[i] = (char *)Malloc(nmlngth+1);
  for (i = 0; i < spp; i++)
    y[i] = (char *)Malloc(sites * sizeof(char));
  for (i = 0; i < spp; i++)
    curtree.nodep[i] = (node *)Malloc(sizeof(node));
  for (i = spp; i < nonodes; i++) {
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
  bestree.nodep = (node **)Malloc(nonodes * sizeof(node *));
  for (i = 0; i < spp; i++)
    bestree.nodep[i] = (node *)Malloc(sizeof(node));
  for (i = spp; i < nonodes; i++) {
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
  bestree2.nodep = (node **)Malloc(nonodes*sizeof(node *));
  for (i = 0; i < (spp); i++)
    bestree2.nodep[i] = (node *)Malloc(sizeof(node));
  for (i = spp; i < (nonodes); i++) {
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
      printf("BAD WEIGHT CHARACTER: %c -- WEIGHTS IN DNAMLK MUST BE 0 OR 1\n",
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

  fprintf(outfile, "\nSites are weighted as follows:\n");
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
    if (cursp != spp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4hd\n",
	     ith);
      exit(-1);
    }
    sites = curst;
  }
  for (i = 0; i < sites; i++)
    weight[i] = 1;
  weightsum = sites;
  extranum = 0;
  while (!(eoln(infile))) {
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'W')
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
    if (ch == 'W')
      inputweights();
    else{
      printf("ERROR: INCORRECT AUXILIARY OPTIONS LINE WHICH STARTS WITH %c\n",
	     ch);
      exit(-1);}
  }
  if (categs > 1) {
    fprintf(outfile, "\nSite category   Rate of change\n\n");
    for (i = 1; i <= categs; i++)
      fprintf(outfile, "%11hd%13.3f\n", i, rate[i - 1]);
    putc('\n', outfile);
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
  fracchange = xi * (2.0 * freqa * freqgr + 2.0 * freqc * freqty) +
      xv*(1.0 - freqa * freqa - freqc * freqc - freqg * freqg - freqt * freqt);
}  /* getbasefreqs */

void getdata()
{
  short i, j, k, l, basesread, basesnew;
  Char ch;
  boolean allread, done;

  if (printdata)
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
    while (i <= spp) {
      if ((interleaved && basesread == 0) || !interleaved) {
	for (j = 0; j < nmlngth; j++) {
	  if (eof(infile) || eoln(infile)){
	    printf("ERROR: END-OF-LINE OR END-OF-FILE");
	    printf(" IN THE MIDDLE OF A SPECIES NAME\n");
	    exit(-1);
	  }
	  naym[i - 1][j] = getc(infile);
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
	  if (ch != 'A' && ch != 'B' && ch != 'C' && ch != 'D' && ch != 'G' &&
	      ch != 'H' && ch != 'K' && ch != 'M' && ch != 'N' && ch != 'R' &&
	      ch != 'S' && ch != 'T' && ch != 'U' && ch != 'V' && ch != 'W' &&
	      ch != 'X' && ch != 'Y' && ch != '?' && ch != 'O' && ch != '-' &&
	      ch != '.') {
	    printf("ERROR: BAD BASE:%c AT POSITION%5ld OF SPECIES %3ld\n",
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
	exit(-1);}
      i++;
    }
    if (interleaved) {
      basesread = basesnew;
      allread = (basesread == sites);
    } else
      allread = (i > spp);
  }
  if (printdata) {
    for (i = 1; i <= ( (sites - 1) / 60 + 1);++i){
      for (j = 1; j <= spp; j++) {
	for (k = 0; k < nmlngth; k++)
	  putc(naym[j - 1][k], outfile);
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
	while (k <= spp && tied) {
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
      while (k <= spp && tied) {
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
  short i,j,k;
  node *p1,*p2,*p3;

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
    ally[alias[i - 1] - 1] = alias[i - 1];
    location[alias[i - 1] - 1] = i;
  }
  sumrates = 0.0;
  for (i = 0; i < categs; i++)
    sumrates += probcat[i] * rate[i];
  for (i = 0; i < categs; i++)
    rate[i] /= sumrates;

  for (i=0; i<spp;i++){
    curtree.nodep[i]->x=(phenotype)Malloc(endsite * sizeof(sitelike *));
    if (!usertree) {
      bestree.nodep[i]->x=(phenotype)Malloc(endsite * sizeof(sitelike *));
      if (njumble > 1)
	bestree2.nodep[i]->x=(phenotype)Malloc(endsite * sizeof(sitelike *));
    }
    for (j=0;j<endsite;++j){
      curtree.nodep[i]->x[j] = (sitelike *)Malloc(categs * sizeof(sitelike));
      if (!usertree) {
        bestree.nodep[i]->x[j] = (sitelike *)Malloc(categs*sizeof(sitelike));
        if (njumble > 1)
          bestree2.nodep[i]->x[j] = (sitelike *)Malloc(categs*sizeof(sitelike));
      }
    }
  }

  for (i=spp; i < nonodes ; i++){
    p1 = curtree.nodep[i];
    if (!usertree) {
      p2 = bestree.nodep[i];
      if (njumble > 1)
        p3 = bestree2.nodep[i];
    }
    for (k=0;k<3;++k) {
      p1->x = (phenotype)Malloc(endsite * sizeof(sitelike *));
      if (!usertree) {
	p2->x = (phenotype)Malloc(endsite * sizeof(sitelike *));
        if (njumble > 1)
          p3->x = (phenotype)Malloc(endsite * sizeof(sitelike *));
      }
      for (j=0;j<endsite;++j){
	p1->x[j] = (sitelike *)Malloc(categs * sizeof(sitelike));
	if (!usertree) {
	  p2->x[j] = (sitelike *)Malloc(categs * sizeof(sitelike));
          if (njumble > 1)
            p3->x[j] = (sitelike *)Malloc(categs * sizeof(sitelike));
        }
      }
      p1=p1->next;
      if (!usertree) {
	p2=p2->next;
        if (njumble > 1)
          p3=p3->next;
      }
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
    for (i = 0; i < spp; i++) {
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
    for (i = 0; i < spp; i++) {
      for (j = 0; j < endsite; j++) {
	w = aliasweight[j];
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
{  /* DNA Maximum Likelihood with molecular clock */
char infilename[100],outfilename[100],trfilename[100];
int i;
#ifdef MAC
  macsetup("Dnamlk","");
  argv[0] = "Dnamlk";
#endif

  openfile(&infile,INFILE,"r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);

  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  datasets = 1;
  mulsets = false;
  firstset = true;
  doinit();

  ttratio0    = ttratio;
  enterorder  = (short *)Malloc(spp*sizeof(short));
  weight      = (short *)Malloc(sites*sizeof(short));
  alias       = (short *)Malloc(sites*sizeof(short));
  aliasweight = (short *)Malloc(sites*sizeof(short));
  ally        = (short *)Malloc(sites*sizeof(short));
  location    = (short *)Malloc(sites*sizeof(short));

  if (trout)
    openfile(&treefile,TREEFILE,"w",argv[0],trfilename);
  for (ith = 1; ith <= datasets; ith++) {
    ttratio = ttratio0;
    if (datasets > 1) {
      fprintf(outfile, "Data set # %hd:\n\n",ith);
      if (progress)
        printf("\nData set # %hd:\n",ith);
    }
    getinput();
    contribution = (double **)Malloc(endsite * sizeof(double *));
    for (i=0;i<endsite;++i)
        contribution[i] = (double *)Malloc(categs * sizeof(double));
    if (ith == 1)
      firstset = false;
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
    for (i=0;i<endsite;++i)
        free (contribution[i]);
    free(contribution);
  }
  FClose(infile);
  FClose(outfile);
  FClose(treefile);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(trfilename);
#endif
  exit(0);
}  /* DNA Maximum Likelihood with molecular clock */

int eof(f)
FILE *f;
{
    register int ch;

    if (feof(f))
        return 1;
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
mem = (MALLOCRETURN *)calloc((size_t)1,x);
if (!mem)
  memerror();
else
  return (MALLOCRETURN *)mem;
}

