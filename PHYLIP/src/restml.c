#include "phylip.h"

/* version 3.52c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define maxcutter       8    /* maximum number of bases in a site           */
#define maxtrees        8    /* maximum number of user trees                */
#define smoothings      10   /* number of passes in smoothing algorithm     */
#define iterations      10   /* number of iterates of EM for each branch    */
#define nmlngth         10   /* number of characters max. in species name   */

#define epsilon         0.00001 /* used in update                           */
#define extrap0         100.0   /* extrapolation factor to speed iteration  */
#define initialv        0.1     /* starting value of branch length          */
#define point           "."

#define ibmpc0          false
#define ansi0           true
#define vt520           false
#define down            2
#define over            60


typedef double sitelike[maxcutter + 1];
typedef sitelike *phenotype;
typedef Char **sequence;
typedef Char naym[nmlngth];
typedef short longer[6];
typedef double **transmatrix;
typedef transmatrix *transptr;

typedef struct node {
  struct node *next, *back;
  boolean tip, iter, initialized;
  short branchnum, number;
  phenotype x;
  naym nayme;
  double v;
  short xcoord, ycoord, ymin, ymax;
} node;

typedef struct tree {
  node **nodep;
  transptr trans, transprod;
  double likelihood;
  node *start;
} tree;

FILE *infile, *outfile, *treefile;
short numsp, numsp1, numsp2, sites, enzymes, endsite, weightsum,
            sitelength, inseed, outgrno, datasets, ith, i, j, l,
            jumb, njumble=0;
double extrapol;
boolean  global, jumble, lengths, outgropt, weights, trout, trunc8,
               usertree, printdata, progress, treeprint, mulsets, firstset,
               interleaved, ibmpc, vt52, ansi;
tree curtree, priortree, bestree, bestree2;
sequence y;
longer seed;
short *enterorder;
short *weight, *alias, *aliasweight;




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
      case 'a':
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


void uppercase(ch)
Char *ch;
{
  /* convert ch to upper case -- either ASCII or EBCDIC */
    *ch = (isupper (*ch) ? (*ch) : toupper(*ch));
}  /* uppercase */


void getnums()
{
  /* read and print out numbers of species and sites */
  fscanf(infile, "%hd%hd%hd", &numsp, &sites, &enzymes);
  if (printdata)
    fprintf(outfile, "%4hd Species, %4hd Sites,%4hd Enzymes\n",
            numsp, sites, enzymes);
  numsp1 = numsp + 1;
  numsp2 = numsp * 2 - 2;
}  /* getnums */

void getoptions()
{
  /* interactively set options */
  short i, inseed0;
  Char ch;
  boolean done, done1;

  fprintf(outfile, "\nRestriction site Maximum Likelihood");
  fprintf(outfile, " method, version %s\n\n",VERSION);
  putchar('\n');
  sitelength = 6;
  trunc8 = true;
  extrapol = extrap0;
  global = false;
  jumble = false;
  njumble = 1;
  lengths = false;
  outgrno = 1;
  outgropt = false;
  trout = true;
  usertree = false;
  weights = false;
  printdata = false;
  progress = true;
  treeprint = true;
  interleaved = true;
  for (;;) {
    printf("%s",           (ansi ? ("\033[2J\033[H") :
                            vt52 ? ("\033E\033H")    :
                            "\n"));
    printf("\nRestriction site Maximum Likelihood");
    printf(" method, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
           (usertree ? "No, use user trees in input file" : "Yes"));
    if (usertree) {
      printf("  N          Use lengths from user trees?  %s\n",
	     (lengths ? "Yes" : "No"));
    }
    printf("  A               Are all sites detected?  %s\n",
           (trunc8 ? "No" : "Yes"));
    if (!usertree) {
      printf("  G                Global rearrangements?  %s\n",
             (global ? "Yes" : "No"));
      printf("  J   Randomize input order of sequences?  ");
      if (jumble)
        printf("Yes (seed =%8hd,%3hd times)\n", inseed0, njumble);
      else
        printf("No. Use input order\n");
    }
    printf("  L                          Site length?%3hd\n",sitelength);
    printf("  O                        Outgroup root?  %s%3hd\n",
           (outgropt ? "Yes, at sequence number" :
                       "No, use as outgroup species"),outgrno);

    printf("  E                  Extrapolation factor%7.1f\n", extrapol);
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2hd sets\n", datasets);
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n",
           (interleaved ? "Yes" : "No, sequential"));
    printf("  0   Terminal type (IBM PC, VT52, ANSI)?  %s\n",
           ibmpc ? "IBM PC" :
            ansi  ? "ANSI"   :
            vt52  ? "vt52"    : "(none)");
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
           (treeprint ? "Yes" : "No"));
    printf("  4       Write out trees onto tree file?  %s\n",
           (trout ? "Yes" : "No"));
    printf("\nAre these settings correct?");
    printf(" (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (strchr("JOUNAGLTLEMI01234",ch) != NULL){
      switch (ch) {
	
      case 'A':
	trunc8 = !trunc8;
	break;
	
      case 'E':
	printf("Extrapolation factor?\n");
	do {
	  scanf("%lf%*[^\n]", &extrapol);
	  getchar();
	  if (extrapol <= 0.0)
	    printf("Must be positive!\n");
	} while (extrapol <= 0.0);
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
	do {
	  printf("New Sitelength?\n");
	  scanf("%hd%*[^\n]", &sitelength);
	  getchar();
	  if (sitelength < 1)
	    printf("BAD RESTRICTION SITE LENGTH: %hd\n", sitelength);
	} while (sitelength < 1);
	break;
	
      case 'N':
        lengths = !lengths;
        break;

      case 'O':
	outgropt = !outgropt;
	if (outgropt) {
	  done1 = true;
	  do {
	    printf("Type number of the outgroup:\n");
	    scanf("%hd%*[^\n]", &outgrno);
	    getchar();
	    done1 = (outgrno >= 1 && outgrno <= numsp);
	    if (!done1) {
	      printf("BAD OUTGROUP NUMBER: %4hd\n", outgrno);
	      printf("  Must be in range 1 -%2hd\n", numsp);
	    }
	  } while (done1 != true);
	} else outgrno = 1;
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
  curtree.trans = (transptr)Malloc(numsp2*sizeof(transmatrix));
  for (i=0;i<numsp2;++i){
    curtree.trans[i] = (transmatrix)Malloc((sitelength + 1) * sizeof (double *));
    for (j=0;j<sitelength+1;++j)
      curtree.trans[i][j] = (double *)Malloc((sitelength + 1) * sizeof(double));
  }
  curtree.transprod = (transptr)Malloc(numsp2*sizeof(transmatrix));
  for (i=0;i<numsp2;++i){
    curtree.transprod[i] = (transmatrix)Malloc((sitelength+ 1) * sizeof (double *));
    for (j=0;j<(sitelength + 1);++j)
      curtree.transprod[i][j] = (double *)Malloc((sitelength + 1)* sizeof(double));
  }
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
  bestree.trans = (transptr)Malloc(numsp2*sizeof(transmatrix));
  for (i=0;i<numsp2;++i){
    bestree.trans[i] = (transmatrix)Malloc((sitelength + 1) * sizeof (double *));
    for (j=0;j<(sitelength + 1);++j){
      bestree.trans[i][j] = (double *)Malloc((sitelength + 1)* sizeof(double));
      }
  }
  bestree.transprod = (transptr)Malloc(numsp2*sizeof(transmatrix));
  for (i=0;i<numsp2;++i){
    bestree.transprod[i] = (transmatrix)Malloc((sitelength + 1)* sizeof (double *));
    for (j=0;j<sitelength+1;++j)
      bestree.transprod[i][j] = (double *)Malloc((sitelength +1)* sizeof(double));
  }
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
  priortree.trans = (transptr)Malloc(numsp2*sizeof(transmatrix));
  for (i=0;i<numsp2;++i){
    priortree.trans[i] = (transmatrix)Malloc((sitelength + 1) * sizeof (double *));
    for (j=0;j<(sitelength + 1);++j){
      priortree.trans[i][j] = (double *)Malloc((sitelength + 1)* sizeof(double));
      }
  }
  priortree.transprod = (transptr)Malloc(numsp2*sizeof(transmatrix));
  for (i=0;i<numsp2;++i){
    priortree.transprod[i] = (transmatrix)Malloc((sitelength + 1)* sizeof (double *));
    for (j=0;j<sitelength+1;++j)
      priortree.transprod[i][j] = (double *)Malloc((sitelength +1)* sizeof(double));
  }
  priortree.nodep = (node **)Malloc(numsp2*sizeof(node *));
  for (i = 0; i < numsp; i++)
    priortree.nodep[i] = (node *)Malloc(sizeof(node));
  for (i = numsp1 - 1; i < numsp2; i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc(sizeof(node));
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    priortree.nodep[i] = p;
  }
  if (njumble == 1) return;
  bestree2.trans = (transptr)Malloc(numsp2*sizeof(transmatrix));
  for (i=0;i<numsp2;++i){
    bestree2.trans[i] = (transmatrix)Malloc((sitelength + 1) * sizeof (double *));
    for (j=0;j<sitelength + 1;++j)
      bestree2.trans[i][j] = (double *)Malloc((sitelength +1)* sizeof(double));
  }
  bestree2.transprod = (transptr)Malloc(numsp2*sizeof(transmatrix));
  for (i=0;i<numsp2;++i){
    bestree2.transprod[i] = (transmatrix)Malloc((sitelength +1)* sizeof (double *));
    for (j=0;j<sitelength+1;++j)
      bestree2.transprod[i][j] = (double *)Malloc((sitelength +1)* sizeof(double));
  }
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
  for (i = 1; i <= sites; i++) {
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
      printf("BAD WEIGHT CHARACTER: %c -- WEIGHTS IN RESTML MUST BE 0 OR 1\n",
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

  fprintf(outfile, "\n   Sites are weighted as follows:\n\n");
  for (i = 1; i <= sites; i++) {
    if ((i - 1) % 60 == 0) {
      putc('\n', outfile);
      for (j = 1; j <= nmlngth + 3; j++)
        putc(' ', outfile);
    }
    fprintf(outfile, "%hd", weight[i]);
    if (i % 10 == 0 && i % 60 != 0)
      putc(' ', outfile);
  }
  fprintf(outfile, "\n\n");
}  /* printweights */

void inputoptions()
{
  /* read the options information */
  Char ch;
  short i, extranum, cursp, curst, curenz;

  if (!firstset) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%hd%hd%hd", &cursp, &curst, &curenz);
    if (cursp != numsp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4hd\n",
             ith);
      exit(-1);
    }
    if (curenz != enzymes) {
      printf("\nERROR: INCONSISTENT NUMBER OF ENZYMES IN DATA SET %4hd\n",
             ith);
      exit(-1);
    }
    sites = curst;
  }
  for (i = 1; i <= sites; i++)
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
    if (ch != 'W') {
      printf("ERROR: INCORRECT AUXILIARY OPTIONS LINE");
      printf(" WHICH STARTS WITH %c\n", ch);
      exit(-1);
    }
    else
      inputweights();
  }
  fprintf(outfile, "\n  Recognition sequences all%2hd bases long\n",
          sitelength);
  if (trunc8)
    fprintf(outfile,
      "\nSites absent from all species are assumed to have been omitted\n\n");
  if (weights)
    printweights();
}  /* inputoptions */

void setuptree(a)
tree *a;
{
  /* set up data structures for a tree */
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
  /* read the species and sites data */
  short i, j, k, l, sitesread, sitesnew;
  Char ch;
  boolean allread, done;

  if (printdata)
    putc('\n', outfile);
  j = nmlngth + (sites + (sites - 1) / 10) / 2 - 5;
  if (j < nmlngth - 1)
    j = nmlngth - 1;
  if (j > 39)
    j = 39;
  if (printdata) {
    fprintf(outfile, "Name");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "Sites\n");
    fprintf(outfile, "----");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "-----\n\n");
  }
  setuptree(&curtree);
  sitesread = 0;
  allread = false;
  while (!(allread)) {
    allread = true;
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    i = 1;
    while (i <= numsp ) {
      if ((interleaved && sitesread == 0) || !interleaved) {
        for (j = 0; j < nmlngth; j++) {
	  if (eof(infile) | eoln(infile)){
	    printf("ERROR: END-OF-LINE OR END-OF-FILE");
	    printf(" IN THE MIDDLE OF A SPECIES NAME\n");
	    exit(-1);
	  }
	  curtree.nodep[i - 1]->nayme[j] = getc(infile);
        }
      }
      if (interleaved)
        j = sitesread;
      else
        j = 0;
      done = false;
      while (!done && !eof(infile)) {
        if (interleaved)
          done = true;
        while (j < sites && !(eoln(infile) || eof(infile))) {
          ch = getc(infile);
          if (ch == '\n')
            ch = ' ';
          if (ch == ' ')
            continue;
          uppercase(&ch);
          if (ch != '1' && ch != '0' && ch != '+' && ch != '-' && ch != '?' &&
              ch != '.') {
            printf(" WARNING -- BAD SYMBOL %c",ch);
	    printf(" AT POSITION %5hd OF SPECIES %3hd\n",j,i);
	    exit(-1);
          }
          if (ch == '1')
            ch = '+';
          if (ch == '0')
            ch = '-';
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
        sitesnew = j;
      fscanf(infile, "%*[^\n]");
      getc(infile);
      if ((interleaved && j != sitesnew ) || ((!interleaved) && j != sites)){
        printf("ERROR: SEQUENCES OUT OF ALIGNMENT\n");
	exit(-1);}
      i++;
    }
    if (interleaved) {
      sitesread = sitesnew;
      allread = (sitesread == sites);
    } else
      allread = (i > numsp);
  }
  if (printdata) {
    for (i = 1; i <= ((sites - 1) / 60 + 1); i++) {
      for (j = 0; j < numsp; j++) {
        for (k = 0; k < nmlngth; k++)
          putc(curtree.nodep[j]->nayme[k], outfile);
        fprintf(outfile, "   ");
        l = i * 60;
        if (l > sites)
          l = sites;
        for (k = (i - 1) * 60 + 1; k <= l; k++) {
	  putc(y[j][k - 1], outfile);
	  if (k % 10 == 0 && k % 60 != 0)
	    putc(' ', outfile);
        }
        putc('\n', outfile);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
  if (!usertree) {
    setuptree(&priortree);
    setuptree(&bestree);
    if (njumble > 1) setuptree(&bestree2);
  }
}  /* getdata */

void sitesort()
{
  /* Shell sort keeping alias, aliasweight in same order */
  short gap, i, j, jj, jg, k, itemp;
  boolean flip, tied;

  gap = sites / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= sites; i++) {
      j = i - gap;
      flip = true;
      while (j > 0 && flip) {
        jj = alias[j];
        jg = alias[j + gap];
        flip = false;
        tied = true;
        k = 1;
        while (k <= numsp && tied) {
          flip = (y[k - 1][jj - 1] > y[k - 1][jg - 1]);
          tied = (tied && y[k - 1][jj - 1] == y[k - 1][jg - 1]);
          k++;
        }
        if (tied) {
          aliasweight[j] += aliasweight[j + gap];
          aliasweight[j + gap] = 0;
        }
        if (!flip)
          break;
        itemp = alias[j];
        alias[j] = alias[j + gap];
        alias[j + gap] = itemp;
        itemp = aliasweight[j];
        aliasweight[j] = aliasweight[j + gap];
        aliasweight[j + gap] = itemp;
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
      k = 1;
      while (k <= numsp && tied) {
        tied = (tied && y[k - 1][alias[i] - 1] == y[k - 1][alias[j] - 1]);
        k++;
      }
      if (tied && aliasweight[j] > 0) {
        aliasweight[i] += aliasweight[j];
        aliasweight[j] = 0;
        alias[j] = alias[i];
      }
      j++;
    }
    i = j - 1;
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
    if (aliasweight[i] > 0)
      i++;
    else {
      if (j <= i)
        j = i + 1;
      if (j <= sites) {
        found = false;
        do {
          found = (aliasweight[j] > 0);
          j++;
        } while (!(found || j > sites));
        if (found) {
          j--;
          itemp = alias[i];
          alias[i] = alias[j];
          alias[j] = itemp;
          itemp = aliasweight[i];
          aliasweight[i] = aliasweight[j];
          aliasweight[j] = itemp;
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
  short i,j;
  node *p;

  for (i = 1; i <= sites; i++) {
    alias[i] = i;
    aliasweight[i] = weight[i];
  }
  sitesort();
  sitecombine();
  sitescrunch();
  for (i = 1; i <= sites; i++) {
    weight[i] = aliasweight[i];
    if (weight[i] > 0)
      endsite = i;
  }
  weight[0] = 1;
  for (i = 0; i < numsp; i++)
    curtree.nodep[i]->x = (phenotype)Malloc((endsite+1)*sizeof(sitelike));
  for (i = numsp; i < numsp2; i++) {
    p = curtree.nodep[i];
    for (j = 1; j <= 3; j++) {
      p->x = (phenotype)Malloc((endsite+1)*sizeof(sitelike));
      p = p->next;
    }
  }
  if (!usertree) {
    for (i = 0; i < numsp; i++)
      priortree.nodep[i]->x = (phenotype)Malloc((endsite+1)*sizeof(sitelike));
    for (i = numsp; i < numsp2; i++) {
      p = priortree.nodep[i];
      for (j = 1; j <= 3; j++) {
        p->x = (phenotype)Malloc((endsite+1)*sizeof(sitelike));
        p = p->next;
      }
    }
    for (i = 0; i < numsp; i++)
      bestree.nodep[i]->x = (phenotype)Malloc((endsite+1)*sizeof(sitelike));
    for (i = numsp; i < numsp2; i++) {
      p = bestree.nodep[i];
      for (j = 1; j <= 3; j++) {
        p->x = (phenotype)Malloc((endsite+1)*sizeof(sitelike));
        p = p->next;
      }
    }
    if (njumble > 1) {
      for (i = 0; i < numsp; i++)
        bestree2.nodep[i]->x = (phenotype)Malloc((endsite+1)*sizeof(sitelike));
      for (i = numsp; i < numsp2; i++) {
        p = bestree2.nodep[i];
        for (j = 1; j <= 3; j++) {
          p->x = (phenotype)Malloc((endsite+1)*sizeof(sitelike));
          p = p->next;
        }
      }
    }
  }
}  /* makeweights */

void makevalues()
{
  /* set up fractional likelihoods at tips */
  short i, j, k, l, m;
  boolean found;

  for (k = 1; k <= endsite; k++) {
    j = alias[k];
    for (i = 0; i < numsp; i++) {
      for (l = 0; l <= sitelength; l++)
        curtree.nodep[i]->x[k][l] = 1.0;
      switch (y[i][j - 1]) {

      case '+':
        for (m = 1; m <= sitelength; m++)
          curtree.nodep[i]->x[k][m] = 0.0;
        break;

      case '-':
        curtree.nodep[i]->x[k][0] = 0.0;
        break;

      case '?':
        /* blank case */
        break;
      }
    }
  }
  for (i = 0; i < numsp; i++) {
    for (k = 1; k <= sitelength; k++)
      curtree.nodep[i]->x[0][k] = 1.0;
    curtree.nodep[i]->x[0][0] = 0.0;
  }
  if (trunc8)
    return;
  found = false;
  i = 1;
  while (!found && i <= endsite) {
    found = true;
    for (k = 0; k < numsp; k++)
      found = (found && y[k][alias[i] - 1] == '-');
    if (!found)
      i++;
  }
  if (found) {
    weightsum += (enzymes - 1) * weight[i];
    weight[i] *= enzymes;
  }
}  /* makevalues */


void getinput()
{
  /* reads the input data */
  inputoptions();
  getdata();
  makeweights();
  makevalues();
}  /* getinput */


main(argc, argv)
int argc;
Char *argv[];
{  /* maximum likelihood phylogenies from restriction sites */
char infilename[100],outfilename[100],trfilename[100];
#ifdef MAC
  macsetup("Restml","");
  argv[0] = "Restml";
#endif
  openfile(&infile,INFILE,"r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);
  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  mulsets = false;
  datasets = 1;
  firstset = true;
  doinit();
  if (trout)
    openfile(&treefile,TREEFILE,"w",argv[0],trfilename);
  enterorder = (short *)Malloc(numsp*sizeof(short));
  weight = (short *)Malloc((sites+1)*sizeof(short));
  alias = (short *)Malloc((sites+1)*sizeof(short));
  aliasweight = (short *)Malloc((sites+1)*sizeof(short));
  for (ith = 1; ith <= datasets; ith++) {
    if (datasets > 1) {
      fprintf(outfile, "Data set # %hd:\n",ith);
      if (progress)
        printf("\nData set # %hd:\n",ith);
    }
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
}  /* maximum likelihood phylogenies from restriction sites */

/*
  misc. support routines.  Some of these should eventually be removed
*/

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
mem = (MALLOCRETURN *)calloc(1,x);
if (!mem)
  memerror();
else
  return (MALLOCRETURN *)mem;
}

