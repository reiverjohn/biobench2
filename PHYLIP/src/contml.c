#include "phylip.h"

/* version 3.56c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#define epsilon1        0.000001   /* small number */
#define epsilon2        0.02   /* not such a small number */
#define point           "."

#define smoothings      4   /* number of passes through smoothing algorithm */
#define namelength      10   /* number of characters max. in species name */
#define maxtrees        10   /* maximum number of user trees in KHT test */

#define ibmpc0          false
#define ansi0           true
#define vt520           false


typedef double *phenotype;
typedef Char naym[namelength];
typedef short longer[6];

typedef struct node {
  struct node *next, *back;
  boolean tip, iter;
  short number;
  phenotype view;
  naym nayme;
  double v, deltav, bigv, dist, xcoord;
  short ycoord, ymin, ymax;
} node;

typedef struct tree {
  node **nodep;
  double likelihood;
  node *start;
} tree;


Static FILE *infile, *outfile, *treefile;
Static short numsp, numsp1, numsp2, loci, inseed, totalleles, df, outgrno,
             col, datasets, ith, i, j, l, jumb, njumble=0;
Static short *alleles, *locus;
Static phenotype *x;
Static naym *nayms;
Static boolean all, contchars, global, jumble, lengths, outgropt, trout,
               usertree, printdata, progress, treeprint, mulsets, ibmpc,
               vt52, ansi, firstset;
Static longer seed;
Static short *enterorder;
Static tree curtree, priortree, bestree, bestree2;
short nextsp,numtrees,which,maxwhich; /* From maketree, propogated to global */
boolean succeeded;
double maxlogl;
double l0gl[maxtrees];
double *l0gf[maxtrees];
Char ch;



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

double randum(seed)
short *seed;
{
  /* random number generator -- slow but machine independent */
  short i, j, k, sum;
  longer mult, newseed;
  double x;

  mult[0] = 13;
  mult[1] = 24;
  mult[2] = 22;
  mult[3] = 6;
  for (i = 0; i <= 5; i++)
    newseed[i] = 0;
  for (i = 0; i <= 5; i++) {
    sum = newseed[i];
    k = i;
    if (i > 3)
      k = 3;
    for (j = 0; j <= k; j++)
      sum += mult[j] * seed[i - j];
    newseed[i] = sum;
    for (j = i; j <= 4; j++) {
      newseed[j + 1] += newseed[j] / 64;
      newseed[j] &= 63;
    }
  }
  memcpy(seed, newseed, sizeof(longer));
  seed[5] &= 3;
  x = 0.0;
  for (i = 0; i <= 5; i++)
    x = x / 64.0 + seed[i];
  x /= 4.0;
  return x;
}  /* randum */


void uppercase(ch)
Char *ch;
{
  /* convert a character to upper case -- either ASCII or EBCDIC */
   *ch = (islower(*ch) ?  toupper(*ch) : (*ch));
}  /* uppercase */


void getnums()
{
  /* read species numbers and number of characters or loci */
  fscanf(infile, "%hd%hd", &numsp, &loci);
  fprintf(outfile, "\n%4hd Populations, %4hd Loci\n", numsp, loci);
  numsp1 = numsp + 1;
  numsp2 = numsp * 2 - 2;
}  /* getnums */

void getoptions()
{
  /* interactively set options */
  short i, inseed0;
  Char ch;
  boolean done, done1;

  fprintf(outfile, "\nContinuous character Maximum Likelihood");
  fprintf(outfile, " method version %s\n\n",VERSION);
  putchar('\n');
  global = false;
  jumble = false;
  njumble = 1;
  lengths = false;
  outgrno = 1;
  outgropt = false;
  all = false;
  contchars = false;
  trout = true;
  usertree = false;
  printdata = false;
  progress = true;
  treeprint = true;
  do {
    printf(ansi ?  "\033[2J\033[H" :
	   vt52 ?  "\033E\033H"    : "\n");
    printf("\nContinuous character Maximum Likelihood");
    printf(" method version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  U                       Search for best tree?  %s\n",
	   (usertree ? "No, use user trees in input" : "Yes"));
    if (usertree) {
      printf("  L                Use lengths from user trees?%s\n",
	     (lengths ? "  Yes" : "  No"));
    }
    printf("  C  Gene frequencies or continuous characters?  %s\n",
	   (contchars ? "Continuous characters" : "Gene frequencies"));
    if (!contchars)
      printf("  A   Input file has all alleles at each locus?  %s\n",
	     (all ? "Yes" : "No, one allele missing at each"));
    printf("  O                              Outgroup root?  %s %hd\n",
	   (outgropt ? "Yes, at species number" :
                       "No, use as outgroup species"),outgrno);
    if (!usertree) {
      printf("  G                      Global rearrangements?  %s\n",
	     (global ? "Yes" : "No"));
      printf("  J           Randomize input order of species?");
      if (jumble)
        printf("  Yes (seed =%8hd,%3hd times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  M                 Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2hd sets\n", datasets);
    else
      printf("  No\n");
    printf("  0         Terminal type (IBM PC, VT52, ANSI)?  %s\n",
	   ibmpc ? "IBM PC" :
	   ansi  ? "ANSI"   :
	   vt52  ? "VT52"   : "(none");
    printf("  1          Print out the data at start of run  %s\n",
	   (printdata ? "Yes" : "No"));
    printf("  2        Print indications of progress of run  %s\n",
	   (progress ? "Yes" : "No"));
    printf("  3                              Print out tree  %s\n",
	   (treeprint ? "Yes" : "No"));
    printf("  4             Write out trees onto tree file?  %s\n",
	   (trout ? "Yes" : "No"));
    printf("\nAre these settings correct?");
    printf(" (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      if (strchr("JLOUGACM12340",ch) != NULL){

        switch (ch) {
        case 'A':
          if (!contchars)
            all = !all;
          break;

        case 'C':
          contchars = !contchars;
          break;

        case 'G':
          global = !global;
          break;

        case 'J':
          jumble = !jumble;
          if (jumble) {
            do {
              printf("Random number seed (must be odd)?\n");
              scanf("%hd%*[^\n]", &inseed);
              getchar();
            } while (!(inseed & 1));
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
              done1 = (outgrno >= 1 || outgrno <= numsp);
              if (!done1) {
                printf("BAD OUTGROUP NUMBER: %4hd\n", outgrno);
                printf("  Must be in range 1 -%2hd\n", numsp);
              }
            } while (done1 != true);
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


void doinit()
{
  /* initializes variables */
  short i, j, k, n;
  node *p, *q;

  getnums();
  getoptions();
  curtree.nodep = (node **)Malloc(numsp2*sizeof(node *));
  for (i = 0; i < numsp; i++)
    curtree.nodep[i] = (node *)Malloc(sizeof(node));

  n = 1;
  if (!usertree) {
    bestree.nodep = (node **)Malloc(numsp2*sizeof(node *));
    for (i = 0; i < numsp; i++)
      bestree.nodep[i] = (node *)Malloc(sizeof(node));
      if (!bestree.nodep[i])
    priortree.nodep = (node **)Malloc(numsp2*sizeof(node *));
    for (i = 0; i < numsp; i++)
      priortree.nodep[i] = (node *)Malloc(sizeof(node));
    n = 3;
    if (njumble > 1) {
      bestree2.nodep = (node **)Malloc(numsp2*sizeof(node *));
      for (i = 0; i < numsp; i++)
        bestree2.nodep[i] = (node *)Malloc(sizeof(node));
      n = 4;
    }
  }
  for (k = 1; k <= n; k++) {
    for (i = numsp1 - 1; i < numsp2; i++) {
      q = NULL;
      for (j = 1; j <= 3; j++) {
        p = (node *)Malloc(sizeof(node));
        p->next = q;
        q = p;
      }
      p->next->next->next = p;
      if (k == 1)
        curtree.nodep[i] = p;
      else if (n > 1) {
        if (k == 2)
          bestree.nodep[i] = p;
        else if (k == 3)
          priortree.nodep[i] = p;
        else
          bestree2.nodep[i] = p;
      }
    }
  }
  alleles = (short *)Malloc(loci*sizeof(short));
  if (contchars)
    locus = (short *)Malloc(loci*sizeof(short));
}  /* doinit */

void getalleles()
{
  /* set up number of alleles at loci */
  short i, j, k, n, cursp, curloc;
  node *p;

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
  if (contchars ) {
    totalleles = loci;
    for (i = 1; i <= loci; i++) {
      locus[i - 1] = i;
      alleles[i - 1] = 2;
    }
    df = loci;
  } else {
    totalleles = 0;
    fscanf(infile, "%*[^\n]");
    getc(infile);
    if (printdata) {
      fprintf(outfile, "\nNumbers of alleles at the loci:\n");
      fprintf(outfile, "------- -- ------- -- --- -----\n\n");
    }
    for (i = 1; i <= loci; i++) {
      if (eoln(infile)) {
        fscanf(infile, "%*[^\n]");
        getc(infile);
      }
      fscanf(infile, "%hd", &alleles[i - 1]);
      if (alleles[i - 1] <= 0) {
        printf("BAD NUMBER OF ALLELES: %4hd AT LOCUS %3hd\n", alleles[i - 1], i);
        exit(-1);
      }
      totalleles += alleles[i - 1];
      if (printdata)
        fprintf(outfile, "%4hd", alleles[i - 1]);
    }
    locus = (short *)Malloc(totalleles*sizeof(short));
    totalleles = 0;
    for (i = 1; i <= loci; i++) {
      for (j = totalleles; j < (totalleles + alleles[i - 1]); j++)
        locus[j] = i;
      totalleles += alleles[i - 1];
    }
    df = totalleles - loci;
  }
  for (i = 0; i < numsp; i++)
    curtree.nodep[i]->view = (phenotype)Malloc(totalleles*sizeof(double));
  n = 1;
  if (!usertree) {
    for (i = 0; i < numsp; i++)
      bestree.nodep[i]->view = (phenotype)Malloc(totalleles*sizeof(double));
    for (i = 0; i < numsp; i++)
      priortree.nodep[i]->view = (phenotype)Malloc(totalleles*sizeof(double));
    n = 3;
    if (njumble > 1) {
      for (i = 0; i < numsp; i++)
        bestree2.nodep[i]->view = (phenotype)Malloc(totalleles*sizeof(double));
      n = 4;
    }
  }
  for (k = 1; k <= n; k++) {
    for (i = numsp1 - 1; i < numsp2; i++) {
      if (k == 1)
        p = curtree.nodep[i];
      else if (n > 1) {
        if (k == 2)
          p = bestree.nodep[i];
        else if (k == 3)
          p = priortree.nodep[i];
        else
          p = bestree2.nodep[i];
        }
      for (j = 1; j <= 3; j++) {
        p->view = (phenotype)Malloc(totalleles*sizeof(double));
        p = p->next;
      }
    }
  }
  for (i = 0; i < numsp; i++)
    x[i] = (phenotype)Malloc(totalleles*sizeof(double));
  if (usertree)
    for (i = 0; i < maxtrees; i++)
      l0gf[i] = (double *)Malloc(totalleles*sizeof(double));
  if (printdata)
    putc('\n', outfile);
}  /* getalleles */

void getdata()
{
  /* read species data */
  short i, j, k, l, m, n, p;
  double sum;

  if (printdata) {
    fprintf(outfile, "\nName");
    if (contchars)
      fprintf(outfile, "                       Phenotypes\n");
    else
      fprintf(outfile, "                 Gene Frequencies\n");
    fprintf(outfile, "----");
    if (contchars)
      fprintf(outfile, "                       ----------\n");
    else
      fprintf(outfile, "                 ---- -----------\n");
    putc('\n', outfile);
    if (!contchars) {
      for (j = 1; j <= namelength - 8; j++)
        putc(' ', outfile);
      fprintf(outfile, "locus:");
      p = 1;
      for (j = 1; j <= loci; j++) {
	if (all)
	  n = alleles[j - 1];
	else
	  n = alleles[j - 1] - 1;
	for (k = 1; k <= n; k++) {
	  fprintf(outfile, "%10hd", j);
	  if (p % 6 == 0 && (all || p < df)) {
	    putc('\n', outfile);
	    for (l = 1; l <= namelength - 2; l++)
	      putc(' ', outfile);
	  }
	  p++;
	}
      }
      fprintf(outfile, "\n\n");
    }
  }
  for (i = 0; i < numsp; i++) {
    if (true) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
      for (j = 0; j < namelength; j++) {
	nayms[i][j] = getc(infile);
	if ( eof(infile) || eoln(infile)){
	  printf("ERROR: END-OF-LINE OR END-OF-FILE IN");
	  printf(" THE MIDDLE OF A SPECIES NAME\n");
	  exit(-1);}
	else if (printdata)
	  putc(nayms[i][j], outfile);
      }
      m = 1;
      p = 1;
      for (j = 1; j <= loci; j++) {
	sum = 0.0;
	if (contchars)
	  n = 1;
	else if (all)
	  n = alleles[j - 1];
	else
	  n = alleles[j - 1] - 1;
	for (k = 1; k <= n; k++) {
	  if (eoln(infile)) {
	    fscanf(infile, "%*[^\n]");
	    getc(infile);
	  }
	  fscanf(infile, "%lf", &x[i][m - 1]);
	  sum += x[i][m - 1];
	  if (!contchars && x[i][m - 1] < 0.0) {
	    printf("\nLOCUS%3hd IN SPECIES%3hd: AN ALLELE", j, i + 1);
	    printf(" FREQUENCY IS NEGATIVE\n");
	    exit(-1);
	  }
	  if (printdata) {
	    fprintf(outfile, "%10.5f", x[i][m - 1]);
	    if (p % 6 == 0 && (all || p < df)) {
	      putc('\n', outfile);
	      for (l = 1; l <= namelength; l++)
		putc(' ', outfile);
	    }
	  }
	  if (!contchars) {
	    if (x[i][m - 1] >= epsilon1)
	      x[i][m - 1] = sqrt(x[i][m - 1]);
	  }
	  p++;
	  m++;
	}
	if (all && fabs(sum - 1.0) > epsilon2) {
	  printf("\nLOCUS%3hd IN SPECIES%3hd: FREQUENCIES DO NOT ADD UP TO 1\n",
                   j, i + 1);
          exit(-1);
          }
	if (!all && !contchars) {
	  x[i][m - 1] = 1.0 - sum;
	  if (x[i][m - 1] >= epsilon1)
	    x[i][m - 1] = sqrt(x[i][m - 1]);
	  if (x[i][m - 1] < 0.0) {
	    if (x[i][m - 1] > -epsilon2) {
	      for (l = 0; l <= m - 2; l++)
		x[i][l] /= sqrt(sum);
	      x[i][m - 1] = 0.0;
	    } else {
	      printf("\n LOCUS%3hd IN SPECIES%3hd: ");
	      printf("FREQUENCIES ADD UP TO MORE THAN 1\n",j, i + 1);
	      exit(-1);
	    }
	  }
	  m++;
	}
      }
      if (printdata)
        putc('\n', outfile);
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  if (printdata)
    putc('\n', outfile);
}  /* getdata */


void getinput()
{
  /* reads the input data */
  getalleles();
  getdata();
}  /* getinput */


#define down            2
#define over            60


void setuptree(a)
tree *a;
{
  /* initialize a tree */
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
  a->likelihood = -99999.0;
  a->start = a->nodep[0];
}  /* setuptree */

void hookup(p, q)
node *p, *q;
{
  /* hook together two nodes */
  p->back = q;
  q->back = p;
}  /* hookup */


void sumlikely(p, q, sum)
node *p, *q;
double *sum;
{
  /* sum contribution to likelihood over forks in tree */
  short i;
  double term, sumsq, vee;
  double TEMP;

  if (!p->tip)
    sumlikely(p->next->back, p->next->next->back, sum);
  if (!q->tip)
    sumlikely(q->next->back, q->next->next->back, sum);
  if (p->back == q)
    vee = p->v;
  else
    vee = p->v + q->v;
  vee += p->deltav + q->deltav;
  if (vee <= 1.0e-10) {
    printf("ERROR: CHECK FOR TWO IDENTICAL  SPECIES ");
    printf("AND ELIMINATE ONE FROM THE DATA\n");
    exit(-1);
  }
  sumsq = 0.0;
  if (usertree && which <= maxtrees) {
    for (i = 0; i < loci; i++)
      l0gf[which - 1]
        [i] += (1 - alleles[i]) * log(vee) / 2.0;
  }
  for (i = 0; i < totalleles; i++) {
    TEMP = p->view[i] - q->view[i];
    term = TEMP * TEMP;
    if (usertree && which <= maxtrees)
      l0gf[which - 1]
        [locus[i] - 1] -= term / (2.0 * vee);
    sumsq += term;
  }
  (*sum) += df * log(vee) / -2.0 - sumsq / (2.0 * vee);
}  /* sumlikely */

double evaluate(t)
tree *t;
{
  /* evaluate likelihood of a tree */
  short i;
  double sum;

  sum = 0.0;
  if (usertree && which <= maxtrees) {
    for (i = 0; i < loci; i++)
      l0gf[which - 1][i] = 0.0;
  }
  sumlikely(t->start, t->start->back, &sum);
  if (usertree) {
    l0gl[which - 1] = sum;
    if (which == 1) {
      maxwhich = 1;
      maxlogl = sum;
    } else if (sum > maxlogl) {
      maxwhich = which;
      maxlogl = sum;
    }
  }
  t->likelihood = sum;
  return sum;
}  /* evaluate */

/* Local variables for update: */

double distance(p, q)
node *p, *q;
{
  /* distance between two nodes */
  short i;
  double sum;
  double TEMP;

  sum = 0.0;
  for (i = 0; i < totalleles; i++) {
    TEMP = p->view[i] - q->view[i];
    sum += TEMP * TEMP;
  }
  return sum;
}  /* distance */

void makedists(p)
node *p;
{
  short i;
  node *q;
  /* compute distances among three neighbors of a node */


  for (i = 1; i <= 3; i++) {
    q = p->next;
    p->dist = distance(p->back, q->back);
    p = q;
  }
}  /* makedists */

void makebigv(p,negatives)
node *p;
boolean *negatives;
{
  /* make new branch length */
  short i;
  node *temp, *q, *r;

  q = p->next;
  r = q->next;
  *negatives = false;
  for (i = 1; i <= 3; i++) {
    p->bigv = p->v + p->back->deltav;
    if (p->iter) {
      p->bigv = (p->dist + r->dist - q->dist) / (df * 2);
      p->back->bigv = p->bigv;
      if (p->bigv < p->back->deltav)
        *negatives = true;
    }
    temp = p;
    p = q;
    q = r;
    r = temp;
  }
}  /* makebigv */

void correctv(p)
node *p;
{
  /* iterate branch lengths if some are to be zero */
  node *q, *r, *temp;
  short i, j;
  double f1, f2, vtot;

  q = p->next;
  r = q->next;
  for (i = 1; i <= smoothings; i++) {
    for (j = 1; j <= 3; j++) {
      vtot = q->bigv + r->bigv;
      if (vtot > 0.0)
        f1 = q->bigv / vtot;
      else
        f1 = 0.5;
      f2 = 1.0 - f1;
      p->bigv = (f1 * r->dist + f2 * p->dist - f1 * f2 * q->dist) / df;
      p->bigv -= vtot * f1 * f2;
      if (p->bigv < p->back->deltav)
        p->bigv = p->back->deltav;
      p->back->bigv = p->bigv;
      temp = p;
      p = q;
      q = r;
      r = temp;
    }
  }
}  /* correctv */

void littlev(p)
node *p;
{
  /* remove part of it that belongs to other branches */
  short i;

  for (i = 1; i <= 3; i++) {
    if (p->iter)
      p->v = p->bigv - p->back->deltav;
    if (p->back->iter)
      p->back->v = p->v;
    p = p->next;
  }
}  /* littlev */

void nuview(p)
node *p;
{
  /* renew information about subtrees */
  short i, j;
  node *q, *r, *a, *b, *temp;
  double v1, v2, vtot, f1, f2;

  q = p->next;
  r = q->next;
  for (i = 1; i <= 3; i++) {
    a = q->back;
    b = r->back;
    v1 = q->bigv;
    v2 = r->bigv;
    vtot = v1 + v2;
    if (vtot > 0.0)
      f1 = v2 / vtot;
    else
      f1 = 0.5;
    f2 = 1.0 - f1;
    for (j = 0; j <totalleles; j++)
      p->view[j] = f1 * a->view[j] + f2 * b->view[j];
    p->deltav = v1 * f1;
    temp = p;
    p = q;
    q = r;
    r = temp;
  }
}  /* nuview */

void update(p)
node *p;
{
  /* update branch lengths around a node */
  boolean negatives;

  if (p->tip)
    return;
  makedists(p);
  makebigv(p,&negatives);
  if (negatives)
    correctv(p);
  littlev(p);
  nuview(p);
}  /* update */

void smooth(p)
node *p;
{
  /* go through tree getting new branch lengths and views */
  if (p->tip)
    return;
  update(p);
  smooth(p->next->back);
  smooth(p->next->next->back);
}  /* smooth */

void insert_(p, q)
node *p, *q;
{
  /* put p and q together and iterate info. on resulting tree */
  short i;

  hookup(p->next->next, q->back);
  hookup(p->next, q);
  for (i = 1; i <= smoothings; i++) {
    smooth(p);
    smooth(p->back);
  }
}  /* insert */

void copynode(c, d)
node *c, *d;
{
  /* make a copy of a node */
  memcpy(d->nayme, c->nayme, sizeof(naym));
  memcpy(d->view, c->view, totalleles*sizeof(double));
  d->v = c->v;
  d->iter = c->iter;
  d->deltav = c->deltav;
  d->bigv = c->bigv;
  d->dist = c->dist;
  d->xcoord = c->xcoord;
  d->ycoord = c->ycoord;
  d->ymin = c->ymin;
  d->ymax = c->ymax;
}  /* copynode */

void copy_(a, b)
tree *a, *b;
{
  /* make a copy of a tree */
  short i, j=0;
  node *p, *q;

  for (i = 0; i < numsp; i++) {
    copynode(a->nodep[i], b->nodep[i]);
    if (a->nodep[i]->back) {
      if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->number - 1])
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->number - 1];
      else if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->number - 1]->next)
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->number - 1]->next;
      else
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->number - 1]->next->next;
    }
    else b->nodep[i]->back = NULL;
  }
  for (i = numsp; i < numsp2; i++) {
    p = a->nodep[i];
    q = b->nodep[i];
    for (j = 1; j <= 3; j++) {
      copynode(p, q);
      if (p->back) {
        if (p->back == a->nodep[p->back->number - 1])
          q->back = b->nodep[p->back->number - 1];
        else if (p->back == a->nodep[p->back->number - 1]->next)
          q->back = b->nodep[p->back->number - 1]->next;
        else
          q->back = b->nodep[p->back->number - 1]->next->next;
      }
      else
        q->back = NULL;
      p = p->next;
      q = q->next;
    }
  }
  b->likelihood = a->likelihood;
  b->start = a->start;
}  /* copy */

void setuptip(m, t)
short m;
tree *t;
{
  /* initialize branch lengths and views in a tip */
  node *tmp;

  tmp = t->nodep[m - 1];
  memcpy(tmp->view, x[m - 1], totalleles*sizeof(double));
  memcpy(tmp->nayme, nayms[m - 1], sizeof(naym));
  tmp->deltav = 0.0;
  tmp->v = 0.0;
}  /* setuptip */

void buildnewtip(m, t,nextsp)
short m;
tree *t;
short nextsp;
{
  /* initialize and hook up a new tip */
  node *p;

  setuptip(m, t);
  p = t->nodep[nextsp + numsp - 3];
  hookup(t->nodep[m - 1], p);
}  /* buildnewtip */

void buildsimpletree(t)
tree *t;
{
  /* make and initialize a three-species tree */
  setuptip(enterorder[0], t);
  setuptip(enterorder[1], t);
  hookup(t->nodep[enterorder[0] - 1], t->nodep[enterorder[1] - 1]);
  buildnewtip(enterorder[2], t, nextsp);
  insert_(t->nodep[enterorder[2] - 1]->back, t->nodep[enterorder[0] - 1]);
}  /* buildsimpletree */

void addtraverse(p, q, contin)
node *p, *q;
boolean contin;
{
  /* traverse through a tree, finding best place to add p */
  insert_(p, q);
  numtrees++;
  if (evaluate(&curtree) > bestree.likelihood)
    copy_(&curtree, &bestree);
  copy_(&priortree, &curtree);
  if (!q->tip && contin) {
    addtraverse(p, q->next->back, contin);
    addtraverse(p, q->next->next->back, contin);
  }
}  /* addtraverse */


void re_move(p, q)
node **p, **q;
{
  /* remove p and record in q where it was */
  *q = (*p)->next->back;
  hookup(*q, (*p)->next->next->back);
  (*p)->next->back = NULL;
  (*p)->next->next->back = NULL;
  update(*q);
  update((*q)->back);
}  /* re_move */

void rearrange(p)
node *p;
{
  /* rearranges the tree, globally or locally */
  node *q, *r;

  if (!p->tip && !p->back->tip) {
    r = p->next->next;
    re_move(&r, &q );
    copy_(&curtree, &priortree);
    addtraverse(r, q->next->back, global && nextsp == numsp);
    addtraverse(r, q->next->next->back, global && nextsp == numsp);
    copy_(&bestree, &curtree);
    if (global && nextsp == numsp && progress)
      putchar('.');
    if (global && nextsp == numsp && !succeeded) {
      if (r->back->tip) {
        r = r->next->next;
        re_move(&r, &q );
        q = q->back;
        copy_(&curtree, &priortree);
        if (!q->tip) {
          addtraverse(r, q->next->back, true);
          addtraverse(r, q->next->next->back, true);
        }
        q = q->back;
        if (!q->tip) {
          addtraverse(r, q->next->back, true);
          addtraverse(r, q->next->next->back, true);
        }
        copy_(&bestree, &curtree);
      }
    }
  }
  if (!p->tip) {
    rearrange(p->next->back);
    rearrange(p->next->next->back);
  }
}  /* rearrange */


void coordinates(p, lengthsum, tipy,tipmax)
node *p;
double lengthsum;
short *tipy;
double *tipmax;
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;

  if (p->tip) {
    p->xcoord = lengthsum;
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
    if (lengthsum > (*tipmax))
      (*tipmax) = lengthsum;
    return;
  }
  q = p->next;
  do {
    coordinates(q->back, lengthsum + q->v, tipy,tipmax);
    q = q->next;
  } while ((p == curtree.start->back || p != q) &&
           (p != curtree.start->back || p->next != q));
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = lengthsum;
  if (p == curtree.start->back)
    p->ycoord = p->next->next->back->ycoord;
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */

void drawline(i, scale)
short i;
double scale;
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q;
  short n, j;
  boolean extra;
  node *r, *first, *last;
  boolean done;

  p = curtree.start->back;
  q = curtree.start->back;
  extra = false;
  if (i == p->ycoord && p == curtree.start->back) {
    if (p->number - numsp >= 10)
      fprintf(outfile, "-%2hd", p->number - numsp);
    else
      fprintf(outfile, "--%hd", p->number - numsp);
    extra = true;
  } else
    fprintf(outfile, "  ");
  do {
    if (!p->tip) {
      r = p->next;
      done = false;
      do {
        if (i >= r->back->ymin && i <= r->back->ymax) {
          q = r->back;
          done = true;
        }
        r = r->next;
      } while (!(done || (p != curtree.start->back && r == p) ||
                 (p == curtree.start->back && r == p->next)));
      first = p->next->back;
      r = p;
      while (r->next != p)
        r = r->next;
      last = r->back;
      if (p == curtree.start->back)
        last = p->back;
    }
    done = (p->tip || p == q);
    n = (long)(scale * (q->xcoord - p->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if (q->ycoord == i && !done) {
      if (p->ycoord != q->ycoord)
        putc('+', outfile);
      else
        putc('-', outfile);
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        if (q->number - numsp >= 10)
          fprintf(outfile, "%2hd", q->number - numsp);
        else
          fprintf(outfile, "-%hd", q->number - numsp);
        extra = true;
      } else {
        for (j = 1; j < n; j++)
          putc('-', outfile);
      }
    } else if (!p->tip) {
      if (last->ycoord > i && first->ycoord < i && i != p->ycoord) {
        putc('!', outfile);
        for (j = 1; j < n; j++)
          putc(' ', outfile);
      } else {
        for (j = 1; j <= n; j++)
          putc(' ', outfile);
      }
    } else {
      for (j = 1; j <= n; j++)
        putc(' ', outfile);
    }
    if (q != p)
      p = q;
  } while (!done);
  if (p->ycoord == i && p->tip) {
    for (j = 0; j < namelength; j++)
      putc(p->nayme[j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */

void printree()
{
  /* prints out diagram of the tree */
  short i;
  short tipy;
  double tipmax,scale;

  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  tipmax = 0.0;
  coordinates(curtree.start->back, 0.0, &tipy,&tipmax);
  scale = over / (tipmax + 0.0001);
  for (i = 1; i <= (tipy - down); i++)
    drawline(i,scale);
  putc('\n', outfile);
}  /* printree */

#undef down
#undef over

void treeout(p)
node *p;
{
  /* write out file with representation of final tree */
  short i, n, w;
  Char c;
  double x;

  if (p->tip) {
    n = 0;
    for (i = 1; i <= namelength; i++) {
      if (p->nayme[i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = p->nayme[i];
      if (c == ' ')
        c = '_';
      putc(c, treefile);
    }
    col += n;
  } else {
    putc('(', treefile);
    col++;
    treeout(p->next->back);
    putc(',', treefile);
    col++;
    if (col > 55) {
      putc('\n', treefile);
      col = 0;
    }
    treeout(p->next->next->back);
    if (p == curtree.start->back) {
      putc(',', treefile);
      col++;
      if (col > 45) {
        putc('\n', treefile);
        col = 0;
      }
      treeout(p->back);
    }
    putc(')', treefile);
    col++;
  }
  x = p->v;
  if (x > 0.0)
    w = (long)(0.43429448222 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.43429448222 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if (p == curtree.start->back)
    fprintf(treefile, ";\n");
  else {
    fprintf(treefile, ":%*.5f", w + 7, x);
    col += w + 8;
  }
}  /* treeout */

/* Local variables for summarize: */

void describe(p, chilow,chihigh)
node *p;
double chilow,chihigh;
{
  /* print out information for one branch */
  short i;
  node *q;
  double bigv, delta;

  q = p->back;
  fprintf(outfile, "%3hd       ", q->number - numsp);
  if (p->tip) {
    for (i = 0; i < namelength; i++)
      putc(p->nayme[i], outfile);
  } else
    fprintf(outfile, "%4hd      ", p->number - numsp);
  fprintf(outfile, "%15.5f", q->v);
  delta = p->deltav + p->back->deltav;
  bigv = p->v + delta;
  if (p->iter)
     fprintf(outfile, "   (%12.5f,%12.5f)",
             chilow * bigv - delta, chihigh * bigv - delta);
  fprintf(outfile, "\n");
  if (!p->tip) {
    describe(p->next->back, chilow,chihigh);
    describe(p->next->next->back, chilow,chihigh);
  }
}  /* describe */

void summarize(numtrees)
short numtrees;
{
  /* print out branch lengths etc. */
  double chilow,chihigh;

  fprintf(outfile, "\nremember: ");
  if (outgropt)
    fprintf(outfile, "(although rooted by outgroup) ");
  fprintf(outfile, "this is an unrooted tree!\n\n");
  fprintf(outfile, "Ln Likelihood = %11.5f\n", curtree.likelihood);
  if (!usertree)
    fprintf(outfile, "\nexamined %4hd trees\n", numtrees);
  if (df == 1) {
    chilow = 0.000982;
    chihigh = 5.02389;
  } else if (df == 2) {
    chilow = 0.05064;
    chihigh = 7.3777;
  } else {
    chilow = 1.0 - 2.0 / (df * 9);
    chihigh = chilow;
    chilow -= 1.95996 * sqrt(2.0 / (df * 9));
    chihigh += 1.95996 * sqrt(2.0 / (df * 9));
    chilow *= chilow * chilow;
    chihigh *= chihigh * chihigh;
  }
  fprintf(outfile, "\nBetween     And             Length");
  fprintf(outfile, "      Approx. Confidence Limits\n");
  fprintf(outfile, "-------     ---             ------");
  fprintf(outfile, "      ------- ---------- ------\n");
  describe(curtree.start->back->next->back, chilow,chihigh);
  describe(curtree.start->back->next->next->back, chilow,chihigh);
  describe(curtree.start, chilow,chihigh);
  fprintf(outfile, "\n\n");
  if (trout) {
    col = 0;
    treeout(curtree.start->back);
  }
}  /* summarize */

void getch(c)
Char *c;
{
  /* get next nonblank character */
  do {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    *c = getc(infile);
    if (*c == '\n')
      *c = ' ';
  } while (*c == ' ');
}  /* getch */

void findch(c, lparens,rparens)
Char c;
short *lparens,*rparens;
{
  /* skip forward in user tree until find character c */
  boolean done;

  done = false;
  while (!done) {
    if (c == ',') {
      if (ch == '(' || ch == ')' || ch == ':' || ch == ';') {
        printf("\nERROR IN USER TREE: ");
	printf("UNMATCHED PARENTHESIS OR MISSING COMMA\n");
        printf(" OR NON-TRIFURCATED BASE\n");
	exit(-1);
      } else if (ch == ',')
        done = true;
    } else if (c == ')') {
      if (ch == '(' || ch == ',' || ch == ':' || ch == ';') {
        printf("\nERROR IN USER TREE: UNMATCHED PARENTHESIS OR NON-BIFURCATED NODE\n");
	exit(-1);
      } else if (ch == ')') {
        (*rparens)++;
        if ((*lparens) > 0 && (*lparens) == (*rparens)) {
          if ((*lparens) == numsp - 2) {
            if (eoln(infile)) {
              fscanf(infile, "%*[^\n]");
              getc(infile);
            }
            ch = getc(infile);
            if (ch != ';') {
              printf( "\nERROR IN USER TREE: ");
	      printf("UNMATCHED PARENTHESIS OR MISSING SEMICOLON\n");
	      exit(-1);
            }
          }
        }
	done = true;
      }
    }
    if ((done && ch == ')') || (!done)) {
      if (eoln(infile)) {
        fscanf(infile, "%*[^\n]");
        getc(infile);
      }
      ch = getc(infile);
    }
  }
}  /* findch */

void processlength(p)
node *p;
{
  long digit, ordzero;
  double valyew, divisor;
  boolean pointread;

  ordzero = '0';
  pointread = false;
  valyew = 0.0;
  divisor = 1.0;
  getch(&ch);
  digit = ch - ordzero;
  while (((unsigned long)digit <= 9) | ch == '.'){
    if (ch == '.')
      pointread = true;
    else {
      valyew = valyew * 10.0 + digit;
      if (pointread)
	divisor *= 10.0;
    }
    getch(&ch);
    digit = ch - ordzero;
}
  if (lengths) {
    p->v = valyew / divisor;
    p->back->v = p->v;
    p->iter = false;
    p->back->iter = false;
  }
}  /* processlength */

#undef point   /*   ??? Why is this necessary */

void addelement(p, lparens,rparens,nextnode,names,nolengths)
node *p;
short *lparens,*rparens,*nextnode;
boolean *names, *nolengths;
{
  /* add one node to the user tree */
  node *q;
  short i, n;
  boolean found;
  Char str[namelength];

  do {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    ch = getc(infile);
  } while (ch == ' ');
 if (ch == '(') {
    (*lparens)++;
    if ((*lparens) > numsp - 2) {
      printf("\nERROR IN USER TREE: TOO MANY LEFT PARENTHESES\n");
      exit(-1);
    }
    (*nextnode)++;
    q = curtree.nodep[(*nextnode) - 1];
    hookup(p, q);
    addelement(q->next,lparens,rparens,nextnode,names,nolengths);
    findch(',', lparens,rparens);
    addelement(q->next->next, lparens,rparens,nextnode,names,nolengths);
    findch(')', lparens,rparens);
  }
  else {
    for (i = 0; i < namelength; i++)
      str[i] = ' ';
    n = 1;
    do {
      if (ch == '_')
        ch = ' ';
      str[n - 1] = ch;
      if (eoln(infile)) {
        fscanf(infile, "%*[^\n]");
        getc(infile);
      }
      ch = getc(infile);
      n++;
    } while (ch != ':' && ch != ',' && ch != ')' &&
             n <= namelength);
    n = 1;
    do {
      found = true;
      for (i = 0; i < namelength; i++)
        found = (found && str[i] == nayms[n - 1][i]);
      if (found) {
        if (names[n - 1] == false)
          names[n - 1] = true;
        else {
          printf("\nERROR IN USER TREE: DUPLICATE NAME FOUND -- ");
          for (i = 0; i < namelength; i++)
            putchar(curtree.nodep[n - 1]->nayme[i]);
          putchar('\n');
	  exit(-1);
        }
      } else
        n++;
    } while (!(n > numsp || found));
    if (n > numsp) {
      printf("CANNOT FIND SPECIES: ");
      for (i = 0; i < namelength; i++)
        putchar(str[i]);
      putchar('\n');
    }
    hookup(curtree.nodep[n - 1], p);
    if (curtree.start->number > n)
      curtree.start = curtree.nodep[n - 1];
  }
  if (ch == ':') {
    processlength(p);
    *nolengths = false;
  }
}  /* addelement */

void treeread()
{
  /* read in a user tree */
  node *p;
  short i;
  short nextnode=0,lparens=0,rparens=0;
  boolean *names, nolengths;

  curtree.start = curtree.nodep[numsp - 1];
  do {
    ch = getc(infile);
  } while (ch == ' ');
  if (ch != '(')
    return;
  names = (boolean *)Malloc(numsp*sizeof(boolean));
  for (i = 0; i < numsp; i++)
    names[i] = false;
  lparens = 1;
  rparens = 0;
  nolengths = true;
  nextnode = numsp + 1;
  p = curtree.nodep[nextnode - 1];
  for (i = 1; i <= 2; i++) {
    addelement(p, &(lparens),&(rparens),&(nextnode),names,&nolengths);
    p = p->next;
    findch(',', &(lparens),&(rparens));
  }
  addelement(p, &lparens,&rparens,&nextnode,names,&nolengths);
  if (nolengths && lengths)
    printf("\nNO LENGTHS FOUND IN INPUT FILE WITH LENGTH OPTION CHOSEN\n");
  findch(')', &lparens,&rparens);
  fscanf(infile, "%*[^\n]");
  getc(infile);
  free(names);
}  /* treeread */

void nodeinit(p)
node *p;
{
  /* initialize a node */
  node *q, *r;
  short i;

  if (p->tip)
    return;
  q = p->next->back;
  r = p->next->next->back;
  nodeinit(q);
  nodeinit(r);
  for (i = 0; i < totalleles; i++)
    p->view[i] = 0.5 * q->view[i] + 0.5 * r->view[i];
  if (p->iter)
    p->v = 0.1;
  if (p->back->iter)
    p->back->v = 0.1;
}  /* nodeinit */

void initrav(p)
node *p;
{
  /* traverse to initialize */
  if (p->tip)
    nodeinit(p->back);
  else {
    initrav(p->next->back);
    initrav(p->next->next->back);
  }
}  /* initrav */

void treevaluate()
{
  /* evaluate user-defined tree, iterating branch lengths */
  short i;
  double dummy;

  initrav(curtree.start);
  initrav(curtree.start->back);
  for (i = 1; i <= smoothings * 4; i++)
    smooth(curtree.start->back);
  dummy = evaluate(&curtree);
}  /* treevaluate */


void maketree()
{
  /* construct the tree */
  short i, j, k, n, num;
  double sum, sum2, sd;
  double TEMP;
  node *p;

  if (usertree) {
    fscanf(infile, "%hd%*[^\n]", &numtrees);
    getc(infile);
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      putc('\n', outfile);
    }
    setuptree(&curtree);
    for (which = 1; which <= numsp; which++)
      setuptip(which, &curtree);
    which = 1;
    while (which <= numtrees) {
      treeread();
      curtree.start = curtree.nodep[outgrno - 1];
      treevaluate();
      printree();
      summarize(numtrees);
      which++;
    }
    if (numtrees > 1 && loci > 1 ) {
      fprintf(outfile, "Tree    Ln L      Diff Ln L     Its S.D.");
      fprintf(outfile, "   Significantly worse?\n\n");
      num = (numtrees > maxtrees) ? maxtrees : numtrees;
      for (i = 1; i <= num; i++) {
        fprintf(outfile, "%3hd%12.5f", i, l0gl[i - 1]);
        if (maxwhich == i)
          fprintf(outfile, "  <------ best\n");
        else {
          sum = 0.0;
          sum2 = 0.0;
          for (j = 0; j < loci; j++) {
            sum += l0gf[maxwhich - 1][j] - l0gf[i - 1][j];
            TEMP = l0gf[maxwhich - 1][j] - l0gf[i - 1][j];
            sum2 += TEMP * TEMP;
          }
          sd = sqrt(loci / (loci - 1.0) * (sum2 - sum * sum / loci));
          fprintf(outfile, "%12.5f%12.4f", l0gl[i - 1] - maxlogl, sd);
          if (sum > 1.95996 * sd)
            fprintf(outfile, "           Yes\n");
          else
            fprintf(outfile, "            No\n");
        }
      }
      fprintf(outfile, "\n\n");
    }
  } else {
    if (jumb == 1) {
      setuptree(&curtree);
      setuptree(&priortree);
      setuptree(&bestree);
      if (njumble > 1) setuptree(&bestree2);
    }
    for (i = 1; i <= numsp; i++)
      enterorder[i - 1] = i;
    if (jumble) {
      for (i = 0; i < numsp; i++) {
        j = (long)(randum(seed) * numsp) + 1;
        k = enterorder[j - 1];
        enterorder[j - 1] = enterorder[i];
        enterorder[i] = k;
      }
    }
    nextsp = 3;
    buildsimpletree(&curtree);
    curtree.start = curtree.nodep[enterorder[0] - 1];
    if (jumb == 1) numtrees = 1;
    nextsp = 4;
    if (progress) {
      printf("\nAdding species:\n");
      printf("   ");
      for (i = 0; i < namelength; i++)
        putchar(nayms[enterorder[0] - 1][i]);
      printf("\n   ");
      for (i = 0; i < namelength; i++)
        putchar(nayms[enterorder[1] - 1][i]);
      printf("\n   ");
      for (i = 0; i < namelength; i++)
        putchar(nayms[enterorder[2] - 1][i]);
      putchar('\n');
    }
    while (nextsp <= numsp) {
      buildnewtip(enterorder[nextsp - 1], &curtree, nextsp);
      copy_(&curtree, &priortree);
      bestree.likelihood = -99999.0;
      addtraverse(curtree.nodep[enterorder[nextsp - 1] - 1]->back,
                  curtree.start->back, true );
      copy_(&bestree, &curtree);
      if (progress) {
        printf("   ");
        for (j = 0; j < namelength; j++)
          putchar(nayms[enterorder[nextsp - 1] - 1][j]);
        putchar('\n');
      }
      if (global && nextsp == numsp) {
        if (progress) {
          printf("\nDoing global rearrangements\n");
          printf("  !");
          for (i = 1; i <= numsp - 2; i++)
            putchar('-');
          printf("!\n");
          printf("   ");
        }
      }
      succeeded = true;
      while (succeeded) {
        succeeded = false;
        rearrange(curtree.start->back);
        if (global && nextsp == numsp)
          putc('\n', outfile);
      }
      if (njumble > 1) {
        if (jumb == 1 && nextsp == numsp)
          copy_(&bestree, &bestree2);
        else if (nextsp == numsp) {
          if (bestree2.likelihood < bestree.likelihood)
            copy_(&bestree, &bestree2);
        }
      }
      if (nextsp == numsp && jumb == njumble) {
        if (njumble > 1) copy_(&bestree2, &curtree);
        curtree.start = curtree.nodep[outgrno - 1];
        printree();
        summarize(numtrees);
      }
      nextsp++;
    }
  }
  if ( jumb < njumble)
    return;
  if (progress) {
    printf("\n\nOutput written to output file\n\n");
    if (trout)
      printf("Tree also written onto file\n\n");
  }
  for (i = 0; i < numsp; i++)
    free(curtree.nodep[i]->view);
  n = 1;
  if (!usertree) {
    for (i = 0; i < numsp; i++)
      free(bestree.nodep[i]->view);
    for (i = 0; i < numsp; i++)
      free(priortree.nodep[i]->view);
    n = 3;
  }
  for (k = 1; k <= n; k++) {
    for (i = numsp1 - 1; i < numsp2; i++) {
      if (k == 1)
        p = curtree.nodep[i];
      else if (n > 1) {
        if (k == 2)
          p = bestree.nodep[i];
        else
          p = priortree.nodep[i];
        }
      for (j = 1; j <= 3; j++) {
        free(p->view);
        p = p->next;
      }
    }
  }
  for (i = 0; i < numsp; i++)
    free(x[i]);
  if (!contchars)
    free(locus);
  if (usertree)
    for (i = 0; i < maxtrees; i++)
      free(l0gf[i]);
}  /* maketree */


main(argc, argv)
int argc;
Char *argv[];
{  /* main program */
  char infilename[100],outfilename[100],trfilename[100];
#ifdef MAC
  macsetup("Contml","");
  argv[0] = "Contml";
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
  if (trout)
    openfile(&treefile,TREEFILE,"w",argv[0],trfilename);
  x = (phenotype *)Malloc(numsp*sizeof(phenotype));
  nayms = (naym *)Malloc(numsp*sizeof(naym));
  enterorder = (short *)Malloc(numsp*sizeof(short));
  for (ith = 1; ith <= datasets; ith++) {
    getinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1) {
      fprintf(outfile, "Data set # %hd:\n\n", ith);
      if (progress)
        printf("\nData set # %hd:\n", ith);
    }
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
  }
  FClose(outfile);
  FClose(treefile);
  FClose(infile);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(trfilename);
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
mem = (MALLOCRETURN *)calloc(1,x);
if (!mem)
  memerror();
else
  return (MALLOCRETURN *)mem;
}
