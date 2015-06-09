#include "phylip.h"

/* version 3.56c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define smoothings      4    /* number of passes through smoothing algorithm */
#define namelength      10   /* number of characters max. in species name    */
#define epsilon         0.000001   /* a very small but not too small number  */

#define ibmpc0          false
#define ansi0           true
#define vt520           false

typedef double *vector;
typedef short *intvector;
typedef Char naym[namelength];
typedef short longer[6];

typedef struct node {
  struct node *next, *back;
  boolean tip, iter;
  short number;
  naym nayme;
  vector d, w;
  double v, dist;
  short xcoord, ycoord, ymin, ymax;
} node;

typedef struct tree {
  node **nodep;
  double likelihood;
  node *start;
} tree;


FILE *infile, *outfile, *treefile;
short numsp,numsp1,numsp2,nums,inseed,outgrno,col,datasets,ith,
             i, j, l, jumb, njumble=0;
vector *x;
intvector *reps;
naym *nayms;
boolean global,jumble,lengths,usertree,lower,upper, negallowed,
               outgropt,replicates, trout,  printdata, progress, treeprint,
               mulsets, ibmpc, vt52, ansi, firstset=false;
double power;
longer seed;
short *enterorder;
tree curtree, priortree, bestree, bestree2;
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
  short i,j,k,sum;
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
  /* read species number */
  fscanf(infile, "%hd", &numsp);
  fprintf(outfile, "\n%4hd Populations\n", numsp);
  numsp1 = numsp + 1;
  numsp2 = numsp * 2 - 2;
}  /* getnums */

void getoptions()
{
  /* interactively set options */
  short i, inseed0=0;
  Char ch;
  boolean done=false, done1=false;

  fprintf(outfile, "\nFitch-Margoliash method version %s\n\n",VERSION);
  putchar('\n');
  global = false;
  jumble = false;
  njumble = 1;
  lengths = false;
  lower = false;
  negallowed = false;
  outgrno = 1;
  outgropt = false;
  power = 2.0;
  replicates = false;
  trout = true;
  upper = false;
  usertree = false;
  printdata = false;
  progress = true;
  treeprint = true;
  do {
    printf(ansi ? "\033[2J\033[H" :
	   vt52 ? "\033E\033H"    : "\n");
    printf("\nFitch-Margoliash method version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
	   (usertree ? "No, use user trees in input file" : "Yes"));
    if (usertree) {
      printf("  N          Use lengths from user trees?  %s\n",
	     (lengths ? "Yes" : "No"));
    }
    printf("  P                                Power?%9.5f\n",power);
    printf("  -      Negative branch lengths allowed?  %s\n",
	   negallowed ? "Yes" : "No");
    printf("  O                        Outgroup root?");
    if (outgropt)
      printf("  Yes, at species number%3hd\n", outgrno);
    else
      printf("  No, use as outgroup species%3hd\n", outgrno);
    printf("  L         Lower-triangular data matrix?");
    if (lower)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  R         Upper-triangular data matrix?");
    if (upper)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  S                        Subreplicates?");
    if (replicates)
      printf("  Yes\n");
    else
      printf("  No\n");
    if (!usertree) {
      printf("  G                Global rearrangements?");
      if (global)
        printf("  Yes\n");
      else
        printf("  No\n");
      printf("  J     Randomize input order of species?");
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
    printf(
   "\nAre these settings correct? (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    done = (ch == 'Y');
   if (!done) {
      if (ch == 'J' || ch == 'O' || ch == 'U' || ch == 'N' || ch == 'P' ||
          ch == 'G' || ch == '-' || ch == 'L' || ch == 'R' || ch == 'S' ||
          ch == 'M' || ch == '0' || ch == '1' || ch == '2' || ch == '3' ||
          ch == '4') {
        switch (ch) {

        case '-':
          negallowed = !negallowed;
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
          lower = !lower;
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
          }
          break;

        case 'P':
          printf("New power?\n");
          scanf("%lf%*[^\n]", &power);
          getchar();
          break;

        case 'R':
          upper = !upper;
          break;

        case 'S':
          replicates = !replicates;
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
  short i=0, j=0, k=0, n=0;
  node *p, *q;

  getnums(); /* set numsp1/numsp2 */
  /* initialize the X structure */
  x = (vector *)Malloc(numsp*sizeof(vector));
  reps = (intvector *)Malloc(numsp*sizeof(intvector));
  for (i=0;i<numsp;++i){
    x[i]=(vector)Malloc(numsp2 * sizeof(double));
    reps[i]=(intvector)Malloc(numsp * sizeof(short));
  }

  getoptions();
  curtree.nodep = (node **)Malloc(numsp2*sizeof(node *));
  for (i = 0; i < numsp; i++){
    curtree.nodep[i] = (node *)Malloc(sizeof(node));
    if (curtree.nodep[i] == NULL)
      memerror();
    curtree.nodep[i]->d = (vector)Malloc(numsp2 * sizeof(double));
    if (!curtree.nodep[i]->d)
      memerror();
    curtree.nodep[i]->w = (vector)Malloc(numsp2 * sizeof(double));
    if (!curtree.nodep[i]->w)
      memerror();  }

  n = 1;
  if (!usertree) {
    bestree.nodep = (node **)Malloc(numsp2*sizeof(node *));
    for (i = 0; i < numsp; i++){
      bestree.nodep[i] = (node *)Malloc(sizeof(node));
      if (bestree.nodep[i]==NULL)
        memerror();
      bestree.nodep[i]->d = (vector)Malloc(numsp2 * sizeof(double));
      if (!bestree.nodep[i]->d)
        memerror();
      bestree.nodep[i]->w = (vector)Malloc(numsp2 * sizeof(double));
      if (!bestree.nodep[i]->w)
        memerror();}

    priortree.nodep = (node **)Malloc(numsp2*sizeof(node *));
    for (i = 0; i < numsp; i++){
      priortree.nodep[i] = (node *)Malloc(sizeof(node));
      if (priortree.nodep[i]==NULL)
        memerror();
      priortree.nodep[i]->d = (double *)Malloc(numsp2 * sizeof(double));
      if (!priortree.nodep[i]->d)
        memerror();
      priortree.nodep[i]->w = (double *)Malloc(numsp2 * sizeof(double));
      if (!priortree.nodep[i]->w)
        memerror(); }
    n = 3;

    if (njumble > 1) {
      bestree2.nodep = (node **)Malloc(numsp2*sizeof(node *));
      for (i = 0; i < numsp; i++){
        bestree2.nodep[i] = (node *)Malloc(sizeof(node));
        if (bestree2.nodep[i]==NULL)
          memerror();
        bestree2.nodep[i]->d = (double *)Malloc(numsp2 * sizeof(double));
        if (!bestree2.nodep[i]->d)
          memerror();
        bestree2.nodep[i]->w = (double *)Malloc(numsp2 * sizeof(double));
        if (!bestree2.nodep[i]->w)
          memerror();}
      n = 4;
    }
  }
  for (k = 1; k <= n; k++) {
    for (i = numsp1 - 1; i < numsp2; i++) {
      q = NULL;
      for (j = 1; j <= 3; j++) {
        p = (node *)Malloc(sizeof(node));
        if (!p)
          memerror();
        p->d = (double *)Malloc(numsp2 * sizeof(double));
        if (!p->d)
          memerror();
        p->w = (double *)Malloc(numsp2 * sizeof(double));
        if (!p->w)
          memerror();
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
}  /* doinit */


void inputoptions()
{
  /* read options information */
  Char ch;
  short cursp=0;

  if (!firstset) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%hd", &cursp);
    if (cursp != numsp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4hd\n",ith);
      exit(-1);
    }
  }
  while (!(eoln(infile))) {
    ch = getc(infile);
    uppercase(&ch);
    if (ch != ' ') {
      printf("BAD OPTION CHARACTER: %c\n", ch);
      exit(-1);
    }
  }
  fprintf(outfile, "                  __ __             2\n");
  fprintf(outfile, "                  \\  \\   (Obs - Exp)\n");
  fprintf(outfile, "Sum of squares =  /_ /_  ------------\n");
  fprintf(outfile, "                               ");
  if (power == (short)power)
    fprintf(outfile, "%2hd\n", (short)power);
  else
    fprintf(outfile, "%4.1f\n", power);
  fprintf(outfile, "                   i  j      Obs\n\n");
  fprintf(outfile, "Negative branch lengths ");
  if (!negallowed)
    fprintf(outfile, "not ");
  fprintf(outfile, "allowed\n\n");
  if (global)
    fprintf(outfile, "global optimization\n\n");
}  /* inputoptions */


Static void getinput()
{
  /* reads the input data */
  inputoptions();
}  /* getinput */



#define down            2
#define over            60


void setuptree(a)
tree *a;
{
  /* initialize a tree */
  short i=0, j=0;
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
  a->likelihood = -1.0;
  a->start = a->nodep[0];
}  /* setuptree */

void getdata()
{
  /* read in distance matrix */
  short i=0, j=0, k=0, columns=0;
  boolean skipit=false, skipother=false;

  if (replicates)
    columns = 4;
  else
    columns = 6;
  if (printdata) {
    fprintf(outfile, "\nName                       Distances");
    if (replicates)
      fprintf(outfile, " (replicates)");
    fprintf(outfile, "\n----                       ---------");
    if (replicates)
      fprintf(outfile, "-------------");
    fprintf(outfile, "\n\n");
  }
  for (i = 0; i < numsp; i++) {
    x[i][i] = 0.0;
    fscanf(infile, "%*[^\n]");
    getc(infile);
    for (j = 0; j < namelength; j++) {
      if (eoln(infile))
	nayms[i][j] = ' ';
      else
	nayms[i][j] = getc(infile);
    }
    for (j = 0; j < numsp; j++) {
      skipit = ((lower && j + 1 >= i + 1) || (upper && j + 1 <= i + 1));
      skipother = ((lower && i + 1 >= j + 1) || (upper && i + 1 <= j + 1));
      if (!skipit) {
	if (eoln(infile)) {
	  fscanf(infile, "%*[^\n]");
	  getc(infile);
	}
	fscanf(infile, "%lf", &x[i][j]);
	if (replicates) {
	  if (eoln(infile)) {
	    fscanf(infile, "%*[^\n]");
	    getc(infile);
	  }
	  fscanf(infile, "%hd", &reps[i][j]);
	} else
	  reps[i][j] = 1;
      }
      if (!skipit && skipother) {
          x[j][i] = x[i][j];
          reps[j][i] = reps[i][j];
        }
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  if (!printdata)
    return;
  for (i = 0; i < numsp; i++) {
    for (j = 0; j < namelength; j++)
      putc(nayms[i][j], outfile);
    putc(' ', outfile);
    for (j = 1; j <= numsp; j++) {
      fprintf(outfile, "%10.5f", x[i][j - 1]);
      if (replicates)
        fprintf(outfile, " (%3hd)", reps[i][j - 1]);
      if (j % columns == 0 && j < numsp) {
        putc('\n', outfile);
        for (k = 1; k <= namelength + 1; k++)
          putc(' ', outfile);
      }
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* getdata */

void hookup(p, q)
node *p, *q;
{
  /* hook together two nodes */
  p->back = q;
  q->back = p;
}  /* hookup */



void secondtraverse(q, y, nx,sum)
node *q;
double y;
short *nx;          /* comes from firsttraverse                              */
double *sum;        /* comes from evaluate via firsttraverse                 */
{
  /* from each of those places go back to all others */
  double z=0.0, TEMP=0.0;

  z = y + q->v;
  if (q->tip) {
    TEMP = q->d[(*nx) - 1] - z;
    *sum += q->w[(*nx) - 1] * (TEMP * TEMP);
  } else {
    secondtraverse(q->next->back, z, nx, sum);
    secondtraverse(q->next->next->back, z, nx,sum);
  }
}  /* secondtraverse */

void firsttraverse(p, nx,sum)
node *p;
short *nx;
double *sum;
{
  /* go through tree calculating branch lengths */
  if (p->tip) {
    *nx = p->number;
  secondtraverse(p->back, 0.0, nx,sum);

  } else {
    firsttraverse(p->next->back, nx,sum);
    firsttraverse(p->next->next->back, nx,sum);
  }
}  /* firsttraverse */

double evaluate(t)
tree *t;
{
  double sum=0.0;
  short nx=0;
  /* evaluate likelihood of a tree */
  firsttraverse(t->start->back->back,&nx,&sum);
  firsttraverse(t->start->back, &nx,&sum);
  if (replicates && (lower || upper))
    sum /= 2;
  t->likelihood = -sum;
  return (-sum);
}  /* evaluate */

void nudists(x, y)
node *x, *y;
{
  /* compute distance between an interior node and tips */
  short nq=0, nr=0, nx=0, ny=0;
  double dil=0, djl=0, wil=0, wjl=0, vi=0, vj=0;
  node *qprime, *rprime;

  qprime = x->next;
  rprime = qprime->next->back;
  qprime = qprime->back;
  ny = y->number;
  dil = qprime->d[ny - 1];
  djl = rprime->d[ny - 1];
  wil = qprime->w[ny - 1];
  wjl = rprime->w[ny - 1];
  vi = qprime->v;
  vj = rprime->v;
  x->w[ny - 1] = wil + wjl;
  if (wil + wjl <= 0.0)
    x->d[ny - 1] = 0.0;
  else
    x->d[ny - 1] = ((dil - vi) * wil + (djl - vj) * wjl) / (wil + wjl);
  nx = x->number;
  nq = qprime->number;
  nr = rprime->number;
  dil = y->d[nq - 1];
  djl = y->d[nr - 1];
  wil = y->w[nq - 1];
  wjl = y->w[nr - 1];
  y->w[nx - 1] = wil + wjl;
  if (wil + wjl <= 0.0)
    y->d[nx - 1] = 0.0;
  else
    y->d[nx - 1] = ((dil - vi) * wil + (djl - vj) * wjl) / (wil + wjl);
}  /* nudists */


void makedists(p)
node *p;
{
  /* compute distances among three neighbors of a node */
  short i=0, nr=0, ns=0;
  node *q, *r, *s;


  r = p->back;
  nr = r->number;
  for (i = 1; i <= 3; i++) {
    q = p->next;
    s = q->back;
    ns = s->number;

    if (s->w[nr - 1] + r->w[ns - 1] <= 0.0)
      p->dist = 0.0;
    else
      p->dist = (s->w[nr - 1] * s->d[nr - 1] + r->w[ns - 1] * r->d[ns - 1]) /
                (s->w[nr - 1] + r->w[ns - 1]);
    p = q;
    r = s;
    nr = ns;
  }
}  /* makedists */

void makebigv(p)
node *p;
{
  /* make new branch length */
  short i=0;
  node *temp, *q, *r;

  q = p->next;
  r = q->next;
  for (i = 1; i <= 3; i++) {
    if (p->iter) {
      p->v = (p->dist + r->dist - q->dist) / 2.0;
      p->back->v = p->v;
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
  short i=0, j=0, n=0, nq=0, nr=0, ntemp=0;
  double wq=0.0, wr=0.0;

  q = p->next;
  r = q->next;
  n = p->back->number;
  nq = q->back->number;
  nr = r->back->number;
  for (i = 1; i <= smoothings; i++) {
    for (j = 1; j <= 3; j++) {
      if (p->iter) {
        wr = r->back->w[n - 1] + p->back->w[nr - 1];
        wq = q->back->w[n - 1] + p->back->w[nq - 1];
        if (wr + wq <= 0.0 && !negallowed)
          p->v = 0.0;
        else
          p->v = ((p->dist - q->v) * wq + (r->dist - r->v) * wr) / (wr + wq);
        if (p->v < 0 && !negallowed)
          p->v = 0.0;
        p->back->v = p->v;
      }
      temp = p;
      p = q;
      q = r;
      r = temp;
      ntemp = n;
      n = nq;
      nq = nr;
      nr = ntemp;
    }
  }
}  /* correctv */

void alter(x, y)
node *x, *y;
{
  /* traverse updating these views */
  nudists(x, y);
  if (!y->tip) {
    alter(x, y->next->back);
    alter(x, y->next->next->back);
  }
}  /* alter */

void nuview(p)
node *p;
{
  /* renew information about subtrees */
  short i=0;
  node *q, *r, *pprime, *temp;

  q = p->next;
  r = q->next;
  for (i = 1; i <= 3; i++) {
    temp = p;
    pprime = p->back;
    alter(p, pprime);
    p = q;
    q = r;
    r = temp;
  }
}  /* nuview */

void update(p)
node *p;
{
  /* update branch lengths around a node */

  if (p->tip)
    return;
  makedists(p);
  if (p->iter || p->next->iter || p->next->next->iter) {
    makebigv(p);
    correctv(p);
  }
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

void filltraverse(pb, qb,contin)
node *pb, *qb;
boolean contin;
{
  if (qb->tip)
    return;
  if (contin) {
    filltraverse(pb, qb->next->back,contin);
    filltraverse(pb, qb->next->next->back,contin);
    nudists(qb, pb);
    return;
  }
  if (!qb->next->back->tip)
    nudists(qb->next->back, pb);
  if (!qb->next->next->back->tip)
    nudists(qb->next->next->back, pb);
}  /* filltraverse */

void fillin(pa, qa,contin)
node *pa, *qa;
boolean contin;
{
  if (!pa->tip) {
    fillin(pa->next->back, qa, contin);
    fillin(pa->next->next->back, qa, contin);
  }
  filltraverse(pa, qa, contin);
}  /* fillin */

void insert_(p, q, contin_)
node *p, *q;
boolean contin_;
{
  /* put p and q together and iterate info. on resulting tree */
  short i=0;
  double x=0.0;
  hookup(p->next->next, q->back);
  hookup(p->next, q);
  x = q->v / 2.0;
  p->v = 0.0;
  p->back->v = 0.0;
  p->next->v = x;
  p->next->back->v = x;
  p->next->next->back->v = x;
  p->next->next->v = x;
  fillin(p->back, p, contin_);
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
  memcpy(d->d, c->d, numsp2*sizeof(double));
  memcpy(d->w, c->w, numsp2*sizeof(double));
  d->v = c->v;
  d->iter = c->iter;
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
      else if (a->nodep[i]->back
                 == a->nodep[a->nodep[i]->back->number - 1]->next)
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->number - 1]->next;
      else
        b->nodep[i]->back
          = b->nodep[a->nodep[i]->back->number - 1]->next->next;
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

void setuptip(m,t)
short m;
tree *t;
{
  /* initialize branch lengths and views in a tip */
  short i=0;
  intvector n=(short *)Malloc(numsp * sizeof(short)); 
  node *WITH;

  WITH = t->nodep[m - 1];
  memcpy(WITH->d, x[m - 1], (numsp2 * sizeof(double)));
  memcpy(n, reps[m - 1], (numsp * sizeof(short)));

  for (i = 0; i < numsp; i++) {
    if (i + 1 != m && n[i] > 0) {
      if (WITH->d[i] < epsilon)
        WITH->d[i] = epsilon;
      WITH->w[i] = n[i] / exp(power * log(WITH->d[i]));
    } else {
      WITH->w[m - 1] = 0.0;
      WITH->d[m - 1] = 0.0;
    }
  }
  for (i = numsp; i < numsp2; i++) {
    WITH->w[i] = 1.0;
    WITH->d[i] = 0.0;
  }
  WITH->number = m;
  memcpy(WITH->nayme, nayms[m - 1], sizeof(naym));
  if (WITH->iter) WITH->v = 0.0;
  free(n);
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

void buildsimpletree(t,nextsp)
tree *t;
short nextsp;
{
  /* make and initialize a three-species tree */
  setuptip(enterorder[0], t);
  setuptip(enterorder[1], t);
  hookup(t->nodep[enterorder[0] - 1], t->nodep[enterorder[1] - 1]);
  buildnewtip(enterorder[2], t, nextsp);
  insert_(t->nodep[enterorder[2] - 1]->back, t->nodep[enterorder[0] - 1],
          false);
}  /* buildsimpletree */

void addtraverse(p, q, contin, numtrees,succeeded)
node *p, *q;
boolean contin,*succeeded;
short *numtrees;
{
 /* traverse through a tree, finding best place to add p */
  insert_(p, q, true);
  (*numtrees)++;
  if (evaluate(&curtree) > bestree.likelihood){
    copy_(&curtree, &bestree);
    (*succeeded)=true;
  }
  copy_(&priortree, &curtree);
  if (!q->tip && contin) {
    addtraverse(p, q->next->back, contin,numtrees,succeeded);
    addtraverse(p, q->next->next->back, contin,numtrees,succeeded);
  }
}  /* addtraverse */


void re_move(p, q)
node **p, **q;
{
  /* re_move p and record in q where it was */
  *q = (*p)->next->back;
  hookup(*q, (*p)->next->next->back);
  (*p)->next->back = NULL;
  (*p)->next->next->back = NULL;
  update(*q);
  update((*q)->back);
}  /* re_move */

void rearrange(p, numtrees,nextsp,succeeded)
node *p;
short *numtrees,*nextsp;
boolean *succeeded;
{
  node *q, *r;
  if (!p->tip && !p->back->tip) {
    r = p->next->next;
    re_move(&r, &q);
    copy_(&curtree, &priortree);
    addtraverse(r, q->next->back, global && ((*nextsp) == numsp),
                numtrees,succeeded);
    addtraverse(r, q->next->next->back, global && ((*nextsp) == numsp),
                numtrees,succeeded);
    copy_(&bestree, &curtree);
    if (global && ((*nextsp) == numsp))
      putchar('.');
    if (global && ((*nextsp) == numsp) && !(*succeeded)) {
      if (r->back->tip) {
        r = r->next->next;
        re_move(&r, &q);
        q = q->back;
        copy_(&curtree, &priortree);
        if (!q->tip) {
          addtraverse(r, q->next->back, true, numtrees,succeeded);
          addtraverse(r, q->next->next->back, true, numtrees,succeeded);
        }
        q = q->back;
        if (!q->tip) {
          addtraverse(r, q->next->back, true, numtrees,succeeded);
          addtraverse(r, q->next->next->back, true, numtrees,succeeded);
        }
        copy_(&bestree, &curtree);
      }
    }
  }
  if (!p->tip) {
    rearrange(p->next->back, numtrees,nextsp,succeeded);
    rearrange(p->next->next->back, numtrees,nextsp,succeeded);
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
    p->xcoord = (short)((double)over * lengthsum + 0.5);
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
    if (lengthsum > *tipmax){
      *tipmax = lengthsum;}
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
  p->xcoord = (short)((double)over * lengthsum + 0.5);
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
  short n=0, j=0;
  boolean extra=false, trif=false;
  node *r, *first, *last;
  boolean done=false;

  p = curtree.start->back;
  q = curtree.start->back;
  extra = false;
  trif = false;
  if (i == p->ycoord && p == curtree.start->back) {
    if (p->number - numsp >= 10)
      fprintf(outfile, "-%2hd", p->number - numsp);
    else
      fprintf(outfile, "--%hd", p->number - numsp);
    extra = true;
    trif = true;
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
    n = (short)(scale * (double)(q->xcoord - p->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if (q->ycoord == i && !done) {
      if (p->ycoord != q->ycoord)
        putc('+', outfile);
      if (trif) {
        n++;
        trif = false;
      }
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
        trif = false;
      }
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
  short i=0;
  short tipy=0;
  double scale=0.0,tipmax=0.0,divisor=0.0;

  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  tipmax = 0.0;
  coordinates(curtree.start->back, 0.0, &tipy,&tipmax);
  divisor = ((short)(tipmax + 1.000));
  scale = 1.0 / (double)divisor;
  for (i = 1; i <= (tipy - down); i++)
    drawline(i, scale);
  putc('\n', outfile);
}  /* printree */


void treeout(p)
node *p;
{
  /* write out file with representation of final tree */
  short i=0, n=0, w=0;
  Char c;
  double x=0.0;

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
      treeout(p->back);
    }
    putc(')', treefile);
    col++;
  }
  x = p->v;
  if (x > 0.0)
    w = (short)(0.43429445222 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (short)(0.43429445222 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if (p == curtree.start->back)
    fprintf(treefile, ";\n");
  else {
    fprintf(treefile, ":%*.5f", w + 7, x);
    col += w + 8;
  }
}  /* treeout */

void describe(p)
node *p;
{
  /* print out information for one branch */
  short i=0;
  node *q;

  q = p->back;
  fprintf(outfile, "%4hd          ", q->number - numsp);
  if (p->tip) {
    for (i = 0; i < namelength; i++)
      putc(p->nayme[i], outfile);
  } else
    fprintf(outfile, "%4hd      ", p->number - numsp);
  fprintf(outfile, "%15.5f\n", q->v);
  if (!p->tip) {
    describe(p->next->back);
    describe(p->next->next->back);
  }
}  /* describe */

void summarize(numtrees)
short numtrees;
{
  /* print out branch lengths etc. */
  short i, j, totalnum;

  fprintf(outfile, "\nremember:");
  if (outgropt)
    fprintf(outfile, " (although rooted by outgroup)");
  fprintf(outfile, " this is an unrooted tree!\n\n");
  fprintf(outfile, "Sum of squares = %11.5f\n\n", -curtree.likelihood);
  if (power == 2.0) {
    totalnum = 0;
    for (i = 1; i <= nums; i++) {
      for (j = 1; j <= nums; j++) {
        if (i != j)
          totalnum += reps[i - 1][j - 1];
      }
    }
    fprintf(outfile, "Average percent standard deviation = ");
    fprintf(outfile, "%11.5f\n\n",
            100 * sqrt(curtree.likelihood / (2 - totalnum)));
  }
  if (!usertree)
    fprintf(outfile, "examined %4hd trees\n\n", numtrees);
  fprintf(outfile, "Between        And            Length\n");
  fprintf(outfile, "-------        ---            ------\n");
  describe(curtree.start->back->next->back);
  describe(curtree.start->back->next->next->back);
  describe(curtree.start);
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
  while (!(done)) {
    if (c == ',') {
      if (ch == '(' || ch == ')' || ch == ':' || ch == ';') {
        printf(
             "\nERROR IN USER TREE: UNMATCHED PARENTHESIS OR MISSING COMMA\n");
        printf(" OR NOT TRIFURCATED BASE\n");
	exit(-1);
      } else if (ch == ',')
        done = true;
    } else if (c == ')') {
      if (ch == '(' || ch == ',' || ch == ':' || ch == ';') {
        printf("\nERROR IN USER TREE:");
	printf(" UNMATCHED PARENTHESIS OR NON-BIFURCATED NODE\n");
	exit(-1);
      } else if (ch == ')') {
        (*rparens)++;
        if (*lparens > 0 && *lparens == *rparens ) {
          if (*lparens == numsp - 2) {
            if (eoln(infile)) {
              fscanf(infile, "%*[^\n]");
              getc(infile);
            }
            ch = getc(infile);
            if (ch != ';') {
              printf("\nERROR IN USER TREE:");
	      printf(" UNMATCHED PARENTHESIS OR MISSING SEMICOLON\n");
	      exit(-1);
            }
          }
        }
	done = true;
      }
    }
    if ((done && ch == ')') || !(done)) {
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
  short digit, ordzero;
  double valyew, divisor;
  boolean pointread;

  ordzero = '0';
  pointread = false;
  valyew = 0.0;
  divisor = 1.0;
  getch(&ch);
  digit = ch - ordzero;
  while (((unsigned short)digit <= 9) || (ch == '.')) {
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

#undef point

void addelement(p,nextnode,lparens,rparens,names,nolengths)
node *p;
short *nextnode,*lparens,*rparens;
boolean *names;                             /* a boolean array               */
boolean *nolengths;
{
                                            /* add one node to the user tree */
  node *q;
  short i=0, n=0;
  boolean found=false;
  Char str[namelength];

  strcpy(str,"");
  do {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    ch = getc(infile);
  } while (ch == ' ');
  if (ch == '(') {
    (*lparens)++;
    (*nextnode)++;
    q = curtree.nodep[(*nextnode) - 1];
    hookup(p, q);
    addelement(q->next,nextnode,lparens,rparens,names,nolengths);
    findch(',', lparens,rparens);
    addelement(q->next->next,nextnode,lparens,rparens,names,nolengths);
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
      printf("Cannot find species: ");
      for (i = 0; i < namelength; i++)
        putchar(str[i]);
      putchar('\n');
    }
    nums++;
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
  short i=0;
  short nextnode=0,lparens=0,rparens=0;
  boolean *names, nolengths;

  nums = 0;
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
    addelement(p, &nextnode,&lparens,&rparens,names,&nolengths);
    p = p->next;
    findch(',', &lparens,&rparens);
  }
  addelement(p, &nextnode,&lparens,&rparens,names,&nolengths);
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
  short i, j;

  for (i = 1; i <= 3; i++) {
    for (j = 0; j < numsp2; j++) {
      p->w[j] = 1.0;
      p->d[j] = 0.0;
    }
    p = p->next;
  }
  if (p->iter)
    p->v = 1.0;
  if (p->back->iter)
    p->back->v = 1.0;
}  /* nodeinit */

void initrav(p)
node *p;
{
  /* traverse to initialize */
  if (p->tip)
    return;
  nodeinit(p);
  initrav(p->next->back);
  initrav(p->next->next->back);
}  /* initrav */

void treevaluate()
{
  /* evaluate user-defined tree, iterating branch lengths */
  short i;
  double dummy;

  for (i = 1; i <= numsp; i++)
    setuptip(i, &curtree);
  initrav(curtree.start);
  if (curtree.start->back != NULL) {
    initrav(curtree.start->back);
    for (i = 1; i <= smoothings * 4; i++)
      smooth(curtree.start->back);
  }
  dummy = evaluate(&curtree);
}  /* treevaluate */


void  maketree()
{
  /* contruct the tree */
  short nextsp,numtrees;
  boolean succeeded=false;
  short i,j,k,which;

  if (usertree) {
    getdata();
    setuptree(&curtree);
    for (which = 1; which <= numsp; which++)
      setuptip(which, &curtree);
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%hd%*[^\n]", &numtrees);
    getc(infile);
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }
    which = 1;
    while (which <= numtrees) {
      treeread();
      curtree.start = curtree.nodep[outgrno - 1];
      treevaluate();
      printree();
      summarize(numtrees);
      which++;
    }
  } else {
    if (jumb == 1) {
      getdata();
      setuptree(&curtree);
      setuptree(&priortree);
      setuptree(&bestree);
      if (njumble > 1) setuptree(&bestree2);
    }
    for (i = 1; i <= numsp; i++)
      enterorder[i - 1] = i;
    if (jumble) {
      for (i = 0; i < numsp; i++) {
        j = (short)(randum(seed) * (double)numsp) + 1;
        k = enterorder[j - 1];
        enterorder[j - 1] = enterorder[i];
        enterorder[i] = k;
      }
    }
    nextsp = 3;
    buildsimpletree(&curtree, nextsp);
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
      nums = nextsp;
      buildnewtip(enterorder[nextsp - 1], &curtree, nextsp);
      copy_(&curtree, &priortree);
      bestree.likelihood = -99999.0;
      addtraverse(curtree.nodep[enterorder[nextsp - 1] - 1]->back,
                  curtree.start->back, true, &numtrees,&succeeded);
      copy_(&bestree, &curtree);
      if (progress) {
        printf("   ");
        for (j = 0; j < namelength; j++)
          putchar(nayms[enterorder[nextsp - 1] - 1][j]);
        putchar('\n');
      }
      if (global && nextsp == numsp) {
        if (progress) {
          printf("Doing global rearrangements\n");
          printf("  !");
          for (j = 1; j <= (numsp - 2); j++)
            putchar('-');
          printf("!\n");
          printf("   ");
        }
      }
      succeeded = true;
      while (succeeded) {
        succeeded = false;
        rearrange(curtree.start->back,
                  &numtrees,&nextsp,&succeeded);
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
  if (jumb == njumble && progress) {
    printf("\nOutput written to output file\n\n");
    if (trout)
      printf("Tree also written onto file\n");
    putchar('\n');
  }
}  /* maketree */


main(argc, argv)
int argc;
Char *argv[];
{
  int i;
  char infilename[100],outfilename[100],trfilename[100];
#ifdef MAC
  macsetup("Fitch","");
  argv[0] = "Fitch";
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
  nayms = (naym *)Malloc(numsp*sizeof(naym));
  enterorder = (short *)Malloc(numsp*sizeof(short));
  for (i=0;i<numsp;++i){
    enterorder[i]=0;}
  for (ith = 1; ith <= datasets; ith++) {
    if (datasets > 1) {
      fprintf(outfile, "Data set # %hd:\n\n",ith);
      if (progress)
        printf("\nData set # %hd:\n",ith);
    }
    getinput();
    for (jumb = 1; jumb <= njumble; jumb++)
        maketree();
    firstset = false;
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
  }
  if (trout)
    FClose(treefile);
  FClose(outfile);
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
long  x;
{
MALLOCRETURN *mem;
mem = (MALLOCRETURN *)calloc(1,x);
if (!mem)
  memerror();
else
  return (MALLOCRETURN *)mem;
}

