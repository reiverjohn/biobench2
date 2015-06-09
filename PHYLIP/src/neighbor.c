#include "phylip.h"

/* version 3.52c. (c) Copyright 1993 by Joseph Felsenstein and Mary Kuhner.
   Written by Mary Kuhner, Jon Yamato, Joseph Felsenstein, Akiko Fuseki,
   Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define namelength      10   /* number of characters max. in species name */
#define down            2
#define over            60

#define ibmpc0          false
#define ansi0           true
#define vt520           false

typedef double *vector;
typedef long *intvector;
typedef Char naym[namelength];
typedef long longer[6];

typedef struct node {
  struct node *next, *back;
  boolean tip;
  long number;
  naym nayme;
  double v;
  long xcoord, ycoord, ymin, ymax;
} node;

typedef struct tree {
  node **nodep;
  node *start;
} tree;


Static FILE *infile, *outfile, *treefile;
Static long numsp, numsp1, numsp2, inseed, outgrno, col, datasets, ith, j;
Static vector *x;
Static intvector *reps;
Static naym *names;
Static boolean jumble, lower, upper, outgropt, replicates, trout,
               printdata, progress, treeprint, mulsets, ibmpc, vt52, ansi,
               njoin;
Static tree curtree;
Static longer seed;
Static long *enterorder;

/* Local variables for maketree, propagated globally for C version: */
node **cluster;


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



double randum(seed)
long *seed;
{
  /* random number generator -- slow but machine independent */
  long i, j, k, sum;
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
    *ch = isupper(*ch) ? (*ch) : toupper(*ch);
}  /* uppercase */


void getnums()
{
  /* read species numbers or number of bootstrap data sets */
  fscanf(infile, "%ld", &numsp);
  fprintf(outfile, "\n%4ld Populations\n", numsp);
  numsp1 = numsp + 1;
  numsp2 = numsp * 2 - 2;
}  /* getnums */

void getoptions()
{
  /* interactively set options */
  long i, inseed0;
  Char ch;
  boolean done1;

  fprintf(outfile, "\nNeighbor-Joining/UPGMA method version %s\n\n",VERSION);
  putchar('\n');
  jumble = false;
  lower = false;
  outgrno = 1;
  outgropt = false;
  replicates = false;
  trout = true;
  upper = false;
  printdata = false;
  progress = true;
  treeprint = true;
  njoin = true;
  for(;;) {
    printf(ansi ? "\033[2J\033[H" :
           vt52 ? "\033E\033H"    : "\n");
    printf("\nNeighbor-Joining/UPGMA method version 3.5\n\n");
    printf("Settings for this run:\n");
    printf("  N       Neighbor-joining or UPGMA tree?  %s\n",
           (njoin ? "Neighbor-joining" : "UPGMA"));
    if (njoin) {
      printf("  O                        Outgroup root?");
      if (outgropt)
        printf("  Yes, at species number%3ld\n", outgrno);
      else
        printf("  No, use as outgroup species%3ld\n", outgrno);
    }
    printf("  L         Lower-triangular data matrix?  %s\n",
           (lower ? "Yes" : "No"));
    printf("  R         Upper-triangular data matrix?  %s\n",
           (upper ? "Yes" : "No"));
    printf("  S                        Subreplicates?  %s\n",
           (replicates ? "Yes" : "No"));
    printf("  J     Randomize input order of species?");
    if (jumble)
      printf("  Yes (random number seed =%8ld)\n", inseed0);
    else
      printf("  No. Use input order\n");
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  0   Terminal type (IBM PC, VT52, ANSI)?  %s\n",
           (ibmpc ? "IBM PC" : ansi ? "ANSI" : vt52 ? "VT52" : "(none"));

    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
           (treeprint ? "Yes" : "No"));
    printf("  4       Write out trees onto tree file?  %s\n",
           (trout ? "Yes" : "No"));
    printf( "\nAre these settings correct?");
    printf(" (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if  (ch == 'Y')
      break;
    if (strchr("NJOULRSM01234",ch) != NULL){
      switch (ch) {
	
      case 'J':
	jumble = !jumble;
	if (jumble) {
	  do {
	    printf("Random number seed (must be odd)?\n");
	    scanf("%ld%*[^\n]", &inseed);
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
	}
	break;
	
      case 'L':
	lower = !lower;
	break;
	
      case 'O':
	outgropt = !outgropt;
	if (outgropt) {
	  done1 = true;
	  do {
	    printf("Type number of the outgroup:\n");
	    scanf("%ld%*[^\n]", &outgrno);
	    getchar();
	    done1 = (outgrno >= 1 && outgrno <= numsp);
	    if (!done1) {
	      printf("BAD OUTGROUP NUMBER: %4ld\n", outgrno);
	      printf("  Must be in range 1 -%2ld\n", numsp);
	    }
	  } while (done1 != true);
	}
        else
          outgrno = 1;
	break;
	
      case 'R':
	upper = !upper;
	break;
	
      case 'S':
	replicates = !replicates;
	break;
	
      case 'N':
	njoin = !njoin;
	break;
	
      case 'M':
	mulsets = !mulsets;
	if (mulsets) {
	  done1 = false;
	  do {
	    printf("How many data sets?\n");
	    scanf("%ld%*[^\n]", &datasets);
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
}  /* getoptions */


void doinit()
{
  /* initializes variables */
  getnums();
  getoptions();
}  /* doinit */


void inputoptions()
{
  /* read options information */
  Char ch;
  long cursp;

  if (ith != 1) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%ld", &cursp);
    if (cursp != numsp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4ld\n",
             ith);
      exit(-1);
    }
  }
  while (!eoln(infile)) {
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch != ' ') {
      printf("BAD OPTION CHARACTER: %c\n", ch);
      exit(-1);
    }
  }
  putc('\n', outfile);
  if (njoin)
    fprintf(outfile, " Neighbor-joining method\n");
  else
    fprintf(outfile, " UPGMA method\n");
  fprintf(outfile, "\n Negative branch lengths allowed\n\n");
}  /* inputoptions */


void getinput()
{
  /* reads the input data */
  inputoptions();
}  /* getinput */




void setuptree(a)
tree *a;
{
  /* initialize a tree */
  long i, j;
  node *p, *q;

  a->nodep = (node **)Malloc((numsp2+1)*sizeof(node *));
  for (i = 0; i < numsp; i++) {
    a->nodep[i] = (node *)Malloc(sizeof(node));
    a->nodep[i]->tip = true;
    a->nodep[i]->number = i + 1;
    memcpy(a->nodep[i]->nayme, names[i], sizeof(naym));
    a->nodep[i]->v = 0.0;
  }
  for (i = numsp1; i <= numsp2; i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc(sizeof(node));
      p->tip = false;
      p->number = i;
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    a->nodep[i - 1] = p;
  }
  q = NULL;
  for (j = 1; j <= 2; j++) {
    p = (node *)Malloc(sizeof(node));
    p->tip = false;
    p->number = numsp2 + 1;
    p->next = q;
    q = p;
  }
  p->next->next = p;
  a->nodep[numsp2] = p;
  a->start = a->nodep[0];
}  /* setuptree */

void getdata()
{
  /* read in distance matrix */
  long i, j, k, columns;
  boolean skipit, skipother;

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
	names[i][j] = ' ';
      else {
	names[i][j] = getc(infile);
	if (names[i][j] == '\n')
	  names[i][j] = ' ';
      }
    }
    for (j = 0; j < numsp; j++) {
      skipit = ((lower && j + 1 >= i + 1) || (upper && j + 1 <= i + 1));
      skipother = ((lower && i + 1 >= j + 1 )|| (upper && i + 1 <= j + 1));
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
	  fscanf(infile, "%ld", &reps[i][j]);
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
      putc(names[i][j], outfile);
    putc(' ', outfile);
    for (j = 1; j <= numsp; j++) {
      fprintf(outfile, "%10.5f", x[i][j - 1]);
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

Local Void hookup(p, q)
node *p, *q;
{
  /* hook together two nodes */
  p->back = q;
  q->back = p;
}  /* hookup */


void coordinates(p, lengthsum,tipy,tipmax)
node *p;
double lengthsum;
long *tipy;
double *tipmax;
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;

  if (p->tip) {
    p->xcoord = (long)(over * lengthsum + 0.5);
    p->ycoord = (*tipy);
    p->ymin = (*tipy);
    p->ymax = (*tipy);
    (*tipy) += down;
    if (lengthsum > (*tipmax))
      (*tipmax) = lengthsum;
    return;
  }
  q = p->next;
  do {
    coordinates(q->back, lengthsum + q->v, tipy,tipmax);
    q = q->next;
  } while ((p == curtree.start || p != q)
        && (p != curtree.start || p->next != q));
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (long)(over * lengthsum + 0.5);
  if (p == curtree.start) {
    if (njoin)
      p->ycoord = p->next->next->back->ycoord;
    else
      p->ycoord = (p->back->ycoord + p->next->back->ycoord) / 2;
  } else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */

void drawline(i, scale)
long i;
double scale;
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q;
  long n, j;
  boolean extra, trif;
  node *r, *first, *last;
  boolean done;

  p = curtree.start;
  q = curtree.start;
  extra = false;
  trif = false;
  if (i == p->ycoord && p == curtree.start) {
    if (p->number - numsp >= 10)
      fprintf(outfile, "-%2ld", p->number - numsp);
    else
      fprintf(outfile, "--%ld", p->number - numsp);
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
      } while (!(done || (p != curtree.start && r == p) ||
                 (p == curtree.start && r == p->next)));
      first = p->next->back;
      r = p;
      while (r->next != p)
        r = r->next;
      last = r->back;
      if (p == curtree.start)
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
      if (trif) {
        n++;
        trif = false;
      }
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        if (q->number - numsp >= 10)
          fprintf(outfile, "%2ld", q->number - numsp);
        else
          fprintf(outfile, "-%ld", q->number - numsp);
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

  long tipy = 1;
  double scale, tipmax = 0.0;
  long i;

  if (!treeprint)
    return;
  putc('\n', outfile);
  coordinates(curtree.start, 0.0, &tipy,&tipmax);
  scale = 1.0 / (long)(tipmax + 1.000);
  for (i = 1; i <= tipy - down; i++)
    drawline(i, scale);
  putc('\n', outfile);
}  /* printree */


void describe(p)
node *p;
{
  /* print out information for one branch */
  long i;
  node *q;

  q = p->back;
  fprintf(outfile, "%4ld          ", q->number - numsp);
  if (p->tip) {
    for (i = 0; i < namelength; i++)
      putc(p->nayme[i], outfile);
  } else
    fprintf(outfile, "%4ld      ", p->number - numsp);
  fprintf(outfile, "%15.5f\n", q->v);
  if (!p->tip) {
    describe(p->next->back);
    describe(p->next->next->back);
  }
}  /* describe */

void summarize()
{
  /* print out branch lengths etc. */
  putc('\n', outfile);
  if (njoin) {
    fprintf(outfile, "remember:");
    if (outgropt)
      fprintf(outfile, " (although rooted by outgroup)");
    fprintf(outfile, " this is an unrooted tree!\n");
  }
  fprintf(outfile, "\nBetween        And            Length\n");
  fprintf(outfile, "-------        ---            ------\n");
  describe(curtree.start->next->back);
  describe(curtree.start->next->next->back);
  if (njoin)
    describe(curtree.start->back);
  fprintf(outfile, "\n\n");
}  /* summarize */

void treeout(p)
node *p;
{
  /* write out file with representation of final tree */
  long i, n, w;
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
    if (p == curtree.start && njoin) {
      putc(',', treefile);
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
  if (p == curtree.start)
    fprintf(treefile, ";\n");
  else {
    fprintf(treefile, ":%*.5f", (int)(w + 7), x);
    col += w + 8;
  }
}  /* treeout */

void nodelabel(isnode)
boolean isnode;
{
  if (isnode)
    printf("NODE");
  else
    printf("OTU ");
}  /* nodelabel */

Local Void jointree()
{
  /* calculate the tree */
  long nc, nextnode, mini, minj, i, j, ia, ja, ii, jj, nude, iter;
  double diq, djq, dij, fotu2, total, tmin, dio, djo, bi, bj, bk, dmin, da;
  long el[3];
  vector av;
  intvector oc;

  for (i = 0; i <= numsp - 2; i++) {
    for (j = i + 1; j < numsp; j++) {
      da = (x[i][j] + x[j][i]) / 2.0;
      x[i][j] = da;
      x[j][i] = da;
    }
  }
  /* First initialization */
  fotu2 = numsp - 2.0;
  nextnode = numsp + 1;
  av = (vector)Malloc(numsp*sizeof(double));
  oc = (intvector)Malloc(numsp*sizeof(long));
  for (i = 0; i < numsp; i++) {
    av[i] = 0.0;
    oc[i] = 1;
  }
  /* Enter the main cycle */
  if (njoin)
    iter = numsp - 3;
  else
    iter = numsp - 1;
  for (nc = 1; nc <= iter; nc++) {
    for (j = 2; j <= numsp; j++) {
      for (i = 0; i <= j - 2; i++)
        x[j - 1][i] = x[i][j - 1];
    }
    tmin = 99999.0;
    /* Compute sij and minimize */
    for (ja = 2; ja <= numsp; ja++) {
      jj = enterorder[ja - 1];
      if (cluster[jj - 1] != NULL) {
        for (ia = 0; ia <= ja - 2; ia++) {
          ii = enterorder[ia];
          if (cluster[ii - 1] != NULL) {
            if (njoin) {
              diq = 0.0;
              djq = 0.0;
              dij = x[ii - 1][jj - 1];
              for (i = 0; i < numsp; i++) {
                diq += x[i][ii - 1];
                djq += x[i][jj - 1];
              }
              total = fotu2 * dij - diq - djq;
            } else
              total = x[ii - 1][jj - 1];
            if (total < tmin) {
              tmin = total;
              mini = ii;
              minj = jj;
            }
          }
        }
      }
    }
    /* compute lengths and print */
    if (njoin) {
      dio = 0.0;
      djo = 0.0;
      for (i = 0; i < numsp; i++) {
        dio += x[i][mini - 1];
        djo += x[i][minj - 1];
      }
      dmin = x[mini - 1][minj - 1];
      dio = (dio - dmin) / fotu2;
      djo = (djo - dmin) / fotu2;
      bi = (dmin + dio - djo) * 0.5;
      bj = dmin - bi;
      bi -= av[mini - 1];
      bj -= av[minj - 1];
    } else {
      bi = x[mini - 1][minj - 1] / 2.0 - av[mini - 1];
      bj = x[mini - 1][minj - 1] / 2.0 - av[minj - 1];
      av[mini - 1] += bi;
    }
    if (progress) {
      printf("CYCLE %3ld: ", iter - nc + 1);
      if (njoin)
        nodelabel(av[mini - 1] > 0.0);
      else
        nodelabel(oc[mini - 1] > 1.0);
      printf("%3ld (%10.5f) JOINS ", mini, bi);
      if (njoin)
        nodelabel(av[minj - 1] > 0.0);
      else
        nodelabel(oc[minj - 1] > 1.0);
      printf("%3ld (%10.5f)\n", minj, bj);
    }
    hookup(curtree.nodep[nextnode - 1]->next, cluster[mini - 1]);
    hookup(curtree.nodep[nextnode - 1]->next->next, cluster[minj - 1]);
    cluster[mini - 1]->v = bi;
    cluster[minj - 1]->v = bj;
    cluster[mini - 1]->back->v = bi;
    cluster[minj - 1]->back->v = bj;
    cluster[mini - 1] = curtree.nodep[nextnode - 1];
    cluster[minj - 1] = NULL;
    nextnode++;
    if (njoin)
      av[mini - 1] = dmin * 0.5;
    /* re-initialization */
    fotu2 -= 1.0;
    for (j = 0; j < numsp; j++) {
      if (cluster[j] != NULL) {
        if (njoin) {
          da = (x[mini - 1][j] + x[minj - 1][j]) * 0.5;
          if (mini - j - 1 < 0)
            x[mini - 1][j] = da;
          if (mini - j - 1 > 0)
            x[j][mini - 1] = da;
        } else {
          da = x[mini - 1][j] * oc[mini - 1] + x[minj - 1][j] * oc[minj - 1];
          da /= oc[mini - 1] + oc[minj - 1];
          x[mini - 1][j] = da;
          x[j][mini - 1] = da;
        }
      }
    }
    for (j = 0; j < numsp; j++) {
      x[minj - 1][j] = 0.0;
      x[j][minj - 1] = 0.0;
    }
    oc[mini - 1] += oc[minj - 1];
  }
  /* the last cycle */
  nude = 1;
  for (i = 1; i <= numsp; i++) {
    if (cluster[i - 1] != NULL) {
      el[nude - 1] = i;
      nude++;
    }
  }
  if (!njoin) {
    curtree.start = cluster[el[0] - 1];
    return;
  }
  bi = (x[el[0] - 1][el[1] - 1] + x[el[0] - 1][el[2] - 1] - x[el[1] - 1]
        [el[2] - 1]) * 0.5;
  bj = x[el[0] - 1][el[1] - 1] - bi;
  bk = x[el[0] - 1][el[2] - 1] - bi;
  bi -= av[el[0] - 1];
  bj -= av[el[1] - 1];
  bk -= av[el[2] - 1];
  if (progress) {
    printf("LAST CYCLE:\n");
    putchar(' ');
    nodelabel(av[el[0] - 1] > 0.0);
    printf("%3ld  (%10.5f) JOINS ", el[0], bi);
    nodelabel(av[el[1] - 1] > 0.0);
    printf("%3ld  (%10.5f) JOINS ", el[1], bj);
    nodelabel(av[el[2] - 1] > 0.0);
    printf("%3ld  (%10.5f)\n", el[2], bk);
  }
  hookup(curtree.nodep[nextnode - 1], cluster[el[0] - 1]);
  hookup(curtree.nodep[nextnode - 1]->next, cluster[el[1] - 1]);
  hookup(curtree.nodep[nextnode - 1]->next->next, cluster[el[2] - 1]);
  cluster[el[0] - 1]->v = bi;
  cluster[el[1] - 1]->v = bj;
  cluster[el[2] - 1]->v = bk;
  cluster[el[0] - 1]->back->v = bi;
  cluster[el[1] - 1]->back->v = bj;
  cluster[el[2] - 1]->back->v = bk;
  curtree.start = cluster[el[0] - 1]->back;
  free(av);
  free(oc);
}  /* jointree */


void maketree()
{
  /* construct the tree */
  long i, j, k;

  getdata();
  if (progress)
    putchar('\n');
  if (ith == 1)
    setuptree(&curtree);
  else
  for (i = 0; i < numsp; i++)
    memcpy(curtree.nodep[i]->nayme, names[i], sizeof(naym));
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
  for (i = 0; i < numsp; i++)
    cluster[i] = curtree.nodep[i];
  jointree();
  if (outgropt)
    curtree.start = curtree.nodep[outgrno - 1]->back;
  printree();
  if (treeprint)
    summarize();
  if (trout) {
    col = 0;
    treeout(curtree.start);
  }
  if (progress) {
    printf("\nOutput written on output file\n\n");
    if (trout)
      printf("Tree written on tree file\n\n");
  }
}  /* maketree */


main(argc, argv)
int argc;
Char *argv[];
{  /* main program */
char infilename[100],outfilename[100],trfilename[100];
#ifdef MAC
  macsetup("Neighbor","");
  argv[0] = "Neighbor";
#endif
  openfile(&infile,INFILE,"r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);
  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  mulsets = false;
  datasets = 1;
  doinit();
  if (trout)
    openfile(&treefile,TREEFILE,"w",argv[0],trfilename);
  x = (vector *)Malloc(numsp*sizeof(vector));
  for (j = 0; j < numsp; j++)
    x[j] = (vector)Malloc(numsp*sizeof(double));
  reps = (intvector *)Malloc(numsp*sizeof(intvector));
  for (j = 0; j < numsp; j++)
    reps[j] = (intvector)Malloc(numsp*sizeof(long));
  names = (naym *)Malloc(numsp*sizeof(naym));
  enterorder = (long *)Malloc(numsp*sizeof(long));
  cluster = (node **)Malloc(numsp*sizeof(node *));
  ith = 1;
  while (ith <= datasets) {
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n",ith);
      if (progress)
        printf("Data set # %ld:\n",ith);
    }
    getinput();
    maketree();
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    ith++;
  }
  FClose(infile);
  FClose(outfile);
  FClose(treefile);
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
mem = (MALLOCRETURN *)malloc(x);
if (!mem)
  memerror();
else
  return (MALLOCRETURN *)mem;
}

