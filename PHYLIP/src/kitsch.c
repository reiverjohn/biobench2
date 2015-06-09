#include "phylip.h"

/* version 3.56c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define namelength      10   /* number of characters in species name    */
#define epsilon         0.000001   /* a very small but not too small number */

#define ibmpc0          false
#define ansi0           true
#define vt520           false
#define down            2
#define over            60

typedef double *vector;       /* nodes will form binary tree           */
typedef Char naym[namelength];

typedef struct node {         /* describes a tip species or an ancestor */
  struct node *next, *back;
  long index;
  boolean tip;                /* present species are tips of tree       */
  vector d, w;                /* distances and weights                  */
  double t;                   /* time                                   */
  boolean sametime;           /* bookkeeps scrunched nodes              */
  double weight;              /* weight of node used by scrunch         */
  boolean processed;          /* used by evaluate                       */
  long xcoord, ycoord, ymin;  /* used by printree                       */
  long ymax;
} node;

typedef node **pointptr;
typedef long longer[6];

Static node *root, *best;
Static FILE *infile, *outfile, *treefile;
Static long numsp, numsp2, inseed, numtrees, col, datasets, ith,
            i, j, l, jumb, njumble;
/* numsp = number of species
   numtrees is used by usertree option part of maketree */
Static pointptr treenode, bestree;   /* pointers to all nodes in tree */
Static naym *nayms;
Static boolean jumble, usertree, lower, upper, negallowed, replicates, trout,
               printdata, progress, treeprint, mulsets, ibmpc, vt52,
               ansi, firstset;
Static double power;
Static longer seed;
Static long *enterorder;
Char ch;
/* Local variables for maketree, propagated globally for C version: */
  long examined;
  double like, bestyet;
  node *there;


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
{/* convert a character to upper case -- either ASCII or EBCDIC */
    *ch = islower(*ch) ? toupper(*ch) :(*ch);
}  /* uppercase */


void getnums()
{
  /* read species number */
  fscanf(infile, "%ld", &numsp);
  fprintf(outfile, "\n%4ld Populations\n", numsp);
  numsp2 = numsp * 2 - 1;
}  /* getnums */

void getoptions()
{
  /* interactively set options */
  long i, inseed0;
  Char ch;
  boolean done, done1;

  fprintf(outfile, "\nFitch-Margoliash method ");
  fprintf(outfile, "with contemporary tips, version %s\n\n",VERSION);
  jumble = false;
  njumble = 1;
  lower = false;
  negallowed = false;
  power = 2.0;
  replicates = false;
  upper = false;
  usertree = false;
  trout = true;
  printdata = false;
  progress = true;
  treeprint = true;
  for(;;) {
   printf( ansi ? "\033[2J\033[H" :
           vt52 ? "\033E\033H"    :
                  "\n");
    printf("\nFitch-Margoliash method ");
    printf("with contemporary tips, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
           usertree ? "No, use user trees in input file" : "Yes");
    printf("  P                                Power?%9.5f\n",power);
    printf("  -      Negative branch lengths allowed?  %s\n",
           (negallowed ? "Yes" : "No"));
    printf("  L         Lower-triangular data matrix?  %s\n",
           (lower ? "Yes" : "No"));
    printf("  R         Upper-triangular data matrix?  %s\n",
           (upper ? "Yes" : "No"));
    printf("  S                        Subreplicates?  %s\n",
           (replicates ? "Yes" : "No"));
    if (!usertree) {
      printf("  J     Randomize input order of species?");
      if (jumble)
            printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  0   Terminal type (IBM PC, VT52, ANSI)?  %s\n",
           (ibmpc ? "IBM PC" :
            ansi  ? "ANSI"   :
            vt52  ? "VT52"   :
                    "(none)"));

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
   if (strchr("JUP-LRSM12340",ch)){
     switch (ch) {

     case '-':
       negallowed = !negallowed;
       break;

     case 'J':
       jumble = !jumble;
       if (jumble) {
	 printf("Random number seed (must be odd)?\n");
	 scanf("%ld%*[^\n]", &inseed);
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
	 scanf("%ld%*[^\n]", &njumble);
	 getchar();
       }
       else njumble = 1;
       break;

     case 'L':
       lower = !lower;
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
  long i, j;
  node *p, *q;

  getnums();
  getoptions();
  treenode = (node **)Malloc(numsp2*sizeof(node *));
  for (i = 0; i < (numsp); i++) {
    treenode[i] = (node *)Malloc(sizeof(node));
    treenode[i]->d = (vector)Malloc(numsp2*(sizeof(double)));
    treenode[i]->w = (vector)Malloc(numsp2*(sizeof(double)));
  }
  for (i = numsp; i < (numsp2); i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc(sizeof(node));
      p->d = (vector)Malloc(numsp2*(sizeof(double)));
      p->w = (vector)Malloc(numsp2*(sizeof(double)));
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    treenode[i] = p;
  }
  if (!usertree && njumble > 1) {
    bestree = (node **)Malloc(numsp2*sizeof(node *));
    for (i = 0; i < (numsp); i++) {
      bestree[i] = (node *)Malloc(sizeof(node));
      bestree[i]->d = (vector)Malloc(numsp2*(sizeof(double)));
      bestree[i]->w = (vector)Malloc(numsp2*(sizeof(double)));
    }
    for (i = numsp; i < (numsp2); i++) {
      q = NULL;
      for (j = 1; j <= 3; j++) {
        p = (node *)Malloc(sizeof(node));
        p->d = (vector)Malloc(numsp2*(sizeof(double)));
        p->w = (vector)Malloc(numsp2*(sizeof(double)));
        p->next = q;
        q = p;
      }
      p->next->next->next = p;
      bestree[i] = p;
    }
  }

}  /* doinit */

void inputoptions()
{
  /* read options information */
  Char ch;
  long cursp;

  if (!firstset) {
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
  fprintf(outfile, "                  __ __             2\n");
  fprintf(outfile, "                  \\  \\   (Obs - Exp)\n");
  fprintf(outfile, "Sum of squares =  /_ /_  ------------\n");
  fprintf(outfile, "                               ");
  if (power == (long)power)
    fprintf(outfile, "%2ld\n", (long)power);
  else
    fprintf(outfile, "%4.1f\n", power);
  fprintf(outfile, "                   i  j      Obs\n\n");
  fprintf(outfile, "negative branch lengths");
  if (!negallowed)
    fprintf(outfile, " not");
  fprintf(outfile, " allowed\n\n");
}  /* inputoptions */

void getinput()
{
  /* reads the input data */
  inputoptions();
}  /* getinput */

void getdata()
{
  /* read in distance matrix */
  long i, j, k, columns, n;
  boolean skipit, skipother;
  double x;
  node *p;

  columns = replicates ? 4 : 6;
  if (printdata) {
    fprintf(outfile, "\nName                       Distances");
    if (replicates)
      fprintf(outfile, " (replicates)");
    fprintf(outfile, "\n----                       ---------");
    if (replicates)
      fprintf(outfile, "-------------");
    fprintf(outfile, "\n\n");
  }
  for (i = 1; i <= (numsp2); i++) {
    treenode[i - 1]->back = NULL;
    treenode[i - 1]->index = i;
    treenode[i - 1]->tip = (i <= numsp);
    treenode[i - 1]->t = 0.0;
    treenode[i - 1]->sametime = false;
    if (i > numsp) {
      p = treenode[i - 1]->next;
      while (p != treenode[i - 1]) {
        p->back = NULL;
        p->tip = false;
        p->index = i;
        p = p->next;
      }
    }
  }
  if (!usertree && njumble > 1)
    for (i = 1; i <= (numsp2); i++) {
      bestree[i - 1]->back = NULL;
      bestree[i - 1]->index = i;
      bestree[i - 1]->tip = (i <= numsp);
      bestree[i - 1]->t = 0.0;
      bestree[i - 1]->sametime = false;
      if (i > numsp) {
        p = bestree[i - 1]->next;
        while (p != bestree[i - 1]) {
          p->back = NULL;
          p->tip = false;
          p->index = i;
          p = p->next;
        }
      }
    }
  for (i = 0; i < (numsp); i++) {
    treenode[i]->d[i] = 0.0;
    treenode[i]->w[i] = 0.0;
    treenode[i]->weight = 0.0;
    fscanf(infile, "%*[^\n]");
    getc(infile);
    for (j = 0; j < namelength; j++) {
      if (eoln(infile))
	nayms[i][j] = ' ';
      else {
	nayms[i][j] = getc(infile);
	if (nayms[i][j] == '\n')
	  nayms[i][j] = ' ';
      }
    }
    for (j = 1; j <= (numsp); j++) {
      skipit = ((lower && j >= i + 1) || (upper && j <= i + 1));
      skipother = ((lower && i + 1 >= j) || (upper && i + 1 <= j));
      if (!skipit) {
	if (eoln(infile)) {
	  fscanf(infile, "%*[^\n]");
	  getc(infile);
	}
	fscanf(infile, "%lf", &x);
	treenode[i]->d[j - 1] = x;
	if (replicates) {
	  if (eoln(infile)) {
	    fscanf(infile, "%*[^\n]");
	    getc(infile);
	  }
	  fscanf(infile, "%ld", &n);
	} else
	  n = 1;
	if (n > 0 && x < 0) {
	  printf("NEGATIVE DISTANCE BETWEEN SPECIES%5ld AND %5ld\n",
		 i + 1, j);
	  exit(-1);
	}
	treenode[i]->w[j - 1] = n;
	if (skipother) {
	  treenode[j - 1]->d[i] = treenode[i]->d[j - 1];
	  treenode[j - 1]->w[i] = treenode[i]->w[j - 1];
	}
      }
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  if (printdata) {
    for (i = 0; i < (numsp); i++) {
      for (j = 0; j < namelength; j++)
        putc(nayms[i][j], outfile);
      putc(' ', outfile);
      for (j = 1; j <= (numsp); j++) {
        fprintf(outfile, "%10.5f", treenode[i]->d[j - 1]);
        if (replicates)
          fprintf(outfile, " (%3ld)", (long)treenode[i]->w[j - 1]);
        if (j % columns == 0 && j < numsp) {
          putc('\n', outfile);
          for (k = 1; k <= namelength + 1; k++)
            putc(' ', outfile);
        }
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  for (i = 0; i < (numsp); i++) {
    for (j = 0; j < (numsp); j++) {
      if (i + 1 != j + 1) {
        if (treenode[i]->d[j] < epsilon)
          treenode[i]->d[j] = epsilon;
        treenode[i]->w[j] /= exp(power * log(treenode[i]->d[j]));
      }
    }
  }
}  /* getdata */

void add(below, newtip, newfork)
node *below, *newtip, *newfork;
{
  /* inserts the nodes newfork and its left descendant, newtip,
     to the tree.  below becomes newfork's right descendant */
  if (below != treenode[below->index - 1])
    below = treenode[below->index - 1];
  if (below->back != NULL)
    below->back->back = newfork;
  newfork->back = below->back;
  below->back = newfork->next->next;
  newfork->next->next->back = below;
  newfork->next->back = newtip;
  newtip->back = newfork->next;
  if (root == below)
    root = newfork;
  root->back = NULL;
}  /* add */

void re_move(item, fork)
node **item, **fork;
{
  /* removes nodes item and its ancestor, fork, from the tree.
     the new descendant of fork's ancestor is made to be
     fork's second descendant (other than item).  Also
     returns pointers to the deleted nodes, item and fork */
  node *p, *q;

  if ((*item)->back == NULL) {
    *fork = NULL;
    return;
  }
  *fork = treenode[(*item)->back->index - 1];
  if (root == *fork) {
    if (*item == (*fork)->next->back)
      root = (*fork)->next->next->back;
    else
      root = (*fork)->next->back;
  }
  p = (*item)->back->next->back;
  q = (*item)->back->next->next->back;
  if (p != NULL)
    p->back = q;
  if (q != NULL)
    q->back = p;
  (*fork)->back = NULL;
  p = (*fork)->next;
  while (p != *fork) {
    p->back = NULL;
    p = p->next;
  }
  (*item)->back = NULL;
}  /* remove */

void scrunchtraverse(u,closest,tmax)
node *u,**closest;
double *tmax;
{
  /* traverse to find closest node to the current one */
  if (!u->sametime) {
    if (u->t > *tmax) {
      *closest = u;
      *tmax = u->t;
    }
    return;
  }
  u->t = treenode[u->back->index - 1]->t;
  if (!u->tip) {
    scrunchtraverse(u->next->back, closest,tmax);
    scrunchtraverse(u->next->next->back, closest,tmax);
  }
}  /* scrunchtraverse */

void combine(a, b)
node *a, *b;
{
  /* put node b into the set having the same time as a */
  if (a->weight + b->weight <= 0.0)
    a->t = 0.0;
  else
    a->t = (a->t * a->weight + b->t * b->weight) / (a->weight + b->weight);
  a->weight += b->weight;
  b->sametime = true;
}  /* combine */

void scrunch(s)
node *s;
{
  /* see if nodes can be combined to prevent negative lengths */
/* Local variables for scrunch: */
  double tmax;
  node *closest;
  boolean found;

  closest = NULL;
  tmax = -1.0;
  do {
    if (!s->tip) {
      scrunchtraverse(s->next->back, &closest,&tmax);
      scrunchtraverse(s->next->next->back, &closest,&tmax);
    }
    found = (tmax > s->t);
    if (found)
      combine(s, closest);
    tmax = -1.0;
  } while (found);
}  /* scrunch */

void secondtraverse(a,q,u,v,i,j,k,sum)
node *a,*q,*u,*v;
long i,j,k;
double *sum;
{
  /* recalculate distances, add to sum */
  long l;
  double wil, wjl, wkl, wli, wlj, wlk, TEMP;

  if (!(a->processed || a->tip)) {
    secondtraverse(a->next->back, q,u,v,i,j,k,sum);
    secondtraverse(a->next->next->back, q,u,v,i,j,k,sum);
    return;
  }
  if (!(a != q && a->processed))
    return;
  l = a->index;
  wil = u->w[l - 1];
  wjl = v->w[l - 1];
  wkl = wil + wjl;
  wli = a->w[i - 1];
  wlj = a->w[j - 1];
  wlk = wli + wlj;
  q->w[l - 1] = wkl;
  a->w[k - 1] = wlk;
  if (wkl <= 0.0)
    q->d[l - 1] = 0.0;
  else
    q->d[l - 1] = (wil * u->d[l - 1] + wjl * v->d[l - 1]) / wkl;
  if (wlk <= 0.0)
    a->d[k - 1] = 0.0;
  else
    a->d[k - 1] = (wli * a->d[i - 1] + wlj * a->d[j - 1]) / wlk;
  if (wkl > 0.0) {
    TEMP = u->d[l - 1] - v->d[l - 1];
    (*sum) += wil * wjl / wkl * (TEMP * TEMP);
  }
  if (wlk > 0.0) {
    TEMP = a->d[i - 1] - a->d[j - 1];
    (*sum) += wli * wlj / wlk * (TEMP * TEMP);
  }
}  /* secondtraverse */

void firstraverse(q_, r,sum)
node *q_,*r;
double *sum;
{  /* firsttraverse                              */
   /* go through tree calculating branch lengths */
  /* Local variables for firstraverse: */
  node *q;
  long i, j, k;
  node *u, *v;

  q = q_;
  if (q == NULL)
    return;
  q->sametime = false;
  if (!q->tip) {
    firstraverse(q->next->back, r,sum);
    firstraverse(q->next->next->back, r,sum);
  }
  q->processed = true;
  if (q->tip)
    return;
  u = q->next->back;
  v = q->next->next->back;
  i = u->index;
  j = v->index;
  k = q->index;
  if (u->w[j - 1] + v->w[i - 1] <= 0.0)
    q->t = 0.0;
  else
    q->t = (u->w[j - 1] * u->d[j - 1] +
              v->w[i - 1] * v->d[i - 1]) /
             (2.0 * (u->w[j - 1] + v->w[i - 1]));
  q->weight = u->weight + v->weight + u->w[j - 1] + v->w[i - 1];
  if (!negallowed)
    scrunch(q);
  secondtraverse(r,q,u,v,i,j,k,sum);
}  /* firstraverse */

void sumtraverse(q, sum)
node *q;
double *sum;
{
  /* traverse to finish computation of sum of squares */
  long i, j;
  node *u, *v;
  double TEMP, TEMP1;

  if (q->tip)
    return;
  sumtraverse(q->next->back, sum);
  sumtraverse(q->next->next->back, sum);
  u = q->next->back;
  v = q->next->next->back;
  i = u->index;
  j = v->index;
  TEMP = u->d[j - 1] - 2.0 * q->t;
  TEMP1 = v->d[i - 1] - 2.0 * q->t;
  (*sum) += u->w[j - 1] * (TEMP * TEMP) + v->w[i - 1] * (TEMP1 * TEMP1);
}  /* sumtraverse */

void evaluate(r)
node *r;
{
  /* fill in times and evaluate sum of squares for tree */
  /* Local variables for evaluate: */
  double sum;
  long i;
  sum = 0.0;
  for (i = 0; i < (numsp2); i++)
    treenode[i]->processed = treenode[i]->tip;
  firstraverse(r, r,&sum);
  sumtraverse(r, &sum);
  examined++;
  if (replicates && (lower || upper))
    sum /= 2;
  like = -sum;
}  /* evaluate */


void tryadd(p, item,nufork)
node *p,**item,**nufork;
{
  /* temporarily adds one fork and one tip to the tree.
     if the location where they are added yields greater
     "likelihood" than other locations tested up to that
     time, then keeps that location as there */
  add(p, *item, *nufork);
  evaluate(root);
  if (like > bestyet) {
    bestyet = like;
    there = p;
  }
  re_move(item, nufork);
}  /* tryadd */

void addpreorder(p, item, nufork)
node *p, *item, *nufork;
{
  /* traverses a binary tree, calling PROCEDURE tryadd
     at a node before calling tryadd at its descendants */
/* Local variables for addpreorder: */
  if (p == NULL)
    return;
  tryadd(p, &item,&nufork);
  if (!p->tip) {
    addpreorder(p->next->back, item, nufork);
    addpreorder(p->next->next->back, item, nufork);
  }
}  /* addpreorder */

void tryrearr(p,r,success)
node *p,**r;
boolean *success;
{
  /* evaluates one rearrangement of the tree.
     if the new tree has greater "likelihood" than the old
     one sets success := TRUE and keeps the new tree.
     otherwise, restores the old tree */
  node *frombelow, *whereto, *forknode;
  double oldlike;

  if (p->back == NULL)
    return;
  forknode = treenode[p->back->index - 1];
  if (forknode->back == NULL)
    return;
  oldlike = like;
  if (p->back->next->next == forknode)
    frombelow = forknode->next->next->back;
  else
    frombelow = forknode->next->back;
  whereto = forknode->back;
  re_move(&p, &forknode);
  add(whereto, p, forknode);
  if ((*r)->back != NULL)
    *r = treenode[(*r)->back->index - 1];
  evaluate(*r);
  if (like > oldlike) {
    bestyet = like;
    *success = true;
    return;
  }
  re_move(&p, &forknode);
  add(frombelow, p, forknode);
  if ((*r)->back != NULL)
    *r = treenode[(*r)->back->index - 1];
  like = oldlike;
}  /* tryrearr */

void repreorder(p,r,success)
node *p,**r;
boolean *success;
{
  /* traverses a binary tree, calling PROCEDURE tryrearr
     at a node before calling tryrearr at its descendants */
  if (p == NULL)
    return;
  tryrearr(p,r,success);
  if (!p->tip) {
    repreorder(p->next->back,r,success);
    repreorder(p->next->next->back,r,success);
  }
}  /* repreorder */

void rearrange(r_)
node **r_;
{
  /* traverses the tree (preorder), finding any local
     rearrangement which decreases the number of steps.
     if traversal succeeds in increasing the tree's
     "likelihood", PROCEDURE rearrange runs traversal again */
/* Local variables for rearrange: */
  node **r;
  boolean success;
  r = r_;
  success = true;
  while (success) {
    success = false;
    repreorder(*r,r,&success);
  }
}  /* rearrange */

void findch(c)
Char c;
{
  /* scan forward until find character c */
  boolean done;

  done = false;
  while (!(done)) {
    if (c == ',') {
      if (ch == '(' || ch == ')' || ch == ';') {
        printf("\nERROR IN USER TREE: ");
	printf("UNMATCHED PARENTHESIS OR MISSING COMMA\n");
	exit(-1);
      } else if (ch == ',')
        done = true;
    } else if (c == ')') {
      if (ch == '(' || ch == ',' || ch == ';') {
        printf("\nERROR IN USER TREE:");
	printf(" UNMATCHED PARENTHESIS OR NOT BIFURCATED NODE\n");
	exit(-1);
      } else {
        if (ch == ')')
          done = true;
      }
    } else if (c == ';') {
      if (ch != ';') {
        printf("\nERROR IN USER TREE:");
	printf(" UNMATCHED PARENTHESIS OR MISSING SEMICOLON\n");
	exit(-1);
      } else
        done = true;
    }
    if (done && ch != ')')   /* was : if (!(done && ch == ')') && (done))  */
      continue;
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
  }
}  /* findch */

void addelement(p, nextnode,lparens,names)
node **p;
long *nextnode,*lparens;
boolean *names;
{
  /* recursive procedure adds nodes to user-defined tree */
  node *q;
  long i, n;
  boolean found;
  Char str[namelength];

  do {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
  } while (ch == ' ');
  if (ch == '(' ) {
    if (*lparens >= numsp - 1) {
      printf("\nERROR IN USER TREE: TOO MANY LEFT PARENTHESES\n");
      exit(-1);
    }
    (*nextnode)++;
    (*lparens)++;
    q = treenode[*nextnode - 1];
    addelement(&q->next->back, nextnode,lparens,names);
    q->next->back->back = q->next;
    findch(',');
    addelement(&q->next->next->back, nextnode,lparens,names);
    q->next->next->back->back = q->next->next;
    findch(')');
    *p = q;
    return;
  }

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
    if (ch == '\n')
      ch = ' ';
    n++;
  } while (ch != ',' && ch != ')' && ch != ':' &&
           n <= namelength);
  n = 1;
  do {
    found = true;
    for (i = 0; i < namelength; i++)
      found = (found && str[i] == nayms[n - 1][i]);
    if (found) {
      if (names[n - 1] == false) {
        *p = treenode[n - 1];
        names[n - 1] = true;
      } else {
        printf("\nERROR IN USER TREE: DUPLICATE NAME FOUND -- ");
        for (i = 0; i < namelength; i++)
          putchar(nayms[n - 1][i]);
        putchar('\n');
	exit(-1);
      }
    } else
      n++;
  } while (!(n > numsp || found ));
  if (n <= numsp)
    return;
  printf("CANNOT FIND SPECIES: ");
  for (i = 0; i < namelength; i++)
    putchar(str[i]);
  putchar('\n');
}  /* addelement */

void treeread()
{
  /* read in user-defined tree and set it up */
  long nextnode, lparens;
  boolean *names;
  long i;

  root = treenode[numsp];
  nextnode = numsp;
  root->back = NULL;
  names = (boolean *)Malloc(numsp*sizeof(boolean));
  for (i = 0; i < (numsp); i++)
    names[i] = false;
  lparens = 0;
  addelement(&root, &nextnode,&lparens,names);
  findch(';');
  fscanf(infile, "%*[^\n]");
  getc(infile);
}  /* treeread */

void coordinates(p, tipy)
node *p;
long *tipy;
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
    return;
  }
  q = p->next;
  do {
    coordinates(q->back, tipy);
    q = q->next;
  } while (p != q);
  first = p->next->back;
  q = p->next;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (long)(over * p->t + 0.5);
  p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */

void drawline(i, scale)
long i;
double scale;
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q, *r, *first, *last;
  long n, j;
  boolean extra, done;

  p = root;
  q = root;
  extra = false;
  if (i == p->ycoord && p == root) {
    if (p->index - numsp >= 10)
      fprintf(outfile, "-%2ld", p->index - numsp);
    else
      fprintf(outfile, "--%ld", p->index - numsp);
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
      } while (!(done || r == p));
      first = p->next->back;
      r = p->next;
      while (r->next != p)
        r = r->next;
      last = r->back;
    }
    done = (p == q);
    n = (long)(scale * (p->xcoord - q->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if (q->ycoord == i && !done) {
      if (p->ycoord != q->ycoord)
        putc('+', outfile);
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        if (q->index - numsp >= 10)
          fprintf(outfile, "%2ld", q->index - numsp);
        else
          fprintf(outfile, "-%ld", q->index - numsp);
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
    }
    if (p != q)
      p = q;
  } while (!done);
  if (p->ycoord == i && p->tip) {
    for (j = 0; j < namelength; j++)
      putc(nayms[p->index - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */

void printree()
{
 /* prints out diagram of the tree */
 /* Local variables for printree: */
  long tipy;
  double scale;
  long i;
  node *p;

  putc('\n', outfile);
  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  coordinates(root, &tipy);
  p = root;
  while (!p->tip)
    p = p->next->back;
  scale = 1.0 / (long)(root->t - p->t + 1.000);
  putc('\n', outfile);
  for (i = 1; i <= (tipy - down); i++)
    drawline(i, scale);
  putc('\n', outfile);
}  /* printree */

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
      if (nayms[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = nayms[p->index - 1][i];
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
    putc(')', treefile);
    col++;
  }
  if (p != root)
    x = treenode[p->back->index - 1]->t - p->t;
  if (x > 0.0)
    w = (long)(0.43429448222 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.43429448222 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if (p == root)
    fprintf(treefile, ";\n");
  else {
    fprintf(treefile, ":%*.5f", (int)(w + 7), x);
    col += w + 8;
  }
}  /* treeout */

void dtraverse(q)
node *q;
{
  /* print table of lengths etc. */
  long i;

  if (!q->tip)
    dtraverse(q->next->back);
  if (q->back != NULL) {
    fprintf(outfile, "%4ld  ", q->back->index - numsp);
    if (q->index <= numsp) {
      for (i = 0; i < namelength; i++)
        putc(nayms[q->index - 1][i], outfile);
    } else
      fprintf(outfile, "%4ld      ", q->index - numsp);
    fprintf(outfile, "%13.5f", treenode[q->back->index - 1]->t - q->t);
    fprintf(outfile, "%15.5f\n", root->t - q->t);
  }
  if (!q->tip)
    dtraverse(q->next->next->back);
}  /* dtraverse */

void describe()
{
  /* prints table of lengths, times, sum of squares, etc. */
  long i, j;
  double totalnum;
  double TEMP;

  fprintf(outfile, "\nSum of squares = %10.3f\n\n", -like);
  if (fabs(power - 2) < 0.01) {
    totalnum = 0.0;
    for (i = 0; i < (numsp); i++) {
      for (j = 0; j < (numsp); j++) {
        if (i + 1 != j + 1 && treenode[i]->d[j] > 0.0) {
          TEMP = treenode[i]->d[j];
          totalnum += treenode[i]->w[j] * (TEMP * TEMP);
        }
      }
    }
    totalnum -= 2;
    if (replicates && (lower || upper))
      totalnum /= 2;
    fprintf(outfile, "Average percent standard deviation =");
    fprintf(outfile, "%10.5f\n\n", 100 * sqrt(-(like / totalnum)));
  }
  if (!usertree)
    fprintf(outfile, "examined %4ld trees\n\n", examined);
  fprintf(outfile, "From    To           Length          Time\n");
  fprintf(outfile, "----    --           ------          ----\n\n");
  dtraverse(root);
  putc('\n', outfile);
  if (trout) {
    col = 0;
    treeout(root);
  }
}  /* describe */


void copynode(c, d)
node *c, *d;
{
  /* make a copy of a node */

  memcpy(d->d, c->d, numsp2*sizeof(double));
  memcpy(d->w, c->w, numsp2*sizeof(double));
  d->t = c->t;
  d->sametime = c->sametime;
  d->weight = c->weight;
  d->processed = c->processed;
  d->xcoord = c->xcoord;
  d->ycoord = c->ycoord;
  d->ymin = c->ymin;
  d->ymax = c->ymax;
}  /* copynode */

void copy_(a, b)
pointptr a, b;
{
  /* make a copy of a tree */
  short i, j=0;
  node *p, *q;

  for (i = 0; i < numsp; i++) {
    copynode(a[i], b[i]);
    if (a[i]->back != NULL) {
      if (a[i]->back == a[a[i]->back->index - 1])
        b[i]->back = b[a[i]->back->index - 1];
      else if (a[i]->back == a[a[i]->back->index - 1]->next)
      b[i]->back = b[a[i]->back->index - 1]->next;
    else
      b[i]->back = b[a[i]->back->index - 1]->next->next;
    }
    else b[i]->back = NULL;
  }
  for (i = numsp; i < numsp2; i++) {
    p = a[i];
    q = b[i];
    for (j = 1; j <= 3; j++) {
      copynode(p, q);
      if (p->back) {
        if (p->back == a[p->back->index - 1])
          q->back = b[p->back->index - 1];
        else if (p->back == a[p->back->index - 1]->next)
          q->back = b[p->back->index - 1]->next;
        else
          q->back = b[p->back->index - 1]->next->next;
      }
      else
        q->back = NULL;
      p = p->next;
      q = q->next;
    }
  }
}  /* copy */


void maketree()
{
  /* constructs a binary tree from the pointers in treenode.
     adds each node at location which yields highest "likelihood"
     then rearranges the tree for greatest "likelihood" */
  long i, j, k;
  double bestlike, bstlike2, gotlike;
  boolean lastrearr;
  node *item, *nufork;

  if (!usertree) {
    if (jumb == 1) {
      getdata();
      examined = 0;
    }
    for (i = 1; i <= (numsp); i++)
      enterorder[i - 1] = i;
    if (jumble) {
      for (i = 0; i < (numsp); i++) {
        j = (long)(randum(seed) * numsp) + 1;
        k = enterorder[j - 1];
        enterorder[j - 1] = enterorder[i];
        enterorder[i] = k;
      }
    }
    root = treenode[enterorder[0] - 1];
    add(treenode[enterorder[0] - 1], treenode[enterorder[1] - 1],
        treenode[numsp]);
    if (progress) {
      printf("\nAdding species:\n");
      printf("   ");
      for (i = 0; i < namelength; i++)
        putchar(nayms[enterorder[0] - 1][i]);
      printf("\n   ");
      for (i = 0; i < namelength; i++)
        putchar(nayms[enterorder[1] - 1][i]);
      putchar('\n');
    }
    lastrearr = false;
    for (i = 3; i <= (numsp); i++) {
      bestyet = -10000000.0;
      item = treenode[enterorder[i - 1] - 1];
      nufork = treenode[numsp + i - 2];
      addpreorder(root, item, nufork);
      add(there, item, nufork);
      like = bestyet;
      rearrange(&root);
      evaluate(root);
      examined--;
      if (progress) {
        printf("   ");
        for (j = 0; j < namelength; j++)
          putchar(nayms[enterorder[i - 1] - 1][j]);
        putchar('\n');
      }
      lastrearr = (i == numsp);
      if (lastrearr) {
        if (progress) {
          printf("\nDoing global rearrangements\n");
          printf("  !");
          for (j = 1; j <= (numsp2); j++)
            putchar('-');
          printf("!\n");
        }
        bestlike = bestyet;
        do {
          gotlike = bestlike;
          if (progress)
            printf("   ");
          for (j = 0; j < (numsp2); j++) {
            there = root;
            bestyet = -32000.0;
            item = treenode[j];
            if (item != root) {
              re_move(&item, &nufork);
              there = root;
              addpreorder(root, item, nufork);
              add(there, item, nufork);
            }
            if (progress){
              putchar('.');
              fflush(stdout);}
          }
          if (progress)
            putchar('\n');
        } while (bestlike > gotlike);
        if (njumble > 1) {
          if (jumb == 1 || (jumb > 1 && bestlike > bstlike2)) {
            copy_(treenode, bestree);
            best = bestree[root->index -1];
            bstlike2 = bestlike;
          }
        }
      }
      if (i == numsp && njumble == jumb) {
        if (njumble > 1) {
          copy_(bestree, treenode);
          root = treenode[best->index - 1];
        }
        evaluate(root);
        printree();
        describe();
      }
    }
  } else {
    getdata();
    fscanf(infile, "%ld%*[^\n]", &numtrees);
    getc(infile);
    if (treeprint)
      fprintf(outfile, "\n\nUser-defined trees:\n\n");
    i = 1;
    while (i <= numtrees ) {
      treeread();
      evaluate(root);
      printree();
      describe();
      i++;
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
{  /* Fitch-Margoliash criterion with contemporary tips */
char infilename[100],outfilename[100],trfilename[100];
#ifdef MAC
  macsetup("Kitsch","");
  argv[0] = "Kitsch"; 
#endif
  /* reads in numsp, options, and the data, then calls maketree to
     construct the tree */
  openfile(&infile,INFILE,"r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);

  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  mulsets = false;
  firstset = true;
  datasets = 1;
  doinit();
  openfile(&treefile,TREEFILE,"w",argv[0],trfilename);
  nayms = (naym *)Malloc(numsp*sizeof(naym));
  enterorder = (long *)Malloc(numsp*sizeof(long));
  for (ith = 1; ith <= datasets; ith++) {
    if (datasets > 1) {
      fprintf(outfile, "\nData set # %ld:\n",ith);
      if (progress)
        printf("\nData set # %ld:\n",ith);
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
  FClose(infile);
  FClose(outfile);
  FClose(treefile);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(trfilename);
#endif
  exit(0);
}  /* Fitch-Margoliash criterion with contemporary tips */


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

