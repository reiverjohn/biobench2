#include "phylip.h"

/* version 3.52c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define nmlngth         10   /* number of characters in species name    */
#define maxtrees        50   /* maximum number of tied trees stored     */
#define maxuser         8    /* maximum number of user-defined trees    */
#define down            2

#define ibmpc0          false
#define ansi0           true
#define vt520           false


typedef long *steptr;
typedef enum {
  ala, arg, asn, asp, cys, gln, glu, gly, his, ileu, leu, lys, met, phe, pro,
  ser1, ser2, thr, trp, tyr, val, del, stop, asx, glx, ser, unk, quest
} aas;
typedef long sitearray[3];
typedef sitearray *seqptr;

/* nodes will form a binary tree */

typedef struct node {        /* describes a tip species or an ancestor */
  struct node *next, *back;  /* pointers to nodes                      */
  long index;                /* number of the node                     */
  boolean tip, bottom;       /* present species are tips of tree       */
  aas *seq;                  /* the sequence                           */
  seqptr siteset;            /* temporary storage for aa's             */
  steptr numsteps;           /* bookkeeps steps                        */
  long xcoord, ycoord, ymin; /* used by printree                       */

  long ymax;
} node;

typedef node **pointptr;
typedef long longer[6];

typedef struct gseq {
  seqptr seq;
  struct gseq *next;
} gseq;


Static node *root;
Static FILE *infile, *outfile, *treefile;
Static long spp, nonodes, chars, inseed, outgrno, col, datasets, ith,
       i, j, l, jumb, njumble;
/* spp = number of species
   nonodes = number of nodes in tree
   chars = number of sites in actual sequences
   outgrno indicates outgroup */
Static boolean jumble, usertree, weights, thresh, trout, outgropt,
               printdata, progress, treeprint, stepbox, ancseq, mulsets,
               interleaved, ibmpc, vt52, ansi, firstset;
Static long fullset, fulldel;
Static steptr weight;
Static pointptr treenode;   /* pointers to all nodes in tree */
Static Char **nayme;   /* names of species */
Static double threshold;
Static steptr threshwt;
Static longer seed;
Static long *enterorder;
Static sitearray translate[(long)quest - (long)ala + 1];
Static long **fsteps;
Static long **bestrees;
Static gseq *garbage;
Static node *temp, *temp1;
Char ch;
aas tmpa;

/* Local variables for maketree, propagated globally for c version: */
long nextree, which, minwhich;
double like, bestyet, bestlike, minsteps, bstlike2;
boolean lastrearr, recompute;
node *there;
double nsteps[maxuser];
long *place;
boolean *names;



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


void gnu(p)
gseq **p;
{
  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */
  if (garbage != NULL) {
    *p = garbage;
    garbage = garbage->next;
  } else {
    *p = (gseq *)Malloc(sizeof(gseq));
    (*p)->seq = (seqptr)Malloc(chars*sizeof(sitearray));
  }
  (*p)->next = NULL;
}  /* gnu */


void chuck(p)
gseq *p;
{
  /* collect garbage on p -- put it on front of garbage list */
  p->next = garbage;
  garbage = p;
}  /* chuck */


void setup()
{
  /* set up set table to get aasets from aas */
  aas a, b;
  long s;

  for (a = ala; (long)a <= (long)stop; a = (aas)((long)a + 1))
    translate[(long)a - (long)ala][0] = 1L << ((long)a);
  translate[0]
    [1] = (1L << ((long)ala)) | (1L << ((long)asp)) | (1L << ((long)glu)) |
          (1L << ((long)gly)) | (1L << ((long)ser1)) | (1L << ((long)pro)) |
          (1L << ((long)thr)) | (1L << ((long)val));
  translate[(long)arg - (long)ala]
    [1] = (1L << ((long)arg)) | (1L << ((long)cys)) | (1L << ((long)gln)) |
          (1L << ((long)gly)) | (1L << ((long)his)) | (1L << ((long)ileu)) |
          (1L << ((long)leu)) | (1L << ((long)lys)) | (1L << ((long)met)) |
          (1L << ((long)pro)) | (1L << ((long)ser2)) | (1L << ((long)thr)) |
          (1L << ((long)trp)) | (1L << ((long)stop));
  translate[(long)asn - (long)ala]
    [1] = (1L << ((long)asn)) | (1L << ((long)asp)) | (1L << ((long)his)) |
          (1L << ((long)ileu)) | (1L << ((long)lys)) | (1L << ((long)ser2)) |
          (1L << ((long)thr)) | (1L << ((long)tyr));
  translate[(long)asp - (long)ala]
    [1] = (1L << ((long)ala)) | (1L << ((long)asp)) | (1L << ((long)asn)) |
          (1L << ((long)glu)) | (1L << ((long)gly)) | (1L << ((long)his)) |
          (1L << ((long)tyr)) | (1L << ((long)val));
  translate[(long)cys - (long)ala]
    [1] = (1L << ((long)arg)) | (1L << ((long)cys)) | (1L << ((long)gly)) |
          (1L << ((long)phe)) | (1L << ((long)ser1)) | (1L << ((long)ser2)) |
          (1L << ((long)trp)) | (1L << ((long)tyr)) | (1L << ((long)stop));
  translate[(long)gln - (long)ala]
    [1] = (1L << ((long)arg)) | (1L << ((long)gln)) | (1L << ((long)glu)) |
          (1L << ((long)his)) | (1L << ((long)leu)) | (1L << ((long)lys)) |
          (1L << ((long)pro)) | (1L << ((long)stop));
  translate[(long)glu - (long)ala]
    [1] = (1L << ((long)ala)) | (1L << ((long)asp)) | (1L << ((long)gln)) |
          (1L << ((long)glu)) | (1L << ((long)gly)) | (1L << ((long)lys)) |
          (1L << ((long)val)) | (1L << ((long)stop));
  translate[(long)gly - (long)ala]
    [1] = (1L << ((long)ala)) | (1L << ((long)arg)) | (1L << ((long)asp)) |
          (1L << ((long)cys)) | (1L << ((long)glu)) | (1L << ((long)gly)) |
          (1L << ((long)ser2)) | (1L << ((long)trp)) | (1L << ((long)val)) |
          (1L << ((long)stop));
  translate[(long)his - (long)ala]
    [1] = (1L << ((long)arg)) | (1L << ((long)asn)) | (1L << ((long)asp)) |
          (1L << ((long)gln)) | (1L << ((long)his)) | (1L << ((long)leu)) |
          (1L << ((long)pro)) | (1L << ((long)tyr));
  translate[(long)ileu - (long)ala]
    [1] = (1L << ((long)arg)) | (1L << ((long)asn)) | (1L << ((long)ileu)) |
          (1L << ((long)leu)) | (1L << ((long)lys)) | (1L << ((long)met)) |
          (1L << ((long)phe)) | (1L << ((long)ser2)) | (1L << ((long)thr)) |
          (1L << ((long)val));
  translate[(long)leu - (long)ala]
    [1] = (1L << ((long)arg)) | (1L << ((long)gln)) | (1L << ((long)his)) |
          (1L << ((long)ileu)) | (1L << ((long)leu)) | (1L << ((long)met)) |
          (1L << ((long)phe)) | (1L << ((long)pro)) | (1L << ((long)ser1)) |
          (1L << ((long)trp)) | (1L << ((long)val)) | (1L << ((long)stop));
  translate[(long)lys - (long)ala]
    [1] = (1L << ((long)arg)) | (1L << ((long)asn)) | (1L << ((long)gln)) |
          (1L << ((long)glu)) | (1L << ((long)ileu)) | (1L << ((long)lys)) |
          (1L << ((long)met)) | (1L << ((long)thr)) | (1L << ((long)stop));
  translate[(long)met - (long)ala]
    [1] = (1L << ((long)arg)) | (1L << ((long)ileu)) | (1L << ((long)leu)) |
          (1L << ((long)lys)) | (1L << ((long)met)) | (1L << ((long)val)) |
          (1L << ((long)thr));
  translate[(long)phe - (long)ala]
    [1] = (1L << ((long)cys)) | (1L << ((long)ileu)) | (1L << ((long)leu)) |
          (1L << ((long)phe)) | (1L << ((long)ser1)) | (1L << ((long)tyr)) |
          (1L << ((long)val));
  translate[(long)pro - (long)ala]
    [1] = (1L << ((long)ala)) | (1L << ((long)arg)) | (1L << ((long)gln)) |
          (1L << ((long)his)) | (1L << ((long)leu)) | (1L << ((long)pro)) |
          (1L << ((long)ser1)) | (1L << ((long)thr));
  translate[(long)ser1 - (long)ala]
    [1] = (1L << ((long)ala)) | (1L << ((long)cys)) | (1L << ((long)leu)) |
          (1L << ((long)phe)) | (1L << ((long)pro)) | (1L << ((long)ser1)) |
          (1L << ((long)thr)) | (1L << ((long)trp)) | (1L << ((long)tyr)) |
          (1L << ((long)stop));
  translate[(long)ser2 - (long)ala]
    [1] = (1L << ((long)arg)) | (1L << ((long)asn)) | (1L << ((long)cys)) |
          (1L << ((long)gly)) | (1L << ((long)ileu)) | (1L << ((long)ser2)) |
          (1L << ((long)thr));
  translate[(long)thr - (long)ala]
    [1] = (1L << ((long)ala)) | (1L << ((long)arg)) | (1L << ((long)asn)) |
          (1L << ((long)ileu)) | (1L << ((long)lys)) | (1L << ((long)met)) |
          (1L << ((long)pro)) | (1L << ((long)ser1)) | (1L << ((long)ser2)) |
          (1L << ((long)thr));
  translate[(long)trp - (long)ala]
    [1] = (1L << ((long)arg)) | (1L << ((long)cys)) | (1L << ((long)gly)) |
          (1L << ((long)leu)) | (1L << ((long)ser1)) | (1L << ((long)stop)) |
          (1L << ((long)trp));
  translate[(long)tyr - (long)ala]
    [1] = (1L << ((long)asn)) | (1L << ((long)asp)) | (1L << ((long)cys)) |
          (1L << ((long)his)) | (1L << ((long)phe)) | (1L << ((long)ser1)) |
          (1L << ((long)stop)) | (1L << ((long)tyr));
  translate[(long)val - (long)ala]
    [1] = (1L << ((long)ala)) | (1L << ((long)asp)) | (1L << ((long)glu)) |
          (1L << ((long)gly)) | (1L << ((long)ileu)) | (1L << ((long)leu)) |
          (1L << ((long)met)) | (1L << ((long)phe)) | (1L << ((long)val));
  translate[(long)stop - (long)ala]
    [1] = (1L << ((long)arg)) | (1L << ((long)cys)) | (1L << ((long)gln)) |
          (1L << ((long)glu)) | (1L << ((long)gly)) | (1L << ((long)leu)) |
          (1L << ((long)lys)) | (1L << ((long)ser1)) | (1L << ((long)trp)) |
          (1L << ((long)tyr)) | (1L << ((long)stop));
  translate[(long)del - (long)ala][1] = 1L << ((long)del);
  fulldel = (1L << ((long)stop + 1)) - (1L << ((long)ala));
  fullset = fulldel & (~(1L << ((long)del)));
  translate[(long)asx - (long)ala]
    [0] = (1L << ((long)asn)) | (1L << ((long)asp));
  translate[(long)glx - (long)ala]
    [0] = (1L << ((long)gln)) | (1L << ((long)glu));
  translate[(long)ser - (long)ala]
    [0] = (1L << ((long)ser1)) | (1L << ((long)ser2));
  translate[(long)unk - (long)ala][0] = fullset;
  translate[(long)quest - (long)ala][0] = fulldel;
  translate[(long)asx - (long)ala]
    [1] = translate[(long)asn - (long)ala]
          [1] | translate[(long)asp - (long)ala][1];
  translate[(long)glx - (long)ala]
    [1] = translate[(long)gln - (long)ala]
          [1] | translate[(long)glu - (long)ala][1];
  translate[(long)ser - (long)ala]
    [1] = translate[(long)ser1 - (long)ala]
          [1] | translate[(long)ser2 - (long)ala][1];
  translate[(long)unk - (long)ala][1] = fullset;
  translate[(long)quest - (long)ala][1] = fulldel;
  for (a = ala; (long)a <= (long)quest; a = (aas)((long)a + 1)) {
    s = 0;
    for (b = ala; (long)b <= (long)stop; b = (aas)((long)b + 1)) {
      if (((1L << ((long)b)) & translate[(long)a - (long)ala][1]) != 0)
        s |= translate[(long)b - (long)ala][1];
    }
    translate[(long)a - (long)ala][2] = s;
  }
}  /* setup */


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
  /* convert ch to upper case -- either ASCII or EBCDIC */
    *ch = (islower (*ch) ? toupper(*ch) : (*ch));
}  /* uppercase */


void inputnumbers()
{
  /* input the numbers of species and of characters */
  fscanf(infile, "%ld%ld", &spp, &chars);
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n", spp, chars);
  if (printdata)
    putc('\n', outfile);
  nonodes = spp * 2 - 1;
}  /* inputnumbers */

void getoptions()
{
  /* interactively set options */
  long i, inseed0;
  Char ch;
  boolean done, done1;

  fprintf(outfile, "\nProtein parsimony algorithm, version %s\n\n",VERSION);
  putchar('\n');
  jumble = false;
  njumble = 1;
  outgrno = 1;
  outgropt = false;
  thresh = false;
  trout = true;
  usertree = false;
  weights = false;
  printdata = false;
  progress = true;
  treeprint = true;
  stepbox = false;
  ancseq = false;
  interleaved = true;
  for (;;) {
    printf(ansi ? "\033[2J\033[H" : vt52 ? "\033E\033H" : "\n");
    printf("\nProtein parsimony algorithm, version %s\n\n",VERSION);
    printf("Setting for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
           (usertree ? "No, use user trees in input file" : "Yes"));
    if (!usertree) {
      printf("  J   Randomize input order of sequences?");
      if (jumble)
        printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  O                        Outgroup root?");
    if (outgropt)
      printf("  Yes, at sequence number%3ld\n", outgrno);
    else
      printf("  No, use as outgroup species%3ld\n", outgrno);
    printf("  T              Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f per site\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n",
           (interleaved ? "Yes" : "No"));
    printf("  0   Terminal type (IBM PC, VT52, ANSI)?  %s\n",
           (ibmpc ? "IBM PC" : ansi ? "ANSI" : vt52 ? "VT52" : "(none)"));
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
           (treeprint ? "Yes" : "No"));
    printf("  4          Print out steps in each site  %s\n",
           (stepbox ? "Yes" : "No"));
    printf("  5  Print sequences at all nodes of tree  %s\n",
           (ancseq ? "Yes" : "No"));
    printf("  6       Write out trees onto tree file?  %s\n",
           (trout ? "Yes" : "No"));
    printf(
   "\nAre these settings correct? (type Y or the letter for one to change)\n");

//    scanf("%c%*[^\n]", &ch);
//    getchar();
    ch = 'y';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (strchr("JOTUMI1234560",ch)){
      switch (ch) {
	
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
	
      case 'O':
	outgropt = !outgropt;
	if (outgropt) {
	  done1 = true;
	  do {
	    printf("Type number of the outgroup:\n");
	    scanf("%ld%*[^\n]", &outgrno);
	    getchar();
	    done1 = (outgrno >= 1 && outgrno <= spp);
	    if (!done1) {
	      printf("BAD OUTGROUP NUMBER: %4ld\n", outgrno);
	      printf("  Must be in range 1 -%2ld\n", spp);
	    }
	  } while (done1 != true);
	} else outgrno = 1;
	break;
	
      case 'T':
	thresh = !thresh;
	if (thresh) {
	  done1 = false;
	  do {
	    printf("What will be the threshold value?\n");
	    scanf("%lf%*[^\n]", &threshold);
	    getchar();
	    done1 = (threshold >= 1.0);
	    if (!done1)
	      printf("BAD THRESHOLD VALUE:  it must be greater than 1\n");
	    else
	      threshold = (long)(threshold * 10.0 + 0.5) / 10.0;
	  } while (done1 != true);
	}
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
	
      case 'I':
	interleaved = !interleaved;
	break;
	
      case 'U':
	usertree = !usertree;
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
	stepbox = !stepbox;
	break;
	
      case '5':
	ancseq = !ancseq;
	break;
	
      case '6':
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
  long i;
  node *p, *q;

  inputnumbers();
  getoptions();
  treenode = (pointptr)Malloc(nonodes*sizeof(node *));
  for (i = 0; i < (spp); i++) {
    treenode[i] = (node *)Malloc(sizeof(node));
    treenode[i]->numsteps = (steptr)Malloc(chars*sizeof(long));
    treenode[i]->siteset = (seqptr)Malloc(chars*sizeof(sitearray));
    treenode[i]->seq = (aas *)Malloc(chars*sizeof(aas));
  }
  for (i = spp; i < (nonodes); i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc(sizeof(node));
      p->numsteps = (steptr)Malloc(chars*sizeof(long));
      p->siteset = (seqptr)Malloc(chars*sizeof(sitearray));
      p->seq = (aas *)Malloc(chars*sizeof(aas));
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    treenode[i] = p;
  }
}  /* doinit*/

void inputweights()
{
  /* input the character weights, 0-9 and A-Z for weights 10 - 35 */
  Char ch;
  long i;

  for (i = 1; i < nmlngth; i++)
    ch = getc(infile);
  for (i = 0; i < (chars); i++) {
    do {
      if (eoln(infile)) {
        fscanf(infile, "%*[^\n]");
        getc(infile);
      }
      ch = getc(infile);
    } while (ch == ' ');
    weight[i] = 1;
    if (isdigit(ch))
      weight[i] = ch - '0';
    else if (isalpha(ch)) {
      uppercase(&ch);
      if (ch >= 'A' && ch <= 'I')
        weight[i] = ch - 55;
      else if (ch >= 'J' && ch <= 'R')
        weight[i] = ch - 55;
      else
        weight[i] = ch - 55;
    } else {
      printf("BAD WEIGHT CHARACTER: %c\n", ch);
      exit(-1);
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  weights = true;
}  /* inputweights */

void printweights()
{
  /* print out the weights of sites */
  long i, j, k;

  fprintf(outfile, "    Sites are weighted as follows:\n");
  fprintf(outfile, "        ");
  for (i = 0; i <= 9; i++)
    fprintf(outfile, "%3ld", i);
  fprintf(outfile, "\n     *---------------------------------\n");
  for (j = 0; j <= (chars / 10); j++) {
    fprintf(outfile, "%5ld!  ", j * 10);
    for (i = 0; i <= 9; i++) {
      k = j * 10 + i;
      if (k > 0 && k <= chars)
        fprintf(outfile, "%3ld", weight[k - 1]);
      else
        fprintf(outfile, "   ");
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* printweights */

void inputoptions()
{
  /* input the information on the options */
  Char ch;
  long extranum, i, cursp, curchs, FORLIM;

  if (!firstset) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%ld%ld", &cursp, &curchs);
    if (cursp != spp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4ld\n",ith);
      exit(-1);
    }
    chars = curchs;
  }
  extranum = 0;
  while (!(eoln(infile))) {
    ch = getc(infile);
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
  for (i = 0; i < (chars); i++)
    weight[i] = 1;
  for (i = 1; i <= extranum; i++) {
    ch = getc(infile);
    uppercase(&ch);
    if (ch != 'W') {
      printf("ERROR: INCORRECT AUXILIARY OPTIONS LINE");
      printf(" WHICH STARTS WITH %c\n", ch);
      exit(-1);
    }
    if (ch == 'W')
        inputweights();
  }
  if (weights)
    printweights();
  if (!thresh)
    threshold = spp * 3.0;
  for (i = 0; i < (chars); i++) {
    weight[i] *= 10;
    threshwt[i] = (long)(threshold * weight[i] + 0.5);
  }
}  /* inputoptions */

void inputdata()
{
  /* input the names and sequences for each species */
  long i, j, k, l, aasread, aasnew;
  Char charstate;
  boolean allread, done;
  aas aa;   /* temporary amino acid for input */
  long FORLIM, FORLIM1;

  if (progress)
    putchar('\n');
  j = nmlngth + (chars + (chars - 1) / 10) / 2 - 5;
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
  aasread = 0;
  allread = false;
  while (!(allread)) {
    allread = true;
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    i = 1;
    while (i <= spp) {
      if ((interleaved && aasread == 0) || !interleaved) {
        for (j = 0; j < nmlngth; j++) {
          if (eof(infile) || eoln(infile)){
            printf("ERROR: END-OF-LINE OR END-OF-FILE");
            printf(" IN THE MIDDLE OF A SPECIES NAME\n");
            exit(-1);
          }
          nayme[i - 1][j] = getc(infile);
        }
      }
      j = interleaved ? aasread : 0;
      done = false;
      while (!done && !eof(infile)) {
        if (interleaved)
          done = true;
        while (j < chars && !(eoln(infile) || eof(infile))) {
          charstate = getc(infile);
          if (charstate == ' ' || (charstate >= '0' && charstate <= '9'))
            continue;
          uppercase(&charstate);
          if ((!isalpha(charstate) && charstate != '.' && charstate != '?' &&
               charstate != '-' && charstate != '*') || charstate == 'J' ||
              charstate == 'O' || charstate == 'U') {
            printf("WARNING -- BAD AMINO ACID:%c",charstate);
	    printf(" AT POSITION%5ld OF SPECIES %3ld\n",j,i);
	    exit(-1);
          }
          j++;
          if (charstate == '.') {
            treenode[i - 1]->seq[j - 1] = treenode[0]->seq[j - 1];
            memcpy(treenode[i - 1]->siteset[j - 1],
                   treenode[0]->siteset[j - 1], sizeof(sitearray));
            continue;
          }
            aa =  (charstate == 'A') ?  ala :
                  (charstate == 'B') ?  asx :
                  (charstate == 'C') ?  cys :
                  (charstate == 'D') ?  asp :
                  (charstate == 'E') ?  glu :
                  (charstate == 'F') ?  phe :
                  (charstate == 'G') ?  gly : aa;
            aa =  (charstate == 'H') ?  his :
                  (charstate == 'I') ? ileu :
                  (charstate == 'K') ?  lys :
                  (charstate == 'L') ?  leu :
                  (charstate == 'M') ?  met :
                  (charstate == 'N') ?  asn :
                  (charstate == 'P') ?  pro :
                  (charstate == 'Q') ?  gln :
                  (charstate == 'R') ?  arg : aa;
             aa = (charstate == 'S') ?  ser :
                  (charstate == 'T') ?  thr :
                  (charstate == 'V') ?  val :
                  (charstate == 'W') ?  trp :
                  (charstate == 'X') ?  unk :
                  (charstate == 'Y') ?  tyr :
                  (charstate == 'Z') ?  glx :
                  (charstate == '*') ? stop :
                  (charstate == '?') ? quest:
                  (charstate == '-') ? del  :  aa;

          treenode[i - 1]->seq[j - 1] = aa;
          memcpy(treenode[i - 1]->siteset[j - 1],
                 translate[(long)aa - (long)ala], sizeof(sitearray));
        }
        if (interleaved)
          continue;
        if (j < chars) {
          fscanf(infile, "%*[^\n]");
          getc(infile);
        } else if (j == chars)
          done = true;
      }
      if (interleaved && i == 1)
        aasnew = j;
      fscanf(infile, "%*[^\n]");
      getc(infile);
      if ((interleaved && j != aasnew) || ((!interleaved) && j != chars)){
        printf("ERROR: SEQUENCES OUT OF ALIGNMENT\n");
	exit(-1);}
      i++;
    }
    if (interleaved) {
      aasread = aasnew;
      allread = (aasread == chars);
    } else
      allread = (i > spp);
  }
  if (printdata) {
    for (i = 1; i <= ((chars - 1) / 60 + 1); i++) {
      for (j = 1; j <= (spp); j++) {
        for (k = 0; k < nmlngth; k++)
          putc(nayme[j - 1][k], outfile);
        fprintf(outfile, "   ");
        l = i * 60;
        if (l > chars)
          l = chars;
        for (k = (i - 1) * 60 + 1; k <= l; k++) {
          if (j > 1 && treenode[j - 1]->seq[k - 1] == treenode[0]->seq[k - 1])
            charstate = '.';
          else {
            tmpa = treenode[j-1]->seq[k-1];
              charstate =  (tmpa == ala) ? 'A' :
                           (tmpa == asx) ? 'B' :
                           (tmpa == cys) ? 'C' :
                           (tmpa == asp) ? 'D' :
                           (tmpa == glu) ? 'E' :
                           (tmpa == phe) ? 'F' :
                           (tmpa == gly) ? 'G' :
                           (tmpa == his) ? 'H' :
                           (tmpa ==ileu) ? 'I' :
                           (tmpa == lys) ? 'K' :
                           (tmpa == leu) ? 'L' : charstate;
              charstate =  (tmpa == met) ? 'M' :
                           (tmpa == asn) ? 'N' :
                           (tmpa == pro) ? 'P' :
                           (tmpa == gln) ? 'Q' :
                           (tmpa == arg) ? 'R' :
                           (tmpa == ser) ? 'S' :
                           (tmpa ==ser1) ? 'S' :
                           (tmpa ==ser2) ? 'S' : charstate;
              charstate =  (tmpa == thr) ? 'T' :
                           (tmpa == val) ? 'V' :
                           (tmpa == trp) ? 'W' :
                           (tmpa == unk) ? 'X' :
                           (tmpa == tyr) ? 'Y' :
                           (tmpa == glx) ? 'Z' :
                           (tmpa == del) ? '-' :
                           (tmpa ==stop) ? '*' :
                           (tmpa==quest) ? '?' : charstate;
        }
          putc(charstate, outfile);
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
}  /* inputdata */

void makevalues()
{
  /* set up fractional likelihoods at tips */
  long i, j;
  node *p;

  for (i = 1; i <= nonodes; i++) {
    treenode[i - 1]->back = NULL;
    treenode[i - 1]->tip = (i <= spp);
    treenode[i - 1]->index = i;
    for (j = 0; j < (chars); j++)
      treenode[i - 1]->numsteps[j] = 0;
    if (i > spp) {
      p = treenode[i - 1]->next;
      while (p != treenode[i - 1]) {
        p->back = NULL;
        p->tip = false;
        p->index = i;
        for (j = 0; j < (chars); j++)
          p->numsteps[j] = 0;
        p = p->next;
      }
    }
  }
}  /* makevalues */

void doinput()
{
  /* reads the input data */
  inputoptions();
  inputdata();
  makevalues();
}  /* doinput */



void fillin(p, left, rt)
node *p, *left, *rt;
{
  /* sets up for each node in the tree the aa set for site m
     at that point and counts the changes.  The program
     spends much of its time in this PROCEDURE */
  boolean counted;
  aas aa;
  long s;
  sitearray ls, rs, qs;
  long i, j, k, m, n;

  for (m = 0; m < chars; m++) {
    k = 0;
    if (left != NULL)
      memcpy(ls, left->siteset[m], sizeof(sitearray));
    if (rt != NULL)
      memcpy(rs, rt->siteset[m], sizeof(sitearray));
    if (left == NULL) {
      n = rt->numsteps[m];
      memcpy(qs, rs, sizeof(sitearray));
    }
    else if (rt == NULL) {
      n = left->numsteps[m];
      memcpy(qs, ls, sizeof(sitearray));
    }
    else {
      n = left->numsteps[m] + rt->numsteps[m];
      counted = false;
      for (i = 0; i <= 5; i++) {
        if (k < 3) {
          switch (i) {

            case 0:
              s = ls[0] & rs[0];
              break;

            case 1:
              s = (ls[0] & rs[1]) | (ls[1] & rs[0]);
              break;

            case 2:
              s = (ls[0] & rs[2]) | (ls[1] & rs[1]) | (ls[2] & rs[0]);
              break;

            case 3:
              s = ls[0] | (ls[1] & rs[2]) | (ls[2] & rs[1]) | rs[0];
              break;

            case 4:
              s = ls[1] | (ls[2] & rs[2]) | rs[1];
              break;

            case 5:
              s = ls[2] | rs[2];
              break;
            }
          if (counted || s != 0) {
            qs[k] = s;
            k++;
            counted = true;
          } else if (!counted)
            n += weight[m];
        }
      }
    }
    for (i = 0; i <= 1; i++) {
      for (aa = ala; (long)aa <= (long)stop; aa = (aas)((long)aa + 1)) {
        if (((1L << ((long)aa)) & qs[i]) != 0) {
          for (j = i + 1; j <= 2; j++)
            qs[j] |= translate[(long)aa - (long)ala][j - i];
        }
      }
    }
    p->numsteps[m] = n;
    memcpy(p->siteset[m], qs, sizeof(sitearray));
  }
}  /* fillin */

void preorder(p)
node *p;
{
  /* recompute number of steps in preorder taking both ancestoral and
     descendent steps into account */
  if (p != NULL && !p->tip) {
    fillin (p->next, p->next->next->back, p->back);
    fillin (p->next->next, p->back, p->next->back);
    preorder (p->next->back);
    preorder (p->next->next->back);
  }
} /* preorder */


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
  if (recompute) {
    fillin (newfork, newfork->next->back, newfork->next->next->back);
    preorder(newfork);
    if (newfork != root)
      preorder(newfork->back);
  }
}  /* add */

void re_move(item, fork)
node **item, **fork;
{
  /* removes nodes item and its ancestor, fork, from the tree.
     the new descendant of fork's ancestor is made to be
     fork's second descendant (other than item).  Also
     returns pointers to the deleted nodes, item and fork */
  node *p, *q, *other;

  if ((*item)->back == NULL) {
    *fork = NULL;
    return;
  }
  *fork = treenode[(*item)->back->index - 1];
  if ((*item) == (*fork)->next->back)
    other = (*fork)->next->next->back;
  else other = (*fork)->next->back;
  if (root == *fork)
    root = other;
  p = (*item)->back->next->back;
  q = (*item)->back->next->next->back;
  if (p != NULL) p->back = q;
  if (q != NULL) q->back = p;
  (*fork)->back = NULL;
  p = (*fork)->next;
  do {
    p->back = NULL;
    p = p->next;
  } while (p != (*fork));
  (*item)->back = NULL;
  if (recompute) {
    preorder(other);
    if (other != root) preorder(other->back);
  }
}  /* re_move */

void evaluate(r)
node *r;
{
  /* determines the number of steps needed for a tree. this is the
     minimum number of steps needed to evolve sequences on this tree */
  long i, steps, term;
  double sum;

  sum = 0.0;
  for (i = 0; i < (chars); i++) {
    steps = r->numsteps[i];
    if (steps <= threshwt[i])
      term = steps;
    else
      term = threshwt[i];
    sum += term;
    if (usertree && which <= maxuser)
      fsteps[which - 1][i] = term;
  }
  if (usertree && which <= maxuser) {
    nsteps[which - 1] = sum;
    if (which == 1) {
      minwhich = 1;
      minsteps = sum;
    } else if (sum < minsteps) {
      minwhich = which;
      minsteps = sum;
    }
  }
  like = -sum;
}  /* evaluate */

void postorder(p)
node *p;
{
  /* traverses a binary tree, calling PROCEDURE fillin at a
     node's descendants before calling fillin at the node */
  if (p->tip)
    return;
  postorder(p->next->back);
  postorder(p->next->next->back);
  fillin(p,    p->next->back, p->next->next->back);
}  /* postorder */

void reroot(outgroup)
node *outgroup;
{
  /* reorients tree, putting outgroup in desired position. */
  node *p, *q;

  if (outgroup->back->index == root->index)
    return;
  p = root->next;
  q = root->next->next;
  p->back->back = q->back;
  q->back->back = p->back;
  p->back = outgroup;
  q->back = outgroup->back;
  outgroup->back->back = q;
  outgroup->back = p;
}  /* reroot */



void savetraverse(p, pos,found)
node *p;
long *pos;
boolean *found;
{
  /* sets BOOLEANs that indicate which way is down */
  p->bottom = true;
  if (p->tip)
    return;
  p->next->bottom = false;
  savetraverse(p->next->back, pos,found);
  p->next->next->bottom = false;
  savetraverse(p->next->next->back, pos,found);
}  /* savetraverse */

void savetree(pos,found)
long *pos;
boolean *found;
{
  /* record in place where each species has to be
     added to reconstruct this tree */
  long i, j;
  node *p;
  boolean done;

  reroot(treenode[outgrno - 1]);
  savetraverse(root, pos,found);
  for (i = 0; i < (nonodes); i++)
    place[i] = 0;
  place[root->index - 1] = 1;
  for (i = 1; i <= (spp); i++) {
    p = treenode[i - 1];
    while (place[p->index - 1] == 0) {
      place[p->index - 1] = i;
      while (!p->bottom)
        p = p->next;
      p = p->back;
    }
    if (i > 1) {
      place[i - 1] = place[p->index - 1];
      j = place[p->index - 1];
      done = false;
      while (!done) {
        place[p->index - 1] = spp + i - 1;
        while (!p->bottom)
          p = p->next;
        p = p->back;
        done = (p == NULL);
        if (!done)
          done = (place[p->index - 1] != j);
      }
    }
  }
}  /* savetree */

void findtree(pos,found)
long *pos;
boolean *found;
{
  /* finds tree given by ARRAY place in ARRAY
     bestrees by binary search */
  long i, lower, upper;
  boolean below, done;

  below = false;
  lower = 1;
  upper = nextree - 1;
  *found = false;
  while (!(*found) && lower <= upper) {
    (*pos) = (lower + upper) / 2;
    i = 3;
    done = false;
    while (!done) {
      done = (i > spp);
      if (!done)
        done = (place[i - 1] != bestrees[(*pos) - 1][i - 1]);
      if (!done)
        i++;
    }
    (*found) = (i > spp);
    below = (place[i - 1] < bestrees[(*pos) - 1][i - 1]);
    if (*found)
      break;
    if (below)
      upper = (*pos) - 1;
    else
      lower = (*pos) + 1;
  }
  if (!(*found) && !below)
    (*pos)++;
}  /* findtree */

void addtree(pos,found)
long *pos;
boolean *found;
{
  /* puts tree from ARRAY place in its proper position
     in ARRAY bestrees */
  long i, FORLIM;

  for (i = nextree - 1; i >= (*pos); i--)
    memcpy(bestrees[i], bestrees[i - 1], spp*sizeof(long));
  for (i = 0; i < (spp); i++)
    bestrees[(*pos) - 1][i] = place[i];
  nextree++;
}  /* addtree */

void tryadd(p, item,nufork)
node *p,**item,**nufork;
{
  /* temporarily adds one fork and one tip to the tree.
     if the location where they are added yields greater
     "likelihood" than other locations tested up to that
     time, then keeps that location as there */
/* Local variables for tryadd: */
  long pos;
  boolean found;
  node *rute, *q;

  if (p == root)
    fillin(temp, *item, p);
  else {
    fillin(temp1, *item, p);
    fillin(temp, temp1, p->back);
  }
  evaluate(temp);
  if (lastrearr) {
    if (like < bestlike) {
      if ((*item) == (*nufork)->next->next->back) {
        q = (*nufork)->next;
        (*nufork)->next = (*nufork)->next->next;
        (*nufork)->next->next = q;
        q->next = (*nufork);
      }
    }
    else if (like >= bstlike2) {
      recompute = false;
      add (p, (*item), (*nufork));
      rute = root->next->back;
      savetree(&pos,&found);
      reroot(rute);
      if (like > bstlike2) {
        bestlike = bstlike2 = like;
        pos = 1;
        nextree = 1;
        addtree(&pos,&found);
      } else {
        pos = 0;
        findtree(&pos,&found);
        if (!found) {
          if (nextree <= maxtrees)
            addtree(&pos,&found);
        }
      }
      re_move (item, nufork);
      recompute = true;
    }
  }
  if (like > bestyet) {
    bestyet = like;
    there = p;
  }
}  /* tryadd */

void addpreorder(p, item, nufork)
node *p, *item, *nufork;
{
  /* traverses a binary tree, calling PROCEDURE tryadd
     at a node before calling tryadd at its descendants */

  if (p == NULL)
    return;
  tryadd(p, &item,&nufork);
  if (!p->tip) {
    addpreorder(p->next->back, item, nufork);
    addpreorder(p->next->next->back, item, nufork);
  }
}  /* addpreorder */

void tryrearr(p, success)
node *p;
boolean *success;
{
  /* evaluates one rearrangement of the tree.
     if the new tree has greater "likelihood" than the old
     one sets success := TRUE and keeps the new tree.
     otherwise, restores the old tree */
  node *frombelow, *whereto, *forknode, *q;
  double oldlike;

  if (p->back == NULL)
    return;
  forknode = treenode[p->back->index - 1];
  if (forknode->back == NULL)
    return;
  oldlike = bestyet;
  if (p->back->next->next == forknode)
    frombelow = forknode->next->next->back;
  else
    frombelow = forknode->next->back;
  whereto = treenode[forknode->back->index - 1];
  if (whereto->next->back == forknode)
    q = whereto->next->next->back;
  else
    q = whereto->next->back;
  fillin(temp1, frombelow, q);
  fillin(temp, temp1, p);
  fillin(temp1, temp, whereto->back);
  evaluate(temp1);
  if (like <= oldlike) {
    if (p == forknode->next->next->back) {
      q = forknode->next;
      forknode->next = forknode->next->next;
      forknode->next->next = q;
      q->next = forknode;
    }
  }
  else {
    recompute = false;
    re_move(&p, &forknode);
    fillin(whereto, whereto->next->back, whereto->next->next->back);
    recompute = true;
    add(whereto, p, forknode);
    *success = true;
    bestyet = like;
  }
}  /* tryrearr */

void repreorder(p,success)
node *p;
boolean *success;
{
  /* traverses a binary tree, calling PROCEDURE tryrearr
     at a node before calling tryrearr at its descendants */
  if (p == NULL)
    return;
  tryrearr(p,success);
  if (!p->tip) {
    repreorder(p->next->back,success);
    repreorder(p->next->next->back,success);
  }
}  /* repreorder */

void rearrange(r)
node **r;
{
  /* traverses the tree (preorder), finding any local
     rearrangement which decreases the number of steps.
     if traversal succeeds in increasing the tree's
     "likelihood", PROCEDURE rearrange runs traversal again */
  boolean success = true;
  while (success) {
    success = false;
    repreorder(*r, &success);
  }
}  /* rearrange */


void findch(c)
Char c;
{
  /* scan forward until find character c */
  boolean done;

  done = false;
  while (!done) {
    if (c == ',') {
      if (ch == '(' || ch == ')' ||ch == ';') {
        printf("\nERROR IN USER TREE:");
	printf(" UNMATCHED PARENTHESIS OR MISSING COMMA\n");
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
    if ((done && ch == ')') || !done) {
      if (eoln(infile)) {
        fscanf(infile, "%*[^\n]");
        getc(infile);
      }
      ch = getc(infile);
    }
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
  Char str[nmlngth];

  do {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    ch = getc(infile);
  } while (ch == ' ');
  if (ch == '(' ) {
    if ((*lparens) >= spp - 1) {
      printf("\nERROR IN USER TREE: TOO MANY LEFT PARENTHESES\n");
      exit(-1);
    }
    (*nextnode)++;
    (*lparens)++;
    q = treenode[(*nextnode) - 1];
    addelement(&q->next->back, nextnode,lparens,names);
    q->next->back->back = q->next;
    findch(',');
    addelement(&q->next->next->back, nextnode,lparens,names);
    q->next->next->back->back = q->next->next;
    findch(')');
    *p = q;
    return;
  }
  for (i = 0; i < nmlngth; i++)
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
  } while (ch != ',' && ch != ')' && ch != ':' && n <= nmlngth);
  n = 1;
  do {
    found = true;
    for (i = 0; i < nmlngth; i++)
      found = (found && str[i] == nayme[n - 1][i]);
    if (found) {
      if (names[n - 1] == false) {
        *p = treenode[n - 1];
        names[n - 1] = true;
      } else {
        printf("\nERROR IN USER TREE: DUPLICATE NAME FOUND -- ");
        for (i = 0; i < nmlngth; i++)
          putchar(nayme[n - 1][i]);
        putchar('\n');
	exit(-1);
      }
    } else
      n++;
  } while (!(n > spp || found));
  if (n <= spp)
    return;
  printf("CANNOT FIND SPECIES: ");
  for (i = 0; i < nmlngth; i++)
    putchar(str[i]);
  putchar('\n');
}  /* addelement */

void treeread()
{
  /* read in user-defined tree and set it up */
  long nextnode, lparens;
  long i;

  root = treenode[spp];
  nextnode = spp;
  root->back = NULL;
  for (i = 0; i < (spp); i++)
    names[i] = false;
  lparens = 0;
  addelement(&root, &nextnode,&lparens,names);
  findch(';');
  if (progress)
    printf("\n\n");
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
    p->ycoord = (*tipy);
    p->ymin = (*tipy);
    p->ymax = (*tipy);
    (*tipy) += down;
    return;
  }
  q = p->next;
  do {
    coordinates(q->back,tipy);
    q = q->next;
  } while (p != q);
  first = p->next->back;
  q = p->next;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = last->ymax - first->ymin;
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
    if (p->index - spp >= 10)
      fprintf(outfile, "-%2ld", p->index - spp);
    else
      fprintf(outfile, "--%ld", p->index - spp);
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
      putc('+', outfile);
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        if (q->index - spp >= 10)
          fprintf(outfile, "%2ld", q->index - spp);
        else
          fprintf(outfile, "-%ld", q->index - spp);
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
    if (p != q)
      p = q;
  } while (!done);
  if (p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], outfile);
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

  putc('\n', outfile);
  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  coordinates(root, &tipy);
  scale = 1.5;
  putc('\n', outfile);
  for (i = 1; i <= (tipy - down); i++)
    drawline(i, scale);
  fprintf(outfile, "\n  remember:");
  if (outgropt)
    fprintf(outfile, " (although rooted by outgroup)");
  fprintf(outfile, " this is an unrooted tree!\n\n");
}  /* printree */



void ancestset(a, b, c, d, k)
long *a, *b, *c, *d;
long *k;
{
  /* sets up the aa set array. */
  aas aa;
  long s, sa, sb;
  long i, j, m, n;
  boolean counted;

  counted = false;
  *k = 0;
  for (i = 0; i <= 5; i++) {
    if (*k < 3) {
      s = 0;
      if (i > 3)
        n = i - 3;
      else
        n = 0;
      for (j = n; j <= (i - n); j++) {
        if (j < 3)
          sa = a[j];
        else
          sa = fullset;
        for (m = n; m <= (i - j - n); m++) {
          if (m < 3)
            sb = sa & b[m];
          else
            sb = sa;
          if (i - j - m < 3)
            sb &= c[i - j - m];
          s |= sb;
        }
      }
      if (counted || s != 0) {
        d[*k] = s;
        (*k)++;
        counted = true;
      }
    }
  }
  for (i = 0; i <= 1; i++) {
    for (aa = ala; (long)aa <= (long)stop; aa = (aas)((long)aa + 1)) {
      if (((1L << ((long)aa)) & d[i]) != 0) {
        for (j = i + 1; j <= 2; j++)
          d[j] |= translate[(long)aa - (long)ala][j - i];
      }
    }
  }
}  /* ancestset */


void hyprint(b1, b2,bottom,r,nonzero,maybe,nothing)
long b1, b2;
boolean *bottom,*nonzero,*maybe;
node *r;
sitearray nothing;
{
  /* print out states in sites b1 through b2 at node */
  long i;
  boolean dot;
  Char ch;
  aas aa;

  if (*bottom) {
    if (!outgropt)
      fprintf(outfile, "      ");
    else
      fprintf(outfile, "root  ");
  } else
    fprintf(outfile, "%3ld   ", r->back->index - spp);
  if (r->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[r->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", r->index - spp);
  if (*bottom)
    fprintf(outfile, "          ");
  else if (*nonzero)
    fprintf(outfile, "   yes    ");
  else if (*maybe)
    fprintf(outfile, "  maybe   ");
  else
    fprintf(outfile, "   no     ");
  for (i = b1 - 1; i < b2; i++) {
    aa = r->seq[i];
    switch (aa) {

    case ala:
      ch = 'A';
      break;

    case asx:
      ch = 'B';
      break;

    case cys:
      ch = 'C';
      break;

    case asp:
      ch = 'D';
      break;

    case glu:
      ch = 'E';
      break;

    case phe:
      ch = 'F';
      break;

    case gly:
      ch = 'G';
      break;

    case his:
      ch = 'H';
      break;

    case ileu:
      ch = 'I';
      break;

    case lys:
      ch = 'K';
      break;

    case leu:
      ch = 'L';
      break;

    case met:
      ch = 'M';
      break;

    case asn:
      ch = 'N';
      break;

    case pro:
      ch = 'P';
      break;

    case gln:
      ch = 'Q';
      break;

    case arg:
      ch = 'R';
      break;

    case ser:
      ch = 'S';
      break;

    case ser1:
      ch = 'S';
      break;

    case ser2:
      ch = 'S';
      break;

    case thr:
      ch = 'T';
      break;

    case trp:
      ch = 'W';
      break;

    case tyr:
      ch = 'Y';
      break;

    case val:
      ch = 'V';
      break;

    case glx:
      ch = 'Z';
      break;

    case del:
      ch = '-';
      break;

    case stop:
      ch = '*';
      break;

    case unk:
      ch = 'X';
      break;

    case quest:
      ch = '?';
      break;
    }
    if (!(*bottom))
      dot = (r->siteset[i] [0] == treenode[r->back->index - 1]->siteset[i][0]
             || ((r->siteset[i][0] &
                  (~((1L << ((long)ser1)) | (1L << ((long)ser2)) |
                                         (1L << ((long)ser))))) == 0 &&
                (treenode[r->back->index - 1]->siteset[i] [0] &
                  (~((1L << ((long)ser1)) | (1L << ((long)ser2)) |
                                         (1L << ((long)ser))))) == 0));
    else
      dot = false;
    if (dot)
      putc('.', outfile);
    else
      putc(ch, outfile);
    if ((i + 1) % 10 == 0)
      putc(' ', outfile);
  }
  putc('\n', outfile);
}  /* hyprint */

void hyptrav(r, hypset, b1, b2, k,bottom,nothing)
node *r;
sitearray *hypset;
long b1, b2,*k;
boolean *bottom;
sitearray nothing;
{
  boolean maybe, nonzero;
  long i;
  aas aa;
  long anc, hset;
  gseq *ancset, *temparray;

  gnu(&ancset);
  gnu(&temparray);
  maybe = false;
  nonzero = false;
  for (i = b1 - 1; i < b2; i++) {
    if (!r->tip) {
      ancestset(hypset[i], r->next->back->siteset[i],
                r->next->next->back->siteset[i], temparray->seq[i], k);
      memcpy(r->siteset[i], temparray->seq[i], sizeof(sitearray));
    }
    if (!(*bottom))
      anc = treenode[r->back->index - 1]->siteset[i][0];
    if (!r->tip) {
      hset = r->siteset[i][0];
      r->seq[i] = quest;
      for (aa = ala; (long)aa <= (long)stop; aa = (aas)((long)aa + 1)) {
        if (hset == 1L << ((long)aa))
          r->seq[i] = aa;
      }
      if (hset == ((1L << ((long)asn)) | (1L << ((long)asp))))
        r->seq[i] = asx;
      if (hset == ((1L << ((long)gln)) | (1L << ((long)gly))))
        r->seq[i] = glx;
      if (hset == ((1L << ((long)ser1)) | (1L << ((long)ser2))))
        r->seq[i] = ser;
      if (hset == fullset)
        r->seq[i] = unk;
    }
    nonzero = (nonzero || (r->siteset[i][0] & anc) == 0);
    maybe = (maybe || r->siteset[i][0] != anc);
  }
  hyprint(b1, b2,bottom,r,&nonzero,&maybe,nothing);
  *bottom = false;
  if (!r->tip) {
    memcpy(temparray->seq, r->next->back->siteset, chars*sizeof(sitearray));
    for (i = b1 - 1; i < b2; i++)
      ancestset(hypset[i], r->next->next->back->siteset[i], nothing,
                ancset->seq[i], k);
    hyptrav(r->next->back, ancset->seq, b1, b2,k,bottom,nothing );
    for (i = b1 - 1; i < b2; i++)
      ancestset(hypset[i], temparray->seq[i], nothing, ancset->seq[i],k);
    hyptrav(r->next->next->back, ancset->seq, b1, b2, k,bottom,nothing);
  }
  chuck(temparray);
  chuck(ancset);
}  /* hyptrav */

void hypstates(k)
long *k;
{
  /* fill in and describe states at interior nodes */
  boolean bottom;
  sitearray nothing;
  long i, n;
  seqptr hypset;

  fprintf(outfile, "\nFrom    To     Any Steps?    State at upper node\n");
  fprintf(outfile, "                             ");
  fprintf(outfile, "( . means same as in the node below it on tree)\n\n");
  memcpy(nothing, translate[(long)quest - (long)ala], sizeof(sitearray));
  hypset = (seqptr)Malloc(chars*sizeof(sitearray));
  for (i = 0; i < (chars); i++)
    memcpy(hypset[i], nothing, sizeof(sitearray));
  bottom = true;
  for (i = 1; i <= ((chars - 1) / 40 + 1); i++) {
    putc('\n', outfile);
    n = i * 40;
    if (n > chars)
      n = chars;
    bottom = true;
    hyptrav(root, hypset, i * 40 - 39, n, k,&bottom,nothing);
  }
  free(hypset);
}  /* hypstates */

void treeout(p)
node *p;
{
  /* write out file with representation of final tree */
  long i, n;
  Char c;

  if (p->tip) {
    n = 0;
    for (i = 1; i <= nmlngth; i++) {
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = nayme[p->index - 1][i];
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
    if (col > 65) {
      putc('\n', treefile);
      col = 0;
    }
    treeout(p->next->next->back);
    putc(')', treefile);
    col++;
  }
  if (p != root)
    return;
  if (nextree > 2)
    fprintf(treefile, "[%6.4f];\n", 1.0 / (nextree - 1));
  else
    fprintf(treefile, ";\n");
}  /* treeout */

void describe()
{
  /* prints ancestors, steps and table of numbers of steps in
     each site */
  long i,j,k;

  if (treeprint)
    fprintf(outfile, "\nrequires a total of %10.3f\n", like / -10);
  if (stepbox) {
    putc('\n', outfile);
    if (weights)
      fprintf(outfile, " weighted");
    fprintf(outfile, "steps in each position:\n");
    fprintf(outfile, "      ");
    for (i = 0; i <= 9; i++)
      fprintf(outfile, "%4ld", i);
    fprintf(outfile, "\n     *-----------------------------------------\n");
    for (i = 0; i <= (chars / 10); i++) {
      fprintf(outfile, "%5ld", i * 10);
      putc('!', outfile);
      for (j = 0; j <= 9; j++) {
        k = i * 10 + j;
        if (k == 0 || k > chars)
          fprintf(outfile, "    ");
        else
          fprintf(outfile, "%4ld", root->numsteps[k - 1] / 10);
      }
      putc('\n', outfile);
    }
  }
  if (ancseq) {
    hypstates(&k);
    putc('\n', outfile);
  }
  putc('\n', outfile);
  if (trout) {
    col = 0;
    treeout(root);
  }
}  /* describe */


void maketree()
{
  /* constructs a binary tree from the pointers in treenode.
     adds each node at location which yields highest "likelihood"
     then rearranges the tree for greatest "likelihood" */
  long i, j, k, numtrees, num;
  double gotlike, wt, sumw, sum, sum2, sd;
  node *item, *nufork, *dummy;
  double TEMP;

  if (!usertree) {
    for (i = 1; i <= (spp); i++)
      enterorder[i - 1] = i;
    if (jumble) {
      for (i = 0; i < (spp); i++) {
        j = (long)(randum(seed) * spp) + 1;
        k = enterorder[j - 1];
        enterorder[j - 1] = enterorder[i];
        enterorder[i] = k;
      }
    }
    root = treenode[enterorder[0] - 1];
    recompute = true;
    add(treenode[enterorder[0] - 1], treenode[enterorder[1] - 1],
        treenode[spp]);
    if (progress) {
      printf("Adding species:\n");
      printf("   ");
      for (i = 0; i < nmlngth; i++)
        putchar(nayme[enterorder[0] - 1][i]);
      printf("\n   ");
      for (i = 0; i < nmlngth; i++)
        putchar(nayme[enterorder[1] - 1][i]);
      putchar('\n');
    }
    lastrearr = false;
    for (i = 3; i <= (spp); i++) {
      bestyet = -1.0e6;
      there = root;
      item = treenode[enterorder[i - 1] - 1];
      nufork = treenode[spp + i - 2];
      addpreorder(root, item, nufork);
      add(there, item, nufork);
      like = bestyet;
      rearrange(&root);
      if (progress) {
        printf("   ");
        for (j = 0; j < nmlngth; j++)
          putchar(nayme[enterorder[i - 1] - 1][j]);
        putchar('\n');
      }
      lastrearr = (i == spp);
      if (lastrearr) {
        if (progress) {
          printf("\nDoing global rearrangements\n");
          printf("  !");
          for (j = 1; j <= nonodes; j++)
            putchar('-');
          printf("!\n");
        }
        bestlike = bestyet;
        if (jumb == 1) {
          bstlike2 = bestlike;
          nextree = 1;
        }
        do {
          if (progress)
            printf("   ");
          gotlike = bestlike;
          for (j = 0; j < (nonodes); j++) {
            bestyet = -1.0e6;
            item = treenode[j];
            if (item != root) {
              nufork = treenode[treenode[j]->back->index - 1];
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
      }
    }
    if (progress)
      putchar('\n');
    for (i = spp - 1; i >= 1; i--)
      re_move(&treenode[i], &dummy);
    if (jumb == njumble) {
      if (treeprint) {
        putc('\n', outfile);
        if (nextree == 2)
          fprintf(outfile, "One most parsimonious tree found:\n");
        else
          fprintf(outfile, "%6ld trees in all found\n", nextree - 1);
      }
      if (nextree > maxtrees + 1) {
        if (treeprint)
          fprintf(outfile, "here are the first%4ld of them\n", (long)maxtrees);
        nextree = maxtrees + 1;
      }
      if (treeprint)
        putc('\n', outfile);
      recompute = false;
      for (i = 0; i <= (nextree - 2); i++) {
        root = treenode[0];
        add(treenode[0], treenode[1], treenode[spp]);
        for (j = 3; j <= (spp); j++)
          add(treenode[bestrees[i][j - 1] - 1], treenode[j - 1],
              treenode[spp + j - 2]);
        reroot(treenode[outgrno - 1]);
        postorder(root);
        evaluate(root);
        printree();
        describe();
        for (j = 1; j < (spp); j++)
          re_move(&treenode[j], &dummy);
      }
    }
  } else {
    fscanf(infile, "%ld%*[^\n]", &numtrees);
    getc(infile);
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n\n\n");
    }
    which = 1;
    while (which <= numtrees) {
      treeread();
      if (outgropt)
        reroot(treenode[outgrno - 1]);
      postorder(root);
      evaluate(root);
      printree();
      describe();
      which++;
    }
    putc('\n', outfile);
    if (numtrees > 1 && chars > 1 ) {
      fprintf(outfile, "Tree    Steps   Diff Steps   Its S.D.");
      fprintf(outfile, "   Significantly worse?\n\n");
      if (numtrees > maxuser)
        num = maxuser;
      else
        num = numtrees;
      which = 1;
      while (which <= num) {
        fprintf(outfile, "%3ld%10.1f", which, nsteps[which - 1] / 10);
        if (minwhich == which)
          fprintf(outfile, "  <------ best\n");
        else {
          sumw = 0.0;
          sum = 0.0;
          sum2 = 0.0;
          for (j = 0; j < (chars); j++) {
            if (weight[j] > 0) {
              wt = weight[j] / 10.0;
              sumw += wt;
              sum += (fsteps[which - 1][j] -
                      fsteps[minwhich - 1][j]) / 10.0;
              TEMP = (fsteps[which - 1][j] -
                      fsteps[minwhich - 1][j]) / 10.0;
              sum2 += TEMP * TEMP / wt;
            }
          }
          TEMP = sum / sumw;
          sd = sqrt(sumw / (sumw - 1.0) * (sum2 - TEMP * TEMP));
          fprintf(outfile, "%10.1f%12.4f",
                  (nsteps[which - 1] - minsteps) / 10, sd);
          if (sum > 1.95996 * sd)
            fprintf(outfile, "           Yes\n");
          else
            fprintf(outfile, "           No\n");
        }
        which++;
      }
      fprintf(outfile, "\n\n");
    }
  }
  if (jumb == njumble && progress) {
    printf("Output written to output file\n\n");
    if (trout)
      printf("Trees also written onto file\n\n");
  }
}  /* maketree */


main(argc, argv)
int argc;
Char *argv[];
{  /* Protein parsimony by uphill search */
char infilename[100],outfilename[100],trfilename[100];
#ifdef MAC
  macsetup("Protpars","");
  argv[0] = "Protpars";
#endif
//   openfile(&infile,INFILE,"r",argv[0],infilename);
// Aamer Modified Input File Line To Take From Argument
  if( argc < 2 ) {
    fprintf( stderr, "Usage: %s <inputfile name>\n", argv[0]);
    exit(0);
  }
  openfile(&infile,argv[1],"r",argv[0],infilename);

  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);

  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  garbage = NULL;
  mulsets = false;
  datasets = 1;
  firstset = true;
  setup();
  doinit();
  if (usertree) {
    fsteps = (long **)Malloc(maxuser*sizeof(long *));
    for (j = 1; j <= maxuser; j++)
      fsteps[j - 1] = (long *)Malloc(chars*sizeof(long));
  }
  bestrees = (long **)Malloc(maxtrees*sizeof(long *));
  for (j = 1; j <= maxtrees; j++)
    bestrees[j - 1] = (long *)Malloc(spp*sizeof(long));
  nayme = (Char **)Malloc(spp*sizeof(Char *));
  for (j = 1; j <= spp; j++)
    nayme[j - 1] = (Char *)Malloc(nmlngth*sizeof(Char));
  enterorder = (long *)Malloc(spp*sizeof(long));
  names = (boolean *)Malloc(spp*sizeof(boolean));
  place = (long *)Malloc(nonodes*sizeof(long));
  weight = (steptr)Malloc(chars*sizeof(long));
  threshwt = (steptr)Malloc(chars*sizeof(long));
  temp = (node *)Malloc(sizeof(node));
  temp->numsteps = (steptr)Malloc(chars*sizeof(long));
  temp->siteset = (seqptr)Malloc(chars*sizeof(sitearray));
  temp->seq = (aas *)Malloc(chars*sizeof(aas));
  temp1 = (node *)Malloc(sizeof(node));
  temp1->numsteps = (steptr)Malloc(chars*sizeof(long));
  temp1->siteset = (seqptr)Malloc(chars*sizeof(sitearray));
  temp1->seq = (aas *)Malloc(chars*sizeof(aas));
  if (trout)
    openfile(&treefile,TREEFILE,"w",argv[0],trfilename);
  for (ith = 1; ith <= datasets; ith++) {
    doinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n\n",ith);
      if (progress)
        printf("Data set # %ld:\n\n",ith);
    }
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
}  /* Protein parsimony by uphill search */


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
