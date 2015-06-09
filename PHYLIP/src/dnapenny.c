#include "phylip.h"

/* version 3.56c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define nmlngth         10    /* number of characters in species name        */
#define maxtrees        100   /* maximum number of trees to be printed out   */
#define often           100   /* how often to notify how many trees examined */
#define many            1000  /* how many multiples of howoften before stop  */

#define ibmpc0          false
#define ansi0           true
#define vt520           false
#define down            2

typedef short *steptr;
typedef short *stepshortptr;
typedef Char **sequence;

typedef enum {
  A, C, G, U, O
} bases;
typedef char *baseptr;

typedef struct gbases {
  baseptr base;
  struct gbases *next;
} gbases;

/* nodes will form a binary tree */

typedef struct node {         /* describes a tip species or an ancestor */
  struct node *next, *back;
  short index;
  boolean tip;                /* present species are tips of tree       */
  baseptr base;             /* the sequence                           */
  stepshortptr numsteps;    /* bookkeeps steps                        */
  short xcoord, ycoord, ymin;  /* used by printree                       */
  short ymax;
} node;

typedef node **pointptr;

/* Local variables for hyptrav: */
struct LOC_hyptrav {
  node *r;
  char *hypset;
  boolean maybe, nonzero;
  short tempset, anc;
} ;

Static node *root, *p;
Static FILE *infile, *outfile, *treefile;
Static short spp, nonodes, chars, endsite, outgrno, howmany, howoften, col,
	    datasets, ith, i, j;
/* spp = number of species
   nonodes = number of nodes in tree
   chars = number of sites in actual sequences
   outgrno indicates outgroup */
Static boolean weights, thresh, simple, trout, outgropt,  printdata,
	       progress, treeprint, stepbox, ancseq, mulsets, firstset,
	       interleaved, ibmpc, vt52, ansi;
Static short *weight, *oldweight;
Static steptr alias, ally, location;
Static pointptr treenode;   /* pointers to all nodes in tree */
Static Char **name;   /* names of species */
Static sequence y;
Static double threshold, fracdone, fracinc;
Static boolean *added;
Static gbases *garbage;
Static Char basechar[]="ACMGRSVTWYHKDBNO???????????????";

typedef short *treenumbers;
typedef double *valptr;
typedef short *placeptr;

/* Variables for maketree, propagated globally for C version: */
short nextree, examined, mults;
boolean firsttime, done, recompute;
double like, bestyet;
treenumbers *bestorders, *bestrees;
treenumbers current, order;
long *threshwt;
baseptr nothing;
node *temp, *temp1;
static short suppno[] =
  { 0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,4};

static char suppset[] =          /* this was previously a function. */
{                                  /* in C, it doesn't need to be.    */
   1 << ((Char)A),
   1 << ((Char)C),
   1 << ((Char)G),
   1 << ((Char)U),
   1 << ((Char)O),
   (1 << ((Char)A)) | (1 << ((Char)C)),
   (1 << ((Char)A))|(1 << ((Char)G)),
   (1 << ((Char)A))|(1 << ((Char)U)),
   (1 << ((Char)A)) | (1 << ((Char)O)),
   (1 << ((Char)C)) | (1 << ((Char)G)),
   (1 << ((Char)C)) | (1 << ((Char)U)),
   (1 << ((Char)C)) | (1 << ((Char)O)),
   (1 << ((Char)G)) | (1 << ((Char)U)),
   (1 << ((Char)G)) | (1 << ((Char)O)),
   (1 << ((Char)U)) | (1 << ((Char)O)),
   (1 << ((Char)A)) | (1 << ((Char)C)) | (1 << ((Char)G)),
   (1 << ((Char)A)) | (1 << ((Char)C)) | (1 << ((Char)U)),
   (1 << ((Char)A)) | (1 << ((Char)C)) | (1 << ((Char)O)),
   (1 << ((Char)A)) | (1 << ((Char)G)) | (1 << ((Char)U)),
   (1 << ((Char)A)) | (1 << ((Char)G)) | (1 << ((Char)O)),
   (1 << ((Char)A)) | (1 << ((Char)U)) | (1 << ((Char)O)),
   (1 << ((Char)C)) | (1 << ((Char)G)) | (1 << ((Char)U)),
   (1 << ((Char)C)) | (1 << ((Char)G)) | (1 << ((Char)O)),
   (1 << ((Char)C)) | (1 << ((Char)U)) | (1 << ((Char)O)),
   (1 << ((Char)G)) | (1 << ((Char)U)) | (1 << ((Char)O)),
   (1 << ((Char)A))|(1 << ((Char)C))|(1 << ((Char)G))|(1 << ((Char)U)),
   (1 << ((Char)A))|(1 << ((Char)C))|(1 << ((Char)G))|(1 << ((Char)O)),
   (1 << ((Char)A))|(1 << ((Char)C))|(1 << ((Char)U))|(1 << ((Char)O)),
   (1 << ((Char)A))|(1 << ((Char)G))|(1 << ((Char)U))|(1 << ((Char)O)),
   (1 << ((Char)C))|(1 << ((Char)G))|(1 << ((Char)U)) | (1 << ((Char)O)),
   (1 << ((Char)A))|(1 << ((Char)C))|(1 << ((Char)G)) | (1 << ((Char)U)) | (1 << ((Char)O))};


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



Static Void gnu(p)
gbases **p;
{
  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */
  if (garbage != NULL) {
    *p = garbage;
    garbage = garbage->next;
  } else {
    *p = (gbases *)Malloc(sizeof(gbases));
    (*p)->base = (baseptr)Malloc(endsite*sizeof(char));
  }
  (*p)->next = NULL;
}  /* gnu */


Static Void chuck(p)
gbases *p;
{
  /* collect garbage on p -- put it on front of garbage list */
  p->next = garbage;
  garbage = p;
}  /* chuck */


Static Void uppercase(ch)
Char *ch;
{
    *ch = (islower (*ch) ? toupper(*ch) : (*ch));
}  /* uppercase */



Local Void inputnumbers()
{
  /* input the numbers of species and of characters */
  fscanf(infile, "%hd%hd", &spp, &chars);
  if (printdata)
    fprintf(outfile, "%2hd species, %3hd  sites\n", spp, chars);
  if (printdata)
    putc('\n', outfile);
  if (progress)
    putchar('\n');
  nonodes = spp * 2 - 1;
}  /* inputnumbers */

Local Void getoptions()
{
  /* interactively set options */
  Char ch;
  boolean  done1;

  fprintf(outfile, "\nPenny algorithm for DNA, version %s\n",VERSION);
  fprintf(outfile, " branch-and-bound to find all");
  fprintf(outfile, " most parsimonious trees\n\n");
  howoften = often;
  howmany = many;
  outgrno = 1;
  outgropt = false;
  simple = true;
  thresh = false;
  threshold = spp;
  trout = true;
  weights = false;
  printdata = false;
  progress = true;
  treeprint = true;
  stepbox = false;
  ancseq = false;
  interleaved = true;
  for (;;) {
    printf(ansi ? "\033[2J\033[H" : vt52 ? "\033E\033H" : "\n");
    printf("\nPenny algorithm for DNA, version %s\n",VERSION);
    printf(" branch-and-bound to find all most parsimonious trees\n\n");
    printf("Settings for this run:\n");
    printf("  H        How many groups of %4hd trees:%6hd\n", howoften, howmany);
    printf("  F        How often to report, in trees:  %4hd\n", howoften);
    printf("  S           Branch and bound is simple?  %s\n",
           (simple ?  "Yes" : "No. reconsiders order of species"));
    printf("  O                        Outgroup root?  %s%3hd\n",
           (outgropt ? "Yes, at sequence number" :
                       "No, use as outgroup species"),outgrno);
    printf("  T              Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f per site\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2hd sets\n", datasets);
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n",
           (interleaved ? "Yes" : "No, sequential"));
    printf("  0   Terminal type (IBM PC, VT52, ANSI)?  %s\n",
           (ibmpc ? "IBM PC" : ansi ? "ANSI" : vt52 ? "VT52" : "(none)"));
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
           (treeprint ? "Yes" : "No"));
    printf("  4          Print out steps in each site  %s\n",
           (stepbox ? "Yes" : "No" ));
    printf("  5  Print sequences at all nodes of tree  %s\n",
           (ancseq ? "Yes" : "No"));
    printf("  6       Write out trees onto tree file?  %s\n",
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
    if ((strchr("HMSOFTI1234560",ch)) != NULL){
      switch (ch) {
	
      case 'H':
	do {
	  printf("How many cycles of %4hd trees?\n", howoften);
	  scanf("%hd%*[^\n]", &howmany);
	  getchar();
	} while (howmany <= 0);
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
	
      case 'F':
	do {
	  printf("How many trees per cycle?\n");
	  scanf("%hd%*[^\n]", &howoften);
	  getchar();
	} while (howoften <= 0);
	break;
	
      case 'S':
	simple = !simple;
	break;
	
      case 'O':
	outgropt = !outgropt;
	if (outgropt) {
	  done1 = true;
	  do {
	    printf("Type number of the outgroup:\n");
	    scanf("%hd%*[^\n]", &outgrno);
	    getchar();
	    done1 = (outgrno >= 1 && outgrno <= spp);
	    if (!done1) {
	      printf("BAD OUTGROUP NUMBER: %4hd\n", outgrno);
	      printf("  Must be in range 1 -%2hd\n", spp);
	    }
	  } while (done1 != true);
	}
	else
	   outgrno = 1;
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
	      threshold = (short)(threshold * 10.0 + 0.5) / 10.0;
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


Static Void doinit()
{
  /* initializes variables */
  short i, j;
  node *p, *q;

  inputnumbers();
  getoptions();
  y = (Char **)Malloc(spp*sizeof(Char *));
  for (i = 0; i < spp; i++)
    y[i] = (Char *)Malloc(chars*sizeof(Char));
  treenode = (pointptr)Malloc(nonodes*sizeof(node *));
  for (i = 0; i < spp; i++)
    treenode[i] = (node *)Malloc(sizeof(node));
  for (i = spp; i < nonodes; i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc(sizeof(node));
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    treenode[i] = p;
  }
}  /* doinit*/



Local Void inputweights()
{
  /* input the character weights, 0-9 and A-Z for weights
     0 - 35 */
  Char ch;
  short i;

  for (i = 1; i < nmlngth; i++) {
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
  }
  for (i = 0; i < chars; i++) {
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

Local Void printweights()
{
  /* print out the weights of sites */
  short i, j, k;

  fprintf(outfile, "    Sites are weighted as follows:\n");
  fprintf(outfile, "        ");
  for (i = 0; i <= 9; i++)
    fprintf(outfile, "%3hd", i);
  fprintf(outfile, "\n     *---------------------------------\n");
  for (j = 0; j <= (chars / 10); j++) {
    fprintf(outfile, "%5hd!  ", j * 10);
    for (i = 0; i <= 9; i++) {
      k = j * 10 + i;
      if (k > 0 && k <= chars)
	fprintf(outfile, "%3hd", weight[k - 1]);
      else
	fprintf(outfile, "   ");
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* printweights */

Local Void inputoptions()
{
  /* input the information on the options */
  Char ch;
  short extranum, i, cursp, curchs;

  if (!firstset) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%hd%hd", &cursp, &curchs);
    if (cursp != spp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4hd\n",ith);
      exit(-1);
    }
    chars = curchs;
  }
  weights = false;
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
  for (i = 0; i < chars; i++)
    weight[i] = 1;
  for (i = 1; i <= extranum; i++) {
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'W')
      inputweights();
    else {
      printf("ERROR: INCORRECT AUXILIARY OPTIONS LINE WHICH STARTS WITH %c\n",
	     ch);
      exit(-1);}
  }
  if (weights)
    printweights();
}  /* inputoptions */

Local Void inputdata()
{
  /* input the names and sequences for each species */
  short i, j, k, l, basesread, basesnew;
  Char charstate;
  boolean allread, done;

  putc('\n', outfile);
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
  basesread = 0;
  allread = false;
  while (!(allread)) {
    allread = true;
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    i = 1;
    while (i <= spp ) {
      if ((interleaved && basesread == 0) || !interleaved) {
	for (j = 0; j < nmlngth; j++) {
	  name[i - 1][j] = getc(infile);
	  if (eof(infile) | eoln(infile)){
	    printf("ERROR: END-OF-LINE OR END-OF-FILE");
	    printf(" IN THE MIDDLE OF A SPECIES NAME\n");
	    exit(-1);}
	}
      }
      j = interleaved ? basesread : 0;
      done = false;
      while (!done && !eof(infile)) {
	if (interleaved)
	  done = true;
	while (j < chars && !(eoln(infile) || eof(infile))) {
	  charstate = getc(infile);
	  if (charstate == ' ' || (charstate >= '0' && charstate <= '9'))
	    continue;
	  uppercase(&charstate);
	  if ((strchr("ABCDGHKMNRSTUVWXY?O-.",charstate)) == NULL){
	    printf("ERROR: BAD BASE:%c AT POSITION%5hd OF SPECIES %3hd\n",
		   charstate, j, i);
	    exit(-1);
	  }
	  j++;
	  if (charstate == '.')
	    charstate = y[0][j - 1];
	  y[i - 1][j - 1] = charstate;
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
	basesnew = j;
      fscanf(infile, "%*[^\n]");
      getc(infile);
      if ((interleaved && j != basesnew) || ((!interleaved) && j != chars)){
	printf("ERROR: SEQUENCES OUT OF ALIGNMENT\n");
	exit(-1);}
      i++;
    }
    if (interleaved) {
      basesread = basesnew;
      allread = (basesread == chars);
    } else
      allread = (i > spp);
  }
  if (!printdata)
    return;
   for (i = 1; i <= ((chars - 1) / 60 + 1); i++) {
     for (j = 1; j <= spp; j++) {
      for (k = 0; k < nmlngth; k++)
	putc(name[j - 1][k], outfile);
      fprintf(outfile, "   ");
      l = i * 60;
      if (l > chars)
	l = chars;
      for (k = (i - 1) * 60 + 1; k <= l; k++) {
	if (j > 1 && y[j - 1][k - 1] == y[0][k - 1])
	  charstate = '.';
	else
	  charstate = y[j - 1][k - 1];
	putc(charstate, outfile);
	if (k % 10 == 0 && k % 60 != 0)
	  putc(' ', outfile);
      }
      putc('\n', outfile);
    }
     putc('\n', outfile);
   }
  putc('\n', outfile);
}  /* inputdata */

Local Void sitesort()
{
  /* Shell sort keeping sites, weights in same order */
  short gap, i, j, jj, jg, k, itemp;
  boolean flip, tied;

  gap = chars / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= chars; i++) {
      j = i - gap;
      flip = true;
      while (j > 0 && flip) {
	jj = alias[j - 1];
	jg = alias[j + gap - 1];
	tied = true;
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
	itemp = weight[j - 1];
	weight[j - 1] = weight[j + gap - 1];
	weight[j + gap - 1] = itemp;
	j -= gap;
      }
    }
    gap /= 2;
  }
}  /* sitesort */

Local Void sitecombine()
{
  /* combine sites that have identical patterns */
  short i, j, k;
  boolean tied;

  i = 1;
  while (i < chars) {
    j = i + 1;
    tied = true;
    while (j <= chars && tied) {
      k = 1;
      while (k <= spp && tied) {
	tied = (tied &&
	    y[k - 1][alias[i - 1] - 1] == y[k - 1][alias[j - 1] - 1]);
	k++;
      }
      if (tied) {
	weight[i - 1] += weight[j - 1];
	weight[j - 1] = 0;
	ally[alias[j - 1] - 1] = alias[i - 1];
      }
      j++;
    }
    i = j - 1;
  }
}  /* sitecombine */

Local Void sitescrunch()
{
  /* move so positively weighted sites come first */
  short i, j, itemp;
  boolean done, found;

  done = false;
  i = 1;
  j = 2;
  while (!done) {
    if (ally[alias[i - 1] - 1] != alias[i - 1]) {
      if (j <= i)
	j = i + 1;
      if (j <= chars) {
	found = false;
	do {
	  found = (ally[alias[j - 1] - 1] == alias[j - 1]);
	  j++;
	} while (!(found || j > chars));
	if (found) {
	  j--;
	  itemp = alias[i - 1];
	  alias[i - 1] = alias[j - 1];
	  alias[j - 1] = itemp;
	  itemp = weight[i - 1];
	  weight[i - 1] = weight[j - 1];
	  weight[j - 1] = itemp;
	} else
	  done = true;
      } else
	done = true;
    }
    i++;
    done = (done || i >= chars);
  }
}  /* sitescrunch */

Local Void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  short i;

  for (i = 1; i <= chars; i++) {
    alias[i - 1] = i;
    oldweight[i - 1] = weight[i - 1];
    ally[i - 1] = i;
  }
  sitesort();
  sitecombine();
  sitescrunch();
  endsite = 0;
  for (i = 1; i <= chars; i++) {
    if (ally[i - 1] == i)
      endsite++;
  }
  for (i = 1; i <= endsite; i++)
    location[alias[i - 1] - 1] = i;
  if (!thresh)
    threshold = spp;
  threshwt = (long *)Malloc(endsite*sizeof(long));
  for (i = 0; i < endsite; i++) {
    weight[i] *= 10;
    threshwt[i] = (long)(threshold * weight[i] + 0.5);
  }
}  /* makeweights */

Local Void makevalues()
{
  /* set up fractional likelihoods at tips */
  short i, j;
  char ns;
  node *p;
  for (i = 1; i <= nonodes; i++) {
    treenode[i - 1]->back = NULL;
    treenode[i - 1]->tip = (i <= spp);
    treenode[i - 1]->index = i;
    if (i > spp) {
      p = treenode[i - 1]->next;
      while (p != treenode[i - 1]) {
	p->back = NULL;
	p->tip = false;
	p->index = i;
	p = p->next;
      }
    }
  }
  for (i = 0; i < spp; i++) {
    treenode[i]->numsteps = (stepshortptr)Malloc(endsite*sizeof(short));
    treenode[i]->base = (baseptr)Malloc(endsite*sizeof(char));
  }
  for (i = spp; i < nonodes; i++) {
    p = treenode[i];
    for (j = 1; j <= 3; j++) {
      p->numsteps = (stepshortptr)Malloc(endsite*sizeof(short));
      p->base = (baseptr)Malloc(endsite*sizeof(char));
      p = p->next;
    }
  }
  for (j = 0; j < endsite; j++) {
    for (i = 0; i < spp; i++) {
      switch (y[i][alias[j] - 1]) {

      case 'A':
	ns = 1 << ((char)A);
	break;

      case 'C':
	ns = 1 << ((char)C);
	break;

      case 'G':
	ns = 1 << ((char)G);
	break;

      case 'U':
	ns = 1 << ((char)U);
	break;

      case 'T':
	ns = 1 << ((char)U);
	break;

      case 'M':
	ns = (1 << ((char)A)) | (1 << ((char)C));
	break;

      case 'R':
	ns = (1 << ((char)A)) | (1 << ((char)G));
	break;

      case 'W':
	ns = (1 << ((char)A)) | (1 << ((char)U));
	break;

      case 'S':
	ns = (1 << ((char)C)) | (1 << ((char)G));
	break;

      case 'Y':
	ns = (1 << ((char)C)) | (1 << ((char)U));
	break;

      case 'K':
	ns = (1 << ((char)G)) | (1 << ((char)U));
	break;

      case 'B':
	ns = (1 << ((char)C)) | (1 << ((char)G)) | (1 << ((char)U));
	break;

      case 'D':
	ns = (1 << ((char)A)) | (1 << ((char)G)) | (1 << ((char)U));
	break;

      case 'H':
	ns = (1 << ((char)A)) | (1 << ((char)C)) | (1 << ((char)U));
	break;

      case 'V':
	ns = (1 << ((char)A)) | (1 << ((char)C)) | (1 << ((char)G));
	break;

      case 'N':
	ns = (1 << ((char)A)) | (1 << ((char)C)) | (1 << ((char)G)) |
	     (1 << ((char)U));
	break;

      case 'X':
	ns = (1 << ((char)A)) | (1 << ((char)C)) | (1 << ((char)G)) |
	     (1 << ((char)U));
	break;

      case '?':
	ns = (1 << ((char)A)) | (1 << ((char)C)) | (1 << ((char)G)) |
	     (1 << ((char)U)) | (1 << ((char)O));
	break;

      case 'O':
	ns = 1 << ((char)O);
	break;

      case '-':
	ns = 1 << ((char)O);
	break;
      }
      treenode[i]->base[j] = ns;
      treenode[i]->numsteps[j] = 0;
    }
  }
}  /* makevalues */


Static Void doinput()
{
  /* reads the input data */
  inputoptions();
  inputdata();
  makeweights();
  makevalues();
}  /* doinput */




Local Void fillin(p, left, rt)
node *p, *left, *rt;
{
  /* sets up for each node in the tree the base sequence
     at that point and counts the changes.  The program
     spends much of its time in this PROCEDURE */
  short i;
  short ns, rs, ls;

  for (i = 0; i < endsite; i++) {
    if (left != NULL)
      ls = left->base[i];
    if (rt != NULL)
      rs = rt->base[i];
    if (left == NULL) {
      ns = rs;
      p->numsteps[i] = rt->numsteps[i];
    } else if (rt == NULL) {
      ns = ls;
      p->numsteps[i] = left->numsteps[i];
    } else {
      ns = ls & rs;
      p->numsteps[i] = left->numsteps[i] + rt->numsteps[i];
    }
    if (ns == 0) {
      if (left != NULL && rt != NULL)
	ns = ls | rs;
      p->numsteps[i] += weight[i];
    }
    p->base[i] = ns;
  }
}  /* fillin */

Local Void preorder(p)
node *p;
{
  /* recompute number of steps in preorder taking both ancestoral and
     descendent steps into account */
  if (p == NULL || p->tip)
    return;
  fillin(p->next, p->next->next->back, p->back);
  fillin(p->next->next, p->back, p->next->back);
  preorder(p->next->back);
  preorder(p->next->next->back);
}  /* preorder */

Local Void add(below, newtip, newfork)
node **below, **newtip, **newfork;

{
  /* inserts the nodes newfork and its left descendant, newtip,
     to the tree.  below becomes newfork's right descendant.
     The global variable root is also updated */
  if (*below != treenode[(*below)->index - 1])
    *below = treenode[(*below)->index - 1];
  if ((*below)->back != NULL)
    (*below)->back->back = *newfork;
  (*newfork)->back = (*below)->back;
  (*below)->back = (*newfork)->next->next;
  (*newfork)->next->next->back = *below;
  (*newfork)->next->back = *newtip;
  (*newtip)->back = (*newfork)->next;
  if (root == *below)
    root = *newfork;
  if (!recompute)
    return;
  fillin(*newfork, (*newfork)->next->back, (*newfork)->next->next->back);
  preorder(*newfork);
  if (*newfork != root)
    preorder((*newfork)->back);
}  /* add */

Local Void re_move(item, fork)
node **item, **fork;
{
  /* removes nodes item and its ancestor, fork, from the tree.
     the new descendant of fork's ancestor is made to be
     fork's second descendant (other than item).  Also
     returns pointers to the deleted nodes, item and fork.
     The global variable root is also updated */
  node *p, *q, *other;

  if ((*item)->back == NULL) {
    *fork = NULL;
    return;
  }
  *fork = treenode[(*item)->back->index - 1];
  if (*item == (*fork)->next->back)
    other = (*fork)->next->next->back;
  else
    other = (*fork)->next->back;
  if (root == *fork)
    root = other;
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
  if (!recompute)
    return;
  preorder(other);
  if (other != root)
    preorder(other->back);
}  /* re_move */


Local Void supplement(r)
node *r;
{
  /* determine minimum number of steps more which will
     be added when rest of species are put in tree */
  short i, j, k, sum, sumall, sumadded;
  boolean doneadded, allhave, addedhave, has;
  char supps;

  for (i = 0; i < endsite; i++) {
    sum = 3;
    j = 1;
    doneadded = false;
    do {
      allhave = true;
      addedhave = true;
      supps = suppset[j-1];
      for (k = 0; k < spp; k++) {
	has = ((treenode[k]->base[i] & supps) != 0);
	if (added[k] && !doneadded)
	  addedhave = (addedhave && has);
	allhave = (allhave && has);
      }
      if (allhave)
	sumall = suppno[j - 1];
      if (addedhave)
	sumadded = suppno[j - 1];
      doneadded = (doneadded || addedhave);
      j++;
    } while (!(j > 31 || (allhave && doneadded)));
    if (addedhave && allhave)
      sum = sumall - sumadded;
    r->numsteps[i] += sum * weight[i];
  }
}  /* supplement */

Local Void evaluate(r)
node *r;
{
  /* determines the number of steps needed for a tree. this is
     the minimum number of steps needed to evolve sequences on
     this tree */
  short i, steps;
  double sum;

  sum = 0.0;
  supplement(r);
  for (i = 0; i < endsite; i++) {
    steps = r->numsteps[i];
    if ((long)steps <= threshwt[i])
      sum += steps;
    else
      sum += threshwt[i];
  }
  if (examined == 0 && mults == 0)
    bestyet = -1.0;
  like = sum;
}  /* evaluate */

Local Void postorder(p)
node *p;
{
  /* traverses a binary tree, calling PROCEDURE fillin at a
     node's descendants before calling fillin at the node */
  if (p->tip)
    return;
  postorder(p->next->back);
  postorder(p->next->next->back);
  fillin(p, p->next->back, p->next->next->back);
}  /* postorder */



Local Void addtraverse(p, item, fork, m,n,valyew,place)
node *p, *item, *fork;
short *m,*n;
valptr valyew;
placeptr place;
{
  /* traverse all places to add item */
  if (done)
    return;
  if (*m <= 2 || (p != root && p != root->next->back)) {
    if (p == root)
      fillin(temp, item, p);
    else {
      fillin(temp1, item, p);
      fillin(temp, temp1, p->back);
    }
    (*n)++;
    evaluate(temp);
    examined++;
    if (examined == howoften) {
      examined = 0;
      mults++;
      if (mults == howmany)
	done = true;
      if (progress) {
	printf("%7hd", mults);
	if (bestyet >= 0)
	  printf("%16.1f", bestyet / 10.0);
	else
	  printf("            -   ");
	printf("%17hd%20.2f\n", nextree - 1, fracdone * 100);
      }
    }
    valyew[(*n) - 1] = like;
    place[(*n) - 1] = p->index;
  }
  if (!p->tip) {
    addtraverse(p->next->back, item, fork, m,n,valyew,place);
    addtraverse(p->next->next->back, item, fork,m,n,valyew,place);
  }
}  /* addtraverse */

Local Void shellsort(a, b, n)
double *a;
short *b;
short n;
{
  /* Shell sort keeping a, b in same order */
  short gap, i, j, itemp;
  double rtemp;

  gap = n / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= n; i++) {
      j = i - gap;
      while (j > 0) {
	if (a[j - 1] > a[j + gap - 1]) {
	  rtemp = a[j - 1];
	  a[j - 1] = a[j + gap - 1];
	  a[j + gap - 1] = rtemp;
	  itemp = b[j - 1];
	  b[j - 1] = b[j + gap - 1];
	  b[j + gap - 1] = itemp;
	}
	j -= gap;
      }
    }
    gap /= 2;
  }
}  /* shellsort */

Local Void addit(m)
short m;
{
  /* adds the species one by one, recursively */
  short n;
  valptr valyew;
  placeptr place;
  short i, j, n1, besttoadd;
  valptr bestval;
  placeptr bestplace;
  double oldfrac, oldfdone, sum, bestsum;

  valyew = (valptr)Malloc(nonodes*sizeof(double));
  bestval = (valptr)Malloc(nonodes*sizeof(double));
  place = (placeptr)Malloc(nonodes*sizeof(short));
  bestplace = (placeptr)Malloc(nonodes*sizeof(short));
  if (simple && !firsttime) {
    n = 0;
    added[order[m - 1] - 1] = true;
    addtraverse(root, treenode[order[m - 1] - 1],
		treenode[spp + m - 2], &m,&n,valyew,place);
    besttoadd = order[m - 1];
    memcpy(bestplace, place, nonodes*sizeof(short));
    memcpy(bestval, valyew, nonodes*sizeof(double));
  } else {
    bestsum = -1.0;
    for (i = 1; i <= spp; i++) {
      if (!added[i - 1]) {
	n = 0;
	added[i - 1] = true;
	addtraverse(root, treenode[i - 1], treenode[spp + m - 2],
		    &m,&n,valyew,place);
	added[i - 1] = false;
	sum = 0.0;
	for (j = 0; j < n; j++)
	  sum += valyew[j];
	if (sum > bestsum) {
	  bestsum = sum;
	  besttoadd = i;
	  memcpy(bestplace, place, nonodes*sizeof(short));
	  memcpy(bestval, valyew, nonodes*sizeof(double));
	}
      }
    }
  }
  order[m - 1] = besttoadd;
  memcpy(place, bestplace, nonodes*sizeof(short));
  memcpy(valyew, bestval, nonodes*sizeof(double));
  shellsort(valyew, place, n);
  oldfrac = fracinc;
  oldfdone = fracdone;
  n1 = 0;
  for (i = 0; i < n; i++) {
    if (valyew[i] <= bestyet || bestyet < 0.0)
      n1++;
  }
  if (n1 > 0)
    fracinc /= n1;
  for (i = 0; i < n; i++) {
    if (valyew[i] <= bestyet || bestyet < 0.0) {
      current[m - 1] = place[i];
      recompute = (m < spp);
      add(&treenode[place[i] - 1], &treenode[besttoadd - 1],
	  &treenode[spp + m - 2]);
      added[besttoadd - 1] = true;
      if (m < spp)
	addit(m + 1);
      else {
	if (valyew[i] < bestyet || bestyet < 0.0) {
	  nextree = 1;
	  bestyet = valyew[i];
	}
	if (nextree <= maxtrees) {
	  memcpy(bestorders[nextree - 1], order,
		 spp*sizeof(short));
	  memcpy(bestrees[nextree - 1], current,
		 spp*sizeof(short));
	}
	nextree++;
	firsttime = false;
      }
      recompute = (m < spp);
      re_move(&treenode[besttoadd - 1], &treenode[spp + m - 2]);
      added[besttoadd - 1] = false;
    }
    fracdone += fracinc;
  }
  fracinc = oldfrac;
  fracdone = oldfdone;
  free(valyew);
  free(bestval);
  free(place);
  free(bestplace);
}  /* addit */

Local Void reroot(outgroup)
node *outgroup;
{
  /* reorients tree, putting outgroup in desired position. */
  node *p, *q, *newbottom, *oldbottom;

  if (outgroup->back->index == root->index)
    return;
  newbottom = outgroup->back;
  p = treenode[newbottom->index - 1]->back;
  while (p->index != root->index) {
    oldbottom = treenode[p->index - 1];
    treenode[p->index - 1] = p;
    p = oldbottom->back;
  }
  p = root->next;
  q = root->next->next;
  p->back->back = q->back;
  q->back->back = p->back;
  p->back = outgroup;
  q->back = outgroup->back;
  outgroup->back->back = root->next->next;
  outgroup->back = root->next;
  treenode[newbottom->index - 1] = newbottom;
}  /* reroot */


Local Void coordinates(p, tipy)
node *p;
short *tipy;
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
    coordinates(q->back, tipy);
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

Local Void drawline(i, scale)
short i;
double scale;
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q, *r, *first, *last;
  short n, j;
  boolean extra, done;

  p = root;
  q = root;
  extra = false;
  if (i == p->ycoord && p == root) {
    if (p->index - spp >= 10)
      fprintf(outfile, "-%2hd", p->index - spp);
    else
      fprintf(outfile, "--%hd", p->index - spp);
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
    n = (short)(scale * (p->xcoord - q->xcoord) + 0.5);
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
	  fprintf(outfile, "%2hd", q->index - spp);
	else
	  fprintf(outfile, "-%hd", q->index - spp);
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
      putc(name[p->index - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */

Local Void printree()
{
  /* prints out diagram of the tree */
/* Local variables for printree: */
  short tipy;
  double scale;
  short i;

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



Local Void ancestset(a, b, c)
char *a, *b, *c;
{
  /* make the set of ancestral states below nodes
     whose base sets are a and b */
  *c = (*a) & (*b);
  if (*c == 0)
    *c = (*a) | (*b);
}  /* ancestset */


Local Void hyprint(b1, b2, bottom,HyptravVars)
short b1, b2;
boolean *bottom;
struct LOC_hyptrav *HyptravVars;
{
  /* print out states in sites b1 through b2 at node */
  short i, j, k, n;
  boolean dot;
  bases b;

  if (*bottom) {
    if (!outgropt)
      fprintf(outfile, "       ");
    else
      fprintf(outfile, "root   ");
  } else
    fprintf(outfile, "%4hd   ", HyptravVars->r->back->index - spp);
  if (HyptravVars->r->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(name[HyptravVars->r->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4hd      ", HyptravVars->r->index - spp);
  if (*bottom)
    fprintf(outfile, "          ");
  else if (HyptravVars->nonzero)
    fprintf(outfile, "   yes    ");
  else if (HyptravVars->maybe)
    fprintf(outfile, "  maybe   ");
  else
    fprintf(outfile, "   no     ");
  for (i = b1; i <= b2; i++) {
    j = location[ally[i - 1] - 1];
    HyptravVars->tempset = HyptravVars->r->base[j - 1];
    HyptravVars->anc = HyptravVars->hypset[j - 1];
    if (!(*bottom))
      HyptravVars->anc = treenode[HyptravVars->r->back->index - 1]->base[j - 1];
    dot = (HyptravVars->tempset == HyptravVars->anc && !(*bottom));
    if (dot)
      putc('.', outfile);
    else if (HyptravVars->tempset == 1 << ((char)A))
      putc('A', outfile);
    else if (HyptravVars->tempset == 1 << ((char)C))
      putc('C', outfile);
    else if (HyptravVars->tempset == 1 << ((char)G))
      putc('G', outfile);
    else if (HyptravVars->tempset == 1 << ((char)U))
      putc('T', outfile);
    else if (HyptravVars->tempset == 1 << ((char)O))
      putc('-', outfile);
    else {
      k = 1;
      n = 0;
      for (b = A; (char)b <= (char)O; b = (bases)((char)b + 1)) {
	if (((1 << ((char)b)) & HyptravVars->tempset) != 0)

	  n += k;
	k += k;
      }
      putc(basechar[n - 1], outfile);
    }
    if (i % 10 == 0)
      putc(' ', outfile);
  }
  putc('\n', outfile);
}  /* hyprint */


Local Void hyptrav(r_, hypset_, b1, b2, bottom)
node *r_;
char *hypset_;
short b1, b2;
boolean *bottom;
{
  /*  compute, print out states at one interior node */
  struct LOC_hyptrav Vars;
  short i, j;
  char left, rt;
  gbases *temparray, *ancset;

  Vars.r = r_;
  Vars.hypset = hypset_;
  gnu(&ancset);
  gnu(&temparray);
  Vars.maybe = false;
  Vars.nonzero = false;
  for (i = b1 - 1; i < b2; i++) {
    j = location[ally[i] - 1];
    Vars.anc = Vars.hypset[j - 1];
    if (!Vars.r->tip) {
      left = Vars.r->next->back->base[j - 1];
      rt = Vars.r->next->next->back->base[j - 1];
      Vars.tempset = left & rt & Vars.anc;
      if (Vars.tempset == 0) {
	Vars.tempset = (left & rt) | (left & Vars.anc) | (rt & Vars.anc);
	if (Vars.tempset == 0)
	  Vars.tempset = left | rt | Vars.anc;
      }
      Vars.r->base[j - 1] = Vars.tempset;
    }
    if (!(*bottom))
      Vars.anc = treenode[Vars.r->back->index - 1]->base[j - 1];
    Vars.nonzero = (Vars.nonzero || (Vars.r->base[j - 1] & Vars.anc) == 0);
    Vars.maybe = (Vars.maybe || Vars.r->base[j - 1] != Vars.anc);
  }
  hyprint(b1, b2, bottom,&Vars);
  (*bottom) = false;
  if (!Vars.r->tip) {
    memcpy(temparray->base, Vars.r->next->back->base, endsite*sizeof(char));
    for (i = b1 - 1; i < b2; i++) {
      j = location[ally[i] - 1];
      ancestset(&Vars.hypset[j - 1], &Vars.r->next->next->back->base[j - 1],
		&ancset->base[j - 1]);
    }
    hyptrav(Vars.r->next->back, ancset->base, b1, b2,bottom);
    for (i = b1 - 1; i < b2; i++) {
      j = location[ally[i] - 1];
      ancestset(&Vars.hypset[j - 1], &temparray->base[j - 1],
		&ancset->base[j - 1]);
    }
    hyptrav(Vars.r->next->next->back, ancset->base, b1, b2,bottom);
  }
  chuck(temparray);
  chuck(ancset);
}  /* hyptraVars */

Local Void hypstates()
{
  /* fill in and describe states at interior nodes */
  boolean bottom;
  short i, n;

  fprintf(outfile, "\nFrom    To     Any Steps?    State at upper node\n");
  fprintf(outfile, "                            ");
  fprintf(outfile, " ( . means same as in the node below it on tree)\n\n");
  nothing = (baseptr)Malloc(endsite*sizeof(char));
  for (i = 0; i < endsite; i++)
    nothing[i] = 0;
  for (i = 1; i <= ((chars - 1 ) / 40 + 1); i++) {
    putc('\n', outfile);
    bottom = true;
    n = i * 40;
    if (n > chars)
      n = chars;
    hyptrav(root, nothing, i * 40 - 39, n, &bottom);
  }
  free(nothing);
}  /* hypstates */

Local Void treeout(p)
node *p;
{
  /* write out file with representation of final tree */
  short i, n;
  Char c;

  if (p->tip) {
    n = 0;
    for (i = 1; i <= nmlngth; i++) {
      if (name[p->index - 1][i - 1] != ' ')
	n = i;
    }
    for (i = 0; i < n; i++) {
      c = name[p->index - 1][i];
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

Local Void describe()
{
  /* prints ancestors, steps and table of numbers of steps in
     each site */
  short i, j, k, l;

  if (stepbox) {
    putc('\n', outfile);
    if (weights)
      fprintf(outfile, " weighted");
    fprintf(outfile, " steps in each site:\n");
    fprintf(outfile, "      ");
    for (i = 0; i <= 9; i++)
      fprintf(outfile, "%4hd", i);
    fprintf(outfile, "\n     *------------------------------------");
    fprintf(outfile, "-----\n");
    for (i = 0; i <= (chars / 10); i++) {
      fprintf(outfile, "%5hd", i * 10);
      putc('!', outfile);
      for (j = 0; j <= 9; j++) {
	k = i * 10 + j;
	if (k == 0 || k > chars)
	  fprintf(outfile, "    ");
	else {
	  l = location[ally[k - 1] - 1];
	  if (oldweight[k - 1] > 0)
	    fprintf(outfile, "%4hd",
		    oldweight[k - 1] *
		    (root->numsteps[l - 1] / weight[l - 1]));
	  else
	    fprintf(outfile, "   0");
	}
      }
      putc('\n', outfile);
    }
  }
  if (ancseq) {
    hypstates();
    putc('\n', outfile);
  }
  putc('\n', outfile);
  if (trout) {
    col = 0;
    treeout(root);
  }
}  /* describe */


Static Void maketree()
{
  /* tree construction recursively by branch and bound */
  short i, j, k;
  node *dummy;

  if (progress) {
    printf("\nHow many\n");
    printf("trees looked                                       Approximate\n");
    printf("at so far      Length of        How many           percentage\n");
    printf("(multiples     shortest tree    trees this long    searched\n");
    printf("of %4hd):       found so far     found so far       so far\n",
	   howoften);
    printf("----------     ------------     ------------       ------------\n");
  }
  done = false;
  mults = 0;
  examined = 0;
  nextree = 1;
  root = treenode[0];
  firsttime = true;
  for (i = 0; i < spp; i++)
    added[i] = false;
  added[0] = true;
  order[0] = 1;
  k = 2;
  fracdone = 0.0;
  fracinc = 1.0;
  bestyet = -1.0;
  recompute = true;
  temp = (node *)Malloc(sizeof(node));
  temp->numsteps = (stepshortptr)Malloc(endsite*sizeof(short));
  temp->base = (baseptr)Malloc(endsite*sizeof(char));
  temp1 = (node *)Malloc(sizeof(node));
  temp1->numsteps = (stepshortptr)Malloc(endsite*sizeof(short));
  temp1->base = (baseptr)Malloc(endsite*sizeof(char));
  addit(k);
  if (done) {
    if (progress) {
      printf("Search broken off!  Not guaranteed to\n");
      printf(" have found the most parsimonious trees.\n");
    }
    if (treeprint) {
      fprintf(outfile, "Search broken off!  Not guaranteed to\n");
      fprintf(outfile, " have found the most parsimonious\n");
      fprintf(outfile, " trees, but here is what we found:\n");
    }
  }
  if (treeprint) {
    fprintf(outfile, "\nrequires a total of %18.3f\n\n", bestyet / 10.0);
    if (nextree == 2)
      fprintf(outfile, "One most parsimonious tree found:\n");
    else
      fprintf(outfile, "%6hd trees in all found\n", nextree - 1);
  }
  if (nextree > maxtrees + 1) {
    if (treeprint)
      fprintf(outfile, "here are the first%4hd of them\n", (short)maxtrees);
    nextree = maxtrees + 1;
  }
  if (treeprint)
    putc('\n', outfile);
  for (i = 0; i < spp; i++)
    added[i] = true;
  for (i = 0; i <= nextree - 2; i++) {
    root = treenode[0];
    for (j = k; j <= spp; j++)
      add(&treenode[bestrees[i][j - 1] - 1],
	  &treenode[bestorders[i][j - 1] - 1], &treenode[spp + j - 2]);
    reroot(treenode[outgrno - 1]);
    postorder(root);
    evaluate(root);
    printree();
    describe();
    for (j = k - 1; j < spp; j++)
      re_move(&treenode[bestorders[i][j] - 1], &dummy);
  }
  if (progress) {
    printf("\nOutput written to output file\n\n");
    if (trout)
      printf("Trees also written onto treefile\n\n");
  }
  free(temp->numsteps);
  free(temp->base);
  free(temp);
  free(temp1->numsteps);
  free(temp1->base);
  free(temp1);
}  /* maketree */


main(argc, argv)
int argc;
Char *argv[];
{  /* Penny's branch-and-bound method for DNA sequences */
char infilename[100],outfilename[100],trfilename[100];
#ifdef MAC
  macsetup("Dnapenny","");
  argv[0] = "Dnapenny";
#endif

  /* Reads in the number of species, number of characters,
     options and data.  Then finds all most parsimonious trees */
  openfile(&infile,INFILE,"r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);

  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  mulsets = false;
  garbage = NULL;
  datasets = 1;
  firstset = true;
  doinit();
  if (trout)
    openfile(&treefile,TREEFILE,"w",argv[0],trfilename);
  weight = (short *)Malloc(chars*sizeof(short));
  oldweight = (short *)Malloc(chars*sizeof(short));
  alias = (steptr)Malloc(chars*sizeof(short));
  ally = (steptr)Malloc(chars*sizeof(short));
  location = (steptr)Malloc(chars*sizeof(short));
  name = (Char **)Malloc(spp*sizeof(Char *));
  for (j = 1; j <= spp; j++)
    name[j - 1] = (Char *)Malloc(nmlngth*sizeof(char));
  bestorders = (treenumbers *)Malloc(maxtrees*sizeof(treenumbers));
  for (j = 1; j <= maxtrees; j++)
    bestorders[j - 1] = (treenumbers)Malloc(spp*sizeof(short));
  bestrees = (treenumbers *)Malloc(maxtrees*sizeof(treenumbers));
  for (j = 1; j <= maxtrees; j++)
    bestrees[j - 1] = (treenumbers)Malloc(spp*sizeof(short));
  current = (treenumbers)Malloc(spp*sizeof(short));
  order = (treenumbers)Malloc(spp*sizeof(short));
  added = (boolean *)Malloc(nonodes*sizeof(boolean));
  for (ith = 1; ith <= datasets; ith++) {
    doinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1) {
      fprintf(outfile, "\nData set # %hd:\n",ith);
      if (progress)
        printf("\nData set # %hd:\n",ith);
    }
    maketree();
    free(threshwt);
    for (i = 0; i < spp; i++) {
      free(treenode[i]->numsteps);
      free(treenode[i]->base);
    }
    for (i = spp; i < nonodes; i++) {
      p = treenode[i];
      for (j = 1; j <= 3; j++) {
	free(p->numsteps);
	free(p->base);
	p = p->next;
      }
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
}  /* Penny's branch-and-bound method for DNA sequences */

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
mem = (MALLOCRETURN *)malloc((size_t)x);
if (!mem)
  memerror();
else
  return (MALLOCRETURN *)mem;
}
