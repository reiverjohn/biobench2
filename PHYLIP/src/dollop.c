#include "phylip.h"

/* version 3.572c. (c) Copyright 1995 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define nmlngth         10   /* number of characters in species name    */
#define maxtrees        100  /* maximum number of tied trees stored     */
#define maxuser         10   /* maximum number of user-defined trees    */

#define ibmpc0          false
#define ansi0           true
#define vt520           false
#define down            2

typedef long *bitptr;

/* nodes will form a binary tree */

typedef struct node {
                              /* describes a tip species or an ancestor */
  struct node *next, *back;   /* pointers to nodes                      */
  long index;                /* number of the node                     */
  boolean tip;                /* present species are tips of tree       */
  bitptr stateone, statezero; /* see in PROCEDURE fillin                */
  long xcoord, ycoord, ymin; /* used by printree                       */
  long ymax;
} node;

typedef long *steptr;
typedef node **pointptr;
typedef long longer[6];

typedef struct gbit {
  bitptr bits_;
  struct gbit *next;
} gbit;

Static node *root;
Static FILE *infile, *outfile, *treefile;
Static long spp, nonodes, chars, wrds, words, inseed, col, datasets, ith, j,
  l, jumb, njumble;
/* spp = number of species
   nonodes = number of nodes in tree
   chars = number of binary characters
   words = number of words needed for characters of one organism */
Static boolean jumble, usertree, weights, thresh, ancvar, questions, dollo,
               trout,  printdata, progress, treeprint, stepbox,
               ancseq, mulsets, firstset, ibmpc, vt52, ansi;
Static steptr extras, weight;
Static boolean *ancone, *anczero, *ancone0, *anczero0;
Static pointptr treenode;   /* pointers to all nodes in tree */
Static Char **nayme;   /* names of species */
Static double threshold;
Static double *threshwt;
Static longer seed;
Static long *enterorder;
Static double **fsteps;
Static steptr numsteps;
Static long **bestrees;
Static Char *guess;
Static gbit *garbage;
Static long bits = 8*sizeof(long) - 1;

/* Local variables for maketree, propagated globally for C version: */
long nextree, which, minwhich;
double like, bestyet, bestlike, bstlike2, minsteps;
boolean lastrearr;
double nsteps[maxuser];
node *there;
long fullset;
bitptr zeroanc, oneanc;
long *place;
Char ch;
boolean *naymes;
steptr numsone, numszero;
bitptr steps;


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
gbit **p;
{
  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */
  if (garbage != NULL) {
    *p = garbage;
    garbage = garbage->next;
  } else {
    *p = (gbit *)Malloc(sizeof(gbit));
    (*p)->bits_ = (bitptr)Malloc(words*sizeof(long));
  }
  (*p)->next = NULL;
}  /* gnu */


void chuck(p)
gbit *p;
{
  /* collect garbage on p -- put it on front of garbage list */
  p->next = garbage;
  garbage = p;
}  /* chuck */

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
{  /* convert a character to upper case -- either ASCII or EBCDIC */
  *ch = (islower (*ch) ? toupper(*ch) : (*ch));
}  /* uppercase */

void newline(i, j, k)
long i, j, k;
{
  /* go to new line if i is a multiple of j, indent k spaces */
  long m;

  if ((i - 1) % j != 0 || i <= 1)
    return;
  putchar('\n');
  for (m = 1; m <= k; m++)
    putchar(' ');
}  /* newline */


void getoptions()
{
  /* interactively set options */
  long i, inseed0;
  Char ch;
  boolean  done1;

  fprintf(outfile,"\nDollo and polymorphism parsimony algorithm,");
  fprintf(outfile," version %s\n\n",VERSION);
  putchar('\n');
  ancvar = false;
  dollo = true;
  jumble = false;
  njumble = 1;
  thresh = false;
  threshold = spp;
  trout = true;
  usertree = false;
  weights = false;
  printdata = false;
  progress = true;
  treeprint = true;
  stepbox = false;
  ancseq = false;
  for (;;) {
    printf(ansi ? "\033[2J\033[H" :
	   vt52 ? "\033E\033H"    : "\n");
    printf("\nDollo and polymorphism parsimony algorithm, version %s\n\n",
           VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
	   (usertree ? "No, use user trees in input file" : "Yes"));
    printf("  P                     Parsimony method?  %s\n",
	   dollo ? "Dollo" : "Polymorphism");
    if (!usertree) {
      printf("  J     Randomize input order of species?");
      if (jumble)
        printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  T              Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f per char.\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  A   Use ancestral states in input file?  %s\n",
	   ancvar ? "Yes" : "No");
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  0   Terminal type (IBM PC, VT52, ANSI)?  %s\n",
	   ibmpc ? "IBM PC" :
	   ansi  ? "ANSI"   :
	   vt52  ? "VT52"   : "(none)");

    printf("  1    Print out the data at start of run  %s\n",
	   printdata ? "Yes" : "No");
    printf("  2  Print indications of progress of run  %s\n",
	   progress ? "Yes" : "No");
    printf("  3                        Print out tree  %s\n",
	   treeprint ? "Yes" : "No");
    printf("  4     Print out steps in each character  %s\n",
	   stepbox ? "Yes" : "No");
    printf("  5     Print states at all nodes of tree  %s\n",
	   ancseq ? "Yes" : "No");
    printf("  6       Write out trees onto tree file?  %s\n",
	   trout ? "Yes" : "No");
    printf("\nAre these settings correct? ");
    printf("(type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (strchr("APJTUM1234560",ch) != NULL){
      switch (ch) {
	
      case 'A':
	ancvar = !ancvar;
	break;
	
      case 'P':
	dollo = !dollo;
	break;
	
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
	  printf("Number of times to jumble?\n");
	  scanf("%ld%*[^\n]", &njumble);
	  getchar();
	}
	else njumble = 1;
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

void inputnumbers()
{
  /* input the numbers of species and of characters */
  fscanf(infile, "%ld%ld", &spp, &chars);
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  characters\n", spp, chars);
  if (printdata)
    putc('\n', outfile);
  words = chars / bits + 1;
  nonodes = spp * 2 - 1;
}  /* inputnumbers */


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
    treenode[i]->stateone = (bitptr)Malloc(words*sizeof(long));
    treenode[i]->statezero = (bitptr)Malloc(words*sizeof(long));
  }
  for (i = spp; i < (nonodes); i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc(sizeof(node));
      p->stateone = (bitptr)Malloc(words*sizeof(long));
      p->statezero = (bitptr)Malloc(words*sizeof(long));
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    treenode[i] = p;
  }
}  /* doinit */

void inputweights()
{
  /* input the character weights, 0-9 and A-Z for weights 0 - 35 */
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
        weight[i] = ch - 'A' + 10;
      else if (ch >= 'J' && ch <= 'R')
        weight[i] = ch - 'J' + 19;
      else
        weight[i] = ch - 'S' + 28;
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
  /* print out the weights of characters */
  long i, j, k;

  fprintf(outfile, "Characters are weighted as follows:\n");
  fprintf(outfile, "       ");
  for (i = 0; i <= 9; i++)
    fprintf(outfile, "%3ld", i);
  fprintf(outfile, "\n*---------------------------------\n");
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

void inputancestors()
{
  /* reads the ancestral states for each character */
  long i;
  Char ch;

  for (i = 1; i < nmlngth; i++)
    ch = getc(infile);
  for (i = 0; i < (chars); i++) {
    anczero0[i] = true;
    ancone0[i] = true;
    do {
      if (eoln(infile)) {
        fscanf(infile, "%*[^\n]");
        getc(infile);
      }
      ch = getc(infile);
    } while (ch == ' ');
    if (ch == 'p')
      ch = 'P';
    if (ch == 'b')
      ch = 'B';
    if (ch == '0' || ch == '1' || ch == 'P' || ch == 'B' || ch == '?') {
      switch (ch) {

      case '1':
        anczero0[i] = false;
        break;

      case '0':
        ancone0[i] = false;
        break;

      case 'P':
        /* blank case */
        break;

      case 'B':
        /* blank case */
        break;

      case '?':
        /* blank case */
        break;
      }
    } else {
      printf("BAD ANCESTOR STATE: %c AT CHARACTER %4ld\n", ch, i + 1);
      exit(-1);
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
}  /* inputancestors */

void printancestors()
{
  /* print out list of ancestral states */
  long i;

  fprintf(outfile, "Ancestral states:\n");
  for (i = 1; i <= nmlngth + 3; i++)
    putc(' ', outfile);
  for (i = 1; i <= (chars); i++) {
    newline(i, 55, (int)(nmlngth + 3));
    if (ancone[i - 1] && anczero[i - 1])
      putc('?', outfile);
    else if (ancone[i - 1])
      putc('1', outfile);
    else
      putc('0', outfile);
    if (i % 5 == 0)
      putc(' ', outfile);
  }
  fprintf(outfile, "\n\n");
}  /* printancestor */

void inputoptions()
{
  /* input the information on the options */
  Char ch;
  long extranum, i, cursp, curchs;
  boolean avar;

  if (!firstset) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%ld%ld", &cursp, &curchs);
    if (cursp != spp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4ld\n", ith);
      exit(-1);
    }
    chars = curchs;
  }
  extranum = 0;
  avar = false;
  while (!(eoln(infile))) {
    ch = getc(infile);
    uppercase(&ch);
    if (ch == 'A' || ch == 'W')
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
    if (ch != 'A' && ch != 'W') {
      printf("ERROR: INCORRECT AUXILIARY OPTIONS LINE");
      printf(" WHICH STARTS WITH %c\n", ch);
      exit(-1);
    }
    if (ch == 'A') {
      avar = true;
      if (!ancvar) {
	printf("ERROR: ANCESTOR OPTION NOT CHOSEN IN MENU");
	printf(" WITH OPTION %c IN INPUT\n", ch);
	exit(-1);
      } else
	inputancestors();
    }
    if (ch == 'W')
      inputweights();
  }
  if (ancvar && !avar) {
    puts("ERROR: ANCESTOR OPTION CHOSEN IN MENU WITH NO OPTION A IN INPUT");
    exit(-1);
  }
  if (dollo)
    fprintf(outfile, "Dollo");
  else
    fprintf(outfile, "Polymorphism");
  fprintf(outfile, " parsimony method\n\n");
  if (weights && printdata)
    printweights();
  for (i = 0; i < (chars); i++) {
    if (!ancvar) {
      anczero[i] = true;
      ancone[i] = false;
    } else {
      anczero[i] = anczero0[i];
      ancone[i] = ancone0[i];
    }
  }
  if (ancvar && avar && printdata)
    printancestors();
  questions = false;
  for (i = 0; i < (chars); i++) {
    questions = (questions || (ancone[i] && anczero[i]));
    threshwt[i] = threshold * weight[i];
  }
}  /* inputoptions */

void inputdata()
{
  /* input the names and character state data for species */
  long i, j, l;
  char k;
  Char charstate;
  /* possible states are
     '0', '1', 'P', 'B', and '?' */
  node *p;

  putc('\n', outfile);
  j = nmlngth + (chars + (chars - 1) / 10) / 2 - 4;
  if (j < nmlngth - 1)
    j = nmlngth - 1;
  if (j > 36)
    j = 36;
  if (printdata) {
    fprintf(outfile, "Name");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "Characters\n");
    fprintf(outfile, "----");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "----------\n\n");
  }
  for (i = 0; i < (chars); i++)
    extras[i] = 0;
  for (i = 1; i <= (nonodes); i++) {
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
    } else {
      for (j = 0; j < nmlngth; j++) {
	nayme[i - 1][j] = getc(infile);
	if (eof(infile) || eoln(infile)){
	  printf("ERROR: END-OF-LINE OR END-OF-FILE");
	  printf(" IN THE MIDDLE OF A SPECIES NAME\n");
	  exit(-1);}
      }
      if (printdata) {
	for (j = 0; j < nmlngth; j++)
	  putc(nayme[i - 1][j], outfile);
      }
      fprintf(outfile, "   ");
      for (j = 0; j < (words); j++) {
	treenode[i - 1]->stateone[j] = 0;
	treenode[i - 1]->statezero[j] = 0;
      }
      for (j = 1; j <= (chars); j++) {
	k = (j - 1) % bits + 1;
	l = (j - 1) / bits + 1;
	do {
	  if (eoln(infile)) {
	    fscanf(infile, "%*[^\n]");
	    getc(infile);
	  }
	  charstate = getc(infile);
	} while (charstate == ' ');
	if (charstate == 'b')
	  charstate = 'B';
	if (charstate == 'p')
	  charstate = 'P';
	if (charstate != '0' && charstate != '1' && charstate != '?' &&
	    charstate != 'P' && charstate != 'B') {
	  printf("ERROR: BAD CHARACTER STATE: %c ",charstate);
	  printf("AT CHARACTER %5ld OF SPECIES %3ld\n",j,i);
	  exit(-1);
	}
	if (printdata) {
	  newline(j, 55, (int)(nmlngth + 3));
	  putc(charstate, outfile);
	  if (j % 5 == 0)
	    putc(' ', outfile);
	}
	if (charstate == '1')
	  treenode[i - 1]->stateone[l - 1] =
	    ((long)treenode[i - 1]->stateone[l - 1]) | (1L << k);
	if (charstate == '0')
	  treenode[i - 1]->statezero[l - 1] =
	    ((long)treenode[i - 1]->statezero[l - 1]) | (1L << k);
	if (charstate == 'P' || charstate == 'B') {
	  if (dollo)
	    extras[j - 1] += weight[j - 1];
	  else {
	    treenode[i - 1]->stateone[l - 1] =
	      ((long)treenode[i - 1]->stateone[l - 1]) | (1L << k);
	    treenode[i - 1]->statezero[l - 1] =
	      ((long)treenode[i - 1]->statezero[l - 1]) | (1L << k);
	  }
	}
      }
      fscanf(infile, "%*[^\n]");
      getc(infile);
      if (printdata)
	putc('\n', outfile);
    }
  }
  fprintf(outfile, "\n\n");
}  /* inputdata */


void doinput()
{
  /* reads the input data */
  inputoptions();
  inputdata();
}  /* doinput */


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
}  /* re_move */

void fillin(p)
node *p;
{
  /* Sets up for each node in the tree two statesets.
     stateone and statezero are the sets of character
     states that must be 1 or must be 0, respectively,
     in a most parsimonious reconstruction, based on the
     information at or above this node.  Note that this
     state assignment may change based on information further
     down the tree.  If a character is in both sets it is in
     state "P".  If in neither, it is "?". */
  long i;

  for (i = 0; i < words; i++) {
    p->stateone[i] = p->next->back->stateone[i] |
                     p->next->next->back->stateone[i];
    p->statezero[i] = p->next->back->statezero[i] |
                      p->next->next->back->statezero[i];
  }
}  /* fillin */

void correct(p)
node *p;
{  /* get final states for intermediate nodes */
  long i;
  long z0, z1, s0, s1, temp;

  if (p->tip)
    return;
  for (i = 0; i < (words); i++) {
    if (p->back == NULL) {
      s0 = zeroanc[i];
      s1 = fullset & (~zeroanc[i]);
    } else {
      s0 = treenode[p->back->index - 1]->statezero[i];
      s1 = treenode[p->back->index - 1]->stateone[i];
    }
    z0 = (s0 & p->statezero[i]) |
         (p->next->back->statezero[i] & p->next->next->back->statezero[i]);
    z1 = (s1 & p->stateone[i]) |
         (p->next->back->stateone[i] & p->next->next->back->stateone[i]);
    if (dollo) {
      temp = z0 & (~(zeroanc[i] & z1));
      z1 &= ~(fullset & (~zeroanc[i]) & z0);
      z0 = temp;
    }
    temp = fullset & (~z0) & (~z1);
    p->statezero[i] = z0 | (temp & s0 & (~s1));
    p->stateone[i] = z1 | (temp & s1 & (~s0));
  }
}  /* correct */


void postorder(p,numsone,numszero)
node *p;
steptr numsone,numszero;
{
  /* traverses a binary tree, calling PROCEDURE fillin at a
     node's descendants before calling fillin at the node */
  if (p->tip)
    return;
  postorder(p->next->back, numsone,numszero);
  postorder(p->next->next->back, numsone,numszero);
  fillin(p);
}  /* postorder */


void count(p, numsone,numszero)
node *p;
steptr numsone,numszero;
{
  /* counts the number of steps in a fork of the tree.
     The program spends much of its time in this PROCEDURE */
  long i, j, l;

  if (dollo) {
    for (i = 0; i < (words); i++)
      steps[i] = (treenode[p->back->index - 1]->stateone[i] &
                  p->statezero[i] & zeroanc[i]) |
          (treenode[p->back->index - 1]->statezero[i] & p->stateone[i] &
           fullset & (~zeroanc[i]));
  } else {
    for (i = 0; i < (words); i++)
      steps[i] = treenode[p->back->index - 1]->stateone[i] &
                 treenode[p->back->index - 1]->statezero[i] & p->stateone[i] &
                 p->statezero[i];
  }
  j = 1;
  l = 0;
  for (i = 0; i < (chars); i++) {
    l++;
    if (l > bits) {
      l = 1;
      j++;
    }
    if (((1L << l) & steps[j - 1]) != 0) {
      if (((1L << l) & zeroanc[j - 1]) != 0)
        numszero[i] += weight[i];
      else
        numsone[i] += weight[i];
    }
  }
}  /* count */

void preorder(p, numsone,numszero)
node *p;
steptr numsone,numszero;
{
  /* go back up tree setting up and counting interior node
     states */

  if (!p->tip) {
    correct(p);
    preorder(p->next->back, numsone,numszero);
    preorder(p->next->next->back, numsone,numszero);
  }
  if (p->back != NULL)
    count(p, numsone,numszero);
}  /* preorder */

void evaluate(r)
node *r;
{
  /* Determines the number of losses or polymorphisms needed
     for a tree. This is the minimum number needed to evolve
     chars on this tree */
  long i, stepnum, smaller;
  double sum, term;

  sum = 0.0;
  for (i = 0; i < (chars); i++) {
    numszero[i] = 0;
    numsone[i] = 0;
  }
  for (i = 0; i < (words); i++)
    zeroanc[i] = fullset;
  postorder(r, numsone,numszero);
  preorder(r, numsone,numszero);
  for (i = 0; i < (words); i++)
    zeroanc[i] = 0;
  postorder(r, numsone,numszero);
  preorder(r, numsone,numszero);
  for (i = 0; i < (chars); i++) {
    smaller = spp * weight[i];
    numsteps[i] = smaller;
    if (anczero[i]) {
      numsteps[i] = numszero[i];
      smaller = numszero[i];
    }
    if (ancone[i] && numsone[i] < smaller)
      numsteps[i] = numsone[i];
    stepnum = numsteps[i] + extras[i];
    if (stepnum <= threshwt[i])
      term = stepnum;
    else
      term = threshwt[i];
    sum += term;
    if (usertree && which <= maxuser)
      fsteps[which - 1][i] = term;
    guess[i] = '?';
    if (!ancone[i] || (anczero[i] && numszero[i] < numsone[i]))
      guess[i] = '0';
    else if (!anczero[i] || (ancone[i] && numsone[i] < numszero[i]))
      guess[i] = '1';
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


void savetree()
{
  /* record in place where each species has to be
     added to reconstruct this tree */
  long i, j;
  node *p;
  boolean done;

  for (i = 0; i < (nonodes); i++)
    place[i] = 0;
  place[root->index - 1] = 1;
  for (i = 1; i <= (spp); i++) {
    p = treenode[i - 1];
    while (place[p->index - 1] == 0) {
      place[p->index - 1] = i;
      p = p->back;
      if (p != NULL)
        p = treenode[p->index - 1];
    }
    if (i > 1) {
      place[i - 1] = place[p->index - 1];
      j = place[p->index - 1];
      done = false;
      while (!done) {
        place[p->index - 1] = spp + i - 1;
        p = treenode[p->index - 1];
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
  (*found) = false;
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
    if ((*found))
      break;
    if (below)
      upper = (*pos) - 1;
    else
      lower = (*pos) + 1;
  }
  if (!(*found) && !below)
    (*pos)++;
}  /* findtree */

void addtree(pos)
long *pos;
{
  /* puts tree from ARRAY place in its proper position
     in ARRAY bestrees */
  long i;

  for (i = nextree - 1; i >=(*pos); i--)
    memcpy(bestrees[i], bestrees[i - 1], spp*sizeof(long));
  for (i = 0; i < (spp); i++)
    bestrees[(*pos) - 1][i] = place[i];
  nextree++;
}  /* addtree */

void tryadd(p, item,nufork)
node *p;
node **item,**nufork;
{
  /* temporarily adds one fork and one tip to the tree.
     if the location where they are added yields greater
     "likelihood" than other locations tested up to that
     time, then keeps that location as there */
/* Local variables for tryadd: */
  long pos;
  boolean found;
  add(p, *item, *nufork);
  evaluate(root);
  if (lastrearr) {
    if (like >= bstlike2) {
      savetree();
      if (like > bstlike2) {
        bestlike = bstlike2 = like;
        pos = 1;
        nextree = 1;
        addtree(&pos);
      } else {
        pos = 0;
        findtree(&pos,&found);
        if (!found) {
          if (nextree <= maxtrees)
            addtree(&pos);
        }
      }
    }
  }
  if (like > bestyet) {
    bestyet = like;
    there = p;
  }
  re_move(item, nufork);
}  /* tryadd */

void addpreorder(p, item_, nufork_)
node *p, *item_, *nufork_;
{
  /* traverses a binary tree, calling PROCEDURE tryadd
     at a node before calling tryadd at its descendants */
  node *item= item_;
  node *nufork = nufork_;

  if (p == NULL)
    return;
  tryadd(p, &item,&nufork);
  if (!p->tip) {
    addpreorder(p->next->back, item, nufork);
    addpreorder(p->next->next->back, item, nufork);
  }
}  /* addpreorder */


void tryrearr(p, r,success)
node *p;
node **r;
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
  oldlike = bestyet;
  if (p->back->next->next == forknode)
    frombelow = forknode->next->next->back;
  else
    frombelow = forknode->next->back;
  whereto = forknode->back;
  re_move(&p, &forknode);
  add(whereto, p, forknode);
  evaluate(*r);
  if (like <= oldlike) {
    re_move(&p, &forknode);
    add(frombelow, p, forknode);
  } else {
    (*success) = true;
    bestyet = like;
  }
}  /* tryrearr */

void repreorder(p,r,success)
node *p;
node **r;
boolean *success;
{
  /* traverses a binary tree, calling PROCEDURE tryrearr
     at a node before calling tryrearr at its descendants */
  if (p == NULL)
    return;
  tryrearr(p, r,success);
  if (!p->tip) {
    repreorder(p->next->back, r,success);
    repreorder(p->next->next->back, r,success);
  }
}  /* repreorder */

void rearrange(r_)
node **r_;
{
  /* traverses the tree (preorder), finding any local
     rearrangement which decreases the number of steps.
     if traversal succeeds in increasing the tree's
     "likelihood", PROCEDURE rearrange runs traversal again */
  node **r         = r_;
  boolean success  = true;

  while (success) {
    success = false;
    repreorder(*r, r,&success);
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
    if ((done && ch == ')') || !(done)) {
      if (eoln(infile)) {
        fscanf(infile, "%*[^\n]");
        getc(infile);
      }
      ch = getc(infile);
    }
  }
}  /* findch */

void addelement(p, nextnode,lparens,naymes)
node **p;
long *nextnode,*lparens;
boolean *naymes;
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
    addelement(&q->next->back, nextnode,lparens,naymes);
    q->next->back->back = q->next;
    findch(',' );
    addelement(&q->next->next->back, nextnode,lparens,naymes);
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
      if (naymes[n - 1] == false) {
        *p = treenode[n - 1];
        naymes[n - 1] = true;
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
  Char ch;
  long nextnode    = spp;
  long lparens     =0;
  long i;

  root = treenode[spp];
  root->back = NULL;
  for (i = 0; i < (spp); i++)
    naymes[i] = false;
  addelement(&root, &nextnode,&lparens,naymes);
  findch(';');
  if (progress)
    printf("\n\n");
  fscanf(infile, "%*[^\n]");
  getc(infile);
}  /* treeread */


void coordinates(p,tipy)
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
  long tipy        = 1;
  double scale     = 1.5;
  long i;

  putc('\n', outfile);
  if (!treeprint)
    return;
  putc('\n', outfile);
  coordinates(root,&tipy);
  scale = 1.5;
  putc('\n', outfile);
  for (i = 1; i <= (tipy - down); i++)
    drawline(i, scale);
  putc('\n', outfile);
}  /* printree*/



void filltrav(r)
node *r;
{
  /* traverse to fill in interior node states */
  if (r->tip)
    return;
  filltrav(r->next->back);
  filltrav(r->next->next->back);
  fillin(r);
}  /* filltrav */

/* Local variables for hyptrav: */
struct LOC_hyptrav {
  struct LOC_hypstates *LINK;
  node *r;
  boolean bottom, nonzero;
  gbit *zerobelow, *onebelow;
} ;

void hyprint(Hyptrav,unknown,dohyp)
struct LOC_hyptrav *Hyptrav;
boolean *unknown;
bitptr dohyp;
{
  /* print out states at node */
  long i, j, k;
  char l;
  boolean dot, a0, a1, s0, s1;

  if (Hyptrav->bottom)
    fprintf(outfile, "root   ");
  else
    fprintf(outfile, "%3ld    ", Hyptrav->r->back->index - spp);
  if (Hyptrav->r->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[Hyptrav->r->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", Hyptrav->r->index - spp);
  if (Hyptrav->nonzero)
    fprintf(outfile, "   yes    ");
  else if (*unknown)
    fprintf(outfile, "    ?     ");
  else
    fprintf(outfile, "   no     ");
  for (j = 1; j <= (chars); j++) {
    newline(j, 40, (int)(nmlngth + 17));
    k = (j - 1) / bits + 1;
    l = (j - 1) % bits + 1;
    dot = (((1L << l) & dohyp[k - 1]) == 0 && guess[j - 1] == '?');
    s0 = (((1L << l) & Hyptrav->r->statezero[k - 1]) != 0);
    s1 = (((1L << l) & Hyptrav->r->stateone[k - 1]) != 0);
    a0 = (((1L << l) & Hyptrav->zerobelow->bits_[k - 1]) != 0);
    a1 = (((1L << l) & Hyptrav->onebelow->bits_[k - 1]) != 0);
    dot = (dot || (a1 == s1 && a0 == s0));
    if (dot)
      putc('.', outfile);
    else {
      if (s0) {
        if (s1)
          putc('P', outfile);
        else
          putc('0', outfile);
      } else if (s1)
        putc('1', outfile);
      else
        putc('?', outfile);
    }
    if (j % 5 == 0)
      putc(' ', outfile);
  }
  putc('\n', outfile);
}  /* hyprint */

void hyptrav(r_,unknown,dohyp)
node *r_;
boolean *unknown;
bitptr dohyp;
{
  /*  compute, print out states at one interior node */
  struct LOC_hyptrav HypVars;
  long i;

  HypVars.r = r_;
  gnu(&HypVars.zerobelow);
  gnu(&HypVars.onebelow);
  if (!HypVars.r->tip)
    correct(HypVars.r);
  HypVars.bottom = (HypVars.r->back == NULL);
  HypVars.nonzero = false;
  if (HypVars.bottom) {
    memcpy(HypVars.zerobelow->bits_, zeroanc, words*sizeof(long));
    memcpy(HypVars.onebelow->bits_, oneanc, words*sizeof(long));
  } else {
    memcpy(HypVars.zerobelow->bits_,
          treenode[HypVars.r->back->index - 1]->statezero, words*sizeof(long));
    memcpy(HypVars.onebelow->bits_,
          treenode[HypVars.r->back->index - 1]->stateone, words*sizeof(long));
  }
  for (i = 0; i < (words); i++)
    HypVars.nonzero = (HypVars.nonzero ||
                      ((HypVars.r->stateone[i] & HypVars.zerobelow->bits_[i])
                      | (HypVars.r->statezero[i]
                        & HypVars.onebelow->bits_[i])) != 0);
  hyprint(&HypVars,unknown,dohyp);
  if (!HypVars.r->tip) {
    hyptrav(HypVars.r->next->back, unknown,dohyp);
    hyptrav(HypVars.r->next->next->back, unknown,dohyp);
  }
  chuck(HypVars.zerobelow);
  chuck(HypVars.onebelow);
}  /* hyptrav */

void hypstates()
{
  boolean unknown  = false;
  bitptr dohyp;
  /* fill in and describe states at interior nodes */
  long i, j, k;

  for (i = 0; i < (words); i++) {
    zeroanc[i] = 0;
    oneanc[i] = 0;
  }
  for (i = 0; i < (chars); i++) {
    j = i / bits + 1;
    k = i % bits + 1;
    if (guess[i] == '0')
      zeroanc[j - 1] = ((long)zeroanc[j - 1]) | (1L << k);
    if (guess[i] == '1')
      oneanc[j - 1] = ((long)oneanc[j - 1]) | (1L << k);
    unknown = (unknown || guess[i] == '?');
  }
  dohyp = (bitptr)Malloc(words*sizeof(long));
  for (i = 0; i < words; i++)
    dohyp[i] = zeroanc[i] | oneanc[i];
  filltrav(root);
  fprintf(outfile, "From    To     Any Steps?");
  fprintf(outfile, "    State at upper node\n");
  fprintf(outfile, "                            ");
  fprintf(outfile, " ( . means same as in the node below it on tree)\n\n");
  hyptrav(root, &unknown,dohyp);
  free(dohyp);
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
     each character */
  long i, j, k;

  if (treeprint)
    fprintf(outfile, "\nrequires a total of %10.3f\n", -like);
  if (stepbox) {
    putc('\n', outfile);
    if (weights)
      fprintf(outfile, " weighted");
    if (dollo)
      fprintf(outfile, " reversions ");
    else
      fprintf(outfile, " polymorphisms ");
    fprintf(outfile, "in each character:\n");
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
          fprintf(outfile, "%4ld", numsteps[k - 1] + extras[k - 1]);
      }
      putc('\n', outfile);
    }
  }
  putc('\n', outfile);
  if (questions) {
    fprintf(outfile, "best guesses of ancestral states:\n");
    fprintf(outfile, "      ");
    for (i = 0; i <= 9; i++)
      fprintf(outfile, "%2ld", i);
    fprintf(outfile, "\n     *--------------------\n");
    for (i = 0; i <= (chars / 10); i++) {
      putc(' ', outfile);
      fprintf(outfile, "%4ld!", i * 10);
      for (j = 0; j <= 9; j++) {
        if (i * 10 + j == 0 || i * 10 + j > chars)
          fprintf(outfile, "  ");
        else
          fprintf(outfile, " %c", guess[i * 10 + j - 1]);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
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


void maketree()
{
  /* constructs a binary tree from the pointers in treenode.
     adds each node at location which yields highest "likelihood"
     then rearranges the tree for greatest "likelihood" */
  long i, j, k, numtrees, num, sumw;
  double gotlike, sum, sum2, sd;
  node *item, *nufork, *dummy;
  double TEMP;

  steps = (bitptr)Malloc(words*sizeof(long));
  fullset = (1L << (bits + 1)) - (1L << 1);
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
    add(treenode[enterorder[0] - 1], treenode[enterorder[1] - 1],
        treenode[spp]);
    if (progress) {
      printf("\nAdding species:\n");
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
      bestyet = -10.0 * spp * chars;
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
          for (j = 1; j <= (nonodes); j++)
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
            bestyet = -10.0 * spp * chars;
            item = treenode[j];
            if (item != root) {
              nufork = treenode[j]->back;
              re_move(&item, &nufork);
              there = root;
              addpreorder(root, item, nufork);
              add(there, item, nufork);
            }
            if (progress)
              putchar('.');
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
      for (i = 0; i <= (nextree - 2); i++) {
        root = treenode[0];
        add(treenode[0], treenode[1], treenode[spp]);
        for (j = 3; j <= spp; j++) {
          add(treenode[bestrees[i][j - 1] - 1], treenode[j - 1],
              treenode[spp + j - 2]);}
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
      fprintf(outfile, ":\n");
    }
    which = 1;
    while (which <= numtrees) {
      treeread();
      evaluate(root);
      printree();
      describe();
      which++;
    }
    fprintf(outfile, "\n\n");
    if (numtrees > 1 && chars > 1) {
      if (numtrees > maxuser) {
        printf("TOO MANY USER-DEFINED TREES");
        printf("  test performed on only the first%4ld of them\n",
               (long)maxuser);
      }
      fprintf(outfile, "Tree    Steps   Diff Steps   Its S.D.");
      fprintf(outfile, "   Significantly worse?\n\n");
      if (numtrees > maxuser)
        num = maxuser;
      else
        num = numtrees;
      for (which = 1; which <= num; which++) {
        fprintf(outfile, "%4ld%10.1f", which, nsteps[which - 1]);
        if (which == minwhich)
          fprintf(outfile, "  <------- best\n");
        else {
          sumw = 0;
          sum = 0.0;
          sum2 = 0.0;
          for (j = 0; j < (chars); j++) {
            if (weight[j] > 0) {
              sumw += weight[j];
              sum += fsteps[which - 1][j] - fsteps[minwhich - 1][j];
              TEMP = fsteps[which - 1][j] - fsteps[minwhich - 1][j];
              sum2 += TEMP * TEMP;
            }
          }
          TEMP = sum / sumw;
          sd = sqrt(sumw / (sumw - 1.0) * (sum2 - sum * sum / sumw));
          fprintf(outfile, "%8.1f%15.5f",
                  nsteps[which - 1] - minsteps, sd);
          if (sum > 1.95996 * sd)
            fprintf(outfile, "           Yes\n");
          else
            fprintf(outfile, "           No\n");
        }
      }
      fprintf(outfile, "\n\n");
    }
  }
  if (jumb == njumble) {
    if (progress) {
      printf("Output written to output file\n\n");
      if (trout)
        printf("Trees also written onto tree file\n\n");
    }
    free(steps);
  }
}  /* maketree */


main(argc, argv)
int argc;
Char *argv[];
{  /* Dollo or polymorphism parsimony by uphill search */
char infilename[100],outfilename[100],trfilename[100];
#ifdef MAC
  macsetup("Dollop","");
  argv[0] = "Dollop";
#endif
  /* reads in spp, chars, and the data. Then calls maketree to
     construct the tree */
  openfile(&infile,INFILE,"r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);

  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  garbage = NULL;
  mulsets = false;
  datasets = 1;
  firstset = true;
  doinit();
  if (trout)
    openfile(&treefile,TREEFILE,"w",argv[0],NULL);
  extras = (steptr)Malloc(chars*sizeof(long));
  weight = (steptr)Malloc(chars*sizeof(long));
  threshwt = (double *)Malloc(chars*sizeof(double));
  if (usertree) {
    fsteps = (double **)Malloc(maxuser*sizeof(double *));
    for (j = 1; j <= maxuser; j++)
      fsteps[j - 1] = (double *)Malloc(chars*sizeof(double));
  }
  bestrees = (long **)Malloc(maxtrees*(sizeof(long *)));
  for (j = 1; j <= maxtrees; j++)
    bestrees[j - 1] = (long *)Malloc(spp*sizeof(long));
  numsteps = (steptr)Malloc(chars*sizeof(long));
  nayme = (Char **)Malloc(spp*sizeof(Char *));
  for (j = 1; j <= spp; j++)
    nayme[j - 1] = (Char *)Malloc(nmlngth*sizeof(Char));
  enterorder = (long *)Malloc(spp*sizeof(long));
  naymes = (boolean *)Malloc(spp*sizeof(boolean));
  place = (long *)Malloc(nonodes*sizeof(long));
  ancone = (boolean *)Malloc(chars*sizeof(boolean));
  anczero = (boolean *)Malloc(chars*sizeof(boolean));
  ancone0 = (boolean *)Malloc(chars*sizeof(boolean));
  anczero0 = (boolean *)Malloc(chars*sizeof(boolean));
  numsone = (steptr)Malloc(chars*sizeof(long));
  numszero = (steptr)Malloc(chars*sizeof(long));
  guess = (Char *)Malloc(chars*sizeof(Char));
  zeroanc = (bitptr)Malloc(words*sizeof(long));
  oneanc = (bitptr)Malloc(words*sizeof(long));
  for (ith = 1; ith <= (datasets); ith++) {
    doinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n\n",ith);
      if (progress)
        printf("\nData set # %ld:\n",ith);
    }
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
  }
  FClose(infile);
  FClose(outfile);
  FClose(treefile);
#ifdef MAC
  fixmacfile(infilename);
  fixmacfile(outfilename);
#endif
  exit(0);
}  /* Dollo or polymorphism parsimony by uphill search */


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

