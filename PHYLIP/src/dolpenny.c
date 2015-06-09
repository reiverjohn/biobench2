#include "phylip.h"

/* version 3.572c. (c) Copyright 1995 by Joseph Felsenstein.
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
typedef long *bitptr;
typedef short *treenumbers;

/* nodes will form a binary tree */

typedef struct node {          /* describes a tip species or an ancestor */
  struct node *next, *back;
  short index;
  boolean tip;                 /* present species are tips of tree       */
  bitptr stateone, statezero;  /* see in PROCEDURE fillin                */
  short xcoord, ycoord, ymin;  /* used by printree                       */
  short ymax;
} node;

typedef node **pointptr;

typedef struct gbit {
  bitptr bits_;
  struct gbit *next;
} gbit;


Static node *root;
Static FILE *infile, *outfile, *treefile;
Static short spp, nonodes, chars, words, howmany, howoften, col, datasets,
	     ith, j;
/* spp = number of species
   nonodes = number of nodes in tree
   chars = number of binary characters
   words = number of words needed
           to represent characters of one organism */
Static boolean weights, thresh, ancvar, questions, dollo, simple, trout,
	       printdata, progress, treeprint, stepbox, ancseq,
	       mulsets, firstset, ibmpc, vt52, ansi;
Static steptr extras, weight;
Static boolean *ancone, *anczero, *ancone0, *anczero0;
Static pointptr treenode;   /* pointers to all nodes in tree */
Static Char **name;   /* names of species */
Static double threshold, fracdone, fracinc;
Static double *threshwt;
Static boolean *added;
Static Char *guess;
Static steptr numsteps, numsone, numszero;
Static gbit *garbage;
Static short **bestorders, **bestrees;
Static short bits = 8*sizeof(long) - 1;

/* Local variables for maketree, propagated globally for C version: */
short nextree, examined, mults;
boolean firsttime, done;
double like, bestyet;
treenumbers current, order;
long fullset;
bitptr zeroanc, oneanc;
bitptr stps;


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



void uppercase(ch)
Char *ch;
{  /* convert a character to upper case -- either ASCII or EBCDIC */
  *ch = (islower (*ch) ? toupper(*ch) : (*ch));
}  /* uppercase */

void newline(i, j, k)
short i, j, k;
{
  /* go to new line if i is a multiple of j, indent k spaces */
  short m;

  if ((i - 1) % j != 0 || i <= 1)
    return;
  putchar('\n');
  for (m = 1; m <= k; m++)
    putchar(' ');
}  /* newline */


void getoptions()
{
  /* interactively set options */
  Char ch;
  boolean done, done1;

  fprintf(outfile, "\nPenny algorithm for Dollo or polymorphism");
  fprintf(outfile, " parsimony, version %s\n",VERSION);
  fprintf(outfile, " branch-and-bound to find all");
  fprintf(outfile, " most parsimonious trees\n\n");
  howoften = often;
  howmany = many;
  simple = true;
  thresh = false;
  threshold = spp;
  trout = true;
  weights = false;
  ancvar = false;
  dollo = true;
  printdata = false;
  progress = true;
  treeprint = true;
  stepbox = false;
  ancseq = false;
  do {
    printf(ansi ? "\033[2J\033[H" :
	   vt52 ? "\033E\033H"    : "\n");
    printf("\nPenny algorithm for Dollo or polymorphism parsimony,");
    printf(" version %s\n",VERSION);
    printf(" branch-and-bound to find all most parsimonious trees\n\n");
    printf("Settings for this run:\n");
    printf("  P                     Parsimony method?  %s\n",
	   (dollo ? "Dollo" : "Polymorphism"));
    printf("  H        How many groups of %4hd trees:%6hd\n",howoften,howmany);
    printf("  F        How often to report, in trees:%5hd\n",howoften);
    printf("  S           Branch and bound is simple?  %s\n",
	   (simple ? "Yes" : "No.  Reconsiders order of species"));
    printf("  T              Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f per char.\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  A   Use ancestral states in input file?  %s\n",
	   (ancvar ? "Yes" : "No"));
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2hd sets\n", datasets);
    else
      printf("  No\n");
    printf("  0   Terminal type (IBM PC, VT52, ANSI)?  %s\n",
	   (ibmpc ? "IBM PC" :
	    ansi  ? "ANSI"   :
	    vt52  ? "VT52"   : "(none)"));

    printf("  1    Print out the data at start of run  %s\n",
	   (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
	   (progress ? "Yes" : "no"));
    printf("  3                        Print out tree  %s\n",
	   (treeprint ? "Yes": "No"));
    printf("  4     Print out steps in each character  %s\n",
	   (stepbox ? "Yes" : "No"));
    printf("  5     Print states at all nodes of tree  %s\n",
	   (ancseq ? "Yes" : "No"));
    printf("  6       Write out trees onto tree file?  %s\n",
	   (trout ? "Yes" : "No"));
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      if (strchr("HMSTAPF1234560",ch)){
	switch (ch) {

	case 'H':
	  do {
	    printf("How many cycles of %4hd trees?\n", howoften);
	    scanf("%hd%*[^\n]", &howmany);
	    getchar();
	  } while (howmany <= 0);
	  break;

	case 'F':
	  do {
	    printf("How trees per cycle?\n");
	    scanf("%hd%*[^\n]", &howoften);
	    getchar();
	  } while (howoften <= 0);
	  break;

	case 'A':
	  ancvar = !ancvar;
	  break;

	case 'P':
	  dollo = !dollo;
	  break;

	case 'S':
	  simple = !simple;
	  break;

	case 'T':
	  thresh = !thresh;
	  if (thresh) {
	    done1 = false;
	    do {
	      printf("What will be the threshold value?\n");
	      scanf("%lg%*[^\n]", &threshold);
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
  } while (!done);
}  /* getoptions */

void inputnumbers()
{
  /* input the numbers of species and of characters */
  fscanf(infile, "%hd%hd", &spp, &chars);
  if (printdata)
    fprintf(outfile, "%2hd species, %3hd  characters\n", spp, chars);
  if (printdata)
    putc('\n', outfile);
  if (progress)
    putchar('\n');
  words = chars / bits + 1;
  nonodes = spp * 2 - 1;
}  /* inputnumbers */


void doinit()
{
  /* initializes variables */
  short i, j;
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
  short i;

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
  /* print out the weights of characters */
  short i, j, k;

  fprintf(outfile, "   Characters are weighted as follows:\n");
  fprintf(outfile, "        ");
  for (i = 0; i <= 9; i++)
    fprintf(outfile, "%3hd", i);
  fprintf(outfile, "\n    *---------------------------------\n");
  for (j = 0; j <= (chars / 10); j++) {
    fprintf(outfile, "%4hd!  ", j * 10);
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

void inputancestors()
{
  /* reads the ancestral states for each character */
  short i;
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
      printf("BAD ANCESTOR STATE: %c AT CHARACTER %4hd\n", ch, i + 1);
      exit(-1);
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
}  /* inputancestors */

void printancestors()
{
  /* print out list of ancestral states */
  short i;

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
  short extranum, i, cursp, curchs;
  boolean avar;

  if (!firstset) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%hd%hd", &cursp, &curchs);
    if (cursp != spp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4hd\n", ith);
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
    if ((ch != 'A' && ch != 'W')){
      printf("ERROR: INCORRECT AUXILIARY OPTIONS LINE WHICH STARTS WITH %c\n",
	     ch);
      exit(-1);}
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
    printf("ERROR: ANCESTOR OPTION CHOSEN IN MENU");
    printf(" WITH NO OPTION A IN INPUT\n");
    exit(-1);
  }
  fprintf(outfile,"%s parsimony method\n\n",dollo ? "Dollo" : "Polymorphism");
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
  if (ancvar && printdata)
    printancestors();
  questions = false;
  if (!thresh)
    threshold = spp;
  for (i = 0; i < (chars); i++) {
    questions = (questions || (ancone[i] && anczero[i]));
    threshwt[i] = threshold * weight[i];
  }
}  /* inputoptions */

void inputdata()
{
  /* input the names and character state data for species */
  short i, j, l;
  char k;
  Char charstate;
  /* possible states are
                        '0', '1', 'P', 'B', and '?' */
  node *p;

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
	name[i - 1][j] = getc(infile);
	if (eof(infile) || eoln(infile)){
	  printf("ERROR: END-OF-LINE OR END-OF-FILE");
	  printf(" IN THE MIDDLE OF A SPECIES NAME\n");
	  exit(-1);}
      }
      if (printdata) {
	for (j = 0; j < nmlngth; j++)
	  putc(name[i - 1][j], outfile);
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
	  printf("WARNING -- BAD CHARACTER STATE: %c ",charstate);
	  printf("AT CHARACTER %5hd OF SPECIES %3hd\n",j,i);
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

typedef double *valptr;
typedef short *placeptr;

void add(below, newtip, newfork)
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
}  /* add */

void re_move(item, fork)
node **item, **fork;
{
  /* removes nodes item and its ancestor, fork, from the tree.
     the new descendant of fork's ancestor is made to be
     fork's second descendant (other than item).  Also
     returns pointers to the deleted nodes, item and fork.
     The global variable root is also updated */
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
  short i;

  for (i = 0; i < (words); i++) {
    p->stateone[i] = p->next->back->stateone[i] |
                     p->next->next->back->stateone[i];
    p->statezero[i] = p->next->back->statezero[i] |
		      p->next->next->back->statezero[i];
  }
}  /* fillin */

void correct(p)
node *p;
{  /* get final states for intermediate nodes */
  short i;
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



void postorder(p)
node *p;
{
  /* traverses a binary tree, calling PROCEDURE fillin at a
     node's descendants before calling fillin at the node */
  if (p->tip)
    return;
  postorder(p->next->back);
  postorder(p->next->next->back);
  fillin(p);
}  /* postorder */

void count(stps)
long *stps;
{
  /* counts the number of steps in a branch of the tree.
     The program spends much of its time in this PROCEDURE */
  short i, j, l;

  j = 1;
  l = 0;
  for (i = 0; i < (chars); i++) {
    l++;
    if (l > bits) {
      l = 1;
      j++;
    }
    if (((1L << l) & stps[j - 1]) != 0) {
      if (((1L << l) & zeroanc[j - 1]) != 0)
	numszero[i] += weight[i];
      else
	numsone[i] += weight[i];
    }
  }
}  /* count */

void preorder(p)
node *p;
{
  /* go back up tree setting up and counting interior node
     states */
  short i;

  if (!p->tip) {
    correct(p);
    preorder(p->next->back);
    preorder(p->next->next->back);
  }
  if (p->back == NULL)
    return;
  if (dollo) {
    for (i = 0; i < (words); i++)
      stps[i] = (treenode[p->back->index - 1]->stateone[i] & p->statezero[i] &
		 zeroanc[i]) |
		(treenode[p->back->index - 1]->statezero[i] & p->stateone[i] &
		 fullset & (~zeroanc[i]));
  } else {
    for (i = 0; i < (words); i++)
      stps[i] = treenode[p->back->index - 1]->stateone[i] &
		treenode[p->back->index - 1]->statezero[i] & p->stateone[i] &
		p->statezero[i];
  }
  count(stps);
}  /* preorder */

void evaluate(r)
node *r;
{
  /* Determines the number of losses or polymorphisms needed
     for a tree. This is the minimum number needed to evolve
     chars on this tree */
  short i, stepnum, smaller;
  double sum;

  sum = 0.0;
  for (i = 0; i < (chars); i++) {
    numszero[i] = 0;
    numsone[i] = 0;
  }
  for (i = 0; i < (words); i++)
    zeroanc[i] = fullset;
  postorder(r);
  preorder(r);
  for (i = 0; i < (words); i++)
    zeroanc[i] = 0;
  postorder(r);
  preorder(r);
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
      sum += stepnum;
    else
      sum += threshwt[i];
    guess[i] = '?';
    if (!ancone[i] || (anczero[i] && numszero[i] < numsone[i]))
      guess[i] = '0';
    else if (!anczero[i] || (ancone[i] && numsone[i] < numszero[i]))
      guess[i] = '1';
  }
  if (examined == 0 && mults == 0)
    bestyet = -1.0;
  like = sum;
}  /* evaluate */



void addtraverse(a, b, c,place,valyew,n)
node *a, *b, *c;
valptr valyew;
placeptr place;
short *n;
{
  /* traverse all places to add b */
  if (done)
    return;
  add(&a, &b, &c);
  (*n)++;
  evaluate(root);
  examined++;
  if (examined == howoften) {
    examined = 0;
    mults++;
    if (mults == howmany)
      done = true;
    if (progress) {
      printf("%6hd",mults);
      if (bestyet >= 0)
	printf("%18.5f", bestyet);
      else
	printf("         -        ");
      printf("%17hd%20.2f\n", nextree - 1, fracdone * 100);
    }
  }
  valyew[*n - 1] = like;
  place[*n - 1] = a->index;
  re_move(&b, &c);
  if (!a->tip) {
    addtraverse(a->next->back, b, c, place,valyew,n);
    addtraverse(a->next->next->back, b, c, place,valyew,n);
  }
}  /* addtraverse */

void shellsort(a, b, n)
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

void addit(m)
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
    addtraverse(root, treenode[order[m - 1] - 1], treenode[spp + m - 2],
		place,valyew,&n);
    besttoadd = order[m - 1];
    memcpy(bestplace, place, nonodes*sizeof(short));
    memcpy(bestval, valyew, nonodes*sizeof(double));
  } else {
    bestsum = -1.0;
    for (i = 1; i <= (spp); i++) {
      if (!added[i - 1]) {
	n = 0;
	added[i - 1] = true;
	addtraverse(root, treenode[i - 1], treenode[spp + m - 2],
		    place,valyew,&n);
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
  for (i = 0; i < (n); i++) {
    if (valyew[i] <= bestyet || bestyet < 0.0)
      n1++;
  }
  if (n1 > 0)
    fracinc /= n1;
  for (i = 0; i < n; i++) {
    if (valyew[i] <=bestyet ||bestyet < 0.0) {
      current[m - 1] = place[i];
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


void coordinates(p, tipy)
node *p;
short *tipy;
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = *tipy;
    p->ymin = *tipy;
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
      } while (!(done ||r == p));
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

void printree()
{
  /* prints out diagram of the tree */
  short i,tipy;
  double scale;

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
  putc('\n', outfile);
}  /* printree */


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

void hyprint(dohyp,unknown,bottom,nonzero,zerobelow,onebelow,r)
bitptr dohyp;
boolean unknown,bottom,nonzero;
gbit *zerobelow,*onebelow;
node *r;
{
  /* print out states at node */
  short i, j, k;
  char l;
  boolean dot, a0, a1, s0, s1;

  if (bottom)
    fprintf(outfile, "root   ");
  else
    fprintf(outfile, "%3hd    ", r->back->index - spp);
  if (r->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(name[r->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4hd      ", r->index - spp);
  if (nonzero)
    fprintf(outfile, "   yes    ");
  else if (unknown)
    fprintf(outfile, "    ?     ");
  else
    fprintf(outfile, "   no     ");
  for (j = 1; j <= (chars); j++) {
    newline(j, 40, (int)(nmlngth + 17));
    k = (j - 1) / bits + 1;
    l = (j - 1) % bits + 1;
    dot = (((1L << l) & dohyp[k - 1]) == 0 && guess[j - 1] == '?');
    s0 = (((1L << l) & r->statezero[k - 1]) != 0);
    s1 = (((1L << l) & r->stateone[k - 1]) != 0);
    a0 = (((1L << l) & zerobelow->bits_[k - 1]) != 0);
    a1 = (((1L << l) & onebelow->bits_[k - 1]) != 0);
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

void hyptrav(r_, unknown,dohyp)
node *r_;
boolean unknown;
bitptr dohyp;
{
  /*  compute, print out states at one interior node */
node *r;
boolean bottom, nonzero;
gbit *zerobelow, *onebelow;

  short i;
  r = r_;
  gnu(&zerobelow);
  gnu(&onebelow);
  if (!r->tip)
    correct(r);
  bottom = (r->back == NULL);
  nonzero = false;
  if (bottom) {
    memcpy(zerobelow->bits_, zeroanc, words*sizeof(long));
    memcpy(onebelow->bits_, oneanc, words*sizeof(long));
  } else {
    memcpy(zerobelow->bits_, treenode[r->back->index - 1]->statezero,
	   words*sizeof(long));
    memcpy(onebelow->bits_, treenode[r->back->index - 1]->stateone,
	   words*sizeof(long));
  }
  for (i = 0; i < (words); i++)
    nonzero = (nonzero || ((r->stateone[i] & zerobelow->bits_[i]) |
			       (r->statezero[i] & onebelow->bits_[i])) != 0);
  hyprint(dohyp,unknown,bottom,nonzero,zerobelow,onebelow,r);
  if (!r->tip) {
    hyptrav(r->next->back, unknown,dohyp);
    hyptrav(r->next->next->back, unknown,dohyp);
  }
  chuck(zerobelow);
  chuck(onebelow);
}  /* hyptrav */

void hypstates()
{
/* fill in and describe states at interior nodes */
/* Local variables for hypstates:                */
  boolean unknown;
  bitptr dohyp;
  short i, j, k;

  for (i = 0; i < (words); i++) {
    zeroanc[i] = 0;
    oneanc[i] = 0;
  }
  unknown = false;
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
  for (i = 0; i < (words); i++)
    dohyp[i] = zeroanc[i] | oneanc[i];
  filltrav(root);
  fprintf(outfile, "From    To     Any Steps?    ");
  fprintf(outfile, "State at upper node\n");
  fprintf(outfile, "                             ");
  fprintf(outfile, "( . means same as in the node below it on tree)\n\n");
  hyptrav(root,unknown,dohyp);
  free(dohyp);
}  /* hypstates */

void treeout(p)
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

void describe()
{
  /* prints ancestors, steps and table of numbers of steps in
     each character */
  short i, j, k;

  if (stepbox) {
    putc('\n', outfile);
    if (weights)
      fprintf(outfile, "weighted");
    if (dollo)
      fprintf(outfile, " reversions ");
    else
      fprintf(outfile, " polymorphisms ");
    fprintf(outfile, "in each character:\n");
    fprintf(outfile, "      ");
    for (i = 0; i <= 9; i++)
      fprintf(outfile, "%4hd", i);
    fprintf(outfile, "\n     *-----------------------------------------\n");
    for (i = 0; i <= (chars / 10); i++) {
      fprintf(outfile, "%5hd", i * 10);
      putc('!', outfile);
      for (j = 0; j <= 9; j++) {
	k = i * 10 + j;
	if (k == 0 || k > chars)
	  fprintf(outfile, "    ");
	else
	  fprintf(outfile, "%4hd", numsteps[k - 1] + extras[k - 1]);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  if (questions) {
    fprintf(outfile, "best guesses of ancestral states:\n");
    fprintf(outfile, "      ");
    for (i = 0; i <= 9; i++)
      fprintf(outfile, "%2hd", i);
    fprintf(outfile, "\n     *--------------------\n");
    for (i = 0; i <= (chars / 10); i++) {
      fprintf(outfile, "%5hd!", i * 10);
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
  /* tree construction recursively by branch and bound */
  short i, j, k;
  node *dummy;

  fullset = (1L << (bits + 1)) - (1L << 1);
  if (progress) {
    printf("\nHow many\n");
    printf("trees looked                                       Approximate\n");
    printf("at so far      Length of        How many           percentage\n");
    printf("(multiples     shortest tree    trees this long    searched\n");
    printf("of %4hd):      found so far     found so far       so far\n",
	   howoften);
    printf("----------     ------------     ------------       ------------\n");
  }
  done = false;
  mults = 0;
  examined = 0;
  nextree = 1;
  root = treenode[0];
  firsttime = true;
  for (i = 0; i < (spp); i++)
    added[i] = false;
  added[0] = true;
  order[0] = 1;
  k = 2;
  fracdone = 0.0;
  fracinc = 1.0;
  bestyet = -1.0;
  free(stps);
  stps = (bitptr)Malloc(words*sizeof(long));
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
    fprintf(outfile, "\nrequires a total of %18.3f\n\n", bestyet);
    if (nextree == 2)
      fprintf(outfile, "One most parsimonious tree found:\n");
    else
      fprintf(outfile, "%5hd trees in all found\n", nextree - 1);
  }
  if (nextree > maxtrees + 1) {
    if (treeprint)
      fprintf(outfile, "here are the first%4ld of them\n", (long)maxtrees);
    nextree = maxtrees + 1;
  }
  if (treeprint)
    putc('\n', outfile);
  for (i = 0; i < (spp); i++)
    added[i] = true;
  for (i = 0; i <= (nextree - 2); i++) {
    for (j = k; j <= (spp); j++)
      add(&treenode[bestrees[i][j - 1] - 1],
	  &treenode[bestorders[i][j - 1] - 1], &treenode[spp + j - 2]);
    evaluate(root);
    printree();
    describe();
    for (j = k - 1; j < (spp); j++)
      re_move(&treenode[bestorders[i][j] - 1], &dummy);
  }
  if (progress) {
    printf("\nOutput written to output file\n\n");
    if (trout)
      printf("Trees also written onto tree file\n\n");
  }
  free(stps);
}  /* maketree */


main(argc, argv)
int argc;
Char *argv[];
{  /* branch-and-bound method for Dollo, polymorphism parsimony */
char infilename[100],outfilename[100],trfilename[100];
  /* Reads in the number of species, number of characters,
     options and data.  Then finds all most parsimonious trees */
#ifdef MAC
  macsetup("Dolpenny","");
  argv[0] = "Dolpenny";
#endif
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
  extras = (short *)Malloc(chars*sizeof(short));
  weight = (short *)Malloc(chars*sizeof(short));
  threshwt = (double *)Malloc(chars*sizeof(double));
  guess = (Char *)Malloc(chars*sizeof(Char));
  numsteps = (short *)Malloc(chars*sizeof(short));
  numszero = (short *)Malloc(chars*sizeof(short));
  numsone = (short *)Malloc(chars*sizeof(short));
  bestorders = (short **)Malloc(maxtrees*sizeof(short *));
  bestrees = (short **)Malloc(maxtrees*sizeof(short *));
  for (j = 1; j <= maxtrees; j++) {
    bestorders[j - 1] = (short *)Malloc(spp*sizeof(short));
    bestrees[j - 1] = (short *)Malloc(spp*sizeof(short));
  }
  current = (treenumbers)Malloc(spp*sizeof(short));
  order = (treenumbers)Malloc(spp*sizeof(short));
  name = (Char **)Malloc(spp*sizeof(Char *));
  for (j = 1; j <= spp; j++)
    name[j - 1] = (Char *)Malloc(nmlngth*sizeof(Char));
  added = (boolean *)Malloc(nonodes*sizeof(boolean));
  ancone = (boolean *)Malloc(chars*sizeof(boolean));
  anczero = (boolean *)Malloc(chars*sizeof(boolean));
  ancone0 = (boolean *)Malloc(chars*sizeof(boolean));
  anczero0 = (boolean *)Malloc(chars*sizeof(boolean));
  zeroanc = (bitptr)Malloc(words*sizeof(long));
  oneanc = (bitptr)Malloc(words*sizeof(long));
  if (trout)
    openfile(&treefile,TREEFILE,"w",argv[0],trfilename);
  for (ith = 1; ith <= datasets; ith++) {
    doinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1) {
      fprintf(outfile, "Data set # %hd:\n\n",ith);
      if (progress)
        printf("\nData set # %hd:\n",ith);
    }
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
}  /* branch-and-bound method for Dollo, polymorphism parsimony */


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

