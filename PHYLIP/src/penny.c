#include "phylip.h"

/* version 3.572c. (c) Copyright 1995 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define nmlngth         10   /* number of characters in species name        */
#define maxtrees        1000  /* maximum number of trees to be printed out   */
#define often           100  /* how often to notify how many trees examined */
#define many            1000 /* how many multiples of howoften before stop  */
#define down            2

#define ibmpc0          false
#define ansi0           true
#define vt520           false


typedef short *steptr;
typedef long *bitptr;
typedef short *treenumbers;

typedef struct gbit {
  bitptr bits_;
  struct gbit *next;
} gbit;

/* nodes will form a binary tree */

typedef struct node {             /* describes a tip species or an ancestor */
  struct node *next,*back;
  short index;
  boolean tip, visited;           /* present species are tips of tree        */
  bitptr fulstte1, empstte1, fulstte0, empstte0, fulsteps, empsteps;
    /* see in PROCEDURE fillin */
  short xcoord, ycoord, ymin;     /* used by printree                        */
  short ymax;
} node;

typedef node **pointptr;
Static node *root;
Static FILE *infile, *outfile, *treefile;
Static short spp, nonodes, chars, words, rno, howmany, howoften, col,
	     datasets, ith, j, outgrno;

/* spp = number of species                                                   *
 * nonodes = number of nodes in tree                                         *
 *  chars = number of binary characters                                      *
 *  words = number of words needed to represent characters of one organism   *
 *  outgrno indicates outgroup                                               */

Static boolean weights, thresh, ancvar, questions, allsokal, allwagner,
	       mixture, simple, trout, noroot, didreroot, outgropt,
	       printdata, progress, treeprint, stepbox, ancseq, mulsets,
	       firstset, ibmpc, vt52, ansi;
Static steptr extras, weight;
Static boolean *ancone, *anczero, *ancone0, *anczero0;
Static pointptr treenode;         /* pointers to all nodes in tree   */
Static Char **name;   /* names of species                */
Static double threshold, fracdone, fracinc;
Static double *threshwt;
Static bitptr wagner, wagner0;
Static boolean *added;
Static Char *guess;
Static steptr numsteps;
Static short **bestorders, **bestrees;
Static steptr numsone, numszero;
Static gbit *garbage;
Static bits = sizeof(long)*8 - 1;

struct hyptrav_vars {
  node *r;
  boolean bottom, maybe, nonzero;
  gbit *zerobelow, *onebelow;
} ;

short nextree, examined, mults;

boolean firsttime, done, full;
double like, bestyet;
treenumbers current, order;
long fullset;
bitptr zeroanc, oneanc;
bitptr suppsteps;


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
     *ch = (islower(*ch) ? toupper(*ch) : (*ch));
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


void inputnumbers()
{
  /* input the numbers of species and of characters */
  fscanf(infile, "%hd%hd", &spp, &chars);
  if (printdata)
    fprintf(outfile, "%2hd species, %3hd  characters\n", spp, chars);
  if (printdata)
    putc('\n', outfile);
  putchar('\n');
  words = chars / bits + 1;
  nonodes = spp * 2 - 1;
}  /* inputnumbers */

void getoptions()
{
  /* interactively set options */
  Char ch;
  boolean done, done1;

  fprintf(outfile, "\nPenny algorithm, version %s\n",VERSION);
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
  ancvar = false;
  allsokal = false;
  allwagner = true;
  mixture = false;
  printdata = false;
  progress = true;
  treeprint = true;
  stepbox = false;
  ancseq = false;
  for(;;) {
    printf(ansi ? "\033[2J\033[H" :
	   vt52 ? "\033E\033H"    : "\n");
    printf("\nPenny algorithm, version %s\n",VERSION);
    printf(" branch-and-bound to find all most parsimonious trees\n\n");
    printf("Settings for this run:\n");
    printf("  X                     Use Mixed method?  %s\n",
	   mixture ? "Yes" : "No");
    printf("  P                     Parsimony method?  %s\n",
	   (allwagner && !mixture)  ? "Wagner"       :
	   (!(allwagner || mixture)) ? "Camin-Sokal"  : "(methods in mixture");
    printf("  F        How often to report, in trees:%5hd\n",howoften);
    printf("  H        How many groups of%5hd trees:%6hd\n",howoften,howmany);
    printf("  O                        Outgroup root?");
    if (outgropt)
      printf("  Yes, at species number%3hd\n", outgrno);
    else
      printf("  No, use as outgroup species%3hd\n", outgrno);
    printf("  S           Branch and bound is simple?  %s\n",
	   simple ? "Yes" : "No. reconsiders order of species");
    printf("  T              Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f per char.\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  A   Use ancestral states in input file?  %s\n",
	   ancvar ? "Yes" : "No");
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2hd sets\n", datasets);
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
    printf("\nAre these settings correct?");
    printf(" (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (strchr("HFSOMPATX1234560",ch)){
      switch (ch) {
	
      case 'X':
	mixture = !mixture;
	break;
	
      case 'P':
	allwagner = !allwagner;
	break;
	
      case 'A':
	ancvar = !ancvar;
	break;
	
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
	    done1 = (outgrno >= 1 || outgrno <= spp);
	    if (!done1) {
	      printf("BAD OUTGROUP NUMBER: %4hd\n", outgrno);
	      printf("  Must be in range 1 -%2hd\n", spp);
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
  allsokal = (!allwagner && !mixture);
}  /* getoptions */


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
    treenode[i]->fulstte1 = (bitptr)Malloc(words*sizeof(long));
    treenode[i]->fulstte0 = (bitptr)Malloc(words*sizeof(long));
    treenode[i]->empstte1 = (bitptr)Malloc(words*sizeof(long));
    treenode[i]->empstte0 = (bitptr)Malloc(words*sizeof(long));
    treenode[i]->fulsteps = (bitptr)Malloc(words*sizeof(long));
    treenode[i]->empsteps = (bitptr)Malloc(words*sizeof(long));
  }
  for (i = spp; i < (nonodes); i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc(sizeof(node));
      p->fulstte1 = (bitptr)Malloc(words*sizeof(long));
      p->fulstte0 = (bitptr)Malloc(words*sizeof(long));
      p->empstte1 = (bitptr)Malloc(words*sizeof(long));
      p->empstte0 = (bitptr)Malloc(words*sizeof(long));
      p->fulsteps = (bitptr)Malloc(words*sizeof(long));
      p->empsteps = (bitptr)Malloc(words*sizeof(long));
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    treenode[i] = p;
  }
}  /* doinit */

void inputmixture()
{
  /* input mixture of methods */
  short i, j, k;
  Char ch;
  boolean wag;

  for (i = 1; i < nmlngth; i++)
    ch = getc(infile);
  for (i = 0; i < (words); i++)
    wagner0[i] = 0;
  j = 0;
  k = 1;
  for (i = 1; i <= (chars); i++) {
    do {
      if (eoln(infile)) {
	fscanf(infile, "%*[^\n]");
	getc(infile);
      }
      ch = getc(infile);
    } while (ch == ' ');
    uppercase(&ch);
    wag = false;
    if (ch == 'W' || ch == '?')
      wag = true;
    else if (ch == 'S' || ch == 'C')
      wag = false;
    else {
      printf("BAD METHOD: %c\n", ch);
      exit(-1);
    }
    j++;
    if (j > bits) {
      j = 1;
      k++;
    }
    if (wag)
      wagner0[k - 1] = ((long)wagner0[k - 1]) | (1L << j);
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
}  /* inputmixture */

void printmixture()
{
  /* print out list of parsimony methods */
  short i, k, l;

  fprintf(outfile, "Parsimony methods:\n");
  l = 0;
  k = 1;
  for (i = 1; i <= nmlngth + 3; i++)
    putc(' ', outfile);
  for (i = 1; i <= (chars); i++) {
    newline(i, 55, (int)(nmlngth + 3));
    l++;
    if (l > bits) {
      l = 1;
      k++;
    }
    if (((1L << l) & wagner[k - 1]) != 0)
      putc('W', outfile);
    else
      putc('S', outfile);
    if (i % 5 == 0)
      putc(' ', outfile);
  }
  fprintf(outfile, "\n\n");
}  /* printmixture */

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
  fprintf(outfile, "       ");
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
  boolean avar, mix;

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
  mix = false;
  while (!(eoln(infile))) {
    ch = getc(infile);
    uppercase(&ch);
    if (ch == 'A' || ch == 'M' || ch == 'W')
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
    if (ch != 'A' && ch != 'M' && ch != 'W') {
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
    if (ch == 'M') {
      mix = true;
      if (!mixture) {
	printf("ERROR: MIXTURE OPTION NOT CHOSEN IN MENU");
	printf(" WITH OPTION %c IN INPUT\n", ch);
	exit(-1);
      } else
	inputmixture();
    }
    if (ch == 'W')
      inputweights();
  }
  if (ancvar && !avar) {
    printf("ERROR: ANCESTOR OPTION CHOSEN IN MENU");
    printf(" WITH NO OPTION A IN INPUT\n");
    exit(-1);
  }
  if (mixture && !mix) {
    printf("ERROR: MIXTURE OPTION CHOSEN IN MENU");
    printf(" WITH NO OPTION M IN INPUT\n");
    exit(-1);
  }
  if (weights && printdata)
    printweights();
  for (i = 0; i < (words); i++) {
    if (mixture && mix)
      wagner[i] = wagner0[i];
    else if (allsokal)
      wagner[i] = 0;
    else
      wagner[i] = (1L << (bits + 1)) - (1L << 1);
  }
  if (allsokal)
    fprintf(outfile, "Camin-Sokal parsimony method\n\n");
  if (allwagner)
    fprintf(outfile, "Wagner parsimony method\n\n");
  if (mixture && mix && printdata)
    printmixture();
  for (i = 0; i < (chars); i++) {
    if (!ancvar) {
      anczero[i] = true;
      ancone[i] = (((1L << (i % bits + 1)) & wagner[i / bits]) != 0);
    } else {
      anczero[i] = anczero0[i];
      ancone[i] = ancone0[i];
    }
  }
  if (ancvar && avar && printdata)
    printancestors();
  noroot = true;
  questions = false;
  for (i = 0; i < (chars); i++) {
    if (weight[i] > 0) {
      noroot = (noroot && ancone[i] && anczero[i] &&
	  ((((1L << (i % bits + 1)) & wagner[i / bits]) != 0)
            || threshold <= 2.0));
    }
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
  /* possible states are     '0', '1', 'P', 'B', and '?' */
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
    treenode[i - 1]->visited = false;
    treenode[i - 1]->index = i;
    if (i > spp) {
      p = treenode[i - 1]->next;
      while (p != treenode[i - 1]) {
	p->back = NULL;
	p->tip = false;
	p->visited = false;
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
	treenode[i - 1]->fulstte1[j] = 0;
	treenode[i - 1]->fulstte0[j] = 0;
	treenode[i - 1]->empstte1[j] = 0;
	treenode[i - 1]->empstte0[j] = 0;
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
	  printf("WARNING -- BAD CHARACTER STATE: %c",charstate);
	  printf(" AT CHARACTER %5hd OF SPECIES %3hd\n",charstate, j, i);
	  exit(-1);
	}
	if (printdata) {
	  newline(j, 55, (int)(nmlngth + 3));
	  putc(charstate, outfile);
	  if (j % 5 == 0)
	    putc(' ', outfile);
	}
	if (charstate == '1') {
	  treenode[i - 1]->fulstte1[l - 1] =
	    ((long)treenode[i - 1]->fulstte1[l - 1]) | (1L << k);
	  treenode[i - 1]->empstte1[l - 1] =
	    treenode[i - 1]->fulstte1[l - 1];
	}
	if (charstate == '0') {
	  treenode[i - 1]->fulstte0[l - 1] =
	    ((long)treenode[i - 1]->fulstte0[l - 1]) | (1L << k);
	  treenode[i - 1]->empstte0[l - 1] =
	    treenode[i - 1]->fulstte0[l - 1];
	}
	if (charstate == 'P' || charstate == 'B')
	  extras[j - 1] += weight[j - 1];
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
  node *p;

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
  root->back = NULL;
  p = *newfork;
  do {
    p->visited = false;
    p = p->back;
    if (p != NULL) p = treenode[p->index - 1];
  } while (p != NULL);
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
  q = (*fork)->back;
  (*fork)->back = NULL;
  p = (*fork)->next;
  while (p != *fork) {
    p->back = NULL;
    p = p->next;
  }
  (*item)->back = NULL;
  if (q != NULL) q = treenode[q->index - 1];
  while (q != NULL) {
    q-> visited = false;
    q = q->back;
    if (q != NULL) q = treenode[q->index - 1];
  }
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
  short i;
  long l0, l1, r0, r1, st, wa, za;

  for (i = 0; i < (words); i++) {
    if (full) {
      l0 = p->next->back->fulstte0[i];
      l1 = p->next->back->fulstte1[i];
      r0 = p->next->next->back->fulstte0[i];
      r1 = p->next->next->back->fulstte1[i];
    }
    else {
      l0 = p->next->back->empstte0[i];
      l1 = p->next->back->empstte1[i];
      r0 = p->next->next->back->empstte0[i];
      r1 = p->next->next->back->empstte1[i];
    }
    st = (l1 & r0) | (l0 & r1);
    wa = wagner[i];
    za = zeroanc[i];
    if (full) {
      p->fulstte1[i] = (l1 | r1) & (~(st & (wa | za)));
      p->fulstte0[i] = (l0 | r0) & (~(st & (wa | (fullset & (~za)))));
      p->fulsteps[i] = st;
    }
    else {
      p->empstte1[i] = (l1 | r1) & (~(st & (wa | za)));
      p->empstte0[i] = (l0 | r0) & (~(st & (wa | (fullset & (~za)))));
      p->empsteps[i] = st;
    }
  }
}  /* fillin */


void count(stps)
long *stps;
{
  /* counts the number of steps in a fork of the tree.
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

void postorder(p)
node *p;
{
  /* traverses a binary tree, calling PROCEDURE fillin at a
     node's descendants before calling fillin at the node */
  if (p->tip)
    return;
  postorder(p->next->back);
  postorder(p->next->next->back);
  if (!p->visited) {
    fillin(p);
    if (!full) p->visited = true;
  }
}  /* postorder */

void cpostorder(p)
node *p;
{
  /* traverses a binary tree, calling PROCEDURE count at a
     node's descendants before calling count at the node */
  if (p->tip)
    return;
  cpostorder(p->next->back);
  cpostorder(p->next->next->back);
  if (full) count (p->fulsteps); else count (p->empsteps);
}  /* cpostorder */

void supplement(suppsteps)
bitptr suppsteps;
{
  /* determine minimum number of steps more which will
     be added when rest of species are put in tree */
  short i, j, k, l;
  long defone, defzero, a;

  k = 0;
  for (i = 0; i < (words); i++) {
    defone = 0;
    defzero = 0;
    a = 0;
    for (l = 1; l <= bits; l++) {
      k++;
      if (k <= chars) {
	if (!ancone[k - 1])
	  defzero = ((long)defzero) | (1L << l);
	if (!anczero[k - 1])
	  defone = ((long)defone) | (1L << l);
      }
    }
    for (j = 0; j < (spp); j++) {
      defone |= treenode[j]->empstte1[i] & (~treenode[j]->empstte0[i]);
      defzero |= treenode[j]->empstte0[i] & (~treenode[j]->empstte1[i]);
      if (added[j])
	a |= defone & defzero;
    }
    suppsteps[i] = defone & defzero & (~a);
  }
}  /* supplement */

void evaluate(r)
node *r;
{
  /* Determines the number of steps needed for a tree.
     This is the minimum number needed to evolve chars on
     this tree */
  /* Local variables for evaluate: */
  short i, stepnum, smaller;
  double sum;

  sum = 0.0;
  for (i = 0; i < (chars); i++) {
    numszero[i] = 0;
    numsone[i] = 0;
  }
  supplement(suppsteps);
  for (i = 0; i < (words); i++)
    zeroanc[i] =fullset;
  full = true;
  postorder(r);
  cpostorder(r);
  count(r->fulstte1);
  count(suppsteps);
  for (i = 0; i < (words); i++)
    zeroanc[i] = 0;
  full = false;
  postorder(r);
  cpostorder(r);
  count(r->empstte0);
  count(suppsteps);
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


void addtraverse(a, b, c, m,n,valyew,place)
node *a, *b, *c;
short *m,*n;
valptr valyew;
placeptr place;
{
  /* traverse all places to add b */
  if (done)
    return;
  if ((*m) <= 2 || !(noroot && (a == root || a == root->next->back))) {
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
	printf("%6hd", mults);
	if (bestyet >= 0)
	  printf("%18.5f", bestyet);
	else
	  printf("         -        ");
	printf("%17hd%20.2f\n", nextree - 1, fracdone * 100);
      }
    }
    valyew[(*n) - 1] = like;
    place[(*n) - 1] = a->index;
    re_move(&b, &c);
  }
  if (!a->tip) {
    addtraverse(a->next->back, b, c, m,n,valyew,place);
    addtraverse(a->next->next->back, b, c, m,n,valyew,place);
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
/* Local variables for addit: */
  short n;
  valptr valyew;
  placeptr place;
  /* adds the species one by one, recursively */
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
    for (i = 1; i <= (spp); i++) {
      if (!added[i - 1]) {
	n = 0;
	added[i - 1] = true;
	addtraverse(root, treenode[i - 1], treenode[spp + m - 2], &m,
		    &n,valyew,place);
	added[i - 1] = false;
	sum = 0.0;
	for (j = 0; j < (n); j++)
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
  for (i = 0; i < (n); i++) {
    if (valyew[i] <= bestyet || bestyet < 0.0) {
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

void reroot(outgroup)
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
  if (noroot) {
    fprintf(outfile, "\n  remember:");
    if (didreroot)
      fprintf(outfile, " (although rooted by outgroup)");
    fprintf(outfile, " this is an unrooted tree!\n");
  }
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


void hyprint(hyptravv,unknown)
struct hyptrav_vars *hyptravv;
boolean unknown;
{
  /* print out states at node */
  short i, j, k;
  char l;
  boolean dot, a0, a1, s0, s1;

  if (hyptravv->bottom) {
    if (noroot && !didreroot)
      fprintf(outfile, "       ");
    else
      fprintf(outfile, "root   ");
  } else
    fprintf(outfile, "%3hd    ", hyptravv->r->back->index - spp);
  if (hyptravv->r->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(name[hyptravv->r->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4hd      ", hyptravv->r->index - spp);
  if (hyptravv->bottom && noroot && !didreroot)
    fprintf(outfile, "          ");
  else if (hyptravv->nonzero)
    fprintf(outfile, "   yes    ");
  else if (unknown)
    fprintf(outfile, "    ?     ");
  else if (hyptravv->maybe)
    fprintf(outfile, "  maybe   ");
  else
    fprintf(outfile, "   no     ");
  for (j = 1; j <= (chars); j++) {
    newline(j, 40, (int)(nmlngth + 17));
    k = (j - 1) / bits + 1;
    l = (j - 1) % bits + 1;
    dot = (((1L << l) & wagner[k - 1]) == 0 && guess[j - 1] == '?');
    s0 = (((1L << l) & hyptravv->r->empstte0[k - 1]) != 0);
    s1 = (((1L << l) & hyptravv->r->empstte1[k - 1]) != 0);
    a0 = (((1L << l) & hyptravv->zerobelow->bits_[k - 1]) != 0);
    a1 = (((1L << l) & hyptravv->onebelow->bits_[k - 1]) != 0);
    dot = (dot || ((!hyptravv->bottom || !noroot || didreroot) && a1 == s1 &&
		   a0 == s0 && s0 != s1));
    if (dot)
      putc('.', outfile);
    else {
      if (s0)
	putc('0', outfile);
      else if (s1)
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
boolean unknown;
bitptr dohyp;
{
  /*  compute, print out states at one interior node */
  struct hyptrav_vars vars;
  short i;
  long l0, l1, r0, r1, s0, s1, a0, a1, temp, dh, wa;

  vars.r = r_;
  gnu(&vars.zerobelow);
  gnu(&vars.onebelow);
  vars.bottom = (vars.r->back == NULL);
  vars.maybe = false;
  vars.nonzero = false;
  if (vars.bottom) {
    memcpy(vars.zerobelow->bits_, zeroanc, words*sizeof(long));
    memcpy(vars.onebelow->bits_, oneanc, words*sizeof(long));
  } else {
    memcpy(vars.zerobelow->bits_, treenode[vars.r->back->index - 1]->empstte0,
	   words*sizeof(long));
    memcpy(vars.onebelow->bits_, treenode[vars.r->back->index - 1]->empstte1,
	   words*sizeof(long));
  }
  for (i = 0; i < (words); i++) {
    dh = dohyp[i];
    s0 = vars.r->empstte0[i];
    s1 = vars.r->empstte1[i];
    a0 = vars.zerobelow->bits_[i];
    a1 = vars.onebelow->bits_[i];
    if (!vars.r->tip) {
      wa = wagner[i];
      l0 = vars.r->next->back->empstte0[i];
      l1 = vars.r->next->back->empstte1[i];
      r0 = vars.r->next->next->back->empstte0[i];
      r1 = vars.r->next->next->back->empstte1[i];
      s0 = (wa & ((a0 & l0) | (a0 & r0) | (l0 & r0))) |
	   (dh & fullset & (~wa) & s0);
      s1 = (wa & ((a1 & l1) | (a1 & r1) | (l1 & r1))) |
	   (dh & fullset & (~wa) & s1);
      temp = fullset & (~(s0 | s1 | l1 | l0 | r1 | r0));
      s0 |= temp & a0;
      s1 |= temp & a1;
      vars.r->empstte0[i] = s0;
      vars.r->empstte1[i] = s1;
    }
    vars.maybe = (vars.maybe || (dh & (s0 | s1)) != (a0 | a1));
    vars.nonzero = (vars.nonzero || ((s1 & a0) | (s0 & a1)) != 0);
  }
  hyprint(&vars,unknown);
  if (!vars.r->tip) {
    hyptrav(vars.r->next->back,unknown,dohyp);
    hyptrav(vars.r->next->next->back, unknown,dohyp);
  }
  chuck(vars.zerobelow);
  chuck(vars.onebelow);
}  /* hyptrav */

void hypstates()
{
  /* fill in and describe states at interior nodes */
/* Local variables for hypstates: */
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
    unknown = (unknown || ((((1L << k) & wagner[j - 1]) == 0) &&
		  guess[i] == '?'));
  }
  dohyp = (bitptr)Malloc(words*sizeof(long));
  for (i = 0; i < (words); i++)
    dohyp[i] = wagner[i] | zeroanc[i] | oneanc[i];
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
      fprintf(outfile, "weighted ");
    fprintf(outfile, "steps in each character:\n");
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
  if (questions && (!noroot || didreroot)) {
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
    if (noroot)
      reroot(treenode[outgrno - 1]);
    didreroot = (outgropt && noroot);
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
}  /* maketree */


main(argc, argv)
int argc;
Char *argv[];
{  /* Penny's branch-and-bound method */
   /* Reads in the number of species, number of characters,
     options and data.  Then finds all most parsimonious trees */
char infilename[100],outfilename[100],trfilename[100];
#ifdef MAC
  macsetup("Penny","");
  argv[0] = "Penny";
#endif
  openfile(&infile,INFILE,"r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);
  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  mulsets = false;
  datasets = 1;
  firstset = true;
  garbage = NULL;
  doinit();
  extras = (steptr)Malloc(chars*sizeof(steptr));
  weight = (steptr)Malloc(chars*sizeof(steptr));
  threshwt = (double *)Malloc(chars*sizeof(double));
  bestorders = (short **)Malloc(maxtrees*sizeof(short *));
  bestrees = (short **)Malloc(maxtrees*sizeof(short *));
  for (j = 1; j <= maxtrees; j++) {
    bestorders[j - 1] = (short *)Malloc(spp*sizeof(short));
    bestrees[j - 1] = (short *)Malloc(spp*sizeof(short));
  }
  numsteps = (steptr)Malloc(chars*sizeof(steptr));
  guess = (Char *)Malloc(chars*sizeof(Char));
  numszero = (steptr)Malloc(chars*sizeof(steptr));
  numsone = (steptr)Malloc(chars*sizeof(steptr));
  current = (treenumbers)Malloc(spp*sizeof(short));
  order = (treenumbers)Malloc(spp*sizeof(short));
  name = (Char **)Malloc(spp*sizeof(Char *));
  added = (boolean *)Malloc(nonodes*sizeof(boolean));
  for (j = 1; j <= spp; j++)
    name[j - 1] = (Char *)Malloc(nmlngth*sizeof(Char));
  ancone = (boolean *)Malloc(chars*sizeof(boolean));
  anczero = (boolean *)Malloc(chars*sizeof(boolean));
  ancone0 = (boolean *)Malloc(chars*sizeof(boolean));
  anczero0 = (boolean *)Malloc(chars*sizeof(boolean));
  wagner = (bitptr)Malloc(words*sizeof(long));
  wagner0 = (bitptr)Malloc(words*sizeof(long));
  zeroanc = (bitptr)Malloc(words*sizeof(long));
  oneanc = (bitptr)Malloc(words*sizeof(long));
  suppsteps = (bitptr)Malloc(words*sizeof(long));
  if (trout)
    openfile(&treefile,TREEFILE,"w",argv[0],trfilename);
  for (ith = 1; ith <= datasets; ith++) {
    doinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1) {
      fprintf(outfile, "\nData set # %hd:\n\n",ith);
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
}  /* Penny's branch-and-bound method */


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

