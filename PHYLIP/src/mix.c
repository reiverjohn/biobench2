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

typedef struct node {           /* describes a tip species or an ancestor */
  struct node *next, *back;     /* pointers to nodes                      */
  long index;                   /* number of the node                     */
  boolean tip, bottom,visited;  /* present species are tips of tree       */
  bitptr fulstte1, fulstte0;  /* see in PROCEDURE fillin                */
  bitptr empstte1, empstte0;  /* see in PROCEDURE fillin                */
  bitptr fulsteps,empsteps;
  long xcoord, ycoord, ymin;    /* used by printree                       */
  long ymax;
} node;

typedef long *steptr;
typedef node **pointptr;
typedef long longer[6];
typedef long *placeptr;
char ch;

typedef struct gbit {
  bitptr bits_;
  struct gbit *next;
} gbit;


node *root;
FILE *infile, *outfile, *treefile;
long spp, nonodes, chars, words, inseed, outgrno,  datasets, ith,
            j, l, jumb, njumble;
/* spp = number of species
  nonodes = number of nodes in tree
  chars = number of binary characters
  words = number of words needed to represent characters of one organism
  outgrno indicates outgroup */
boolean jumble, usertree, weights, thresh, ancvar, questions, allsokal,
               allwagner, mixture, trout, noroot, outgropt, didreroot,
                printdata, progress, treeprint, stepbox, ancseq,
               mulsets, firstset, ibmpc, vt52, ansi;
steptr extras, weight;
boolean *ancone, *anczero, *ancone0, *anczero0;
pointptr treenode;   /* pointers to all nodes in tree */
char **nayme;   /* names of species */
double threshold;
double *threshwt;
bitptr wagner, wagner0;
longer seed;
long *enterorder;
double **fsteps;
char *guess;
long **bestrees;
steptr numsteps, numsone, numszero;
gbit *garbage;
long bits =  (8*sizeof(long) - 1);

void openfile(fp,filename,mode,application,perm)
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
  putc('\n', outfile);
  for (m = 1; m <= k; m++)
    putc(' ', outfile);
}  /* newline */


void getoptions()
{
  /* interactively set options */
  long i, inseed0;
  Char ch;
  boolean  done1;

  fprintf(outfile, "\nMixed parsimony algorithm, version %s\n\n",VERSION);
  putchar('\n');
  jumble = false;
  njumble = 1;
  outgrno = 1;
  outgropt = false;
  thresh = false;
  threshold = spp;
  trout = true;
  usertree = false;
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
  for (;;) {
    printf(ansi ?  "\033[2J\033[H" :
           vt52 ?  "\033E\033H"    : "\n");
    printf("\nMixed parsimony algorithm, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
           (usertree ? "No, use user trees in input file" : "Yes"));
    printf("  X                     Use Mixed method?  %s\n",
           (mixture ? "Yes" : "No"));
    printf("  P                     Parsimony method?");
    if (!mixture) {
      printf("  %s\n",(allwagner ? "Wagner" : "Camin-Sokal"));
    } else
      printf("  (methods in mixture)\n");
    if (!usertree) {
      printf("  J     Randomize input order of species?");
      if (jumble)
        printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  O                        Outgroup root?");
    if (outgropt)
      printf("  Yes, at species number%3ld\n", outgrno);
    else
      printf("  No, use as outgroup species%3ld\n", outgrno);
    printf("  T              Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f per char.\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  A   Use ancestral states in input file?  %s\n",
    (ancvar ? "Yes" : "No"));
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  0   Terminal type (IBM PC, VT52, ANSI)?  %s\n",
           (ibmpc ? "IBM PC" :
            ansi  ? "ANSI"   :
            vt52  ? "VT52"   : "(none)"));

    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
          (treeprint ? "Yes" : "No"));
    printf("  4     Print out steps in each character  %s\n",
           (stepbox ? "Yes" : "No"));
    printf("  5     Print states at all nodes of tree  %s\n",
           (ancseq ? "Yes" : "No"));
    printf("  6       Write out trees onto tree file?  %s\n",
           (trout ? "Yes" : "No"));
    printf("\nAre these settings correct? ");
    printf("(type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (strchr("JOTUMPAX1234560",ch)){
      switch (ch) {
	
      case 'U':
	usertree = !usertree;
	break;
	
      case 'X':
	mixture = !mixture;
	break;
	
      case 'P':
	allwagner = !allwagner;
	break;
	
      case 'A':
	ancvar = !ancvar;
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
	}
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


Static Void doinit()
{
  /* initializes variables */
  long i;
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
  long i, j, k;
  Char ch;
  boolean wag;

  for (i = 1; i < nmlngth; i++) {
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
  }
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
      if (ch == '\n')
        ch = ' ';
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
      wagner0[k - 1] = ((long)wagner0[k - 1]) | (1L << ((int)j));
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
}  /* inputmixture */

void printmixture()
{
  /* print out list of parsimony methods */
  long i, k, l;

  fprintf(outfile, "Parsimony methods:\n");
  l = 0;
  k = 1;
  for (i = 1; i <= nmlngth + 3; i++)
    putc(' ', outfile);
  for (i = 1; i <= (chars); i++) {
    newline(i, 55L, nmlngth + 3L);
    l++;
    if (l > bits) {
      l = 1;
      k++;
    }
    if ((unsigned long)l < 32 && ((1L << l) & wagner[k - 1]) != 0)
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
  long i;

  for (i = 1; i < nmlngth; i++) {
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
  }
  for (i = 0; i < (chars); i++) {
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

void printweights()
{
  /* print out the weights of characters */
  long i, j, k;

  fprintf(outfile, "Characters are weighted as follows:\n");
  fprintf(outfile, "       ");
  for (i = 0; i <= 9; i++)
    fprintf(outfile, "%3ld", i);
  fprintf(outfile, "\n      *---------------------------------\n");
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

  for (i = 1; i < nmlngth; i++) {
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
  }
  for (i = 0; i < (chars); i++) {
    anczero0[i] = true;
    ancone0[i] = true;
    do {
      if (eoln(infile)) {
        fscanf(infile, "%*[^\n]");
        getc(infile);
      }
      ch = getc(infile);
      if (ch == '\n')
        ch = ' ';
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
    newline(i, 55L, nmlngth + 3L);
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
  boolean avar, mix;

  if (!firstset) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%ld%ld", &cursp, &curchs);
    if (cursp != spp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4ld\n",
             ith);
      exit(-1);
    }
    chars = curchs;
  }
  extranum = 0;
  avar = false;
  mix = false;
  while (!eoln(infile)) {
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
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
    if (ch == '\n')
      ch = ' ';
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
    printf("ERROR: MIXTURE OPTION CHOSEN IN MENU WITH NO OPTION M IN INPUT\n");
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
  long i, j, l;
  char k;
  Char charstate;
  /* possible states are '0', '1', 'P', 'B', and '?' */
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
	nayme[i - 1][j] = getc(infile);
	if (nayme[i - 1][j] == '\n')
	  nayme[i - 1][j] = ' ';
	if ( eof(infile) || eoln(infile)){
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
/*	  anerror |= eof(infile); */
	  charstate = getc(infile);
	  if (charstate == '\n')
	    charstate = ' ';
	} while (charstate == ' ');
	if (charstate == 'b')	  charstate = 'B';
	if (charstate == 'p')	  charstate = 'P';
	if (charstate != '0' && charstate != '1' && charstate != '?' &&
	    charstate != 'P' && charstate != 'B') {
	  printf("ERROR: BAD CHARACTER STATE: %c ",charstate);
	  printf("AT CHARACTER %5ld OF SPECIES %3ld\n",j,i);
	  exit(-1);
	}
	if (printdata) {
	  newline(j, 55L, nmlngth + 3L);
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

main(argc, argv)
int argc;
Char *argv[];
{  /* Mixed parsimony by uphill search */
char infilename[100],outfilename[100],trfilename[100];
#ifdef MAC
  macsetup("Mix","");
  argv[0] = "Mix";
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
  if (usertree) {
    fsteps = (double **)Malloc(maxuser*sizeof(double *));
    for (j = 1; j <= maxuser; j++)
      fsteps[j - 1] = (double *)Malloc(chars*sizeof(double));
  }
  bestrees = (long **)Malloc(maxtrees*sizeof(long *));
  for (j = 1; j <= maxtrees; j++)
    bestrees[j - 1] = (long *)Malloc(spp*sizeof(long));
  extras = (steptr)Malloc(chars*sizeof(long));
  weight = (steptr)Malloc(chars*sizeof(long));
  threshwt = (double *)Malloc(chars*sizeof(double));
  numsteps = (steptr)Malloc(chars*sizeof(long));
  numszero = (steptr)Malloc(chars*sizeof(long));
  numsone = (steptr)Malloc(chars*sizeof(long));
  guess = (Char *)Malloc(chars*sizeof(Char));
  nayme = (Char **)Malloc(spp*sizeof(Char *));
  for (j = 1; j <= spp; j++)
    nayme[j - 1] = (Char *)Malloc(nmlngth*sizeof(Char));
  enterorder = (long *)Malloc(spp*sizeof(long));
  ancone = (boolean *)Malloc(chars*sizeof(boolean));
  anczero = (boolean *)Malloc(chars*sizeof(boolean));
  ancone0 = (boolean *)Malloc(chars*sizeof(boolean));
  anczero0 = (boolean *)Malloc(chars*sizeof(boolean));
  wagner = (bitptr)Malloc(words*sizeof(long));
  wagner0 = (bitptr)Malloc(words*sizeof(long));
  if (trout)
    openfile(&treefile,TREEFILE,"w",argv[0],trfilename);
  for (ith = 1; ith <= datasets; ith++) {
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
  FClose(outfile);
  FClose(infile);
  FClose(treefile);
#ifdef MAC
  fixmacfile(trfilename);
  fixmacfile(outfilename);
#endif
  exit(0);
}  /* Mixed parsimony by uphill search */


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
