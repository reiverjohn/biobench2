#include "phylip.h"

/* version 3.56c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define nmlngth         10   /* number of characters in species name    */
#define maxtrees        100   /* maximum number of tied trees stored     */
#define maxuser         8   /* maximum number of user-defined trees */

#define ibmpc0          false
#define ansi0           true
#define vt520           false


typedef short *steptr;
typedef short *stepshortptr;
typedef boolean *boolptr;
typedef enum {
  A, C, G, U, O
} bases;
typedef short *baseptr;

typedef struct gbases {
  baseptr base;
  struct gbases *next;
} gbases;

/* nodes will form a binary tree */

typedef struct node {         /* describes a tip species or an ancestor */
  struct node *next, *back;   /* pointers to nodes                      */
  short index;                 /* number of the node                     */
  boolean tip, bottom;        /* present species are tips of tree       */
  baseptr base;               /* the sequence                           */
  stepshortptr numsteps;      /* bookkeeps steps                        */
  short xcoord, ycoord, ymin;  /* used by printree                       */
  short ymax;
} node;

typedef Char **sequence;
typedef node **pointptr;
typedef short longer[6];

/* Local variables for hyptrav, in original pascal: */
struct LOC_hyptrav {
  boolean bottom;
  node *r;
  short *hypset;
  boolean maybe, nonzero;
  short tempset, anc;
} ;


Static node *root, *p;
Static FILE *infile, *outfile, *treefile;
Static short spp, nonodes, chars, endsite, inseed, outgrno, col, datasets,
             ith, i, j, l, jumb, njumble;
/*                                                                    *
 * spp    = number of species                                         *
 * nonodes = number of nodes in tree                                  *
 * Bchars   = number of sites in actual sequences                     *
 * outgrno indicates outgroup                                         */

Static boolean jumble, usertree, weights, trout, outgropt,  printdata,
               progress, treeprint, stepbox, ancseq, firstset, mulsets,
               interleaved, ibmpc, vt52, ansi;
Static steptr weight, oldweight;
Static steptr necsteps, alias, ally, location;
Static pointptr treenode;   /* pointers to all nodes in tree */
Static Char **nayme;   /* names of species */
Static sequence y;
Static longer seed;
Static short *enterorder;
Static boolean *names;
Static Char basechar[32]="ACMGRSVTWYHKDBNO???????????????";
Static gbases *garbage;
Static short **bestrees;
Char ch;
/* Local variables for maketree, propogated globally for C version: */
  short nextree, which, maxwhich;
  double like, maxsteps, bestyet, bestlike, bstlike2;
  boolean lastrearr, recompute;
  double nsteps[maxuser];
  short **fsteps;
  node *there;
  short *place;
  boolptr intree;
  baseptr nothing;
  node *temp, *temp1;



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
gbases **p;
{
  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */
  if (garbage != NULL) {
    *p = garbage;
    garbage = garbage->next;
  } else {
    *p = (gbases *)Malloc(sizeof(gbases));
    (*p)->base = (baseptr)Malloc(endsite*sizeof(short));
  }
  (*p)->next = NULL;
}  /* gnu */


void chuck(p)
gbases *p;
{
  /* collect garbage on p -- put it on front of garbage list */
  p->next = garbage;
  garbage = p;
}  /* chuck */


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
  /* convert ch to upper case -- either ASCII or EBCDIC */
    *ch = (islower (*ch) ? toupper(*ch) : (*ch));
}  /* uppercase */


void getoptions()
{
  /* interactively set options */
  short i, inseed0;
  Char ch;
  boolean  done1;

  fprintf(outfile, "\nDNA compatibility algorithm, version %s\n\n",VERSION);
  putchar('\n');
  jumble = false;
  njumble = 1;
  outgrno = 1;
  outgropt = false;
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
    if (ansi)
      printf("\033[2J\033[H");
    else if (vt52)
      printf("\033E\033H");
    else
      putchar('\n');
    printf("\nDNA compatibility algorithm, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
	   (usertree ? "No, use user trees in input file" : "Yes"));
    if (!usertree) {
      printf("  J   Randomize input order of sequences?");
      if (jumble) {
        printf(
         "  Yes (seed =%8hd,%3hd times)\n", inseed0, njumble);
      }
      else
        printf("  No. Use input order\n");
    }
    printf("  O                        Outgroup root?");
    if (outgropt)
      printf("  Yes, at sequence number%3hd\n", outgrno);
    else
      printf("  No, use as outgroup species%3hd\n", outgrno);
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2hd sets\n", datasets);
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n",
	   (interleaved ? "Yes" : "No, sequential"));
    printf("  0   Terminal type (IBM PC, VT52, ANSI)?  %s\n",
	   ibmpc ? "IBM PC" :
	   ansi  ? "ANSI"   :
	   vt52  ? "VT52"   : "(none)");
    printf("  1    Print out the data at start of run  %s\n",
	   (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
	   (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
	   (treeprint ? "Yes" : "No"));
    printf("  4  Print steps & compatibility at sites  %s\n",
	   (stepbox ? "Yes" : "No"));
    printf("  5  Print sequences at all nodes of tree  %s\n",
	   (ancseq ? "Yes" : "No"));
    printf("  6       Write out trees onto tree file?  %s\n",
	   (trout ? "Yes" : "No"));
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (strchr("JOTUMI1234560",ch) != NULL){
      switch (ch) {
	
      case 'J':
	jumble = !jumble;
	if (jumble) {
	  printf("Random number seed (must be odd)?\n");
	  scanf("%hd%*[^\n]", &inseed);
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
	  scanf("%hd%*[^\n]", &njumble);
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
	    scanf("%hd%*[^\n]", &outgrno);
	    getchar();
	    done1 = (outgrno >= 1 && outgrno <= spp);
	    if (!done1) {
	      printf("BAD OUTGROUP NUMBER: %4hd\n", outgrno);
	      printf("  Must be in range 1 -%2hd\n", spp);
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

void inputnumbers()
{
  /* input the numbers of species and of characters */
  fscanf(infile, "%hd%hd", &spp, &chars);
  if (printdata)
    fprintf(outfile, "%2hd species, %3hd sites\n", spp, chars);
  putc('\n', outfile);
  nonodes = spp * 2 - 1;
}  /* inputnumbers */


void doinit()
{
  /* initializes variables */
  short i;
  node *p, *q;

  inputnumbers();
  getoptions();
  y = (Char **)Malloc(spp*sizeof(Char *));
  for (i = 0; i < spp; i++)
    y[i] = (Char *)Malloc(chars*sizeof(Char));
  treenode = (node **)Malloc((spp*2-1)*sizeof(node *));
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
}  /* doinit */

void inputweights()
{
  /* input the character weights, 0-9 and A-Z for weights
    0 - 35 */
  Char ch;
  short i;

  for (i = 1; i < nmlngth; i++)
    ch = getc(infile);
  for (i = 0; i < chars; i++) {
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

void inputoptions()
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
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4hd\n", ith);
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
  for (i = 0; i < chars; i++)
    weight[i] = 1;
  for (i = 1; i <= extranum; i++) {
    ch = getc(infile);
    uppercase(&ch);
    if (ch == 'W')
      inputweights();
    else{
      printf("ERROR: INCORRECT AUXILIARY OPTIONS LINE WHICH STARTS WITH %c\n",
	     ch);
    exit(-1);}
  }
  if (weights)
    printweights();
}  /* inputoptions */

void inputdata()
{
  /* input the names and sequences for each species */
  short i, j, k, l, basesread, basesnew;
  Char charstate;
  boolean allread, done;

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
	  if ( eof(infile) || eoln(infile)){
              printf("ERROR: END-OF-LINE OR END-OF-FILE IN");
	      printf(" THE MIDDLE OF A SPECIES NAME\n");
	      exit(-1);}
	  nayme[i - 1][j] = getc(infile);
        }
      }
      if (interleaved)
        j = basesread;
      else
        j = 0;
      done = false;
      while (!done && !eof(infile)){
        if (interleaved)
          done = true;
        while (j < chars && !(eoln(infile) || eof(infile))) {
          charstate = getc(infile);
          if (charstate == ' ' || (charstate >= '0' && charstate <= '9'))
            continue;
          uppercase(&charstate);
	  if (strchr("ABCDGHKMNRSTUVWXY?O-.",charstate) == NULL){
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
      if ((interleaved && j != basesnew) || ((!interleaved) && (j != chars))){
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
  if ( !printdata)
    return;
  for (i = 1; i <= ((chars - 1) / 60 + 1); i++) {
    for (j = 1; j <= spp; j++) {
      for (k = 0; k < nmlngth; k++)
        putc(nayme[j - 1][k], outfile);
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

void sitesort()
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

void sitecombine()
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

void sitescrunch()
{
  /* move so one representative of each pattern of
    sites comes first */
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

void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  short i;
  node *p;

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
}  /* makeweights */

void makevalues()
{
  /* set up fractional likelihoods at tips */
  short i, j;
  short ns;
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
    treenode[i]->base = (baseptr)Malloc(endsite*sizeof(short));
  }
  for (i = spp; i < nonodes; i++) {
    p = treenode[i];
    for (j = 1; j <= 3; j++) {
      p->numsteps = (stepshortptr)Malloc(endsite*sizeof(short));
      p->base = (baseptr)Malloc(endsite*sizeof(short));
      p = p->next;
    }
  }
  for (j = 0; j < endsite; j++) {
    for (i = 0; i < spp; i++) {
      switch (y[i][alias[j] - 1]) {

      case 'A':
        ns = 1L << ((short)A);
        break;

      case 'C':
        ns = 1L << ((short)C);
        break;

      case 'G':
        ns = 1L << ((short)G);
        break;

      case 'U':
        ns = 1L << ((short)U);
        break;

      case 'T':
        ns = 1L << ((short)U);
        break;

      case 'M':
        ns = (1L << ((short)A)) | (1L << ((short)C));
        break;

      case 'R':
        ns = (1L << ((short)A)) | (1L << ((short)G));
        break;

      case 'W':
        ns = (1L << ((short)A)) | (1L << ((short)U));
        break;

      case 'S':
        ns = (1L << ((short)C)) | (1L << ((short)G));
        break;

      case 'Y':
        ns = (1L << ((short)C)) | (1L << ((short)U));
        break;

      case 'K':
        ns = (1L << ((short)G)) | (1L << ((short)U));
        break;

      case 'B':
        ns = (1L << ((short)C)) | (1L << ((short)G)) | (1L << ((short)U));
        break;

      case 'D':
        ns = (1L << ((short)A)) | (1L << ((short)G)) | (1L << ((short)U));
        break;

      case 'H':
        ns = (1L << ((short)A)) | (1L << ((short)C)) | (1L << ((short)U));
        break;

      case 'V':
        ns = (1L << ((short)A)) | (1L << ((short)C)) | (1L << ((short)G));
        break;

      case 'N':
        ns = (1L << ((short)A)) | (1L << ((short)C)) | (1L << ((short)G)) |
             (1L << ((short)U));
        break;

      case 'X':
        ns = (1L << ((short)A)) | (1L << ((short)C)) | (1L << ((short)G)) |
             (1L << ((short)U));
        break;

      case '?':
        ns = (1L << ((short)A)) | (1L << ((short)C)) | (1L << ((short)G)) |
             (1L << ((short)U)) | (1L << ((short)O));
        break;

      case 'O':
        ns = 1L << ((short)O);
        break;

      case '-':
        ns = 1L << ((short)O);
        break;
      }
      treenode[i]->base[j] = ns;
      treenode[i]->numsteps[j] = 0;
    }
  }
}  /* makevalues */


void doinput()
{
  /* reads the input data */
  inputoptions();
  inputdata();
  makeweights();
  makevalues();
}  /* doinput */


#define down            2




void mincomp(n)
short n;
{
  /* computes for each site the minimum number of steps
    necessary to accomodate those species already
    in the analysis, adding in species n */
  short i, j, k, l, m;
  bases b;
  short s;
  boolean allowable, deleted;

  intree[n - 1] = true;
  for (i = 0; i < endsite; i++)
    necsteps[i] = 3;
  for (m = 0; m <= 31; m++) {
    s = 0;
    l = -1;
    k = m;
    for (b = A; (short)b <= (short)O; b = (bases)((short)b + 1)) {
      if ((k & 1) == 1) {
        s |= 1L << ((short)b);
        l++;
      }
      k /= 2;
    }
    for (j = 0; j < endsite; j++) {
      allowable = true;
      i = 1;
      while (allowable && i <= spp) {
        if (intree[i - 1] && treenode[i - 1]->base[j] != 0) {
          if ((treenode[i - 1]->base[j] & s) == 0)
            allowable = false;
        }
        i++;
      }
      if (allowable) {
        if (l < necsteps[j])
          necsteps[j] = l;
      }
    }
  }
  for (j = 0; j < endsite; j++) {
    deleted = false;
    for (i = 0; i < spp; i++) {
      if (intree[i] && treenode[i]->base[j] == 0)
        deleted = true;
    }
    if (deleted)
      necsteps[j]++;
  }
  for (i = 0; i < endsite; i++)
    necsteps[i] *= weight[i];
}  /* mincomp */

void fillin(p, left, rt)
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

void preorder(p)
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
  if (!recompute)
    return;
  fillin(newfork, newfork->next->back, newfork->next->next->back);
  preorder(newfork);
  if (newfork != root)
    preorder(newfork->back);
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
    (*item)->back = NULL;
  }
  if (!recompute)
    return;
  preorder(other);
  if (other != root)
    preorder(other->back);
}  /* re_move */

void evaluate(r)
node *r;
{
  /* determines the number of steps needed for a tree. this is
    the minimum number of steps needed to evolve sequences on
    this tree */
  short i, term;
  double sum;

   sum = 0.0;
   for (i = 0; i < endsite; i++) {
    if (r->numsteps[i] == necsteps[i])
      term = weight[i];
    else
      term = 0;
    sum += term;
    if (usertree && which <= maxuser)
      fsteps[which - 1][i] = term;
  }
  if (usertree && which <= maxuser) {
    nsteps[which - 1] = sum;
    if (which == 1) {
      maxwhich = 1;
      maxsteps = sum;
    } else if (sum > maxsteps) {
      maxwhich = which;
      maxsteps = sum;
    }
  }
  like = sum;
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
  fillin(p, p->next->back, p->next->next->back);
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
  outgroup->back->back = root->next->next;
  outgroup->back = root->next;
}  /* reroot */



void savetraverse(p)
node *p;
{
  /* sets BOOLEANs that indicate which way is down */
  p->bottom = true;
  if (p->tip)
    return;
  p->next->bottom = false;
  savetraverse(p->next->back);
  p->next->next->bottom = false;
  savetraverse(p->next->next->back);
}  /* savetraverse */

void savetree()
{
  /* record in place where each species has to be
    added to reconstruct this tree */
  short i, j;
  node *p;
  boolean done;

  reroot(treenode[outgrno - 1]);
  savetraverse(root);
  for (i = 0; i < nonodes; i++)
    place[i] = 0;
  place[root->index - 1] = 1;
  for (i = 1; i <= spp; i++) {
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

void findtree(found,pos)
boolean *found;
short *pos;
{
  /* finds tree given by ARRAY place in ARRAY
    bestrees by binary search */
  short i, lower, upper;
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

void addtree(pos)
short pos;
{
  /* puts tree from ARRAY place in its proper position
    in ARRAY bestrees */
  short i;

  for (i = nextree - 1; i >= pos; i--)
    memcpy(bestrees[i], bestrees[i - 1], nonodes*sizeof(short));
  for (i = 0; i < spp; i++)
    bestrees[pos - 1][i] = place[i];
  nextree++;
}  /* addtree */

void tryadd(p,item,nufork)
node *p,**item,**nufork;
{
  /* temporarily adds one fork and one tip to the tree.
    if the location where they are added yields greater
    "likelihood" than other locations tested up to that
    time, then keeps that location as there */
  short pos;
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
    } else if (like >= bstlike2) {
      recompute = false;
      add(p, *item,*nufork);
      rute = root->next->back;
      savetree();
      reroot(rute);
      if (like > bstlike2) {
        bestlike = bstlike2 = like;
        pos = 1;
        nextree = 1;
        addtree(pos);
      } else {
        pos = 0;
        findtree(&found,&pos);
        if (!found) {
          if (nextree <= maxtrees)
            addtree(pos);
        }
      }
      re_move(item, nufork);
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


void tryrearr(p,success)
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
    if (p != forknode->next->next->back)
      return;
    q = forknode->next;
    forknode->next = forknode->next->next;
    forknode->next->next = q;
    q->next = forknode;
    return;
  }
  recompute = false;
  re_move(&p, &forknode);
  fillin(whereto, whereto->next->back, whereto->next->next->back);
  recompute = true;
  add(whereto, p, forknode);
  *success = true;
  bestyet = like;
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
  boolean success=true;

  while (success) {
    success = false;
    repreorder(*r,&success);
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
      if (ch == '(' || ch == ')' || ch == ';') {
        printf("\nERROR IN USER TREE:");
	printf(" UNMATCHED PARENTHESIS OR MISSING COMMA\n");
	exit(-1);
      } else if (ch == ',')
        done = true;
    } else if (c == ')') {
      if (ch == '(' || ch == ',' ||ch == ';') {
        printf("\nERROR IN USER TREE:");
	printf(" UNMATCHED PARENTHESIS OR NON-BIFURCATED NODE\n");
	exit(-1);
      } else if (ch == ')')
          done = true;
    } else if (c == ';') {
      if (ch != ';') {
        printf("\nERROR IN USER TREE:");
	printf(" UNMATCHED PARENTHESIS OR MISSING SEMICOLON\n");
	exit(-1);
      } else
        done = true;
    }
    if (ch != ')' && done)
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
short *nextnode,*lparens;
boolean *names;
{
  /* recursive procedure adds nodes to user-defined tree */
  node *q;
  short i, n;
  boolean found;
  Char str[nmlngth];

  do {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
  } while (ch == ' ' );
  if (ch == '(' ) {
    if ((*lparens) >= spp - 1) {
      printf("\nERROR IN USER TREE: TOO MANY LEFT PARENTHESES\n");
      exit(-1);
    }
    (*nextnode)++;
    (*lparens)++;
    q = treenode[(*nextnode) - 1];
    addelement(&q->next->back,nextnode,lparens,names);
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
    if (ch == '\n')
      ch = ' ';
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
  if (n <= spp) {
    intree[n - 1] = true;
    return;
  }
  printf("CANNOT FIND SPECIES: ");
  for (i = 0; i < nmlngth; i++)
    putchar(str[i]);
  putchar('\n');
}  /* addelement */

void treeread()
{
  /* read in user-defined tree and set it up */
  short nextnode, lparens,i;

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
      putc(nayme[p->index - 1][j], outfile);
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
  fprintf(outfile, "\n  remember:");
  if (outgropt)
    fprintf(outfile, " (although rooted by outgroup)");
  fprintf(outfile, " this is an unrooted tree!\n\n");
}  /* printree */


void ancestset(a, b, c)
short *a, *b, *c;
{
  /* make the set of ancestral states below nodes
    whose base sets are a and b */
  *c = (*a) & (*b);
  if (*c == 0)
    *c = (*a) | (*b);
}  /* ancestset */



void hyprint(b1, b2, hyptrav_vars)
short b1, b2;
struct LOC_hyptrav *hyptrav_vars;
{
  /* print out states in sites b1 through b2 at node */
  short i, j, k, n;
  boolean dot;
  bases b;

  if (hyptrav_vars->bottom) {
    if (!outgropt)
      fprintf(outfile, "       ");
    else
      fprintf(outfile, "root   ");
  } else
    fprintf(outfile, "%4hd   ", hyptrav_vars->r->back->index - spp);
  if (hyptrav_vars->r->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[hyptrav_vars->r->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4hd      ", hyptrav_vars->r->index - spp);
  if (hyptrav_vars->bottom)
    fprintf(outfile, "          ");
  else if (hyptrav_vars->nonzero)
    fprintf(outfile, "   yes    ");
  else if (hyptrav_vars->maybe)
    fprintf(outfile, "  maybe   ");
  else
    fprintf(outfile, "   no     ");
  for (i = b1; i <= b2; i++) {
    j = location[ally[i - 1] - 1];
    hyptrav_vars->tempset = hyptrav_vars->r->base[j - 1];
    hyptrav_vars->anc = hyptrav_vars->hypset[j - 1];
    if (!hyptrav_vars->bottom)
      hyptrav_vars->anc
           = treenode[hyptrav_vars->r->back->index - 1]->base[j - 1];
    dot
      = (hyptrav_vars->tempset == hyptrav_vars->anc && !hyptrav_vars->bottom);
    if (dot)
      putc('.', outfile);
    else if (hyptrav_vars->tempset == 1L << ((short)A))
      putc('A', outfile);
    else if (hyptrav_vars->tempset == 1L << ((short)C))
      putc('C', outfile);
    else if (hyptrav_vars->tempset == 1L << ((short)G))
      putc('G', outfile);
    else if (hyptrav_vars->tempset == 1L << ((short)U))
      putc('T', outfile);
    else if (hyptrav_vars->tempset == 1L << ((short)O))
      putc('-', outfile);
    else {
      k = 1;
      n = 0;
      for (b = A; (short)b <= (short)O; b = (bases)((short)b + 1)) {
        if (((1L << ((short)b)) & hyptrav_vars->tempset) != 0)
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

void hyptrav(r_, hypset_, b1, b2, bottom_)
node *r_;
short *hypset_;
short b1, b2;
boolean bottom_;
{
  /*  compute, print out states at one interior node */
  struct LOC_hyptrav Vars;
  short i, j;
  short left, rt;
  gbases *temparray, *ancset;

  Vars.bottom = bottom_;
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
    if (!Vars.bottom)
      Vars.anc = treenode[Vars.r->back->index - 1]->base[j - 1];
    Vars.nonzero = (Vars.nonzero || (Vars.r->base[j - 1] & Vars.anc) == 0);
    Vars.maybe = (Vars.maybe || Vars.r->base[j - 1] != Vars.anc);
  }
  hyprint(b1, b2, &Vars);
  Vars.bottom = false;
  if (!Vars.r->tip) {
    memcpy(temparray->base, Vars.r->next->back->base, endsite*sizeof(short));
    for (i = b1 - 1; i < b2; i++) {
      j = location[ally[i] - 1];
      ancestset(&Vars.hypset[j - 1], &Vars.r->next->next->back->base[j - 1],
                &ancset->base[j - 1]);
    }
    hyptrav(Vars.r->next->back, ancset->base, b1, b2, Vars.bottom);
    for (i = b1 - 1; i < b2; i++) {
      j = location[ally[i] - 1];
      ancestset(&Vars.hypset[j - 1], &temparray->base[j - 1],
                &ancset->base[j - 1]);
    }
    hyptrav(Vars.r->next->next->back, ancset->base, b1, b2, Vars.bottom);
  }
  chuck(temparray);
  chuck(ancset);
}  /* hyptrav */

void hypstates()
{
  /* fill in and describe states at interior nodes */
  short i, n;

  fprintf(outfile, "\nFrom    To     Any Steps?    State at upper node\n");
  fprintf(outfile, "                            ");
  fprintf(outfile, " ( . means same as in the node below it on tree)\n");
  nothing = (baseptr)Malloc(endsite*sizeof(short));
  for (i = 0; i < endsite; i++)
    nothing[i] = 0;
  for (i = 1; i <= ((chars - 1) / 40 + 1); i++) {
    putc('\n', outfile);
    n = i * 40;
    if (n > chars)
      n = chars;
    hyptrav(root, nothing, i * 40 - 39, n, true);
  }
  free(nothing);
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
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = nayme[p->index - 1][i];
      if (c == ' ')
        c = '#';
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
    each site and table of compatibilities */
  short i, j, k, l;

  if (treeprint) {
    fprintf(outfile, "\ntotal number of compatible sites is ");
    fprintf(outfile, "%10.1f\n", like);
  }
  if (stepbox) {
    putc('\n', outfile);
    if (weights)
      fprintf(outfile, "weighted ");
    fprintf(outfile, "\n steps in each site:\n");
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
    fprintf(outfile,
            "\n compatibility (Y or N) of each site with this tree:\n\n");
    fprintf(outfile, "      ");
    for (i = 0; i <= 9; i++)
      fprintf(outfile, "%hd", i);
    fprintf(outfile, "\n     *----------\n");
    for (i = 0; i <= (chars / 10); i++) {
      putc(' ', outfile);
      fprintf(outfile, "%3hd !", i * 10);
      for (j = 0; j <= 9; j++) {
        k = i * 10 + j;
        if (k > 0 && k <= chars) {
          if (root->numsteps[location[ally[k - 1] - 1] - 1] ==
              necsteps[location[ally[k - 1] - 1] - 1]) {
            if (oldweight[k - 1] > 0)
              putc('Y', outfile);
            else
              putc('y', outfile);
          } else {
            if (oldweight[k - 1] > 0)
              putc('N', outfile);
            else
              putc('n', outfile);
          }
        } else
          putc(' ', outfile);
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



void maketree()
{
  /* constructs a binary tree from the pointers in treenode.
    adds each node at location which yields highest "likelihood"
    then rearranges the tree for greatest "likelihood" */
  short i, j, k, numtrees, num;
  double gotlike, wt, sumw, sum, sum2, sd;
  node *item, *nufork, *dummy;
  short TEMP2;

  temp = (node *)Malloc(sizeof(node));
  temp->numsteps = (stepshortptr)Malloc(endsite*sizeof(short));
  temp->base = (baseptr)Malloc(endsite*sizeof(short));
  temp1 = (node *)Malloc(sizeof(node));
  temp1->numsteps = (stepshortptr)Malloc(endsite*sizeof(short));
  temp1->base = (baseptr)Malloc(endsite*sizeof(short));
  if (!usertree) {
    recompute = true;
    for (i = 0; i < spp; i++)
      intree[i] = false;
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble) {
      for (i = 0; i < spp; i++) {
        j = (short)(randum(seed) * spp) + 1;
        k = enterorder[j - 1];
        enterorder[j - 1] = enterorder[i];
        enterorder[i] = k;
      }
    }
    root = treenode[enterorder[0] - 1];
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
    intree[0] = true;
    intree[1] = true;
    lastrearr = false;
    for (i = 3; i <= spp; i++) {
      mincomp(i);
      bestyet = -10.0 * spp * chars;
      item = treenode[enterorder[i - 1] - 1];
      nufork = treenode[spp + i - 2];
      there = root;
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
          for (j = 0; j < nonodes; j++) {
            bestyet = -10.0 * spp * chars;
            item = treenode[j];
            there = root;
            if (item != root) {
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
          fprintf(outfile, "%6hd trees in all found\n", nextree - 1);
      }
      if (nextree > maxtrees + 1) {
        if (treeprint)
          fprintf(outfile, "here are the first%4hd of them\n", (short)maxtrees);
        nextree = maxtrees + 1;
      }
      if (treeprint)
        putc('\n', outfile);
      recompute = false;
      for (i = 0; i <= (nextree - 2); i++) {
        root = treenode[0];
        add(treenode[0], treenode[1], treenode[spp]);
        for (j = 3; j <= spp; j++)
          add(treenode[bestrees[i][j - 1] - 1], treenode[j - 1],
            treenode[spp + j - 2]);
        reroot(treenode[outgrno - 1]);
        postorder(root);
        evaluate(root);
        printree();
        describe();
        for (j = 1; j < spp; j++)
          re_move(&treenode[j], &dummy);
      }
    }
  } else {
    fscanf(infile, "%hd%*[^\n]", &numtrees);
    getc(infile);
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n\n\n");
    }
    fsteps = (short **)Malloc(maxuser*sizeof(short *));
    for (j = 1; j <= maxuser; j++)
      fsteps[j - 1] = (short *)Malloc(endsite*sizeof(short));
    maxsteps = 0.0;
    which = 1;
    while (which <= numtrees) {
      for (j = 0; j < spp; j++)
        intree[j] = false;
      treeread();
      j = 1;
      while (!intree[j - 1])
	j++;
      mincomp(j);
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
      fprintf(outfile, "Tree    Compatible  Difference  Its S.D.");
      fprintf(outfile, "   Significantly worse?\n\n");
      if (numtrees > maxuser)
        num = maxuser;
      else
        num = numtrees;
      for (which = 1; which <= num; which++) {
        fprintf(outfile, "%3hd%15.1f", which, nsteps[which - 1]);
        if (maxwhich == which)
          fprintf(outfile, "  <------ best\n");
        else {
          sumw = 0.0;
          sum = 0.0;
          sum2 = 0.0;
          for (j = 0; j < endsite; j++) {
            if (weight[j] > 0) {
              wt = weight[j];
              sumw += wt;
              sum += fsteps[maxwhich - 1][j] - fsteps[which - 1][j];
              TEMP2 = fsteps[which - 1][j] - fsteps[maxwhich - 1][j];
              sum2 += TEMP2 * TEMP2 / wt;
            }
          }
          sd = sqrt(sumw / (sumw - 1.0) * (sum2 - (sum /sumw) * (sum / sumw)));
          fprintf(outfile, "%9.1f%12.4f",
                  maxsteps - nsteps[which - 1], sd);
          if (sum > 1.95996 * sd)
            fprintf(outfile, "           Yes\n");
          else
            fprintf(outfile, "            No\n");
        }
      }
      fprintf(outfile, "\n\n");
    }
    for (j = 1; j <= maxuser; j++)
      free(fsteps[j - 1]);
    free(fsteps);
  }
  if (jumb == njumble) {
    if (progress) {
      printf("Output written to output file\n\n");
      if (trout)
        printf("Trees also written onto file\n\n");
    }
    free(temp->numsteps);
    free(temp->base);
    free(temp);
    free(temp1->numsteps);
    free(temp1->base);
    free(temp1);
  }
}  /* maketree */


main(argc, argv)
int argc;
Char *argv[];
{  /* DNA compatibility by uphill search */
char infilename[100],outfilename[100],trfilename[100];
  /* reads in spp, chars, and the data. Then calls maketree to
    construct the tree */
#ifdef MAC
  macsetup("Dnacomp","");
  argv[0] = "Dnacomp";
#endif
  openfile(&infile,INFILE,"r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);
  mulsets = false;
  garbage = NULL;
  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  datasets = 1;
  firstset = true;
  doinit();
  bestrees = (short **) Malloc(maxtrees*sizeof(short *));
  for (j = 1; j <= maxtrees; j++)
    bestrees[j - 1] = (short *)Malloc(nonodes*sizeof(short));
  nayme = (Char **)Malloc(spp*sizeof(Char *));
  for (j = 1; j <= spp; j++)
    nayme[j - 1] = (Char *)Malloc(nmlngth*sizeof(Char));
  weight = (steptr)Malloc(chars*sizeof(short));
  oldweight = (steptr)Malloc(chars*sizeof(short));
  enterorder = (short *)Malloc(spp*sizeof(short));
  necsteps = (steptr)Malloc(chars*sizeof(short));
  alias = (steptr)Malloc(chars*sizeof(short));
  ally = (steptr)Malloc(chars*sizeof(short));
  location = (steptr)Malloc(chars*sizeof(short));
  place = (short *)Malloc((2*spp-1)*sizeof(short));
  intree = (boolptr)Malloc(chars*sizeof(boolean));
  names = (boolean *)Malloc(spp*sizeof(boolean));
  if (trout)
    openfile(&treefile,TREEFILE,"w",argv[0],trfilename);
  for (ith = 1; ith <= datasets; ith++) {
    doinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1) {
      fprintf(outfile, "Data set # %hd:\n\n", ith);
      if (progress)
        printf("Data set # %hd:\n\n", ith);
    }
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
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
}  /* DNA compatibility by uphill search */

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

