#include "phylip.h"

/* version 3.572c. (c) Copyright 1995 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define nmlngth         10   /* number of characters in species name    */

#define ibmpc0          false
#define ansi0           true
#define vt520           false
#define downn           2
#define overr           4

typedef short *steptr;

typedef Char naym[nmlngth];

typedef enum {
  horiz, vert, up, over, upcorner, downcorner, one, zero, quest, polym
} chartype;
typedef long *bitptr;
/* nodes will form a binary tree */

typedef struct node {         /* describes a tip species or an ancestor */
  struct node *next, *back;   /* pointers to nodes                      */
  short index;                /* number of the node                     */
  boolean tip;                /* present species are tips of tree       */
  bitptr stateone, statezero;/* see in PROCEDURE fillin               */
  short xcoord, ycoord, ymin, ymax;/* used by printree                  */
  /* stores state for printing */
  Char state;
} node;
typedef enum {  rearr, flipp, reroott, none } rearrtype;

typedef node **pointptr;
typedef enum {
  arb, use, spec
} howtree;

char infilename[100],trfilename[100];
Static node *root;
Static FILE *infile, *treefile;
Static short spp, nonodes, chars, words, outgrno, screenlines, col;
/* spp = number of species
  nonodes = number of nodes in tree
  chars = number of binary characters
  words = number of words needed to represent characters of one organism
  outgrno indicates outgroup */
Static boolean weights, thresh, ancvar, questions, dollo, factors,
	       waswritten, ibmpc, vt52, ansi;
Static steptr extras, weight;
Static boolean *ancone, *anczero, *ancone0, *anczero0;
Static Char *factor;
Static pointptr treenode;   /* pointers to all nodes in tree */
Static naym *nayme;         /* names of species              */
Static double threshold;
Static double *threshwt;
Static Char cha[10];
Static boolean reversed[10];
Static boolean graphic[10];
Static howtree how;
char progname[25];
char ch;
Static short bits = 8*sizeof(long) - 1;

/* Local variables for treeconstruct, propagated globally for c version: */
short dispchar, dispword, dispbit, atwhat, what, fromwhere, towhere,
  oldoutgrno, compatible;
double like, bestyet, gotlike;
Char *guess;
boolean display, newtree, changed, subtree, written, oldwritten, restoring,
  wasleft, oldleft, earlytree;
boolean *intree;
steptr numsteps;
long fullset;
bitptr zeroanc, oneanc;
node *nuroot;
rearrtype lastop;
steptr numsone, numszero;


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

boolean digit(ch)
Char ch;
{
  /* tests whether ch is a digit */
  return isdigit(ch);
}  /* digit */


void uppercase(ch)
Char *ch;
{
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

void inpnum(n, success)
short *n;
boolean *success;
{
  int fields;
  char line[100];
  gets(line);
  *n = atof(line);
  fields = sscanf(line,"%hd",n);
  *success = (fields == 1);

}  /* inpnum */



void inputnumbers()
{
  /* input the numbers of species and of characters */
  fscanf(infile, "%hd%hd", &spp, &chars);
  printf("%3hd Species, %3hd Characters\n\n", spp, chars);
  words = chars / bits + 1;
  nonodes = spp * 2 - 1;
}  /* inputnumbers */

void getoptions()
{
  /* interactively set options */
  Char ch;
  boolean done, done1, gotopt;

  how = arb;
  thresh = false;
  threshold = spp;
  weights = false;
  ancvar = false;
  factors = false;
  dollo = true;
  do {
    printf((ansi || ibmpc) ? "\033[2J\033[H" :
	   vt52 ? "\033E\033H"    : "\n");
    printf("\nInteractive Dollo or polymorphism parsimony,");
    printf(" version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  P                        Parsimony method?");
    printf("  %s\n",(dollo ? "Dollo" : "Polymorphism"));
    printf("  T                 Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  A      Use ancestral states in input file?");
    printf("  %s\n",(ancvar ? "Yes" : "No"));
    printf("  U Initial tree (arbitrary, user, specify)?");
    printf("  %s\n",(how == arb) ? "Arbitrary" :
	            (how == use) ? "User tree from tree file" :
	                           "Tree you specify");
    printf("  0      Graphics type (IBM PC, VT52, ANSI)?  %s\n",
	   ibmpc ? "IBM PC" :
	   ansi  ? "ANSI"   :
	   vt52  ? "VT52"   :
	          "(none)");
    printf("  L               Number of lines on screen?%4hd",screenlines);
    do {
      printf(
	"\n\nAre these settings correct? (type Y or the letter for one to change)\n");
      scanf("%c%*[^\n]", &ch);
      getchar();
      uppercase(&ch);
      done = (ch == 'Y');
      gotopt = (strchr("TPULA0",ch) != NULL) ? true : false;
      if (gotopt) {
	switch (ch) {

	case 'A':
	  ancvar = !ancvar;
	  break;

	case 'P':
	  dollo = !dollo;
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
	  if (how == arb)
	    how = use;
	  else if (how == use)
	    how = spec;
	  else
	    how = arb;
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

	case 'L':
	  do {
	    printf("Number of lines on screen?\n");
	    scanf("%hd", &screenlines);
	  } while (screenlines <= 12);
	  break;
	}
      }
      if (!(gotopt || done))
	printf("Not a possible option!\n");
    } while (!(gotopt || done));
  } while (!done);
}  /* getoptions */

void inputfactors()
{
  /* reads the factor symbols */
  short i;
  Char ch;

  for (i = 1; i < nmlngth; i++)
    ch = getc(infile);
  for (i = 0; i < (chars); i++) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    factor[i] = getc(infile);
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  factors = true;
}  /* inputfactors */


void printfactors()
{
  /* print out list of factor symbols */
  short i;
  printf("Factors:");
  for (i = 1; i <= nmlngth - 5; i++)
    putchar(' ');
  for (i = 1; i <= (chars); i++) {
    newline(i, 55, (int)(nmlngth + 3));
    putchar(factor[i - 1]);
    if (i % 5 == 0)
      putchar(' ');
  }
  printf("\n\n");
}  /* printfactors */

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
    if (digit(ch))
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

  printf("   Characters are weighted as follows:\n");
  printf("        ");
  for (i = 0; i <= 9; i++)
    printf("%3hd", i);
  printf("\n     *---------------------------------\n");
  for (j = 0; j <= (chars / 10); j++) {
    printf("%5hd!  ", j * 10);
    for (i = 0; i <= 9; i++) {
      k = j * 10 + i;
      if (k > 0 && k <= chars)
	printf("%3hd", weight[k - 1]);
      else
	printf("   ");
    }
    putchar('\n');
  }
  putchar('\n');
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
    if (strchr("10PB?",ch)){
      anczero0[i] = (ch == '1') ? false : anczero0[i];
      ancone0[i] = (ch == '0') ? false : ancone0[i];
    } else {
      printf("BAD ANCESTOR STATE: %cAT CHARACTER %4hd\n", ch, i + 1);
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

  printf("Ancestral states:\n");
  for (i = 1; i <= nmlngth + 3; i++)
    putchar(' ');
  for (i = 1; i <= (chars); i++) {
    newline(i, 55, (int)(nmlngth + 3));
    if (ancone[i - 1] && anczero[i - 1])
      putchar('?');
    else if (ancone[i - 1])
      putchar('1');
    else
      putchar('0');
    if (i % 5 == 0)
      putchar(' ');
  }
  printf("\n\n");
}  /* printancestor */

void inputoptions()
{
  /* input the information on the options */
  Char ch;
  short extranum, i;
  boolean avar;

  extranum = 0;
  avar = false;
  while (!eoln(infile)) {
    ch = getc(infile);
    uppercase(&ch);
    if (ch == 'A' || ch == 'F' || ch == 'W')
      extranum++;
    else if (ch != ' ') {
      printf("BAD OPTION CHARACTER: %c\n\n", ch);
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
    if (ch != 'A' && ch != 'F' && ch != 'W'){
      printf("ERROR: INCORRECT AUXILIARY OPTIONS LINE WHICH STARTS WITH %c\n",
	     ch);
      exit(-1);}
    if (ch == 'A') {
      avar = true;
      if (!ancvar) {
	printf("ERROR: ANCESTOR OPTION NOT CHOSEN IN MENU WITH OPTION %c IN INPUT\n",
	       ch);
	exit(-1);
      } else
	inputancestors();
    }
    if (ch == 'F')
      inputfactors();
    if (ch == 'W')
      inputweights();
  }
  if (ancvar && !avar) {
    puts("ERROR: ANCESTOR OPTION CHOSEN IN MENU WITH NO OPTION A IN INPUT");
    exit(-1);
  }
  putchar('\n');
  if (weights)
    printweights();
  if (factors)
    printfactors();
  for (i = 0; i < (chars); i++) {
    if (!ancvar) {
      anczero[i] = true;
      ancone[i] = false;
    } else {
      anczero[i] = anczero0[i];
      ancone[i] = ancone0[i];
    }
  }
  if (ancvar && avar)
    printancestors();
  if (!thresh)
    threshold = spp;
  questions = false;
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
  node *p, *q;

  for (i = 0; i < (chars); i++)
    extras[i] = 0;
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
      for (j = 0; j < (words); j++) {
	treenode[i - 1]->stateone[j] = 0;
	treenode[i - 1]->statezero[j] = 0;
      }
      for (j = 0; j < (chars); j++) {
	k = j % bits + 1;
	l = j / bits + 1;
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
	  printf("INPUT ERROR -- BAD CHARACTER STATE: %c ",charstate);
	  printf("AT CHARACTER %5hd OF SPECIES %3hd\n",j+1,i);
	  exit(-1);
	}
	if (charstate == '1')
	  treenode[i - 1]->stateone[l - 1] =
	    ((long)treenode[i - 1]->stateone[l - 1]) | (1L << k);
	if (charstate == '0')
	  treenode[i - 1]->statezero[l - 1] =
	    ((long)treenode[i - 1]->statezero[l - 1]) | (1L << k);
	if (charstate == 'P' || charstate == 'B') {
	  if (dollo)
	    extras[j] += weight[j];
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
    }
  }
  root = NULL;
  printf("\n\n");
}  /* inputdata */


void doinput()
{
  /* reads the input data */
  inputnumbers();
  printf("\nReading input file ...\n\n");
  getoptions();
  nayme = (naym *)Malloc(spp*sizeof(naym));
  intree = (boolean *)Malloc(nonodes*sizeof(boolean));
  extras = (steptr)Malloc(chars*sizeof(short));
  weight = (steptr)Malloc(chars*sizeof(short));
  numszero = (steptr)Malloc(chars*sizeof(short));
  numsone = (steptr)Malloc(chars*sizeof(short));
  threshwt = (double *)Malloc(chars*sizeof(double));
  factor = (Char *)Malloc(chars*sizeof(Char));
  ancone = (boolean *)Malloc(chars*sizeof(boolean));
  anczero = (boolean *)Malloc(chars*sizeof(boolean));
  ancone0 = (boolean *)Malloc(chars*sizeof(boolean));
  anczero0 = (boolean *)Malloc(chars*sizeof(boolean));
  zeroanc = (bitptr)Malloc(words*sizeof(long));
  oneanc = (bitptr)Malloc(words*sizeof(long));
  inputoptions();
  inputdata();
}  /* doinput */


void configure()
{
  /* configure to machine -- set up special characters */
  chartype a;

  for (a = horiz; (long)a <= (long)polym; a = (chartype)((long)a + 1))
    reversed[(long)a] = false;
  for (a = horiz; (long)a <= (long)polym; a = (chartype)((long)a + 1))
    graphic[(long)a] = false;
  if (ibmpc) {
    cha[(long)horiz] = 205;
    graphic[(long)horiz] = true;
    cha[(long)vert] = 186;
    graphic[(long)vert] = true;
    cha[(long)up] = 186;
    graphic[(long)up] = true;
    cha[(long)over] = 205;
    graphic[(long)over] = true;
    cha[(long)one] = 219;
    reversed[(long)one] = true;
    cha[(long)zero] = 176;
    graphic[(long)zero] = true;
    cha[(long)quest] = 178;   /* or try CHR(177) */
    cha[(long)polym] = '\001';
    reversed[(long)polym] = true;
    cha[(long)upcorner] = 200;
    graphic[(long)upcorner] = true;
    cha[(long)downcorner] = 201;
    graphic[(long)downcorner] = true;
    graphic[(long)quest] = true;
    return;
  }
  if (vt52) {
    cha[(long)one] = ' ';
    reversed[(long)one] = true;
    cha[(long)horiz] = cha[(long)one];
    reversed[(long)horiz] = true;
    cha[(long)vert] = cha[(long)one];
    reversed[(long)vert] = true;
    cha[(long)up] = '`';
    graphic[(long)up] = true;
    cha[(long)over] = 'a';
    graphic[(long)over] = true;
    cha[(long)zero] = 'i';
    graphic[(long)zero] = true;
    cha[(long)quest] = '?';
    reversed[(long)quest] = true;
    cha[(long)polym] = '%';
    reversed[(long)polym] = true;
    cha[(long)upcorner] = 'e';
    graphic[(long)upcorner] = true;
    cha[(long)downcorner] = 'f';
    graphic[(long)downcorner] = true;
    return;
  }
  if (ansi) {
    cha[(long)one] = ' ';
    reversed[(long)one] = true;
    cha[(long)horiz] = cha[(long)one];
    reversed[(long)horiz] = true;
    cha[(long)vert] = cha[(long)one];
    reversed[(long)vert] = true;
    cha[(long)up] = 'x';
    graphic[(long)up] = true;
    cha[(long)over] = 'q';
    graphic[(long)over] = true;
    cha[(long)zero] = 'a';
    graphic[(long)zero] = true;
    reversed[(long)zero] = true;
    cha[(long)quest] = '?';
    cha[(long)quest] = '?';
    reversed[(long)quest] = true;
    cha[(long)polym] = '%';
    reversed[(long)polym] = true;
    cha[(long)upcorner] = 'm';
    graphic[(long)upcorner] = true;
    cha[(long)downcorner] = 'l';
    graphic[(long)downcorner] = true;
    return;
  }
  cha[(long)horiz] = '=';
  cha[(long)vert] = ' ';
  cha[(long)up] = '!';
  cha[(long)over] = '-';
  cha[(long)one] = '*';
  cha[(long)zero] = '=';
  cha[(long)quest] = '.';
  cha[(long)polym] = '%';
  cha[(long)upcorner] = '`';
  cha[(long)downcorner] = ',';
}  /* configure */


void pregraph()
{
  /* turn on graphic characters */
  if (vt52)
    printf("\033F");
  if (ansi)
    printf("\033(0");
}  /* pregraph */

void prereverse()
{
  /* turn on reverse video */
  if (vt52)
    printf("\033p");
  if (ansi)
    printf("\033[7m");
}  /* prereverse */


void prefix(a)
chartype a;
{
  /* give prefix appropriate for this character */
  if (reversed[(long)a])
    prereverse();
  if (graphic[(long)a])
    pregraph();
}  /* prefix */


void postgraph()
{
  /* turn off graphic characters */
  if (vt52)
    printf("\033G");
  if (ansi)
    printf("\033(B");
}  /* postgraph */

void postreverse()
{
  /* turn off reverse video */
  if (vt52)
    printf("\033q");
  if (ansi)
    printf("\033[0m");
}  /* postreverse */


void postfix(a)
chartype a;
{
  /* give postfix appropriate for this character */
  if (reversed[(long)a])
    postreverse();
  if (graphic[(long)a])
    postgraph();
}  /* postfix */


void makechar(a)
chartype a;
{
  /* print out a character with appropriate prefix or postfix */
  prefix(a);
  putchar(cha[(long)a]);
  postfix(a);
}  /* makechar */



void add(below, newtip, newfork)
node *below, *newtip, *newfork;
{
  /* inserts the nodes newfork and its left descendant, newtip,
    to the tree.  below becomes newfork's right descendant */
  boolean putleft;
  node *leftdesc, *rtdesc;

  if (below != treenode[below->index - 1])
    below = treenode[below->index - 1];
  if (below->back != NULL)
    below->back->back = newfork;
  newfork->back = below->back;
  putleft = true;
  if (restoring)
    putleft =wasleft;
  if (putleft) {
    leftdesc = newtip;
    rtdesc = below;
  } else {
    leftdesc = below;
    rtdesc = newtip;
  }
  rtdesc->back = newfork->next->next;
  newfork->next->next->back = rtdesc;
  newfork->next->back = leftdesc;
  leftdesc->back = newfork->next;
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
  if (*item == (*fork)->next->back) {
    if (root == *fork)
      root = (*fork)->next->next->back;
    wasleft = true;
  } else {
    if (root == *fork)
      root = (*fork)->next->back;
    wasleft = false;
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
    in a most parnsimonious reconstruction, based on the
    information at or above this node.  Note that this
    state assignment may change based on information further
    down the tree.  If a character is in both sets it is in
    state "P".  If in neither, it is "?". */
  short i;

  for (i = 0; i < (words); i++) {
    p->stateone[i] = p->next->back->stateone[i] | p->next->next->back->stateone[i];
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
      s1 = oneanc[i];
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
      z1 &= ~(oneanc[i] & z0);
      z0 = temp;
    }
    temp = fullset & (~z0) & (~z1);
    p->statezero[i] = z0 | (temp & s0 & (~s1));
    p->stateone[i] = z1 | (temp & s1 & (~s0));
  }
}  /* correct */

/* Local variables for evaluate: */


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


void count(p)
node *p;
{
  /* counts the number of steps in a fork of the tree.
    The program spends much of its time in this PROCEDURE */
  short i, j, l;
  bitptr steps;

  steps = (bitptr)Malloc(words*sizeof(long));
  if (dollo) {
    for (i = 0; i < (words); i++)
      steps[i] = (treenode[p->back->index - 1]->stateone[i] &
		  p->statezero[i] & zeroanc[i]) |
		 (treenode[p->back->index - 1]->statezero[i] &
		  p->stateone[i] & oneanc[i]);
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
  free(steps);
}  /* count */

void preorder(p)
node *p;
{
  /* go back up tree setting up and counting interior node
    states */

  if (!p->tip) {
    correct(p);
    preorder(p->next->back);
    preorder(p->next->next->back);
  }
  if (p->back != NULL)
    count(p);
}  /* preorder */

void evaluate(r)
node *r;
{
  /* Determines the number of losses or polymorphisms needed for a tree.
     This is the minimum number needed to evolve chars on this tree */
  short i, stepnum, smaller;
  double sum;
  boolean nextcompat, thiscompat, done;

  sum = 0.0;
  for (i = 0; i < (chars); i++) {
    numszero[i] = 0;
    numsone[i] = 0;
  }
  for (i = 0; i < (words); i++) {
    zeroanc[i] = fullset;
    oneanc[i] = 0;
  }
  compatible = 0;
  nextcompat = true;
  postorder(r);
  preorder(r);
  for (i = 0; i < (words); i++) {
    zeroanc[i] = 0;
    oneanc[i] = fullset;
  }
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
    thiscompat = (stepnum <= weight[i]);
    if (factors) {
      done = (i + 1 == chars);
      if (!done)
	done = (factor[i + 1] != factor[i]);
      nextcompat = (nextcompat && thiscompat);
      if (done) {
	if (nextcompat)
	  compatible += weight[i];
	nextcompat = true;
      }
    } else if (thiscompat)
      compatible += weight[i];
    guess[i] = '?';
    if (!ancone[i] ||
	(anczero[i] && numszero[i] < numsone[i]))
      guess[i] = '0';
    else if (!anczero[i] ||
	     (ancone[i] && numsone[i] < numszero[i]))
      guess[i] = '1';
  }
  like = -sum;
}  /* evaluate */

void reroot(outgroup)
node *outgroup;
{
  /* reorients tree, putting outgroup in desired position. */
  node *p, *q, *newbottom, *oldbottom;
  boolean onleft;

  if (outgroup->back->index == root->index)
    return;
  newbottom = outgroup->back;
  p = treenode[newbottom->index - 1]->back;
  while (p->index != root->index) {
    oldbottom = treenode[p->index - 1];
    treenode[p->index - 1] = p;
    p = oldbottom->back;
  }
  onleft = (p == root->next);
  if (restoring)
    if (!onleft && wasleft){
      p = root->next->next;
      q = root->next;
    } else {
      p = root->next;
      q = root->next->next;
    }
  else {
    if (onleft)
      oldoutgrno = root->next->next->back->index;
    else
      oldoutgrno = root->next->back->index;
    wasleft = onleft;
    p = root->next;
    q = root->next->next;
  }
  p->back->back = q->back;
  q->back->back = p->back;
  p->back = outgroup;
  q->back = outgroup->back;
  if (restoring) {
    if (!onleft && wasleft) {
      outgroup->back->back = root->next;
      outgroup->back = root->next->next;
    } else {
      outgroup->back->back = root->next->next;
      outgroup->back = root->next;
    }
  } else {
    outgroup->back->back = root->next->next;
    outgroup->back = root->next;
  }
  treenode[newbottom->index - 1] = newbottom;
}  /* reroot */

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

void hyptrav(r)
node *r;
{
  /* compute states at interior nodes for one character */
  if (!r->tip)
    correct(r);
  if (((1L << dispbit) & r->stateone[dispword - 1]) != 0) {
    if (((1L << dispbit) & r->statezero[dispword - 1]) != 0) {
      if (dollo)
	r->state = '?';
      else
	r->state = 'P';
    } else
      r->state = '1';
  } else {
    if (((1L << dispbit) & r->statezero[dispword - 1]) != 0)
      r->state = '0';
    else
      r->state = '?';
  }
  if (!r->tip) {
    hyptrav(r->next->back);
    hyptrav(r->next->next->back);
  }
}  /* hyptrav */

void hypstates()
{
  /* fill in and describe states at interior nodes */
  short i, j, k;

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
  }
  filltrav(root);
  hyptrav(root);
}  /* hypstates */

/* Local variables for printree: */


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
    (*tipy) += downn;
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
  p->xcoord = (last->ymax - first->ymin) * 3 / 2;
  p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */

void drawline(i, lastline)
short i, lastline;
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q, *r, *first, *last;
  short n, j;
  boolean extra, done;
  Char s, cc;
  chartype c, d;

  p = nuroot;
  q = nuroot;
  extra = false;
  if (i == p->ycoord && (p == root || subtree)) {
    c = over;
    if (p == root)
      cc = guess[dispchar - 1];
    else
      cc = p->state;
    if (display) {
      switch (cc) {

      case '1':
	c = one;
	break;

      case '0':
	c = zero;
	break;

      case '?':
	c = quest;
	break;

      case 'P':
	c = polym;
	break;
      }
    }
    if (p->index >= 10) {
      makechar(c);
      printf("%2hd", p->index);
    } else {
      makechar(c);
      makechar(c);
      printf("%hd", p->index);
    }
    extra = true;
  } else
    printf("  ");
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
    n = p->xcoord - q->xcoord;
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if (q->ycoord == i && !done) {
      if (q->ycoord > p->ycoord)
	d = upcorner;
      else
	d = downcorner;
      c = over;
      s = q->state;
      if (s == 'P' && p->state != 'P')
	s = p->state;
      if (display) {
	switch (s) {

	case '1':
	  c = one;
	  break;

	case '0':
	  c = zero;
	  break;

	case '?':
	  c = quest;
	  break;

	case 'P':
	  c = polym;
	  break;
	}
	d = c;
      }
      if (n > 1) {
	makechar(d);
	prefix(c);
	for (j = 1; j <= n - 2; j++)
	  putchar(cha[(long)c]);
	postfix(c);
      }
      if (q->index >= 10)
	printf("%2hd", q->index);
      else {
	makechar(c);
	printf("%hd", q->index);
      }
      extra = true;
    } else if (!q->tip) {
      if (last->ycoord > i && first->ycoord < i && i != p->ycoord) {
	c = up;
	if (i < p->ycoord)
	  s = p->next->back->state;
	else
	  s = p->next->next->back->state;
	if (s == 'P' && p->state != 'P')
	  s = p->state;
	if (display) {
	  switch (s) {

	  case '1':
	    c = one;
	    break;

	  case '0':
	    c = zero;
	    break;

	  case '?':
	    c = quest;
	    break;

	  case 'P':
	    c = polym;
	    break;
	  }
	}
	makechar(c);
	for (j = 1; j < n; j++)
	  putchar(' ');
      } else {
	for (j = 1; j <= n; j++)
	  putchar(' ');
      }
    } else {
      for (j = 1; j <= n; j++)
	putchar(' ');
    }
    if (p != q)
      p = q;
  } while (!done);
  if (p->ycoord == i && p->tip) {
    putchar(':');
    for (j = 0; j < nmlngth; j++)
      putchar(nayme[p->index - 1][j]);
  }
  if (i == 1) {
    if (display)
      printf("   CHARACTER%4hd", dispchar);
    else
      printf("            ");
    if (subtree)
      printf(" Subtree ");
    else
      printf("         ");
    if (dollo)
      printf("Dollo");
    else
      printf("Polymorphism ");
  }
  if (i == 3 && display) {
    printf("   ");
    makechar(one);
    printf(":1 ");
    makechar(quest);
    printf(":? ");
    makechar(zero);
    printf(":0 ");
    if (!dollo) {
      putchar(' ');
      makechar(polym);
      printf(":0/1");
    }
  }
  if ((i == 5 || (lastline < 5 && i == lastline)) && !earlytree) {
    if (lastline < 5)
      printf("\n                   ");
    printf("   Total:%10.5f", -like);
    gotlike = -like;
  }
  if ((i == 7 || (lastline < 7 && i == lastline)) && changed && !earlytree) {
    if (lastline < 7)
      printf("\n                     ");
    if (-like < bestyet) {
      printf("     BEST YET!");
      bestyet = -like;
    } else if (fabs(-like - bestyet) < 0.000001)
      printf("     (as good as best)");
    else {
      if (-like < gotlike)
        printf("     better");
      else if (-like > gotlike)
        printf("     worse!");
    }
  }
  if ((i == 9 || (lastline < 9 && i == lastline)) && !earlytree) {
    if (lastline < 9)
      printf("\n                       ");
    printf("     %3hd characters compatible", compatible);
  }
  putchar('\n');
}  /* drawline */

void printree()
{
  /* prints out diagram of the tree */
  short tipy;
  short i, rover, dow;

  if (!subtree)
    nuroot = root;
  if (changed || newtree)
    evaluate(root);
  if (display)
    hypstates();
  if (ansi || ibmpc)
    printf("\033[2J\033[H");
  else if (vt52)
    printf("\033E\033H");
  else
    putchar('\n');
  tipy = 1;
  rover = 100 / spp;
  if (rover > overr)
    rover = overr;
  dow = downn;
  if (spp * dow > screenlines && !subtree) {
    dow--;
    rover--;
  }
  coordinates(nuroot, &tipy);
  for (i = 1; i <= (tipy - dow); i++)
    drawline(i, tipy - dow);
  if (spp <= screenlines / 2 ||subtree)
    putchar('\n');
  changed = false;
}  /* printree */



void arbitree()
{
  short i;

  root = treenode[0];
  add(treenode[0], treenode[1], treenode[spp]);
  for (i = 3; i <= (spp); i++)
    add(treenode[spp + i - 3], treenode[i - 1], treenode[spp + i - 2]);
  for (i = 0; i < (nonodes); i++)
    intree[i] = true;
}  /* arbitree */

void yourtree()
{
  short i, j;
  boolean ok;

  root = treenode[0];
  add(treenode[0], treenode[1], treenode[spp]);
  i = 2;
  do {
    i++;
    printree();
    printf("Add species%3hd: ", i);
    for (j = 0; j < nmlngth; j++)
      putchar(nayme[i - 1][j]);
    do {
      printf("\nbefore node (type number): ");
      inpnum(&j, &ok);
      ok = (ok && ((j >= 1 && j < i) || (j > spp && j < spp + i - 1)));
      if (!ok)
	printf("Impossible number. Please try again:\n");
    } while (!ok);
    add(treenode[j - 1], treenode[i - 1], treenode[spp + i - 2]);
  } while (i != spp);
  for (i = 0; i < (nonodes); i++)
    intree[i] = true;
}  /* yourtree */


void findch(c)
Char c;
{
  /* scan forward until find character c */
  boolean done;

  done = false;
  while (!(done)) {
    if (c == ',') {
      if (ch == '(' || ch == ')' || ch == ';') {
	printf("\nERROR IN USER TREE: UNMATCHED PARENTHESIS OR MISSING COMMA\n");
	exit(-1);
      } else if (ch == ',')
	done = true;
    } else if (c == ')') {
      if (ch == '(' || ch == ',' || ch == ';') {
	printf("\nERROR IN USER TREE: ");
	printf("UNMATCHED PARENTHESIS OR NOT BIFURCATED NODE\n");
	exit(-1);
      } else {
	if (ch == ')')
	  done = true;
      }
    } else if (c == ';') {
      if (ch != ';') {
	printf("\nERROR IN USER TREE: UNMATCHED PARENTHESIS");
	printf(" OR MISSING SEMICOLON\n");
	exit(-1);
      } else
	done = true;
    }
    if ((done && ch == ')') || !(done)) {
      if (eoln(treefile)) {
	fscanf(treefile, "%*[^\n]");
	getc(treefile);
      }
      ch = getc(treefile);
    }
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
    if (eoln(treefile)) {
      fscanf(treefile, "%*[^\n]");
      getc(treefile);
    }
    ch = getc(treefile);
  } while (!(ch != ' '));
  if (ch == '(' ) {
    if ((*lparens) >= spp - 1) {
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
    intree[*nextnode - 1] = true;
    return;
  }
  for (i = 0; i < nmlngth; i++)
    str[i] = ' ';
  n = 1;
  do {
    if (ch == '_')
      ch = ' ';
    str[n - 1] = ch;
    if (eoln(treefile)) {
      fscanf(treefile, "%*[^\n]");
      getc(treefile);
    }
    ch = getc(treefile);
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
	intree[n - 1] = true;
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
  /* Local variables for treeread:           */
  short nextnode, lparens;
  boolean *names;
  short i;

  openfile(&treefile,TREEFILE,"r",progname,trfilename);

  root = treenode[spp];
  nextnode = spp;
  root->back = NULL;
  for (i = 0; i < nonodes; i++)
    intree[i] = false;
  names = (boolean *)Malloc(spp*sizeof(boolean));
  for (i = 0; i < spp; i++)
    names[i] = false;
  lparens = 0;
  addelement(&root, &nextnode,&lparens,names);
  findch(';');
  FClose(treefile);
#ifdef MAC
   fixmacfile(trfilename);
#endif
  free(names);
}  /* treeread */

void buildtree()
{

  changed = false;
  newtree = false;
  switch (how) {

  case arb:
    arbitree();
    break;

  case use:
    treeread();
    break;

  case spec:
    yourtree();
    break;
  }
  outgrno = root->next->back->index;
  if (intree[outgrno - 1])
    reroot(treenode[outgrno - 1]);
}  /* buildtree */


void help()
{
  /* display help information */
  printf("\n\nR Rearrange a tree by moving a node or group\n");
  printf("# Show the states of the next char. that doesn't fit tree\n");
  printf("+ Show the states of the next character\n");
  printf("-         ...     of the previous character\n");
  printf("S Show the states of a given character\n");
  printf(". redisplay the same tree again\n");
  printf("T Try all possible positions of a node or group\n");
  printf("U Undo the most recent rearrangement\n");
  printf("W Write tree to a file\n");
  printf("O select an Outgroup for the tree\n");
  printf("F Flip (rotate) branches at a node\n");
  printf("C show only one Clade (subtree) (useful if tree is too big)\n");
  printf("H Help (this screen)\n");
  printf("?  \"     \"      \"   \n");
  printf("Q (Quit) Exit from program\n");
  printf("X Exit from program\n\n\n");
  printf("TO CONTINUE, PRESS ON THE Return OR Enter KEY");
  scanf("%c%*[^\n]", &ch);
  getchar();
  printree();
}  /* help */

void rearrange()
{
  short i, j;
  boolean ok1, ok2;
  node *p, *q;

  printf("Remove everything to the right of which node? ");
  inpnum(&i, &ok1);
  ok1 = (ok1 && i >= 1 && i < spp * 2 && i != root->index);
  if (ok1) {
    printf("Add before which node? ");
    inpnum(&j, &ok2);
    ok2 = (ok2 && j >= 1 && j < spp * 2);
    if (ok2) {
      ok2 = (treenode[j - 1] != treenode[treenode[i - 1]->back->index - 1]);
      p = treenode[j - 1];
      while (p != root) {
	ok2 = (ok2 && p != treenode[i - 1]);
	p = treenode[p->back->index - 1];
      }
      if (ok1 && ok2) {
	what = i;
	q = treenode[treenode[i - 1]->back->index - 1];
	if (q->next->back->index == i)
	  fromwhere = q->next->next->back->index;
	else
	  fromwhere = q->next->back->index;
	towhere = j;
	re_move(&treenode[i - 1], &q);
	add(treenode[j - 1], treenode[i - 1], q);
      }
      lastop = rearr;
    }
  }
  changed = (ok1 && ok2);
  printree();
  if (!(ok1 && ok2))
    printf("Not a possible rearrangement.  Try again: \n");
  else {
    oldwritten = written;
    written = false;
  }
}  /* rearrange */

void nextinc()
{
  /* show next incompatible character */
  short disp0;
  boolean done;

  display = true;
  disp0 = dispchar;
  done = false;
  do {
    dispchar++;
    if (dispchar > chars) {
      dispchar = 1;
      done = (disp0 == 0);
    }
  } while (!(numsteps[dispchar - 1] >
	     weight[dispchar - 1] ||
	     dispchar == disp0 || done));
  dispword = (dispchar - 1) / bits + 1;
  dispbit = (dispchar - 1) % bits + 1;
  printree();
}  /* nextinc */

void nextchar()
{
  /* show next character */
  display = true;
  dispchar++;
  if (dispchar > chars)
    dispchar = 1;
  dispword = (dispchar - 1) / bits + 1;
  dispbit = (dispchar - 1) % bits + 1;
  printree();
}  /* nextchar */

void prevchar()
{
  /* show previous character */
  display = true;
  dispchar--;
  if (dispchar < 1)
    dispchar = chars;
  dispword = (dispchar - 1) / bits + 1;
  dispbit = (dispchar - 1) % bits + 1;
  printree();
}  /* prevchar */

void show()
{
  short i;
  boolean ok;

  do {
    printf("SHOW: (Character number or 0 to see none)? ");
    inpnum(&i, &ok);
    ok = (ok && (i == 0 || (i >= 1 && i <= chars)));
    if (ok && i != 0) {
      display = true;
      dispchar = i;
      dispword = (i - 1) / bits + 1;
      dispbit = (i - 1) % bits + 1;
    }
    if (ok && i == 0)
      display = false;
  } while (!ok);
  printree();
}  /* show */

void tryadd(p, item,nufork,place)
node *p,**item,**nufork;
double *place;
{
  /* temporarily adds one fork and one tip to the tree.
    Records scores in ARRAY place */
  add(p, *item, *nufork);
  evaluate(root);
  place[p->index - 1] = -like;
  re_move(item, nufork);
}  /* tryadd */

void addpreorder(p, item_, nufork_, place)
node *p, *item_, *nufork_;
double *place;
{
  /* traverses a binary tree, calling PROCEDURE tryadd
    at a node before calling tryadd at its descendants */
  node *item, *nufork;


  item = item_;
  nufork = nufork_;
  if (p == NULL)
    return;
  tryadd(p, &item,&nufork,place);
  if (!p->tip) {
    addpreorder(p->next->back, item,nufork,place);
    addpreorder(p->next->next->back,item,nufork,place);
  }
}  /* addpreorder */

void try()
{
  /* Remove node, try it in all possible places */
/* Local variables for try: */
  double *place;
  short i, j, oldcompat;
  double current;
  node *q, *dummy, *rute;
  boolean tied, better, ok;

  printf("Try other positions for which node? ");
  inpnum(&i, &ok);
  if (!(ok && i >= 1 && i <= nonodes && i != root->index)) {
    printf("Not a possible choice! ");
    return;
  }
  printf("WAIT ...\n");
  place = (double *)Malloc(nonodes*sizeof(double));
  for (j = 0; j < (nonodes); j++)
    place[j] = -1.0;
  evaluate(root);
  current = -like;
  oldcompat = compatible;
  what = i;
  q = treenode[treenode[i - 1]->back->index - 1];
  if (q->next->back->index == i)
    fromwhere = q->next->next->back->index;
  else
    fromwhere = q->next->back->index;
  rute = root;
  if (root->index == treenode[i - 1]->back->index) {
    if (treenode[treenode[i - 1]->back->index - 1]->next->back == treenode[i - 1])
      rute = treenode[treenode[i - 1]->back->index - 1]->next->next->back;
    else
      rute = treenode[treenode[i - 1]->back->index - 1]->next->back;
  }
  re_move(&treenode[i - 1], &dummy);
  oldleft = wasleft;
  root = rute;
  addpreorder(root, treenode[i - 1], dummy, place);
  wasleft = oldleft;
  restoring = true;
  add(treenode[fromwhere - 1], treenode[what - 1],
      dummy);
  like = -current;
  compatible = oldcompat;
  restoring = false;
  better = false;
  printf("       BETTER: ");
  for (j = 1; j <= (nonodes); j++) {
    if (place[j - 1] < current && place[j - 1] >= 0.0) {
      printf("%3hd:%6.2f", j, place[j - 1]);
      better = true;
    }
  }
  if (!better)
    printf(" NONE");
  printf("\n       TIED:    ");
  tied = false;
  for (j = 1; j <= (nonodes); j++) {
    if (fabs(place[j - 1] - current) < 1.0e-6 && j != fromwhere) {
      if (j < 10)
	printf("%2hd", j);
      else
	printf("%3hd", j);
      tied = true;
    }
  }
  if (tied)
    printf(":%6.2f\n", current);
  else
    printf("NONE\n");
  changed = true;
  free(place);
}  /* try */

void undo()
{
  /* restore to tree before last rearrangement */
  short temp;
  boolean btemp;
  node *q;

  switch (lastop) {

  case rearr:
    restoring = true;
    oldleft = wasleft;
    re_move(&treenode[what - 1], &q);
    btemp = wasleft;
    wasleft = oldleft;
    add(treenode[fromwhere - 1], treenode[what - 1],q);
    wasleft = btemp;
    restoring = false;
    temp = fromwhere;
    fromwhere = towhere;
    towhere = temp;
    changed = true;
    break;

  case flipp:
    q = treenode[atwhat - 1]->next->back;
    treenode[atwhat - 1]->next->back =
      treenode[atwhat - 1]->next->next->back;
    treenode[atwhat - 1]->next->next->back = q;
    treenode[atwhat - 1]->next->back->back = treenode[atwhat - 1]->next;
    treenode[atwhat - 1]->next->next->back->back =
      treenode[atwhat - 1]->next->next;
    break;

  case reroott:
    restoring = true;
    temp = oldoutgrno;
    oldoutgrno = outgrno;
    outgrno = temp;
    reroot(treenode[outgrno - 1]);
    restoring = false;
    break;

  case none:
    /* blank case */
    break;
  }
  printree();
  if (lastop == none) {
    printf("No operation to undo! ");
    return;
  }
  btemp = oldwritten;
  oldwritten = written;
  written = btemp;
}  /* undo */

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
  if (p == root)
    fprintf(treefile, ";\n");
}  /* treeout */

void treewrite(done)
boolean done;
{
  /* write out tree to a file */
  Char ch;

  if (waswritten) {   /*MIPS*/
    printf("Tree file already was open.\n");
    printf("   A   Add to this tree to tree file\n");
    printf("   R   Replace tree file contents by this tree\n");
    printf("   F   Write out tree to a different tree file\n");
    printf("   N   Do Not write out this tree\n");
    do {
      printf("Which should we do? ");
      scanf("%c%*[^\n]", &ch);
      getchar();
      uppercase(&ch);
    } while (ch != 'A' && ch != 'R' && ch != 'N' && ch != 'F');
  }
  if (ch == 'F'){
    trfilename[0] = '\0';
    while (trfilename[0] =='\0'){
      printf("Please enter a tree file name>");
      gets(trfilename);}
  }
  if (ch == 'R' || ch == 'A' || ch == 'F' || !waswritten){
    openfile(&treefile,trfilename,(ch == 'A' && waswritten) ? "a" : "w",
	     progname,trfilename);
  }

  if (!done)
    printree();
  if (waswritten && ch == 'N')
    return;
  col = 0;
  treeout(root);
  printf("\nTree written to file\n\n");
  waswritten = true;
  written = true;
  FClose(treefile);
#ifdef MAC
  fixmacfile(trfilename);
#endif
}  /* treewrite */

void clade()
{
  /* pick a subtree and show only that on screen */
  short i;
  boolean ok;

  printf("Select subtree rooted at which node (0 for whole tree)? ");
  inpnum(&i, &ok);
  ok = (ok && (unsigned)i <= nonodes);
  if (ok) {
    subtree = (i > 0);
    if (subtree)
      nuroot = treenode[i - 1];
    else
      nuroot = root;
  }
  printree();
  if (!ok)
    printf("Not possible to use this node. ");
}  /* clade */

void flip()
{
  /* flip at a node left-right */
  short i;
  boolean ok;
  node *p;

  printf("Flip branches at which node? ");
  inpnum(&i, &ok);
  ok = (ok && i > spp && i <= nonodes);
  if (ok) {
    p = treenode[i - 1]->next->back;
    treenode[i - 1]->next->back = treenode[i - 1]->next->next->back;
    treenode[i - 1]->next->next->back = p;
    treenode[i - 1]->next->back->back = treenode[i - 1]->next;
    treenode[i - 1]->next->next->back->back = treenode[i - 1]->next->next;
    atwhat = i;
    lastop = flipp;
  }
  printree();
  if (ok) {
    oldwritten = written;
    written = false;
    return;
  }
  if (i >= 1 && i <= spp)
    printf("Can't flip there. ");
  else
    printf("No such node. ");
}  /* flip */

void changeoutgroup()
{
  short i;
  boolean ok;

  oldoutgrno = outgrno;
  do {
    printf("Which node should be the new outgroup? ");
    inpnum(&i, &ok);
    ok = (ok && intree[i - 1] && i >= 1 && i <= nonodes &&
	  i != root->index);
    if (ok)
      outgrno = i;
  } while (!ok);
  if (intree[outgrno - 1])
    reroot(treenode[outgrno - 1]);
  changed = true;
  lastop = reroott;
  printree();
  oldwritten = written;
  written = false;
}  /* changeoutgroup */

void redisplay()
{
/* Local variables for redisplay: */
  boolean done;

  done = false;
  waswritten = false;
  do {
    printf("NEXT? (Options: R # + - S . T U W O F C H ? X Q) ");
    printf("(H or ? for Help) ");
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    if (strchr("RSWH#.O?+TFX-UCQ",ch)){
      switch (ch) {

      case 'R':
	rearrange();
	break;

      case '#':
	nextinc();
	break;

      case '+':
	nextchar();
	break;

      case '-':
	prevchar();
	break;

      case 'S':
	show();
	break;

      case '.':
	printree();
	break;

      case 'T':
	try();
	break;

      case 'U':
	undo();
	break;

      case 'W':
	treewrite(done);
	break;

      case 'O':
	changeoutgroup();
	break;

      case 'F':
	flip();
	break;

      case 'C':
	clade();
	break;

      case 'H':
	help();
	break;

      case '?':
	help();
	break;

      case 'X':
	done = true;
	break;

      case 'Q':
	done = true;
	break;
      }
    }
  } while (!done);
  if (!written) {
    do {
      printf("Do you want to write out the tree to a file? (Y or N) ");
      scanf("%c%*[^\n]", &ch);
      getchar();
    } while (ch != 'Y' && ch != 'y' && ch != 'N' && ch != 'n');
  }
  if (ch == 'Y' || ch == 'y')
    treewrite(done);
}  /* redisplay */


void treeconstruct()
{
  /* constructs a binary tree from the pointers in treenode. */

  restoring = false;
  subtree = false;
  display = false;
  dispchar = 0;
  fullset = (1L << (bits + 1)) - (1L << 1);
  guess = (Char *)Malloc(chars*sizeof(Char));
  numsteps = (steptr)Malloc(chars*sizeof(short));
  earlytree = true;
  buildtree();
  waswritten = false;
  printf("\nComputing steps needed for compatibility in sites ...\n\n");
  newtree = true;
  earlytree = false;
  printree();
  bestyet = -like;
  gotlike = -like;
  lastop = none;
  newtree = false;
  written = false;
  lastop = none;
  redisplay();
}  /* treeconstruct */


main(argc, argv)
int argc;
Char *argv[];
{  /* Interactive Dollo/polymorphism parsimony */
  /* reads in spp, chars, and the data. Then calls treeconstruct to
    construct the tree and query the user */
#ifdef MAC
  macsetup("Dolmove","");
  argv[0] = "Dolmove";
#endif
  strcpy(progname,argv[0]);
  strcpy(infilename,INFILE);
  strcpy(trfilename,TREEFILE);

  openfile(&infile,infilename,"r",argv[0],infilename);

  screenlines = 24;
  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  doinput();
  configure();
  treeconstruct();
  if (waswritten) {
    FClose(treefile);
#ifdef MAC
    fixmacfile(trfilename);
#endif
  }
  FClose(infile);
  exit(0);
}  /* Interactive Dollo/polymorphism parsimony */



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

