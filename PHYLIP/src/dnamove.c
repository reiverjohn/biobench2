#include "phylip.h"

/* version 3.56c. (c) Copyright 1993 by Joseph Felsenstein.
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
typedef enum {  A, C, G, U, O } bases;
typedef long *baseptr;
typedef enum {
  horiz, vert, up, over, upcorner, downcorner, aa, cc, gg, tt, quest
  } chartype;
/* nodes will form a binary tree */
typedef enum {
  rearr, flipp, reroott, none
  } rearrtype;

typedef struct node {         /* describes a tip species or an ancestor */
  struct node *next, *back;   /* pointers to nodes                      */
  short index;                /* number of the node                     */
  boolean tip;                /* present species are tips of tree       */
  baseptr base;             /* the sequence                           */
  Char state;                 /* stores state for printing              */
  short xcoord, ycoord;       /* used by printree                       */
  short ymin, ymax;
} node;

typedef node **pointptr;
typedef Char naym[nmlngth];
typedef struct gbase {
  baseptr base;
  struct gbase *next;
} gbase;

typedef enum {
  arb, use, spec
  } howtree;

char infilename[100],outfilename[100],trfilename[100];
Static node *root;
Static FILE *infile, *outfile, *treefile;
Static short spp, nonodes, chars, outgrno, screenlines, col;
/* spp = number of species
   nonodes = number of nodes in tree
   chars = number of sites in actual sequences
   outgrno indicates outgroup */
Static boolean outgropt, weights, thresh,  waswritten, interleaved,
  ibmpc, vt52, ansi;
Static steptr weight;
Static pointptr treenode;   /* pointers to all nodes in tree */
Static naym *nayme;   /* names of species */
Static double threshold;
Static double *threshwt;
Static boolean reversed[11];
Static boolean graphic[11];
Static Char chh[11];
Static howtree how;
Static gbase *garbage;
Char progname[20];
Char  ch;
/* Local variables for treeconstruct, propogated global for C version: */

short dispchar, atwhat, what, fromwhere, towhere, oldoutgrno, compatible;
double like, bestyet, gotlike;
boolean display, newtree, changed, subtree, written, oldwritten, restoring,
  wasleft, oldleft, earlytree;
steptr necsteps;
boolean *intree;
long sett[31];
steptr numsteps;
node *nuroot;
rearrtype lastop;


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
     gbase **p;
{
  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */
  if (garbage != NULL) {
    *p = garbage;
    garbage = garbage->next;
  } else {
    *p = (gbase *)Malloc(sizeof(gbase));
    (*p)->base = (baseptr)Malloc(chars*sizeof(long));
  }
  (*p)->next = NULL;
}  /* gnu */


void chuck(p)
     gbase *p;
{
  /* collect garbage on p -- put it on front of garbage list */
  p->next = garbage;
  garbage = p;
}  /* chuck */

void uppercase(ch)
     Char *ch;
{  /* convert a character to upper case -- either ASCII or EBCDIC */
  *ch = (islower(*ch) ?  toupper(*ch) : (*ch));
}  /* uppercase */



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
  /* input the numbers of species and of sites */
  fscanf(infile, "%hd%hd", &spp, &chars);
  printf("\n%2hd species, %3hd  sites\n", spp, chars);
  putc('\n', outfile);
  nonodes = spp * 2 - 1;
}  /* inputnumbers */

void getoptions()
{
  /* interactively set options */
  Char ch;
  boolean done, done1, gotopt;

  how = arb;
  outgrno = 1;
  outgropt = false;
  thresh = false;
  weights = false;
  interleaved = true;
  do {
    printf((ansi || ibmpc) ?  "\033[2J\033[H" :
	   vt52 ?  "\033E\033H"    : "\n");
    printf("\nInteractive DNA parsimony, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  O                           Outgroup root?");
    if (outgropt)
      printf("  Yes, at sequence number%3hd\n", outgrno);
    else
      printf("  No, use as outgroup species%3hd\n", outgrno);
    printf("  T                 Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count up to%4.1f per site\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  I             Input sequences interleaved?  %s\n",
	   (interleaved ? "Yes" : "No, sequential"));

    printf("  U Initial tree (arbitrary, user, specify)?  %s\n",
	   (how == arb) ? "Arbitrary"                :
	   (how == use) ? "User tree from tree file" : "Tree you specify");
    printf("  0      Graphics type (IBM PC, VT52, ANSI)?  %s\n",
	   ibmpc ? "IBM PC\n" :
	   ansi  ? "ANSI"     :
	   vt52  ? "VT52"     : "(none)");
    printf("  L               Number of lines on screen?%4hd\n",screenlines);
    do {
      printf("\nAre these settings correct? ");
      printf("(type Y or the letter for one to change)\n");
      scanf("%c%*[^\n]", &ch);
      getchar();
      uppercase(&ch);
      done = (ch == 'Y');
      gotopt = (strchr("OTIU0L",ch) != NULL) ? true : false;
      if (gotopt) {
        switch (ch) {
	
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
	
        case 'I':
          interleaved = !interleaved;
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
            scanf("%hd%*[^\n]", &screenlines);
            getchar();
          } while (screenlines <= 12);
          break;
        }
      }
      if (!(gotopt || done))
        printf("Not a possible option!\n");
    } while (!(gotopt || done));
  } while (!done);
}  /* getoptions */

void inputweights()
{
  /* input the site weights, 0-9 and A-Z for weights 0 - 35 */
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
  /* print out the weights of sites */
  short i, j, k;

  printf("Sites are weighted as follows:\n");
  printf("        ");
  for (i = 0; i <= 9; i++)
    printf("%3hd", i);
  putc('\n', outfile);
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

void inputoptions()
{
  /* input the information on the options */
  Char ch;
  short extranum, i;

  extranum = 0;
  while (!eoln(infile)) {
    ch = getc(infile);
    uppercase(&ch);
    if (ch == 'W')
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
    if (ch == 'W')
      inputweights();
    if (ch != 'W'){
      printf("ERROR: INCORRECT AUXILIARY OPTIONS LINE WHICH STARTS WITH %c\n",
	     ch);
      exit(-1);}
  }
  if (weights)
    printweights();
  if (!thresh)
    threshold = spp;
  for (i = 0; i < (chars); i++)
    threshwt[i] = threshold * weight[i];
}  /* inputoptions */

void inputdata()
{
  /* input the names and sequences for each species */
  short i, j, basesread, basesnew;
  Char charstate;
  boolean allread, done;
  long ns;   /* temporary base set for input */
  node *p, *q;

  treenode = (node **)Malloc(nonodes*sizeof(node *));
  for (i = 0; i < (spp); i++) {
    treenode[i] = (node *)Malloc(sizeof(node));
    treenode[i]->base = (baseptr)Malloc(chars*sizeof(long));
  }
  for (i = spp; i < (nonodes); i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc(sizeof(node));
      p->base = (baseptr)Malloc(chars*sizeof(long));
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
    }
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
	  if (eof(infile) || eoln(infile)){
	    printf("ERROR: END-OF-LINE OR END-OF-FILE");
	    printf(" IN THE MIDDLE OF A SPECIES NAME\n");
            exit(-1);}
	  nayme[i - 1][j] = getc(infile);
        }
      }
      if (interleaved)
        j = basesread;
      else
        j = 0;
      done = false;
      while (!done && !eof(infile)) {
        if (interleaved)
          done = true;
        while (j < chars && !(eoln(infile) || eof(infile))) {
          charstate = getc(infile);
	  if (ch == '\n')
	    ch = ' ';
          if (charstate == ' ' || (charstate >= '0' && charstate <= '9'))
            continue;
          uppercase(&charstate);
	  if (strchr("ABCDGHKMNRSTUVWXY?O-.",charstate) == NULL){
            printf("ERROR: BAD BASE:%c AT POSITION%5hd OF SPECIES %3hd\n",
                   charstate, j, i);
	    exit(-1);
          }
          j++;
          switch (charstate) {
	
          case 'A':
            ns = 1L << ((long)A);
            break;
	
          case 'C':
            ns = 1L << ((long)C);
            break;
	
          case 'G':
            ns = 1L << ((long)G);
            break;
	
          case 'U':
            ns = 1L << ((long)U);
            break;
	
          case 'T':
            ns = 1L << ((long)U);
            break;
	
          case 'M':
            ns = (1L << ((long)A)) | (1L << ((long)C));
            break;
	
          case 'R':
            ns = (1L << ((long)A)) | (1L << ((long)G));
            break;
	
          case 'W':
            ns = (1L << ((long)A)) | (1L << ((long)U));
            break;
	
          case 'S':
            ns = (1L << ((long)C)) | (1L << ((long)G));
            break;
	
          case 'Y':
            ns = (1L << ((long)C)) | (1L << ((long)U));
            break;
	
          case 'K':
            ns = (1L << ((long)G)) | (1L << ((long)U));
            break;
	
          case 'B':
            ns = (1L << ((long)C)) | (1L << ((long)G)) | (1L << ((long)U));
            break;
	
          case 'D':
            ns = (1L << ((long)A)) | (1L << ((long)G)) | (1L << ((long)U));
            break;
	
          case 'H':
            ns = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)U));
            break;
	
          case 'V':
            ns = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G));
            break;
	
          case 'N':
            ns = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G)) |
	      (1L << ((long)U));
            break;
	
          case 'X':
            ns = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G)) |
	      (1L << ((long)U));
            break;
	
          case '?':
            ns = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G)) |
   	         (1L << ((long)U)) | (1L << ((long)O));
            break;
	
          case 'O':
            ns = 1L << ((long)O);
            break;
	
          case '.':
            ns = treenode[0]->base[j - 1];
            break;
	
          case '-':
            ns = 1L << ((long)O);
            break;
          }
          treenode[i - 1]->base[j - 1] = ns;
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
  root = NULL;
  printf("\n\n");
}  /* inputdata */


void doinput()
{
  /* reads the input data */
  short j;

  inputnumbers();
  getoptions();
  printf("\nReading input file ...\n\n");
  nayme = (naym *)Malloc(spp*sizeof(naym));
  intree = (boolean *)Malloc(nonodes*sizeof(boolean));
  weight = (steptr)Malloc(chars*sizeof(short));
  numsteps = (steptr)Malloc(chars*sizeof(short));
  necsteps = (steptr)Malloc(chars*sizeof(short));
  threshwt = (double *)Malloc(chars*sizeof(double));
  inputoptions();
  inputdata();
}  /* doinput */

void configure()
{
  /* configure to machine -- set up special characters */
  chartype a;

  for (a = horiz; (long)a <= (long)quest; a = (chartype)((long)a + 1))
    reversed[(long)a] = false;
  for (a = horiz; (long)a <= (long)quest; a = (chartype)((long)a + 1))
    graphic[(long)a] = false;
  if (ibmpc) {
    chh[(long)horiz] = 205;
    graphic[(long)horiz] = true;
    chh[(long)vert] = 186;
    graphic[(long)vert] = true;
    chh[(long)up] = 186;
    graphic[(long)up] = true;
    chh[(long)over] = 205;
    graphic[(long)over] = true;
    chh[(long)upcorner] = 200;
    graphic[(long)upcorner] = true;
    chh[(long)downcorner] = 201;
    graphic[(long)downcorner] = true;
    chh[(long)aa] = 176;
    chh[(long)cc] = 178;
    chh[(long)gg] = 177;
    chh[(long)tt] = 219;
    chh[(long)quest] = '\001';
    return;
  }
  if (vt52) {
    chh[(long)horiz] = ' ';
    reversed[(long)horiz] = true;
    chh[(long)vert] = chh[(long)horiz];
    reversed[(long)vert] = true;
    chh[(long)up] = '`';
    graphic[(long)up] = true;
    chh[(long)over] = 'a';
    graphic[(long)over] = true;
    chh[(long)upcorner] = 'e';
    graphic[(long)upcorner] = true;
    chh[(long)downcorner] = 'f';
    graphic[(long)downcorner] = true;
    chh[(long)aa] = 'w';
    graphic[(long)aa] = true;
    chh[(long)cc] = ' ';
    reversed[(long)cc] = true;
    chh[(long)gg] = 'j';
    graphic[(long)gg] = true;
    chh[(long)tt] = 'i';
    graphic[(long)tt] = true;
    chh[(long)quest] = '?';
    reversed[(long)quest] = true;
    return;
  }
  if (ansi) {
    chh[(long)horiz] = ' ';
    reversed[(long)horiz] = true;
    chh[(long)vert] = chh[(long)horiz];
    reversed[(long)vert] = true;
    chh[(long)up] = 'x';
    graphic[(long)up] = true;
    chh[(long)over] = 'q';
    graphic[(long)over] = true;
    chh[(long)upcorner] = 'm';
    graphic[(long)upcorner] = true;
    chh[(long)downcorner] = 'l';
    graphic[(long)downcorner] = true;
    chh[(long)aa] = 'a';
    reversed[(long)aa] = true;
    chh[(long)cc] = 'c';
    reversed[(long)cc] = true;
    chh[(long)gg] = 'g';
    reversed[(long)gg] = true;
    chh[(long)tt] = 't';
    reversed[(long)tt] = true;
    chh[(long)quest] = '?';
    reversed[(long)quest] = true;
    return;
  }
  chh[(long)horiz] = '=';
  chh[(long)vert] = ' ';
  chh[(long)up] = '!';
  chh[(long)upcorner] = '`';
  chh[(long)downcorner] = ',';
  chh[(long)over] = '-';
  chh[(long)aa] = 'a';
  chh[(long)cc] = 'c';
  chh[(long)gg] = 'g';
  chh[(long)tt] = 't';
  chh[(long)quest] = '.';
}  /* configure */


void pregraph()
{
  /* turn on graphic characters */
  if (vt52)
    printf("\033F");
  if (ansi) {
    printf("\033(0");
    printf("\033[10m");
  }
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
  if (ansi) {
    printf("\033[11m");
    printf("\033(B");
  }
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
  putchar(chh[(long)a]);
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
    putleft = wasleft;
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
  /* sets up for each node in the tree the base sequence
     at that point and counts the changes.  The program
     spends much of its time in this PROCEDURE */
  short i;
  long ns, rs, ls;

  for (i = 0; i < (chars); i++) {
    ls = p->next->back->base[i];
    rs = p->next->next->back->base[i];
    ns = ls & rs;
    if (ns == 0) {
      ns = ls | rs;
      numsteps[i] += weight[i];
    }
    p->base[i] = ns;
  }
}  /* fillin */

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

void evaluate(r)
     node *r;
{
  /* determines the number of steps needed for a tree. this is
     the minimum number of steps needed to evolve sequences on
     this tree */
  short i, steps;
  double sum;

  compatible = 0;
  sum = 0.0;
  for (i = 0; i < (chars); i++)
    numsteps[i] = 0;
  postorder(r);
  for (i = 0; i < (chars); i++) {
    steps = numsteps[i];
    if (steps <= threshwt[i])
      sum += steps;
    else
      sum += threshwt[i];
    if (steps <= necsteps[i] && !earlytree)
      compatible += weight[i];
  }
  like = -sum;
  /*printf("like: %f\n",like);*/
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


void ancestset(a, b, c)
     long a, b, *c;
{
  /* make the set of ancestral states below nodes
     whose base sets are a and b */
  *c = a & b;
  if (*c == 0)
    *c = a | b;
}  /* ancestset */

void firstrav(r,i)
     node *r;
     short i;
{
  /* initial traverse for hypothetical states */
  if (r->tip)
    return;
  firstrav(r->next->back, i);
  firstrav(r->next->next->back, i);
  ancestset(r->next->back->base[i - 1],
            r->next->next->back->base[i - 1], &r->base[i - 1]);
}  /* firstrav */

void hyptrav(r, hypset, i,bottom)
     node *r;
     long *hypset;
     short i;
     boolean *bottom;
{
  /*  compute, print out state at one interior node */
  long tempset, left, rt, anc;
  gbase *temparray, *ancset;

  gnu(&ancset);
  gnu(&temparray);
  anc = hypset[i - 1];
  if (!r->tip) {
    left = r->next->back->base[i - 1];
    rt = r->next->next->back->base[i - 1];
    tempset = left & rt & anc;
    if (tempset == 0) {
      tempset = (left & rt) | (left & anc) | (rt & anc);
      if (tempset == 0)
        tempset = left | rt | anc;
    }
    r->base[i - 1] = tempset;
  }
  if (!(*bottom))
    anc = treenode[r->back->index - 1]->base[i - 1];
  r->state = '?';
  if (r->base[dispchar - 1] == 1L << ((long)A))
    r->state = 'A';
  if (r->base[dispchar - 1] == 1L << ((long)C))
    r->state = 'C';
  if (r->base[dispchar - 1] == 1L << ((long)G))
    r->state = 'G';
  if (r->base[dispchar - 1] == 1L << ((long)U))
    r->state = 'T';
  /*printf("%c",r->state);*/
  *bottom = false;
  if (!r->tip) {
    memcpy(temparray->base, r->next->back->base, chars*sizeof(long));
    ancestset(hypset[i - 1], r->next->next->back->base[i - 1],
              &ancset->base[i - 1]);
    hyptrav(r->next->back, ancset->base, i,bottom);
    ancestset(hypset[i - 1], temparray->base[i - 1],
              &ancset->base[i - 1]);
    hyptrav(r->next->next->back, ancset->base, i,bottom);
  }
  chuck(temparray);
  chuck(ancset);
}  /* hyptrav */

void hypstates()
{
  /* fill in and describe states at interior nodes */
  /* Local variables for hypstates: */
  short i;
  boolean bottom;
  baseptr nothing;

  i = dispchar;
  nothing = (baseptr)Malloc(chars*sizeof(long));
  nothing[i - 1] = 0;
  bottom = true;
  firstrav(root, i);
  hyptrav(root, nothing, i,&bottom);
  free(nothing);
}  /* hypstates */


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
  Char st;
  chartype c, d;

  p = nuroot;
  q = nuroot;
  extra = false;
  if (i == p->ycoord && (p == root || subtree)) {
    c = over;
    if (display) {
      switch (p->state) {
	
      case 'A':
        c = aa;
        break;
	
      case 'C':
        c = cc;
        break;
	
      case 'G':
        c = gg;
        break;
	
      case 'T':
        c = tt;
        break;
	
      case '?':
        c = quest;
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
      if (display) {
        switch (q->state) {
	
        case 'A':
          c = aa;
          break;
	
        case 'C':
          c = cc;
          break;
	
        case 'G':
          c = gg;
          break;
	
        case 'T':
          c = tt;
          break;
	
        case '?':
          c = quest;
          break;
        }
        d = c;
      }
      if (n > 1) {
        makechar(d);
        prefix(c);
        for (j = 1; j <= n - 2; j++)
          putchar(chh[(long)c]);
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
          st = p->next->back->state;
        else
          st = p->next->next->back->state;
        if (display) {
          switch (st) {
	
          case 'A':
            c = aa;
            break;
	
          case 'C':
            c = cc;
            break;
	
          case 'G':
            c = gg;
            break;
	
          case 'T':
            c = tt;
            break;
	
          case '?':
            c = quest;
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
      printf("  SITE%4hd", dispchar);
    else
      printf("           ");
    if (subtree)
      printf(" Subtree ");
    else
      printf("         ");
    printf("(unrooted)");
  }
  if (i == 3 && display) {
    printf("   ");
    makechar(aa);
    printf(":A, ");
    makechar(cc);
    printf(":C, ");
    makechar(gg);
    printf(":G, ");
    makechar(tt);
    printf(":T, ");
    makechar(quest);
    printf(":?");
  }
  if ((i == 5 || (lastline < 5 && i == lastline)) && !earlytree) {
    if (lastline < 5)
      printf("\n                   ");
    printf("   Steps:%10.5f", - like);
    gotlike = -like;
  }
  if ((i == 7 || (lastline < 7 && i == lastline)) && changed && !earlytree) {
    if (lastline < 7)
      printf("\n                     ");
    if (-like <bestyet) {
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
    printf("     %3hd sites compatible", compatible);
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
  printf((ansi || ibmpc) ? "\033[2J\033[H" :
	 vt52 ? "\033E\033H"    : "\n");
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
    drawline(i, tipy-dow);
  if (spp <= screenlines / 2 || subtree)
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
      printf("\n before node (type number): ");
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
        printf("\nERROR IN USER TREE: UNMATCHED PARENTHESIS OR NOT BIFURCATED NODE\n");
	exit(-1);
      } else {
        if (ch == ')')
          done = true;
      }
    } else if (c == ';') {
      if (ch != ';') {
        printf("\nERROR IN USER TREE: UNMATCHED PARENTHESIS OR MISSING SEMICOLON\n");
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

void addelement(p, nextnode,lparens)
     node **p;
     short *nextnode,*lparens;
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
  } while (ch == ' ' );
  if (ch == '(' ) {
    if (*lparens>= spp - 1) {
      printf("\nERROR IN USER TREE: TOO MANY LEFT PARENTHESES\n");
      exit(-1);
    }
    (*nextnode)++;
    (*lparens)++;
    q = treenode[*nextnode - 1];
    addelement(&q->next->back, nextnode,lparens);
    q->next->back->back = q->next;
    findch(',');
    addelement(&q->next->next->back, nextnode,lparens);
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
      if (intree[n - 1] == false) {
        *p = treenode[n - 1];
        intree[n - 1] = true;
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
  short nextnode, lparens;
  short i;
  openfile(&treefile,trfilename,"r",progname,trfilename);
  root = treenode[spp];
  nextnode = spp;
  root->back = NULL;
  for (i = 0; i < (spp); i++)
    intree[i] = false;
  lparens = 0;
  addelement(&root, &nextnode,&lparens);
  findch(';');
  printf("\n\n");
  FClose(treefile);
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
  if (!outgropt)
    outgrno = root->next->back->index;
  if (outgropt && intree[outgrno - 1])
    reroot(treenode[outgrno - 1]);
}  /* buildtree */

void setorder()
{
  /* sets in order of number of members */
  sett[0] = 1L << ((long)A);
  sett[1] = 1L << ((long)C);
  sett[2] = 1L << ((long)G);
  sett[3] = 1L << ((long)U);
  sett[4] = 1L << ((long)O);
  sett[5] = (1L << ((long)A)) | (1L << ((long)C));
  sett[6] = (1L << ((long)A)) | (1L << ((long)G));
  sett[7] = (1L << ((long)A)) | (1L << ((long)U));
  sett[8] = (1L << ((long)A)) | (1L << ((long)O));
  sett[9] = (1L << ((long)C)) | (1L << ((long)G));
  sett[10] = (1L << ((long)C)) | (1L << ((long)U));
  sett[11] = (1L << ((long)C)) | (1L << ((long)O));
  sett[12] = (1L << ((long)G)) | (1L << ((long)U));
  sett[13] = (1L << ((long)G)) | (1L << ((long)O));
  sett[14] = (1L << ((long)U)) | (1L << ((long)O));
  sett[15] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G));
  sett[16] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)U));
  sett[17] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)O));
  sett[18] = (1L << ((long)A)) | (1L << ((long)G)) | (1L << ((long)U));
  sett[19] = (1L << ((long)A)) | (1L << ((long)G)) | (1L << ((long)O));
  sett[20] = (1L << ((long)A)) | (1L << ((long)U)) | (1L << ((long)O));
  sett[21] = (1L << ((long)C)) | (1L << ((long)G)) | (1L << ((long)U));
  sett[22] = (1L << ((long)C)) | (1L << ((long)G)) | (1L << ((long)O));
  sett[23] = (1L << ((long)C)) | (1L << ((long)U)) | (1L << ((long)O));
  sett[24] = (1L << ((long)G)) | (1L << ((long)U)) | (1L << ((long)O));
  sett[25] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G)) |
    (1L << ((long)U));
  sett[26] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G)) |
    (1L << ((long)O));
  sett[27] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)U)) |
    (1L << ((long)O));
  sett[28] = (1L << ((long)A)) | (1L << ((long)G)) | (1L << ((long)U)) |
    (1L << ((long)O));
  sett[29] = (1L << ((long)C)) | (1L << ((long)G)) | (1L << ((long)U)) |
    (1L << ((long)O));
  sett[30] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G)) |
    (1L << ((long)U)) | (1L << ((long)O));
}  /* setorder */

void mincomp()
{
  /* computes for each site the minimum number of steps
     necessary to accomodate those species already
     in the analysis */
  short i, j, k;
  boolean done;

  for (i = 0; i < (chars); i++) {
    done = false;
    j = 0;
    while (!done) {
      j++;
      done = true;
      k = 1;
      do {
        if (intree[k - 1])
          done = (done && (treenode[k - 1]->base[i] & sett[j - 1]) != 0);
        k++;
      } while (k <= spp && done);
    }
    if (j == 31)
      necsteps[i] = 4;
    if (j <= 30)
      necsteps[i] = 3;
    if (j <= 25)
      necsteps[i] = 2;
    if (j <= 15)
      necsteps[i] = 1;
    if (j <= 5)
      necsteps[i] = 0;
    necsteps[i] *= weight[i];
  }
}  /* mincomp */


void help()
{
  /* display help information */
  printf("\n\nR Rearrange a tree by moving a node or group\n");
  printf("# Show the states of the next site that doesn't fit tree\n");
  printf("+             ... of the next site\n");
  printf("-             ... of the previous site\n");
  printf("S Show the states of a given site\n");
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

Local Void nextinc()
{
  /* show next incompatible site */
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
  } while (!(necsteps[dispchar - 1] != numsteps[dispchar - 1] ||
             dispchar == disp0 || done));
  printree();
}  /* nextinc */

void nextchar()
{
  /* show next site */
  display = true;
  dispchar++;
  if (dispchar > chars)
    dispchar = 1;
  printree();
}  /* nextchar */

void prevchar()
{
  /* show previous site */
  display = true;
  dispchar--;
  if (dispchar < 1)
    dispchar = chars;
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
    }
    if (ok && i == 0)
      display = false;
  } while (!ok);
  printree();
}  /* show */

/* Local variables for addpreorder: */

void tryadd(p,item,nufork,place)
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

void addpreorder(p, item_, nufork_,place)
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
  tryadd(p,&item,&nufork,place);
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
  if (root == treenode[treenode[i - 1]->back->index - 1]) {
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
  add(treenode[fromwhere - 1], treenode[what - 1], dummy);
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
    treenode[atwhat - 1]->next->back =treenode[atwhat - 1]->next->next->back;
    treenode[atwhat - 1]->next->next->back = q;
    treenode[atwhat - 1]->next->back->back = treenode[atwhat - 1]->next;
    treenode[atwhat - 1]->next->next->back->back =
      treenode[atwhat - 1]->next->next;
    break;

  case reroott:
    restoring = true;
    temp =oldoutgrno;
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

  if (waswritten) {
    printf("\nTree file already was open.\n");
    printf("   A   Add to this tree to tree file\n");
    printf("   R   Replace tree file contents by this tree\n");
    printf("   F   Write out tree to a different tree file\n");
    printf("   N   Do Not write out this tree\n");
    do {
      printf("Which should we do? ");
      scanf("%c%*[^\n]", &ch);
      getchar();
      uppercase(&ch);
    } while (ch != 'A' && ch != 'R' && ch != 'N'  && ch != 'F');
  }
  if (ch == 'F'){
    trfilename[0] = '\0';
    while (trfilename[0] =='\0'){
      printf("Please enter a tree file name>");
      gets(trfilename);}
  }

  if (ch == 'R' || ch == 'A' || !waswritten){
    openfile(&treefile,trfilename, (ch == 'A' && waswritten) ? "a" : "w",
	     progname,trfilename);
  }
  if (!done)
    printree();
  if (waswritten && ch != 'A' && ch != 'R')
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
}
/* treewrite */

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
  boolean done = false;
  waswritten = false;
  do {
    printf("NEXT? (Options: R # + - S . T U W O F C H ? X Q) ");
    printf("(H or ? for Help) ");
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    if (strchr("CFORSTUXQ+#-.WH?",ch)){
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
  if (written)
    return;
  do {
    printf("Do you want to write out the tree to a file? (Y or N) ");
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == 'Y' || ch == 'y')
      treewrite(done);
  } while (ch != 'Y' && ch != 'y' && ch != 'N' && ch != 'n');
}  /* redisplay */

void treeconstruct()
{
  /* constructs a binary tree from the pointers in treenode. */

  restoring = false;
  subtree = false;
  display = false;
  dispchar = 0;
  earlytree = true;
  waswritten = false;
  buildtree();
  printf("\nComputing steps needed for compatibility in sites ...\n\n");
  setorder();
  mincomp();
  newtree = true;
  earlytree = false;
  printree();
  bestyet = -like;
  gotlike = -like;
  lastop = none;
  newtree = false;
  written = false;
  redisplay();
}  /* treeconstruct */

main(argc, argv)
     int argc;
     Char *argv[];

{  /* Interactive DNA parsimony */
  /* reads in spp, chars, and the data. Then calls treeconstruct to
     construct the tree and query the user */
#ifdef MAC
  macsetup("Dnamove","");
  argv[0] = "Dnamove";
#endif
  strcpy(progname,argv[0]);
  strcpy(infilename,INFILE);
  strcpy(outfilename,OUTFILE);
  strcpy(trfilename,TREEFILE);

  openfile(&infile,infilename,"r",argv[0],infilename);
  openfile(&outfile,outfilename,"w",argv[0],outfilename);

  garbage = NULL;
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
  FClose(outfile);
#ifdef MAC
  fixmacfile(outfilename);
#endif
  exit(0);
}  /* Interactive DNA parsimony */

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

