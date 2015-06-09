#include "phylip.h"

/* version 3.52c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define namelength      10   /* maximum number of characters in species name */

#define ibmpc0          false
#define ansi0           true
#define vt520           false

#define point           "."

typedef double *phenotype;
typedef Char naym[namelength];
typedef struct node {
  struct node *next, *back;
  boolean tip;
  long number;
  phenotype view;
  naym nayme;
  double v, deltav;
} node;

typedef struct tree {
  node **nodep;
  node *start;
} tree;

Static FILE *infile, *outfile, *treefile;
Static long numsp, numsp1, numsp2, chars, numtrees, j;
Static phenotype *x, *cntrast;
Static naym *nayms;
Static boolean printdata, progress, reg, mulsets, ibmpc, vt52, ansi,
       bifurcating;
Static Char ch;

/* Local variables for maketree, propagated globally for c version: */
tree curtree;
long contno;


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
        printf("Please enter a new filename>");
        gets(file);
        break;
      case 'w':
        printf("%s: can't write %s\n",application,file);
        printf("Please enter a new filename>");
        gets(file);
        break;
      }
    }
  }
  *fp=of;
  if (perm != NULL)
    strcpy(perm,file);
}



void uppercase(ch)
Char *ch;
{
  /* convert a character to upper case -- either ASCII or EBCDIC */
    *ch = isupper(*ch) ? (*ch) :toupper(*ch);
}  /* uppercase */


void getnums()
{
  /* read species numbers and number of characters */
  fscanf(infile, "%ld%ld", &numsp, &chars);
  numsp1 = numsp + 1;
  numsp2 = numsp * 2 - 1;
}  /* getnums */

void getoptions()
{
  /* interactively set options */
  Char ch;
  boolean done, done1;

  printdata = false;
  progress = true;
  do {
    printf(ansi ? "\033[2J\033[H" :
           vt52 ? "\033E\033H"    : "\n");
    printf("\nContinuous character Contrasts, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  R     Print out correlations and regressions?  %s\n",
	   (reg ? "Yes" : "No"));
    printf("  M                     Analyze multiple trees?");
    if (mulsets)
      printf("  Yes, %2ld trees\n", numtrees);
    else
      printf("  No\n");
    printf("  0         Terminal type (IBM PC, VT52, ANSI)?  %s\n",
	   ibmpc ? "IBM PC"  :
	   ansi  ? "ANSI"    :
	   vt52  ? "VT52"    : "(none)");

    printf("  1          Print out the data at start of run  %s\n",
	   (printdata ? "Yes" : "No"));
    printf("  2        Print indications of progress of run  %s\n",
	   (progress ? "Yes" : "No"));

    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      if (ch == 'R' || ch == 'M' || ch == '1' || ch == '2' || ch == '0') {
	switch (ch) {

	case 'R':
	  reg = !reg;
	  break;

	case 'M':
	  mulsets = !mulsets;
	  if (mulsets) {
	    done1 = false;
	    do {
	      printf("How many trees?\n");
	      scanf("%ld%*[^\n]", &numtrees);
	      getchar();
	      done1 = (numtrees >= 1);
	      if (!done1)
		printf("BAD TREES NUMBER:  it must be greater than 1\n");
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
	}
      } else
	printf("Not a possible option!\n");
    }
  } while (!done);
}  /* getoptions */

void getdata()
{
  /* read species data */
  long i, j, l;

  if (printdata) {
    fprintf(outfile, "\nContinuous character Contrasts, version %s\n\n",VERSION);
    fprintf(outfile, "%4ld Populations, %4ld Characters\n\n", numsp, chars);
    fprintf(outfile, "Name");
    fprintf(outfile, "                       Phenotypes\n");
    fprintf(outfile, "----");
    fprintf(outfile, "                       ----------\n\n");
  }
  for (i = 0; i < numsp; i++) {
    fscanf(infile, "%*[^\n]");
    getc(infile);
    for (j = 0; j < namelength; j++) {
      nayms[i][j] = getc(infile);
      if (nayms[i][j] == '\n')
	nayms[i][j] = ' ';
      if ( eof(infile) | eoln(infile)){
	printf("ERROR: END-OF-LINE OR END-OF-FILE");
        printf(" IN THE MIDDLE OF A SPECIES NAME\n");
	exit(-1);}
      else if (printdata)
	putc(nayms[i][j], outfile);
    }
    for (j = 1; j <= chars; j++) {
      if (eoln(infile)) {
	fscanf(infile, "%*[^\n]");
	getc(infile);
      }
      fscanf(infile, "%lf", &x[i][j - 1]);
      if (printdata) {
	fprintf(outfile, "%10.5lf", x[i][j - 1]);
	if (j % 6 == 0) {
	  putc('\n', outfile);
	  for (l = 1; l <= namelength; l++)
	    putc(' ', outfile);
	}
      }
    }
    if (printdata)
      putc('\n', outfile);
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  if (printdata)
    putc('\n', outfile);
}  /* getdata */


void doinit()
{
  /* initializes variables */
  getnums();
  getoptions();
}  /* doinit */


void setuptree(a)
tree *a;
{
  /* initialize a tree */
  long i, j;
  node *p, *q;

  a->nodep = (node **)Malloc((long)numsp2*sizeof(node *));
  for (i = 1; i <= numsp; i++) {
    a->nodep[i - 1] = (node *)Malloc((long)sizeof(node));
    a->nodep[i - 1]->view = (phenotype)Malloc((long)chars*sizeof(double));
    a->nodep[i - 1]->tip = true;
    a->nodep[i - 1]->number = i;
    a->nodep[i - 1]->back = NULL;
  }
  for (i = numsp1; i <= numsp2; i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc((long)sizeof(node));
      p->view = (phenotype)Malloc((long)chars*sizeof(double));
      p->tip = false;
      p->number = i;
      p->next = q;
      p->back = NULL;
      q = p;
    }
    p->next->next->next = p;
    a->nodep[i - 1] = p;
  }
}  /* setuptree */

void hookup(p, q)
node *p, *q;
{
  /* hook together two nodes */
  p->back = q;
  q->back = p;
}  /* hookup */

void setuptip(m, t)
long m;
tree *t;
{
  /* initialize branch lengths and views in a tip */
  node *with;

  with = t->nodep[m - 1];
  memcpy(with->view, x[m - 1], chars*sizeof(double));
  memcpy(with->nayme, nayms[m - 1], sizeof(naym));
  with->deltav = 0.0;
  with->v = 0.0;
}  /* setuptip */


void getch(c)
Char *c;
{
  /* get next nonblank character */
  do {
    if (eoln(treefile)) {
      fscanf(treefile, "%*[^\n]");
      getc(treefile);
    }
    *c = getc(treefile);
    if (*c == '\n')
      *c = ' ';
  } while (*c == ' ');
}  /* getch */

void findch(c, lparens,rparens)
Char c;
long *lparens,*rparens;
{
  /* skip forward in user tree until find character c */
  boolean done;

  done = false;
  while (!done) {
    if (c == ',') {
      if (ch == '(' || ch == ')' || ch == ':' || ch == ';') {
	printf("\nERROR IN USER TREE: UNMATCHED PARENTHESIS OR MISSING");
        printf(" COMMA\n OR NOT TRIFURCATED BASE\n");
	exit(-1);
      } else if (ch == ',')
	done = true;
    } else if (c == ')') {
      if (ch == '(' || ch == ',' || ch == ':' || ch == ';') {
	printf("\nERROR IN USER TREE:");
	printf(" UNMATCHED PARENTHESIS OR NOT BIFURCATED NODE\n");
	exit(-1);
      } else if (ch == ')') {
	(*rparens)++;
	if ((*lparens) > 0 && (*lparens) == (*rparens)) {
	  if ((*lparens) < numsp - 2) {
	    printf("\nERROR: UNMATCHED PARENTHESIS OR TOO FEW SPECIES\n");
	    exit(-1);
	  } else if ((*lparens) == numsp - 1) {
	    getch(&ch);
	    if (ch != ';') {
	      printf("\nERROR IN USER TREE: ");
	      printf("UNMATCHED PARENTHESIS OR MISSING SEMICOLON\n");
	      exit(-1);
	    }
	  }
	}
	done = true;
      }
    }
    if ((done && ch == ')') || (!done))
      getch(&ch);
  }
}  /* findch */

void findeither(lparens,rparens,rtparen)
long *lparens,*rparens;
boolean *rtparen;
{
  /* find either a rt paren or a comma */
  boolean done;

  done = false;
  while (!(done)) {
    if (ch == '(' || ch == ':' || ch == ';') {
      printf("\nERROR IN USER TREE: UNMATCHED PARENTHESIS OR MISSING COMMA\n");
      printf(" OR NOT TRIFURCATED BASE\n");
      exit(-1);
    } else if (ch == ',' || ch == ')')
      done = true;
    else
      getch(&ch);
  }
  if (ch == ')') {
    (*rparens)++;
    if ((*lparens) > 0 && (*lparens) == (*rparens)) {
      if ((*lparens) < numsp - 2) {
	printf("\nERROR: UNMATCHED PARENTHESIS OR TOO FEW SPECIES\n");
	exit(-1);
      } else if ((*lparens) == numsp - 2) {
	getch(&ch);
	if (ch != ';') {
	  printf("\nERROR IN USER TREE:");
          printf(" UNMATCHED PARENTHESIS OR MISSING SEMICOLON\n");
	  exit(-1);
	}
      }
    }
  }
  (*rtparen) = (ch == ')');
  if ((done && ch == ')') || !(done))
    getch(&ch);
}  /* findeither */


void processlength(p)
node *p;
{
  long digit, ordzero;
  double valyew, divisor;
  boolean pointread;
  Char STR1[256], STR2[256];

  ordzero = '0';
  pointread = false;
  valyew = 0.0;
  divisor = 1.0;
  getch(&ch);
  digit = ch - ordzero;
  while (((unsigned long)digit <= 9) |
	 (strcmp((sprintf(STR1, "%c", ch), STR1), point) == 0)) {
    sprintf(STR2, "%c", ch);
    if (!strcmp(STR2, point))
      pointread = true;
    else {
      valyew = valyew * 10.0 + digit;
      if (pointread)
	divisor *= 10.0;
    }
    getch(&ch);
    digit = ch - ordzero;
  }
  p->v = valyew / divisor;
  p->back->v = p->v;
}  /* processlength */


void addelement(p, nextnode,lparens,rparens,names)
node *p;
long *nextnode,*lparens,*rparens;
boolean *names;
{
  /* add one node to the user tree */
  node *q;
  long i, n;
  boolean found;
  Char str[namelength];

  getch(&ch);
  if (ch == '(') {
    (*lparens)++;
    if ((*lparens) >= numsp) {
      printf("\nERROR IN USER TREE: TOO MANY LEFT PARENTHESES\n");
      exit(-1);
    } else {
      (*nextnode)++;
      q = curtree.nodep[(*nextnode) - 1];
      hookup(p, q);
      addelement(q->next, nextnode,lparens,rparens,names);
      findch(',', lparens,rparens);
      addelement(q->next->next, nextnode,lparens,rparens,names);
      findch(')',lparens,rparens);
    }
  } else {
    for (i = 0; i < namelength; i++)
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
      if (ch == '\n')
	ch = ' ';
      n++;
    } while (ch != ':' && ch != ',' && ch != ')' && n <= namelength);
    n = 1;
    do {
      found = true;
      for (i = 0; i < namelength; i++)
	found = (found && str[i] == nayms[n - 1][i]);
      if (found) {
	if (names[n - 1] == false)
	  names[n - 1] = true;
	else {
	  printf("\nERROR IN USER TREE: DUPLICATE NAME FOUND -- ");
	  for (i = 0; i < namelength; i++)
	    putchar(curtree.nodep[n - 1]->nayme[i]);
	  putchar('\n');
	  exit(-1);
	}
      } else
	n++;
    } while (!(n > numsp || found));
    if (n > numsp) {
      printf("Cannot find species: ");
      for (i = 0; i < namelength; i++)
	putchar(str[i]);
      putchar('\n');
      exit(-1);
    }
    hookup(curtree.nodep[n - 1], p);
  }
  if (ch == ':')
    processlength(p);
}  /* addelement */

void treeread()
{
  /* read in a user tree */
  /* Local variables for treeread: */
  boolean rtparen;
  long nextnode, lparens, rparens;
  boolean *names;
  node *p;
  long i;

  getch(&ch);
  if (ch != '(')
    return;
  names = (boolean *)Malloc((long)numsp*sizeof(boolean));
  for (i = 0; i < numsp; i++)
    names[i] = false;
  lparens = 1;
  rparens = 0;
  nextnode = numsp + 1;
  p = curtree.nodep[nextnode - 1];
  p->back = NULL;
  curtree.start = p;
  p = p->next;
  i = 1;
  rtparen = false;
  while (i <= 3 && !rtparen) {
    addelement(p, &nextnode,&lparens,&rparens,names);
    p = p->next;
    findeither(&lparens,&rparens,&rtparen);
    i++;
  }
  bifurcating = (i == 3);
  free(names);
}  /* treeread */

void initcontrasts()
{
  /* compute the contrasts */
  long i, j;
  for (i = 0; i <= numsp - 2 ; i++) {
    for (j = 0; j < chars; j++)
      cntrast[i][j] = 0.0;
  }
  contno = 1;
}  /* initcontrasts */

void contbetween(p, q)
node *p, *q;
{
  /* compute one contrast */
  long i;
  double sqbl;

  if (p->back != q)
    sqbl = sqrt(p->v + q->v + p->deltav + q->deltav);
  else
    sqbl = sqrt(p->v + p->deltav + q->deltav);
  for (i = 0; i < chars; i++)
    cntrast[contno - 1][i] = (p->view[i] - q->view[i]) / sqbl;
  contno++;
}  /* contbetween */

void nuview(p)
node *p;
{
  /* renew information about subtrees */
  long j;
  node *q, *r;
  double v1, v2, vtot, f1, f2;

  q = p->next->back;
  r = p->next->next->back;
  v1 = q->v + q->deltav;
  v2 = r->v + r->deltav;
  vtot = v1 + v2;
  if (vtot > 0.0)
    f1 = v2 / vtot;
  else
    f1 = 0.5;
  f2 = 1.0 - f1;
  for (j = 0; j < chars; j++)
    p->view[j] = f1 * q->view[j] + f2 * r->view[j];
  p->deltav = v1 * f1;
}  /* nuview */

void makecontrasts(p)
node *p;
{
  /* compute the contrasts, recursively */
  if (p->tip)
    return;
  makecontrasts(p->next->back);
  makecontrasts(p->next->next->back);
  nuview(p);
  contbetween(p->next->back, p->next->next->back);
}  /* makecontrasts */

void writecontrasts()
{
  /* write out the contrasts */
  long i, j;

  if (printdata || reg) {
    fprintf(outfile, "\nContrasts (columns are different characters)\n");
    fprintf(outfile, "--------- -------- --- --------- -----------\n\n");
  }
  for (i = 0; i <= contno - 2; i++) {
    for (j = 0; j < chars; j++)
      fprintf(outfile, "%10.5lf", cntrast[i][j]);
    putc('\n', outfile);
  }
}  /* writecontrasts */

void regressions()
{
  /* compute regressions and correlations among contrasts */
  long i, j, k;
  double **sumprod;

  sumprod = (double **)Malloc((long)chars*sizeof(double *));
  for (i = 0; i < chars; i++) {
    sumprod[i] = (double *)Malloc((long)chars*sizeof(double));
    for (j = 0; j < chars; j++)
      sumprod[i][j] = 0.0;
  }
    for (i = 0; i <= contno - 2; i++) {
    for (j = 0; j < chars; j++) {
      for (k = 0; k < chars; k++)
	sumprod[j][k] += cntrast[i][j] * cntrast[i][k];
    }
  }
  fprintf(outfile, "\nCovariance matrix\n");
  fprintf(outfile, "---------- ------\n\n");
  for (i = 0; i < chars; i++) {
    for (j = 0; j < chars; j++)
      sumprod[i][j] /= contno - 1;
  }
  for (i = 0; i < chars; i++) {
    for (j = 0; j < chars; j++)
      fprintf(outfile, "%10.4lf", sumprod[i][j]);
    putc('\n', outfile);
  }
  fprintf(outfile, "\nRegressions (columns on rows)\n");
  fprintf(outfile, "----------- -------- -- -----\n\n");
  for (i = 0; i < chars; i++) {
    for (j = 0; j < chars; j++)
      fprintf(outfile, "%10.4lf", sumprod[i][j] / sumprod[i][i]);
    putc('\n', outfile);
  }
  fprintf(outfile, "\nCorrelations\n");
  fprintf(outfile, "------------\n\n");
  for (i = 0; i < chars; i++) {
    for (j = 0; j < chars; j++)
      fprintf(outfile, "%10.4lf",
	      sumprod[i][j] / sqrt(sumprod[i][i] * sumprod[j][j]));
    putc('\n', outfile);
  }
  for (i = 0; i < chars; i++)
    free(sumprod[i]);
  free(sumprod);

}  /* regressions */


void maketree()
{
  /* set up the tree and use it */
  long which;

  setuptree(&curtree);
  for (which = 1; which <= numsp; which++)
    setuptip(which, &curtree);
  which = 1;
  while (which <= numtrees) {
    if ((printdata || reg) && numtrees > 1) {
      fprintf(outfile, "Tree number%4ld\n", which);
      fprintf(outfile, "==== ====== ====\n\n");
    }
    bifurcating = false;
    treeread();
    initcontrasts();
    makecontrasts(curtree.start);
    if (!bifurcating)
      contbetween(curtree.start, curtree.start->back);
    writecontrasts();
    if (reg)
      regressions();
    putc('\n', outfile);
    which++;
  }
  if (progress)
    printf("\nOutput written on output file\n\n");
}  /* maketree */


main(argc, argv)
int argc;
Char *argv[];
{  /* main program */
char infilename[100],outfilename[100],trfilename[100];
#ifdef MAC
  macsetup("Contrast","");
  argv[0] = "Contrast";
#endif
  openfile(&infile,INFILE,"r",argv[0],infilename);
  openfile(&treefile,TREEFILE,"r",argv[0],trfilename);
  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);
  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  mulsets = false;
  reg = true;
  numtrees = 1;
  doinit();
  x = (phenotype *)Malloc((long)numsp*sizeof(phenotype));
  for (j = 0; j < numsp; j++)
    x[j] = (phenotype)Malloc((long)chars*sizeof(double));
  cntrast = (phenotype *)Malloc((long)numsp*sizeof(phenotype));
  for (j = 0; j < numsp; j++)
    cntrast[j] = (phenotype)Malloc((long)chars*sizeof(double));
  nayms = (naym *)Malloc((long)numsp*sizeof(naym));
  getdata();
  maketree();
  FClose(infile);
  FClose(outfile);
  FClose(treefile);
  exit(0);
}

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
mem = (MALLOCRETURN *)calloc(1,(size_t)x);
if (!mem)
  memerror();
else
  return (MALLOCRETURN *)mem;

}
