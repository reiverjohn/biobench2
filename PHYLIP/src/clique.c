#include "phylip.h"

/* version 3.52c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Jerry Shurman, Hisashi Horino,
   Akiko Fuseki, Sean Lamont, and Andrew Keeffe.  Permission is granted
   to copy and use this program provided no fee is charged for it and
   provided that this copyright notice is not removed. */

#define NmLngth         10   /* number of characters in a name */
#define FormWide        80   /* width of outfile page */

#define ibmpc0          false
#define ansi0           true
#define vt520           false
#define down            2
typedef boolean *aPtr;
typedef long *SpPtr;
typedef long *ChPtr;

typedef struct vecrec {
  aPtr vec;
  struct vecrec *next;
} vecrec;

typedef Char **NamePtr;
typedef vecrec **aDataPtr;
typedef vecrec **Matrix;      /* nodes will form a binary tree         */
typedef struct node {         /* describes a tip species or an ancestor*/
  struct node *next, *back;   /* pointers to nodes                     */
  long index;                 /* number of the node                    */
  long maxpos;                /*                                       */
  boolean tip;                /* present species are tips of tree      */
  long nodeset;               /* used by accumulate                    */
  long xcoord, ycoord, ymin;  /* used by printree                      */
  long ymax;
} node;

typedef node **pointptr;
char infilename[100],outfilename[100],trfilename[100];

Static long NumSpp, NumChars, ActualChars, Cliqmin, outgrno, col, datasets,
	     ith, j, setsz;
Static boolean ancvar, Clmin, Factors, outgropt, trout, weights,
	       noroot, printdata, printcomp, progress, treeprint, mulsets,
               ibmpc, vt52, ansi, firstset;
Static NamePtr Nayme;
Static aPtr ancone;
Static Char *Factor;
Static long *ActChar, *oldweight, *weight;
Static aDataPtr Data;
Static Matrix Comp;            /* the character compatibility matrix      */
Static FILE *infile, *outfile, *treefile;
Static node *root;
Static long **grouping;
Static pointptr treenode;   /* pointers to all nodes in tree              */
Static vecrec *garbage;

/* these variables are local to DoAll in the pascal Version. */
Static  aPtr aChars;
Static  boolean *Rarer;
Static  long n, MaxChars;
Static  SpPtr SpOrder;
Static  ChPtr ChOrder;

/* Local variables for GetMaxCliques: */
Static  vecrec **Comp2;
Static  long tcount;
Static  aPtr Temp, Processed, Rarer2;

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

Static Void gnu(p)
vecrec **p;
{  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */

 if (garbage != NULL) {
    *p = garbage;
    garbage = garbage->next;
  } else {
    *p = (vecrec *)Malloc((long)sizeof(vecrec));
    (*p)->vec = (aPtr)Malloc((long)NumChars*sizeof(boolean));
  }
  (*p)->next = NULL;
}  /* gnu */


Static Void chuck(p)
vecrec *p;
{  /* collect garbage on p -- put it on front of garbage list */

 p->next = garbage;
  garbage = p;
}  /* chuck */


void nunode(p)
node **p;
{  /* replacement for NEW */
  *p = (node *)Malloc((long)sizeof(node));
  (*p)->next = NULL;
  (*p)->tip = false;
}  /* nunode */


void NewLine(i, j)
long i, j;
{
  /* goes to new line if i MOD j is zero */
  long k;

  if (i % j != 0)
    return;
  putc('\n', outfile);
  for (k = 1; k <= NmLngth + 1; k++)
    putc(' ', outfile);
}  /* NewLine */


void uppercase(ch)
Char *ch;
{  /* convert a character to upper case -- either ASCII or EBCDIC */
   *ch = (islower(*ch) ?  toupper(*ch) : (*ch));
}  /* uppercase */

void inputnumbers()
{
  /* set variables */
  fscanf(infile, "%ld%ld", &NumSpp, &NumChars);
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld character states\n\n",
	    NumSpp, NumChars);
}  /* inputnumbers */

void getoptions()
{
  /* interactively set options */
  Char ch;
  boolean done, done1;

  fprintf(outfile, "\nLargest clique program, version %s\n\n",VERSION);
  putchar('\n');
  ancvar = false;
  Clmin = false;
  Factors = false;
  outgrno = 1;
  outgropt = false;
  trout = true;
  weights = false;
  printdata = false;
  printcomp = false;
  progress = true;
  treeprint = true;
  do {
    printf((ansi) ? "\033[2J\033[H" :
           (vt52) ? "\033E\033H"    : "\n");
    printf("\nLargest clique program, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  A   Use ancestral states in input file?  %s\n",
	   (ancvar ? "Yes" : "No"));
    printf("  C          Specify minimum clique size?");
    if (Clmin)
      printf("  Yes, at size%3ld\n", Cliqmin);
    else
      printf("  No\n");
    printf("  O                        Outgroup root?  %s%3ld\n",
	   (outgropt ? "Yes, at species number" :
                       "No, use as outgroup species"),outgrno);
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
	   (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
	   (progress ? "Yes" : "No"));
    printf("  3        Print out compatibility matrix  %s\n",
	   (printcomp ? "Yes" : "No"));
    printf("  4                        Print out tree  %s\n",
	   (treeprint ? "Yes" : "No"));
    printf("  5       Write out trees onto tree file?  %s\n",
	   (trout ? "Yes" : "No"));
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      if (strchr("OACM012345",ch)){
	switch (ch) {

	case 'A':
	  ancvar = !ancvar;
	  break;

	case 'C':
	  Clmin = !Clmin;
	  if (Clmin) {
	    do {
	      printf("Minimum clique size:\n");
	      scanf("%ld%*[^\n]", &Cliqmin);
	      getchar();
	    } while (Cliqmin < 0);
	  }
	  break;

	case 'O':
	  outgropt = !outgropt;
	  if (outgropt) {
	    done1 = true;
	    do {
	      printf("Type number of the outgroup:\n");
	      scanf("%ld%*[^\n]", &outgrno);
	      getchar();
	      done1 = (outgrno >= 1 && outgrno <= NumSpp);
	      if (!done1) {
		printf("BAD OUTGROUP NUMBER: %4ld\n", outgrno);
		printf("  Must be in range 1 -%2ld\n", NumSpp);
	      }
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
	  printcomp = !printcomp;
	  break;

	case '4':
	  treeprint = !treeprint;
	  break;

	case '5':
	  trout = !trout;
	  break;
	}
      } else
	printf("Not a possible option!\n");
    }
  } while (!done);
}  /* getoptions */


Static Void setuptree()
{
  /* initialization of tree pointers, variables */
  long i;

  treenode = (pointptr)Malloc((long)NumSpp*sizeof(node *));
  for (i = 0; i < NumSpp; i++) {
    treenode[i] = (node *)Malloc((long)sizeof(node));
    treenode[i]->next = NULL;
    treenode[i]->back = NULL;
    treenode[i]->index = i + 1;
    treenode[i]->tip = false;
  }
}  /* setuptree */

Static Void doinit()
{
  /* initializes variables */
  long i;

  inputnumbers();
  getoptions();
  setuptree();
  Data = (aDataPtr)Malloc((long)NumSpp*sizeof(vecrec *));
  for (i = 0; i < (NumSpp); i++)
    gnu(&Data[i]);
  Comp = (Matrix)Malloc((long)NumChars*sizeof(vecrec *));
  for (i = 0; i < (NumChars); i++)
    gnu(&Comp[i]);
}  /* doinit */


Local Void inputancestors()
{
  /* reads the ancestral states for each character */
  long i;
  Char ch;

  for (i = 1; i < NmLngth; i++)
    ch = getc(infile);
  for (i = 0; i < (NumChars); i++) {
    do {
      if (eoln(infile)) {
	fscanf(infile, "%*[^\n]");
	getc(infile);
      }
      ch = getc(infile);
    } while (ch == ' ');
    if (ch == '0' || ch == '1') {
      switch (ch) {

      case '1':
	ancone[i] = true;
	break;

      case '0':
	ancone[i] = false;
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

Local Void printancestors()
{
  /* print out list of ancestral states */
  long i;

  fprintf(outfile, "Ancestral states:\n");
  for (i = 1; i <= NmLngth + 2; i++)
    putc(' ', outfile);
  for (i = 1; i <= (NumChars); i++) {
    NewLine(i, 55);
    if (ancone[i - 1])
      putc('1', outfile);
    else
      putc('0', outfile);
    if (i % 5 == 0)
      putc(' ', outfile);
  }
  fprintf(outfile, "\n\n");
}  /* printancestor */

Local Void inputfactors()
{
  /* reads the factor symbols */
  long i;
  Char ch;

  ActualChars = 1;
  for (i = 1; i < NmLngth; i++)
    ch = getc(infile);
  for (i = 1; i <= (NumChars); i++) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    Factor[i - 1] = getc(infile);
    if (i > 1) {
      if (Factor[i - 1] != Factor[i - 2])
	ActualChars++;
    }
    ActChar[i - 1] = ActualChars;
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  Factors = true;
}  /* inputfactors */

Local Void printfactors()
{
  /* print out list of factor symbols */
  long i;

  fprintf(outfile, "Factors:");
  for (i = 1; i <= NmLngth - 6; i++)
    putc(' ', outfile);
  for (i = 1; i <= (NumChars); i++) {
    NewLine(i, 50);
    putc(Factor[i - 1], outfile);
    if (i % 5 == 0)
      putc(' ', outfile);
  }
  fprintf(outfile, "\n\n");

}  /* printfactors */

Local Void inputweights()
{
  /* input the character weights, 0-9 and A-Z for weights 0 - 35 */
  Char ch;
  long i;

  for (i = 1; i < NmLngth; i++)
    ch = getc(infile);
  for (i = 0; i < (NumChars); i++) {
    do {
      if (eoln(infile)) {
	fscanf(infile, "%*[^\n]");
	getc(infile);
      }
      ch = getc(infile);
    } while (ch == ' ');
    oldweight[i] = 1;
    if (isdigit(ch))
      oldweight[i] = ch - '0';
    else if (isalpha(ch)) {
      uppercase(&ch);
      if (ch >= 'A' && ch <= 'I')
	oldweight[i] = ch - 55;
      else if (ch >= 'J' && ch <= 'R')
	oldweight[i] = ch - 55;
      else
	oldweight[i] = ch - 55;
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
  /* print out the weights of characters */
  long i, j, k;

  fprintf(outfile, "   Characters are weighted as follows:\n");
  fprintf(outfile, "       ");
  for (i = 0; i <= 9; i++)
    fprintf(outfile, "%3ld", i);
  fprintf(outfile, "\n     *---------------------------------\n");
  for (j = 0; j <= (ActualChars / 10); j++) {
    fprintf(outfile, "%5ld! ", j * 10);
    for (i = 0; i <= 9; i++) {
      k = j * 10 + i;
      if (k > 0 && k <= ActualChars)
	fprintf(outfile, "%3ld", oldweight[k - 1]);
      else
	fprintf(outfile, "   ");
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* printweights */


Static Void ReadData()
{
  /* reads the species names and character data */
  long i, j, Extranum;
  Char ch;
  boolean avar;
  long cursp, curchs;

  if (!firstset) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%ld%ld", &cursp, &curchs);
    if (cursp != NumSpp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4ld\n", ith);
      exit(-1);
    }
    NumChars = curchs;
  }
  avar = false;
  ActualChars = NumChars;
  for (i = 1; i <= (NumChars); i++)
    ActChar[i - 1] = i;
  for (i = 0; i < (NumChars); i++)
    oldweight[i] = 1;
  Extranum = 0;
  while (!(eoln(infile))) {
    ch = getc(infile);
    uppercase(&ch);
    if (ch == 'A' || ch == 'F' || ch == 'W')
      Extranum++;
    else if (ch != ' ') {
      printf("BAD OPTION CHARACTER: %c\n", ch);
      putc('\n', outfile);
      exit(-1);
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  for (i = 1; i <= Extranum; i++) {
    ch = getc(infile);
    uppercase(&ch);
    if (ch != 'A' && ch != 'F' && ch != 'W') {
      printf("ERROR: INCORRECT AUXILIARY OPTIONS LINE");
      printf(" WHICH STARTS WITH %c\n", ch);
      }
    if (ch == 'A') {
      avar = true;
      if (!ancvar) {
	printf("ERROR: ANCESTOR OPTION NOT CHOSEN IN MENU");
        printf(" WITH OPTION %c IN INPUT\n",ch);
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
    printf("ERROR: ANCESTOR OPTION CHOSEN IN MENU");
    printf(" WITH NO OPTION A IN INPUT\n");
    exit(-1);
  }
  if (weights && printdata)
    printweights();
  if (Factors)
    printfactors();
  if (ancvar && avar && printdata)
    printancestors();
  noroot = !(outgropt || (ancvar && avar));
  j = NumChars / 2 + (NumChars / 5 - 1) / 2 - 5;
  if (j < 0)
    j = 0;
  if (j > 27)
    j = 27;
  if (printdata) {
    fprintf(outfile, "Species  ");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "Character states\n");
    fprintf(outfile, "-------  ");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "--------- ------\n\n");
  }
  for (i = 0; i < (NumSpp); i++) {
    for (j = 0; j < NmLngth; j++) {
      if (eof(infile) || eoln(infile)){
	printf("ERROR: END-OF-LINE OR END-OF-FILE IN");
	printf(" THE MIDDLE OF A SPECIES NAME\n");
	exit(-1);}
      Nayme[i][j] = getc(infile);
      if (printdata)
	putc(Nayme[i][j], outfile);
    }
    if (printdata)
      fprintf(outfile, "  ");
    for (j = 1; j <= (NumChars); j++) {
      do {
	if (eoln(infile)) {
	  fscanf(infile, "%*[^\n]");
	  getc(infile);
	}
	ch = getc(infile);
      } while (ch == ' ');
      if (printdata) {
	putc(ch, outfile);
	NewLine(j, 55);
	if (j % 5 == 0)
	  putc(' ', outfile);
      }
      if (ch != '0' && ch != '1') {
	printf("BAD CHARACTER STATE: %c (NOT 0 OR 1)");
	printf("AT CHARACTER%3ld OF SPECIES%3ld\n", ch, j, i + 1);
	exit(-1);
      }
      Data[i]->vec[j - 1] = (ch == '1');
    }
    fscanf(infile, "%*[^\n]");
    getc(infile);
    if (printdata)
      putc('\n', outfile);
  }
  putc('\n', outfile);
  for (i = 0; i < (NumChars); i++) {
    if (i + 1 == 1 || !Factors)
      weight[i] = oldweight[i];
    else if (Factor[i] != Factor[i - 1])
      weight[ActChar[i] - 1] = oldweight[i];
  }
}  /* ReadData */


Local boolean Compatible(ch1, ch2)
long ch1, ch2;
{
  /* TRUE if two characters ch1 < ch2 are compatible */
  long i, j, k;
  boolean Compt, Done1, Done2;
  boolean Info[4];

  Compt = true;
  j = 1;
  while (ch1 > ActChar[j - 1])
    j++;
  Done1 = (ch1 != ActChar[j - 1]);
  while (!Done1) {
    k = j;
    while (ch2 > ActChar[k - 1])
      k++;
    Done2 = (ch2 != ActChar[k - 1]);
    while (!Done2) {
      for (i = 0; i <= 3; i++)
	Info[i] = false;
      if (ancvar) {
	if (ancone[j - 1] && ancone[k - 1])
	  Info[0] = true;
	else if (ancone[j - 1] && !ancone[k - 1])
	  Info[1] = true;
	else if (!ancone[j - 1] && ancone[k - 1])
	  Info[2] = true;
	else
	  Info[3] = true;
      }
      for (i = 0; i < (NumSpp); i++) {
	if (Data[i]->vec[j - 1] && Data[i]->vec[k - 1])
	  Info[0] = true;
	else if (Data[i]->vec[j - 1] && !Data[i]->vec[k - 1])
	  Info[1] = true;
	else if (!Data[i]->vec[j - 1] && Data[i]->vec[k - 1])
	  Info[2] = true;
	else
	  Info[3] = true;
      }
      Compt = (Compt && !(Info[0] && Info[1] && Info[2] && Info[3]));
      k++;
      Done2 = (k > NumChars);
      if (!Done2)
	Done2 = (ch2 != ActChar[k - 1]);
    }
    j++;
    Done1 = (j > NumChars);
    if (!Done1)
      Done1 = (ch1 != ActChar[j - 1]);
  }
  return Compt;
}  /* Compatible */


Static Void SetUp(Comp)
vecrec **Comp;
{
  /* sets up the compatibility matrix */
  long i, j;

  if (printcomp) {
    if (Factors)
      fprintf(outfile, "      (For original multistate characters)\n");
    fprintf(outfile, "Character Compatibility Matrix (1 if compatible)\n");
    fprintf(outfile, "--------- ------------- ------ -- -- -----------\n\n");
  }
  for (i = 0; i < (ActualChars); i++) {
    if (printcomp) {
      for (j = 1; j <= ((48 - ActualChars) / 2); j++)
	putc(' ', outfile);
      for (j = 1; j < i + 1; j++) {
	if (Comp[i]->vec[j - 1])
	  putc('1', outfile);
	else
	  putc('.', outfile);
	NewLine(j, 70);
      }
    }
    Comp[i]->vec[i] = true;
    if (printcomp)
      putc('1', outfile);
    for (j = i + 1; j < (ActualChars); j++) {
      Comp[i]->vec[j] = Compatible(i + 1, j + 1);
      if (printcomp) {
	if (Comp[i]->vec[j])
	  putc('1', outfile);
	else
	  putc('.', outfile);
      }
      Comp[j]->vec[i] = Comp[i]->vec[j];
    }
    if (printcomp)
      putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* Setup */





Local Void Intersect(V1, V2, V3)
boolean *V1, *V2, *V3;
{
  /* takes the logical intersection V1 AND V2 */
  long i;

  for (i = 0; i < (ActualChars); i++)
    V3[i] = (V1[i] && V2[i]);
}  /* Intersect */

Local long CountStates(V)
boolean *V;
{
  /* counts the 1's in V */
  long i, TempCount;

  TempCount = 0;
  for (i = 0; i < (ActualChars); i++) {
    if (V[i])
      TempCount += weight[i];
  }
  return TempCount;
}  /* CountStates */

Local Void Gen1(i, CurSize, aChars, Candidates, Excluded)
long i, CurSize;
boolean *aChars, *Candidates, *Excluded;
{
  /* finds largest size cliques and prints them out */
  long CurSize2, j, k, Actual, Possible;
  boolean Futile;
  vecrec *Chars2, *Cands2, *Excl2, *Cprime, *Exprime;

  gnu(&Chars2);
  gnu(&Cands2);
  gnu(&Excl2);
  gnu(&Cprime);
  gnu(&Exprime);
  CurSize2 = CurSize;
  memcpy(Chars2->vec, aChars, NumChars*sizeof(boolean));
  memcpy(Cands2->vec, Candidates, NumChars*sizeof(boolean));
  memcpy(Excl2->vec, Excluded, NumChars*sizeof(boolean));
  j = i;
  while (j <= ActualChars) {
    if (Cands2->vec[j - 1]) {
      Chars2->vec[j - 1] = true;
      Cands2->vec[j - 1] = false;
      CurSize2 += weight[j - 1];
      Possible = CountStates(Cands2->vec);
      Intersect(Cands2->vec, Comp2[j - 1]->vec, Cprime->vec);
      Actual = CountStates(Cprime->vec);
      Intersect(Excl2->vec, Comp2[j - 1]->vec, Exprime->vec);
      Futile = false;
      for (k = 0; k <= j - 2; k++) {
	if (Exprime->vec[k] && !Futile) {
	  Intersect(Cprime->vec, Comp2[k]->vec, Temp);
	  Futile = (CountStates(Temp) == Actual);
	}
      }
      if (CurSize2 + Actual >= Cliqmin && !Futile) {
	if (Actual > 0)
	  Gen1(j + 1,CurSize2,Chars2->vec,Cprime->vec,Exprime->vec);
	else if (CurSize2 > Cliqmin) {
	  Cliqmin = CurSize2;
	  if (tcount >= 0)
	    tcount = 1;
	} else if (CurSize2 == Cliqmin)
	  tcount++;
      }
      if (Possible > Actual) {
	Chars2->vec[j - 1] = false;
	Excl2->vec[j - 1] = true;
	CurSize2 -= weight[j - 1];
      } else
	j = ActualChars;
    }
    j++;
  }
  chuck(Chars2);
  chuck(Cands2);
  chuck(Excl2);
  chuck(Cprime);
  chuck(Exprime);
}  /* Gen1 */


Local boolean Ingroupstate(i)
long i;
{
  /* the ingroup state for the i-th character */
  boolean outstate;

  if (noroot) {
    outstate = Data[0]->vec[i - 1];
    return (!outstate);
  }
  if (ancvar)
    outstate = ancone[i - 1];
  else
    outstate = Data[outgrno - 1]->vec[i - 1];
  return (!outstate);
}  /* Ingroupstate */

Local Void makeset()
{
  /* make up set of species for given set of characters */
  long i, j, k, l, m;
  boolean instate;
  long *st;

  st = (long *)Malloc(setsz*sizeof(long));
  n = 0;
  for (i = 0; i < (MaxChars); i++) {
    for (j = 0; j < setsz; j++)
      st[j] = 0;
    instate = Ingroupstate(ChOrder[i]);
    for (j = 0; j < (NumSpp); j++) {
      k = (long)((j+1)/SETBITS);
      if (Data[SpOrder[j] - 1]->vec[ChOrder[i] - 1] == instate) {
        m = (long)(SpOrder[j]/SETBITS);
	st[m] = ((long)st[m]) | (1L << (SpOrder[j] % SETBITS));
      }
    }
    memcpy(grouping[++n - 1], st, setsz*sizeof(long));
  }
  for (i = 0; i < (NumSpp); i++) {
    k = (long)(SpOrder[i]/SETBITS);
    grouping[++n - 1][k] = 1L << (SpOrder[i] % SETBITS);
  }
  free(st);
}  /* makeset */

Local Void Init(ChOrder, Count, MaxChars,aChars)
long *ChOrder, *Count;
long *MaxChars;
aPtr aChars;
{
  /* initialize vectors and character count */
  long i, j, temp;
  boolean instate;

  *MaxChars = 0;
  for (i = 1; i <= (NumChars); i++) {
    if (aChars[ActChar[i - 1] - 1]) {
      (*MaxChars)++;
      ChOrder[*MaxChars - 1] = i;
      instate = Ingroupstate(i);
      temp = 0;
      for (j = 0; j < (NumSpp); j++) {
	if (Data[j]->vec[i - 1] == instate)
	  temp++;
      }
      Count[i - 1] = temp;
    }
  }
}  /*Init */

Local Void ChSort(ChOrder, Count, MaxChars)
long *ChOrder, *Count;
long MaxChars;
{
  /* sorts the characters by number of ingroup states */
  long j, temp;
  boolean ordered;

  ordered = false;
  while (!ordered) {
    ordered = true;
    for (j = 1; j < MaxChars; j++) {
      if (Count[ChOrder[j - 1] - 1] < Count[ChOrder[j] - 1]) {
	ordered = false;
	temp = ChOrder[j - 1];
	ChOrder[j - 1] = ChOrder[j];
	ChOrder[j] = temp;
      }
    }
  }
}  /* ChSort */

Local Void PrintClique(aChars)
boolean *aChars;
{
  /* prints the characters in a clique */
  long i, j;
  fprintf(outfile, "\n\n");
  if (Factors) {
    fprintf(outfile, "Actual Characters: (");
    j = 0;
    for (i = 1; i <= (ActualChars); i++) {
      if (aChars[i - 1]) {
	fprintf(outfile, "%3ld", i);
	j++;
	NewLine(j, (long)((FormWide - 22) / 3));
      }
    }
    fprintf(outfile, ")\n");
  }
  if (Factors)
    fprintf(outfile, "Binary ");
  fprintf(outfile, "Characters: (");
  j = 0;
  for (i = 1; i <= (NumChars); i++) {
    if (aChars[ActChar[i - 1] - 1]) {
      fprintf(outfile, "%3ld", i);
      j++;
      if (Factors)
	NewLine(j, (long)((FormWide - 22) / 3));
      else
	NewLine(j, (long)((FormWide - 15) / 3));
    }
  }
  fprintf(outfile, ")\n\n");
}  /* PrintClique */


Local Void bigsubset(st, n)
long *st;
long n;
{
  /* find a maximal subset of st among the groupings */
  long i, j;
  long *su;
  boolean max, same;

  su = (long *)Malloc(setsz*sizeof(long));
  for (i = 0; i < setsz; i++)
    su[i] = 0;
  for (i = 0; i < n; i++) {
    max = true;
    for (j = 0; j < setsz; j++)
      if ((grouping[i][j] & ~st[j]) != 0)
       max = false;
    if (max) {
      same = true;
      for (j = 0; j < setsz; j++)
        if (grouping[i][j] != st[j])
          same = false;
      if (!same) {
        for (j = 0; j < setsz; j++)
          if ((su[j] & ~grouping[i][j]) != 0)
            max = false;
        if (max) {
          same = true;
          for (j = 0; j < setsz; j++)
            if (grouping[i][j] != su[j])
              same = false;
          if (!same)
            memcpy(su, grouping[i], setsz*sizeof(long));
        }
      }
    }
  }
  memcpy(st, su, setsz*sizeof(long));
  free(su);
}  /* bigsubset */

Local Void recontraverse(p, st, n, MaxChars)
node **p;
long *st;
long n;
long MaxChars;
{
  /* traverse to reconstruct the tree from the characters */
  long i, j, k, maxpos;
  long *tempset, *st2;
  boolean found, zero, zero2, same;
  node *q;

  j = k = 0;
  for (i = 1; i <= (NumSpp); i++) {
    if (((1L << (i % SETBITS)) & st[(long)(i / SETBITS)]) != 0) {
      k++;
      j = i;
    }
  }
  if (k == 1) {
    *p = treenode[j - 1];
    (*p)->tip = true;
    (*p)->index = j;
    return;
  }
  nunode(p);
  (*p)->index = 0;
  tempset = (long*)Malloc(setsz*sizeof(long));
  memcpy(tempset, st, setsz*sizeof(long));
  q = *p;
  zero = true;
  for (i = 0; i < setsz; i++)
    if (tempset[i] != 0)
      zero = false;
  if (!zero)
    bigsubset(tempset, n);
  zero = true;
  zero2 = true;
  for (i = 0; i < setsz; i++)
    if (st[i] != 0)
      zero = false;
  if (!zero) {
    for (i = 0; i < setsz; i++)
      if (tempset[i] != 0)
        zero2 = false;
  }
  st2 = (long *)Malloc(setsz*sizeof(long));
  memcpy(st2, st, setsz*sizeof(long));
  while (!zero2) {
    nunode(&q->next);
    q = q->next;
    recontraverse(&q->back, tempset, n,MaxChars);
    i = 1;
    maxpos = 0;
    while (i <= MaxChars) {
      same = true;
      for (j = 0; j < setsz; j++)
        if (grouping[i - 1][j] != tempset[j])
          same = false;
      if (same)
	maxpos = i;
      i++;
    }
    q->back->maxpos = maxpos;
    q->back->back = q;
    for (j = 0; j < setsz; j++)
      st2[j] &= ~tempset[j];
    memcpy(tempset, st2, setsz*sizeof(long));
    found = false;
    i = 1;
    while (!found && i <= n) {
      same = true;
      for (j = 0; j < setsz; j++)
        if (grouping[i - 1][j] != tempset[j])
          same = false;
      if (same)
        found = true;
      else
	i++;
    }
    zero = true;
    for (j = 0; j < setsz; j++)
      if (tempset[j] != 0)
        zero = false;
    if (!zero && !found)
      bigsubset(tempset, n);
    zero = true;
    zero2 = true;
    for (j = 0; j < setsz; j++)
      if (st2[j] != 0)
        zero = false;
    if (!zero)
      for (j = 0; j < setsz; j++)
        if (tempset[j] != 0)
          zero2 = false;
}
  q->next = *p;
  free(tempset);
  free(st2);
}  /* recontraverse */

Local Void reconstruct(n,MaxChars)
long n,MaxChars;
{  /* reconstruct tree from the subsets */
  long i;
  long *s;
  s = (long *)Malloc(setsz*sizeof(long));
  for (i = 0; i < setsz; i++) {
    if (i+1 == setsz) {
      s[i] = 1L << ((NumSpp % SETBITS) + 1);
      if (setsz > 1)
        s[i] -= 1;
      else s[i] -= 1L << 1;
    }
    else if (i == 0) {
      if (setsz > 1)
        s[i] = ~0L - 1;
    }
    else {
      if (setsz > 2)
        s[i] = ~0L;
    }
  }
  recontraverse(&root,s,n,MaxChars);
  free(s);
}  /* reconstruct */

Local Void reroot(outgroup)
node *outgroup;
{
  /* reorients tree, putting outgroup in desired position. */
  long i;
  boolean nroot;
  node *p, *q;

  nroot = false;
  p = root->next;
  while (p != root) {
    if (outgroup->back == p) {
      nroot = true;
      p = root;
    } else
      p = p->next;
  }
  if (nroot)
    return;
  p = root;
  i = 0;
  while (p->next != root) {
    p = p->next;
    i++;
  }
  if (i == 2) {
    root->next->back->back = p->back;
    p->back->back = root->next->back;
    q = root->next;
  } else {
    p->next = root->next;
    nunode(&root->next);
    q = root->next;
    nunode(&q->next);
    p = q->next;
    p->next = root;
    q->tip = false;
    p->tip = false;
  }
  q->back = outgroup;
  p->back = outgroup->back;
  outgroup->back->back = p;
  outgroup->back = q;
}  /* reroot */

Local Void coordinates(p, tipy,MaxChars)
node *p;
long *tipy;
long MaxChars;
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;
  long maxx;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
    return;
  }
  q = p->next;
  maxx = 0;
  while (q != p) {
    coordinates(q->back,tipy,MaxChars);
    if (!q->back->tip) {
      if (q->back->xcoord > maxx)
	maxx = q->back->xcoord;
    }
    q = q->next;
  }
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (MaxChars - p->maxpos) * 3 - 2;
  if (p == root)
    p->xcoord += 2;
  p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */

Local Void drawline(i)
long i;
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q;
  long n, m, j, k, l, sumlocpos, size, locpos, branchpos;
  long *poslist;
  boolean extra, done, plus, found, same;
  node *r, *first, *last;

  poslist = (long *)Malloc(((long)NumSpp + MaxChars)*sizeof(long));
  branchpos = 0;
  p = root;
  q = root;
  fprintf(outfile, "  ");
  extra = false;
  plus = false;
  do {
    if (!p->tip) {
      found = false;
      r = p->next;
      while (r != p && !found) {
	if (i >= r->back->ymin && i <= r->back->ymax) {
	  q = r->back;
	  found = true;
	} else
	  r = r->next;
      }
      first = p->next->back;
      r = p;
      while (r->next != p)
	r = r->next;
      last = r->back;
    }
    done = (p->tip || p == q);
    n = p->xcoord - q->xcoord;
    m = n;
    if (extra) {
      n--;
      extra = false;
    }
    if (q->ycoord == i && !done) {
      if (!q->tip) {
	putc('+', outfile);
	plus = true;
	j = 1;
	for (k = 1; k <= (q->maxpos); k++) {
          same = true;
          for (l = 0; l < setsz; l++)
            if (grouping[k - 1][l] != grouping[q->maxpos - 1][l])
              same = false;
          if (same) {
	    poslist[j - 1] = k;
	    j++;
	  }
	}
	size = j - 1;
	if (size == 0) {
	  for (k = 1; k < n; k++)
	    putc('-', outfile);
	  sumlocpos = n;
	} else {
	  sumlocpos = 0;
	  j = 1;
	  while (j <= size) {
	    locpos = poslist[j - 1] * 3;
	    if (j != 1)
	      locpos -= poslist[j - 2] * 3;
	    else
	      locpos -= branchpos;
	    for (k = 1; k < locpos; k++)
	      putc('-', outfile);
	    if (Rarer[ChOrder[poslist[j - 1] - 1] - 1])
	      putc('1', outfile);
	    else
	      putc('0', outfile);
	    sumlocpos += locpos;
	    j++;
	  }
	  for (j = sumlocpos + 1; j < n; j++)
	    putc('-', outfile);
	  putc('+', outfile);
	  if (m > 0)
	    branchpos += m;
	  extra = true;
	}
      } else {
	if (!plus) {
	  putc('+', outfile);
	  plus = false;
	} else
	  n++;
	j = 1;
	for (k = 1; k <= (q->maxpos); k++) {
          same = true;
          for (l = 0; l < setsz; l++)
 	    if (grouping[k - 1][l] != grouping[q->maxpos - 1][l])
              same = false;
          if (same) {
	    poslist[j - 1] = k;
	    j++;
	  }
	}
	size = j - 1;
	if (size == 0) {
	  for (k = 1; k <= n; k++)
	    putc('-', outfile);
	  sumlocpos = n;
	} else {
	  sumlocpos = 0;
	  j = 1;
	  while (j <= size) {
	    locpos = poslist[j - 1] * 3;
	    if (j != 1)
	      locpos -= poslist[j - 2] * 3;
	    else
	      locpos -= branchpos;
	    for (k = 1; k < locpos; k++)
	      putc('-', outfile);
	    if (Rarer[ChOrder[poslist[j - 1] - 1] - 1])
	      putc('1', outfile);
	    else
	      putc('0', outfile);
	    sumlocpos += locpos;
	    j++;
	  }
	  for (j = sumlocpos + 1; j <= n; j++)
	    putc('-', outfile);
	  if (m > 0)
	    branchpos += m;
	}
	putc('-', outfile);
      }
    } else if (!p->tip && last->ycoord > i && first->ycoord < i &&
	       (i != p->ycoord || p == root)) {
      putc('!', outfile);
      for (j = 1; j < n; j++)
	putc(' ', outfile);
      plus = false;
      if (m > 0)
	branchpos += m;
    } else {
      for (j = 1; j <= n; j++)
	putc(' ', outfile);
      plus = false;
      if (m > 0)
	branchpos += m;
    }
    if (q != p)
      p = q;
  } while (!done);
  if (p->ycoord == i && p->tip) {
    for (j = 0; j < NmLngth; j++)
      putc(Nayme[p->index - 1][j], outfile);
}
  putc('\n', outfile);
  free(poslist);
}  /* drawline */

Local Void printree()
{
  /* prints out diagram of the tree */
  long tipy;
  long i;

  if (!treeprint)
    return;
  tipy = 1;
  coordinates(root, &tipy,MaxChars);
  fprintf(outfile, "\n  Tree and");
  if (Factors)
    fprintf(outfile, " binary");
  fprintf(outfile, " characters:\n\n");
  fprintf(outfile, "   ");
  for (i = 0; i < (MaxChars); i++)
    fprintf(outfile, "%3ld", ChOrder[i]);
  fprintf(outfile, "\n   ");
  for (i = 0; i < (MaxChars); i++) {
    if (Rarer[ChOrder[i] - 1])
      fprintf(outfile, "%3c", '1');
    else
      fprintf(outfile, "%3c", '0');
  }
  fprintf(outfile, "\n\n");
  for (i = 1; i <= (tipy - down); i++)
    drawline(i);
  fprintf(outfile, "\nremember: this is an unrooted tree!\n\n");
}  /* printree */

Local Void treeout(p, tcount)
node *p;
long tcount;
{
  /* write out file with representation of final tree */
  long i, n;
  Char c;
  node *q;

  if (p->tip) {
    n = 0;
    for (i = 1; i <= NmLngth; i++) {
      if (Nayme[p->index - 1][i - 1] != ' ')
	n = i;
    }
    for (i = 0; i < n; i++) {
      c = Nayme[p->index - 1][i];
      if (c == ' ')
	c = '_';
      putc(c, treefile);
    }
    col += n;
  } else {
    q = p->next;
    putc('(', treefile);
    col++;
    while (q != p) {
      treeout(q->back, tcount);
      q = q->next;
      if (q == p)
	break;
      putc(',', treefile);
      col++;
      if (col > 72) {
	col = 0;
	putc('\n', treefile);
      }
    }
    putc(')', treefile);
    col++;
  }
  if (p != root)
    return;
  putc(';', treefile);
  if (tcount <= 1)
    putc('\n', treefile);
  else
    fprintf(treefile, "%7.4llf\n", 1.0 / tcount);
}  /* treeout */

Local Void DoAll(Chars_, Processed, Rarer_, tcount)
boolean *Chars_, *Processed, *Rarer_;
long tcount;
{
  /* print out a clique and its tree */
  long i, j;
  ChPtr Count;

  aChars = (aPtr)Malloc((long)NumChars*sizeof(boolean));
  SpOrder = (SpPtr)Malloc((long)NumSpp*sizeof(long));
  ChOrder = (ChPtr)Malloc((long)NumChars*sizeof(long));
  Count = (ChPtr)Malloc((long)NumChars*sizeof(long));
  memcpy(aChars, Chars_, NumChars*sizeof(boolean));
  Rarer = Rarer_;
  Init(ChOrder, Count, &MaxChars, aChars);
  ChSort(ChOrder, Count, MaxChars);
  for (i = 1; i <= (NumSpp); i++)
    SpOrder[i - 1] = i;
  for (i = 1; i <= (NumChars); i++) {
    if (aChars[ActChar[i - 1] - 1]) {
      if (!Processed[ActChar[i - 1] - 1]) {
	Rarer[i - 1] = Ingroupstate(i);
	Processed[ActChar[i - 1] - 1] = true;
      }
    }
  }
  PrintClique(aChars);
  grouping = (long **)Malloc(((long)(NumSpp + MaxChars))*sizeof(long *));
  for (i = 0; i < NumSpp + MaxChars; i++) {
    grouping[i] = (long *)Malloc(setsz*sizeof(long));
    for (j = 0; j < setsz; j++)
      grouping[i][j] = 0;
  }
  makeset();
  setuptree();
  reconstruct(n,MaxChars);
  if (noroot)
    reroot(treenode[outgrno - 1]);
  printree();
  if (trout) {
    col = 0;
    treeout(root, tcount);
  }
  free(SpOrder);
  free(ChOrder);
  free(Count);
  for (i = 0; i < NumSpp + MaxChars; i++)
    free(grouping[i]);
  free(grouping);
}  /* DoAll */

Local Void Gen2(i, CurSize, aChars, Candidates, Excluded)
long i, CurSize;
boolean *aChars, *Candidates, *Excluded;
{
  /* finds largest size cliques and prints them out */
  long CurSize2, j, k, Actual, Possible;
  boolean Futile;
  vecrec *Chars2, *Cands2, *Excl2, *Cprime, *Exprime;

  gnu(&Chars2);
  gnu(&Cands2);
  gnu(&Excl2);
  gnu(&Cprime);
  gnu(&Exprime);
  CurSize2 = CurSize;
  memcpy(Chars2->vec, aChars, NumChars*sizeof(boolean));
  memcpy(Cands2->vec, Candidates, NumChars*sizeof(boolean));
  memcpy(Excl2->vec, Excluded, NumChars*sizeof(boolean));
  j = i;
  while (j <= ActualChars) {
    if (Cands2->vec[j - 1]) {
      Chars2->vec[j - 1] = true;
      Cands2->vec[j - 1] = false;
      CurSize2 += weight[j - 1];
      Possible = CountStates(Cands2->vec);
      Intersect(Cands2->vec, Comp2[j - 1]->vec, Cprime->vec);
      Actual = CountStates(Cprime->vec);
      Intersect(Excl2->vec, Comp2[j - 1]->vec, Exprime->vec);
      Futile = false;
      for (k = 0; k <= j - 2; k++) {
	if (Exprime->vec[k] && !Futile) {
	  Intersect(Cprime->vec, Comp2[k]->vec, Temp);
	  Futile = (CountStates(Temp) == Actual);
	}
      }
      if (CurSize2 + Actual >= Cliqmin && !Futile) {
	if (Actual > 0)
	  Gen2(j + 1,CurSize2,Chars2->vec,Cprime->vec,Exprime->vec);
	else
	  DoAll(Chars2->vec,Processed,Rarer2,tcount);
      }
      if (Possible > Actual) {
	Chars2->vec[j - 1] = false;
	Excl2->vec[j - 1] = true;
	CurSize2 -= weight[j - 1];
      } else
	j = ActualChars;
    }
    j++;
  }
  chuck(Chars2);
  chuck(Cands2);
  chuck(Excl2);
  chuck(Cprime);
  chuck(Exprime);
}  /* Gen2 */


Static Void GetMaxCliques(Comp_)
vecrec **Comp_;
{
  /* recursively generates the largest cliques */
  long i;
  aPtr aChars, Candidates, Excluded;

  Temp = (aPtr)Malloc((long)NumChars*sizeof(boolean));
  Processed = (aPtr)Malloc((long)NumChars*sizeof(boolean));
  Rarer2 = (aPtr)Malloc((long)NumChars*sizeof(boolean));
  aChars = (aPtr)Malloc((long)NumChars*sizeof(boolean));
  Candidates = (aPtr)Malloc((long)NumChars*sizeof(boolean));
  Excluded = (aPtr)Malloc((long)NumChars*sizeof(boolean));
  Comp2 = Comp_;
  putc('\n', outfile);
  if (Clmin) {
    fprintf(outfile, "Cliques with at least%3ld characters\n", Cliqmin);
    fprintf(outfile, "------- ---- -- ----- -- ----------\n");
  } else {
    Cliqmin = 0;
    fprintf(outfile, "Largest Cliques\n");
    fprintf(outfile, "------- -------\n");
    for (i = 0; i < (ActualChars); i++) {
      aChars[i] = false;
      Excluded[i] = false;
      Candidates[i] = true;
    }
    tcount = 0;
    Gen1(1, 0, aChars, Candidates, Excluded);
  }
  for (i = 0; i < (ActualChars); i++) {
    aChars[i] = false;
    Candidates[i] = true;
    Processed[i] = false;
    Excluded[i] = false;
  }
  Gen2(1, 0, aChars, Candidates, Excluded);
  putc('\n', outfile);
  free(Temp);
  free(Processed);
  free(Rarer2);
  free(aChars);
  free(Candidates);
  free(Excluded);
}  /* GetMaxCliques */


main(argc, argv)
int argc;
Char *argv[];
{  /* Main Program */
char infilename[100],outfilename[100],trfilename[100];
#ifdef MAC
   macsetup("Clique","");
   argv[0] = "Clique";
#endif
  openfile(&infile,INFILE,"r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);
  openfile(&treefile,TREEFILE,"w",argv[0],trfilename);
  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  mulsets = false;
  firstset = true;
  datasets = 1;
  doinit();
  setsz = (short)ceil(((double)NumSpp+1.0)/(double)SETBITS);
  ancone = (aPtr)Malloc((long)NumChars*sizeof(boolean));
  Factor = (Char *)Malloc((long)NumChars*sizeof(Char));
  ActChar = (long *)Malloc((long)NumChars*sizeof(long));
  oldweight = (long *)Malloc((long)NumChars*sizeof(long));
  weight = (long *)Malloc((long)NumChars*sizeof(long));
  Nayme = (Char **)Malloc((long)NumSpp*sizeof(Char *));
  for (j = 0; j < NumSpp; j++)
    Nayme[j] = (Char *)Malloc((long)NmLngth*sizeof(Char));
  for (ith = 1; ith <= (datasets); ith++) {
    ReadData();
    firstset = false;
    SetUp(Comp);
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n\n",ith);
      if (progress)
        printf("\nData set # %ld:\n",ith);
    }
    GetMaxCliques(Comp);
    if (progress) {
      printf("\nOutput written to output file\n");
      if (trout)
        printf("\nTree(s) written on tree file\n\n");
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
}


int eof(f)
FILE *f;
{
    register long ch;

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
    register long ch;

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

