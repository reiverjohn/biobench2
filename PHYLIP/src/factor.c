#include "phylip.h"

/* version 3.56c. (c) Copyright 1988-1993 by Joseph Felsenstein and
   Christopher Meacham.
   A program to factor multistate character trees.
   Originally programmed 29 May 1983 by C. A. Meacham, Botany Department,
     University of Georgia
   Additional code by Joe Felsenstein, 1988-1991
   C version code by Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define nmlngth         10   /* Length of species name                       */
#define maxstates       20   /* Maximum number of states in multi chars      */
#define maxoutput       80   /* Maximum length of output line                */
#define sizearray       5000 /* Size of symbarray; must be >= the sum of     */
                             /* squares of the number of states in each multi*/
                             /* char to be factored                          */
#define factchar        ':'  /* character to indicate state connections      */
#define unkchar         '?'  /* input character to indicate state unknown    */

#define ibmpc0          false /* true if terminal is an IBM PC screen        */
#define ansi0           true  /* true if terminal is an ANSI terminal        */
#define vt520           false /* true if terminal is a VT52  terminal        */

typedef struct node {     /* Node of multifurcating tree */
  struct node *ancstr, *sibling, *descendant;
  Char state;             /* Symbol of character state   */
  short edge;             /* Number of subtending edge   */
} node;

Static FILE *infile, *outfile;
Static short neus, nchars, charindex, lastindex;
Static Char ch;
Static boolean ancstrrequest, factorrequest, rooted,  ibmpc, vt52,
	       ansi, progress;
Static Char symbarray[sizearray];
 /* Holds multi symbols and their factored equivs        */
Static short *charnum;     /* Multis           */
Static short *chstart;     /* Position of each */
Static short *numstates;   /* Number of states */
Static Char  *ancsymbol;   /* Ancestral state  */

/*  local variables for dotrees, propagated to global level. */
  short npairs, offset, charnumber, nstates;
  node *root;
  Char pair[maxstates][2];
  node *nodes[maxstates];


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

void uppercase(ch)
Char *ch;
{  /* convert a character to upper case -- either ASCII or EBCDIC */
   *ch = (islower(*ch) ?  toupper(*ch) : (*ch));
}  /* uppercase */

void getoptions()
{
  /* interactively set options */
  Char ch;
  boolean done;

  ibmpc = ibmpc0;
  vt52 = vt520;
  ansi = ansi0;
  progress = true;
  factorrequest = false;
  ancstrrequest = false;
  putchar('\n');
  for (;;){
    printf(ansi ? "\033[2J\033[H" :
	   vt52 ? "\033E\033H"    : "\n");
    printf("\nFactor -- multistate to binary recoding program, version %s\n\n"
	   ,VERSION);
    printf("Settings for this run:\n");
    printf("  A      put ancestral states in output file?  %s\n",
	   ancstrrequest ? "Yes" : "No");
    printf("  F   put factors information in output file?  %s\n",
	   factorrequest ? "Yes" : "No");
    printf("  0       Terminal type (IBM PC, VT52, ANSI)?  %s\n",
	   ibmpc ? "IBM PC" :
	   ansi  ? "ANSI"   :
	   vt52  ? "VT52"   : "(none)");
    printf("  1      Print indications of progress of run  %s\n",
	   (progress ? "Yes" : "No"));
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (ch == 'A' || ch == 'F' || ch == '0' || ch == '1') {
      switch (ch) {
	
      case 'A':
	ancstrrequest = !ancstrrequest;
	break;
	
      case 'F':
	factorrequest = !factorrequest;
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
        progress = !progress;
      break;
      }
    } else
      printf("Not a possible option!\n");
  }
}  /* getoptions */

void nextch(ch)
Char *ch;
{
  *ch = ' ';
  while (*ch == ' ' && !eoln(infile))
    *ch = getc(infile);
}  /* nextch */

void readtree()
{
  /* Reads a single character-state tree; puts adjacent symbol
     pairs into array 'pairs' */

  npairs = 0;
  while (!eoln(infile)) {
    nextch(&ch);
    if (eoln(infile))
      break;
    npairs++;
    pair[npairs - 1][0] = ch;
    nextch(&ch);
    if (eoln(infile) || (ch != factchar)) {
      printf("CHARACTER%4hd:  ERROR IN CHARACTER STATE TREE FORMAT\n",
	     charnumber);
      exit(-1);}

    nextch(&pair[npairs - 1][1]);
    if (eoln(infile) && pair[npairs - 1][1] == ' '){
      printf("CHARACTER%4hd:  ERROR IN CHARACTER STATE TREE FORMAT\n",
	     charnumber);
      exit(-1);}
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);

}  /* readtree */

void attachnodes(poynter,otherone)
node *poynter;
Char *otherone;
{
  /* Makes linked list of all nodes to which passed node is
     ancestral.  First such node is 'descendant'; second
     such node is 'sibling' of first; third such node is
     sibling of second; etc.  */
  node *linker, *ptr;
  short i, j, k;

  linker = poynter;
  for (i = 0; i < (npairs); i++) {
    for (j = 1; j <= 2; j++) {
      if (poynter->state == pair[i][j - 1]) {
	if (j == 1)
	  *otherone = pair[i][1];
	else
	  *otherone = pair[i][0];
	if (*otherone != '.' && *otherone != poynter->ancstr->state) {
	  k = offset + 1;
	  while (*otherone != symbarray[k - 1])
	    k++;
	  if (nodes[k - offset - 1] != NULL)
	    exit(-1);
	  ptr = (node *)Malloc(sizeof(node));
	  ptr->ancstr = poynter;
	  ptr->descendant = NULL;
	  ptr->sibling = NULL;
	  ptr->state = *otherone;
	  if (linker == poynter)   /* If not first */
	    poynter->descendant = ptr;   /* If first */
	  else
	    linker->sibling = ptr;
	  nodes[k - offset - 1] = ptr;
	      /* Save pntr to node */
	  linker = ptr;
	}
      }
    }
  }
}  /* attachnodes */

void maketree(poynter, otherone)
node *poynter;
Char *otherone;
{
  /* Recursively attach nodes */
  if (poynter == NULL )
    return;
  attachnodes(poynter, otherone);
  maketree(poynter->descendant, otherone);
  maketree(poynter->sibling, otherone);
}  /* maketree */

void construct()
{
  /* Puts tree together from array 'pairs' */
  Char rootstate;
  short i, j, k;
  boolean done;
  node *poynter;
  char otherone;

  rooted = false;
  ancsymbol[charindex - 1] = '?';
  rootstate = pair[0][0];
  nstates = 0;
  for (i = 0; i < (npairs); i++) {
    for (j = 1; j <= 2; j++) {
      k = 1;
      done = false;
      while (!done) {
	if (k > nstates) {
	  done = true;
	  break;
	}
	if (pair[i][j - 1] == symbarray[offset + k - 1])
	  done = true;
	else
	  k++;
      }
      if (k > nstates) {
	if (pair[i][j - 1] == '.') {
	  if (rooted)
	    exit(-1);
	  rooted = true;
	  ancsymbol[charindex - 1] = '0';
	  if (j == 1)
	    rootstate = pair[i][1];
	  else
	    rootstate = pair[i][0];
	} else {
	  nstates++;
	  symbarray[offset + nstates - 1] = pair[i][j - 1];
	}
      }
    }
  }
  if ((rooted && nstates != npairs) ||
      (!rooted && nstates != npairs + 1))
    exit(-1);
  root = (node *)Malloc(sizeof(node));
  root->state = ' ';
  root->descendant = (node *)Malloc(sizeof(node));
  root->descendant->ancstr = root;
  root = root->descendant;
  root->descendant = NULL;
  root->sibling = NULL;
  root->state = rootstate;
  for (i = 0; i < (nstates); i++)
    nodes[i] = NULL;
  i = 1;
  while (symbarray[offset + i - 1] != rootstate)
    i++;
  nodes[i - 1] = root;
  maketree(root, &otherone);
  for (i = 0; i < (nstates); i++) {
    if (nodes[i] != root) {
      if (nodes[i] == NULL){
	printf("CHARACTER%4hd:  INVALID CHARACTER STATE TREE DESCRIPTION\n",
	       charnumber);
	exit(-1);}
      else {
	poynter = nodes[i]->ancstr;
	while (poynter != root && poynter != nodes[i])
	  poynter = poynter->ancstr;
	if (poynter != root){
	  printf("CHARACTER%4hd:  INVALID CHARACTER STATE TREE DESCRIPTION\n",
		 charnumber);
	  exit(-1);}
      }
    }
  }
}  /* construct */


void numberedges(poynter,edgenum)
node *poynter;
short *edgenum;
{
  /* Assign to each node a number for the edge below it.
     The root is zero */
  if (poynter == NULL)
    return;
  poynter->edge = *edgenum;
  (*edgenum)++;
  numberedges(poynter->descendant, edgenum);
  numberedges(poynter->sibling, edgenum);
}  /* numberedges */

void factortree()
{
  /* Generate the string of 0's and 1's that will be
     substituted for each symbol of the multistate char. */
  short i, j, place, factoroffset;
  node *poynter;
  short edgenum=0;

  numberedges(root, &edgenum);
  factoroffset = offset + nstates;
  for (i = 0; i < (nstates); i++) {
    place = factoroffset + (nstates - 1) * i;
    for (j = place; j <= (place + nstates - 2); j++)
      symbarray[j] = '0';
    poynter = nodes[i];
    while (poynter != root) {
      symbarray[place + poynter->edge - 1] = '1';
      poynter = poynter->ancstr;
    }
  }
}  /* factortree */


void dotrees()
{
  /* Process character-state trees */
  short lastchar;

  charindex = 0;
  lastchar = 0;
  offset = 0;
  charnumber = 0;
  fscanf(infile, "%hd", &charnumber);
  while (charnumber < 999) {
    if (charnumber < lastchar) {
      printf("CHARACTER%4hd:  OUT OF ORDER", charnumber);
      exit(-1);
    }
    charindex++;
    lastindex = charindex;
    readtree();   /* Process character-state tree  */
    if (npairs > 0) {
      construct();   /* Link tree together  */
      factortree();
    } else {
      nstates = 0;
      ancsymbol[charindex - 1] = '?';
    }
    lastchar = charnumber;
    charnum[charindex - 1] = charnumber;
    chstart[charindex - 1] = offset;
    numstates[charindex - 1] = nstates;
    offset += nstates * nstates;
    fscanf(infile, "%hd", &charnumber);
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);

  /*    each multistate character */
  /*    symbol  */
}  /* dotrees */



void writech(ch, chposition)
Char ch;
short *chposition;
{
  /* Writes a single character to output */
  if (*chposition > maxoutput) {
    putc('\n', outfile);
    *chposition = 1;
  }
  putc(ch, outfile);
  (*chposition)++;
}  /* writech */

Local Void writefactors(chposition)
short *chposition;
{  /* Writes 'FACTORS' line */

  short i, charindex;
  Char symbol;

  fprintf(outfile, "FACTORS   ");
  *chposition = 11;
  symbol = '-';
  for (charindex = 0; charindex < (lastindex); charindex++) {
    if (symbol == '-')
      symbol = '+';
    else
      symbol = '-';
    if (numstates[charindex] == 0)
      writech(symbol,chposition);
    else {
      for (i = 1; i < (numstates[charindex]); i++)
	writech(symbol, chposition);
    }
  }
  putc('\n', outfile);
}  /* writefactors */

void writeancestor(chposition)
short *chposition;
{
  /* Writes 'ANCESTOR' line */
  short i, charindex;

  charindex = 1;
  while (ancsymbol[charindex - 1] == '?')
    charindex++;
  if (charindex > lastindex)
    return;
  fprintf(outfile, "ANCESTOR  ");
  *chposition = 11;
  for (charindex = 0; charindex < (lastindex); charindex++) {
    if (numstates[charindex] == 0)
      writech(ancsymbol[charindex], chposition);
    else {
      for (i = 1; i < (numstates[charindex]); i++)
	writech(ancsymbol[charindex], chposition);
    }
  }
  putc('\n', outfile);
}  /* writeancestor */

void doeu(chposition,eu)
short *chposition,eu;
{
  /* Writes factored data for a single species  */
  short i, charindex, place;
  Char *multichar;

  for (i = 1; i <= nmlngth; i++) {
    ch = getc(infile);
    putc(ch, outfile);
  }
  multichar = (Char *)Malloc(nchars*sizeof(Char));
  *chposition = 11;
  for (i = 0; i < (nchars); i++) {
    do {
      if (eoln(infile)) {
	fscanf(infile, "%*[^\n]");
	getc(infile);
      }
      ch = getc(infile);
    } while (ch == ' ');
    multichar[i] = ch;
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  for (charindex = 0; charindex < (lastindex); charindex++) {
    if (numstates[charindex] == 0)
      writech(multichar[charnum[charindex] - 1], chposition);
    else {
      i = 1;
      while (symbarray[chstart[charindex] + i - 1] !=
	     multichar[charnum[charindex] - 1] && i <= numstates[charindex])
	i++;
      if (i > numstates[charindex]) {
	if( multichar[charnum[charindex] - 1] == unkchar){
	  for (i = 1; i < (numstates[charindex]); i++)
	    writech('?', chposition);
	} else {
	  putc('\n', outfile);
	  printf("IN SPECIES %3hd, MULTISTATE CHARACTER%4hd:  ",
		 eu,charnum[charindex]);
	  printf("'%c' IS NOT A DOCUMENTED STATE\n",
		 multichar[charnum[charindex] - 1]);
	  exit(-1);
	}
      } else {
	place = chstart[charindex] + numstates[charindex] +
		(numstates[charindex] - 1) * (i - 1);
	for (i = 0; i <= (numstates[charindex] - 2); i++)
	  writech(symbarray[place + i], chposition);
      }
    }
  }
  putc('\n', outfile);
  free(multichar);
}  /* doeu */


void dodatamatrix()
{
  /* Reads species information and write factored data set */
  short charindex, totalfactors,eu,chposition;
  totalfactors = 0;
  for (charindex = 0; charindex < (lastindex); charindex++) {
    if (numstates[charindex] == 0)
      totalfactors++;
    else
      totalfactors += numstates[charindex] - 1;
  }
  if (rooted && ancstrrequest)
    fprintf(outfile, "%5hd %5hd    A", neus + 1, totalfactors);
  else
    fprintf(outfile, "%5hd %5hd", neus, totalfactors);
  if (factorrequest)
    fprintf(outfile, " F\n");
  else
    fprintf(outfile, "\n");
  if (factorrequest)
    writefactors(&chposition);
  if (ancstrrequest)
    writeancestor(&chposition);
  eu = 1;
  while (eu <= neus) {
    eu++;
    doeu(&chposition,eu);
  }
  if (progress)
    printf("\nData matrix written on output file\n\n");
}  /* dodatamatrix */


main(argc, argv)
int argc;
Char *argv[];
{
char infilename[100],outfilename[100];
#ifdef MAC
  macsetup("Factor","");
  argv[0] = "Factor";
#endif
  openfile(&infile,INFILE,"r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);

  getoptions();
  fscanf(infile, "%hd%hd%*[^\n]", &neus, &nchars);
  getc(infile);
  charnum = (short *)Malloc(nchars*sizeof(short));
  chstart = (short *)Malloc(nchars*sizeof(short));
  numstates = (short *)Malloc(nchars*sizeof(short));
  ancsymbol = (Char *)Malloc(nchars*sizeof(Char));
  dotrees();   /* Read and factor character-state trees */
  dodatamatrix();
  FClose(infile);
  FClose(outfile);
#ifdef MAC
  fixmacfile(outfilename);
#endif
  exit(0);
}  /* factor */


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

