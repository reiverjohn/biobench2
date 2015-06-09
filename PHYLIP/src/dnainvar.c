#include "phylip.h"

/* version 3.52c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define maxsp           4   /* maximum number of species -- must be 4 */
#define nmlngth         10   /* max. number characters in species name */

#define ibmpc0          false
#define ansi0           true
#define vt520           false


typedef enum {
  xx, yy, zz, ww
} simbol;
typedef Char **sequence;
typedef Char naym[nmlngth];


Static FILE *infile, *outfile;
Static long numsp, sites, endsite, datasets, ith;
Static boolean weights, anerror, printdata, progress, prntpat, printinv,
               mulsets, interleaved, ibmpc, vt52, ansi, firstset;
Static naym nayme[maxsp];
Static sequence y;
Static long *weight,*alias,*aliasweight;

long f[(long)ww - (long)xx + 1][(long)ww - (long)xx + 1]
       [(long)ww - (long)xx + 1]; /* made global from being local to makeinv */

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
{
   *ch = (islower(*ch) ?  toupper(*ch) : (*ch));
}  /* uppercase */


void getnums()
{
  /* input number of species, number of sites */
  fscanf(infile, "%ld%ld", &numsp, &sites);
  if (numsp > maxsp){
    printf("TOO MANY SPECIES: only 4 allowed\n");
    exit(-1);}
  if (printdata)
    fprintf(outfile, "%4ld Species, %4ld Sites\n", numsp, sites);
}  /* getnums */

void getoptions()
{
  /* interactively set options */
  Char ch;
  boolean done, done1;

  fprintf(outfile, "\nNucleic acid sequence Invariants ");
  fprintf(outfile, "method, version %s\n\n",VERSION);
  putchar('\n');
  printdata = false;
  progress = true;
  prntpat = true;
  printinv = true;
  interleaved = true;
  do {
    if (ansi)
      printf("\033[2J\033[H");
    else if (vt52)
      printf("\033E\033H");
    else
      putchar('\n');
    printf("\nNucleic acid sequence Invariants ");
    printf("method, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?");
    if (interleaved)
      printf("  Yes\n");
    else
      printf("  No, sequential\n");
    printf("  0   Terminal type (IBM PC, VT52, ANSI)?");
    if (ibmpc)
      printf("  IBM PC\n");
    if (ansi)
      printf("  ANSI\n");
    if (vt52)
      printf("  VT52\n");
    if (!(ibmpc || vt52 || ansi))
      printf("  (none)\n");
    printf("  1    Print out the data at start of run");
    if (printdata)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  2  Print indications of progress of run  %s\n",
	   (progress ? "Yes" : "No"));
    printf("  3      Print out the counts of patterns");
    if (prntpat)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  4              Print out the invariants");
    if (printinv)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      if (ch == 'M' || ch == 'I' || ch == '1' || ch == '2' || ch == '3' ||
	  ch == '4' || ch == '0') {
	switch (ch) {

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
	  prntpat = !prntpat;
	  break;

	case '4':
	  printinv = !printinv;
	  break;
	}
      } else
	printf("Not a possible option!\n");
    }
  } while (!done);
}  /* getoptions */


void doinit()
{
  /* initializes variables */
  long i;

  getnums();
  if (!anerror)
    getoptions();
  y       = (Char **)Malloc(numsp*sizeof(Char *));
  for (i = 0; i < numsp; i++)
    y[i] = (Char *)Malloc(sites*sizeof(Char));

  weight       = (long *)Malloc(sites * sizeof(long));
  alias        = (long *)Malloc(sites * sizeof(long));
  aliasweight  = (long *)Malloc(sites * sizeof(long));
}  /* doinit*/


void inputweights()
{
  /* input the character weights, which must be 0 or 1 */
  Char ch;
  long i;

  for (i = 1; i < nmlngth; i++)
    ch = getc(infile);
  for (i = 0; i < sites; i++) {
    do {
      if (eoln(infile)) {
	fscanf(infile, "%*[^\n]");
	getc(infile);
      }
      ch = getc(infile);
    } while (ch == ' ');
    weight[i] = 1;
    if (ch == '0')
      weight[i] = 0;
    else if (ch != '1') {
      printf("BAD WEIGHT CHARACTER: %c\n", ch);
      anerror = true;
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  weights = true;
}  /* inputweights */

void printweights()
{
  /* print out the weights of sites */
  long i, j, k;

  fprintf(outfile, "\n\n   Sites are weighted as follows:\n");
  fprintf(outfile, "        ");
  for (i = 0; i <= 9; i++)
    fprintf(outfile, "%3ld", i);
  fprintf(outfile, "\n     *---------------------------------\n");
  for (j = 0; j <= (sites/10); j++) {
    fprintf(outfile, "%5ld!  ", j * 10);
    for (i = 0; i <= 9; i++) {
      k = j * 10 + i;
      if (k > 0 && k <= sites)
	fprintf(outfile, "%3ld", weight[k - 1]);
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
  long extranum, i, cursp, curst;

  if (!firstset) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%ld%ld", &cursp, &curst);
    if (cursp != numsp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4ld\n", ith);
      anerror = true;
    }
    sites = curst;
  }
  weights = false;
  extranum = 0;
  while (!(eoln(infile) || anerror)) {
    ch = getc(infile);
    uppercase(&ch);
    if (ch == 'W')
      extranum++;
    else if (ch != ' ') {
      printf("BAD OPTION CHARACTER: %c\n", ch);
      anerror = true;
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  for (i = 0; i < sites; i++)
    weight[i] = 1;
  for (i = 1; i <= extranum; i++) {
    if (!anerror) {
      ch = getc(infile);
      uppercase(&ch);
      if (ch == 'W')
	inputweights();
      anerror = (ch != 'W');
      if (anerror)
	printf("ERROR: INCORRECT AUXILIARY OPTIONS LINE WHICH STARTS WITH %c\n",
	       ch);
    }
  }
  if (weights)
    printweights();
}  /* inputoptions */

void inputdata()
{
  /* Input the names and sequences for each species */
  long i, j, k, l, basesread, basesnew;
  Char charstate;
  boolean allread, done;

  putc('\n', outfile);
  j = nmlngth + (sites + (sites - 1) / 10) / 2 - 5;
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
  while (!(allread || anerror)) {
    allread = true;
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    i = 1;
    while (i <= numsp && !anerror) {
      if ((interleaved && basesread == 0) || !interleaved) {
	for (j = 0; j < nmlngth; j++) {
	  if (!anerror) {
	    nayme[i - 1][j] = getc(infile);
	    anerror = (anerror || eof(infile) || eoln(infile));
	    if (anerror)
	      printf(
		"ERROR: END-OF-LINE OR END-OF-FILE IN THE MIDDLE OF A SPECIES NAME\n");
	  }
	}
      }
      if (interleaved)
	j = basesread;
      else
	j = 0;
      done = false;
      while (!done && !eof(infile) && !anerror) {
	if (interleaved)
	  done = true;
	while (j < sites && !(eoln(infile) || eof(infile)) && !anerror) {
	  charstate = getc(infile);
	  if (charstate == ' ' || (charstate >= '0' && charstate <= '9'))
	    continue;
	  uppercase(&charstate);
	  if (charstate != 'A' && charstate != 'B' && charstate != 'C' &&
	      charstate != 'D' && charstate != 'G' && charstate != 'H' &&
	      charstate != 'K' && charstate != 'M' && charstate != 'N' &&
	      charstate != 'R' && charstate != 'S' && charstate != 'T' &&
	      charstate != 'U' && charstate != 'V' && charstate != 'W' &&
	      charstate != 'X' && charstate != 'Y' && charstate != '?' &&
	      charstate != 'O' && charstate != '-' && charstate != '.') {

	    printf("ERROR: BAD BASE:%c AT POSITION%5ld OF SPECIES %3ld\n",
		   charstate, j, i);
	    anerror = true;
	  }
	  j++;
	  if (charstate == '.')
	    charstate = y[0][j - 1];
	  y[i - 1][j - 1] = charstate;
	}
	if (interleaved)
	  continue;
	if (j < sites) {
	  fscanf(infile, "%*[^\n]");
	  getc(infile);
	} else if (j == sites)
	  done = true;
      }
      if (interleaved && i == 1)
	basesnew = j;
      fscanf(infile, "%*[^\n]");
      getc(infile);
      if (interleaved)
	anerror = (anerror || j != basesnew);
      else
	anerror = (anerror || j != sites);
      if (anerror)
	printf("ERROR: SEQUENCES OUT OF ALIGNMENT\n");
      i++;
    }
    if (interleaved) {
      basesread = basesnew;
      allread = (basesread == sites);
    } else
      allread = (i > numsp);
  }
  if (!printdata || anerror)
    return;
  for (i = 1; i <= ((sites - 1) / 60 + 1); i++) {
    for (j = 1; j <= numsp; j++) {
      for (k = 0; k < nmlngth; k++)
	putc(nayme[j - 1][k], outfile);
      fprintf(outfile, "   ");
      l = i * 60;
      if (l > sites)
	l = sites;
      for (k = (i - 1) * 60 + 1; k <= l; k++) {
	if (!anerror) {
	  if (j > 1 && y[j - 1][k - 1] == y[0][k - 1])
	    charstate = '.';
	  else
	    charstate = y[j - 1][k - 1];
	  putc(charstate, outfile);
	  if (k % 10 == 0 && k % 60 != 0)
	    putc(' ', outfile);
	}
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
  long gap, i, j, jj, jg, k, itemp;
  boolean flip, tied;

  gap = sites / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= sites; i++) {
      j = i - gap;
      flip = true;
      while (j > 0 && flip) {
	jj = alias[j - 1];
	jg = alias[j + gap - 1];
	flip = false;
	k = 1;
	tied = true;
	while (k <= numsp && tied) {
	  flip = (y[k - 1][jj - 1] > y[k - 1][jg - 1]);
	  tied = (tied && y[k - 1][jj - 1] == y[k - 1][jg - 1]);
	  k++;
	}
	if (!flip)
	  break;
	itemp = alias[j - 1];
	alias[j - 1] = alias[j + gap - 1];
	alias[j + gap - 1] = itemp;
	itemp = aliasweight[j - 1];
	aliasweight[j - 1] = aliasweight[j + gap - 1];
	aliasweight[j + gap - 1] = itemp;
	j -= gap;
      }
    }
    gap /= 2;
  }
}  /* sitesort */

void sitecombine()
{
  /* combine sites that have identical patterns */
  long i, j, k;
  boolean tied;

  i = 1;
  while (i < sites) {
    j = i + 1;
    tied = true;
    while (j <= sites && tied) {
      k = 1;
      while (k <= numsp && tied) {
	tied = (tied &&
	    y[k - 1][alias[i - 1] - 1] == y[k - 1][alias[j - 1] - 1]);
	k++;
      }
      if (tied && aliasweight[j - 1] > 0) {
	aliasweight[i - 1] += aliasweight[j - 1];
	aliasweight[j - 1] = 0;
      }
      j++;
    }
    i = j - 1;
  }
}  /* sitecombine */

void sitescrunch()
{
  /* Bubble sort so positively weighted sites come first */
  long i, j, itemp;
  boolean done, found;

  done = false;
  i = 1;
  j = 2;
  while (!done) {
    found = false;
    if (aliasweight[i - 1] > 0)
      i++;
    else {
      if (j <= i)
	j = i + 1;
      if (j <= sites) {
	found = false;
	do {
	  found = (aliasweight[j - 1] > 0);
	  j++;
	} while (!(found || j > sites));
	if (found) {
	  j--;
	  itemp = alias[i - 1];
	  alias[i - 1] = alias[j - 1];
	  alias[j - 1] = itemp;
	  itemp = aliasweight[i - 1];
	  aliasweight[i - 1] = aliasweight[j - 1];
	  aliasweight[j - 1] = itemp;
	} else
	  done = true;
      }
    }
    done = (done || i >= sites);
  }
}  /* sitescrunch */

void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  long i;

  for (i = 1; i <= sites; i++) {
    alias[i - 1] = i;
    aliasweight[i - 1] = weight[i - 1];
  }
  sitesort();
  sitecombine();
  sitescrunch();
  for (i = 1; i <= sites; i++) {
    weight[i - 1] = aliasweight[i - 1];
    if (weight[i - 1] > 0)
      endsite = i;
  }
}  /* makeweights */


void doinput()
{  /* getinput */
  /* reads the input data */
  if (!anerror)
    inputoptions();
  if (!anerror)
    inputdata();
  if (!anerror)
    makeweights();
}  /* getinput */



void prntpatterns()
{
  /* print out patterns */
  long i, j;

  fprintf(outfile, "\n   Pattern");
  if (prntpat)
    fprintf(outfile, "   Number of times");
  fprintf(outfile, "\n\n");
  for (i = 0; i < endsite; i++) {
    fprintf(outfile, "     ");
    for (j = 0; j < numsp; j++)
      putc(y[j][alias[i] - 1], outfile);
    if (prntpat)
      fprintf(outfile, "  %8ld", weight[i]);
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* prntpatterns */

void makesymmetries()
{
  /* get frequencies of symmetrized patterns */
  long i, j;
  boolean drop, usedz;
  Char ch, ch1, zchar;
  simbol s1, s2, s3;
  simbol t[maxsp - 1];

  for (s1 = xx; (long)s1 <= (long)ww; s1 = (simbol)((long)s1 + 1)) {
    for (s2 = xx; (long)s2 <= (long)ww; s2 = (simbol)((long)s2 + 1)) {
      for (s3 = xx; (long)s3 <= (long)ww; s3 = (simbol)((long)s3 + 1))
	f[(long)s1 - (long)xx][(long)s2 - (long)xx]
	  [(long)s3 - (long)xx] = 0;
    }
  }
  for (i = 0; i < endsite; i++) {
    drop = false;
    for (j = 0; j < numsp; j++) {
      ch = y[j][alias[i] - 1];
      drop = (drop ||
	      (ch != 'A' && ch != 'C' && ch != 'G' && ch != 'T' && ch != 'U'));
    }
    ch1 = y[0][alias[i] - 1];
    if (!drop) {
      usedz = false;
      zchar = ' ';
      for (j = 2; j <= numsp; j++) {
	ch = y[j - 1][alias[i] - 1];
	if (ch == ch1)
	  t[j - 2] = xx;
	else if ((ch1 == 'A' && ch == 'G') || (ch1 == 'G' && ch == 'A') ||
		 (ch1 == 'C' && (ch == 'T' || ch == 'U')) ||
		 ((ch1 == 'T' || ch1 == 'U') && ch == 'C'))
	  t[j - 2] = yy;
	else if (!usedz) {
	  t[j - 2] = zz;
	  usedz = true;
	  zchar = ch;
	} else if (usedz && ch == zchar)
	  t[j - 2] = zz;
	else if (usedz && ch != zchar)
	  t[j - 2] = ww;
      }
      f[(long)t[0] - (long)xx][(long)t[1] - (long)xx]
	[(long)t[2] - (long)xx] += weight[i];
    }
  }
}  /* makesymmetries */

void prntsymbol(s)
simbol s;
{
  /* print 1, 2, 3, 4 as appropriate */
  switch (s) {

  case xx:
    putc('1', outfile);
    break;

  case yy:
    putc('2', outfile);
    break;

  case zz:
    putc('3', outfile);
    break;

  case ww:
    putc('4', outfile);
    break;
  }
}  /* prntsymbol */

void prntsymmetries()
{
  /* print out symmetrized pattern numbers */
  simbol s1, s2, s3;

  fprintf(outfile, "\nSymmetrized patterns (1, 2 = the two purines  ");
  fprintf(outfile, "and  3, 4 = the two pyrimidines\n");
  fprintf(outfile, "                  or  1, 2 = the two pyrimidines  ");
  fprintf(outfile, "and  3, 4 = the two purines)\n\n");
  for (s1 = xx; (long)s1 <= (long)ww; s1 = (simbol)((long)s1 + 1)) {
    for (s2 = xx; (long)s2 <= (long)ww; s2 = (simbol)((long)s2 + 1)) {
      for (s3 = xx; (long)s3 <= (long)ww; s3 = (simbol)((long)s3 + 1)) {
	if (f[(long)s1 - (long)xx][(long)s2 - (long)xx]
	    [(long)s3 - (long)xx] > 0) {
	  fprintf(outfile, "     1");
	  prntsymbol(s1);
	  prntsymbol(s2);
	  prntsymbol(s3);
	  if (prntpat)
	    fprintf(outfile, "   %7ld",
		    f[(long)s1 - (long)xx][(long)s2 - (long)xx]
		    [(long)s3 - (long)xx]);
	  putc('\n', outfile);
	}
      }
    }
  }
}  /* prntsymmetries */


void tabulate(mm, nn, pp, qq, mr,nr,pr,qr)
long mm, nn, pp, qq;
double *mr,*nr,*pr,*qr;
{
  /* make quadratic invariant, table, chi-square */
  long total;
  double k, TEMP;

  fprintf(outfile, "\n   Contingency Table\n\n");
  fprintf(outfile, "%7ld%6ld\n", mm, nn);
  fprintf(outfile, "%7ld%6ld\n\n", pp, qq);
  *mr = (long)(mm);
  *nr = (long)(nn);
  *pr = (long)pp;
  *qr = (long)qq;
  total = mm + nn + pp + qq;
  if (printinv)
    fprintf(outfile, "   Quadratic invariant = %15.1f\n\n",
	    (*nr) * (*pr) - (*mr) * (*qr));
  fprintf(outfile, "   Chi-square = ");
  TEMP = (*mr) * (*qr) - (*nr) * (*pr);
  k = total * (TEMP * TEMP) / (((*mr) + (*nr)) * ((*mr) + (*pr)) *
			       ((*nr) + (*qr)) * ((*pr) + (*qr)));
  fprintf(outfile, "%10.5f", k);
  if ((*mr) * (*qr) > (*nr) * (*pr) && k > 2.71)
    fprintf(outfile, " (P < 0.05)\n");
  else
    fprintf(outfile, " (not significant)\n");
  fprintf(outfile, "\n\n");
}  /* tabulate */


void writename(m)
long m;
{
  /* write out a species name */
  long i, n;

  n = nmlngth;
  while (nayme[m - 1][n - 1] == ' ')
    n--;
  if (n == 0)
    n = 1;
  for (i = 0; i < n; i++)
    putc(nayme[m - 1][i], outfile);
}  /* writename */

void writetree(i, j, k, l)
long i, j, k, l;
{
  /* write out tree topology ((i,j),(k,l)) using names */
  fprintf(outfile, "((");
  writename(i);
  putc(',', outfile);
  writename(j);
  fprintf(outfile, "),(");
  writename(k);
  putc(',', outfile);
  writename(l);
  fprintf(outfile, "))\n");
}  /* writetree */

void exacttest(m, n)
long m, n;
{
  /* exact binomial test that m <= n */
  long i;
  double p, sum;

  p = 1.0;
  for (i = 1; i <= m + n; i++)
    p /= 2.0;
  sum = p;
  for (i = 1; i <= n; i++) {
    p = p * (m + n - i + 1) / i;
    sum += p;
  }
  fprintf(outfile, "      %7.4f", sum);
  if (sum <= 0.05)
    fprintf(outfile, "              yes\n");
  else
    fprintf(outfile, "               no\n");
}  /* exacttest */

void invariants()
{
  /* compute invariants */
  long  m, n, p, q;
  double L1, L2, L3;
  double mr,nr,pr,qr;

  fprintf(outfile, "\nTree topologies (unrooted): \n\n");
  fprintf(outfile, "    I:  ");
  writetree(1, 2, 3, 4);
  fprintf(outfile, "   II:  ");
  writetree(1, 3, 2, 4);
  fprintf(outfile, "  III:  ");
  writetree(1, 4, 2, 3);
  fprintf(outfile, "\n\nLake's linear invariants\n");
  fprintf(outfile,
    " (these are expected to be zero for the two incorrect tree topologies.\n");
  fprintf(outfile,
	  "  This is tested by testing the equality of the two parts\n");
  fprintf(outfile,
	  "  of each expression using a one-sided exact binomial test.\n");
  fprintf(outfile,
    "  The null hypothesis is that the first part is no larger than the second.)\n\n");
  fprintf(outfile, " Tree                           ");
  fprintf(outfile, "  Exact test P value    Significant?\n\n");
  m = f[(long)yy - (long)xx][(long)zz - (long)xx]
      [(long)ww - (long)xx] + f[0][(long)zz - (long)xx]
      [(long)zz - (long)xx];
  n = f[(long)yy - (long)xx][(long)zz - (long)xx]
      [(long)zz - (long)xx] + f[0][(long)zz - (long)xx]
      [(long)ww - (long)xx];
  fprintf(outfile, "   I  %5ld    - %5ld   = %5ld", m, n, m - n);
  exacttest(m, n);
  m = f[(long)zz - (long)xx][(long)yy - (long)xx]
      [(long)ww - (long)xx] + f[(long)zz - (long)xx][0]
      [(long)zz - (long)xx];
  n = f[(long)zz - (long)xx][(long)yy - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx][0]
      [(long)ww - (long)xx];
  fprintf(outfile, "   II %5ld    - %5ld   = %5ld", m, n, m - n);
  exacttest(m, n);
  m = f[(long)zz - (long)xx][(long)ww - (long)xx]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx]
      [(long)zz - (long)xx][0];
  n = f[(long)zz - (long)xx][(long)zz - (long)xx]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx]
      [(long)ww - (long)xx][0];
  fprintf(outfile, "   III%5ld    - %5ld   = %5ld", m, n, m - n);
  exacttest(m, n);
  fprintf(outfile, "\n\nCavender's quadratic invariants (type L)");
  fprintf(outfile, " using purines vs. pyrimidines\n");
  fprintf(outfile,
	  " (these are expected to be zero, and thus have a nonsignificant\n");
  fprintf(outfile, "  chi-square, for the correct tree topology)\n");
  fprintf(outfile, "They will be misled if there are substantially\n");
  fprintf(outfile, "different evolutionary rate between sites, or\n");
  fprintf(outfile, "different purine:pyrimidine ratios from 1:1.\n\n");
  fprintf(outfile, "  Tree I:\n");
  m = f[0][0][0] + f[0][(long)yy - (long)xx]
      [(long)yy - (long)xx] + f[0][(long)zz - (long)xx]
      [(long)zz - (long)xx];
  n = f[0][0][(long)yy - (long)xx] + f[0][0]
      [(long)zz - (long)xx] + f[0][(long)yy - (long)xx][0] + f[0]
      [(long)yy - (long)xx][(long)zz - (long)xx] + f[0]
      [(long)zz - (long)xx][0] + f[0][(long)zz - (long)xx]
      [(long)yy - (long)xx] + f[0][(long)zz - (long)xx]
      [(long)ww - (long)xx];
  p = f[(long)yy - (long)xx][0][0] + f[(long)yy - (long)xx]
      [(long)yy - (long)xx]
      [(long)yy - (long)xx] + f[(long)yy - (long)xx]
      [(long)zz - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx][0]
      [0] + f[(long)zz - (long)xx][(long)yy - (long)xx]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx]
      [(long)zz - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx]
      [(long)ww - (long)xx][(long)ww - (long)xx];
  q = f[(long)yy - (long)xx][0][(long)yy - (long)xx] +
      f[(long)yy - (long)xx][0][(long)zz - (long)xx] +
      f[(long)yy - (long)xx][(long)yy - (long)xx][0] +
      f[(long)yy - (long)xx][(long)yy - (long)xx][(long)zz - (long)xx] +
      f[(long)yy - (long)xx][(long)zz - (long)xx][0] +
      f[(long)yy - (long)xx][(long)zz - (long)xx][(long)yy - (long)xx] +
      f[(long)yy - (long)xx][(long)zz - (long)xx][(long)ww - (long)xx] +
      f[(long)zz - (long)xx][0][(long)yy - (long)xx] +
      f[(long)zz - (long)xx][0][(long)zz - (long)xx] +
      f[(long)zz - (long)xx][0][(long)ww - (long)xx] +
      f[(long)zz - (long)xx][(long)yy - (long)xx][0] +
      f[(long)zz - (long)xx][(long)yy - (long)xx][(long)zz - (long)xx] +
      f[(long)zz - (long)xx][(long)yy - (long)xx][(long)ww - (long)xx] +
      f[(long)zz - (long)xx][(long)zz - (long)xx][0] +
      f[(long)zz - (long)xx][(long)zz - (long)xx][(long)yy - (long)xx] +
      f[(long)zz - (long)xx][(long)zz - (long)xx][(long)ww - (long)xx] +
      f[(long)zz - (long)xx][(long)ww - (long)xx][0] +
      f[(long)zz - (long)xx][(long)ww - (long)xx][(long)yy - (long)xx] +
      f[(long)zz - (long)xx][(long)ww - (long)xx][(long)zz - (long)xx];

  nr = n;
  pr = p;
  mr = m;
  qr = q;
  L1 = nr * pr - mr * qr;
  tabulate(m, n, p, q, &mr,&nr,&pr,&qr);
  fprintf(outfile, "  Tree II:\n");
  m = f[0][0][0] + f[(long)yy - (long)xx][0]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx][0]
      [(long)zz - (long)xx];
  n = f[0][0][(long)yy - (long)xx] + f[0][0]
      [(long)zz - (long)xx] + f[(long)yy - (long)xx][0]
      [0] + f[(long)yy - (long)xx][0]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx][0]
      [0] + f[(long)zz - (long)xx][0]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx][0]
      [(long)ww - (long)xx];
  p = f[0][(long)yy - (long)xx][0] + f[(long)yy - (long)xx]
      [(long)yy - (long)xx]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx]
      [(long)yy - (long)xx][(long)zz - (long)xx] + f[0]
      [(long)zz - (long)xx][0] + f[(long)yy - (long)xx]
      [(long)zz - (long)xx]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx]
      [(long)zz - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx]
      [(long)ww - (long)xx][(long)zz - (long)xx];
  q = f[0][(long)yy - (long)xx][(long)yy - (long)xx] + f[0]
      [(long)yy - (long)xx][(long)zz - (long)xx] +
      f[(long)yy - (long)xx][(long)yy - (long)xx][0] +
      f[(long)yy - (long)xx][(long)yy - (long)xx][(long)zz - (long)xx] +
      f[(long)zz - (long)xx][(long)yy - (long)xx][0] +
      f[(long)zz - (long)xx][(long)yy - (long)xx][(long)yy - (long)xx] +
      f[(long)zz - (long)xx][(long)yy - (long)xx][(long)ww - (long)xx] +
      f[0][(long)zz - (long)xx][(long)yy - (long)xx] + f[0]
      [(long)zz - (long)xx][(long)zz - (long)xx] + f[0]
      [(long)zz - (long)xx][(long)ww - (long)xx] +
      f[(long)yy - (long)xx][(long)zz - (long)xx][0] +
      f[(long)yy - (long)xx][(long)zz - (long)xx][(long)zz - (long)xx] +
      f[(long)yy - (long)xx][(long)zz - (long)xx][(long)ww - (long)xx] +
      f[(long)zz - (long)xx][(long)zz - (long)xx][0] +
      f[(long)zz - (long)xx][(long)zz - (long)xx][(long)yy - (long)xx] +
      f[(long)zz - (long)xx][(long)zz - (long)xx][(long)ww - (long)xx] +
      f[(long)zz - (long)xx][(long)ww - (long)xx][0] +
      f[(long)zz - (long)xx][(long)ww - (long)xx][(long)yy - (long)xx] +
      f[(long)zz - (long)xx][(long)ww - (long)xx][(long)ww - (long)xx];
  nr = n;
  pr = p;
  mr = m;
  qr = q;
  L2 = nr * pr - mr * qr;
  tabulate(m, n, p, q, &mr,&nr,&pr,&qr);
  fprintf(outfile, "  Tree III:\n");
  m = f[0][0][0] + f[(long)yy - (long)xx][(long)yy - (long)xx]
      [0] + f[(long)zz - (long)xx][(long)zz - (long)xx][0];
  n = f[(long)yy - (long)xx][0][0] + f[(long)zz - (long)xx][0]
      [0] + f[0][(long)yy - (long)xx][0] + f[(long)zz - (long)xx]
      [(long)yy - (long)xx][0] + f[0][(long)zz - (long)xx]
      [0] + f[(long)yy - (long)xx][(long)zz - (long)xx]
      [0] + f[(long)zz - (long)xx][(long)ww - (long)xx][0];
  p = f[0][0][(long)yy - (long)xx] + f[(long)yy - (long)xx]
      [(long)yy - (long)xx]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx]
      [(long)zz - (long)xx][(long)yy - (long)xx] + f[0][0]
      [(long)zz - (long)xx] + f[(long)yy - (long)xx]
      [(long)yy - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx]
      [(long)zz - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx]
      [(long)zz - (long)xx][(long)ww - (long)xx];
  q = f[(long)yy - (long)xx][0][(long)yy - (long)xx] + f[(long)zz - (long)xx]
      [0][(long)yy - (long)xx] + f[0][(long)yy - (long)xx][(long)yy - (long)xx] +
      f[(long)zz - (long)xx][(long)yy - (long)xx][(long)yy - (long)xx] +
      f[0][(long)zz - (long)xx]
      [(long)yy - (long)xx] + f[(long)yy - (long)xx][(long)zz - (long)xx]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx][(long)ww - (long)xx]
      [(long)yy - (long)xx] + f[(long)yy - (long)xx][0]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx][0]
      [(long)zz - (long)xx] + f[0][(long)zz - (long)xx]
      [(long)ww - (long)xx] + f[0][(long)yy - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx]
      [(long)yy - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx]
      [(long)yy - (long)xx][(long)ww - (long)xx] + f[0]
      [(long)zz - (long)xx]
      [(long)zz - (long)xx] + f[(long)yy - (long)xx]
      [(long)zz - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx]
      [(long)ww - (long)xx]
      [(long)ww - (long)xx] + f[(long)zz - (long)xx][0]
      [(long)ww - (long)xx] + f[(long)yy - (long)xx]
      [(long)zz - (long)xx][(long)ww - (long)xx] +
      f[(long)zz - (long)xx][(long)ww - (long)xx][(long)zz - (long)xx];
  nr = n;
  pr = p;
  mr = m;
  qr = q;
  L3 = nr * pr - mr * qr;
  tabulate(m, n, p, q, &mr,&nr,&pr,&qr);
  fprintf(outfile, "\n\nCavender's quadratic invariants (type K)");
  fprintf(outfile, " using purines vs. pyrimidines\n");
  fprintf(outfile,
	  " (these are expected to be zero for the correct tree topology)\n");
  fprintf(outfile, "They will be misled if there are substantially\n");
  fprintf(outfile, "different evolutionary rate between sites, or\n");
  fprintf(outfile, "different purine:pyrimidine ratios from 1:1.\n");
  fprintf(outfile, "No statistical test is done on them here.\n\n");
  fprintf(outfile, "  Tree I:   %15.1f\n", L2 - L3);
  fprintf(outfile, "  Tree II:  %15.1f\n", L3 - L1);
  fprintf(outfile, "  Tree III: %15.1f\n\n", L1 - L2);
}  /* invariants */

void makeinv()
{
  /* print out patterns and compute invariants */

  prntpatterns();
  makesymmetries();
  prntsymmetries();
  if (printinv)
    invariants();
}  /* makeinv */


main(argc, argv)
int argc;
Char *argv[];
{  /* DNA Invariants */
char infilename[100],outfilename[100];
#ifdef MAC
  macsetup("Dnainvar","");
  argv[0] = "Dnainvar";
#endif
openfile(&infile,INFILE,"r",argv[0],infilename);
openfile(&outfile,OUTFILE,"w",argv[0],outfilename);

  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  mulsets = false;
  firstset = true;
  datasets = 1;
  anerror = false;
  doinit();
  if (!anerror) {
    for (ith = 1; ith <= datasets; ith++) {
      if (!anerror)
	doinput();
      if (ith == 1)
	firstset = false;
      if (!anerror) {
        if (datasets > 1) {
          if (progress)
            printf("\nData set # %ld:\n",ith);
          fprintf(outfile, "Data set # %ld:\n\n",ith);
        }
	makeinv();
      }
      if (progress) {
        putchar('\n');
        if (!anerror)
	  printf("Output written to output file\n");
        putchar('\n');
      }
    }
  }
  FClose(outfile);
  FClose(infile);
#ifdef MAC
  fixmacfile(outfilename);
#endif
  exit(0);
}  /* DNA Invariants */

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

