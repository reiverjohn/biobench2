
#include "phylip.h"

/* version 3.57c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define nmlngth         10   /* number of characters in species name   */

#define ibmpc0          false
#define ansi0           true
#define vt520           false

typedef Char naym[nmlngth];
typedef long longer[6];
typedef long *stepp;
typedef struct steptr {
  stepp steps;
  double *dummy;
} steptr;
typedef enum {
  seqs, morphology, restsites, genefreqs
} datatype;


Static FILE *infile, *outfile;
Static long spp, sites, loci, maxalleles, groups, newsites, newersites,
             newgroups, newergroups, nenzymes, reps, ws, maxfactor;
Static boolean permute, jackknife, weights, factors, enzymes, all,
	       printdata, progress, interleaved, ibmpc, vt52, ansi;
Static longer seed;
Static datatype data;
Static steptr oldweight, weight, where, how_many, newwhere, newhowmany,
             newerwhere, newerhowmany, factor, newerfactor;
Static naym *nayme;   /* names of species */
Static long *alleles;
Static Char **nodep;
Static double **nodef;
Static long **sppord;



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

void uppercase(ch)
Char *ch;
{  /* convert ch to upper case -- either ASCII or EBCDIC */
   *ch = (islower(*ch) ?  toupper(*ch) : (*ch));
}  /* uppercase */

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



void getoptions()
{
  /* interactively set options */
  long i, inseed, reps0;
  Char ch;
  boolean done1;

  data = seqs;
  jackknife = false;
  all = false;
  reps = 100;
  printdata = false;
  progress = true;
  interleaved = true;
  do {
    printf("\nRandom number seed (must be odd)?\n");
    scanf("%ld%*[^\n]", &inseed);
    getchar();
  } while (!(inseed > 0 && (inseed & 1)));
  for (i = 0; i <= 5; i++)
    seed[i] = 0;
  i = 0;
  do {
    seed[i] = inseed & 63;
    inseed /= 64;
    i++;
  } while (inseed != 0);
  for (;;) {
    printf(ansi ? "\033[2J\033[H" :
	   vt52 ? "\033E\033H"    : "\n");
    printf("\nBootstrapped sequences algorithm, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  D   Sequence, Morph, Rest., Gene Freqs?  %s\n",
	   (data == seqs       ) ? "Molecular sequences"      :
	   (data == morphology ) ? "Discrete Morphology"      :
	   (data == restsites)   ? "Restriction Sites"        :
	   (data == genefreqs)   ? "Gene Frequencies" : "");
    if (data == restsites)
      printf("  E                    Number of enzymes?  %s\n",
	     enzymes ? "Present in input file" :
	               "Not present in input file");
    if (data == genefreqs)
      printf("  A    All alleles present at each locus?  %s\n",
	     all ? "Yes" : "No, one absent at each locus");

    printf("  J     Bootstrap, Jackknife, or Permute?  %s\n",
	   jackknife ? "Delete-half jackknife"                 :
	   permute   ? "Permute species for each character"    :
	               "Bootstrap");
    printf("  R                  How many replicates?%5ld\n", reps);
    if (data == seqs || data == restsites) {
      printf("  I          Input sequences interleaved?  %s\n",
	     interleaved ? "Yes" : "No, sequential");
    }
    printf("  0   Terminal type (IBM PC, VT52, ANSI)?  %s\n",
	   ibmpc ? "IBM PC" :
	   ansi  ? "ANSI"   :
	   vt52  ? "VT52"   : "(none)");
    printf("  1    Print out the data at start of run  %s\n",
	   printdata ? "Yes" : "No");
    printf("  2  Print indications of progress of run  %s\n",
	   progress ? "Yes" : "No");
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    if (ch == 'Y')
	   break;
    if (strchr("ADEJRI120",ch)){
      switch (ch) {
	
      case 'D':
	if (data == genefreqs)
	  data = seqs;
	else
	  data = (datatype)((long)data + 1);
	break;
	
      case 'A':
	all = !all;
	break;
	
      case 'E':
	enzymes = !enzymes;
	break;
	
      case 'J':
	if (permute)
	  permute = false;
	else if (jackknife) {
	  jackknife = false;
	  permute = true;
	} else
	  jackknife = true;
	break;
	
      case 'R':
	done1 = true;
	reps0 = reps;
	do {
	  printf("Number of replicates?\n");
	  scanf("%ld%*[^\n]", &reps);
	  getchar();
	  done1 = (reps > 0);
	  if (!done1) {
	    printf("BAD NUMBER: must be positive\n");
	    reps = reps0;
	  }
	} while (done1 != true);
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
      }
    } else
      printf("Not a possible option!\n");
  }
  if (printdata)
    fprintf(outfile, "\nBootstrapped sequences algorithm, version %s\n\n\n",
	    VERSION);
}  /* getoptions */

Local Void getnums()
{
  /* read numbers of species and of sites */
  long i;

  fscanf(infile, "%ld%ld", &spp, &sites);
  loci = sites;
  maxalleles = 1;
  if (data == restsites && enzymes)
    fscanf(infile, "%ld", &nenzymes);
  if (data == genefreqs) {
    alleles = (long *)Malloc(sites*sizeof(long));
    fscanf(infile, "%*[^\n]");
    getc(infile);
    sites = 0;
    for (i = 0; i < (loci); i++) {
      if (eoln(infile)) {
	fscanf(infile, "%*[^\n]");
	getc(infile);
      }
      fscanf(infile, "%ld", &alleles[i]);
      if (alleles[i] > maxalleles)
         maxalleles = alleles[i];
      if (all)
         sites += alleles[i];
      else
         sites += alleles[i] - 1;
    }
    if (!all)
       maxalleles--;
  }
  if (printdata) {
    if (data == genefreqs)
       fprintf(outfile, "%3ld species, %3ld  loci\n", spp, loci);
    else
       fprintf(outfile, "%3ld species, %3ld  sites\n", spp, sites);
  }
}  /* getnums */

Local Void inputfactors()
{
  long i, j, factor_run=0;
  Char ch, prevch;
  

  for (i = 2; i <= nmlngth; i++)
    ch = getc(infile);
  prevch = ' ';
  j = 0;
  for (i = 0; i < (sites); i++) {
    
    do {    /* Get rid of any white space, */
      if (eoln(infile)) {
	fscanf(infile, "%*[^\n]");
	getc(infile);
      }
      ch = getc(infile);
    } while (ch == ' ');
    /* Then check to see if the factor changed, */
    if (ch != prevch) {
      j++;
      factor_run=1 ;
    } else {
      factor_run++ ;
      if (factor_run > maxfactor)
	maxfactor = factor_run ;
    }

    prevch = ch;
    factor.steps[i] = j;
  }

  fscanf(infile, "%*[^\n]");
  getc(infile);
}  /* inputfactors */

Local Void printfactors()
{
  long i, j;

  fprintf(outfile, "Factors (least significant digit)\n\n");
  for (i = 1; i <= nmlngth + 3; i++)
    putc(' ', outfile);
  for (i = 1; i <= (sites); i++) {
    if (i % 55 == 1 && i != 1) {
      fprintf(outfile, " \n");
      for (j = 1; j <= nmlngth + 3; j++)
      putc(' ', outfile);
    }
    putc(factor.steps[i - 1], outfile);
  }
  fprintf(outfile, "\n\n");
}  /* printfactors */

Local Void inputweights()
{
  /* input the character weights, 0 or 1 */
  Char ch;
  long i;

  for (i = 1; i < nmlngth; i++)
    ch = getc(infile);
  for (i = 0; i < (sites); i++) {
    do {
      if (eoln(infile)) {
	fscanf(infile, "%*[^\n]");
	getc(infile);
      }
      ch = getc(infile);
    } while (ch == ' ');
    oldweight.steps[i] = 1;
    if (ch == '0' || ch == '1')
      oldweight.steps[i] = ch - '0';
    else {
      printf("BAD WEIGHT CHARACTER: %c -- WEIGHTS IN DNABOOT MUST BE 0 OR 1\n",
	     ch);
      exit(-1);
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
}  /* inputweights */

Local Void printweights()
{
  /* print out the weights of sites */
  long i, j;

  fprintf(outfile, "\n       Sites are weighted as follows:\n");
  for (i = 1; i <= (sites); i++) {
    if ((i - 1) % 60 == 0) {
      putc('\n', outfile);
      for (j = 1; j <= nmlngth + 3; j++)
	putc(' ', outfile);
    }
    fprintf(outfile, "%ld", oldweight.steps[i - 1]);
    if (i % 10 == 0 && i % 60 != 0)
      putc(' ', outfile);

  }
  fprintf(outfile, "\n\n");
}  /* printweights */

Local Void inputoptions()
{
  /* input the information on the options */
  Char ch;
  long extranum, i, j, k, l, m;

  factors = false;
  weights = false;
  maxfactor = maxalleles;
  if (data == genefreqs) {
    k = 0;
    l = 0;
    for (i = 0; i < (loci); i++) {
      if (all)
	m = alleles[i];
      else
	m = alleles[i] - 1;
      k++;
      for (j = 1; j <= m; j++) {
	l++;
	factor.steps[l - 1] = k;
      }
    }
  } else {
    for (i = 1; i <= (sites); i++)
      factor.steps[i - 1] = i;
  }
  for (i = 0; i < (sites); i++)
    oldweight.steps[i] = 1;
  extranum = 0;
  while (!eoln(infile)) {
    ch = getc(infile);
    uppercase(&ch);
    if (ch != 'W' && ch != 'F') {
      if (ch != ' ') {
	printf("BAD OPTION CHARACTER: %c\n", ch);
	exit(-1);
      }
      continue;
    }
    switch (ch) {

    case 'F':
      factors = true;
      extranum++;
      break;

    case 'W':
      weights = true;
      extranum++;
      break;
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  for (i = 1; i <= extranum; i++) {
    ch = getc(infile);
    uppercase(&ch);
    if (ch == 'F')
      inputfactors();
    if (ch == 'W')
      inputweights();
    if (ch != 'W' && ch != 'F'){
      printf("ERROR: INCORRECT AUXILIARY OPTIONS LINE WHICH STARTS WITH %c\n",
	       ch);
      exit(-1);}
  }
  if (factors && printdata)
    printfactors();
  if (weights && printdata)
    printweights();
  for (i = 0; i < (loci); i++)
    how_many.steps[i] = 0;
  for (i = 0; i < (loci); i++)
    where.steps[i] = 0;
  for (i = 1; i <= (sites); i++) {
    how_many.steps[factor.steps[i - 1] - 1]++;
    if (where.steps[factor.steps[i - 1] - 1] == 0)
      where.steps[factor.steps[i - 1] - 1] = i;
  }
  groups = factor.steps[sites - 1];
  newgroups = 0;
  newsites = 0;
  for (i = 0; i < (groups); i++) {
    if (oldweight.steps[where.steps[i] - 1] > 0) {
      newgroups++;
      newsites += how_many.steps[i];
      newwhere.steps[newgroups - 1] = where.steps[i];
      newhowmany.steps[newgroups - 1] = how_many.steps[i];
    }
  }
}  /* inputoptions */


Local Void inputdata()
{
  /* input the names and sequences for each species */
  long i, j, k, l, m, n, basesread, basesnew;
  double x;
  Char charstate;
  boolean allread, done;

  if (data == genefreqs) {
    nodef = (double **)Malloc(spp*sizeof(double *));
    for (i = 0; i < (spp); i++)
      nodef[i] = (double *)Malloc(sites*sizeof(double));
  } else {
    nodep = (Char **)Malloc(spp*sizeof(Char *));
    for (i = 0; i < (spp); i++)
      nodep[i] = (Char *)Malloc(sites*sizeof(Char));
  }
  j = nmlngth + (sites + (sites - 1) / 10) / 2 - 5;
  if (j < nmlngth - 1)
    j = nmlngth - 1;
  if (j > 37)
    j = 37;
  if (printdata) {
    fprintf(outfile, "Name");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "Data\n");
    fprintf(outfile, "----");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "----\n\n");
  }
  interleaved = (interleaved && ((data == seqs) || (data == restsites)));
  if (data == genefreqs) {
    for (i = 1; i <= (spp); i++) {
      for (j = 0; j < nmlngth; j++) {
	if (eof(infile) || eoln(infile)){
	  printf("ERROR: END-OF-LINE OR END-OF-FILE");
	  printf(" IN THE MIDDLE OF A SPECIES NAME\n");
	  exit(-1);
	}
	nayme[i - 1][j] = getc(infile);
      }
      j = 1;
      while (j <= sites && !eof(infile)) {
	if (eoln(infile)) {
	  fscanf(infile, "%*[^\n]");
	  getc(infile);
	}
	fscanf(infile, "%lf", &x);
	if ((unsigned)x > 1.0) {
	  printf("GENE FREQ OUTSIDE [0,1], species%3ld\n", i);
	  exit(-1);
	} else {
	  nodef[i - 1][j - 1] = x;
	  j++;
	}
      }
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    return;
  }
  basesread = 0;
  allread = false;
  while (!allread) {
    allread = true;
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    i = 1;
    while (i <= spp) {
      if ((interleaved && basesread == 0) || !interleaved) {
	for (j = 0; j < nmlngth; j++) {
	  if (eof(infile) || eoln(infile)){
	    printf("ERROR: END-OF-LINE OR END-OF-FILE");
	    printf(" IN THE MIDDLE OF A SPECIES NAME\n");
	    exit(-1);
	  }
	  nayme[i - 1][j] = getc(infile);
	}
      }
      j = interleaved ? basesread : 0;
      done = false;
      while (!done && !eof(infile)) {
	if (interleaved)
	  done = true;
	while (j < sites && !(eoln(infile) ||eof(infile))) {
	  charstate = getc(infile);
	  if (charstate == ' ' ||
	       (data == seqs && charstate >= '0' && charstate <= '9'))
	    continue;
	  uppercase(&charstate);
	  j++;
	  if (charstate == '.')
	    charstate = nodep[0][j - 1];
	  nodep[i - 1][j - 1] = charstate;
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
      if ((interleaved && j != basesnew) || ((!interleaved) && j != sites)){
	printf("ERROR: SEQUENCES OUT OF ALIGNMENT\n");
	exit(-1);}
      i++;
    }
    if (interleaved) {
      basesread = basesnew;
      allread = (basesread == sites);
    } else
      allread = (i > spp);
  }
  if (!printdata)
    return;
  if (data == genefreqs)
    m = (sites - 1) / 8 + 1;
  else
    m = (sites - 1) / 60 + 1;
  for (i = 1; i <= m; i++) {
    for (j = 0; j < (spp); j++) {
      for (k = 0; k < nmlngth; k++)
	putc(nayme[j][k], outfile);
      fprintf(outfile, "   ");
      if (data == genefreqs)
	l = i * 8;
      else
	l = i * 60;
      if (l > sites)
	l = sites;
      if (data == genefreqs)
	n = (i - 1) * 8;
      else
	n = (i - 1) * 60;
      for (k = n; k < l; k++) {
	if (data == genefreqs)
	  fprintf(outfile, "%8.5f", nodef[j][k]);
	else {
	  if (j + 1 > 1 && nodep[j][k] == nodep[0][k])
	    charstate = '.';
	  else
	    charstate = nodep[j][k];
	  putc(charstate, outfile);
	  if ((k + 1) % 10 == 0 && (k + 1) % 60 != 0)
	    putc(' ', outfile);
	
	}
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* inputdata */


Static Void doinput()
{
  /* reads the input data */
  getoptions();
  getnums();
  oldweight.steps    = (stepp)   Malloc(sites*sizeof(long));
  oldweight.dummy    = (double *)Malloc(sizeof(double));
  weight.steps       = (stepp)   Malloc(sites*sizeof(long));
  weight.dummy       = (double *)Malloc(sizeof(double));
  where.steps        = (stepp)   Malloc(loci*sizeof(long));
  where.dummy        = (double *)Malloc(sizeof(double));
  how_many.steps     = (stepp)   Malloc(loci*sizeof(long));
  how_many.dummy     = (double *)Malloc(sizeof(double));
  factor.steps       = (stepp)   Malloc(sites*sizeof(long));
  factor.dummy       = (double *)Malloc(sizeof(double));
  newwhere.steps     = (stepp)   Malloc(loci*sizeof(long));
  newwhere.dummy     = (double *)Malloc(sizeof(double));
  newhowmany.steps   = (stepp)   Malloc(loci*sizeof(long));
  newhowmany.dummy   = (double *)Malloc(sizeof(double));
  newerwhere.steps   = (stepp)   Malloc(loci*sizeof(long));
  newerwhere.dummy   = (double *)Malloc(sizeof(double));
  newerhowmany.steps = (stepp)   Malloc(loci*sizeof(long));
  newerhowmany.dummy = (double *)Malloc(sizeof(double));
  inputoptions();
  newerfactor.steps  = (stepp)   Malloc(loci*maxfactor*sizeof(long));
  newerfactor.dummy  = (double *)Malloc(sizeof(double));
  nayme              = (naym *)  Malloc(spp*sizeof(naym));
  inputdata();
}  /* doinput */

Local Void bootweights()
{
  /* sets up weights by resampling data */
  long i, j, k;
  double p, q, r;

  ws = newgroups;
  for (i = 0; i < (ws); i++)
    weight.steps[i] = 0;
  if (jackknife) {
    if (newgroups & 1) {
      if (randum(seed) < 0.5)
	q = (newgroups - 1.0) / 2;
      else
	q = (newgroups + 1.0) / 2;
    } else
      q = newgroups / 2.0;
    r = newgroups;
    p = q / r;
    ws = 0;
    for (i = 0; i < (newgroups); i++) {
      if (randum(seed) < p) {
	weight.steps[i]++;
	ws++;
	q--;
      }
      r--;
      if (i + 1 < newgroups)
	p = q / r;
    }
  } else if (permute) {
    for (i = 0; i < (newgroups); i++)
      weight.steps[i] = 1;
  } else {
    for (i = 1; i <= (newgroups); i++) {
      j = (long)(newgroups * randum(seed)) + 1;
      weight.steps[j - 1]++;
    }
  }

  for (i = 0; i < (newgroups); i++)
    newerwhere.steps[i] = 0;
  for (i = 0; i < (newgroups); i++)
    newerhowmany.steps[i] = 0;
  newergroups = 0;
  newersites  = 0;

  for (i = 0; i < (newgroups); i++) {
    for (j = 1; j <= (weight.steps[i]); j++) {
      newergroups++;
      for (k = 1; k <= (newhowmany.steps[i]); k++) {
        newersites++;
        newerfactor.steps[newersites - 1] = newergroups;
      }
      newerwhere.steps[newergroups - 1] = newwhere.steps[i];
      newerhowmany.steps[newergroups - 1] = newhowmany.steps[i];
    }
  }
}  /* bootweights */

void sppermute(n)
long n;
{ 

  long i, j, k;
  for (i = 1; i <= (spp - 1); i++) {
    k = (long)((i+1)* randum(seed));
    j = sppord[n - 1][i];
    sppord[n - 1][i] = sppord[n - 1][k];
    sppord[n - 1][k] = j;
  }
}  /* sppermute */


void writedata()
{
  /* write out one set of bootstrapped sequences */
  long i, j, k, l, m, n, n2;
  double x;
  Char charstate;

  sppord = (long **)Malloc(newergroups*sizeof(long *));
  for (i = 0; i < (newergroups); i++)
    sppord[i] = (long *)Malloc(spp*sizeof(long));
  for (i = 0; i < (newergroups); i++) {
    for (j = 1; j <= (spp); j++)
      sppord[i][j - 1] = j;
  }
  if (data == restsites && enzymes)
    fprintf(outfile, "%5ld %5ld %4ld\n", spp, newergroups, nenzymes);
  else if (data == genefreqs)
    fprintf(outfile, "%5ld %5ld\n", spp, newergroups);
  else
    fprintf(outfile, "%5ld %5ld\n", spp, newersites);
  if (data == genefreqs) {
    for (i = 0; i < (newergroups); i++)
      fprintf(outfile, "%3ld", alleles[factor.steps[newerwhere.steps[i] - 1] - 1]);
    putc('\n', outfile);
  }
  l = 1;
  if (interleaved)
    m = 60;
  else
    m = newergroups;
  do {
    if (m > newergroups)
      m = newergroups;
    for (j = 0; j < (spp); j++) {
      n = 0;
      if (l == 1) {
	for (k = 0; k < nmlngth; k++)
	  putc(nayme[j][k], outfile);
      } else {
	for (k = 1; k <= nmlngth; k++)
	  putc(' ', outfile);
      }
      fprintf(outfile, "   ");
      for (k = l - 1; k < m; k++) {
	if (permute && j + 1 == 1)
	  sppermute(newerfactor.steps[n]);
	for (n2 = -1; n2 <= (newerhowmany.steps[k] - 2); n2++) {
	  n++;
	  if (data == genefreqs) {
	    if (n > 1 && (n & 7) == 1)
	      fprintf(outfile, "\n              ");
	    x = nodef[sppord[newerfactor.steps[n - 1] - 1][j] - 1][newerwhere.steps[k] + n2];
	    fprintf(outfile, "%8.5f", x);
	  } else {
	    if (!interleaved && n > 1 && n % 60 == 1)
	      fprintf(outfile, "\n             ");
	    charstate =
	      nodep[sppord[newerfactor.steps[n - 1] - 1][j] - 1][newerwhere.steps[k] + n2];
	    putc(charstate, outfile);
	    if (n % 10 == 0 && n % 60 != 0)
	      putc(' ', outfile);
	  }
	}
      }
      putc('\n', outfile);
    }
    if (interleaved && i <= (newsites - 1) / 60)
      putc('\n', outfile);
    if (interleaved) {
      if (m < newersites)
	putc('\n', outfile);
      l += 60;
      m += 60;
    }
  } while (interleaved && l <= newersites);
  for (i = 0; i < (newergroups); i++)
    free(sppord[i]);
  free(sppord);
}  /* writedata*/


Static Void bootwrite()
{
  /* does bootstrapping and writes out data sets */
  long rr, repdiv10;

  repdiv10 = reps / 10;
  if (repdiv10 < 1)
    repdiv10 = 1;
  if (progress)
    putchar('\n');
  for (rr = 1; rr <= (reps); rr++) {
    bootweights();
    writedata();
    if (progress && rr % repdiv10 == 0)
      printf("completed replicate number %4ld\n", rr);
  }
  if (progress)
    printf("\nOutput written to output file\n\n");
}  /* bootwrite */


main(argc, argv)
int argc;
Char *argv[];
{  /* Read in sequences or frequencies and bootstrap or jackknife them */
char infilename[100],outfilename[100];
#ifdef MAC
  macsetup("Seqboot","");
  argv[0] = "Seqboot";
#endif
  openfile(&infile,INFILE,"r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);

  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  doinput();
  bootwrite();
  FClose(infile);
  FClose(outfile);
#ifdef MAC
  fixmacfile(outfilename);
#endif
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
  mem = (MALLOCRETURN *)calloc(1,x);
  if (!mem)
    {
      memerror();
      /* Just to make -Wall happy, */
      return (MALLOCRETURN *)mem;
    }
  else
    return (MALLOCRETURN *)mem;
}


