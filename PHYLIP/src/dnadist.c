#include "phylip.h"

/* version 3.56c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define maxcategs       9    /* maximum number of site types */
#define iterations      100  /* number of iterates of EM for each distance */
#define nmlngth         10   /* number of characters max. in species name */

#define epsilon         0.00001   /* used in makedist */

#define ibmpc0          false
#define ansi0           true
#define vt520           false


typedef enum {
  A, C, G, T
} base;
typedef double sitelike[(long)T - (long)A + 1];
typedef sitelike *phenotype;
typedef Char naym[nmlngth];

typedef Char** sequence;
typedef struct node {
  phenotype x;
  naym nayme;
} node;

typedef struct valrec {
  double rat, ratxv, z1, y1, z1zz, z1yy, z1xv;
} valrec;


Static FILE *infile, *outfile;
Static short j, numsp, sites, endsite, categs, weightsum, datasets, ith;
Static boolean freqsfrom, jukes, jinnei, kimura, lower, ml, weights,
               printdata, progress, ctgry, mulsets, firstset, interleaved,
               ibmpc, vt52, ansi;
Static node **nodep;
Static double xi, xv, ttratio, ttratio0, freqa, freqc, freqg, freqt, freqr,
	      freqy, freqar, freqcy, freqgr, freqty, fracchange, sumrates, cvi;
Static short *category, *oldweight, *weight, *alias, *ally, *location;
Static double rate[maxcategs];
Static double **d;
Static sequence y;
double sumweightrat;                  /* these values were propogated  */
double *weightrat;                    /* to global values from the     */
valrec tbl[maxcategs];                /* procedure makedists.          */


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


void uppercase(ch)
Char *ch;
{  /* convert ch to upper case -- either ASCII or EBCDIC */
   *ch = (islower(*ch) ?  toupper(*ch) : (*ch));
}  /* uppercase */


void getnums()
{
  /* input number of species, number of sites */
  fscanf(infile, "%hd%hd", &numsp, &sites);
}  /* getnums */

void getoptions()
{
  /* interactively set options */
  short j,i;
  Char ch;
  boolean  done1, ttr;
  char line[256];
  char rest[256];
  int   scanned;

  ctgry = false;
  categs = 1;
  rate[0] = 1.0;
  freqsfrom = false;
  jinnei = false;
  jukes = false;
  kimura = true;
  ml = false;
  lower = false;
  ttratio = 2.0;
  ttr = false;
  weights = false;
  printdata = false;
  progress = true;
  interleaved = true;
  for (;;) {
    printf(ansi ? "\033[2J\033[H" :
	   vt52 ? "\033E\033H"    : "\n");
    printf("\nNucleic acid sequence Distance Matrix program,");
    printf(" version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  D  Distance (Kimura, Jin/Nei, ML, J-C)?  %s\n",
	   kimura ? "Kimura 2-parameter" :
           jinnei ? "Jin and Nei"        :
           jukes  ? "Jukes-Cantor"       : "Maximum Likelihood");
    if (kimura || jinnei || ml) {
      printf("  T        Transition/transversion ratio?");
      if (!ttr)
        printf("  2.0\n");
      else
        printf("  %8.4f\n", ttratio);
    }
    printf("  C   One category of substitution rates?");
    if (!ctgry || categs == 1)
      printf("  Yes\n");
    else
      printf("  %hd categories\n", categs);
    if (ml) {
      printf("  F       Use empirical base frequencies?  %s\n",
	     (freqsfrom ? "Yes" : "No"));
    }
    printf("  L              Form of distance matrix?  %s\n",
	   (lower ? "Lower-triangular" : "Square"));
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
    printf("\nAre these settings correct? (type Y or letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    if  (ch == 'Y')
	   break;
   if (ch == 'C' || ch == 'F' || ch == 'L' || ch == 'I' || ch == 'D' ||
       ch == 'T' || ch == 'M' || ch == '0' || ch == '1' || ch == '2') {
     switch (ch) {

     case 'D':
       if (kimura) {
	 kimura = false;
         jinnei = true;
         freqsfrom = false;
       } else if (jinnei) {
         jinnei = false;
	 ml = true;
	 freqsfrom = true;
       } else if (ml) {
	 ml = false;
	 jukes = true;
	 freqsfrom = false;
       } else {
	 jukes = false;
	 kimura = true;
	 freqsfrom = false;
       }
       break;

     case 'C':
       ctgry = !ctgry;
       if (ctgry) {
	 do {
	   printf("Number of categories?\n");
	   gets(line);
	   categs = (short)atoi(line);
	 } while (categs < 1 || categs > maxcategs);
	 for (;;){
	   printf("Rate for each category? (use a space to separate)\n");
	   gets(line);
	   done1 = true;
	   for (i = 0; i < categs; i++){
	     scanned = sscanf(line,"%lf %[^\n]", &rate[i],rest);
	     if ((scanned != 2 && i < (categs - 1)) ||
		 (scanned != 1 && i == (categs - 1))){
	       printf("Please enter exactly %hd values.\n",categs);
	       done1 = false;
	       break;
	     }
	     strcpy(line,rest);
	   }
	   if (done1)
	     break;
	 }
       }
       break;

     case 'F':
       freqsfrom = !freqsfrom;
       if (!freqsfrom) {
	 printf(
		"Base frequencies for A, C, G, T/U (use blanks to separate)?\n");
	 scanf("%lf%lf%lf%lf%*[^\n]", &freqa, &freqc, &freqg, &freqt);
	 getchar();
       }
       break;

     case 'L':
       lower = !lower;
       break;

     case 'T':
       if (jukes) {
	 printf("WARNING: CANNOT SET TRANSITION/");
	 printf("TRANSVERSION RATIO FOR");
	 printf(" JUKES-CANTOR DISTANCE\n");
       } else {
	 ttr = !ttr;
	 if (ttr) {
	   do {
	     printf("Transition/transversion ratio?\n");
	     scanf("%lf%*[^\n]", &ttratio);
	     getchar();
	   } while (ttratio < 0.0);
	 }
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
  if (jinnei) {
    do {
      printf(
	"\nCoefficient of variation of substitution rate among sites (must be positive)\n");
      scanf("%lf%*[^\n]", &cvi);
      getchar();
    } while (cvi <= 0.0);
    cvi = 1.0 / (cvi * cvi);
  }
    if (!printdata)
    return;
    fprintf(outfile, "\nNucleic acid sequence Distance Matrix program,");
    fprintf(outfile, " version %s\n\n",VERSION);
    fprintf(outfile, "%3hd species, %4hd sites\n", numsp, sites);
}  /* getoptions */


void doinit()
{
  /* initializes variables */
  short i;

  getnums();
  getoptions();
  y = (Char **)Malloc(numsp*sizeof(Char *));
  nodep = (node **)Malloc(numsp*sizeof(node *));
  for (i = 0; i < numsp; i++) {
    y[i] = (Char *)Malloc(sites*sizeof(Char));
    nodep[i] = (node *)Malloc(sizeof(node));
  }
}  /* doinit */

void inputcategories()
{
  /* reads the categories for each site */
  short i;
  Char ch;

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
    category[i] = ch - '0';
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
}  /* inputcategories */

void printcategories()
{
  /* print out list of categories of sites */
  short i, j;

  fprintf(outfile, "Rate categories\n\n");
  for (i = 1; i <= nmlngth + 3; i++)
    putc(' ', outfile);
  for (i = 1; i <= sites; i++) {
    fprintf(outfile, "%hd", category[i - 1]);
    if (i % 60 == 0) {
      putc('\n', outfile);
      for (j = 1; j <= nmlngth + 3; j++)
        putc(' ', outfile);
    } else if (i % 10 == 0)
      putc(' ', outfile);
  }
  fprintf(outfile, "\n\n");

}  /* printcategories */

void inputweights()
{
  /* input the character weights, 0 or 1 */
  Char ch;
  short i;

  for (i = 1; i < nmlngth; i++)
    ch = getc(infile);
  weightsum = 0;
  for (i = 0; i < sites; i++) {
    do {
      if (eoln(infile)) {
        fscanf(infile, "%*[^\n]");
        getc(infile);
      }
      ch = getc(infile);
    } while (ch == ' ');
    oldweight[i] = 1;
    if (ch == '0' || ch == '1')
      oldweight[i] = ch - '0';
    else {
      printf("BAD WEIGHT CHARACTER: %c -- ");
      printf("WEIGHTS IN DNADIST MUST BE 0 OR 1\n",ch);
      exit(-1);
    }
    weightsum += oldweight[i];
  }
  weights = true;
  fscanf(infile, "%*[^\n]");
  getc(infile);
}  /* inputweights */

void printweights()
{
  /* print out the weights of sites */
  short i, j;

  fprintf(outfile, "\nSites are weighted as follows:\n");
  for (i = 1; i <= sites; i++) {
    if ((i - 1) % 60 == 0) {
      putc('\n', outfile);
      for (j = 1; j <= nmlngth + 3; j++)
        putc(' ', outfile);
    }
    fprintf(outfile, "%hd", oldweight[i - 1]);
    if (i % 10 == 0 && i % 60 != 0)
      putc(' ', outfile);
  }
  fprintf(outfile, "\n\n");
}  /* printweights */

void inputoptions()
{
  /* read options information */
  Char ch;
  boolean ctg;
  short i, extranum, cursp, cursts;

  if (!firstset) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%hd%hd", &cursp, &cursts);
    if (cursp != numsp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4hd\n", ith);
      exit(-1);
    }
    sites = cursts;
  }
  for (i = 0; i < sites; i++){
    category[i] = 1;
    oldweight[i]= 1;
  }
  weightsum = sites;
  ctg = false;
  extranum = 0;
  while (!(eoln(infile))) {
    ch = getc(infile);
    uppercase(&ch);
    if (ch == 'C' || ch == 'W')
      extranum++;
    else if (ch != ' ') {
      putc('\n', outfile);
      printf("MAKE SURE YOU WANTED OPTION CHARACTER: %c\n", ch);
      printf("(It was passed along to output file)\n");
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  if (printdata)
    putc('\n', outfile);
  if (jukes && printdata)
    fprintf(outfile, "  Jukes-Cantor Distance\n");
  if (kimura && printdata)
    fprintf(outfile, "  Kimura 2-parameter Distance\n");
  if (jinnei && printdata)
    fprintf(outfile, "  Jin and Nei Distance\n");
  if (ml && printdata)
    fprintf(outfile, "  Maximum Likelihood Distance\n");
  for (i = 1; i <= extranum; i++) {
    ch = getc(infile);
    uppercase(&ch);
    if (ch != 'C' && ch != 'W') {
      printf("ERROR: INCORRECT AUXILIARY OPTIONS");
      printf(" LINE WHICH STARTS WITH %c\n", ch);
    }
    if (ch == 'C') {
      ctg = true;
      if (!ctgry || categs <= 1) {
        printf("ERROR: CATEGORY OPTION NOT CHOSEN IN MENU");
	printf(" WITH OPTION %c IN INPUT\n",ch);
	exit(-1);
      } else
        inputcategories();
    }
    if (ch == 'W')
      inputweights();
  }
  if ((categs > 1 || ctgry) && !ctg) {
    printf("ERROR: CATEGORY OPTION CHOSEN IN MENU");
    printf(" WITH NO OPTION C IN INPUT\n");
    exit(-1);
  } else if (printdata && (categs > 1)) {
    fprintf(outfile, "\nSite category   Rate of change\n\n");
    for (i = 1; i <= categs; i++)
      fprintf(outfile, "%12hd%13.3f\n", i, rate[i - 1]);
    putc('\n', outfile);
    printcategories();
  }
  if ((jukes || kimura || jinnei) && freqsfrom) {
    printf("WARNING: CANNOT USE EMPIRICAL BASE FREQUENCIES");
    printf(" WITH JUKES-CANTOR, KIMURA OR JIN/NEI DISTANCES\n");
    exit(-1);
  }
  if (jukes)
    ttratio = 0.5;
  if (weights && printdata)
    printweights();
}  /* inputoptions */

void getbasefreqs()
{
  /* compute or read in base frequencies */
  double aa, bb;

  if (printdata) {
    putc('\n', outfile);
    if (freqsfrom)
      fprintf(outfile, "Empirical");
    fprintf(outfile, " Base Frequencies:\n");
  }
  if (progress)
    putchar('\n');
  if (!freqsfrom) {
    if (kimura || jukes || jinnei) {
      freqa = 0.25;
      freqc = 0.25;
      freqg = 0.25;
      freqt = 0.25;
    }
  }
  if (printdata) {
    fprintf(outfile, "    A    %10.5f\n", freqa);
    fprintf(outfile, "    C    %10.5f\n", freqc);
    fprintf(outfile, "    G    %10.5f\n", freqg);
    fprintf(outfile, "   T(U)  %10.5f\n", freqt);
  }
  freqr = freqa + freqg;
  freqy = freqc + freqt;
  freqar = freqa / freqr;
  freqcy = freqc / freqy;
  freqgr = freqg / freqr;
  freqty = freqt / freqy;
  if (printdata)
    fprintf(outfile, "\nTransition/transversion ratio = %10.6f\n\n", ttratio);
  aa = ttratio * freqr * freqy - freqa * freqg - freqc * freqt;
  bb = freqa * freqgr + freqc * freqty;
  xi = aa / (aa + bb);
  xv = 1.0 - xi;
  ttratio = xi / xv;
  if (xi <= 0.0 && xi >= -epsilon)
    xi = 0.0;
  if (xi < 0.0) {
    printf("WARNING: This transition/transversion ratio\n");
    printf("is impossible with these base frequencies!\n");
    xi = 3.0 / 5;
    xv = 2.0 / 5;
    printf(" Transition/transversion parameter reset\n\n");
  }
  if (printdata)
    fprintf(outfile, "(Transition/transversion parameter = %10.6f)\n\n",
            xi / xv);
  fracchange = xi * (2 * freqa * freqgr + 2 * freqc * freqty) +
      xv * (1.0 - freqa * freqa - freqc * freqc - freqg * freqg - freqt * freqt);
}  /* getbasefreqs */

void inputdata()
{
  /* Input the names and sequences for each species */
  short i, j, k, l, basesread, basesnew;
  Char charstate;
  boolean allread, done;

  if (printdata)
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
  while (!(allread)) {
    allread = true;
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    i = 1;
    while (i <= numsp) {
      if ((interleaved && basesread == 0 )|| !interleaved) {
        for (j = 0; j < nmlngth; j++) {
	  if (eof(infile) || eoln(infile)){
	    printf("ERROR: END-OF-LINE OR END-OF-FILE IN");
            printf(" THE MIDDLE OF A SPECIES NAME\n");
	    exit(-1);}
	  nodep[i - 1]->nayme[j] = getc(infile);
        }
      }
      j = (interleaved) ? basesread : 0;
      done = false;
      while (!done && !eof(infile)) {
        if (interleaved)
          done = true;
        while (j < sites && !(eoln(infile) || eof(infile))) {
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
      allread = (i > numsp);
  }
  if (!printdata)
    return;
  for (i = 1; i <= ((sites -1) / 60 + 1); i++) {
    for (j = 1; j <= numsp; j++) {
      for (k = 0; k < nmlngth; k++)
        putc(nodep[j - 1]->nayme[k], outfile);
      fprintf(outfile, "   ");
      l = i * 60;
      if (l > sites)
        l = sites;
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
  /* Shell sort of sites lexicographically */
  short gap, i, j, jj, jg, k, itemp;
  boolean flip, tied;

  gap = sites / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= sites; i++) {
      j = i - gap;
      flip = true;
      while (j > 0 && flip) {
        jj = alias[j - 1];
        jg = alias[j + gap - 1];
        tied = (oldweight[jj - 1] == oldweight[jg - 1]);
        flip = (oldweight[jj - 1] < oldweight[jg - 1] ||
                (tied && category[jj - 1] > category[jg - 1]));
        tied = (tied && category[jj - 1] == category[jg - 1]);
        k = 1;
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
  while (i < sites) {
    j = i + 1;
    tied = true;
    while (j <= sites && tied) {
      tied = (oldweight[alias[i - 1] - 1] == oldweight[alias[j - 1] - 1] &&
              category[alias[i - 1] - 1] == category[alias[j - 1] - 1]);
      k = 1;
      while (k <= numsp && tied) {
        tied = (tied &&
            y[k - 1][alias[i - 1] - 1] == y[k - 1][alias[j - 1] - 1]);
        k++;
      }
      if (!tied)
        break;
      ally[alias[j - 1] - 1] = alias[i - 1];
      j++;
    }
    i = j;
  }
}  /* sitecombine */

void sitescrunch()
{
  /* move so one representative of each pattern of
     sites comes first */
  short i, j, itemp;
  boolean done, found, completed;

  done = false;
  i = 1;
  j = 2;
  while (!done) {
    if (ally[alias[i - 1] - 1] != alias[i - 1]) {
      if (j <= i)
        j = i + 1;
      if (j <= sites) {
        found = false;
        do {
          found = (ally[alias[j - 1] - 1] == alias[j - 1]);
          j++;
          completed = (j > sites);
          if (j <= sites)
            completed = (oldweight[alias[j - 1] - 1] == 0);
        } while (!(found || completed));
        if (found) {
          j--;
          itemp = alias[i - 1];
          alias[i - 1] = alias[j - 1];
          alias[j - 1] = itemp;
        } else
          done = true;
      } else
        done = true;
    }
    i++;
    done = (done || i >= sites);
  }
}  /* sitescrunch */

void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  short i;

  for (i = 1; i <= sites; i++) {
    alias[i - 1] = i;
    ally[i - 1] = i;
    weight[i - 1] = 0;
  }
  sitesort();
  sitecombine();
  sitescrunch();
  endsite = 0;
  for (i = 1; i <= sites; i++) {
    if (ally[i - 1] == i && oldweight[i - 1] > 0)
      endsite++;
  }
  for (i = 1; i <= endsite; i++)
    location[alias[i - 1] - 1] = i;
  sumrates = 0.0;
  for (i = 0; i < sites; i++)
    sumrates += oldweight[i] * rate[category[i] - 1];
  for (i = 0; i < categs; i++)
    rate[i] *= weightsum / sumrates;
  sumrates = weightsum;
  for (i = 0; i < sites; i++)
    weight[location[ally[i] - 1] - 1] += oldweight[i];
}  /* makeweights */

void makevalues()
{
  /* set up fractional likelihoods at tips */
  short i, j, k;
  base b;

  for (i = 0; i < numsp; i++)
    nodep[i]->x = (phenotype)Malloc(endsite*sizeof(sitelike));
  for (k = 0; k < endsite; k++) {
    j = alias[k];
    for (i = 0; i < numsp; i++) {
      for (b = A; (long)b <= (long)T; b = (base)((long)b + 1))
        nodep[i]->x[k][(long)b - (long)A] = 0.0;
      switch (y[i][j - 1]) {

      case 'A':
        nodep[i]->x[k][0] = 1.0;
        break;

      case 'C':
        nodep[i]->x[k][(long)C - (long)A] = 1.0;
        break;

      case 'G':
        nodep[i]->x[k][(long)G - (long)A] = 1.0;
        break;

      case 'T':
        nodep[i]->x[k][(long)T - (long)A] = 1.0;
        break;

      case 'U':
        nodep[i]->x[k][(long)T - (long)A] = 1.0;
        break;

      case 'M':
        nodep[i]->x[k][0] = 1.0;
        nodep[i]->x[k][(long)C - (long)A] = 1.0;
        break;

      case 'R':
        nodep[i]->x[k][0] = 1.0;
        nodep[i]->x[k][(long)G - (long)A] = 1.0;
        break;

      case 'W':
        nodep[i]->x[k][0] = 1.0;
        nodep[i]->x[k][(long)T - (long)A] = 1.0;
        break;

      case 'S':
        nodep[i]->x[k][(long)C - (long)A] = 1.0;
        nodep[i]->x[k][(long)G - (long)A] = 1.0;
        break;

      case 'Y':
        nodep[i]->x[k][(long)C - (long)A] = 1.0;
        nodep[i]->x[k][(long)T - (long)A] = 1.0;
        break;

      case 'K':
        nodep[i]->x[k][(long)G - (long)A] = 1.0;
        nodep[i]->x[k][(long)T - (long)A] = 1.0;
        break;

      case 'B':
        nodep[i]->x[k][(long)C - (long)A] = 1.0;
        nodep[i]->x[k][(long)G - (long)A] = 1.0;
        nodep[i]->x[k][(long)T - (long)A] = 1.0;
        break;

      case 'D':
        nodep[i]->x[k][0] = 1.0;
        nodep[i]->x[k][(long)G - (long)A] = 1.0;
        nodep[i]->x[k][(long)T - (long)A] = 1.0;
        break;

      case 'H':
        nodep[i]->x[k][0] = 1.0;
        nodep[i]->x[k][(long)C - (long)A] = 1.0;
        nodep[i]->x[k][(long)T - (long)A] = 1.0;
        break;

      case 'V':
        nodep[i]->x[k][0] = 1.0;
        nodep[i]->x[k][(long)C - (long)A] = 1.0;
        nodep[i]->x[k][(long)G - (long)A] = 1.0;
        break;

      case 'N':
        for (b = A; (long)b <= (long)T; b = (base)((long)b + 1))
          nodep[i]->x[k][(long)b - (long)A] = 1.0;
        break;

      case 'X':
        for (b = A; (long)b <= (long)T; b = (base)((long)b + 1))
          nodep[i]->x[k][(long)b - (long)A] = 1.0;
        break;

      case '?':
        for (b = A; (long)b <= (long)T; b = (base)((long)b + 1))
          nodep[i]->x[k][(long)b - (long)A] = 1.0;
        break;

      case 'O':
        for (b = A; (long)b <= (long)T; b = (base)((long)b + 1))
          nodep[i]->x[k][(long)b - (long)A] = 1.0;
        break;

      case '-':
        for (b = A; (long)b <= (long)T; b = (base)((long)b + 1))
          nodep[i]->x[k][(long)b - (long)A] = 1.0;
        break;
      }
    }
  }
}  /* makevalues */

void empiricalfreqs()
{
  /* Get empirical base frequencies from the data */
  short i, j, k;
  double sum, suma, sumc, sumg, sumt, w;

  freqa = 0.25;
  freqc = 0.25;
  freqg = 0.25;
  freqt = 0.25;
  for (k = 1; k <= 8; k++) {
    suma = 0.0;
    sumc = 0.0;
    sumg = 0.0;
    sumt = 0.0;
    for (i = 0; i < numsp; i++) {
      for (j = 0; j < endsite; j++) {
        w = weight[j];
        sum = freqa * nodep[i]->x[j][0];
        sum += freqc * nodep[i]->x[j][(long)C - (long)A];
        sum += freqg * nodep[i]->x[j][(long)G - (long)A];
        sum += freqt * nodep[i]->x[j][(long)T - (long)A];
        suma += w * freqa * nodep[i]->x[j][0] / sum;
        sumc += w * freqc * nodep[i]->x[j][(long)C - (long)A] / sum;
        sumg += w * freqg * nodep[i]->x[j][(long)G - (long)A] / sum;
        sumt += w * freqt * nodep[i]->x[j][(long)T - (long)A] / sum;
      }
    }
    sum = suma + sumc + sumg + sumt;
    freqa = suma / sum;
    freqc = sumc / sum;
    freqg = sumg / sum;
    freqt = sumt / sum;
  }
}  /* empiricalfreqs */


void getinput()
{
  /* reads the input data */
  inputoptions();
  if (!freqsfrom )
    getbasefreqs();
  inputdata();
  makeweights();
  makevalues();
  if (freqsfrom) {
    empiricalfreqs();
    getbasefreqs();
  }
}  /* getinput */



void inittable()
{
  /* Define a lookup table. Precompute values and store in a table */
  short i;

  for (i = 0; i < categs; i++) {
    tbl[i].rat = rate[i];
    tbl[i].ratxv = rate[i] * xv;
  }
}  /* inittable */


void makev(m, n, v)
short m, n;
double *v;
{
  /* compute one distance */
  short i, it, numerator, denominator, num1, num2, idx;
  double sum, sum1, sum2, sumyr, lz, aa, bb, cc, vv, p1, p2, p3, q1, q2, q3,
	 tt, delta, slope, xx1freqa, xx1freqc, xx1freqg, xx1freqt;
  double *prod, *prod2, *prod3;
  boolean quick, jukesquick, kimquick, jinneiquick;
  base b;
  node *p, *q;
  sitelike xx1, xx2;

  p = nodep[m - 1];
  q = nodep[n - 1];
  quick = (!ctgry || categs == 1);
  if (jukes || kimura || jinnei) {
    numerator = 0;
    denominator = 0;
    for (i = 0; i < endsite; i++) {
      memcpy(xx1, p->x[i], sizeof(sitelike));
      memcpy(xx2, q->x[i], sizeof(sitelike));
      sum = 0.0;
      sum1 = 0.0;
      sum2 = 0.0;
      for (b = A; (long)b <= (long)T; b = (base)((long)b + 1)) {
        sum1 += xx1[(long)b - (long)A];
        sum2 += xx2[(long)b - (long)A];
        sum += xx1[(long)b - (long)A] * xx2[(long)b - (long)A];
      }
      quick = (quick && (sum1 == 1.0 || sum1 == 4.0) &&
               (sum2 == 1.0 || sum2 == 4.0));
      if (sum1 == 1.0 && sum2 == 1.0) {
        numerator += (long)(weight[i] * sum);
        denominator += weight[i];
      }
    }
  }
  jukesquick = (jukes && quick);
  kimquick = (kimura && quick);
  jinneiquick = (jinnei && quick);
  if (jinnei && !jinneiquick) {
    printf("WARNING: CANNOT CALCULATE JIN/NEI DISTANCE\n");
    printf(" WITH PRESENT PROGRAM IF PARTIALLY AMBIGUOUS NUCLEOTIDES\n");
    exit(-1);
  }
  if (jukesquick && numerator * 4 <= denominator) {
    printf(" WARNING: INFINITE DISTANCE BETWEEN ");
    printf(" SPECIES %3hd AND %3hd\n", m, n);
    exit(-1);
  }
  if (jukesquick)
    vv = -0.75 * log((4.0 * ((double)numerator / denominator) - 1.0) / 3.0);
  if (kimquick || jinneiquick) {
    num1 = 0;
    num2 = 0;
    denominator = 0;
    for (i = 0; i < endsite; i++) {
      memcpy(xx1, p->x[i], sizeof(sitelike));
      memcpy(xx2, q->x[i], sizeof(sitelike));
      sum = 0.0;
      sum1 = 0.0;
      sum2 = 0.0;
      for (b = A; (long)b <= (long)T; b = (base)((long)b + 1)) {
        sum1 += xx1[(long)b - (long)A];
        sum2 += xx2[(long)b - (long)A];
        sum += xx1[(long)b - (long)A] * xx2[(long)b - (long)A];
      }
      sumyr = (xx1[0] + xx1[(long)G - (long)A])
            * (xx2[0] + xx2[(long)G - (long)A]) +
              (xx1[(long)C - (long)A] + xx1[(long)T - (long)A]) *
              (xx2[(long)C - (long)A] + xx2[(long)T - (long)A]);
      if (sum1 == 1.0 && sum2 == 1.0) {
        num1 += (long)(weight[i] * sum);
        num2 += (long)(weight[i] * (sumyr - sum));
        denominator += weight[i];
      }
    }
    tt = 1.0 - (double)num1 / denominator;
    if (tt > 0.0) {
      delta = 0.1;
      tt = delta;
      it = 0;
      while (fabs(delta) > 0.00002 && it < iterations) {
	it++;
	if (kimura) {
          p1 = exp(-tt);
          p2 = exp(-xv * tt) - exp(-tt);
          p3 = 1.0 - exp(-xv * tt);
	} else {
	  p1 = exp(-cvi * log(1 + tt / cvi));
	  p2 = exp(-cvi * log(1 + xv * tt / cvi))
              - exp(-cvi * log(1 + tt / cvi));
	  p3 = 1.0 - exp(-cvi * log(1 + xv * tt / cvi));
	}
        q1 = p1 + p2 / 2.0 + p3 / 4.0;
        q2 = p2 / 2.0 + p3 / 4.0;
	q3 = p3 / 2.0;
        if (kimura)
	  slope = 0.5 * exp(-tt) * (num2 / q2 - num1 / q1) +
		  0.25 * xv * exp(-xv * tt) *
		 ((denominator - num1 - num2) * 2 / q3 - num2 / q2 - num1 / q1);
	else
	  slope = 0.5 * (1 / (1 + tt / cvi)) * exp(-cvi * log(1 + tt / cvi)) *
		  (num2 / q2 - num1 / q1) + 0.25 * (xv / (1 + xv * tt / cvi)) *
		    exp(-cvi * log(1 + xv * tt / cvi)) *
		 ((denominator - num1 - num2) * 2 / q3 - num2 / q2 - num1 / q1);
	if (slope < 0.0)
	  delta = fabs(delta) / -2.0;
	else
	  delta = fabs(delta);
	tt += delta;
      }
    }
    vv = fracchange * tt;
  }
  if (!(jukesquick || kimquick || jinneiquick)) {
    prod = (double *)Malloc(sites*sizeof(double));
    prod2 = (double *)Malloc(sites*sizeof(double));
    prod3 = (double *)Malloc(sites*sizeof(double));
    for (i = 0; i < endsite; i++) {
      memcpy(xx1, p->x[i], sizeof(sitelike));
      memcpy(xx2, q->x[i], sizeof(sitelike));
      xx1freqa = xx1[0] * freqa;
      xx1freqc = xx1[(long)C - (long)A] * freqc;
      xx1freqg = xx1[(long)G - (long)A] * freqg;
      xx1freqt = xx1[(long)T - (long)A] * freqt;
      sum1 = xx1freqa + xx1freqc + xx1freqg + xx1freqt;
      sum2 = freqa * xx2[0] + freqc * xx2[(long)C - (long)A] +
             freqg * xx2[(long)G - (long)A] + freqt * xx2[(long)T - (long)A];
      prod[i] = sum1 * sum2;
      prod2[i] = (xx1freqa + xx1freqg) *
                 (xx2[0] * freqar + xx2[(long)G - (long)A] * freqgr) +
          (xx1freqc + xx1freqt) *
          (xx2[(long)C - (long)A] * freqcy + xx2[(long)T - (long)A] * freqty);
      prod3[i] = xx1freqa * xx2[0] + xx1freqc * xx2[(long)C - (long)A] +
         xx1freqg * xx2[(long)G - (long)A] + xx1freqt * xx2[(long)T - (long)A];
    }
    tt = 0.1;
    delta = 0.1;
    it = 1;
    while (it < iterations && fabs(delta) > 0.00002) {
      slope = 0.0;
      if (tt > 0.0) {
	lz = -tt;
	for (i = 0; i < categs; i++) {
	  tbl[i].z1 = exp(tbl[i].ratxv * lz);
	  tbl[i].y1 = 1.0 - tbl[i].z1;
	  tbl[i].z1zz = exp(tbl[i].rat * lz);
	  tbl[i].z1yy = tbl[i].z1 - tbl[i].z1zz;
	  tbl[i].z1xv = tbl[i].z1 * xv;
        }
        for (i = 0; i < endsite; i++) {
          idx = category[alias[i] - 1];
          cc = prod[i];
          bb = prod2[i];
          aa = prod3[i];
	  slope += weightrat[i] * (tbl[idx - 1].z1zz * (bb - aa) +
		tbl[idx - 1].z1xv * (cc - bb)) /
	      (aa * tbl[idx - 1].z1zz + bb * tbl[idx - 1].z1yy +
	       cc * tbl[idx - 1].y1);
        }
      }
      if (slope < 0.0)
	delta = fabs(delta) / -2.0;
      else
	delta = fabs(delta);
      tt += delta;
      it++;
    }
    vv = tt * fracchange;
    free(prod);
    free(prod2);
    free(prod3);
  }
  *v = vv;
}  /* makev */


void makedists()
{
  /* compute distance matrix */
  short i, j;
  double v;

  inittable();
  for (i = 0; i < endsite; i++)
    weightrat[i] = weight[i] * rate[category[alias[i] - 1] - 1];
  if (progress)
    printf("Distances calculated for species\n");
  for (i = 0; i < numsp; i++)
    d[i][i] = 0.0;
  for (i = 1; i < numsp; i++) {
    if (progress) {
      printf("    ");
      for (j = 0; j < nmlngth; j++)
        putchar(nodep[i - 1]->nayme[j]);
      printf("   ");
    }
    for (j = i + 1; j <= numsp; j++) {
      makev(i, j, &v);
      d[i - 1][j - 1] = v;
      d[j - 1][i - 1] = v;
      if (progress)
	putchar('.');
    }
    if (progress)
      putchar('\n');
  }
  if (progress) {
    printf("    ");
    for (j = 0; j < nmlngth; j++)
      putchar(nodep[numsp - 1]->nayme[j]);
    putchar('\n');
  }
  for (i = 0; i < numsp; i++)
    free(nodep[i]->x);
}  /* makedists */

void writedists()
{
  /* write out distances */
  short i, j, k;

  fprintf(outfile, "%5hd\n", numsp);
  for (i = 0; i < numsp; i++) {
    for (j = 0; j < nmlngth; j++)
      putc(nodep[i]->nayme[j], outfile);
    if (lower)
      k = i;
    else
      k = numsp;
    for (j = 1; j <= k; j++) {
      fprintf(outfile, "%8.4f", d[i][j - 1]);
      if ((j + 1) % 9 == 0 && j < k)
        putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  if (progress)
    printf("\nDistances written to file\n\n");
}  /* writedists */


main(argc, argv)
int argc;
Char *argv[];
{  /* DNA Distances by Maximum Likelihood */
char infilename[100],outfilename[100];
#ifdef MAC
  macsetup("Dnadist","");
  argv[0] = "Dnadist";
#endif
  openfile(&infile,INFILE,"r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);

  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  mulsets = false;
  datasets = 1;
  firstset = true;
  doinit();
  d = (double **)Malloc(numsp*sizeof(double *));
  for (j = 0; j < numsp; j++)
    d[j] = (double*)Malloc(numsp*sizeof(double));
  category = (short *)Malloc(sites*sizeof(short));
  oldweight = (short *)Malloc(sites*sizeof(short));
  weight = (short *)Malloc(sites*sizeof(short));
  alias = (short *)Malloc(sites*sizeof(short));
  ally = (short *)Malloc(sites*sizeof(short));
  location = (short *)Malloc(sites*sizeof(short));
  weightrat = (double *)Malloc(sites*sizeof(double));
  ttratio0 = ttratio;
  for (ith = 1; ith <= datasets; ith++) {
    ttratio = ttratio0;
    getinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1 && progress)
      printf("Data set # %hd:\n\n",ith);
    makedists();
    writedists();
  }
  FClose(infile);
  FClose(outfile);
#ifdef MAC
  fixmacfile(outfilename);
#endif
  exit(0);
}  /* DNA Distances by Maximum Likelihood */

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

