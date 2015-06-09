#include "phylip.h"

/* version 3.572c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define nmlngth         10    /* number of characters in species name     */
#define epsilon         0.000001   /* a small number */

#define ibmpc0          false
#define ansi0           true
#define vt520           false


typedef long *steparray;
typedef enum {
  ala, arg, asn, asp, cys, gln, glu, gly, his, ileu, leu, lys, met, phe, pro,
  ser, thr, trp, tyr, val, del, stop, asx, glx, unk, quest
} aas;
typedef enum {
  universal, ciliate, mito, vertmito, flymito, yeastmito
} codetype;
typedef enum {
  chemical, hall, george
} cattype;

typedef double matrix[20][20];


FILE *infile, *outfile;
long spp, chars, datasets,ith;
/* spp = number of species
   chars = number of sites in actual sequences */
double freqa, freqc, freqg, freqt, ttratio, xi, xv, ease, fracchange;
boolean weights, printdata, progress, mulsets, interleaved,
       ibmpc, vt52, ansi, basesequal, usepam, kimura, firstset;
codetype whichcode;
cattype whichcat;
steparray weight;
Char **nayme;
aas **node;
aas trans[4][4][4];
double pie[20];
long cat[(long)val - (long)ala + 1], numaa[(long)val - (long)ala + 1];
double eig[20];
matrix prob, eigvecs;
double **d;

/* Local variables for makedists, propagated globally for c version: */
  double tt, p, dp, d2p, q, elambdat;


static double pameigs[] = {-0.022091252,-0.019297602, 0.000004760,-0.017477817,
                           -0.016575549,-0.015504543,-0.002112213,-0.002685727,
                           -0.002976402,-0.013440755,-0.012926992,-0.004293227,
                           -0.005356688,-0.011064786,-0.010480731,-0.008760449,
                           -0.007142318,-0.007381851,-0.007806557,-0.008127024

			   };
static double pamprobs[20][20] =
{
      {-0.01522976,-0.00746819,-0.13934468, 0.11755315,-0.00212101,
       0.01558456,-0.07408235,-0.00322387, 0.01375826, 0.00448826,
       0.00154174, 0.02013313,-0.00159183,-0.00069275,-0.00399898,
       0.08414055,-0.01188178,-0.00029870, 0.00220371, 0.00042546},
      {-0.07765582,-0.00712634,-0.03683209,-0.08065755,-0.00462872,
       -0.03791039,0.10642147,-0.00912185,0.01436308,-0.00133243,
        0.00166346,0.00624657,-0.00003363,-0.00128729,-0.00690319,
        0.17442028,-0.05373336,-0.00078751,-0.00038151,0.01382718},
      {-0.08810973,-0.04081786,-0.04066004,-0.04736004,-0.03275406,
      -0.03761164,-0.05047487, -0.09086213,-0.03269598,-0.03558015,
      -0.08407966,-0.07970977, -0.01504743,-0.04011920,-0.05182232,
      -0.07026991,-0.05846931, -0.01016998,-0.03047472,-0.06280511},
      { 0.02513756,-0.00578333, 0.09865453, 0.01322314,-0.00310665,
        0.05880899,-0.09252443,-0.02986539,-0.03127460, 0.01007539,
       -0.00360119,-0.01995024, 0.00094940,-0.00145868,-0.01388816,
        0.11358341,-0.12127513,-0.00054696,-0.00055627, 0.00417284},
      { 0.16517316,-0.00254742,-0.03318745,-0.01984173, 0.00031890,
       -0.02817810, 0.02661678,-0.01761215, 0.01665112, 0.10513343,
       -0.00545026, 0.01827470,-0.00207616,-0.00763758,-0.01322808,
       -0.02202576,-0.07434204, 0.00020593, 0.00119979,-0.10827873},
      { 0.16088826, 0.00056313,-0.02579303,-0.00319655, 0.00037228,
       -0.03193150, 0.01655305,-0.03028640, 0.01367746,-0.11248153,
        0.00778371, 0.02675579, 0.00243718, 0.00895470,-0.01729803,
       -0.02686964,-0.08262584, 0.00011794,-0.00225134, 0.09415650},
      {-0.01739295, 0.00572017,-0.00712592,-0.01100922,-0.00870113,
       -0.00663461,-0.01153857,-0.02248432,-0.00382264,-0.00358612,
       -0.00139345,-0.00971460,-0.00133312, 0.01927783,-0.01053838,
       -0.00911362,-0.01010908, 0.09417598, 0.01763850,-0.00955454},
      { 0.01728888, 0.01344211, 0.01200836, 0.01857259,-0.17088517,
        0.01457592, 0.01997839, 0.02844884, 0.00839403, 0.00196862,
        0.01391984, 0.03270465, 0.00347173,-0.01940984, 0.01233979,
        0.00542887, 0.01008836, 0.00126491,-0.02863042, 0.00449764},
      {-0.02881366,-0.02184155,-0.01566086,-0.02593764,-0.04050907,
       -0.01539603,-0.02576729,-0.05089606,-0.00597430, 0.02181643,
        0.09835597,-0.04040940, 0.00873512, 0.12139434,-0.02427882,
       -0.02945238,-0.01566867,-0.01606503, 0.09475319, 0.02238670},
      { 0.04080274,-0.02869626,-0.05191093,-0.08435843, 0.00021141,
        0.13043842, 0.00871530, 0.00496058,-0.02797641,-0.00636933,
        0.02243277, 0.03640362,-0.05735517, 0.00196918,-0.02218934,
       -0.00608972, 0.02872922, 0.00047619, 0.00151285, 0.00883489},
      {-0.02623824, 0.00331152, 0.03640692, 0.04260231,-0.00038223,
       -0.07480340,-0.01022492,-0.00426473, 0.01448116, 0.01456847,
        0.05786680, 0.03368691,-0.10126924,-0.00147454, 0.01275395,
        0.00017574,-0.01585206,-0.00015767,-0.00231848, 0.02310137},
      {-0.00846258,-0.01508106,-0.01967505,-0.02772004, 0.01248253,
       -0.01331243,-0.02569382,-0.04461524,-0.02207075, 0.04663443,
        0.19347923,-0.02745691, 0.02288515,-0.04883849,-0.01084597,
       -0.01947187,-0.00081675, 0.00516540,-0.07815919, 0.08035585},
      {-0.06553111, 0.09756831, 0.00524326,-0.00885098, 0.00756653,
        0.02783099,-0.00427042,-0.16680359, 0.03951331,-0.00490540,
        0.01719610, 0.15018204, 0.00882722,-0.00423197,-0.01919217,
       -0.02963619,-0.01831342,-0.00524338, 0.00011379,-0.02566864},
      {-0.07494341,-0.11348850, 0.00241343,-0.00803016, 0.00492438,
        0.00711909,-0.00829147, 0.05793337, 0.02734209, 0.02059759,
       -0.02770280, 0.14128338, 0.01532479, 0.00364307, 0.05968116,
       -0.06497960,-0.08113941, 0.00319445,-0.00104222, 0.03553497},
      { 0.05948223,-0.08959930, 0.03269977,-0.03272374,-0.00365667,
       -0.03423294,-0.06418925,-0.05902138, 0.05746317,-0.02580596,
        0.01259572, 0.05848832, 0.00672666, 0.00233355,-0.05145149,
        0.07348503, 0.11427955, 0.00142592,-0.01030651,-0.04862799},
      {-0.01606880, 0.05200845,-0.01212967,-0.06824429,-0.00234304,
        0.01094203,-0.07375538, 0.08808629, 0.12394822, 0.02231351,
       -0.03608265,-0.06978045,-0.00618360, 0.00274747,-0.01921876,
       -0.01541969,-0.02223856,-0.00107603,-0.01251777, 0.05412534},
      { 0.01688843, 0.05784728,-0.02256966,-0.07072251,-0.00422551,
       -0.06261233,-0.08502830, 0.08925346,-0.08529597, 0.01519343,
       -0.05008258, 0.10931873, 0.00521033, 0.02593305,-0.00717855,
        0.02291527, 0.02527388,-0.00266188,-0.00871160, 0.02708135},
      {-0.04233344, 0.00076379, 0.01571257, 0.04003092, 0.00901468,
        0.00670577, 0.03459487, 0.12420216,-0.00067366,-0.01515094,
        0.05306642, 0.04338407, 0.00511287, 0.01036639,-0.17867462,
       -0.02289440,-0.03213205, 0.00017924,-0.01187362,-0.03933874},
      { 0.01284817,-0.01685622, 0.00724363, 0.01687952,-0.00882070,
       -0.00555957, 0.01676246,-0.05560456,-0.00966893, 0.06197684,
       -0.09058758, 0.00880607, 0.00108629,-0.08308956,-0.08056832,
       -0.00413297, 0.02973107, 0.00092948, 0.07010111, 0.13007418},
      { 0.00700223,-0.01347574, 0.00691332, 0.03122905, 0.00310308,
        0.00946862, 0.03455040,-0.06712536,-0.00304506, 0.04267941,
       -0.10422292,-0.01127831,-0.00549798, 0.11680505,-0.03352701,
       -0.00084536, 0.01631369, 0.00095063,-0.09570217, 0.06480321}
    };

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
        file[0]='\0';
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
{
 (*ch) = (isupper(*ch) ? (*ch) : toupper(*ch));
}  /* uppercase */


void inputnumbers()
{
  /* input the numbers of species and of characters */
  long i;

  fscanf(infile, "%ld%ld", &spp, &chars);

  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n\n", spp, chars);
  node = (aas **)Malloc(spp * sizeof(aas *));
  if (firstset) {
    for (i = 0; i < spp; i++)
      node[i] = (aas *)Malloc(chars * sizeof(aas ));
  }
  weight = (steparray)Malloc(chars*sizeof(long));
  d      = (double **)Malloc(spp*sizeof(double *));
  nayme  = (char **)Malloc(spp*sizeof(char *));

  for (i=0;i<spp;++i){
    d[i] = (double *)Malloc(spp*sizeof(double));
    nayme[i] = (char *)malloc(nmlngth*sizeof(char));
  }

}  /* inputnumbers */

void getoptions()
{
  /* interactively set options */
  Char ch;
  Char in[100];
  boolean done, done1;

  if (printdata)
    fprintf(outfile, "\nProtein distance algorithm, version 3.5\n\n");
  putchar('\n');
  weights = false;
  printdata = false;
  progress = true;
  interleaved = true;
  ttratio = 2.0;
  whichcode = universal;
  whichcat = george;
  basesequal = true;
  freqa = 0.25;
  freqc = 0.25;
  freqg = 0.25;
  freqt = 0.25;
  usepam = true;
  kimura = false;
  ease = 0.457;
  do {
    printf(ansi ? "\033[2J\033[H" :
	   vt52 ? "\033E\033H" : "\n");
    printf("\nProtein distance algorithm, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  P  Use PAM, Kimura or categories model?  %s\n",
	   usepam ? "Dayhoff PAM matrix" :
	   kimura ? "Kimura formula" : "Categories model");
    if (!(usepam || kimura)) {
      printf("  C               Use which genetic code?  %s\n",
	     (whichcode == universal) ? "Universal"                  :
	     (whichcode == ciliate)   ? "Ciliate"                    :
	     (whichcode == mito)      ? "Universal mitochondrial"    :
             (whichcode == vertmito)  ? "Vertebrate mitochondrial"   :
	     (whichcode == flymito)   ? "Fly mitochondrial\n"        :
	     (whichcode == yeastmito) ? "Yeast mitochondrial"        : "");
      printf("  A  Which categorization of amino acids?  %s\n",
	     (whichcat == chemical) ? "Chemical"              :
	     (whichcat == george)   ? "George/Hunt/Barker"    : "Hall");
	
      printf("  E      Prob change category (1.0=easy):%8.4f\n",ease);
      printf("  T        Transition/transversion ratio:%7.3f\n",ttratio);
      printf("  F                     Base Frequencies:");
      if (basesequal)
	printf("  Equal\n");
      else
	printf("%7.3f%6.3f%6.3f%6.3f\n", freqa, freqc, freqg, freqt);
    }
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
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
	   progress ? "Yes" : "No");
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    in[0] = '\0';
//    gets(in);
    in[0] = 'y';
    ch=in[0];
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      if (ch == 'C' || ch == 'A' || ch == 'P' || ch == 'E' || ch == 'T' ||
	  ch == 'F' || ch == 'M' || ch == 'I' || ch == '1' || ch == '2' ||
	  ch == '0') {
	switch (ch) {

	case 'C':
	  printf("Which genetic code?\n");
	  printf(" type         for\n\n");
	  printf("   U           Universal\n");
	  printf("   M           Mitochondrial\n");
	  printf("   V           Vertebrate mitochondrial\n");
	  printf("   F           Fly mitochondrial\n");
	  printf("   Y           Yeast mitochondrial\n\n");
	  do {
	    printf("type U, M, V, F, or Y\n");
	    scanf("%c%*[^\n]", &ch);
	    getchar();
	    if (ch == '\n')
	      ch = ' ';
	    uppercase(&ch);
	  } while (ch != 'U' && ch != 'M' && ch != 'V' && ch != 'F' && ch != 'Y');
	  switch (ch) {

	  case 'U':
	    whichcode = universal;
	    break;

	  case 'M':
	    whichcode = mito;
	    break;

	  case 'V':
	    whichcode = vertmito;
	    break;

	  case 'F':
	    whichcode = flymito;
	    break;

	  case 'Y':
	    whichcode = yeastmito;
	    break;
	  }
	  break;

	case 'A':
	  printf(
	    "Which of these categorizations of amino acids do you want to use:\n\n");
	  printf(
	    " all have groups: (Glu Gln Asp Asn), (Lys Arg His), (Phe Tyr Trp)\n");
	  printf(" plus:\n");
	  printf("George/Hunt/Barker:");
	  printf(" (Cys), (Met   Val  Leu  Ileu), (Gly  Ala  Ser  Thr    Pro)\n");
	  printf("Chemical:          ");
	  printf(" (Cys   Met), (Val  Leu  Ileu    Gly  Ala  Ser  Thr), (Pro)\n");
	  printf("Hall:              ");
	  printf(" (Cys), (Met   Val  Leu  Ileu), (Gly  Ala  Ser  Thr), (Pro)\n\n");
	  printf("Which do you want to use (type C, H, or G)\n");
	  do {
	    scanf("%c%*[^\n]", &ch);
	    getchar();
	    if (ch == '\n')
	      ch = ' ';
	    uppercase(&ch);
	  } while (ch != 'C' && ch != 'H' && ch != 'G');
	  switch (ch) {

	  case 'C':
	    whichcat = chemical;
	    break;

	  case 'H':
	    whichcat = hall;
	    break;

	  case 'G':
	    whichcat = george;
	    break;
	  }
	  break;

	case 'P':
	  if (usepam) {
	    usepam = false;
	    kimura = true;
	  } else {
	    if (kimura)
	      kimura = false;
	    else
	      usepam = true;
	  }
	  break;

	case 'E':
	  printf("Ease of changing category of amino acid?\n");
	  do {
	    printf(" (1.0 if no difficulty of changing,\n");
	    printf(" less if less easy. Can't be negative\n");
	    scanf("%lf%*[^\n]", &ease);
	    getchar();
	  } while (ease > 1.0 || ease < 0.0);
	  break;

	case 'T':
	  do {
	    printf("Transition/transversion ratio?\n");
	    scanf("%lf%*[^\n]", &ttratio);
	    getchar();
	  } while (ttratio < 0.0);
	  break;

	case 'F':
	  do {
	    basesequal = false;
	    printf("Frequencies of bases A,C,G,T ?\n");
	    scanf("%lf%lf%lf%lf%*[^\n]", &freqa, &freqc, &freqg, &freqt);
	    getchar();
	    if (fabs(freqa + freqc + freqg + freqt - 1.0) >= 1.0e-3)
	      printf("FREQUENCIES MUST SUM TO 1\n");
	  } while (fabs(freqa + freqc + freqg + freqt - 1.0) >= 1.0e-3);
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
		printf("BAD NUMBER OF DATA SETS:  it must be greater than 1\n");
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
  } while (!done);
}  /* getoptions */

void transition()
{
  /* calculations related to transition-transversion ratio */
  double aa, bb, freqr, freqy, freqgr, freqty;

  freqr = freqa + freqg;
  freqy = freqc + freqt;
  freqgr = freqg / freqr;
  freqty = freqt / freqy;
  aa = ttratio * freqr * freqy - freqa * freqg - freqc * freqt;
  bb = freqa * freqgr + freqc * freqty;
  xi = aa / (aa + bb);
  xv = 1.0 - xi;
  if (xi <= 0.0 && xi >= -epsilon)
    xi = 0.0;
  if (xi < 0.0){
    printf("THIS TRANSITION-TRANSVERSION RATIO IS IMPOSSIBLE WITH");
    printf(" THESE BASE FREQUENCIES\n");
    exit(-1);}
}  /* transition */


void doinit()
{
  /* initializes variables */
  inputnumbers();
  getoptions();
  transition();
}  /* doinit*/



void inputweights()
{
  /* input the character weights, 0-9 and A-Z for weights 10 - 35 */
  Char ch;
  long i;

  for (i = 1; i < nmlngth; i++) {
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
  }
  for (i = 0; i < chars; i++) {
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
  /* print out the weights of sites */
  long i, j, k;

  fprintf(outfile, "    Sites are weighted as follows:\n");
  fprintf(outfile, "        ");
  for (i = 0; i <= 9; i++)
    fprintf(outfile, "%3ld", i);
  fprintf(outfile, "\n     *---------------------------------\n");
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

void inputoptions()
{
  /* input the information on the options */
  Char ch;
  long extranum, i, cursp, curchs;

  if (!firstset) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%ld%ld", &cursp, &curchs);
    if (cursp != spp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4ld\n",
	     ith+1);
      exit(-1);
    }
    chars = curchs;
  }
  extranum = 0;
  while (!(eoln(infile))) {
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
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
    if (ch == '\n')
      ch = ' ';
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
  long i, j, k, l, aasread, aasnew;
  Char charstate;
  boolean allread, done;
  aas aa;   /* temporary amino acid for input */

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
  aasread = 0;
  allread = false;
  while (!(allread)) {
    allread = true;
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    i = 1;
    while (i <= spp) {
      if ((interleaved && aasread == 0) || !interleaved) {
	for (j = 0; j < nmlngth; j++) {
	  nayme[i - 1][j] = getc(infile);
	  if (nayme[i - 1][j] == '\n')
	    nayme[i - 1][j] = ' ';
	  if (eof (infile) || eoln(infile)){
	    printf("ERROR: END-OF-LINE OR END-OF-FILE");
	    printf(" IN THE MIDDLE OF A SPECIES NAME\n");
	    exit(-1);}
	}
      }
      if (interleaved)
	j = aasread;
      else
	j = 0;
      done = false;
      while (((!done) && (!(eoln(infile) | eof(infile))))) {
	if (interleaved)
	  done = true;
	while (((j < chars) & (!(eoln(infile) | eof(infile))))) {
	  charstate = getc(infile);
	  if (charstate == '\n')
	    charstate = ' ';
          if (charstate == ' ' || (charstate >= '0' && charstate <= '9'))
            continue;
	  uppercase(&charstate);
	  if ((!isalpha(charstate) && charstate != '.' && charstate != '?' &&
	       charstate != '-' && charstate != '*') || charstate == 'J' ||
	      charstate == 'O' || charstate == 'U') {
	    printf("WARNING -- BAD AMINO ACID:%c AT POSITION%5ld OF SPECIES %3ld\n",
		   charstate, j+1, i);
	    exit(-1);
	  }
	  j++;
	  if (charstate == '.') {
	    node[i - 1][j - 1] = node[0][j - 1];
	    continue;
	  }
	  switch (charstate) {

	  case 'A':
	    aa = ala;
	    break;

	  case 'B':
	    aa = asx;
	    break;

	  case 'C':
	    aa = cys;
	    break;

	  case 'D':
	    aa = asp;
	    break;

	  case 'E':
	    aa = glu;
	    break;

	  case 'F':
	    aa = phe;
	    break;

	  case 'G':
	    aa = gly;
	    break;

	  case 'H':
	    aa = his;
	    break;

	  case 'I':
	    aa = ileu;
	    break;

	  case 'K':
	    aa = lys;
	    break;

	  case 'L':
	    aa = leu;
	    break;

	  case 'M':
	    aa = met;
	    break;

	  case 'N':
	    aa = asn;
	    break;

	  case 'P':
	    aa = pro;
	    break;

	  case 'Q':
	    aa = gln;
	    break;

	  case 'R':
	    aa = arg;
	    break;

	  case 'S':
	    aa = ser;
	    break;

	  case 'T':
	    aa = thr;
	    break;

	  case 'V':
	    aa = val;
	    break;

	  case 'W':
	    aa = trp;
	    break;

	  case 'X':
	    aa = unk;
	    break;

	  case 'Y':
	    aa = tyr;
	    break;

	  case 'Z':
	    aa = glx;
	    break;

	  case '*':
	    aa = stop;
	    break;

	  case '?':
	    aa = quest;
	    break;

	  case '-':
	    aa = del;
	    break;
	  }
	  node[i - 1][j - 1] = aa;
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
	aasnew = j;
      fscanf(infile, "%*[^\n]");
      getc(infile);
      if ((interleaved && j != aasnew) || ((!interleaved) && j != chars)){
	printf("ERROR: SEQUENCES OUT OF ALIGNMENT\n");
	exit(-1);}
      i++;
    }
    if (interleaved) {
      aasread = aasnew;
      allread = (aasread == chars);
    } else
      allread = (i > spp);
  }
  if ( printdata) {
    for (i = 1; i <= ((chars - 1) / 60 + 1); i++) {
      for (j = 1; j <= spp; j++) {
	for (k = 0; k < nmlngth; k++)
	  putc(nayme[j - 1][k], outfile);
	fprintf(outfile, "   ");
	l = i * 60;
	if (l > chars)
	  l = chars;
	for (k = (i - 1) * 60 + 1; k <= l; k++) {
	  if (j > 1 && node[j - 1][k - 1] == node[0][k - 1])
	    charstate = '.';
	  else {
	    switch (node[j - 1][k - 1]) {

	    case ala:
	      charstate = 'A';
	      break;

	    case asx:
	      charstate = 'B';
	      break;

	    case cys:
	      charstate = 'C';
	      break;

	    case asp:
	      charstate = 'D';
	      break;

	    case glu:
	      charstate = 'E';
	      break;

	    case phe:
	      charstate = 'F';
	      break;

	    case gly:
	      charstate = 'G';
	      break;

	    case his:
	      charstate = 'H';
	      break;

	    case ileu:
	      charstate = 'I';
	      break;

	    case lys:
	      charstate = 'K';
	      break;

	    case leu:
	      charstate = 'L';
	      break;

	    case met:
	      charstate = 'M';
	      break;

	    case asn:
	      charstate = 'N';
	      break;

	    case pro:
	      charstate = 'P';
	      break;

	    case gln:
	      charstate = 'Q';
	      break;

	    case arg:
	      charstate = 'R';
	      break;

	    case ser:
	      charstate = 'S';
	      break;

	    case thr:
	      charstate = 'T';
	      break;

	    case val:
	      charstate = 'V';
	      break;

	    case trp:
	      charstate = 'W';
	      break;

	    case tyr:
	      charstate = 'Y';
	      break;

	    case glx:
	      charstate = 'Z';
	      break;

	    case del:
	      charstate = '-';
	      break;

	    case stop:
	      charstate = '*';
	      break;

	    case unk:
	      charstate = 'X';
	      break;

	    case quest:
	      charstate = '?';
	      break;
	    }
	  }
	  putc(charstate, outfile);
	  if (k % 10 == 0 && k % 60 != 0)
	    putc(' ', outfile);
	}
	putc('\n', outfile);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  if (printdata)
    putc('\n', outfile);
}  /* inputdata */


void doinput()
{
  /* reads the input data */
  inputoptions();
  inputdata();
}  /* doinput */


void code()
{
  /* make up table of the code 1 = u, 2 = c, 3 = a, 4 = g */
  long n;
  aas b;

  trans[0][0][0] = phe;
  trans[0][0][1] = phe;
  trans[0][0][2] = leu;
  trans[0][0][3] = leu;
  trans[0][1][0] = ser;
  trans[0][1][1] = ser;
  trans[0][1][2] = ser;
  trans[0][1][3] = ser;
  trans[0][2][0] = tyr;
  trans[0][2][1] = tyr;
  trans[0][2][2] = stop;
  trans[0][2][3] = stop;
  trans[0][3][0] = cys;
  trans[0][3][1] = cys;
  trans[0][3][2] = stop;
  trans[0][3][3] = trp;
  trans[1][0][0] = leu;
  trans[1][0][1] = leu;
  trans[1][0][2] = leu;
  trans[1][0][3] = leu;
  trans[1][1][0] = pro;
  trans[1][1][1] = pro;
  trans[1][1][2] = pro;
  trans[1][1][3] = pro;
  trans[1][2][0] = his;
  trans[1][2][1] = his;
  trans[1][2][2] = gln;
  trans[1][2][3] = gln;
  trans[1][3][0] = arg;
  trans[1][3][1] = arg;
  trans[1][3][2] = arg;
  trans[1][3][3] = arg;
  trans[2][0][0] = ileu;
  trans[2][0][1] = ileu;
  trans[2][0][2] = ileu;
  trans[2][0][3] = met;
  trans[2][1][0] = thr;
  trans[2][1][1] = thr;
  trans[2][1][2] = thr;
  trans[2][1][3] = thr;
  trans[2][2][0] = asn;
  trans[2][2][1] = asn;
  trans[2][2][2] = lys;
  trans[2][2][3] = lys;
  trans[2][3][0] = ser;
  trans[2][3][1] = ser;
  trans[2][3][2] = arg;
  trans[2][3][3] = arg;
  trans[3][0][0] = val;
  trans[3][0][1] = val;
  trans[3][0][2] = val;
  trans[3][0][3] = val;
  trans[3][1][0] = ala;
  trans[3][1][1] = ala;
  trans[3][1][2] = ala;
  trans[3][1][3] = ala;
  trans[3][2][0] = asp;
  trans[3][2][1] = asp;
  trans[3][2][2] = glu;
  trans[3][2][3] = glu;
  trans[3][3][0] = gly;
  trans[3][3][1] = gly;
  trans[3][3][2] = gly;
  trans[3][3][3] = gly;
  if (whichcode == mito)
    trans[0][3][2] = trp;
  if (whichcode == vertmito) {
    trans[0][3][2] = trp;
    trans[2][3][2] = stop;
    trans[2][3][3] = stop;
    trans[2][0][2] = met;
  }
  if (whichcode == flymito) {
    trans[0][3][2] = trp;
    trans[2][0][2] = met;
    trans[2][3][2] = ser;
  }
  if (whichcode == yeastmito) {
    trans[0][3][2] = trp;
    trans[1][0][2] = thr;
    trans[2][0][2] = met;
  }
  n = 0;
  for (b = ala; (long)b <= (long)val; b = (aas)((long)b + 1)) {
    n++;
    numaa[(long)b - (long)ala] = n;
  }
}  /* code */


Static Void cats()
{
  /* define categories of amino acids */
  aas b;

  /* fundamental subgroups */
  cat[(long)cys - (long)ala] = 1;
  cat[(long)met - (long)ala] = 2;
  cat[(long)val - (long)ala] = 3;
  cat[(long)leu - (long)ala] = 3;
  cat[(long)ileu - (long)ala] = 3;
  cat[(long)gly - (long)ala] = 4;
  cat[0] = 4;
  cat[(long)ser - (long)ala] = 4;
  cat[(long)thr - (long)ala] = 4;
  cat[(long)pro - (long)ala] = 5;
  cat[(long)phe - (long)ala] = 6;
  cat[(long)tyr - (long)ala] = 6;
  cat[(long)trp - (long)ala] = 6;
  cat[(long)glu - (long)ala] = 7;
  cat[(long)gln - (long)ala] = 7;
  cat[(long)asp - (long)ala] = 7;
  cat[(long)asn - (long)ala] = 7;
  cat[(long)lys - (long)ala] = 8;
  cat[(long)arg - (long)ala] = 8;
  cat[(long)his - (long)ala] = 8;
  if (whichcat == george) {
    /* George, Hunt and Barker: sulfhydryl, small hydrophobic, small hydrophilic,
                              aromatic, acid/acid-amide/hydrophilic, basic */
    for (b = ala; (long)b <= (long)val; b = (aas)((long)b + 1)) {
      if (cat[(long)b - (long)ala] == 3)
	cat[(long)b - (long)ala] = 2;
      if (cat[(long)b - (long)ala] == 5)
	cat[(long)b - (long)ala] = 4;
    }
  }
  if (whichcat == chemical) {
    /* Conn and Stumpf:  monoamino, aliphatic, heterocyclic, aromatic, dicarboxylic,
                       basic */
    for (b = ala; (long)b <= (long)val; b = (aas)((long)b + 1)) {
      if (cat[(long)b - (long)ala] == 2)
	cat[(long)b - (long)ala] = 1;
      if (cat[(long)b - (long)ala] == 4)
	cat[(long)b - (long)ala] = 3;
    }
  }
  /* Ben Hall's personal opinion */
  if (whichcat != hall)
    return;
  for (b = ala; (long)b <= (long)val; b = (aas)((long)b + 1)) {
    if (cat[(long)b - (long)ala] == 3)
      cat[(long)b - (long)ala] = 2;
  }
}  /* cats */


void maketrans()
{
  /* Make up transition probability matrix from code and category tables */
  long i, j, k, m, n, s, nb1, nb2;
  double x, sum;
  long sub[3], newsub[3];
  double f[4], g[4];
  aas b1, b2;
  double TEMP, TEMP1, TEMP2, TEMP3;

  for (i = 0; i <= 19; i++) {
    pie[i] = 0.0;
    for (j = 0; j <= 19; j++)
      prob[i][j] = 0.0;
  }
  f[0] = freqt;
  f[1] = freqc;
  f[2] = freqa;
  f[3] = freqg;
  g[0] = freqc + freqt;
  g[1] = freqc + freqt;
  g[2] = freqa + freqg;
  g[3] = freqa + freqg;
  TEMP = f[0];
  TEMP1 = f[1];
  TEMP2 = f[2];
  TEMP3 = f[3];
  fracchange = xi * (2 * f[0] * f[1] / g[0] + 2 * f[2] * f[3] / g[2]) +
      xv * (1 - TEMP * TEMP - TEMP1 * TEMP1 - TEMP2 * TEMP2 - TEMP3 * TEMP3);
  sum = 0.0;
  for (i = 0; i <= 3; i++) {
    for (j = 0; j <= 3; j++) {
      for (k = 0; k <= 3; k++) {
	if (trans[i][j][k] != stop)
	  sum += f[i] * f[j] * f[k];
      }
    }
  }
  for (i = 0; i <= 3; i++) {
    sub[0] = i + 1;
    for (j = 0; j <= 3; j++) {
      sub[1] = j + 1;
      for (k = 0; k <= 3; k++) {
	sub[2] = k + 1;
	b1 = trans[i][j][k];
	for (m = 0; m <= 2; m++) {
	  s = sub[m];
	  for (n = 1; n <= 4; n++) {
	    memcpy(newsub, sub, sizeof(long) * 3L);
	    newsub[m] = n;
	    x = f[i] * f[j] * f[k] / (3.0 * sum);
	    if (((s == 1 || s == 2) && (n == 3 || n == 4)) ||
		((n == 1 || n == 2) && (s == 3 || s == 4)))
	      x *= xv * f[n - 1];
	    else
	      x *= xi * f[n - 1] / g[n - 1] + xv * f[n - 1];
	    b2 = trans[newsub[0] - 1][newsub[1] - 1][newsub[2] - 1];
	    if (b1 != stop) {
	      nb1 = numaa[(long)b1 - (long)ala];
	      pie[nb1 - 1] += x;
	      if (b2 != stop) {
		nb2 = numaa[(long)b2 - (long)ala];
		if (cat[(long)b1 - (long)ala] != cat[(long)b2 - (long)ala]) {
		  prob[nb1 - 1][nb2 - 1] += x * ease;
		  prob[nb1 - 1][nb1 - 1] += x * (1.0 - ease);
		} else
		  prob[nb1 - 1][nb2 - 1] += x;
	      } else
		prob[nb1 - 1][nb1 - 1] += x;
	    }
	  }
	}
      }
    }
  }
  for (i = 0; i <= 19; i++)
    prob[i][i] -= pie[i];
  for (i = 0; i <= 19; i++) {
    for (j = 0; j <= 19; j++)
      prob[i][j] /= sqrt(pie[i] * pie[j]);
  }
  /* computes pi^(1/2)*B*pi^(-1/2)  */
}  /* maketrans */



void givens(a, i, j, n, ctheta, stheta, left)
double (*a)[20];
long i, j, n;
double ctheta, stheta;
boolean left;
{
  /* Givens transform at i,j for 1..n with angle theta */
  long k;
  double d;

  for (k = 0; k < n; k++) {
    if (left) {
      d = ctheta * a[i - 1][k] + stheta * a[j - 1][k];
      a[j - 1][k] = ctheta * a[j - 1][k] - stheta * a[i - 1][k];
      a[i - 1][k] = d;
    } else {
      d = ctheta * a[k][i - 1] + stheta * a[k][j - 1];
      a[k][j - 1] = ctheta * a[k][j - 1] - stheta * a[k][i - 1];
      a[k][i - 1] = d;
    }
  }
}  /* givens */

void coeffs(x, y, c, s,accuracy)
double x, y, *c, *s;
double accuracy;
{
  /* compute cosine and sine of theta */
  double root;

  root = sqrt(x * x + y * y);
  if (root < accuracy) {
    *c = 1.0;
    *s = 0.0;
  } else {
    *c = x / root;
    *s = y / root;
  }
}  /* coeffs */

void tridiag(a, n,accuracy)
double (*a)[20];
long n;
double accuracy;
{
  /* Givens tridiagonalization */
  long i, j;
  double s, c;

  for (i = 2; i < n; i++) {
    for (j = i + 1; j <= n; j++) {
      coeffs(a[i - 2][i - 1], a[i - 2][j - 1], &c, &s,accuracy);
      givens(a, i, j, n, c, s, true);
      givens(a, i, j, n, c, s, false);
      givens(eigvecs, i, j, n, c, s, true);
    }
  }
}  /* tridiag */

void shiftqr(a, n, accuracy)
double (*a)[20];
long n;
double accuracy;
{
  /* QR eigenvalue-finder */
  long i, j;
  double approx, s, c, d, TEMP, TEMP1;

  for (i = n; i >= 2; i--) {
    do {
      TEMP = a[i - 2][i - 2] - a[i - 1][i - 1];
      TEMP1 = a[i - 1][i - 2];
      d = sqrt(TEMP * TEMP + TEMP1 * TEMP1);
      approx = a[i - 2][i - 2] + a[i - 1][i - 1];
      if (a[i - 1][i - 1] < a[i - 2][i - 2])
	approx = (approx - d) / 2.0;
      else
	approx = (approx + d) / 2.0;
      for (j = 0; j < i; j++)
	a[j][j] -= approx;
      for (j = 1; j < i; j++) {
	coeffs(a[j - 1][j - 1], a[j][j - 1], &c, &s, accuracy);
	givens(a, j, j + 1, i, c, s, true);
	givens(a, j, j + 1, i, c, s, false);
	givens(eigvecs, j, j + 1, n, c, s, true);
      }
      for (j = 0; j < i; j++)
	a[j][j] += approx;
    } while (fabs(a[i - 1][i - 2]) > accuracy);
  }
}  /* shiftqr */


void qreigen(prob, n)
double (*prob)[20];
long n;
{
  /* QR eigenvector/eigenvalue method for symmetric matrix */
  double accuracy;
  long i, j;

  accuracy = 1.0e-6;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      eigvecs[i][j] = 0.0;
    eigvecs[i][i] = 1.0;
  }
  tridiag(prob, n, accuracy);
  shiftqr(prob, n, accuracy);
  for (i = 0; i < n; i++)
    eig[i] = prob[i][i];
  for (i = 0; i <= 19; i++) {
    for (j = 0; j <= 19; j++)
      prob[i][j] = sqrt(pie[j]) * eigvecs[i][j];
  }
  /* prob[i][j] is the value of U' times pi^(1/2) */
}  /* qreigen */


void pameigen()
{
  /* eigenanalysis for PAM matrix, precomputed */
  memcpy(prob,pamprobs,sizeof(pamprobs));
  memcpy(eig,pameigs,sizeof(pameigs));
  fracchange = 0.01;
}  /* pameigen */


void predict(nb1, nb2)
long nb1, nb2;
{
  /* make contribution to prediction of this aa pair */
  long m;
  double TEMP;

  for (m = 0; m <= 19; m++) {
    elambdat = exp(tt * eig[m]);
    q = prob[m][nb1 - 1] * prob[m][nb2 - 1] * elambdat;
    p += q;
    dp += eig[m] * q;
    TEMP = eig[m];
    d2p += TEMP * TEMP * q;
  }
}  /* predict */


void makedists()
{
  /* compute the distances */
  long i, j, k, m, n, iterations, nb1, nb2;
  double delta, lnlike, slope, curv;
  boolean neginfinity, inf;
  aas b1, b2;

  fprintf(outfile, "%5ld\n", spp);
  if (progress)
    printf("Computing distances:\n");
  for (i = 1; i <= spp; i++) {
    if (progress)
      printf("  ");
    if (progress) {
      for (j = 0; j < nmlngth; j++)
	putchar(nayme[i - 1][j]);
    }
    if (progress)
      printf("   ");
    d[i - 1][i - 1] = 0.0;
    for (j = 0; j <= i - 2; j++) {
      if (!kimura) {
	if (usepam)
	  tt = 10.0;
	else
	  tt = 1.0;
	delta = tt / 2.0;
	iterations = 0;
        inf = false;
	do {
	  lnlike = 0.0;
	  slope = 0.0;
	  curv = 0.0;
	  neginfinity = false;
	  for (k = 0; k < chars; k++) {
	    b1 = node[i - 1][k];
	    b2 = node[j][k];
	    if (b1 != stop && b1 != del && b1 != quest && b1 != unk &&
		b2 != stop && b2 != del && b2 != quest && b2 != unk) {
	      p = 0.0;
	      dp = 0.0;
	      d2p = 0.0;
	      nb1 = numaa[(long)b1 - (long)ala];
	      nb2 = numaa[(long)b2 - (long)ala];
	      if (b1 != asx && b1 != glx && b2 != asx && b2 != glx)
		predict(nb1, nb2);
	      else {
		if (b1 == asx) {
		  if (b2 == asx) {
		    predict(3L, 3L);
		    predict(3L, 4L);
		    predict(4L, 3L);
		    predict(4L, 4L);
		  } else {
		    if (b2 == glx) {
		      predict(3L, 6L);
		      predict(3L, 7L);
		      predict(4L, 6L);
		      predict(4L, 7L);
		    } else {
		      predict(3L, nb2);
		      predict(4L, nb2);
		    }
		  }
		} else {
		  if (b1 == glx) {
		    if (b2 == asx) {
		      predict(6L, 3L);
		      predict(6L, 4L);
		      predict(7L, 3L);
		      predict(7L, 4L);
		    } else {
		      if (b2 == glx) {
			predict(6L, 6L);
			predict(6L, 7L);
			predict(7L, 6L);
			predict(7L, 7L);
		      } else {
			predict(6L, nb2);
			predict(7L, nb2);
		      }
		    }
		  } else {
		    if (b2 == asx) {
		      predict(nb1, 3L);
		      predict(nb1, 4L);
		      predict(nb1, 3L);
		      predict(nb1, 4L);
		    } else if (b2 == glx) {
		      predict(nb1, 6L);
		      predict(nb1, 7L);
		      predict(nb1, 6L);
		      predict(nb1, 7L);
		    }
		  }
		}
	      }
	      if (p <= 0.0)
		neginfinity = true;
	      else {
		lnlike += log(p);
		slope += dp / p;
		curv += d2p / p - dp * dp / (p * p);
	      }
	    }
	  }
	  iterations++;
	  if (!neginfinity) {
	    if (curv < 0.0) {
	      tt -= slope / curv;
              if (tt > 10000.0) {
                printf("\nWARNING: INFINITE DISTANCE BETWEEN SPECIES %ld AND %ld; -1.0 WAS WRITTEN\n", i, j);
                tt = -1.0/fracchange;
                inf = true;
                iterations = 20;
              }
            }
	    else {
	      if ((slope > 0.0 && delta < 0.0) || (slope < 0.0 && delta > 0.0))
		delta /= -2;
	      tt += delta;
	    }
	  } else {
	    delta /= -2;
	    tt += delta;
	  }
	  if (tt < epsilon && !inf)
	    tt = epsilon;
	} while (iterations != 20);
      } else {
	m = 0;
	n = 0;
	for (k = 0; k < chars; k++) {
	  b1 = node[i - 1][k];
	  b2 = node[j][k];
	  if ((long)b1 <= (long)val && (long)b2 <= (long)val) {
	    if (b1 == b2)
	      m++;
	    n++;
	  }
	}
	p = 1 - (double)m / n;
	dp = 1.0 - p - 0.2 * p * p;
	if (dp < 0.0) {
	  printf(
	    "\nDISTANCE BETWEEN SEQUENCES %3ld AND %3ld IS TOO LARGE FOR KIMURA FORMULA\n",
	    i, j + 1);
	  tt = -1.0;
	} else
	  tt = -log(dp);
      }
      d[i - 1][j] = fracchange * tt;
      d[j][i - 1] = d[i - 1][j];
      if (progress)
	putchar('.');
    }
    if (progress)
      putchar('\n');
  }
  for (i = 0; i < spp; i++) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[i][j], outfile);
    fprintf(outfile, "   ");
    for (j = 0; j < spp; j++)
      fprintf(outfile, "%9.5f", d[i][j]);
    putc('\n', outfile);
  }
  if (progress)
    printf("\nOutput written to output file\n\n");
}  /* makedists */


main(argc, argv)
int argc;
Char *argv[];
{  /* ML Protein distances by PAM or categories model */
  char infilename[100],outfilename[100];
#ifdef MAC
   macsetup("Protdist","");
   argv[0] = "Protdist";
#endif

//   openfile(&infile,INFILE,"r",argv[0],infilename);

// Aamer Modified Input File Line To Take From Argument
  if( argc < 2 ) {
    fprintf( stderr, "Usage: %s <inputfile name>\n", argv[0]);
    exit(0);
  }
  openfile(&infile,argv[1],"r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"w",argv[0],outfilename);

  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  mulsets = false;
  datasets = 1;
  firstset = true;
  doinit();
  if (!kimura)
    code();
  if (!(usepam || kimura)) {
    cats();
    maketrans();
    qreigen(prob, 20L);
  } else {
    if (kimura)
      fracchange = 1.0;
    else
      pameigen();
  }
  for (ith = 1; ith <= datasets; ith++) {
    doinput();
    if (ith == 1)
      firstset = false;
    if ((datasets > 1) && progress)
      printf("\nData set # %ld:\n\n", ith);
    makedists();
  }
  FClose(outfile);
  FClose(infile);
#ifdef MAC
  fixmacfile(outfilename);
#endif
  exit(0);
}  /* Protein distances */

int eof(f)
FILE *f;
{
    register int ch;

    if (feof(f))
        return 1;
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

