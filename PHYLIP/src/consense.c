#include "phylip.h"

/* version 3.56c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Hisashi Horino,
   Akiko Fuseki, Sean Lamont, and Andrew Keeffe.  Permission is granted
   to copy and use this program provided no fee is charged for it and
   provided that this copyright notice is not removed. */

#define nmlngth         20     /* max. no. of characters in species name  */

#define ibmpc0          false
#define ansi0           true
#define vt520           false

/* nodes will form a binary tree */

typedef Char naym[nmlngth];

typedef struct node {
  /* describes a tip species or an ancestor */
  struct node *next, *back;              /* pointers to nodes                */
  double index;                          /* number of the node               */
  boolean tip;                           /* present species are tips of tree */
  naym naime;
  long *nodeset;                         /* used by accumulate               */
  short xcoord, ycoord, ymin;            /* used by printree                 */
  short ymax;
} node;

typedef node **pointarray;

Static node *root;
Static double trweight, ntrees;
Static FILE *infile, *outfile, *treefile;
Static short spp, numopts, outgrno, i, j, col, setsz;
Static long maxgrp;                /* max. no. of groups in all trees found  */
Static boolean anerror, trout, first, noroot, outgropt, didreroot, prntsets,
	       progress, treeprint, ibmpc, vt52, ansi, goteof;
Static pointarray treenode;                 /* pointers to all nodes in tree */
Static naym *nayme;                 /* names of species       */
Static long **grouping, **grping2, **group2;/* to store groups found  */
Static long **order, **order2, lasti;
Static double **timesseen, **tmseen2, **times2; /* how often they're seen */
Static long *fullset;
Static node *garbage;

Static short tipy;

openfile(fp, filename, mode, application, perm)
FILE **fp;
char *filename;
char *mode;
char *application;
char *perm;
{
  FILE *of;
  char file[100];
  strcpy(file, filename);
  while (1){
    of = fopen(file, mode);
    if (of)
      break;
    else {
      switch (*mode){
      case 'r':
        printf("%s:  can't read %s\n", application, file);
        printf("Please enter a new filename>");
        gets(file);
        break;
      case 'w':
        printf("%s: can't write %s\n", application, file);
        printf("Please enter a new filename>");
        gets(file);
        break;
      }
    }
  }
  *fp=of;
  if (perm != NULL)
    strcpy(perm, file);
}


Static Void gnu(p)
node **p;
{
  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */
  if (garbage != NULL) {
    *p = garbage;
    garbage = garbage->next;
  } else
    *p = (node *)Malloc(sizeof(node));
  (*p)->next = NULL;
  (*p)->tip = false;
}  /* gnu */


Static Void chuck(p)
node *p;
{
  /* collect garbage on p -- put it on front of garbage list */
  p->next = garbage;
  garbage = p;
}  /* chuck */


Static Void gdispose(p)
node *p;
{
  /* go through tree throwing away nodes */
  node *q, *r;

  if (p->tip)
    return;
  q = p->next;
  while (q != p) {
    gdispose(q->back);
    r = q;
    q = q->next;
    chuck(r);
  }
  chuck(q);
}  /*  gdispose */


Static Void uppercase(ch)
Char *ch;
{  /* convert ch to upper case -- either ASCII or EBCDIC */
    *ch = (islower(*ch) ?  toupper(*ch) : (*ch));
}  /* uppercase */


Static Void getoptions()
{
  /* interactively set options */
  Char ch;
  boolean done, done1;

  fprintf(outfile, "\nMajority-rule and strict consensus tree");
  fprintf(outfile, " program, version %s\n\n", VERSION);
  putchar('\n');
  anerror = false;
  noroot = true;
  numopts = 0;
  outgrno = 1;
  outgropt = false;
  trout = true;
  prntsets = true;
  progress = true;
  treeprint = true;
  do {
    if (ansi)
      printf("\033[2J\033[H");
    else if (vt52)
      printf("\033E\033H");
    else
      putchar('\n');
    printf("\nMajority-rule and strict consensus tree");
    printf(" program, version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    if (noroot) {
      printf("  O                        Outgroup root?");
      if (outgropt)
        printf("  Yes, at species number%3hd\n", outgrno);
      else
        printf("  No, use as outgroup species%3hd\n", outgrno);
      }
    printf("  R        Trees to be treated as Rooted?");
    if (noroot)
      printf("  No\n");
    else
      printf("  Yes\n");
    printf("  0   Terminal type (IBM PC, VT52, ANSI)?");
    if (ibmpc)
      printf("  IBM PC\n");
    if (ansi)
      printf("  ANSI\n");
    if (vt52)
      printf("  VT52\n");
    if (!(ibmpc || vt52 || ansi))
      printf("  (none)\n");
    printf("  1         Print out the sets of species");
    if (prntsets)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  2  Print indications of progress of run  %s\n",
	   (progress ? "Yes" : "No"));
    printf("  3                        Print out tree");
    if (treeprint)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  4       Write out trees onto tree file?");
    if (trout)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      if ((noroot && (ch == 'O')) || ch == 'R' || ch == '0' || ch == '1' ||
          ch == '2' || ch == '3' || ch == '4') {
	switch (ch) {

	case 'O':
	  outgropt = !outgropt;
	  if (outgropt) {
	    numopts++;
	    done1 = true;
	    do {
	      printf("Type number of the outgroup:\n");
	      scanf("%hd%*[^\n]", &outgrno);
	      getchar();
	      done1 = (outgrno >= 1);
	      if (!done1) {
		printf("BAD OUTGROUP NUMBER: %4hd\n", outgrno);
		printf("  Must be greater than zero\n");
	      }
	    } while (done1 != true);
	  }
	  break;

	case 'R':
	  noroot = !noroot;
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
	  prntsets = !prntsets;
	  break;

        case '2':
	  progress = !progress;
          break;
	
	case '3':
	  treeprint = !treeprint;
	  break;

	case '4':
	  trout = !trout;
	  break;
	}
      } else
	printf("Not a possible option!\n");
    }
  } while (!done);
}  /* getoptions */

void getch(c, parens)
Char *c;
short *parens;
{
  /* get next nonblank character */
  do {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    *c = getc(infile);
  } while (!(*c != ' ' || eof(infile)));
  if (*c == '(')
   (*parens)++;
  if (*c == ')')
    (*parens)--;
}  /* getch */

Static Void setupnode(p, i)
node *p;
short i;
{
  /* initialization of node pointers, variables */

  p->next = NULL;
  p->back = NULL;
  p->index = i * 1.0;
  p->tip = false;
}  /* setupnode */

void dupname(name, p, found)
Char name[nmlngth];
node *p;
boolean *found;
{
  /* search for a duplicate name recursively */
  if (p) {
    if (p->next) {
      dupname(name, p->next->back, found);
      if (p->next->next)
        dupname(name, p->next->next->back, found);
    }
    if (p->tip) {
      for (i = 0; i < nmlngth; i++) {
        if (name[i] != p->naime[i])
          *found = true;
      }
    }
  }
}

void addelement(p, q, ch, parens, nextnode)
node **p, *q;
Char *ch;
short *nextnode, *parens;
{
  /* recursive procedure adds nodes to user-defined tree */
  node *r, *pfirst;
  short i, n, len;
  boolean found, notlast;
  Char str[nmlngth];

  if ((*ch) == '(') {
    (*nextnode)++;
    gnu(p);
    pfirst = *p;
    notlast = true;
    while (notlast && !anerror) {
      gnu(&(*p)->next);
      r = (*p)->next;
      getch(ch, parens);
      addelement(&r->back, r, ch, parens, nextnode);
      if (anerror)
	break;
      if (r->back != NULL)
	*p = r;
      if ((*ch) == ')') {
	notlast = false;
	do {
	  getch(ch, parens);
	} while ((*ch) != ',' && (*ch) != ')' && (*ch) != '[' && (*ch) != ';');
      }
    }
    if (!anerror) {
      (*p)->next = pfirst;
      *p = pfirst;
    }
  } else if ((*ch) != ')') {
    if (!anerror) {
      for (i = 0; i < nmlngth; i++)
	str[i] = ' ';
      n = 1;
      do {
	if ((*ch) == '_')
	  (*ch) = ' ';
	str[n - 1] = (*ch);
	if (eoln(infile)) {
	  fscanf(infile, "%*[^\n]");
	  getc(infile);
	}
	(*ch) = getc(infile);
	n++;
      } while ((*ch) != ':' && (*ch) != ',' && (*ch) != ')' &&
	       (*ch) != '[' && (*ch) != ';' && n <= nmlngth);
      if ((*ch) == ')')
	(*parens)--;
      len = n - 1;
      if (first && !anerror) {
	spp++;
        if (!anerror) {
          found = false;
	  dupname(str, root, &found);
	  if (found) {
            printf("\nERROR IN USER TREE: DUPLICATE NAME FOUND: ");
	    for (i = 0; i < nmlngth; i++)
	      putchar(str[i]);
            putchar('\n');
	    anerror = true;
          }
	}
	if (!anerror) {
	  gnu(p);
          setupnode(*p, spp);
	  (*p)->tip = true;
	  for (i = 0; i < nmlngth; i++)
	    (*p)->naime[i] = ' ';
	  for (i = 0; i < len; i++)
	    (*p)->naime[i] = str[i];
	  if (prntsets) {
	    fprintf(outfile, "  ");
	    for (i = 0; i < len; i++)
	      putc(str[i], outfile);
            putc('\n', outfile);
	  }
	}
      } else {
	n = 1;
	do {
	  found = true;
	  for (i = 0; i < nmlngth; i++)
	    found = (found && str[i] == nayme[n - 1][i]);
	  if (found)
	    *p = treenode[n - 1];
	  else
	    n++;
	} while (!(n > spp || found));
	if (n > spp) {
	  anerror = true;
	  printf("ERROR: CANNOT FIND SPECIES: ");
	  for (i = 0; i < nmlngth; i++)
	    putchar(str[i]);
	  putchar('\n');
	}
      }
    }
  } else
    getch(ch, parens);
  if ((*ch) == ':') {
    do {
      getch(ch, parens);
    } while ((*ch) != ',' && (*ch) != ')' && (*ch) != '[' && (*ch) != ';');
  }
  if ((*ch) == '[' && !anerror) {
    if (!eoln(infile)) {
      fscanf(infile, "%lf", &trweight);
      getch(ch, parens);
      if (*ch != ']') {
        printf("ERROR: MISSING RIGHT SQUARE BRACKET\n");
        anerror = true;
      }
      else {
        getch(ch, parens);
        if (*ch != ';') {
          printf("ERROR: MISSING SEMICOLON AFTER SQUARE BRACKETS\n");
          anerror = true;
        }
      }
    }
  }
  else if (*ch == ';' && !anerror) {
    trweight = 1.0;
    if (!eoln(infile))
      printf("WARNING: TREE WEIGHT SET TO 1.0\n");
  }
  if (*p != NULL && !anerror)
    (*p)->back = q;
}  /* addelement */


void initreenode(p)
node *p;
{
  /* traverse tree and assign tips to treenode */
  node *q;

  if (p) {
    q = p->next;
    while (q && q != p) {
      initreenode(q->back);
      q = q->next;
    }
    if (p->tip) {
      treenode[(int)p->index - 1] = p;
      memcpy(nayme[(int)p->index - 1], p->naime, nmlngth);
    }
  }
} /* initreenode */


Static Void treeread()
{
  /* read in user-defined tree and set it up */
  char ch;
  short parens, nextnode;

  goteof = false;
  parens = 0;
  getch(&ch, &parens);
  while (eoln(infile) && !eof(infile)) {
    fscanf(infile, "%*[^\n]");
    getc(infile);
  }
  if (eof(infile)) {
    goteof = true;
    return;
  } 
  addelement(&root, NULL, &ch, &parens, &nextnode);
    fscanf(infile, "%*[^\n]");
  getc(infile);
  if (parens != 0 && !anerror) {
    printf("\nERROR IN TREE FILE:  UNMATCHED PARENTHESES\n");
    anerror = true;
  }
  if (outgrno > spp && !anerror) {
    anerror = true;
    printf("ERROR IN OUTGROUP OPTION: SPECIES NUMBER%3hd DOES NOT EXIST\n",
	   outgrno);
  }
}  /* treeread */


Static Void reroot(outgroup)
node *outgroup;
{
  /* reorients tree, putting outgroup in desired position. */
  short i;
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
    gnu(&root->next);
    q = root->next;
    gnu(&q->next);
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

void rehash()
{
  long *s, k;
  short i, j;
  double temp, ss;
  boolean done;

  s = (long *)Malloc(setsz * sizeof(long));
  for (i = 0; i < maxgrp/2; i++) {
    k = *order[i];
    memcpy(s, grouping[k], setsz * sizeof(long));
    ss = 0.0;
    for (j = 0; j < setsz; j++)
      ss += s[j] /* pow(2, SETBITS*j)*/;
    temp = ss * ((sqrt(5.0) - 1) / 2);
    j = (long)(maxgrp * (temp - floor(temp)));
    done = false;
    while (!done) {
      if (!grping2[j]) {
        grping2[j] = (long *)Malloc(setsz * sizeof(long));
        order2[i] = (long *)Malloc(sizeof(long));
        tmseen2[j] = (double *)Malloc(sizeof(double));
        memcpy(grping2[j], grouping[k], setsz * sizeof(long));
        *tmseen2[j] = *timesseen[k];
        *order2[i] = j;
        grouping[k] = NULL;
        timesseen[k] = NULL;
        order[i] = NULL;
        done = true;
      } else {
        j++;
        if (j >= maxgrp) j -= maxgrp;
      }
    }
  }
  free(s);
}  /* rehash */


void enterset(r)
node *r;
{
  long *s, i, j, start;
  double times, ss, n;
  boolean done, same;

  s = (long *)Malloc(setsz * sizeof(long));
  memcpy(s, r->nodeset, setsz * sizeof(long));
  same = true;
  for (i = 0; i < setsz; i++)
    if (s[i] != fullset[i])
      same = false;
  if (same) return;
  times = trweight;
  ss = 0.0;
  n = ((sqrt(5.0) - 1) / 2);
  for (i = 0; i < setsz; i++)
    ss += s[i] * n;
  i = (long)(maxgrp * (ss - floor(ss))) + 1;
  start = i;
  done = false;
  while (!done) {
    if (grouping[i - 1]) {
      same = true;
      for (j = 0; j < setsz; j++) {
        if (s[j] != grouping[i - 1][j])
          same = false;
      }
    }
    if (grouping[i - 1] && same) {
      *timesseen[i - 1] += times;
      done = true;
    } else if (!grouping[i - 1]) {
      grouping[i - 1] = (long *)Malloc(setsz * sizeof(long));
      lasti++;
      order[lasti] = (long *)Malloc(sizeof(long));
      timesseen[i - 1] = (double *)Malloc(sizeof(double));
      memcpy(grouping[i - 1], s, setsz * sizeof(long));
      *timesseen[i - 1] = times;
      *order[lasti] = i - 1;
      done = true;
    } else {
      i++;
      if (i > maxgrp) i -= maxgrp;
    }
    if (!done && i == start) {
      maxgrp = maxgrp*2;
      tmseen2 = (double **)Malloc(maxgrp*sizeof(double *));
      grping2 = (long **)Malloc(maxgrp*sizeof(long *));
      order2 = (long **)Malloc(maxgrp*sizeof(long *));
      rehash();
      free(timesseen);
      free(grouping);
      free(order);
      timesseen = tmseen2;
      grouping = grping2;
      order = order2;
      done = true;
      lasti = maxgrp/2 - 1;
      enterset(r);
    }
  }
  free(s);
}  /* enterset */


Static Void accumulate(r_)
node *r_;
{
  node *r;
  node *q;
  short i;
  boolean done;

  r = r_;
  if (r->tip) {
    if (!r->nodeset)
      r->nodeset = (long *)Malloc(setsz * sizeof(long));
    for (i = 0; i < setsz; i++)
      r->nodeset[i] = 0L;
    done = false;
    i = 0;
    while (!done && i < setsz) {
      if ((int)r->index < (i+1)*SETBITS) {
        r->nodeset[i] = 1L << ((int)r->index - i*SETBITS);
        done = true;
      }
      else i++;
    }
  }
  else {
    q = r->next;
    while (q != r) {
      accumulate(q->back);
      q = q->next;
    }
    q = r->next;
    if (!r->nodeset)
      r->nodeset = (long *)Malloc(setsz * sizeof(long));
    for (i = 0; i < setsz; i++)
      r->nodeset[i] = 0;
    while (q != r) {
      for (i = 0; i < setsz; i++)
        r->nodeset[i] |= q->back->nodeset[i];
      q = q->next;
    }
  }
  if (!r->tip) {
    if (r->next->next != r)
      enterset(r);
  } else
    enterset(r);
}  /* accumulate */

#define down            2
#define over            5


void compress(n)
short *n;
{
  /* push all the nonempty subsets to the front end of their array */
  short i, j;

  i = 1;
  j = 1;
  do {
    while (grouping[i - 1])
      i++;
    if (j <= i)
      j = i + 1;
    while (!grouping[j - 1] && j < maxgrp)
      j++;
    if (j < maxgrp) {
      grouping[i - 1] = (long *)Malloc(setsz * sizeof(long));
      timesseen[i - 1] = (double *)Malloc(sizeof(double));
      memcpy(grouping[i - 1], grouping[j - 1], setsz * sizeof(long));
      *timesseen[i - 1] = *timesseen[j - 1];
      grouping[j - 1] = NULL;
      timesseen[j - 1] = NULL;
    }
  } while (j != maxgrp);
  (*n) = i - 1;
}  /* compress */


void sort(n)
short n;
{
  /* Shell sort keeping grouping, timesseen in same order */
  short gap, i, j;
  double rtemp;
  long *stemp;

  gap = n / 2;
  stemp = (long *)Malloc(setsz * sizeof(long));
  while (gap > 0) {
    for (i = gap + 1; i <= n; i++) {
      j = i - gap;
      while (j > 0) {
	if (*timesseen[j - 1] < *timesseen[j + gap - 1]) {
	  memcpy(stemp, grouping[j - 1], setsz * sizeof(long));
	  memcpy(grouping[j - 1], grouping[j + gap - 1], setsz * sizeof(long));
	  memcpy(grouping[j + gap - 1], stemp, setsz * sizeof(long));
	  rtemp = *timesseen[j - 1];
	  *timesseen[j - 1] = *timesseen[j + gap - 1];
	  *timesseen[j + gap - 1] = rtemp;
        }
	j -= gap;
      }
    }
    gap /= 2;
  }
  free(stemp);
}  /* sort */


void eliminate(n, n2, fullset)
short *n, *n2;
long *fullset;
{
  /* eliminate groups incompatible with preceding ones */
  short i, j, k;
  boolean comp;

  for (i = 2; i <= (*n); i++) {
    for (j = 0; j <= i - 2; j++) {
      if (timesseen[j] && *timesseen[j] > 0) {
        comp = true;
        for (k = 0; k < setsz; k++)
          if ((grouping[i - 1][k] & grouping[j][k]) != 0)
            comp = false;
        if (!comp) {
          comp = true;
          for (k = 0; k < setsz; k++)
            if ((grouping[i - 1][k] & ~grouping[j][k]) != 0)
              comp = false;
          if (!comp) {
            comp = true;
            for (k = 0; k < setsz; k++)
              if ((grouping[j][k] & ~grouping[i - 1][k]) != 0)
                comp = false;
            if (!comp) {
              comp = noroot;
              if (comp) {
                for (k = 0; k < setsz; k++)
                  if ((fullset[k] & ~grouping[i - 1][k] & ~grouping[j][k]) != 0)
                    comp = false;
              }
            }
          }
        }
        if (!comp) {
	  (*n2)++;
          times2[(*n2) - 1] = (double *)Malloc(sizeof(double));
          group2[(*n2) - 1] = (long *)Malloc(setsz * sizeof(long));
	  *times2[(*n2) - 1] = *timesseen[i - 1];
	  memcpy(group2[(*n2) - 1], grouping[i - 1], setsz * sizeof(long));
	  *timesseen[i - 1] = 0.0;
          for (k = 0; k < setsz; k++)
            grouping[i - 1][k] = 0;
	}
      }
    }
    if (*timesseen[i - 1] == 0.0) {
      timesseen[i - 1] = NULL;
      grouping[i - 1] = NULL;
    }
  }
}  /* eliminate */


void printset(n)
short n;
{
  /* print out a set of species */
  short i, j, k, size;

  fprintf(outfile, "\nSet (species in order)   ");
  for (i = 1; i <= spp - 25; i++)
    putc(' ', outfile);
  fprintf(outfile, "  How many times out of%7.2f\n\n", ntrees);
  for (i = 0; i < n; i++) {
    if (*timesseen[i] > 0) {
      size = 0;
      k = 0;
      for (j = 1; j <= spp; j++) {
        if (j == (k+1)*SETBITS) k++;
	if ((1L << (j - k*SETBITS) & grouping[i][k]) != 0)
	  size++;
      }
      if (size != 1 && !(noroot && size == spp - 1)) {
        k = 0;
	for (j = 1; j <= spp; j++) {
          if (j == (k+1)*SETBITS) k++;
	  if ((1L << (j - k*SETBITS) & grouping[i][k]) != 0)
	    putc('*', outfile);
	  else
	    putc('.', outfile);
	  if (j % 10 == 0)
	    putc(' ', outfile);
	}
	for (j = 1; j <= 23 - spp; j++)
	  putc(' ', outfile);
	fprintf(outfile, "    %5.2f\n", *timesseen[i]);
      }
    }
  }
}  /* printset */


void bigsubset(st, n)
long *st;
short n;
{
  /* find a maximal subset of st among the groupings */
  short i, j, k;
  long *su;
  boolean max, same;

  su = (long *)Malloc(setsz * sizeof(long));
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
      max = !same;
    }
    if (max) {
      for (j = 0; j < setsz; j ++)
        if ((su[j] & ~grouping[i][j]) != 0)	 
          max = false;
      if (max) {
        same = true;
        for (j = 0; j < setsz; j ++)
          if (su[j] != grouping[i][j])
            same = false;
        max = !same;
      }
      if (max)
        memcpy(su, grouping[i], setsz * sizeof(long));
    }
  }
  memcpy(st, su, setsz * sizeof(long));
  free(su);
}  /* bigsubset */


void recontraverse(p, st, n, nextnode)
node **p;
long *st;
short n;
short *nextnode;
{
  /* traverse to add next node of consensus tree */
  short i, j, k, l;
  boolean found, same, zero, zero2;
  long *tempset, *st2;
  node *q;

  k = l = 0;
  for (i = 1; i <= spp; i++) {
    if (i == (l+1)*SETBITS) l++;
    if (((1L << (i-l*SETBITS)) & st[l]) != 0) {
      k++;
      j = i;
    }
  }
  if (k == 1) {
    *p = treenode[j - 1];
    return;
  }
  gnu(p);
  (*nextnode)++;
  (*p)->index = 0.0;
  for (i = 0; i < n; i++) {
    same = true;
    for (j = 0; j < setsz; j++)
      if (grouping[i][j] != st[j])
        same = false;
    if (same)
      (*p)->index = *timesseen[i];
  }
  tempset = (long *)Malloc(setsz * sizeof(long));
  memcpy(tempset, st, setsz * sizeof(long));
  q = *p;
  st2 = (long *)Malloc(setsz * sizeof(long));
  memcpy(st2, st, setsz * sizeof(long));
  zero = true;
  for (j = 0; j < setsz; j++)
    if (tempset[j] != 0)
      zero = false;
  if (!zero)
    bigsubset(tempset, n);
  zero = zero2 = false;
  while (!zero && !zero2) {
    zero = zero2 = true;
    for (j = 0; j < setsz; j++) {
      if (st2[j] != 0)
        zero = false;
      if (tempset[j] != 0)
        zero2 = false;
    }
    if (!zero && !zero2) {
      gnu(&q->next);
      q = q->next;
      recontraverse(&q->back, tempset, n, nextnode);
      q->back->back = q;
      for (j = 0; j < setsz; j++)
        st2[j] &= ~tempset[j];
      memcpy(tempset, st2, setsz * sizeof(long));
      found = false;
      i = 1;
      while (!found && i <= n) {
        if (grouping[i - 1]) {
          same = true;
          for (j = 0; j < setsz; j++)
            if (grouping[i - 1][j] != tempset[j])
              same = false;
	}
        if (grouping[i - 1] && same)
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
    }
  }
  q->next = *p;
  free(tempset);
  free(st2);
}  /* recontraverse */


void reconstruct(n)
short n;
{
  /* reconstruct tree from the subsets */
  short nextnode, i;
  long *s;

  nextnode = spp + 1;
  s = (long *)Malloc(setsz * sizeof(long));
  memcpy(s, fullset, setsz * sizeof(long));
  recontraverse(&root, s, n, &nextnode);
  free(s);
}  /* reconstruct */


void coordinates(p, tipy)
node *p;
short *tipy;
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;
  short maxx;

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
    coordinates(q->back, tipy);
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
  p->xcoord = maxx + over;
  p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */


void drawline(i)
short i;
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q;
  short n, j;
  boolean extra, done, trif;
  node *r, *first, *last;
  boolean found;

  p = root;
  q = root;
  fprintf(outfile, "  ");
  extra = false;
  trif = false;
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
    if (extra) {
      n--;
      extra = false;
    }
    if (q->ycoord == i && !done) {
      if (trif)
	putc('-', outfile);
      else
	putc('+', outfile);
      trif = false;
      if (!q->tip) {
	for (j = 1; j <= n - 5; j++)
	  putc('-', outfile);
	if (q->index >= 100.0)
	  fprintf(outfile, "%5.1f", q->index);
	else if (q->index >= 10.0)
	  fprintf(outfile, "-%4.1f", q->index);
	else
	  fprintf(outfile, "--%3.1f", q->index);
	extra = true;
	trif = true;
      } else {
	for (j = 1; j < n; j++)
	  putc('-', outfile);
      }
    } else if (!p->tip && last->ycoord > i && first->ycoord < i &&
	       (i != p->ycoord || p == root)) {
      putc('!', outfile);
      for (j = 1; j < n; j++)
	putc(' ', outfile);
    } else {
      for (j = 1; j <= n; j++)
	putc(' ', outfile);
      if (trif)
	trif = false;
    }
    if (q != p)
      p = q;
  } while (!done);
  if (p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[(int)((long)p->index) - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */


void printree()
{
  /* prints out diagram of the tree */
  short i;
  short tipy;

  if (treeprint) {
    fprintf(outfile, "\nCONSENSUS TREE:\n");
    fprintf(outfile, "the numbers at the forks indicate the number\n");
    fprintf(outfile, "of times the group consisting of the species\n");
    fprintf(outfile, "which are to the right of that fork occurred\n");
    fprintf(outfile, "among the trees, out of %6.2f trees\n", ntrees);
    tipy = 1;
    coordinates(root, &tipy);
    putc('\n', outfile);
    for (i = 1; i <= tipy - down; i++)
      drawline(i);
    putc('\n', outfile);
  }
  if (noroot) {
    fprintf(outfile, "\n  remember:");
    if (didreroot)
      fprintf(outfile, " (though rerooted by outgroup)");
    fprintf(outfile, " this is an unrooted tree!\n");
  }
  putc('\n', outfile);
}  /* printree */


#undef down
#undef over

void treeout(p)
node *p;
{
  /* write out file with representation of final tree */
  short i, n;
  double x;
  Char c;
  node *q;

  if (p->tip) {
    for (i = 1; i <= nmlngth; i++) {
      if (nayme[(int)((long)p->index) - 1][i - 1] != ' ')
	n = i;
    }
    for (i = 0; i < n; i++) {
      c = nayme[(int)((long)p->index) - 1][i];
      if (c == ' ')
	c = '_';
      putc(c, treefile);
    }
    col += n;
  } else {
    putc('(', treefile);
    col++;
    q = p->next;
    while (q != p) {
      treeout(q->back);
      q = q->next;
      if (q == p)
	break;
      putc(',', treefile);
      col++;
      if (col > 60) {
	putc('\n', treefile);
	col = 0;
      }
    }
    putc(')', treefile);
    col++;
  }
  if (p->tip)
    x = ntrees;
  else
    x = p->index;
  if (p == root) {
    fprintf(treefile, ";\n");
    return;
  }
  if (x >= 100.0) {
    fprintf(treefile, ":%5.1f", x);
    col += 4;
    return;
  }
  if (x >= 10.0) {
    fprintf(treefile, ":%4.1f", x);
    col += 3;
  } else {
    fprintf(treefile, ":%3.1f", x);
    col += 2;
  }
}  /* treeout */


void consensus()
{
  short i;
  short n, n2;
short tipy;

  n2 = 0;
  compress(&n);
  sort(n);
  eliminate(&n, &n2, fullset);
  compress(&n);
  reconstruct(n);
    tipy = 1;
    coordinates(root, &tipy);
  if (prntsets) {
    fprintf(outfile, "\nSets included in the consensus tree\n");
    printset(n);
    for (i = 0; i < n2; i++) {
      if (!grouping[i]) {
        grouping[i] = (long *)Malloc(setsz * sizeof(long));
        timesseen[i] = (double *)Malloc(sizeof(double));
      }
      memcpy(grouping[i], group2[i], setsz * sizeof(long));
      *timesseen[i] = *times2[i];
    }
    n = n2;
    fprintf(outfile, "\n\nSets NOT included in consensus tree:");
    if (n2 == 0)
      fprintf(outfile, " NONE\n");
    else {
      putc('\n', outfile);
      printset(n);
    }
    putc('\n', outfile);
  }
  printree();
  if (progress)
    printf("\nOutput written to output file\n\n");
  if (trout) {
    treeout(root);
    if (progress)
      printf("Tree also written onto file\n\n");
  }
}  /* consensus */


main(argc, argv)
int argc;
Char *argv[];
{  
char infilename[100], outfilename[100], trfilename[100];
#ifdef MAC
  macsetup("Consense", "");
  argv[0] = "Consense";
#endif
  openfile(&infile, INFILE, "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "w", argv[0], outfilename);

  ibmpc = ibmpc0;
  ansi = ansi0;
  vt52 = vt520;
  didreroot = false;
  first = true;
  spp = 0;
  garbage = NULL;
  col = 0;
  getoptions();
  if (trout) 
    openfile(&treefile, TREEFILE, "w", argv[0], trfilename);
  if (prntsets)
    fprintf(outfile, "Species in order: \n\n");
  i = 1;
  ntrees = 0.0;
  maxgrp = 10000;
  lasti = -1;
  grouping = (long **)Malloc(maxgrp*sizeof(long *));
  order = (long **)Malloc(maxgrp*sizeof(long *));
  timesseen = (double **)Malloc(maxgrp*sizeof(double *));
  while (!anerror && !eof(infile)) {
    treeread();
    if (first) {
      treenode = (pointarray)Malloc(spp*sizeof(pointarray *));
      nayme = (naym *)Malloc(spp*sizeof(naym));
      initreenode(root);
      setsz = (short)ceil(((double)spp+1.0)/(double)SETBITS);
      fullset = (long *)Malloc(setsz * sizeof(long));
      for (j = 0; j < setsz; j++)
        fullset[j] = 0L;
      for (j = 0; j < setsz; j++) {
        if (spp + 1 < (j+1)*SETBITS) 
            fullset[j] = (1L << ((spp + 1) - j*SETBITS)) - 1;
        else
            fullset[j] = ~0L;
      }
      fullset[0] -= 1;
    }
    if (goteof)
      continue;
    ntrees += trweight;
    if (noroot) {
      reroot(treenode[outgrno - 1]);
      didreroot = outgropt;
    }
    accumulate(root);
    first = false;
    gdispose(root);
    i++;
  }
  putc('\n', outfile);
  group2 = (long **)Malloc(maxgrp*sizeof(long *));
  times2 = (double **)Malloc(maxgrp*sizeof(double *));
  consensus();

  FClose(treefile);
  FClose(infile);
  FClose(outfile);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(trfilename);
#endif
  exit(0);
}  /* main */


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
long  x;
{
MALLOCRETURN *mem;
mem = (MALLOCRETURN *)calloc(1, x);
if (!mem)
  memerror();
else
  return (MALLOCRETURN *)mem;
}

