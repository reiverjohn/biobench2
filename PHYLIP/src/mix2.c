#include "phylip.h"

/* version 3.572c. (c) Copyright 1993-1995 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define nmlngth         10   /* number of characters in species name    */
#define maxtrees        100  /* maximum number of tied trees stored     */
#define maxuser         10   /* maximum number of user-defined trees    */

#define ibmpc0          false
#define ansi0           true
#define vt520           false
#define down            2

typedef long *bitptr;
/* nodes will form a binary tree */

typedef struct node {           /* describes a tip species or an ancestor */
  struct node *next, *back;     /* pointers to nodes                      */
  long index;                   /* number of the node                     */
  boolean tip, bottom,visited;  /* present species are tips of tree       */
  bitptr fulstte1, fulstte0;  /* see in PROCEDURE fillin                */
  bitptr empstte1, empstte0;  /* see in PROCEDURE fillin                */
  bitptr fulsteps,empsteps;
  long xcoord, ycoord, ymin;    /* used by printree                       */
  long ymax;
} node;

typedef long *steptr;
typedef node **pointptr;
typedef long longer[6];
typedef long *placeptr;
typedef struct gbit {
  bitptr bits_;
  struct gbit *next;
} gbit;


/* Local variables for maketree: */
long nextree, which, minwhich;
double like, bestyet, bestlike, bstlike2, minsteps;
boolean lastrearr,full;
double nsteps[maxuser];
node *there;
long fullset;
bitptr steps, zeroanc, oneanc,fulzeroanc,empzeroanc;
long *place,col;

/* defined in main object module */
extern FILE     *infile,*outfile,*treefile;
extern node     *root;
extern steptr    numsteps,numsone,numszero,extras,weight;
extern long    **bestrees;
extern char     *guess;
extern double  **fsteps;
extern long     *enterorder;
extern longer    seed;
extern bitptr    wagner;
extern double   *threshwt;
extern char    **nayme;
extern boolean  *ancone,*anczero;
extern pointptr treenode;
extern boolean  jumble,usertree,questions,trout,noroot,outgropt,didreroot,
                progress,treeprint,stepbox,ancseq,weights;
extern char     ch;
extern long     chars,words,spp,nonodes,jumb,njumble,outgrno;
extern long      bits;

double randum();

void add(below, newtip, newfork)
node *below, *newtip, *newfork;
{
  /* inserts the nodes newfork and its left descendant, newtip,
    to the tree.  below becomes newfork's right descendant */
  node *p;

  if (below != treenode[below->index - 1])
    below = treenode[below->index - 1];
  if (below->back != NULL)
    below->back->back = newfork;
  newfork->back = below->back;
  below->back = newfork->next->next;
  newfork->next->next->back = below;
  newfork->next->back = newtip;
  newtip->back = newfork->next;
  if (root == below)
    root = newfork;
  root->back = NULL;
  p = newfork;
  while (p != NULL) {
    p->visited = false;
    p = p->back;
    if (p != NULL)
      p = treenode[p->index - 1];
  }
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
  if (root == *fork) {
    if (*item == (*fork)->next->back)
      root = (*fork)->next->next->back;
    else
      root = (*fork)->next->back;
  }
  p = (*item)->back->next->back;
  q = (*item)->back->next->next->back;
  if (p != NULL)
    p->back = q;
  if (q != NULL)
    q->back = p;
  q = (*fork)->back;
  (*fork)->back = NULL;
  p = (*fork)->next;
  while (p != *fork) {
    p->back = NULL;
    p = p->next;
  }
  (*item)->back = NULL;
  if (q != NULL)
    q = treenode[q->index - 1];
  while (q != NULL) {
    q->visited = false;
    q = q->back;
    if (q != NULL)
      q = treenode[q->index - 1];
  }
}  /* re_move */

void fillin(p)
node *p;
{
  /* Sets up for each node in the tree two statesets.
    state1 and state0 are the sets of character
    states that must be 1 or must be 0, respectively,
    in a most parsimonious reconstruction, based on the
    information at or above this node.  Note that this
    state assignment may change based on information further
    down the tree.  If a character is in both sets it is in
    state "P".  If in neither, it is "?". */
  long i;
  long l0, l1, r0, r1, st, wa, za;
  long FORLIM;

  for (i = 0; i < words; i++) {
    if (full) {
      l0 = p->next->back->fulstte0[i];
      l1 = p->next->back->fulstte1[i];
      r0 = p->next->next->back->fulstte0[i];
      r1 = p->next->next->back->fulstte1[i];
    }
    else {
      l0 = p->next->back->empstte0[i];
      l1 = p->next->back->empstte1[i];
      r0 = p->next->next->back->empstte0[i];
      r1 = p->next->next->back->empstte1[i];
    }
    wa = wagner[i];
    st = (l1 & r0) | (l0 & r1);
    if (full) {
      za = fulzeroanc[i];
      p->fulstte1[i] = (l1 | r1) & (~(st & (wa | za)));
      p->fulstte0[i] = (l0 | r0) & (~(st & (wa | (fullset & (~za)))));
      p->fulsteps[i] = st;
    }
    else {
      za = empzeroanc[i];
      p->empstte1[i] = (l1 | r1) & (~(st & (wa | za)));
      p->empstte0[i] = (l0 | r0) & (~(st & (wa | (fullset & (~za)))));
      p->empsteps[i] = st;
    }
  }
}  /* fillin */


void count(stps)
long *stps;
{
  /* counts the number of steps in a fork of the tree.
    The program spends much of its time in this PROCEDURE */
  long i, j, l;

  j = 1;
  l = 0;
  for (i = 0; i < (chars); i++) {
    l++;
    if (l > bits) {
      l = 1;
      j++;
    }
    if (((1L << l) & stps[j - 1]) != 0) {
      if (((1L << l) & zeroanc[j - 1]) != 0)
        numszero[i] += weight[i];
      else
        numsone[i] += weight[i];
    }
  }
}  /* count */

void postorder(p)
node *p;
{
  /* traverses a binary tree, calling PROCEDURE fillin at a
    node's descendants before calling fillin at the node */
  if (p->tip)
    return;
  postorder(p->next->back);
  postorder(p->next->next->back);
  if (p->visited)
    return;
  fillin(p);
  if (!full)
    p->visited = true;
}  /* postorder */


void cpostorder(p)
node *p;
{
  /* traverses a binary tree, calling PROCEDURE count at a
    node's descendants before calling count at the node */
  if (p->tip)
    return;
  cpostorder(p->next->back);
  cpostorder(p->next->next->back);
  if (full)
    count(p->fulsteps);
  else
    count(p->empsteps);
}  /* cpostorder */

void evaluate(r)
node *r;
{
  /* Determines the number of steps needed for a tree.
    This is the minimum number needed to evolve chars on
    this tree */
  long i, stepnum, smaller;
  double sum, term;

  sum = 0.0;
  for (i = 0; i < (chars); i++) {
    numszero[i] = 0;
    numsone[i] = 0;
  }
  full = true;
  for (i = 0; i < (words); i++)
    zeroanc[i] = fullset;
  postorder(r);
  cpostorder(r);
  count(r->fulstte1);
  for (i = 0; i < (words); i++)
    zeroanc[i] = 0;
  full = false;
  postorder(r);
  cpostorder(r);
  count(r->empstte0);
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
      term = stepnum;
    else
      term = threshwt[i];
    sum += term;
    if (usertree && which <= maxuser)
      fsteps[which - 1][i] = term;
    guess[i] = '?';
    if (!ancone[i] || (anczero[i] && numszero[i] < numsone[i]))
      guess[i] = '0';
    else if (!anczero[i] || (ancone[i] && numsone[i] < numszero[i]))
      guess[i] = '1';
  }
  if (usertree && which <= maxuser) {
    nsteps[which - 1] = sum;
    if (which == 1) {
      minwhich = 1;
      minsteps = sum;
    } else if (sum < minsteps) {
      minwhich = which;
      minsteps = sum;
    }
  }
  like = -sum;
}  /* evaluate */


void reroot(outgroup)
node *outgroup;
{
  /* reorients tree, putting outgroup in desired position. */
  node *p, *q;

  if (outgroup->back->index == root->index)
    return;
  p = root->next;
  q = root->next->next;
  p->back->back = q->back;
  q->back->back = p->back;
  p->back = outgroup;
  q->back = outgroup->back;
  outgroup->back->back = q;
  outgroup->back = p;
}  /* reroot */

void savetraverse(p)
node *p;
{
  /* sets BOOLEANs that indicate which way is down */
  p->bottom = true;
  if (p->tip)
    return;
  p->next->bottom = false;
  savetraverse(p->next->back);
  p->next->next->bottom = false;
  savetraverse(p->next->next->back);
}  /* savetraverse */

void savetree()
{
  /* record in place where each species has to be
    added to reconstruct this tree */
  long i, j;
  node *p;
  boolean done;

  if (noroot)
    reroot(treenode[outgrno - 1]);
  savetraverse(root);
  for (i = 0; i < (nonodes); i++)
   place[i] = 0;
  place[root->index - 1] = 1;
  for (i = 1; i <= (spp); i++) {
    p = treenode[i - 1];
    while (place[p->index - 1] == 0) {
      place[p->index - 1] = i;
      while (!p->bottom)
        p = p->next;
      p = p->back;
    }
    if (i > 1) {
      place[i - 1] = place[p->index - 1];
      j = place[p->index - 1];
      done = false;
      while (!done) {
        place[p->index - 1] = spp + i - 1;
        while (!p->bottom)
          p = p->next;
        p = p->back;
        done = (p == NULL);
        if (!done)
          done = (place[p->index - 1] != j);
      }
    }
  }
}  /* savetree */

void findtree(pos,found)
long *pos;
boolean *found;
{
  /* finds tree given by ARRAY place in ARRAY
    bestrees by binary search */
  long i, lower, upper;
  boolean below, done;

  below = false;
  lower = 1;
  upper = nextree - 1;
  *found = false;
  while (!(*found) && lower <= upper) {
    *pos = (lower + upper) / 2;
    i = 3;
    done = false;
    while (!done) {
      done = (i > spp);
      if (!done)
        done = (place[i - 1] != bestrees[*pos - 1][i - 1]);
      if (!done)
        i++;
    }
    *found = (i > spp);
    if (*found)
      break;
    else
      below = (place[i - 1] < bestrees[*pos - 1][i - 1]);
    if (below)
      upper = (*pos) - 1;
    else
      lower = (*pos) + 1;
  }
  if (!(*found) && !below)
    (*pos)++;
}  /* findtree */

void addtree(pos)
long *pos;
{
  /* puts tree from ARRAY place in its proper position
    in ARRAY bestrees */
  long i;
  for (i =nextree - 1; i >= (*pos); i--)
    memcpy(bestrees[i], bestrees[i - 1], spp*sizeof(long));
  for (i = 0; i < (spp); i++)
    bestrees[(*pos) - 1][i] = place[i];
  nextree++;
}  /* addtree */

void tryadd(p, item,nufork)
node *p,**item,**nufork;
{
  /* temporarily adds one fork and one tip to the tree.
    if the location where they are added yields greater
    "likelihood" than other locations tested up to that
    time, then keeps that location as there */

  long pos;
  boolean found;
  node *rute;

  add(p, *item, *nufork);
  evaluate(root);
  if (lastrearr) {
    if (like >= bstlike2) {
      rute = root->next->back;
      savetree();
      reroot(rute);
      if (like > bstlike2) {
        bestlike = bstlike2 = like;
        pos = 1;
        nextree = 1;
        addtree(&pos);
      } else {
        pos = 0;
        findtree(&pos,&found);
        if (!found) {
          if (nextree <= maxtrees)
            addtree(&pos);
        }
      }
    }
  }
  if (like > bestyet) {
    bestyet = like;
    there = p;
  }
  re_move(item, nufork);
}  /* tryadd */

void addpreorder(p, item, nufork)
node *p, *item, *nufork;
{
  /* traverses a binary tree, calling PROCEDURE tryadd
    at a node before calling tryadd at its descendants */
  /* Local variables for addpreorder: */

  if (p == NULL)
    return;
  tryadd(p, &item,&nufork);
  if (!p->tip) {
    addpreorder(p->next->back, item, nufork);
    addpreorder(p->next->next->back, item, nufork);
  }
}  /* addpreorder */

void tryrearr(p, r,success)
node *p,**r;
boolean *success;
{
  /* evaluates one rearrangement of the tree.
    if the new tree has greater "likelihood" than the old
    one sets success := TRUE and keeps the new tree.
    otherwise, restores the old tree */
  node *frombelow, *whereto, *forknode;
  double oldlike;

  if (p->back == NULL)
    return;
  forknode = treenode[p->back->index - 1];
  if (forknode->back == NULL)
    return;
  oldlike = bestyet;
  if (p->back->next->next == forknode)
    frombelow = forknode->next->next->back;
  else
    frombelow = forknode->next->back;
  whereto = treenode[forknode->back->index - 1];
  re_move(&p, &forknode);
  add(whereto, p, forknode);
  evaluate(*r);
  if (like <= oldlike) {
    re_move(&p, &forknode);
    add(frombelow, p, forknode);
  } else {
    *success = true;
    bestyet = like;
  }
}  /* tryrearr */

void repreorder(p, r,success)
node *p,**r;
boolean *success;
{
  /* traverses a binary tree, calling PROCEDURE tryrearr
    at a node before calling tryrearr at its descendants */
  if (p == NULL)
    return;
  tryrearr(p, r,success);
  if (!p->tip) {
    repreorder(p->next->back, r,success);
    repreorder(p->next->next->back, r,success);
  }
}  /* repreorder */

void rearrange(r)
node **r;
{
  /* traverses the tree (preorder), finding any local
    rearrangement which decreases the number of steps.
    if traversal succeeds in increasing the tree's
    "likelihood", PROCEDURE rearrange runs traversal again */
  boolean success=true;
  while (success) {
    success = false;
    repreorder(*r,r,&success);
  }
}  /* rearrange */


void findch(c)
Char c;
{
  /* scan forward until find character c */
  boolean done;

  done = false;
  while (!done) {
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
	printf("\nERROR IN USER TREE:");
	printf(" UNMATCHED PARENTHESIS OR MISSING SEMICOLON\n");
	exit(-1);
      } else
        done = true;
    }
    if (done && ch != ')')
      continue;
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
  }
}  /* findch */

void addelement(p, nextnode, lparens,naymes)
node **p;
long *nextnode,*lparens;
boolean *naymes;
{
  /* recursive procedure adds nodes to user-defined tree */
  node *q;
  long i, n;
  boolean found;
  Char str[nmlngth];

  do {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
  } while (ch == ' ');
  if (ch == '(' ) {
    if ((*lparens) >= spp - 1) {
      printf("\nERROR IN USER TREE: TOO MANY LEFT PARENTHESES\n");
      exit(-1);
    }
    (*nextnode)++;
    (*lparens)++;
    q = treenode[(*nextnode) - 1];
    addelement(&q->next->back, nextnode,lparens,naymes);
    q->next->back->back = q->next;
    findch(',');
    addelement(&q->next->next->back, nextnode,lparens,naymes);
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
    str[n - 1] =ch;
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
    n++;
  } while (ch != ',' && ch != ')' && ch != ':' && n <= nmlngth);
  n = 1;
  do {
    found = true;
    for (i = 0; i < nmlngth; i++)
      found = (found && str[i] == nayme[n - 1][i]);
    if (found) {
      if (naymes[n - 1] == false) {
        *p = treenode[n - 1];
        naymes[n - 1] = true;
      } else {
        printf("\nERROR IN USER TREE: DUPLICATE NAME FOUND -- ");
        for (i = 0; i < nmlngth; i++)
          putchar(nayme[n - 1][i]);
        putchar('\n');
	exit(-1);
      }
    } else
      n++;
  } while (!(n > spp || found ));
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
/* Local variables for treeread: */
  long nextnode, lparens;
  boolean *naymes;
  long i;

  root = treenode[spp];
  nextnode = spp;
  root->back = NULL;
  naymes = (boolean *)Malloc(spp*sizeof(boolean));
  for (i = 0; i < (spp); i++)
    naymes[i] = false;
  lparens = 0;
  addelement(&root, &nextnode,&lparens,naymes);
  for (i = 0; i <= nextnode-1; i++) {
    treenode[i]->visited = false;
    if (i >= spp) {
      treenode[i]->next->visited = false;
      treenode[i]->next->next->visited = false;
    } 
  }
  findch(';');
  if (progress)
    printf("\n\n");
  fscanf(infile, "%*[^\n]");
  getc(infile);
  free(naymes);
}  /* treeread */


void coordinates(p, tipy)
node *p;
long *tipy;
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
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
  p->xcoord = last->ymax - first->ymin;
  p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */

void drawline(i, scale)
long i;
double scale;
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q, *r, *first, *last;
  long n, j;
  boolean extra, done;

  p = root;
  q = root;
  extra = false;
  if (i == p->ycoord && p == root) {
    if (p->index - spp >= 10)
      fprintf(outfile, "-%2ld", p->index - spp);
    else
      fprintf(outfile, "--%ld", p->index - spp);
    extra = true;
  } else
    fprintf(outfile, "  ");
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
    n = (long)(scale * (p->xcoord - q->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if (q->ycoord == i && !done) {
      putc('+', outfile);
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        if (q->index - spp >= 10)
          fprintf(outfile, "%2ld", q->index - spp);
        else
          fprintf(outfile, "-%ld", q->index - spp);
        extra = true;
      } else {
        for (j = 1; j < n; j++)
          putc('-', outfile);
      }
    } else if (!p->tip) {
      if (last->ycoord > i && first->ycoord < i && i != p->ycoord) {
        putc('!', outfile);
        for (j = 1; j < n; j++)
          putc(' ', outfile);
      } else {
        for (j = 1; j <= n; j++)
          putc(' ', outfile);
      }
    } else {
      for (j = 1; j <= n; j++)
        putc(' ', outfile);
    }
    if (p != q)
      p = q;
  } while (!done);
  if (p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */

void printree()
{
  /* prints out diagram of the tree */
  long tipy;
  double scale;
  long i;

  putc('\n', outfile);
  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  coordinates(root, &tipy);
  scale = 1.5;
  putc('\n', outfile);
  for (i = 1; i <= (tipy - down); i++)
    drawline(i, scale);
  if (noroot) {
    fprintf(outfile, "\n  remember:");
    if (didreroot)
      fprintf(outfile, " (although rooted by outgroup)");
    fprintf(outfile, " this is an unrooted tree!\n");
  }
  putc('\n', outfile);
}  /* printree */


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


void hyprint(r,bottom,nonzero,unknown,maybe,zerobelow,onebelow)
node *r;
boolean bottom,nonzero,unknown,maybe;
gbit *zerobelow,*onebelow;
{
  /* print out states at node */
  long i, j, k;
  char l;
  boolean dot, a0, a1, s0, s1;

  if (bottom) {
    if (noroot && !didreroot)
      fprintf(outfile, "       ");
    else
      fprintf(outfile, "root   ");
  } else
    fprintf(outfile, "%3ld    ", r->back->index - spp);
  if (r->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[r->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", r->index - spp);
  if (bottom && noroot && !didreroot)
    fprintf(outfile, "          ");
  else if (nonzero)
    fprintf(outfile, "   yes    ");
  else if (unknown)
    fprintf(outfile, "    ?     ");
  else if (maybe)
    fprintf(outfile, "  maybe   ");
  else
    fprintf(outfile, "   no     ");
  for (j = 1; j <= (chars); j++) {
    newline(j, 40L, nmlngth + 17L);
    k = (j - 1) / bits + 1;
    l = (j - 1) % bits + 1;
    dot = (((1L << l) & wagner[k - 1]) == 0 && guess[j - 1] == '?');
    s0 = (((1L << l) & r->fulstte0[k - 1]) != 0);
    s1 = (((1L << l) & r->fulstte1[k - 1]) != 0);
    a0 = (((1L << l) & zerobelow->bits_[k - 1]) != 0);
    a1 = (((1L << l) & onebelow->bits_[k - 1]) != 0);
    dot = (dot || ((!bottom || !noroot || didreroot) && a1 == s1 &&
                   a0 == s0 && s0 != s1));
    if (dot)
      putc('.', outfile);
    else {
      if (s0)
        putc('0', outfile);
      else if (s1)
        putc('1', outfile);
      else
        putc('?', outfile);
    }
    if (j % 5 == 0)
      putc(' ', outfile);
  }
  putc('\n', outfile);
}  /* hyprint */

void hyptrav(r, dohyp,unknown)
node *r;
bitptr dohyp;
boolean unknown;
{
  /*  compute, print out states at one interior node */
  /* Local variables for hyptrav: */
  boolean bottom, maybe, nonzero;
  gbit *zerobelow, *onebelow;
  long i;
  long l0, l1, r0, r1, s0, s1, a0, a1, temp, dh, wa;

  gnu(&zerobelow);
  gnu(&onebelow);
  bottom = (r->back == NULL);
  maybe = false;
  nonzero = false;
  if (bottom) {
    memcpy(zerobelow->bits_, zeroanc, words*sizeof(long));
    memcpy(onebelow->bits_, oneanc, words*sizeof(long));
  } else {
    memcpy(zerobelow->bits_, treenode[r->back->index - 1]->fulstte0,
           words*sizeof(long));
    memcpy(onebelow->bits_, treenode[r->back->index - 1]->fulstte1,
           words*sizeof(long));
  }
  for (i = 0; i < (words); i++) {
    dh = dohyp[i];
    s0 = r->fulstte0[i];
    s1 = r->fulstte1[i];
    a0 = zerobelow->bits_[i];
    a1 = onebelow->bits_[i];
    if (!r->tip) {
      wa = wagner[i];
      l0 = r->next->back->fulstte0[i];
      l1 = r->next->back->fulstte1[i];
      r0 = r->next->next->back->fulstte0[i];
      r1 = r->next->next->back->fulstte1[i];
      s0 = (wa & ((a0 & l0) | (a0 & r0) | (l0 & r0))) |
           (dh & fullset & (~wa) & s0);
      s1 = (wa & ((a1 & l1) | (a1 & r1) | (l1 & r1))) |
           (dh & fullset & (~wa) & s1);
      temp = fullset & (~(s0 | s1 | l1 | l0 | r1 | r0));
      s0 |= temp & a0;
      s1 |= temp & a1;
      r->fulstte0[i] = s0;
      r->fulstte1[i] = s1;
    }
    maybe = (maybe || (dh & (s0 | s1)) != (a0 | a1));
    nonzero = (nonzero || ((s1 & a0) | (s0 & a1)) != 0);
  }
  hyprint(r,bottom,nonzero,unknown,maybe,zerobelow,onebelow);
  if (!r->tip) {
    hyptrav(r->next->back, dohyp,unknown);
    hyptrav(r->next->next->back, dohyp,unknown);
  }
  chuck(zerobelow);
  chuck(onebelow);
}  /* hyptrav */

void hypstates()
{
  /* fill in and describe states at interior nodes */
/* Local variables for hypstates: */
  boolean unknown;
  bitptr dohyp;
  long i, j, k;


  for (i = 0; i < (words); i++) {
    zeroanc[i] = 0;
    oneanc[i] = 0;
  }
  unknown = false;
  for (i = 0; i < (chars); i++) {
    j = i / bits + 1;
    k = i % bits + 1;
    if (guess[i] == '0')
      zeroanc[j - 1] = ((long)zeroanc[j - 1]) |
                                   (1L << ((int)k));
    if (guess[i] == '1')
      oneanc[j - 1] = ((long)oneanc[j - 1]) | (1L << ((int)k));
    unknown = (unknown ||
        ((((1L << k) & wagner[j - 1]) == 0) && guess[i] == '?'));
  }
  dohyp = (bitptr)Malloc(words*sizeof(long));
  for (i = 0; i < (words); i++)
    dohyp[i] = wagner[i] | zeroanc[i] | oneanc[i];
  filltrav(root);
  fprintf(outfile, "From    To     Any Steps?");
  fprintf(outfile, "    State at upper node\n");
  fprintf(outfile, "                            ");
  fprintf(outfile, " ( . means same as in the node below it on tree)\n\n");
  hyptrav(root, dohyp,unknown);
  free(dohyp);
}  /* hypstates */

void treeout(p)
node *p;
{
  /* write out file with representation of final tree */
  long i, n;
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
  if (p != root)
    return;
  if (nextree > 2)
    fprintf(treefile, "[%6.4f];\n", 1.0 / (nextree - 1));
  else
    fprintf(treefile, ";\n");
}  /* treeout */

void describe()
{
  /* prints ancestors, steps and table of numbers of steps in
    each character */
  long i, j, k;

  if (treeprint)
    fprintf(outfile, "\nrequires a total of %10.3f\n", -like);
  if (stepbox) {
    putc('\n', outfile);
    if (weights)
      fprintf(outfile, " weighted");
    fprintf(outfile, " steps in each character:\n");
    fprintf(outfile, "      ");
    for (i = 0; i <= 9; i++)
      fprintf(outfile, "%4ld", i);
    fprintf(outfile, "\n     *-----------------------------------------\n");
    for (i = 0; i <= (chars / 10); i++) {
      fprintf(outfile, "%5ld", i * 10);
      putc('!', outfile);
      for (j = 0; j <= 9; j++) {
        k = i * 10 + j;
        if (k == 0 || k > chars)
          fprintf(outfile, "    ");
        else
          fprintf(outfile, "%4ld", numsteps[k - 1] + extras[k - 1]);
      }
      putc('\n', outfile);
    }
  }
  putc('\n', outfile);
  if (questions && (!noroot || didreroot)) {
    fprintf(outfile, "best guesses of ancestral states:\n");
    fprintf(outfile, "      ");
    for (i = 0; i <= 9; i++)
      fprintf(outfile, "%2ld", i);
    fprintf(outfile, "\n     *--------------------\n");
    for (i = 0; i <= (chars / 10); i++) {
      putc(' ', outfile);
      fprintf(outfile, "%4ld!", i * 10);
      for (j = 0; j <= 9; j++) {
        if (i * 10 + j == 0 || i * 10 + j > chars)
          fprintf(outfile, "  ");
        else
          fprintf(outfile, " %c", guess[i * 10 + j - 1]);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  if (ancseq) {
    hypstates();
    putc('\n', outfile);
  }
  putc('\n', outfile);
  if (trout) {
    col = 0;
    treeout(root);
  }
}  /* describe */


void maketree()
{
  /* constructs a binary tree from the pointers in treenode.
    adds each node at location which yields highest "likelihood"
    then rearranges the tree for greatest "likelihood" */
  long i, j, k, numtrees, num, sumw;
  double gotlike, sum, sum2, sd;
  node *item, *nufork, *dummy;
  double TEMP;
  place = (long *)Malloc(nonodes*sizeof(long));
  steps = (bitptr)Malloc(words*sizeof(long));
  zeroanc = (bitptr)Malloc(words*sizeof(long));
  oneanc = (bitptr)Malloc(words*sizeof(long));
  fulzeroanc = (bitptr)Malloc(words*sizeof(long));
  empzeroanc = (bitptr)Malloc(words*sizeof(long));

  fullset = (1L << (bits + 1)) - (1L << 1);
  for (i=0 ; i<words ; ++i){
    fulzeroanc[i]=fullset;
    empzeroanc[i]= 0;}
  if (!usertree) {
    for (i = 1; i <= (spp); i++)
      enterorder[i - 1] = i;
    if (jumble) {
      for (i = 0; i < (spp); i++) {
        j = (long)(randum(seed) * spp) + 1;
        k = enterorder[j - 1];
        enterorder[j - 1] = enterorder[i];
        enterorder[i] = k;
      }
    }
    root = treenode[enterorder[0] - 1];
    add(treenode[enterorder[0] - 1], treenode[enterorder[1] - 1],
        treenode[spp]);
    if (progress) {
      printf("\nAdding species:\n");
      printf("   ");
      for (i = 0; i < nmlngth; i++)
        putchar(nayme[enterorder[0] - 1][i]);
      printf("\n   ");
      for (i = 0; i < nmlngth; i++)
        putchar(nayme[enterorder[1] - 1][i]);
      putchar('\n');
    }
    lastrearr = false;
    for (i = 3; i <= (spp); i++) {
      bestyet = -10.0*spp*chars;
      item = treenode[enterorder[i - 1] - 1];
      nufork = treenode[spp + i - 2];
      addpreorder(root, item, nufork);
      add(there, item, nufork);
      like = bestyet;
      rearrange(&root);
      if (progress) {
        printf("   ");
        for (j = 0; j < nmlngth; j++)
          putchar(nayme[enterorder[i - 1] - 1][j]);
        putchar('\n');
      }
      lastrearr = (i == spp);
      if (lastrearr) {
        if (progress) {
          printf("\nDoing global rearrangements\n");
          printf("  !");
          for (j = 1; j <= (nonodes); j++)
            putchar('-');
          printf("!\n");
        }
        bestlike = bestyet;
        if (jumb == 1) {
          bstlike2 = bestlike;
          nextree = 1;
        }
        do {
          if (progress)
            printf("   ");
          gotlike = bestlike;
          for (j = 0; j < (nonodes); j++) {
            bestyet = -10.0 * spp * chars;
            item = treenode[j];
            if (item != root) {
              re_move(&item, &nufork);
              there = root;
              addpreorder(root, item, nufork);
              add(there, item, nufork);
            }
            if (progress)
              putchar('.');
          }
          if (progress)
            putchar('\n');
        } while (bestlike > gotlike);
      }
    }
    if (progress)
      putchar('\n');
    for (i = spp - 1; i >= 1; i--)
      re_move(&treenode[i], &dummy);
    if (jumb == njumble) {
      if (treeprint) {
        putc('\n', outfile);
        if (nextree == 2)
          fprintf(outfile, "One most parsimonious tree found:\n");
        else
          fprintf(outfile, "%6ld trees in all found\n", nextree - 1);
      }
      if (nextree > maxtrees + 1) {
        if (treeprint)
          fprintf(outfile, "here are the first%4ld of them\n", (long)maxtrees);
        nextree = maxtrees + 1;
      }
      if (treeprint)
        putc('\n', outfile);
      for (i = 0; i <= (nextree - 2); i++) {
        root = treenode[0];
        add(treenode[0], treenode[1], treenode[spp]);
        for (j = 3; j <= (spp); j++)
            add(treenode[bestrees[i][j - 1] - 1], treenode[j - 1],
            treenode[spp + j - 2]);
        if (noroot)
          reroot(treenode[outgrno - 1]);
        didreroot = (outgropt && noroot);
        evaluate(root);
        printree();
        describe();
        for (j = 1; j < (spp); j++)
          re_move(&treenode[j], &dummy);
      }
    }
  } else {
    fscanf(infile, "%ld%*[^\n]", &numtrees);
    getc(infile);
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }
    which = 1;
    while (which <= numtrees ) {
      treeread();
      didreroot = (outgropt && noroot);
      if (noroot)
        reroot(treenode[outgrno - 1]);
      evaluate(root);
      printree();
      describe();
      which++;
    }
    fprintf(outfile, "\n\n");
    if (numtrees > 1 && chars > 1 ) {
      if (numtrees > maxuser) {
        printf("TOO MANY USER-DEFINED TREES");
        printf("  test only performed in the first%4ld of them\n",
               (long)maxuser);
        num = maxuser;
      } else
        num = numtrees;
      fprintf(outfile, "Tree    Steps   Diff Steps   Its S.D.");
      fprintf(outfile, "   Significantly worse?\n\n");
      for (which = 1; which <= num; which++) {
        fprintf(outfile, "%4ld%10.1f", which, nsteps[which - 1]);
        if (which == minwhich)
          fprintf(outfile, "  <------- best\n");
        else {
          sumw = 0;
          sum = 0.0;
          sum2 = 0.0;
          for (j = 0; j < (chars); j++) {
            if (weight[j] > 0) {
              sumw += weight[j];
              sum += fsteps[which - 1][j] - fsteps[minwhich - 1][j];
              TEMP = fsteps[which - 1][j] - fsteps[minwhich - 1][j];
              sum2 += TEMP * TEMP;
            }
          }
          sd = sqrt(sumw / (sumw - 1.0) * (sum2 - sum * sum / sumw));
          fprintf(outfile, "%8.1f%15.5f", nsteps[which - 1] - minsteps, sd);
          if (sum > 1.95996 * sd)
            fprintf(outfile, "           Yes\n");
          else
            fprintf(outfile, "           No\n");
        }
      }
      fprintf(outfile, "\n\n");
    }
  }
  if (jumb == njumble) {
    if (progress) {
      printf("\nOutput written to output file\n\n");
      if (trout)
        printf("Trees also written onto tree file\n");
      putchar('\n');
    }
  }
free(place);
free(steps);
free(zeroanc);
free(oneanc);
free(fulzeroanc);
free(empzeroanc);
}  /* maketree */

