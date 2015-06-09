#include "phylip.h"

/* version 3.57c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define maxcutter       8    /* maximum number of bases in a site           */
#define maxtrees        8    /* maximum number of user trees                */
#define smoothings      10   /* number of passes in smoothing algorithm     */
#define iterations      10   /* number of iterates of EM for each branch    */
#define nmlngth         10   /* number of characters max. in species name   */

#define epsilon         0.00001 /* used in update                           */
#define extrap0         100.0   /* extrapolation factor to speed iteration  */
#define initialv        0.1     /* starting value of branch length          */

#define ibmpc0          false
#define ansi0           true
#define vt520           false
#define down            2
#define over            60
#define numsp2          (numsp * 2 - 2)


typedef double sitelike[maxcutter + 1];
typedef sitelike *phenotype;
typedef Char **sequence;
typedef Char naym[nmlngth];
typedef short longer[6];
typedef double **transmatrix;
typedef transmatrix *transptr;

typedef struct node {
  struct node *next, *back;
  boolean tip, iter, initialized;
  short branchnum, number;
  phenotype x;
  naym nayme;
  double v;
  short xcoord, ycoord, ymin, ymax;
} node;

typedef struct tree {
  node **nodep;
  transptr trans, transprod;
  double likelihood;
  node *start;
} tree;

/* Local variables for maketree, propagated globally for C version: */
short        nextsp, numtrees, which, maxwhich, col;
double      maxlogl;
boolean     succeeded, smoothed;
transmatrix tempmatrix, tempprod;
sitelike    pie;
double      l0gl[maxtrees];
double     *l0gf[maxtrees];
Char ch;

/* variables declared in the other segment */
extern boolean     trunc8,usertree,trout,lengths,treeprint,usertree,global,
                   progress,outgropt,jumble;
extern short       endsite,enzymes,weightsum,numsp,outgrno,sitelength,
                   jumb,njumble;
extern tree        curtree,priortree,bestree,bestree2;
extern short       *weight, *alias, *aliasweight, *enterorder;
extern FILE       *infile, *outfile, *treefile;
extern double      extrapol;
extern longer      seed;


double randum(seed)
short *seed;
{
  /* random number generator -- slow but machine independent */
  short i, j, k, sum;
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


void copymatrix(tomat,frommat)
     transmatrix frommat,tomat;
{
  int i,j;
  for (i=0;i<=sitelength;++i){
    for (j=0;j<=sitelength;++j)
      tomat[i][j] = frommat[i][j];
  }
}

void setuppi()
{
  /* set up equilibrium probabilities of being a given
     number of bases away from a restriction site */
  short i;
  double sum;

  pie[0] = 1.0;
  sum = pie[0];
  for (i = 1; i <= sitelength; i++) {
   pie[i] = 3 * pie[i - 1] * (sitelength - i + 1) / i;
    sum += pie[i];
  }
  for (i = 0; i <= sitelength; i++)
    pie[i] /= sum;
}  /* setuppi */

void maketrans(p)
double p;
{
  /* make transition matrix, product matrix with change
     probability p.  Put the results in tempmatrix, tempprod */
  short i, j, k, m1, m2;
  double sump, sumn, pover3, pijk, nijk, f;
  double binom1[maxcutter + 1], binom2[maxcutter + 1];

  pover3 = p / 3;
  f = 2 * pover3 / (1 - pover3);
  for (i = 0; i <= sitelength; i++) {
    if (p > 1.0 - epsilon)
      p = 1.0 - epsilon;
    binom1[0] = exp((sitelength - i) * log(1 - p));
    for (k = 1; k <= sitelength - i; k++)
      binom1[k] = binom1[k - 1] * (p / (1 - p)) * (sitelength - i - k + 1) / k;
    binom2[0] = exp(i * log(1 - pover3));
    for (k = 1; k <= i; k++)
      binom2[k] = binom2[k - 1] * (pover3 / (1 - pover3)) * (i - k + 1) / k;
    for (j = 0; j <= sitelength; ++j) {
      sump = 0.0;
      sumn = 0.0;
      if (i - j > 0)
        m1 = i - j;
      else
        m1 = 0;
      if (sitelength - j < i)
        m2 = sitelength - j;
      else
        m2 = i;
      for (k = m1; k <= m2; k++) {
        nijk = j - i + k * 2 + f * (i - k);
        pijk = binom1[j - i + k] * binom2[k];
        sump += pijk;
        sumn += pijk * nijk;
      }
      tempmatrix[i][j] = sump;
      tempprod[i][j] = sumn;
    }
  }
}  /* maketrans */

void branchtrans(i, p)
short i;
double p;
{
  /* make branch transition matrix, product matrices for branch
     i with probability of change p*/
  maketrans(p);
  copymatrix(curtree.trans[i - 1], tempmatrix);
  copymatrix(curtree.transprod[i - 1], tempprod);
}  /* branchtrans */

double evaluate(tr)
tree *tr;
{
  /* evaluates the likelihood, using info. at one branch */
  double sum, sum2, y, liketerm, like0, lnlike0, term;
  short i, j, k;
  node *p, *q;
  sitelike x1, x2;

  sum = 0.0;
  p = tr->start;
  q = p->back;
  y = p->v;
  maketrans(y);
  memcpy(x1, p->x[0], sizeof(sitelike));
  memcpy(x2, q->x[0], sizeof(sitelike));
  if (trunc8) {
    like0 = 0.0;
    for (j = 0; j <= sitelength; j++) {
      liketerm = pie[j] * x1[j];
      for (k = 0; k <= sitelength; k++)
        like0 += liketerm * tempmatrix[j][k] * x2[k];
    }
    lnlike0 = log(enzymes * (1.0 - like0));
  }
  for (i = 1; i <= endsite; i++) {
    memcpy(x1, p->x[i], sizeof(sitelike));
    memcpy(x2, q->x[i], sizeof(sitelike));
    sum2 = 0.0;
    for (j = 0; j <= sitelength; j++) {
      liketerm = pie[j] * x1[j];
      for (k = 0; k <= sitelength; k++)
        sum2 += liketerm * tempmatrix[j][k] * x2[k];
    }
    term = log(sum2);
    if (trunc8)
      term -= lnlike0;
    if (usertree && which <= maxtrees)
      l0gf[which - 1][i - 1] = term;
    sum += weight[i] * term;
  }
  if (usertree && which <= maxtrees) {
    l0gl[which - 1] = sum;
    if (which == 1) {
      maxwhich = 1;
      maxlogl = sum;
    } else if (sum > maxlogl) {
      maxwhich = which;
      maxlogl = sum;
    }
  }
  tr->likelihood = sum;
  return sum;
}  /* evaluate */

void nuview(p)
node *p;
{
  /* recompute fractional likelihoods for one part of tree */
  short i, j, k, lowlim;
  double sumq, sumr;
  node *q, *r;
  sitelike xq, xr, xp;
  transmatrix tempq, tempr;
  tempq = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  tempr = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  for (i=0;i<=sitelength;++i){
    tempq[i] = (double *)Malloc((sitelength+1) * sizeof (double));
    tempr[i] = (double *)Malloc((sitelength+1) * sizeof (double));
  }
  if (trunc8)
    lowlim = 0;
  else
    lowlim = 1;
  q = p->next->back;
  r = p->next->next->back;
  copymatrix(tempq,curtree.trans[q->branchnum - 1]);
  copymatrix(tempr,curtree.trans[r->branchnum - 1]);
  for (i = lowlim; i <= endsite; i++) {
    memcpy(xq, q->x[i], sizeof(sitelike));
    memcpy(xr, r->x[i], sizeof(sitelike));
    for (j = 0; j <= sitelength; j++) {
      sumq = 0.0;
      sumr = 0.0;
      for (k = 0; k <= sitelength; k++) {
        sumq += tempq[j][k] * xq[k];
        sumr += tempr[j][k] * xr[k];
      }
      xp[j] = sumq * sumr;
    }
    memcpy(p->x[i], xp, sizeof(sitelike));
  }
  for (i=0;i<=sitelength;++i){
    free(tempq[i]);
    free(tempr[i]);
  }
free(tempq);
free(tempr);
}  /* nuview */


void makenewv(p)
node *p;
{
  /* EM algorithm improvement of a branch length */
  short i, j, k, lowlim, it, ite;
  double sum, sumlike, sumprod, liketerm, liket, y, yold, yorig, like0,
         prod0, like, oldlike, olderlike, extrap;
  boolean done;
  node *q;
  sitelike xx1, xx2;

  extrap = extrapol;
  q = p->back;
  y = p->v;
  yorig = y;
  if (trunc8)
    lowlim = 0;
  else
    lowlim = 1;
  done = false;
  oldlike = 0.0;
  olderlike = 0.0;
  it = 1;
  ite = 1;
  while ((it < iterations) && (!done)) {
    like = 0.0;
    maketrans(y);
    yold = y;
    sumlike = 0.0;
    for (i = lowlim; i <= endsite; i++) {
      memcpy(xx1, p->x[i], sizeof(sitelike));
      memcpy(xx2, q->x[i], sizeof(sitelike));
      sum = 0.0;
      sumprod = 0.0;
      for (j = 0; j <= sitelength; j++) {
        liket = xx1[j] * pie[j];
        for (k = 0; k <= sitelength; k++) {
          liketerm = liket * xx2[k];
          sumprod += tempprod[j][k] * liketerm;
          sum += tempmatrix[j][k] * liketerm;
        }
      }
      if (i == 0) {
        like0 = sum;
        prod0 = sumprod;
      } else {
        sumlike += weight[i] * sumprod / sum;
        like += weight[i] * log(sum);
      }
    }
    if (trunc8 && fabs(like0 - 1.0) > 1.0e-10)
      like -= weightsum * log(enzymes * (1.0 - like0));
    y = sumlike / (sitelength * weightsum);
    if (trunc8)
      y = (1.0 - like0) * y + prod0 / sitelength;
    if (ite >= 3 && like > oldlike && like - oldlike < oldlike - olderlike) {
      extrap = (oldlike - olderlike) / (2.0 * oldlike - olderlike - like) *
               extrapol;
      ite = 1;
    } else
      extrap = extrapol;
    if (extrap < 0.0)
      extrap = extrapol;
    olderlike = oldlike;
    oldlike = like;
    y = extrap * (y - yold) + yold;
    if (y < epsilon)
      y = 10.0 * epsilon;
    if (y >= 0.75)
      y = 0.75;
    done = fabs(y-yorig) < epsilon;
    it++;
    ite++;
  }
  smoothed = (smoothed && done);
  p->v = y;
  q->v = y;
  branchtrans(p->branchnum, y);
}  /* makenewv */

void update(p)
node *p;
{
  /* improve branch length and views for one branch */

  if (!p->tip)
    nuview(p);
  if (!p->back->tip)
    nuview(p->back);
  if (p->iter)
    makenewv(p);
}  /* update */

void smooth(p)
node *p;
{
  /* update nodes throughout the tree, recursively */
  update(p);
  if (!p->tip) {
    smooth(p->next->back);
    smooth(p->next->next->back);
    update (p);
  }
}  /* smooth */

void hookup(p, q)
node *p, *q;
{
  /* connect two nodes */
  p->back = q;
  q->back = p;
}  /* hookup */

void insert_(p, q)
node *p, *q;
{
  /* insert a subtree into a branch, improve lengths in tree */
  short i, m, n;
  node *r;

  r = p->next->next;
  hookup(r, q->back);
  hookup(p->next, q);
  if (q->v >= 0.75)
    q->v = 0.75;
  else
    q->v = 0.75 * (1 - sqrt(1 - 1.333333 * q->v));
  q->back->v = q->v;
  r->v = q->v;
  r->back->v = r->v;
  if (q->branchnum == q->number) {
    m = q->branchnum;
    n = r->number;
  } else {
    m = r->number;
    n = q->branchnum;
  }
  q->branchnum = m;
  q->back->branchnum = m;
  r->branchnum = n;
  r->back->branchnum = n;
  branchtrans(q->branchnum, q->v);
  branchtrans(r->branchnum, r->v);
  smoothed = false;
  i = 1;
  while (i < smoothings && !smoothed) {
    smoothed = true;
    smooth(p);
    smooth(p->back);
    i++;
  }
}  /* insert */

void re_move(p, q)
node **p, **q;
{
  /* remove p and record in q where it was */
  *q = (*p)->next->back;
  hookup(*q, (*p)->next->next->back);
  (*q)->back->branchnum = (*q)->branchnum;
  branchtrans((*q)->branchnum, 0.75*(1 - (1 - 1.333333*(*q)->v)
                                        * (1 - 1.333333*(*p)->next->v)));
  (*p)->next->back = NULL;
  (*p)->next->next->back = NULL;
  update(*q);
  update((*q)->back);
}  /* re_move */

void copynode(c, d)
node *c, *d;
{
  /* copy a node */

  d->branchnum = c->branchnum;
  memcpy(d->x, c->x, (endsite+1)*sizeof(sitelike));
  memcpy(d->nayme, c->nayme, sizeof(naym));
  d->v = c->v;
  d->iter = c->iter;
  d->xcoord = c->xcoord;
  d->ycoord = c->ycoord;
  d->ymin = c->ymin;
  d->ymax = c->ymax;
}  /* copynode */

void copy_(a, b)
tree *a, *b;
{
  /* copy a tree */
  short i,j;
  node *p, *q;

  for (i = 0; i < numsp; i++) {
    copynode(a->nodep[i], b->nodep[i]);
    if (a->nodep[i]->back) {
      if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->number - 1])
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->number - 1];
      else if (a->nodep[i]->back
                == a->nodep[a->nodep[i]->back->number - 1]->next)
          b->nodep[i]->back = b->nodep[a->nodep[i]->back->number - 1]->next;
        else
          b->nodep[i]->back
                         = b->nodep[a->nodep[i]->back->number - 1]->next->next;
    }
    else b->nodep[i]->back = NULL;
  }
  for (i = numsp; i < numsp2; i++) {
    p = a->nodep[i];
    q = b->nodep[i];
    for (j = 1; j <= 3; j++) {
      copynode(p, q);
      if (p->back) {
        if (p->back == a->nodep[p->back->number - 1])
          q->back = b->nodep[p->back->number - 1];
        else if (p->back == a->nodep[p->back->number - 1]->next)
          q->back = b->nodep[p->back->number - 1]->next;
        else
          q->back = b->nodep[p->back->number - 1]->next->next;
      }
      else
        q->back = NULL;
      p = p->next;
      q = q->next;
    }
  }
  b->likelihood = a->likelihood;
  for (i=0;i<numsp2;++i){
       copymatrix(b->trans[i],a->trans[i]);
       copymatrix(b->transprod[i],a->transprod[i]);
       }
  b->start = a->start;
}  /* copy */

void buildnewtip(m, tr)
short m;
tree *tr;
{
  /* set up a new tip and interior node it is connected to */
  node *p;

  p = tr->nodep[nextsp + numsp - 3];
  hookup(tr->nodep[m - 1], p);
  p->v = initialv;
  p->back->v = initialv;
  branchtrans(m, initialv);
  p->branchnum = m;
  p->next->branchnum = p->number;
  p->next->next->branchnum = p->number;
  p->back->branchnum = m;
}  /* buildnewtip */

void buildsimpletree(tr)
tree *tr;
{
  /* set up and adjust branch lengths of a three-species tree */
  hookup(tr->nodep[enterorder[0] - 1], tr->nodep[enterorder[1] - 1]);
  tr->nodep[enterorder[0] - 1]->v = initialv;
  tr->nodep[enterorder[1] - 1]->v = initialv;
  branchtrans(enterorder[1], initialv);
  tr->nodep[enterorder[0] - 1]->branchnum = 2;
  tr->nodep[enterorder[1] - 1]->branchnum = 2;
  buildnewtip(enterorder[2], tr);
  insert_(tr->nodep[enterorder[2] - 1]->back, tr->nodep[enterorder[1] - 1]);
}  /* buildsimpletree */

void addtraverse(p, q, contin)
node *p, *q;
boolean contin;
{
int i;
double like;
  /* try adding p at q, proceed recursively through tree */
  insert_(p, q);
  numtrees++;
  like = evaluate(&curtree);
  if (like > bestree.likelihood) {
    copy_(&curtree, &bestree);
    bestree.likelihood = like;
    succeeded = true;
  }
  copy_(&priortree, &curtree);
  if (!q->tip && contin) {
    addtraverse(p, q->next->back, contin);
    addtraverse(p, q->next->next->back, contin);
  }
}  /* addtraverse */

void rearrange(p)
node *p;
{
  /* rearranges the tree, globally or locally */
  node *q, *r;

  if (!p->tip && !p->back->tip) {
    r = p->next->next;
    re_move(&r, &q);
    copy_(&curtree, &priortree);
    addtraverse(r, q->next->back, global && nextsp == numsp);
    addtraverse(r, q->next->next->back, global && nextsp == numsp);
    copy_(&bestree, &curtree);
    if (global && nextsp == numsp && progress)
      putchar('.');
    if (global && nextsp == numsp && !succeeded) {
      if (r->back->tip) {
        r = r->next->next;
        re_move(&r, &q);
        q = q->back;
        if (!q->tip) {
          copy_(&curtree, &priortree);
          addtraverse(r, q->next->back, true);
          copy_(&curtree, &priortree);
          addtraverse(r, q->next->next->back, true);
        }
        q = q->back;
        if (!q->tip) {
          copy_(&curtree, &priortree);
          addtraverse(r, q->next->back, true);
          copy_(&curtree, &priortree);
          addtraverse(r, q->next->next->back, true);
        }
        copy_(&bestree, &curtree);
      }
    }
  }
  if (!p->tip) {
    rearrange(p->next->back);
    rearrange(p->next->next->back);
  }
}  /* rearrange */


void coordinates(p, lengthsum, tipy,tipmax,x)
node *p;
double lengthsum;
short *tipy;
double *tipmax, *x;

{
  /* establishes coordinates of nodes */
  node *q, *first, *last;

  if (p->tip) {
    p->xcoord = (short)(over * lengthsum + 0.5);
    p->ycoord = (*tipy);
    p->ymin = (*tipy);
    p->ymax = (*tipy);
    (*tipy) += down;
    if (lengthsum > (*tipmax))
      (*tipmax) = lengthsum;
    return;
  }
  q = p->next;
  do {
    (*x) = -0.75 * log(1.0 - 1.333333 * q->v);
    coordinates(q->back, lengthsum + (*x),tipy,tipmax,x);
    q = q->next;
  } while ((p == curtree.start->back || p != q) &&
           (p != curtree.start->back || p->next != q));
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (short)(over * lengthsum + 0.5);
  if (p == curtree.start->back)
    p->ycoord = p->next->next->back->ycoord;
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */

void drawline(i, scale)
short i;
double scale;
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q;
  short n, j;
  boolean extra;
  node *r, *first, *last;
  boolean done;

  p = curtree.start->back;
  q = curtree.start->back;
  extra = false;
  if (i == p->ycoord && p == curtree.start->back) {
    if (p->number - numsp >= 10)
      fprintf(outfile, "-%2hd", p->number - numsp);
    else
      fprintf(outfile, "--%hd", p->number - numsp);
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
      } while (!(done || p != curtree.start->back && r == p ||
                 p == curtree.start->back && r == p->next));
      first = p->next->back;
      r = p;
      while (r->next != p)
        r = r->next;
      last = r->back;
      if (p == curtree.start->back)
        last = p->back;
    }
    done = (p->tip || p == q);
    n = (short)(scale * (q->xcoord - p->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if (q->ycoord == i && !done) {
      if (p->ycoord != q->ycoord)
        putc('+', outfile);
      else
        putc('-', outfile);
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        if (q->number - numsp >= 10)
          fprintf(outfile, "%2hd", q->number - numsp);
        else
          fprintf(outfile, "-%hd", q->number - numsp);
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
    if (q != p)
      p = q;
  } while (!done);
  if (p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(p->nayme[j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */

void printree()
{
  /* prints out diagram of the tree */
  short tipy,i;
  double scale, tipmax, x;

  putc('\n', outfile);
  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  tipmax = 0.0;
  coordinates(curtree.start->back, 0.0, &tipy,&tipmax,&x);
  scale = 1.0 / (tipmax + 1.000);
  for (i = 1; i <= tipy - down; i++)
    drawline(i, scale);
  putc('\n', outfile);
}  /* printree */


double sigma(q,sumlr)
node *q;
double *sumlr;
{
  /* get 1.95996 * approximate standard error of branch length */
double sump, sumr, sums, sumc, p, pover3, pijk, Qjk, liketerm, f;
double  slopef,curvef;
  short i, j, k, m1, m2;
  double binom1[maxcutter + 1], binom2[maxcutter + 1];
  transmatrix Prob, slopeP, curveP;
  node *r;
  sitelike x1, x2;
  double TEMP, TEMP1;
  Prob   = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  slopeP = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  curveP = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  for (i=0;i<=sitelength;++i){
    Prob[i]   = (double *)Malloc((sitelength+1) * sizeof(double));
    slopeP[i] = (double *)Malloc((sitelength+1) * sizeof(double));
    curveP[i] = (double *)Malloc((sitelength+1) * sizeof(double));
  }
  p = q->v;
  pover3 = p / 3;
  for (i = 0; i <= sitelength; i++) {
    binom1[0] = exp((sitelength - i) * log(1 - p));
    for (k = 1; k <= (sitelength - i); k++)
      binom1[k] = binom1[k - 1] * (p / (1 - p)) * (sitelength - i - k + 1) / k;
    binom2[0] = exp(i * log(1 - pover3));
    for (k = 1; k <= i; k++)
      binom2[k] = binom2[k - 1] * (pover3 / (1 - pover3)) * (i - k + 1) / k;
    for (j = 0; j <= sitelength; j++) {
      sump = 0.0;
      sums = 0.0;
      sumc = 0.0;
      if (i - j > 0)
        m1 = i - j;
      else
        m1 = 0;
      if (sitelength - j < i)
        m2 = sitelength - j;
      else
        m2 = i;
      for (k = m1; k <= m2; k++) {
        pijk = binom1[j - i + k] * binom2[k];
        sump += pijk;
        sums += pijk * ((j - i + k * 2) / p + (j - sitelength + k) / (1 - p) +
                        (k - i) / (3 - p));
        TEMP = 1 - p;
        TEMP1 = 3 - p;
        sumc += pijk * ((i - j - k * 2) / (p * p) + (j - sitelength + k) /
                          (TEMP * TEMP) + (k - i) / (TEMP1 * TEMP1));
      }
      Prob[i][j] = sump;
      slopeP[i][j] = sums;
      curveP[i][j] = sumc;
    }
  }
  (*sumlr) = 0.0;
  sumc = 0.0;
  sums = 0.0;
  r = q->back;
  for (i = 1; i <= endsite; i++) {
    f = 0.0;
    slopef = 0.0;
    curvef = 0.0;
    sumr = 0.0;
    memcpy(x1, q->x[i], sizeof(sitelike));
    memcpy(x2, r->x[i], sizeof(sitelike));
    for (j = 0; j <= sitelength; j++) {
      liketerm = pie[j] * x1[j];
      sumr += liketerm * x2[j];
      for (k = 0; k <= sitelength; k++) {
        Qjk = liketerm * x2[k];
        f += Qjk * Prob[j][k];
        slopef += Qjk * slopeP[j][k];
        curvef += Qjk * curveP[j][k];
      }
    }
    (*sumlr) += weight[i] * log(f / sumr);
    sums += weight[i] * slopef / f;
    TEMP = slopef / f;
    sumc += weight[i] * (curvef / f - TEMP * TEMP);
  }
  if (trunc8) {
    f = 0.0;
    slopef = 0.0;
    curvef = 0.0;
    sumr = 0.0;
    memcpy(x1, q->x[0], sizeof(sitelike));
    memcpy(x2, r->x[0], sizeof(sitelike));
    for (j = 0; j <= sitelength; j++) {
      liketerm = pie[j] * x1[j];
      sumr += liketerm * x2[j];
      for (k = 0; k <= sitelength; k++) {
        Qjk = liketerm * x2[k];
        f += Qjk * Prob[j][k];
        slopef += Qjk * slopeP[j][k];
        curvef += Qjk * curveP[j][k];
      }
    }
    (*sumlr) += weightsum * log((1.0 - sumr) / (1.0 - f));
    sums += weightsum * slopef / (1.0 - f);
    TEMP = slopef / (1.0 - f);
    sumc += weightsum * (curvef / (1.0 - f) + TEMP * TEMP);
  }
for (i=0;i<sitelength;++i){
  free(Prob[i]);
  free(slopeP[i]);
  free(curveP[i]);
}
free(Prob);
free(slopeP);
free(curveP);
  if (sumc < -1.0e-6)
    return ((-sums - sqrt(sums * sums - 3.841 * sumc)) / sumc);
  else
    return -1.0;
}  /* sigma */

void describe(p)
node *p;
{
  /* print out information on one branch */
  double sumlr;
  short i;
  node *q;
  double s;

  q = p->back;
  fprintf(outfile, "%4hd      ", q->number - numsp);
  fprintf(outfile, "    ");
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(p->nayme[i], outfile);
  } else
    fprintf(outfile, "%4hd      ", p->number - numsp);
  if (q->v >= 0.75)
    fprintf(outfile, "     infinity");
  else
    fprintf(outfile, "%13.5f", -0.75 * log(1 - 1.333333 * q->v));
  if (p->iter) {
    s = sigma(q, &sumlr);
    if (s < 0.0)
      fprintf(outfile, "     (     zero,    infinity)");
    else {
      fprintf(outfile, "     (");
      if (q->v - s <= 0.0)
        fprintf(outfile, "     zero");
      else
        fprintf(outfile, "%9.5f", -0.75 * log(1 - 1.333333 * (q->v - s)));
      putc(',', outfile);
      if (q->v + s >= 0.75)
        fprintf(outfile, "    infinity");
      else
        fprintf(outfile, "%12.5f", -0.75 * log(1 - 1.333333 * (q->v + s)));
      putc(')', outfile);
      }
    if (sumlr > 1.9205)
      fprintf(outfile, " *");
    if (sumlr > 2.995)
      putc('*', outfile);
    }
  else
    fprintf(outfile, "            (not varied)");
  putc('\n', outfile);
  if (!p->tip) {
    describe(p->next->back);
    describe(p->next->next->back);
  }
}  /* describe */

void summarize()
{
  /* print out information on branches of tree */

  fprintf(outfile, "\nremember: ");
  if (outgropt)
    fprintf(outfile, "(although rooted by outgroup) ");
  fprintf(outfile, "this is an unrooted tree!\n\n");
  fprintf(outfile, "Ln Likelihood = %11.5f\n\n", curtree.likelihood);
  if (!usertree)
    fprintf(outfile, "Examined %4hd trees\n", numtrees);
  fprintf(outfile, " \n");
  fprintf(outfile, " Between        And            Length");
  fprintf(outfile, "      Approx. Confidence Limits\n");
  fprintf(outfile, " -------        ---            ------");
  fprintf(outfile, "      ------- ---------- ------\n");
  describe(curtree.start->back->next->back);
  describe(curtree.start->back->next->next->back);
  describe(curtree.start);
  fprintf(outfile, "\n     *  = significantly positive, P < 0.05\n");
  fprintf(outfile, "     ** = significantly positive, P < 0.01\n\n\n");
}  /* summarize */

void treeout(p)
node *p;
{
  /* write out file with representation of final tree */
  short i, n, w;
  Char c;
  double x;

  if (p->tip) {
    n = 0;
    for (i = 1; i <= nmlngth; i++) {
      if (p->nayme[i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = p->nayme[i];
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
    if (col > 45) {
      putc('\n', treefile);
      col = 0;
    }
    treeout(p->next->next->back);
    if (p == curtree.start->back) {
      putc(',', treefile);
      col++;
      if (col > 45) {
        putc('\n', treefile);
        col = 0;
      }
      treeout(p->back);
    }
    putc(')', treefile);
    col++;
  }
  if (p->v >= 0.75)
    x = -1.0;
  else
    x = -0.75 * log(1 - 1.333333 * p->v);
  if (x > 0.0)
    w = (short)(0.43429448222 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (short)(0.43429448222 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if (p == curtree.start->back)
    fprintf(treefile, ";\n");
  else {
    fprintf(treefile, ":%*.5f", (int)(w + 7), x);
    col += w + 8;
  }
}  /* treeout */


void getch(c)
Char *c;
{
  /* get next nonblank character */
  do {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    *c = getc(infile);
    if (*c == '\n')
      *c = ' ';
  } while (*c == ' ');
}  /* getch */

void findch(c,lparens,rparens)
Char c;
short *lparens,*rparens;
{
  /* read forward to next occurrence of character c */
  boolean done;

  done = false;
  while (!done) {
    if (c == ',') {
      if (ch == '(' || ch == ')' || ch == ':' || ch == ';') {
        printf("\nERROR IN USER TREE: UNMATCHED PARENTHESIS");
	printf("OR MISSING COMMA\n OR NOT TRIFURCATED BASE\n");
	exit(-1);
      } else if (ch == ',')
        done = true;
    } else if (c == ')') {
      if (ch == '(' || ch == ',' || ch == ':' || ch == ';') {
        printf( "\nERROR IN USER TREE:");
	printf(" UNMATCHED PARENTHESIS OR NOT BIFURCATED NODE\n");
	exit(-1);
      } else if (ch == ')') {
        (*rparens)++;
        if ((*lparens) > 0 && (*lparens) == (*rparens)) {
          if ((*lparens) == numsp - 2) {
            if (eoln(infile)) {
              fscanf(infile, "%*[^\n]");
              getc(infile);
            }
            ch = getc(infile);
            if (ch == '\n')
              ch = ' ';
            if (ch != ';') {
              printf("\nERROR IN USER TREE:");
	      printf(" UNMATCHED PARENTHESIS OR MISSING SEMICOLON\n");
	      exit(-1);
            }
          }
        }
	done = true;
      }
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

void processlength(p)
node *p;
{
  short digit, ordzero;
  double valyew, divisor;
  boolean pointread;

  ordzero = '0';
  pointread = false;
  valyew = 0.0;
  divisor = 1.0;
  getch(&ch);
  digit = ch - ordzero;
  while (((unsigned short)digit <= 9) | ch == '.') {
    if (ch == '.')
      pointread = true;
    else {
      valyew = valyew * 10.0 + digit;
      if (pointread)
	divisor *= 10.0;
    }
    getch(&ch);
    digit = ch - ordzero;
  }
  if (lengths) {
    p->v = 0.75*(1.0-exp(-1.333333 * valyew / divisor));
    p->back->v = p->v;
    p->iter = false;
    p->back->iter = false;
  }
}  /* processlength */


void addelement(p, nextnode,lparens,rparens,names,nolengths)
node *p;
short *nextnode,*lparens,*rparens;
boolean *names, *nolengths;

{
  /* add one node to user-defined tree */
  node *q;
  short i, j, n;
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
  if (ch == '(') {
    (*lparens)++;
    if ((*lparens) > ( numsp - 2 )) {
      printf("\nERROR IN USER TREE: TOO MANY LEFT PARENTHESES\n");
      exit(-1);
    } else {
      (*nextnode)++;
      q = curtree.nodep[(*nextnode) - 1];
      n = (*nextnode);
      hookup(p, q);
      addelement(q->next, nextnode,lparens,rparens,names,nolengths);
      findch(',',lparens,rparens);
      addelement(q->next->next, nextnode,lparens,rparens,names,nolengths);
      if (*nolengths && lengths)
        printf("NO LENGTHS FOUND IN INPUT FILE WITH LENGTH OPTION CHOSEN");
      findch(')',lparens,rparens);
      for (i = 0; i <= endsite; i++) {
	for (j = 0; j <= sitelength; j++) {
	  q->x[i][j] = 1.0;
	  q->next->x[i][j] = 1.0;
	  q->next->next->x[i][j] = 1.0;
	}
      }
    }
  } else {
    for (i = 0; i < nmlngth; i++)
      str[i] = ' ';
    n = 1;
    do {
      if (ch == '_')
        ch = ' ';
      str[n - 1] = ch;
      if (eoln(infile)) {
        fscanf(infile, "%*[^\n]");
        getc(infile);
      }
      ch = getc(infile);
      if (ch == '\n')
        ch = ' ';
      n++;
    } while (ch != ':' && ch != ',' && ch != ')' && n <= nmlngth);
    n = 1;
    do {
      found = true;
      for (i = 0; i < nmlngth; i++)
        found = (found && str[i] == curtree.nodep[n - 1]->nayme[i]);
      if (found) {
        q = curtree.nodep[n - 1];
        if (names[n - 1] == false)
          names[n - 1] = true;
        else {
          printf("\nERROR IN USER TREE: DUPLICATE NAME FOUND -- ");
          for (i = 0; i < nmlngth; i++)
            putchar(curtree.nodep[n - 1]->nayme[i]);
          putchar('\n');
	  exit(-1);
        }
      } else
        n++;
    } while (!(n > numsp || found));
    if (n > numsp) {
      printf("Cannot find species: ");
      for (i = 0; i < nmlngth; i++)
        putchar(str[i]);
      putchar('\n');
    }
    hookup(p, q);
    if (curtree.start->number > n)
      curtree.start = curtree.nodep[n - 1];
  }
  p->branchnum = q->number;
  q->branchnum = q->number;
  q->v = initialv;
  p->v = initialv;
  branchtrans(q->branchnum, q->v);
  if (ch == ':') {
    processlength(p);
    *nolengths = false;
  }
}  /* addelement */

void treeread()
{
  /* read a user-defined tree */
/* Local variables for treeread: */
  short nextnode, lparens, rparens;
  boolean *names, nolengths;
  node *p;
  short i, j;

  curtree.start = curtree.nodep[numsp - 1];
  do {
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
  } while (ch == ' ');
  if (ch != '(')
    return;
  names = (boolean *)Malloc(numsp*sizeof(boolean));
  for (i = 0; i < numsp; i++)
    names[i] = false;
  lparens = 1;
  rparens = 0;
  nolengths = true;
  nextnode = numsp + 1;
  p = curtree.nodep[nextnode - 1];
  for (i = 1; i <= 2; i++) {
    addelement(p, &nextnode,&lparens,&rparens,names,&nolengths);
    p = p->next;
    findch(',',&lparens,&rparens);
  }
  addelement(p, &nextnode,&lparens,&rparens,names,&nolengths);
  findch(')',&lparens,&rparens);
  for (i = 0; i <= endsite; i++) {
    for (j = 0; j <= sitelength; j++) {
      p->x[i][j] = 1.0;
      p->next->x[i][j] = 1.0;
      p->next->next->x[i][j] = 1.0;
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  free(names);
}  /* treeread */


void travinit(p)
node *p;
{
  /* traverse to set up initial values */
  if (p == NULL)
    return;
  if (p->tip)
    return;
  if (p->initialized)
    return;
  travinit(p->next->back);
  travinit(p->next->next->back);
  nuview(p);
  p->initialized = true;
}  /* travinit */


void travsp(p)
node *p;
{
  /* traverse to find tips */
  if (p == curtree.start)
    travinit(p);
  if (p->tip)
    travinit(p->back);
  else {
    travsp(p->next->back);
    travsp(p->next->next->back);
  }
}  /* travsp */


void treevaluate()
{
  /* find maximum likelihood branch lengths of user tree */
  short i;
  double dummy;

  for (i = 1; i <= numsp; i++)
    curtree.nodep[i-1]->initialized = false;
  for (i = numsp+1; i <= numsp2; i++) {
    curtree.nodep[i-1]->initialized = false;
    curtree.nodep[i-1]->next->initialized = false;
    curtree.nodep[i-1]->next->next->initialized = false;
  }
  for (i = 1; i <= smoothings * 4; i++)
    smooth(curtree.start->back);
  dummy = evaluate(&curtree);
}  /* treevaluate */


void maketree()
{
  /* construct and rearrange tree */
  short i, j, k, num;
  double sum, sum2, sd;
  double TEMP;
  node *p;

  tempmatrix = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  for (i=0;i<=sitelength;++i)
    tempmatrix[i] = (double *)Malloc((sitelength+1) * sizeof (double));
  tempprod = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  for (i=0;i<=sitelength;++i)
    tempprod[i] = (double *)Malloc((sitelength+1) * sizeof (double));

  for (i = 0; i < maxtrees; i++)
    l0gf[i] = (double *)Malloc(endsite*sizeof(double));
  setuppi();
  if (usertree) {
    fscanf(infile, "%hd%*[^\n]", &numtrees);
    getc(infile);
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }
    which = 1;
    while (which <= numtrees) {
      treeread();
      treevaluate();
      printree();
      summarize();
      if (trout) {
        col = 0;
        treeout(curtree.start->back);
      }
      which++;
    }
    putc('\n', outfile);
    if (numtrees > 1 && weightsum > 1) {
      fprintf(outfile, "Tree    Ln L      Diff Ln L     Its S.D.");
      fprintf(outfile, "   Significantly worse?\n\n");
      if (numtrees > maxtrees)
        num = maxtrees;
      else
        num = numtrees;
      for (i = 1; i <= num; i++) {
        fprintf(outfile, "%3hd%12.5f", i, l0gl[i - 1]);
        if (maxwhich == i)
          fprintf(outfile, "  <------ best\n");
        else {
          sum = 0.0;
          sum2 = 0.0;
          for (j = 1; j <= endsite; j++) {
            sum += weight[j] * (l0gf[maxwhich - 1][j - 1]-l0gf[i - 1][j - 1]);
            TEMP = l0gf[maxwhich - 1][j - 1] - l0gf[i - 1][j - 1];
            sum2 += weight[j] * (TEMP * TEMP);
          }
          sd = sqrt(weightsum / (weightsum - 1.0)
                    * (sum2 - sum * sum / weightsum));
          fprintf(outfile, "%12.5f%12.4f", l0gl[i - 1] - maxlogl, sd);
          if (sum > 1.95996 * sd)
            fprintf(outfile, "           Yes\n");
          else
            fprintf(outfile, "           No\n");
        }
      }
      fprintf(outfile, "\n\n");
    }
  } else {
    for (i = 1; i <= numsp; i++)
      enterorder[i - 1] = i;
    if (jumble) {
      for (i = 0; i < numsp; i++) {
        j = (short)(randum(seed) * numsp) + 1;
        k = enterorder[j - 1];
        enterorder[j - 1] = enterorder[i];
        enterorder[i] = k;
      }
    }
    if (progress) {
      printf("\nAdding species:\n");
      printf("   ");
      for (i = 0; i < nmlngth; i++)
        putchar(curtree.nodep[enterorder[0] - 1]->nayme[i]);
      printf("\n   ");
      for (i = 0; i < nmlngth; i++)
        putchar(curtree.nodep[enterorder[1] - 1]->nayme[i]);
      printf("\n   ");
      for (i = 0; i < nmlngth; i++)
        putchar(curtree.nodep[enterorder[2] - 1]->nayme[i]);
      putchar('\n');
    }
    nextsp = 3;
    buildsimpletree(&curtree);
    curtree.start = curtree.nodep[enterorder[0] - 1];
    if (jumb == 1) numtrees = 1;
    nextsp = 4;
    while (nextsp <= numsp) {
      buildnewtip(enterorder[nextsp - 1], &curtree);
      bestree.likelihood = -999999.0;
      copy_(&curtree, &priortree);
      addtraverse(curtree.nodep[enterorder[nextsp - 1] - 1]->back,
                  curtree.start->back, true);
      copy_(&bestree, &curtree);
      if (progress) {
        printf("   ");
        for (j = 0; j < nmlngth; j++)
          putchar(curtree.nodep[enterorder[nextsp - 1] - 1]->nayme[j]);
        putchar('\n');
      }
      if (global && nextsp == numsp) {
        if (progress) {
          printf("Doing global rearrangements\n");
          printf("   ");
        }
      }
      succeeded = true;
      while (succeeded) {
        succeeded = false;
        rearrange(curtree.start->back);
      }
      if (njumble > 1) {
        if (jumb == 1 && nextsp == numsp)
          copy_(&bestree, &bestree2);
        else if (nextsp == numsp) {
          if (bestree2.likelihood < bestree.likelihood)
            copy_(&bestree, &bestree2);
        }
      }
      if (nextsp == numsp && jumb == njumble) {
        if (njumble > 1) copy_(&bestree2, &curtree);
        curtree.start = curtree.nodep[outgrno - 1];
        printree();
        summarize();
        if (trout && nextsp == numsp) {
          col = 0;
          treeout(curtree.start->back);
        }
      }
      nextsp++;
    }
  }
  if (jumb == njumble) {
    if (progress) {
      printf("\nOutput written to output file\n\n");
      if (trout)
        printf("Tree also written onto file\n");
      putchar('\n');
    }
    for (i = 0; i < numsp; i++)
      free(curtree.nodep[i]->x);
    for (i = numsp; i < numsp2; i++) {
      p = curtree.nodep[i];
      for (j = 1; j <= 3; j++) {
        free(p->x);
        p = p->next;
      }
    }
    if (!usertree) {
      for (i = 0; i < numsp; i++)
        free(priortree.nodep[i]->x);
      for (i = numsp; i < numsp2; i++) {
        p = priortree.nodep[i];
        for (j = 1; j <= 3; j++) {
          free(p->x);
          p = p->next;
        }
      }
      for (i = 0; i < numsp; i++)
        free(bestree.nodep[i]->x);
      for (i = numsp; i < numsp2; i++) {
        p = bestree.nodep[i];
        for (j = 1; j <= 3; j++) {
          free(p->x);
          p = p->next;
        }
      }
      if (njumble > 1) {
        for (i = 0; i < numsp; i++)
          free(bestree2.nodep[i]->x);
        for (i = numsp; i < numsp2; i++) {
          p = bestree2.nodep[i];
          for (j = 1; j <= 3; j++) {
            free(p->x);
            p = p->next;
          }
        }
      }
    }
    for (i = 0; i < maxtrees; i++)
      free(l0gf[i]);
  }
}  /* maketree */


