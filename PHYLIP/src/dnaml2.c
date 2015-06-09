#include "phylip.h"

#define maxcategs       9
#define maxtrees        10   /* maximum number of user trees for test        */
#define smoothings      7    /* number of passes through smoothing algorithm */
#define iterations      9    /* number of iterates for each branch           */
#define nmlngth         10   /* max. number of characters in species name    */
#define epsilon         0.0001   /* used in getthree, maeknewv */
#define ibmpc0          false
#define ansi0           true
#define vt520           false
#define down            2
#define over            60
#define point           "."

/* shared defined types and such: */
typedef enum {  A, C, G, T} base;
typedef double sitelike[(short)T - (short)A + 1];
typedef sitelike *ratelike;
typedef ratelike *phenotype;
typedef Char **sequence;
typedef double contribarr[maxcategs];
typedef Char naym[nmlngth];
typedef short longer[6];

typedef struct node {
  struct node *next, *back;
  boolean tip, iter;
  short number;
  phenotype x;
  naym nayme;
  double v;
  short xcoord, ycoord, ymin, ymax;
} node;

typedef struct tree {
  node **nodep;
  double likelihood;
  node *start;
} tree;

typedef double *lf[maxtrees];

typedef struct valrec {
  double rat, ratxi, ratxv, zz, z1, y1, ww1, zz1, ww2, zz2, z1zz, z1yy, xiz1,
	 xiy1xv, ww1zz1, vv1zz1, ww2zz2, vv2zz2;
} valrec;

typedef short val[maxcategs];

/* shared variables between dnaml.c and dnaml2.c                             */
extern short categs,endsite,sites,numsp,numsp2,jumb,lengths,njumble,weightsum,
            outgrno;
extern double *probcat;
extern double *rate;
extern contribarr *contribution;
extern short   *category,*weight,*alias,*ally,*location,*aliasweight;
extern double xi,xv,ttratio,freqa,freqc,freqg,freqt,freqr,freqy,freqar,
              freqcy,freqgr,freqty,lambda,fracchange;
extern short   *enterorder;
extern boolean auto_,usertree,global,progress,treeprint,outgropt,trout,
               ctgry,jumble,lngths;
extern lf l0gf;
extern tree curtree,bestree,bestree2;
extern FILE *infile, *outfile, *treefile;
extern longer seed;

/* Local variables for maketree, propagated globally for c version: */
short k, nextsp, numtrees, which, maxwhich;
double dummy, tdelta, lnlike, slope, curv, maxlogl;
boolean succeeded, smoothed;
double x[3], lnl[3];
double l0gl[maxtrees];
valrec **tbl;
Char ch;
short col;

double randum(seed)
short *seed;
{  /* randum -- slow but machine independent */
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

void inittable()
{
  /* Define a lookup table. Precompute values and store them in a table */
  short i;

  tbl = (valrec **)Malloc(categs * sizeof(valrec *));
  for (i = 0; i < categs; i++)
    tbl[i] = (valrec *)Malloc(sizeof(valrec));
  for (i = 0; i < categs; i++) {
    tbl[i]->ratxi = rate[i] * xi;
    tbl[i]->ratxv = rate[i] * xv;
    tbl[i]->rat = rate[i];
  }
}  /* inittable */

double evaluate(p, saveit)
node *p;
boolean saveit;
{
  static contribarr like,nulike,term,clai;
  double sum, sum2, sumc, y, lz, y1, z1zz, z1yy, prod12, prod1, prod2, prod3,
	 sumterm, lterm;
  short i, j, lai;
  node *q;
  sitelike x1, x2;

  sum = 0.0;
  q = p->back;
  y = p->v;
  lz = -y;
  for (i = 0; i < categs; i++) {
    tbl[i]->zz = exp(tbl[i]->ratxi * lz);
    tbl[i]->z1 = exp(tbl[i]->ratxv * lz);
    tbl[i]->z1zz = tbl[i]->z1 * tbl[i]->zz;
    tbl[i]->z1yy = tbl[i]->z1 * (1.0 - tbl[i]->zz);
  }
  for (i = 0; i < endsite; i++) {
    for (j = 0; j < categs; j++) {
      if (y > 0.0) {
	y1 = 1.0 - tbl[j]->z1;
	z1zz = tbl[j]->z1zz;
	z1yy = tbl[j]->z1yy;
      } else {
	y1 = 0.0;
	z1zz = 1.0;
	z1yy = 0.0;
      }
      memcpy(x1, p->x[i][j], sizeof(sitelike));
      prod1 = freqa * x1[0] + freqc * x1[(short)C - (short)A] +
	      freqg * x1[(short)G - (short)A] + freqt * x1[(short)T - (short)A];
      memcpy(x2, q->x[i][j], sizeof(sitelike));
      prod2 = freqa * x2[0] + freqc * x2[(short)C - (short)A] +
	      freqg * x2[(short)G - (short)A] + freqt * x2[(short)T - (short)A];
      prod3 = (x1[0] * freqa + x1[(short)G - (short)A] * freqg) *
	      (x2[0] * freqar + x2[(short)G - (short)A] * freqgr) +
	  (x1[(short)C - (short)A] * freqc + x1[(short)T - (short)A] * freqt) *
	  (x2[(short)C - (short)A] * freqcy + x2[(short)T - (short)A] * freqty);
      prod12 = freqa * x1[0] * x2[0] +
	       freqc * x1[(short)C - (short)A] * x2[(short)C - (short)A] +
	       freqg * x1[(short)G - (short)A] * x2[(short)G - (short)A] +
	       freqt * x1[(short)T - (short)A] * x2[(short)T - (short)A];
      term[j] = z1zz * prod12 + z1yy * prod3 + y1 * prod1 * prod2;
    }
    sumterm = 0.0;
    for (j = 0; j < categs; j++)
      sumterm += probcat[j] * term[j];
    lterm = log(sumterm);
    for (j = 0; j < categs; j++)
      clai[j] = term[j] / sumterm;
    memcpy(contribution[i], clai, sizeof(contribarr));
    if (saveit && !auto_ && usertree && which <= maxtrees)
      l0gf[which - 1][i] = lterm;
    sum += weight[i] * lterm;
  }
  for (j = 0; j < categs; j++)
    like[j] = 1.0;
  for (i = 0; i < sites; i++) {
    if ((ally[i] > 0) && (location[ally[i]-1] > 0)) {
      sumc = 0.0;
      for (k = 1; k <= categs; k++)
        sumc += probcat[k - 1] * like[k - 1];
      sumc *= lambda;
      lai = location[ally[i] - 1];
      memcpy(clai, contribution[lai - 1], sizeof(contribarr));
      for (j = 0; j < categs; j++)
        nulike[j] = ((1.0 - lambda) * like[j] + sumc) * clai[j];
    } else {
      for (j = 0; j < categs; j++)
        nulike[j] = ((1.0 - lambda) * like[j] + lambda);
    }
    memcpy(like, nulike, sizeof(contribarr));
  }
  sum2 = 0.0;
  for (i = 0; i < categs; i++)
    sum2 += probcat[i] * like[i];
  sum += log(sum2);
  curtree.likelihood = sum;
  if (!saveit || auto_ || !usertree || which > maxtrees)
    return sum;
  l0gl[which - 1] = sum;
  if (which == 1) {
    maxwhich = 1;
    maxlogl = sum;
    return sum;
  }
  if (sum > maxlogl) {
    maxwhich = which;
    maxlogl = sum;
  }
  return sum;
}  /* evaluate */

void nuview(p)
node *p;
{
  short i, j;
  double lw1, lw2, yy1, yy2, ww1zz1, vv1zz1, ww2zz2, vv2zz2, vzsumr1,
	 vzsumr2, vzsumy1, vzsumy2, sum1, sum2, sumr1, sumr2,
	 sumy1, sumy2;
  node *q, *r;
  sitelike xx1, xx2, xx3;

  q = p->next->back;
  r = p->next->next->back;
  lw1 = -q->v;
  for (i = 0; i < categs; i++) {
    tbl[i]->ww1 = exp(tbl[i]->ratxi * lw1);
    tbl[i]->zz1 = exp(tbl[i]->ratxv * lw1);
    tbl[i]->ww1zz1 = tbl[i]->ww1 * tbl[i]->zz1;
    tbl[i]->vv1zz1 = (1.0 - tbl[i]->ww1) * tbl[i]->zz1;
  }
  lw2 = -r->v;
  for (i = 0; i < categs; i++) {
    tbl[i]->ww2 = exp(tbl[i]->ratxi * lw2);
    tbl[i]->zz2 = exp(tbl[i]->ratxv * lw2);
    tbl[i]->ww2zz2 = tbl[i]->ww2 * tbl[i]->zz2;
    tbl[i]->vv2zz2 = (1.0 - tbl[i]->ww2) * tbl[i]->zz2;
  }
  for (i = 0; i < endsite; i++) {
    for (j = 0; j < categs; j++) {
      ww1zz1 = tbl[j]->ww1zz1;
      vv1zz1 = tbl[j]->vv1zz1;
      yy1 = 1.0 - tbl[j]->zz1;
      ww2zz2 = tbl[j]->ww2zz2;
      vv2zz2 = tbl[j]->vv2zz2;
      yy2 = 1.0 - tbl[j]->zz2;

      memcpy(xx1, q->x[i][j], sizeof(sitelike));
      memcpy(xx2, r->x[i][j], sizeof(sitelike));
      sum1 = yy1 * (freqa * xx1[0] + freqc * xx1[(short)C - (short)A] +
	    freqg * xx1[(short)G - (short)A] + freqt * xx1[(short)T - (short)A]);
      sum2 = yy2 * (freqa * xx2[0] + freqc * xx2[(short)C - (short)A] +
	    freqg * xx2[(short)G - (short)A] + freqt * xx2[(short)T - (short)A]);
      sumr1 = freqar * xx1[0] + freqgr * xx1[(short)G - (short)A];
      sumr2 = freqar * xx2[0] + freqgr * xx2[(short)G - (short)A];
      sumy1 = freqcy * xx1[(short)C - (short)A] + freqty * xx1[(short)T - (short)A];
      sumy2 = freqcy * xx2[(short)C - (short)A] + freqty * xx2[(short)T - (short)A];
      vzsumr1 = vv1zz1 * sumr1;
      vzsumr2 = vv2zz2 * sumr2;
      vzsumy1 = vv1zz1 * sumy1;
      vzsumy2 = vv2zz2 * sumy2;
      xx3[0] = (sum1 + ww1zz1 * xx1[0] + vzsumr1) *
	       (sum2 + ww2zz2 * xx2[0] + vzsumr2);
      xx3[(short)C - (short)A] =
	(sum1 + ww1zz1 * xx1[(short)C - (short)A] + vzsumy1) *
	(sum2 + ww2zz2 * xx2[(short)C - (short)A] + vzsumy2);
      xx3[(short)G - (short)A] =
	(sum1 + ww1zz1 * xx1[(short)G - (short)A] + vzsumr1) *
	(sum2 + ww2zz2 * xx2[(short)G - (short)A] + vzsumr2);
      xx3[(short)T - (short)A] =
	(sum1 + ww1zz1 * xx1[(short)T - (short)A] + vzsumy1) *
	(sum2 + ww2zz2 * xx2[(short)T - (short)A] + vzsumy2);
      memcpy(p->x[i][j], xx3, sizeof(sitelike));
    }
  }
}  /* nuview */

void getthree(p)
node *p;
{
  /* compute likelihood, slope, curvature at a new triple of points */
  double tt, td;

  tt = p->v;
  if (tt < epsilon)
    tt = 2 * epsilon;
  td = tdelta;
  if (fabs(td) < epsilon) {
    if (td > 0)
      td = epsilon;
    else
      td = -epsilon;
  }
  if (fabs(td) >= tt - epsilon) {
    if (td > 0)
      td = tt / 2.0;
    else
      td = tt / -2.0;
  }
  p->v = tt + td;
  p->back->v = tt + td;
  x[0] = tt + td;
  lnl[0] = evaluate(p, false);
  p->v = tt - td;
  p->back->v = tt - td;
  x[2] = tt - td;
  lnl[2] = evaluate(p, false);
  p->v = tt;
  p->back->v = tt;
  x[1] = tt;
  lnl[1] = evaluate(p, false);
}  /* getthree */


void makenewv(p)
node *p;
{
  short it, imin, imax, i;
  double tt, oldlike, yold, yorig, ymin, ymax, s32, s21;
  boolean done, already;

  done = false;
  it = 1;
  tt = p->v;
  yorig = tt;
  tdelta = p->v / 10.0;
  getthree(p);
  while (it < iterations && !done) {
    ymax = lnl[0];
    imax = 1;
    for (i = 2; i <= 3; i++) {
      if (lnl[i - 1] > ymax) {
	ymax = lnl[i - 1];
	imax = i;
      }
    }
    if (imax != 2) {
      ymax = x[1];
      x[1] = x[imax - 1];
      x[imax - 1] = ymax;
      ymax = lnl[1];
      lnl[1] = lnl[imax - 1];
      lnl[imax - 1] = ymax;
    }
    tt = x[1];
    p->v = x[1];
    p->back->v = x[1];
    oldlike = lnl[1];
    yold = tt;
    s32 = (lnl[2] - lnl[1]) / (x[2] - x[1]);
    s21 = (lnl[1] - lnl[0]) / (x[1] - x[0]);
    curv = (s32 - s21) / ((x[2] - x[0]) / 2);
    slope = (s32 + s21) / 2 - curv * (x[2] - 2 * x[1] + x[0]) / 2;
    if (curv >= 0.0) {
      if (slope > 0)
	tdelta = fabs(tdelta);
      else
	tdelta = -fabs(tdelta);
    } else
      tdelta = -(slope / curv);
    if (tt + tdelta < epsilon) {
      if (tt / 2.0 > epsilon)
	tdelta = tt / -2.0;
      else
	tdelta = epsilon - tt;
    }
    tt += tdelta;
    p->v = tt;
    p->back->v = tt;
    done = (fabs(yold - p->v) < epsilon || fabs(tdelta) < epsilon);
    lnlike = evaluate(p, false);
    ymin = lnl[0];
    imin = 1;
    for (i = 2; i <= 3; i++) {
      if (lnl[i - 1] < ymin) {
	ymin = lnl[i - 1];
	imin = i;
      }
    }
    already = (tt == x[0]) || (tt == x[1]) || (tt == x[2]);
    if (!already && ymin < lnlike) {
      x[imin - 1] = tt;
      lnl[imin - 1] = lnlike;
    }
    if (already || lnlike < oldlike) {
      p->v = x[1];
      p->back->v = x[1];
      tdelta /= 2;
      curtree.likelihood = oldlike;
      lnlike = oldlike;
      tt = p->v;
    }
    it++;
  }
  smoothed = smoothed && fabs(tt-yorig) < epsilon;
  p->v = x[1];
  p->back->v = p->v;
}  /* makenewv */

void update(p)
node *p;
{
  if (!p->tip)
    nuview(p);
  if (!p->back->tip)
    nuview(p->back);
  if (p->iter)
    makenewv(p);
}  /* update */

void smooth(p, doupdate)
node *p;
boolean doupdate;
{
  if (doupdate)
    update (p);
  if (p->tip)
    return;
  smooth(p->next->back, doupdate);
  smooth(p->next->next->back, doupdate);
  nuview(p);
}  /* smooth */

void hookup(p, q)
node *p, *q;
{
  p->back = q;
  q->back = p;
}  /* hookup */

void in_sert(p, q)
node *p, *q;
{
  short i;
  node *r;

  r = p->next->next;
  hookup(r, q->back);
  hookup(p->next, q);
  q->v = 0.5 * q->v;
  q->back->v = q->v;
  r->v = q->v;
  r->back->v = r->v;
  smoothed = false;
  i = 1;
  while (i < smoothings && !smoothed) {
    smoothed = true;
    smooth(p, true);
    smooth(p->back, p->back->tip);
    i++;
  }
}  /* in_sert */

void re_move(p, q)
node **p, **q;
{
  /* remove p and record in q where it was */
  *q = (*p)->next->back;
  hookup(*q, (*p)->next->next->back);
  (*p)->next->back = NULL;
  (*p)->next->next->back = NULL;
  update(*q);
  update((*q)->back);
}  /* re_move */

void copynode(c, d)
node *c, *d;
{
  int i, j;

  memcpy(d->nayme, c->nayme, sizeof(naym));
  for (i=0;i<endsite;i++)
    for (j = 0; j < categs; j++)
      memcpy(d->x[i][j],c->x[i][j], sizeof(sitelike));
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
  short i, j=0;
  node *p, *q;

  for (i = 0; i < numsp; i++) {
    copynode(a->nodep[i], b->nodep[i]);
    if (a->nodep[i]->back) {
      if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->number - 1])
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->number - 1];
      else if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->number - 1]->next)
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->number - 1]->next;
      else
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->number - 1]->next->next;
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
  b->start = a->start;
}  /* copy */

void buildnewtip(m, tr)
short m;
tree *tr;
{
  node *p;

  p = tr->nodep[nextsp + numsp - 3];
  hookup(tr->nodep[m - 1], p);
  p->v = 0.1;
  p->back->v = 0.1;
}  /* buildnewtip */

void buildsimpletree(tr)
tree *tr;
{
  hookup(tr->nodep[enterorder[0] - 1], tr->nodep[enterorder[1] - 1]);
  tr->nodep[enterorder[0] - 1]->v = 0.1;
  tr->nodep[enterorder[0] - 1]->back->v = 0.1;
  tr->nodep[enterorder[1] - 1]->v = 0.1;
  tr->nodep[enterorder[1] - 1]->back->v = 0.1;
  buildnewtip(enterorder[2], tr);
  in_sert(tr->nodep[enterorder[2] - 1]->back, tr->nodep[enterorder[0] - 1]);
}  /* buildsimpletree */


void addtraverse(p, q, contin)
node *p, *q;
boolean contin;
{
  in_sert(p, q);
  numtrees++;
  if (evaluate(curtree.start, false) > bestree.likelihood)
    copy_(&curtree, &bestree);
  re_move(&p, &q);
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
	  addtraverse(r, q->next->back, true);
	  addtraverse(r, q->next->next->back, true);
	}
	q = q->back;
	if (!q->tip) {
	  addtraverse(r, q->next->back, true);
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


void coordinates(p, lengthsum,tipy,tipmax)
node *p;
double lengthsum;
short *tipy;
double *tipmax;
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;
  double xx;

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
    xx = fracchange * q->v;
    if (xx > 100.0)
      xx = 100.0;
    coordinates(q->back, lengthsum + xx, tipy,tipmax);
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
      } while (!(done || (p != curtree.start->back && r == p) ||
		 (p == curtree.start->back && r == p->next)));
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
      if (last->ycoord > i && first->ycoord < i &&
	  (i != p->ycoord || p == curtree.start)) {
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
/* Local variables for printree: */
  short tipy;
  double scale, tipmax;
  short i;

  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  tipmax = 0.0;
  coordinates(curtree.start->back, 0.0, &tipy,&tipmax);
  scale = 1.0 / (short)(tipmax + 1.000);
  for (i = 1; i <= (tipy - down); i++)
    drawline(i, scale);
  putc('\n', outfile);
}  /* printree */

#undef down
#undef over


/* Local variables for describe: */

void sigma(p, sumlr, s1, s2)
node *p;
double *sumlr, *s1, *s2;
{
  /* compute standard deviation */
  double tt, aa, s32, s21;

  tdelta = p->v / 10.0;
  getthree(p);
  s32 = (lnl[2] - lnl[1]) /  (x[2] - x[1]);
  s21 = (lnl[1] - lnl[0]) /  (x[1] - x[0]);
  curv = (s32 - s21) / ((x[2] - x[0]) / 2);
  slope = (s32 + s21) / 2 -   curv * (x[2] -2 * x[1] + x[0]) / 2;
  tt = p->v;
  p->v = 0.0;
  p->back->v = 0.0;
  aa = evaluate(p, false);
  p->v = tt;
  p->back->v = tt;
  (*sumlr) = evaluate(p, false) - aa;
  if (curv < -epsilon) {
    (*s1) = p->v + (-slope - sqrt(slope * slope -  3.841 * curv)) / curv;
    (*s2) = p->v + (-slope + sqrt(slope * slope -  3.841 * curv)) / curv;
  }
  else {
    (*s1) = -1.0;
    (*s2) = -1.0;
  }
}  /* sigma */

void describe(p)
node *p;
{
  /* print out information for one branch */
  short i;
  node *q;
  double s, sumlr, sigma1, sigma2;

  q = p->back;
  fprintf(outfile, "%4hd          ", q->number - numsp);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(p->nayme[i], outfile);
  } else
    fprintf(outfile, "%4hd      ", p->number - numsp);
  fprintf(outfile, "%15.5f", q->v * fracchange);
  if (p->iter) {
    sigma(q, &sumlr, &sigma1, &sigma2);
    if (sigma1 <= sigma2)
      fprintf(outfile, "     (     zero,    infinity)");
    else {
      fprintf(outfile, "     (");
      if (sigma2 <= 0.0)
        fprintf(outfile, "     zero");
      else
        fprintf(outfile, "%9.5f", sigma2 * fracchange);
      fprintf(outfile, ",%12.5f", sigma1 * fracchange);
      putc(')', outfile);
      }
    if (sumlr > 1.9205)
      fprintf(outfile, " *");
    if (sumlr > 2.995)
      putc('*', outfile);
    }
  putc('\n', outfile);
  if (!p->tip) {
    describe(p->next->back);
    describe(p->next->next->back);
  }
}  /* describe */

void summarize()
{
  /* print out branch length information and node numbers */
  short i, j, mx;
  double mode, sum;
  double like[maxcategs],nulike[maxcategs];
  val *mp;

  fprintf(outfile, "\nremember: ");
  if (outgropt)
    fprintf(outfile, "(although rooted by outgroup) ");
  fprintf(outfile, "this is an unrooted tree!\n\n");
  fprintf(outfile, "Ln Likelihood = %11.5f\n", curtree.likelihood);
  if (!usertree)
    fprintf(outfile, "\nExamined %4hd trees\n", numtrees);
  fprintf(outfile, "\n Between        And            Length");
  if (!(usertree && lngths))
    fprintf(outfile, "      Approx. Confidence Limits");
  fprintf(outfile, "\n");
  fprintf(outfile, " -------        ---            ------");
  if (!(usertree && lngths))
    fprintf(outfile, "      ------- ---------- ------");
  fprintf(outfile, "\n\n");
  describe(curtree.start->back->next->back);
  describe(curtree.start->back->next->next->back);
  describe(curtree.start);
  fprintf(outfile, "\n");
  if (!(usertree && lngths)) {
    fprintf(outfile, "     *  = significantly positive, P < 0.05\n");
    fprintf(outfile, "     ** = significantly positive, P < 0.01\n\n");
  }
  if (ctgry && categs > 1) {
    for (i = 0; i < categs; i++)
      like[i] = 1.0;
    mp = (val *)Malloc(sites*sizeof(val));
    for (i = sites - 1; i >= 0; i--) {
      sum = 0.0;
      for (j = 0; j < categs; j++) {
	nulike[j] = (1.0 - lambda + lambda * probcat[j]) * like[j];
	mp[i][j] = j + 1;
	for (k = 1; k <= categs; k++) {
	  if (k != j + 1) {
	    if (lambda * probcat[k - 1] * like[k - 1] > nulike[j]) {
	      nulike[j] = lambda * probcat[k - 1] * like[k - 1];
	      mp[i][j] = k;
	    }
	  }
	}
        if ((ally[i] > 0) && (location[ally[i]-1] > 0))
	  nulike[j] *= contribution[location[ally[i] - 1] - 1][j];
	sum += nulike[j];
      }
      for (j = 0; j < categs; j++)
	nulike[j] /= sum;
      memcpy(like, nulike, categs * sizeof(double));
    }
    mode = 0.0;
    mx = 1;
    for (i = 1; i <= categs; i++) {
      if (probcat[i - 1] * like[i - 1] > mode) {
	mx = i;
	mode = probcat[i - 1] * like[i - 1];
      }
    }
    fprintf(outfile,
      "Combination of categories that contributes the most to the likelihood:\n\n");
    for (i = 1; i <= nmlngth + 3; i++)
      putc(' ', outfile);
    for (i = 1; i <= sites; i++) {
      fprintf(outfile, "%hd", mx);
      if (i % 10 == 0)
	putc(' ', outfile);
      if (i % 60 == 0 && i != sites) {
	putc('\n', outfile);
	for (j = 1; j <= nmlngth + 3; j++)
	  putc(' ', outfile);
      }
      mx = mp[i - 1][mx - 1];
    }
    putc('\n', outfile);
    free(mp);
  }
  putc('\n', outfile);
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
  x = p->v * fracchange;
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
  /* scan forward until find character c */
  boolean done;

  done = false;
  while (!(done)) {
    if (c == ',') {
      if (ch == '(' || ch == ')' || ch == ':' ||
	  ch == ';') {
	printf("\nERROR IN USER TREE: UNMATCHED PARENTHESIS");
        printf(" OR MISSING COMMA\n OR NOT TRIFURCATED BASE\n");
	exit(-1);
      } else if (ch == ',')
	done = true;
    } else if (c == ')') {
      if (ch == '(' || ch == ',' || ch == ':' || ch == ';') {
	printf("\nERROR IN USER TREE: UNMATCHED PARENTHESIS OR NOT BIFURCATED NODE\n");
	exit(-1);
      } else if (ch == ')') {
	(*rparens)++;
	if ((*lparens) > 0 && (*lparens) == (*rparens)) {
	  if ((*lparens) == numsp - 2) {
	    getch(&ch);
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
    if (ch == ')')
      getch(&ch);
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
  while ( (unsigned short)digit <=9 || ch == '.') {
    if (ch == '.' )
      pointread = true;
    else {
      valyew = valyew * 10.0 + digit;
      if (pointread)
	divisor *= 10.0;
    }
    getch(&ch);
    digit = ch - ordzero;
  }
  if (lngths) {
    p->v = valyew / divisor / fracchange;
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
  node *q;
  short i, n;
  boolean found;
  Char str[nmlngth];

  getch(&ch);
  if (ch == '(') {
    (*lparens)++;
    if ((*lparens) > numsp - 2) {
      printf("\nERROR IN USER TREE: TOO MANY LEFT PARENTHESES\n");
      exit(-1);
    } else {
      (*nextnode)++;
      q = curtree.nodep[(*nextnode) - 1];
      hookup(p, q);
      addelement(q->next,nextnode,lparens,rparens,names,nolengths);
      findch(',',lparens,rparens);
	addelement(q->next->next,nextnode,lparens,rparens,names,nolengths);
	findch(')',lparens,rparens);
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
    hookup(curtree.nodep[n - 1], p);
    if (curtree.start->number > n)
      curtree.start = curtree.nodep[n - 1];
  }
  if (ch == ':') {
    processlength(p);
    *nolengths = false;
  }
}  /* addelement */

void treeread()
{
  short nextnode, lparens, rparens;
  boolean *names, nolengths;
  node *p;
  short i;

  curtree.start = curtree.nodep[numsp - 1];
  getch(&ch);
  if (ch != '(')
    return;
  nextnode = numsp + 1;
  p = curtree.nodep[nextnode - 1];
  names = (boolean *)Malloc(numsp*sizeof(boolean));
  for (i = 0; i < numsp; i++)
    names[i] = false;
  lparens = 1;
  rparens = 0;
  nolengths = true;
  for (i = 1; i <= 2; i++) {
    addelement(p, &nextnode,&lparens,&rparens,names,&nolengths);
    p = p->next;
    findch(',',&lparens,&rparens);
  }
  addelement(p, &nextnode,&lparens,&rparens,names,&nolengths);
  if (nolengths && lngths)
    printf("\nNO LENGTHS FOUND IN INPUT FILE WITH LENGTH OPTION CHOSEN\n");
  findch(')',&lparens,&rparens);
  fscanf(infile, "%*[^\n]");
  getc(infile);
  free(names);
}  /* treeread */

void nodeinit(p)
node *p;
{
  node *q, *r;
  short i, j;
  base b;

  if (p->tip)
    return;
  q = p->next->back;
  r = p->next->next->back;
  nodeinit(q);
  nodeinit(r);
  for (i = 0; i < endsite; i++) {
    for (j = 0; j < categs; j++) {
      for (b = A; (short)b <= (short)T; b = (base)((short)b + 1))
	p->x[i][j]
	  [(short)b - (short)A] = 0.5 * (q->x[i][j][(short)b - (short)A] + r->x[i]
				       [j][(short)b - (short)A]);
    }
  }
  if (p->iter)
    p->v = 0.1;
  if (p->back->iter)
    p->back->v = 0.1;
}  /* nodeinit */

void initrav(p)
node *p;
{
  if (p->tip)
    nodeinit(p->back);
  else {
    initrav(p->next->back);
    initrav(p->next->next->back);
  }
}  /* initrav */

void treevaluate()
{
  short i;

  initrav(curtree.start);
  initrav(curtree.start->back);
  for (i = 1; i <= smoothings * 4; i++)
    smooth(curtree.start->back, true);
  dummy = evaluate(curtree.start, true);
}  /* treevaluate */


void maketree()
{
  short i, j, k, num;
  double sum, sum2, sd;
  double TEMP;
  node *p;

  inittable();
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
    if (!auto_ && numtrees > 1 && weightsum > 1 ) {
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
	  for (j = 0; j < endsite; j++) {
	    sum += weight[j] * (l0gf[maxwhich - 1][j] - l0gf[i - 1][j]);
	    TEMP = l0gf[maxwhich - 1][j] - l0gf[i - 1][j];
	    sum2 += weight[j] * (TEMP * TEMP);
	  }
	  sd = sqrt(weightsum / (weightsum - 1.0) * (sum2 - sum * sum / weightsum));
	  fprintf(outfile, "%12.5f%12.4f", l0gl[i - 1] - maxlogl, sd);
	  if (sum > 1.95996 * sd)
	    fprintf(outfile, "           Yes\n");
	  else
	    fprintf(outfile, "            No\n");
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
      if (njumble > 1 && nextsp == numsp) {
        if (jumb == 1 || bestree2.likelihood < bestree.likelihood)
          copy_(&bestree, &bestree2);
      }
      if (nextsp == numsp && jumb == njumble) {
        if (njumble > 1) copy_(&bestree2, &curtree);
	curtree.start = curtree.nodep[outgrno - 1];
	printree();
        dummy = evaluate(curtree.start, false);
	summarize();
	if (trout && nextsp == numsp) {
	  col = 0;
	  treeout(curtree.start->back);
	}
      }
      nextsp++;
    }
  }
  if ( jumb < njumble)
    return;
  if (progress) {
    printf("\n\nOutput written to output file\n\n");
    if (trout)
      printf("Tree also written onto file\n");
    putchar('\n');
  }
  if (usertree)
    for (i = 0; i < maxtrees; i++)
      free(l0gf[i]);
  free(contribution);
  for (i = 0; i < categs; i++)
    free(tbl[i]);
  free(tbl);
  for (i = 0; i < numsp; i++) {
    for (j = 0; j < endsite; j++)
      free(curtree.nodep[i]->x[j]);
    free(curtree.nodep[i]->x);
  }
  for (i = numsp; i < numsp2; i++) {
    p = curtree.nodep[i];
    for (j = 1; j <= 3; j++) {
      for (k = 0; k < endsite; k++)
        free(p->x[k]);
      free(p->x);
      p = p->next;
    }
  }
  if (usertree)
    return;
  for (i = 0; i < numsp; i++) {
    for (j = 0; j < endsite; j++)
      free(bestree.nodep[i]->x[j]);
    free(bestree.nodep[i]->x);
  }
  for (i = numsp; i < numsp2; i++) {
    p = bestree.nodep[i];
    for (j = 1; j <= 3; j++) {
      for (k = 0; k < endsite; k++)
        free(p->x[k]);
      free(p->x);
      p = p->next;
    }
  }
  if (njumble <= 1)
    return;
  for (i = 0; i < numsp; i++) {
    for (j = 0; j < endsite; j++)
      free(bestree2.nodep[i]->x[j]);
    free(bestree2.nodep[i]->x);
  }
  for (i = numsp; i < numsp2; i++) {
    p = bestree2.nodep[i];
    for (j = 1; j <= 3; j++) {
      for (k = 0; k < endsite; k++)
        free(p->x[k]);
      free(p->x);
      p = p->next;
    }
  }
}  /* maketree */

