#include "phylip.h"

#define maxcategories   5    /* maximum number of site types                 */
#define smoothings      2    /* number of passes through smoothing algorithm */
#define iterations      10   /* number of iterates for each branch           */
#define nmlngth         10   /* number of characters max. in species name    */

#define epsilon         0.0001   /* used in makenewv, getthree, update */

#define point           '.'

#define ibmpc0          false
#define ansi0           true
#define vt520           false


typedef enum {
  A, C, G, T
} base;
typedef double sitelike[(short)T - (short)A + 1];
typedef sitelike **phenotype;
typedef double *contribarr;
typedef char **sequence;


typedef short longer[6];

typedef struct node {
  struct node *next, *back;
  phenotype x;
  double tyme, v;
  short number, xcoord, ycoord, ymin, ymax;
  boolean tip, haslength, initialized;
} node;


typedef struct tree {
  node **nodep;
  double likelihood;
  node *root;
} tree;


extern short      categs,spp,endsite,sites,nonodes,weightsum,datasets,
                  njumble,jumb;
extern double     rate[maxcategories];
extern double     xi,xv,freqa,freqc,freqg,freqt,freqr,freqy,freqar,freqcy,
                  freqgr,freqty,lambda,lambda1,fracchange;
extern boolean    usertree,globle,jumble,lengths,trout,weights,ctgry,auto_,
                  printdata,progress,treeprint;
extern tree       curtree,bestree,bestree2;
extern contribarr probcat;
extern double      **contribution;
extern short       *alias,*ally,*location,*aliasweight;
extern FILE       *infile,*outfile,*treefile;
extern            Char **naym;
extern            short *enterorder;
extern            longer seed;

/* end pseudo-include file */

#define down            2
#define over            60


/* Local variables for maketree, propagated globally for C version: */
short    k, which, maxwhich,col;
double  like, bestyet, tdelta, lnlike, slope, curv, maxlogl;
boolean lastrearr, smoothed;
double  *l0gl;
double  x[3], lnl[3];
double  expon1i[maxcategories], expon1v[maxcategories],
        expon2i[maxcategories], expon2v[maxcategories];
node    *there;
double  **l0gf;
Char ch;

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


void exmake(lz, n)
double lz;
short n;
{
  /* pretabulate tables of exponentials so need not do for each site */
  short i;
  double rat;

  for (i = 0; i < categs; i++) {
    rat = rate[i];
    switch (n) {

    case 1:
      expon1i[i] = exp(rat * xi * lz);
      expon1v[i] = exp(rat * xv * lz);
      break;

    case 2:
      expon2i[i] = exp(rat * xi * lz);
      expon2v[i] = exp(rat * xv * lz);
      break;
    }
  }
}  /* exmake */

double evaluate(p)
node *p;
{
  double sum, sum2, sumc, y, lz, y1, z1zz, z1yy, prod12, prod1, prod2, prod3,
	 sumterm, lterm;
  static contribarr like,term,clai;
  short i, j, lai;
  boolean recorduser;
  node *q, *r;
  sitelike x1, x2;

  if (!like){
    like   = (double *)Malloc(categs * sizeof(double));
    term   = (double *)Malloc(categs * sizeof(double));
    clai   = (double *)Malloc(categs * sizeof(double));
  }
  recorduser = (!auto_ && usertree);
  sum = 0.0;
  if (p != curtree.root) {
    r = p;
    q = p->back;
    y = fabs(r->tyme - q->tyme);
  } else {
    r = p->next->back;
    q = p->next->next->back;
    y = r->tyme + q->tyme - 2 * p->tyme;
  }
  lz = -y;
  exmake(lz, 1L);
  for (i = 0; i < endsite; i++) {
    for (j = 0; j < categs; j++) {
      y1 = 1.0 - expon1v[j];
      z1zz = expon1i[j] * expon1v[j];
      z1yy = (1.0 - expon1i[j]) * expon1v[j];
      memcpy(x1, r->x[i][j], sizeof(sitelike));
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
    memcpy(contribution[i], clai, categs * sizeof(double ));
    if (recorduser)
      l0gf[which - 1][i] = lterm;
    sum += aliasweight[i] * lterm;
  }
  if (auto_) {
    for (j = 0; j < categs; j++)
      like[j] = 1.0;
    for (i = 0; i < sites; i++) {
      if ((ally[i] > 0) && (location[ally[i]-1] > 0)) {
        sumc = 0.0;
        for (k = 1; k <= categs; k++)
          sumc += probcat[k - 1] * like[k - 1];
        sumc *= lambda;
        lai = location[ally[i] - 1];
        memcpy(clai, contribution[lai - 1], categs * sizeof(double));
        for (j = 0; j < categs; j++)
	  like[j] = (lambda1 * like[j] + sumc) * clai[j];
      } else {
        for (j = 0; j < categs; j++)
	  like[j] = lambda1 * like[j] + lambda;
      }
    }
    sum2 = 0.0;
    for (i = 0; i < categs; i++)
      sum2 += probcat[i] * like[i];
    sum += log(sum2);
  }
  curtree.likelihood = sum;
  if (!recorduser)
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
  double lw, ww1, ww2, zz1, zz2, vv1, vv2, yy1, yy2, sum1, sum2, sumr1, sumr2,
	 sumy1, sumy2;
  node *q, *r;
  sitelike xx1, xx2, xx3;

  q = p->next->back;
  r = p->next->next->back;
  if (q != NULL)
    lw = -fabs(p->tyme - q->tyme);
  else
    lw = 0.0;
  exmake(lw, 1L);
  if (r != NULL)
    lw = -fabs(p->tyme - r->tyme);
  else
    lw = 0.0;
  exmake(lw, 2L);
  for (i = 0; i < endsite; i++) {
    for (j = 0; j < categs; j++) {
      ww1 = expon1i[j];
      zz1 = expon1v[j];
      vv1 = 1.0 - ww1;
      yy1 = 1.0 - zz1;
      ww2 = expon2i[j];
      zz2 = expon2v[j];
      vv2 = 1.0 - ww2;
      yy2 = 1.0 - zz2;
      if (q != NULL)
	memcpy(xx1, q->x[i][j], sizeof(sitelike));
      if (r != NULL)
	memcpy(xx2, r->x[i][j], sizeof(sitelike));
      if (q == NULL) {
	xx1[0] = 1.0;
	xx1[(short)C - (short)A] = 1.0;
	xx1[(short)G - (short)A] = 1.0;
	xx1[(short)T - (short)A] = 1.0;
      }
      if (r == NULL) {
	xx2[0] = 1.0;
	xx2[(short)C - (short)A] = 1.0;
	xx2[(short)G - (short)A] = 1.0;
	xx2[(short)T - (short)A] = 1.0;
      }
      sum1 = yy1 * (freqa * xx1[0] + freqc * xx1[(short)C - (short)A] +
	    freqg * xx1[(short)G - (short)A] + freqt * xx1[(short)T - (short)A]);
      sum2 = yy2 * (freqa * xx2[0] + freqc * xx2[(short)C - (short)A] +
	    freqg * xx2[(short)G - (short)A] + freqt * xx2[(short)T - (short)A]);
      sumr1 = vv1 * (freqar * xx1[0] + freqgr * xx1[(short)G - (short)A]);
      sumr2 = vv2 * (freqar * xx2[0] + freqgr * xx2[(short)G - (short)A]);
      sumy1 = vv1 * (freqcy * xx1[(short)C - (short)A] +
		     freqty * xx1[(short)T - (short)A]);
      sumy2 = vv2 * (freqcy * xx2[(short)C - (short)A] +
		     freqty * xx2[(short)T - (short)A]);
      xx3[0] = (sum1 + zz1 * (ww1 * xx1[0] + sumr1)) *
	       (sum2 + zz2 * (ww2 * xx2[0] + sumr2));
      xx3[(short)C - (short)A] =
	(sum1 + zz1 * (ww1 * xx1[(short)C - (short)A] + sumy1)) *
	(sum2 + zz2 * (ww2 * xx2[(short)C - (short)A] + sumy2));
      xx3[(short)G - (short)A] =
	(sum1 + zz1 * (ww1 * xx1[(short)G - (short)A] + sumr1)) *
	(sum2 + zz2 * (ww2 * xx2[(short)G - (short)A] + sumr2));
      xx3[(short)T - (short)A] =
	(sum1 + zz1 * (ww1 * xx1[(short)T - (short)A] + sumy1)) *
	(sum2 + zz2 * (ww2 * xx2[(short)T - (short)A] + sumy2));
      memcpy(p->x[i][j], xx3, sizeof(sitelike));
    }
  }
}  /* nuview */

void getthree(p)
node *p;
{
  /* compute likelihood, slope, curvature at a new triple of points */
  double tt, td;
  node *q, *r, *sdown;

  sdown = curtree.nodep[p->number - 1];
  q = sdown->next->back;
  r = sdown->next->next->back;
  sdown = sdown->back;
  td = fabs(tdelta);
  if (td < epsilon)
    td = epsilon;
  tt = p->tyme;
  if (q->tyme - tt < td)
    td = (q->tyme - tt) / 2.0;
  if (r->tyme - tt < td)
    td = (r->tyme - tt) / 2.0;
  if (sdown != NULL) {
    if (tt - sdown->tyme < td)
      td = (tt - sdown->tyme) / 2.0;
  }
  if (-tt > epsilon && td < epsilon)
    td = epsilon;
  if (tt == -td)
    td += epsilon;
  p->tyme = tt + td;
  x[0] = tt + td;
  nuview(p);
  lnl[0] = evaluate(p);
  p->tyme = tt - td;
  x[2] = tt - td;
  nuview(p);
  lnl[2] = evaluate(p);
  p->tyme = tt;
  x[1] = tt;
  nuview(p);
  lnl[1] = evaluate(p);
}  /* getthree */


void makenewv(p)
node *p;
{
  /* improve a node time */
  short it, imin, imax, i;
  double tt, tfactor, tlow, thigh, oldlike, ymin, ymax, s32, s21, yold;
  boolean done, already;
  node *q, *r, *s, *sdown;

  s = curtree.nodep[p->number - 1];
  q = s->next->back;
  r = s->next->next->back;
  sdown = s->back;
  if (s == curtree.root)
    tlow = -10.0;
  else
    tlow = sdown->tyme;
  thigh = q->tyme;
  if (r->tyme < thigh)
    thigh = r->tyme;
  done = (thigh - tlow < epsilon);
  it = 1;
  tdelta = (thigh - tlow) / 10.0;
  if (s->tyme - tlow < tdelta)
    tdelta = (s->tyme - tdelta)/2.0;
  if (thigh - s->tyme < tdelta)
    tdelta = (thigh - s->tyme)/2.0;
  getthree(s);
  tfactor = 1.0;
  while (it < iterations && !done) {
    tt = s->tyme;
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
    oldlike = lnl[1];
    yold = tt;
    s32 = (lnl[2] - lnl[1]) / (x[2] - x[1]);
    s21 = (lnl[1] - lnl[0]) / (x[1] - x[0]);
    if (fabs(x[2] - x[0]) > epsilon)
      curv = (s32 - s21) / ((x[2] - x[0]) / 2);
    else
      curv = 0.0;
    slope = (s32 + s21) / 2 - curv *	  (x[2] - 2 * x[1] + x[0]) / 4;
    if (curv >= 0.0) {
      if (slope < 0)
	tdelta = -fabs(tdelta);
      else
	tdelta = fabs(tdelta);
    } else
      tdelta = -(tfactor * slope / curv);
    if (tt + tdelta <= tlow + epsilon)
      tdelta = (tlow - tt) / 2;
    if (tt + tdelta >= thigh - epsilon)
      tdelta = (thigh - tt) / 2;
    tt += tdelta;
    done = (fabs(yold - tt) < epsilon || fabs(tdelta) < epsilon);
    s->tyme = tt;
    nuview(s);
    lnlike = evaluate(s);
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
      tt = x[1];
      tfactor /= 2;
      tdelta /= 2;
      curtree.likelihood = oldlike;
      lnlike = oldlike;
    } else
      tfactor = 1.0;
    p->tyme = tt;
    p->next->tyme = tt;
    p->next->next->tyme = tt;
    nuview(p);
    nuview(p->next);
    nuview(p->next->next);
    it++;
  }
  smoothed = smoothed && done;
}  /* makenewv */

void update(p)
node *p;
{
  /* improve time and recompute views at a node */
  if (p == NULL)
    return;
  if (p->back != NULL) {
    if (!p->back->tip)
      nuview(p->back);
  }
  if (p->next->back != NULL) {
    if (!p->next->back->tip)
      nuview(p->next->back);
  }
  if (p->next->next->back != NULL) {
    if (!p->next->next->back->tip)
      nuview(p->next->next->back);
  }
  if (!lengths) {
    makenewv(p);
    return;
  }
  nuview(p);
  nuview(p->next);
  nuview(p->next->next);
}  /* update */

void smooth(p)
node *p;
{
  if (p == NULL)
    return;
  if (p->tip)
    return;
  update(p);
  smooth(p->next->back);
  update(p);
  smooth(p->next->next->back);
  update(p);
}  /* smooth */

void add(below, newtip, newfork)
node *below, *newtip, *newfork;
{
  /* inserts the nodes newfork and its descendant, newtip, in the tree. */
  short i;
  boolean done;
  node *p;

  below = curtree.nodep[below->number - 1];
  if (below->back != NULL)
    below->back->back = newfork;
  newfork->back = below->back;
  below->back = newfork->next->next;
  newfork->next->next->back = below;
  newfork->next->back = newtip;
  newtip->back = newfork->next;
  if (newtip->tyme < below->tyme)
    p = newtip;
  else p = below;
  newfork->tyme = p->tyme;
  if (curtree.root == below)
    curtree.root = newfork;
  if (newfork->back != NULL) {
    if (p->tyme > newfork->back->tyme)
      newfork->tyme = (p->tyme + newfork->back->tyme) / 2.0;
    else newfork->tyme = p->tyme - epsilon;
    newfork->next->tyme = newfork->tyme;
    newfork->next->next->tyme = newfork->tyme;
    do {
      p = curtree.nodep[p->back->number - 1];
      done = (p == curtree.root);
      if (!done)
	done = (curtree.nodep[p->back->number - 1]->tyme < p->tyme - epsilon);
      if (!done) {
	curtree.nodep[p->back->number - 1]->tyme = p->tyme - epsilon;
        curtree.nodep[p->back->number - 1]->next->tyme = p->tyme - epsilon;
        curtree.nodep[p->back->number - 1]->next->next->tyme = p->tyme - epsilon;
      }
    } while (!done);
  } else {
      newfork->tyme = newfork->tyme - 0.1;
      newfork->next->tyme = newfork->tyme;
      newfork->next->next->tyme = newfork->tyme;
    }
  nuview(newfork);
  nuview(newfork->next);
  nuview(newfork->next->next);
  smoothed = false;
  i = 1;
  while (i < smoothings && !smoothed) {
    smoothed = true;
    smooth(newfork);
    smooth(newfork->back);
    i++;
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
  *fork = curtree.nodep[(*item)->back->number - 1];
  if (curtree.root == *fork) {
    if (*item == (*fork)->next->back)
      curtree.root = (*fork)->next->next->back;
    else
      curtree.root = (*fork)->next->back;
  }
  p = (*item)->back->next->back;
  q = (*item)->back->next->next->back;
  if (p != NULL)
    p->back = q;
  if (q != NULL)
    q->back = p;
  (*fork)->back = NULL;
  p = (*fork)->next;
  while (p != *fork) {
    p->back = NULL;
    p = p->next;
  }
  (*item)->back = NULL;
}  /* remove */

void copynode(c, d)
node *c, *d;
{
  short i, j;

  for (i=0;i<endsite;i++)
    for (j = 0; j < categs; j++)
      memcpy(d->x[i][j],c->x[i][j], sizeof(sitelike));
  d->tyme = c->tyme;
  d->v = c->v;
  d->xcoord = c->xcoord;
  d->ycoord = c->ycoord;
  d->ymin = c->ymin;
  d->ymax = c->ymax;
  d->haslength = c->haslength;
  d->initialized = c->initialized;
}  /* copynode */

void copy_(a, b)
tree *a, *b;
{
  short i, j=0;
  node *p, *q;

  for (i = 0; i < spp; i++) {
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
  for (i = spp; i < nonodes; i++) {
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
  b->root = a->root;
}  /* copy */

void tryadd(p, item,nufork)
node *p,**item,**nufork;

{
  /* temporarily adds one fork and one tip to the tree.
    if the location where they are added yields greater
    likelihood than other locations tested up to that
    time, then keeps that location as there */
  add(p, *item, *nufork);
  like = evaluate(p);
  if (lastrearr) {
    if (like >= bestree.likelihood)
      copy_(&curtree, &bestree);
  }
  if (like > bestyet) {
    bestyet = like;
    there = p;
  }
  re_move(item,nufork);
}  /* tryadd */

void addpreorder(p, item_, nufork_)
node *p, *item_, *nufork_;
{
  /* traverses a binary tree, calling PROCEDURE tryadd
    at a node before calling tryadd at its descendants */
  node *item, *nufork;

  item = item_;
  nufork = nufork_;
  if (p == NULL)
    return;
  tryadd(p, &item,&nufork);
  if (!p->tip) {
    addpreorder(p->next->back, item, nufork);
    addpreorder(p->next->next->back,item,nufork);
  }
}  /* addpreorder */


void tryrearr(p, success)
node *p;
boolean *success;
{
  /* evaluates one rearrangement of the tree.
    if the new tree has greater likelihood than the old
    one sets success = TRUE and keeps the new tree.
    otherwise, restores the old tree */
  node *frombelow, *whereto, *forknode;
  double oldlike;

  if (p == curtree.root)
    return;
  forknode = curtree.nodep[p->back->number - 1];
  if (forknode == curtree.root)
    return;
  oldlike = bestyet;
  if (forknode->next->back == p)
    frombelow = forknode->next->next->back;
  else
    frombelow = forknode->next->back;
  whereto = curtree.nodep[forknode->back->number - 1];
  re_move(&p, &forknode);
  nuview(whereto);
  add(whereto, p, forknode);
  like = evaluate(p);
  if (like <= oldlike) {
    re_move(&p, &forknode);
    add(frombelow, p, forknode);
  } else {
    (*success) = true;
    bestyet = like;
  }
}  /* tryrearr */

void repreorder(p, success)
node *p;
boolean *success;
{
  /* traverses a binary tree, calling PROCEDURE tryrearr
    at a node before calling tryrearr at its descendants */
  if (p == NULL)
    return;
  tryrearr(p, success);
  if (p->tip)
    return;
  if (!(*success))
    repreorder(p->next->back, success);
  if (!(*success))
    repreorder(p->next->next->back,success);
}  /* repreorder */

Local Void rearrange(r)
node **r;
{
  /* traverses the tree (preorder), finding any local
    rearrangement which decreases the number of steps.
    if traversal succeeds in increasing the tree's
    likelihood, PROCEDURE rearrange runs traversal again */
  boolean success;
  success = true;
  while (success) {
    success = false;
    repreorder(*r, &success);
  }
}  /* rearrange */

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

void findch(c,rparens)
Char c;
short *rparens;
{
  /* scan forward until find character c */
  boolean done;

  done = false;
  while (!(done)) {
    if (c == ',') {
      if (ch == '(' || ch == ')' || ch == ';') {
	printf("\nERROR IN USER TREE: UNMATCHED PARENTHESIS OR MISSING COMMA\n");
	exit(-1);
      } else if (ch == ',')
	done = true;
    } else if (c == ')') {
      if (ch == '(' || ch == ',' || ch == ';') {
	printf("\nERROR IN USER TREE: UNMATCHED PARENTHESIS OR NOT BIFURCATED NODE\n");
	exit(-1);
      } else {
	if (ch == ')') {
	  done = true;
	  (*rparens)++;
	}
      }
    } else if (c == ';') {
      if (ch != ';') {
	printf("\nERROR IN USER TREE: UNMATCHED PARENTHESIS OR MISSING SEMICOLON\n");
	exit(-1);
      } else
	done = true;
    }
    if ((done && ch == ')') || !(done))
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
  while (((unsigned short)digit <= 9) ||	 ch == point){
    if (ch == point)
      pointread = true;
    else {
      valyew = valyew * 10.0 + digit;
      if (pointread)
	divisor *= 10.0;
    }
    getch(&ch);
    digit = ch - ordzero;
  }
  p->v = valyew / divisor / fracchange;
}  /* processlength */

void addelement(p, nextnode,lparens,rparens,names)
node **p;
short *nextnode,*lparens,*rparens;
boolean *names;
{
  /* recursive procedure adds nodes to user-defined tree */
  node *q;
  short i, n;
  boolean found;
  Char str[nmlngth];

  getch(&ch);
  if (ch == '(') {
    if ((*lparens) >= spp - 1) {
      printf("\nERROR IN USER TREE: TOO MANY LEFT PARENTHESES\n");
      exit(-1);
    } else {
      (*lparens)++;
      (*nextnode)++;
      q = curtree.nodep[(*nextnode) - 1];
      addelement(&q->next->back, nextnode,lparens,rparens,names);
      q->next->back->back = q->next;
      findch(',',rparens);
      addelement(&q->next->next->back, nextnode,lparens,rparens,names);
      q->next->next->back->back = q->next->next;
      findch(')',rparens);
      *p = q;
    }
  } else  {
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
    if (ch == ')')
      (*rparens)++;
    n = 1;
    do {
      found = true;
      for (i = 0; i < nmlngth; i++)
	found = (found && str[i] == naym[n - 1][i]);
      if (found) {
	if (names[n - 1] == false) {
	  *p = curtree.nodep[n - 1];
	  names[n - 1] = true;
	} else {
	  printf("\nERROR IN USER TREE: DUPLICATE NAME FOUND -- ");
	  for (i = 0; i < nmlngth; i++)
	    putchar(naym[n - 1][i]);
	  putchar('\n');
	  exit(-1);
	}
      } else
	n++;
    } while (!(n > spp || found));
    if (n > spp) {
      printf("Cannot find species: ");
      for (i = 0; i < nmlngth; i++)
	putchar(str[i]);
      putchar('\n');
      exit(-1);
    }
  }
  if (!lengths)
    return;
  if (ch == ':' || ch == ';') {
    if ((*rparens) < spp && ch == ':')
      processlength(*p);
    return;
  }
  printf("ERROR IN USER TREE: UNMATCHED PARENTHESES\n");
  printf(" OR BRANCH LENGTHS OPTION CHOSEN IN MENU\n");
  printf(" WITH NO INPUT BRANCH LENGTHS\n");
  exit(-1);
}  /* addelement */

void tymetrav(p, x)
node *p;
double *x;
{
  /* set up times of nodes */
  if (!p->tip) {
    tymetrav(p->next->back, x);
    tymetrav(p->next->next->back,x);
  } else
    (*x)     = 0.0;
    p->tyme  = (*x);
    (*x)    -= p->v;
}  /* tymetrav */

void treeread()
{
  /* read in user-defined tree and set it up */
  short nextnode, lparens, rparens;
  boolean *names;
  double x;
  short i;

  nextnode = spp;
  names = (boolean *)Malloc(spp*sizeof(boolean));
  for (i = 0; i < spp; i++)
    names[i] = false;
  lparens = 0;
  rparens = 0;
  addelement(&curtree.root, &nextnode,&lparens,&rparens,names);
  curtree.root->back = NULL;
  findch(';',&rparens);
  if (progress)
    printf("\n\n");
  if (lengths)
    tymetrav(curtree.root, &x);
  fscanf(infile, "%*[^\n]");
  getc(infile);
  free(names);
}  /* treeread */

void coordinates(p, tipy)
node *p;
short *tipy;
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = (*tipy);
    p->ymin   = (*tipy);
    p->ymax   = (*tipy);
    (*tipy)  += down;
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
  p->xcoord = (short)(0.5 - over * p->tyme);
  p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */

void drawline(i, scale)
short i;
double scale;
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q, *r, *first, *last;
  short n, j;
  boolean extra, done;

  p = curtree.root;
  q = curtree.root;
  extra = false;
  if (p->ycoord == i) {
    if (p->number - spp >= 10)
      fprintf(outfile, "-%2hd", p->number - spp);
    else
      fprintf(outfile, "--%hd", p->number - spp);
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
    n = (short)(scale * (p->xcoord - q->xcoord) + 0.5);
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
	if (q->number - spp >= 10)
	  fprintf(outfile, "%2hd", q->number - spp);
	else
	  fprintf(outfile, "-%hd", q->number - spp);
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
      putc(naym[p->number - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */

void printree()
{
  /* prints out diagram of the tree */
 /* Local variables for printree: */
  short tipy;
  double scale;
  short i;
  node *p;

  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  coordinates(curtree.root, &tipy);
  p = curtree.root;
  while (!p->tip)
    p = p->next->back;
  scale = 1.0 / (short)(p->tyme - curtree.root->tyme + 1.000);
  putc('\n', outfile);
  for (i = 1; i <= tipy - down; i++)
    drawline(i, scale);
  putc('\n', outfile);
}  /* printree */

#undef down
#undef over

void describe(p)
node *p;
{
  short i;
  double v;

  if (p == curtree.root)
    fprintf(outfile, " root         ");
  else
    fprintf(outfile, "%4hd          ", p->back->number - spp);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(naym[p->number - 1][i], outfile);
  } else
    fprintf(outfile, "%4hd      ", p->number - spp);
  if (p != curtree.root) {
    fprintf(outfile, "%10.5f", fracchange * (p->tyme - curtree.root->tyme));
    v = fracchange * (p->tyme - curtree.nodep[p->back->number - 1]->tyme);
    fprintf(outfile, "%12.5f", v);
  }
  putc('\n', outfile);
  if (!p->tip) {
    describe(p->next->back);
    describe(p->next->next->back);
  }
}  /* describe */

void summarize()
{
  short i, j, mx;
  double mode, sum;
  double *like,*nulike;
  short **mp;

  like    = (double *)Malloc(categs * sizeof(double));
  nulike  = (double *)Malloc(categs * sizeof(double));
  mp      = (short **)Malloc(sites * sizeof(short *));
  for (i=0; i<=sites-1; i++)
     mp[i] = (short *)Malloc(sizeof(short)*categs);
  fprintf(outfile, "\nLn Likelihood = %11.5f\n\n", curtree.likelihood);
  fprintf(outfile, " Ancestor      Node      Node Time     Length\n");
  fprintf(outfile, " --------      ----      ---- ----     ------\n");
  describe(curtree.root);
  putc('\n', outfile);
  if (ctgry && categs > 1) {
    for (i = 0; i < categs; i++)
      like[i] = 1.0;
    for (i = sites - 1; i >= 0; i--) {
      sum = 0.0;
      for (j = 0; j < categs; j++) {
	nulike[j] = (lambda1 + lambda * probcat[j]) * like[j]; 
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
 "\nCombination of categories that contributes the most to the likelihood:\n");
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
  }
  putc('\n', outfile);
  free(like);
  free(nulike);
  for (i=0;i<sites;++i)
    free(mp[i]);
  free(mp);
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
      if (naym[p->number - 1][i - 1] != ' ')
	n = i;
    }
    for (i = 0; i < n; i++) {
      c = naym[p->number - 1][i];
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
    if (col > 55) {
      putc('\n', treefile);
      col = 0;
    }
    treeout(p->next->next->back);
    putc(')', treefile);
    col++;
  }
  if (p == curtree.root) {
    fprintf(treefile, ";\n");
    return;
  }
  x = fracchange * (p->tyme - curtree.nodep[p->back->number - 1]->tyme);
  if (x > 0.0)
    w = (short)(0.4342944822 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (short)(0.4342944822 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  fprintf(treefile, ":%*.5f", (int)(w + 7), x);
  col += w + 8;
}  /* treeout */


void nodeinit(p)
node *p;
{
  /* set up times at one node */
  node *q, *r;
  double lowertyme;

  q = p->next->back;
  r = p->next->next->back;
  lowertyme = q->tyme;
  if (r->tyme < q->tyme)
    lowertyme = r->tyme;
  p->tyme = lowertyme - 0.1;
  p->next->tyme = p->tyme;
  p->next->next->tyme = p->tyme;
  q->v = q->tyme - p->tyme;
  p->next->v = q->v;
  r->v = r->tyme - p->tyme;
  p->next->next->v = r->v;
}  /* nodeinit */

void initrav(p)
node *p;
{
  /* traverse to set up times throughout tree */
  if (p->tip)
    return;
  initrav(p->next->back);
  initrav(p->next->next->back);
  nodeinit(p);
}  /* initrav */


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
  if (p == curtree.root)
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
  /* evaluate likelihood of tree, after iterating branch lengths */
  short i;
  double dummy;

  for (i = 0; i < spp; i++)
    curtree.nodep[i]->initialized = false;
  for (i = spp; i < nonodes; i++) {
    curtree.nodep[i]->initialized = false;
    curtree.nodep[i]->next->initialized = false;
    curtree.nodep[i]->next->next->initialized = false;
  }
  if (!lengths)
    initrav(curtree.root);
  travsp(curtree.root);
  for (i = 1; i <= smoothings * 4; i++)
    smooth(curtree.root);
  dummy = evaluate(curtree.root);
}  /* treevaluate */


void maketree()
{
  /* constructs a binary tree from the pointers in curtree.nodep,
    adds each node at location which yields highest likelihood
    then rearranges the tree for greatest likelihood */

  short i, j, numtrees, num;
  double bestlike, gotlike, sum, sum2, sd;
  node *item, *nufork, *dummy, *p1, *p2, *p3;
  double TEMP;

  if (!usertree) {
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble) {
      for (i = 0; i < spp; i++) {
	j = (short)(randum(seed) * spp) + 1;
	k = enterorder[j - 1];
	enterorder[j - 1] = enterorder[i];
	enterorder[i] = k;
      }
    }
    curtree.root = curtree.nodep[spp];
    if (progress) {
      printf("\nAdding species:\n");
      printf("   ");
    }
    add(curtree.nodep[enterorder[0] - 1], curtree.nodep[enterorder[1] - 1],
	curtree.nodep[spp]);
    if (progress) {
      for (i = 0; i < nmlngth; i++)
	putchar(naym[enterorder[0] - 1][i]);
      printf("\n   ");
      for (i = 0; i < nmlngth; i++)
	putchar(naym[enterorder[1] - 1][i]);
      putchar('\n');
    }
    lastrearr = false;
    for (i = 3; i <= spp; i++) {
      bestree.likelihood = -999999.0;
      bestyet = -999999.0;
      there = curtree.root;
      item = curtree.nodep[enterorder[i - 1] - 1];
      nufork = curtree.nodep[spp + i - 2];
      addpreorder(curtree.root, item, nufork);
      add(there, item, nufork);
      like = bestyet;
      rearrange(&curtree.root);
      if (curtree.likelihood > bestree.likelihood)
	copy_(&curtree, &bestree);
      if (progress) {
	printf("   ");
	for (j = 0; j < nmlngth; j++)
	  putchar(naym[enterorder[i - 1] - 1][j]);
	putchar('\n');
      }
      lastrearr = (i == spp);
      if (lastrearr && globle) {
	if (progress) {
	  printf("Doing global rearrangements\n");
	  printf("  !");
	  for (j = 1; j <= nonodes; j++)
	    putchar('-');
	  printf("!\n");
	}
	bestlike = bestyet;
	do {
	  if (progress)
	    printf("   ");
	  gotlike = bestlike;
	  for (j = 0; j < nonodes; j++) {
	    bestyet = -999999.00;
	    item = curtree.nodep[j];
	    if (item != curtree.root) {
	      nufork = curtree.nodep[curtree.nodep[j]->back->number - 1];
	      re_move(&item, &nufork);
	      there = curtree.root;
	      addpreorder(curtree.root, item, nufork);
	      add(there, item, nufork);
	    }
	    if (progress)
	      putchar('.');
	  }
	  if (progress)
	    putchar('\n');
	} while (bestlike < gotlike);
      }
    }
    if (njumble > 1 && lastrearr) {
      for (i = 0; i < spp; i++)
        re_move(&curtree.nodep[i], &dummy);
      if (jumb == 1 || bestree2.likelihood < bestree.likelihood)
        copy_(&bestree, &bestree2);
    }
    if (jumb == njumble) {
      if (njumble > 1) copy_(&bestree2, &curtree);
      else copy_(&bestree, &curtree);
      fprintf(outfile, "\n\n");
      curtree.likelihood = evaluate(curtree.root);
      printree();
      summarize();
      if (trout) {
        col = 0;
        treeout(curtree.root);
      }
    }
  } else {
    fscanf(infile, "%hd%*[^\n]", &numtrees);
    l0gl = (double *)Malloc(numtrees * sizeof(double));
    l0gf = (double **)Malloc(numtrees * sizeof(double *));
    for (i=0;i<numtrees;++i)
      l0gf[i] = (double *)Malloc(endsite * sizeof(double));
    getc(infile);
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
	putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }
    fprintf(outfile, "\n\n");
    which = 1;
    while (which <= numtrees) {
      setuptree(&curtree);
      treeread();
      treevaluate();
      printree();
      summarize();
      if (trout) {
	col = 0;
	treeout(curtree.root);
      }
      which++;
    }
    putc('\n', outfile);
    if (!auto_ && numtrees > 1 && weightsum > 1 ) {
      fprintf(outfile, "Tree    Ln L      Diff Ln L     Its S.D.");
      fprintf(outfile, "   Significantly worse?\n\n");
      num = numtrees;
      for (which = 1; which <= num; which++) {
	fprintf(outfile, "%3hd%12.5f", which, l0gl[which - 1]);
	if (maxwhich == which)
	  fprintf(outfile, "  <------ best\n");
	else {
	  sum = 0.0;
	  sum2 = 0.0;
	  for (j = 0; j < endsite; j++) {
	    sum += aliasweight[j] *
		   (l0gf[maxwhich - 1][j] - l0gf[which - 1][j]);
	    TEMP = l0gf[maxwhich - 1][j] - l0gf[which - 1][j];
	    sum2 += aliasweight[j] * (TEMP * TEMP);
	  }
	  sd = sqrt(weightsum / (weightsum - 1.0) * (sum2 - sum * sum / weightsum));
	  fprintf(outfile, "%12.5f%12.4f",
		  l0gl[which - 1] - maxlogl, sd);
	  if (sum > 1.95996 * sd)
	    fprintf(outfile, "           Yes\n");
	  else
	    fprintf(outfile, "            No\n");
	}
      }
      fprintf(outfile, "\n\n");
    }
    free(l0gl);
    for (i=0;i<numtrees;++i)
      free(l0gf[i]);
    free(l0gf);
  }
  if (jumb == njumble) {
    if (progress) {
      printf("\nOutput written to output file\n\n");
      if (trout)
        printf("Tree(s) also written onto file\n\n");
    }
    for (i=0; i<spp;i++){
      for (j=0;j<endsite;++j){
        free(curtree.nodep[i]->x[j]);
        if (!usertree) {
          free(bestree.nodep[i]->x[j]);
          if (njumble > 1)
            free(bestree2.nodep[i]->x[j]);
        }
      }
      free(curtree.nodep[i]->x);
      if (!usertree) {
        free(bestree.nodep[i]->x);
        if (njumble > 1)
	  free(bestree2.nodep[i]->x);
      }
    }
    for (i=spp; i < nonodes ; i++){
      p1 = curtree.nodep[i];
      if (!usertree) {
        p2 = bestree.nodep[i];
        if (njumble > 1)
          p3 = bestree2.nodep[i];
      }
      for (k=0;k<3;++k) {
        for (j=0;j<endsite;++j){
	  free(p1->x[j]);
          if (!usertree) {
	    free(p2->x[j]);
            if (njumble > 1)
              free(p3->x[j]);
          }
        }
        free(p1->x);
        if (!usertree) {
	  free(p2->x);
          if (njumble > 1)
            free(p3->x);
        }
        p1=p1->next;
        if (!usertree) {
	  p2=p2->next;
          if (njumble > 1)
          p3=p3->next;
        }
      }
    }
  }
}  /* maketree */


