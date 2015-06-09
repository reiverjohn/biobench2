#include "drawgraphics.h"

/* Version 3.52c.  Copyright (c) 1986-1993 by Joseph Felsenstein and
  Christopher A. Meacham.  Additional code written by Hisashi Horino,
  Sean Lamont, Andrew Keefe, and Akiko Fuseki.
  Permission is granted to copy, distribute,
  and modify this program provided that (1) this copyright message is
  not removed and (2) no fee is charged for this program. */

#define maxnch          30
#define point           '.'
#define fontsize        3800
#define pi              3.141592653
#define epsilon         0.00001
#define xstart          10
#define ystart          35
#define gap             0.5
#define iterations      5

#ifdef MAC
#undef maxnodes
#define maxnodes 250
#endif

typedef enum {  fixed, radial, along} labelorient;
FILE *treefile, *plotfile;
char        pltfilename[100];
long        ntips, nextnode,  strpwide, strpdeep,
	    strptop, strpbottom,  payge, numlines;
double       xmargin, ymargin, topoflabels, rightoflabels, leftoflabels,
              bottomoflabels, ark, maxx, maxy, minx, miny, scale, xscale,
	      yscale, xoffset, yoffset, charht, xnow, ynow, xunitspercm,
	      yunitspercm, xsize, ysize, xcorner, ycorner,labelheight,
	      labelrotation, treeangle,  expand, bscale,enthusiasm;
boolean        canbeplotted, preview, previewing, dotmatrix,haslengths,
	       uselengths, regular, didreroot, rotate, empty, rescaled,
               notfirst, improve;
       double textlength[maxnodes], firstlet[maxnodes];
       striptype stripe;
       plottertype plotter, oldplotter, previewer;
       growth grows;
       labelorient labeldirec;
       node *root, *where;
      node *nodep[maxnodes];
       fonttype font;
       enum {  yes, no } penchange, oldpenchange;
char ch;
char fontname[64];
long filesize;
long strpdiv;

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
  /* make ch upper-case */
   *ch = (islower(*ch) ?  toupper(*ch) : (*ch));
}  /* uppercase */

/* Local variables for treeread: */

void getch(c)
Char *c;
{
  /* get next nonblank character */
  do {*c = getc(treefile);
  } while ((*c == ' ')||(*c == '\n')||(*c == '\t'));
}  /* getch */

void processlength(p)
node *p;
{
  long  digit, ordzero;
  double valyew, divisor;
  boolean pointread, minusread;

  ordzero = '0';
  pointread = false;
  minusread = false;
  valyew = 0.0;
  divisor = 1.0;
  getch(&ch);
  digit = ch - ordzero;
  while (((unsigned long)digit <= 9) || (ch == '.') || (ch == '-')){
        if (ch == '.')   pointread = true;
   else if (ch == '-')   minusread = true;
   else {
        valyew = valyew * 10.0 + digit;
        if (pointread)
          divisor *= 10.0;
      }
        getch(&ch);
        digit = ch - ordzero;
      }
 if (!minusread)
   p->oldlen = valyew / divisor;
 else
   p->oldlen = 0.0;
  /* processlength */
}

void addelement(p, q)
node **p, *q;
{
  /* read in and add next part of tree, it will be node p
    and will be hooked to pointer q */
  node *pp;
  long n,indx;
  boolean notlast;

  nextnode++;
  *p = (node *)Malloc((long)sizeof(node));
  nodep[nextnode - 1] = *p;
  indx = nextnode;
  (*p)->index = indx;
  if (ch == '(') {
    (*p)->tip = false;
    pp = *p;
    notlast = true;
    while (notlast) {
      pp->next = (node *)Malloc((long)sizeof(node));
      pp->next->tip = false;
      nextnode++;
      nodep[nextnode - 1] = pp->next;
      pp->next->index = indx;
      pp = pp->next;
      getch(&ch);
      addelement(&pp->back, pp);
      if (ch == ')') {
	notlast = false;
	do {
	  getch(&ch);
	} while (ch != ':' && ch != ',' && ch != ')' && ch != '[' &&
		 ch != ';');
      }
    }
    pp->next = *p;
  } else {
    (*p)->tip = true;
    ntips++;
    n = 1;
    do {
      if (!ebcdic && (ch & 255) == 255)
	ch = '\'';
      if (!ebcdic && (ch & 255) > 175)
	ch -= 48;
      if (!ebcdic && (ch & (~127)) != 0)
	ch -= 64;
      if (ch == '_')
	ch = ' ';
      if (n < maxnch)
	(*p)->nayme[n - 1] = ch;
      if (eoln(treefile)) {
	fscanf(treefile, "%*[^\n]");
	(void)getc(treefile);
      }
      ch = getc(treefile);
      if (ch == '\n')
	ch = ' ';
      n++;
    } while (ch != ':' && ch != ',' && ch != ')');
    if (n > maxnch)
      n = maxnch + 1;
    (*p)->naymlength = n - 1;
  }
  if (ch == ':')
    processlength(*p);
  else
    haslengths = (haslengths && q == NULL);
  (*p)->back = q;
  if (haslengths && q != NULL)
    (*p)->back->oldlen = (*p)->oldlen;
  (*p)->r = 0.0;
  (*p)->theta = 0.0;
}  /* addelement */

void treeread()
{
  /* read a tree from the treefile and set up nodes and pointers */
  haslengths = true;
  ntips = 0;
  nextnode = 0;
  getch(&ch);
  addelement(&root, NULL);
  root->oldlen = 0.0;
  fscanf(treefile, "%*[^\n]");
  (void)getc(treefile);
  uselengths = haslengths;
  where = root;
  rotate = true;
}  /* treeread */

void initialparms()
{
  /* initialize parameters */

  getplotter();
  plotrparms();
  xmargin = 0.08 * xsize;
  ymargin = 0.08 * ysize;
  xscale = xunitspercm;
  yscale = yunitspercm;
  grows = vertical;
  treeangle = pi / 2.0;
  ark = 2 * pi;
  improve = true;
  regular = false;
  rescaled = true;
  bscale = 1.0;
  labeldirec = fixed;
  labelrotation = 0.0;
  charht = 0.3333;
  numlines = dotmatrix ? ((long)floor(yunitspercm * ysize + 0.5) / strpdeep):1;
  enthusiasm = 1.0/(double)ntips;
}  /* initialparms */


long showparms()
{
  long i;
  long numtochange;
  Char ch,input[64];
  double treea;

  putchar('\n');
  if (previewer == tek)
    printf("%c\f", escape);
  else {
    for (i = 1; i <= 24; i++)
      putchar('\n');
  }
  printf("Here are the settings: \n\n");
  printf(" (1)        Use branch lengths:  ");
  if (haslengths)
    printf("%s\n",uselengths ? "Yes" : "No");
   else
    printf("(no branch lengths available)\n");
  printf(" (2)           Angle of labels:");
  if (labeldirec == fixed) {
    printf("  Fixed angle of");
    if (labelrotation >= 10.0)
      printf("%6.1f", labelrotation);
    else if (labelrotation <= -10.0)
      printf("%7.1f", labelrotation);
    else if (labelrotation < 0.0)
      printf("%6.1f", labelrotation);
    else
      printf("%5.1f", labelrotation);
    printf(" degrees\n");
  } else if (labeldirec == radial)
    printf("  Radial\n");
  else
    printf("  Along branches\n");
  printf(" (3)          Rotation of tree:");
  treea = treeangle * 180 / pi;
  if (treea >= 100.0)
    printf("%7.1f\n", treea);
  else if (treea >= 10.0)
    printf("%6.1f\n", treea);
  else if (treea <= -100.0)
    printf("%8.1f\n", treea);
  else if (treea <= -10.0)
    printf("%7.1f\n", treea);
  else if (treea < 0.0)
    printf("%6.1f\n", treea);
  else
    printf("%5.1f\n", treea);
  printf(" (4)     Angle of arc for tree:");
  treea = 180 * ark / pi;
  if (treea >= 100.0)
    printf("%7.1f\n", treea);
  else if (treea >= 10.0)
    printf("%6.1f\n", treea);
  else if (treea <= -100.0)
    printf("%8.1f\n", treea);
  else if (treea <= -10.0)
    printf("%7.1f\n", treea);
  else if (treea < 0.0)
    printf("%6.1f\n", treea);
  else
    printf("%5.1f\n", treea);
  printf(" (5)   Iterate to improve tree:  %s\n",
	 (improve ? "Yes" : "No"));
  printf(" (6)    Scale of branch length:");
  if (rescaled)
    printf("  Automatically rescaled\n");
  else
    printf("  Fixed:%6.2f cm per unit branch length\n", bscale);
  didreroot = false;
  printf(" (7)        Horizontal margins:%6.2f cm\n", xmargin);
  printf(" (7)          Vertical margins:%6.2f cm\n", ymargin);
  printf(" (8) Relative character height:%8.4f\n", charht);
  if (!improve)
    printf(" (9)     Regularize the angles:  %s\n",(regular ? "Yes" : "No"));
  else
    printf(" (9)       Enthusiasm constant:%9.5f\n",enthusiasm);
  if (plotter == lw)
    printf(" (10)           Font          :  %s\n",fontname);

  printf("\n\n Do you want to accept these? (Yes or No)\n");
  for (;;) {
    printf(" Type Y or N or the number (1-%2ld) of the one to change: \n",
            (plotter == lw) ? 10L : 9L);
    gets(input);
    uppercase(&input[0]);
    ch=input[0];
    numtochange = atoi(input);
    if ((ch == 'Y' || ch == 'N') || (numtochange >= 1 && numtochange <=
       ((plotter == lw) ? 10 : 9)))
      break;
  }
  return (ch == 'Y') ? -1 : numtochange;
}  /* showparms */

void rerootit(p)
node *p;
{
  node *q;

  q = root;
  while (q->next != root)
    q = q->next;
  q->next = root->next;
  q = p;
  while (q->next != p)
    q = q->next;
  root->next = p;
  q->next = root;
}  /* rerootit */

void reroot()
{
  /* reroot the tree, traversing tree */
  where = root->next;
  if (!rotate)
    rotate = where->back->tip;
  if (rotate) {
    where = where->next;
    rotate = false;
  } else {
    where = where->back;
    rotate = true;
  }
  rerootit(where);
}  /* reroot */


void getparms(numtochange)
long numtochange;
{
  /* get from user the relevant parameters for the plotter and diagram */
  Char ch;
  boolean ok;

  if (numtochange == 0) {
    do {
      printf(" Type the number of one that you want to change (1-%2ld):\n",
              (plotter == lw) ? 10L : 9L);
      scanf("%ld%*[^\n]", &numtochange);
      (void)getchar();
    } while (numtochange < 1 || numtochange > ((plotter == lw) ? 10 : 9));
  }
  switch (numtochange) {

  case 1:
    if (haslengths)
      uselengths = !uselengths;
    else {
      printf("Cannot use lengths since not all of them exist\n");
      uselengths = false;
    }
    break;

  case 2:
    printf("\nDo you want labels to be Fixed angle, Radial,");
    printf(" or Along branches?\n");
    do {
      printf(" Type F, R, or A\n");
      scanf("%c%*[^\n]", &ch);
      (void)getchar();
      if (ch == '\n')
	ch = ' ';
      uppercase(&ch);
    } while (ch != 'F' && ch != 'R' && ch != 'A');
    switch (ch) {

    case 'A':
      labeldirec = along;
      break;

    case 'F':
      labeldirec = fixed;
      break;

    case 'R':
      labeldirec = radial;
      break;
    }
    if (labeldirec == fixed) {
      printf("Are the labels to be plotted vertically (90),\n");
      printf(" horizontally (0), or downwards (-90) ?\n");
      do {
	printf(" Choose an angle in degrees from 90 to -90: \n");
	scanf("%lf%*[^\n]", &labelrotation);
	(void)getchar();
      } while ((labelrotation < -90.0 || labelrotation > 90.0) &&
	       labelrotation != -99.0);
    }
    break;

  case 3:
    printf("\n At what angle is the tree to be plotted?\n");
    do {
      printf(" Choose an angle in degrees from 360 to -360: \n");
      scanf("%lf%*[^\n]", &treeangle);
      (void)getchar();
      uppercase(&ch);
    } while (treeangle < -360.0 && treeangle > 360.0);
    treeangle = treeangle * pi / 180;
    break;

  case 4:
    printf(" How many degrees (up to 360) of arc\n");
    printf("  should the tree occupy? (Currently it is %5.1f)\n",
	   180 * ark / pi);
    do {
      printf("Enter a number of degrees from 0 up to 360)\n");
      scanf("%lf%*[^\n]", &ark);
      (void)getchar();
    } while (ark <= 0.0 || ark > 360.0);
    ark = ark * pi / 180;
    break;

  case 5:
    improve = !improve;
    break;

  case 6:
    rescaled = !rescaled;
    if (!rescaled) {
      printf("Centimeters per unit branch length?\n");
      scanf("%lf%*[^\n]", &bscale);
      (void)getchar();
    }
    break;

  case 7:
    printf("\nThe tree will be drawn to fit in a rectangle which has \n");
    printf(" margins in the horizontal and vertical directions of:\n");
    printf("%6.2f cm (horizontal margin) and%6.2f cm (vertical margin)\n\n",
	   xmargin, ymargin);
    do {
      printf(" New value (in cm) of horizontal margin?\n");
      scanf("%lf%*[^\n]", &xmargin);
      (void)getchar();
      ok = ((unsigned)xmargin < xsize / 2.0);
      if (!ok)
	printf(" Impossible value.  Please retype it.\n");
    } while (!ok);
    do {
      printf(" New value (in cm) of vertical margin?\n");
      scanf("%lf%*[^\n]", &ymargin);
      (void)getchar();
      ok = ((unsigned)ymargin < ysize / 2.0);
      if (!ok)
	printf(" Impossible value.  Please retype it.\n");
    } while (!ok);
    break;

  case 8:
    printf("New value of character height?\n");
    scanf("%lf%*[^\n]", &charht);
    (void)getchar();
    break;
  case 9:
    if (improve) {
      do {
	printf("Enthusiasm constant? (must be between 0 and 1)\n");
	scanf("%lf%*[^\n]", &enthusiasm);
	getchar();
      } while (enthusiasm <= 0.0 || enthusiasm > 1.0);}
    else
      regular = !regular;
    break;
  case 10:
    printf("Enter font name or \"Hershey\" for default font\n");
    gets(fontname);
    break;
  }
}  /* getparms */


Local void getwidth(p)
node *p;
{
  /* get width and depth beyond each node */
  double nw, nd;
  node *pp, *qq;

  nd = 0.0;
  if (p->tip)
    nw = 1.0;
  else {
    nw = 0.0;
    qq = p;
    pp = p->next;
    do {
      getwidth(pp->back);
      nw += pp->back->width;
      if (pp->back->depth > nd)
	nd = pp->back->depth;
      pp = pp->next;
    } while (pp != qq);
  }
  p->depth = nd + p->length;
  p->width = nw;
}  /* getwidth */

void plrtrans(p, theta, lower, upper)
node *p;
double theta, lower, upper;
{
  /* polar coordinates of a node relative to start */
  long num;
  double nn, pr, ptheta, angle, angle2, subangle, len;
  node *pp, *qq;

  nn = p->width;
  angle = theta;
  subangle = (upper - lower) / nn;
  qq = p;
  pp = p->next;
  if (p->tip)
    return;
  angle = upper;
  do {
    angle -= pp->back->width / 2.0 * subangle;
    pr = p->r;
    ptheta = p->theta;
    if (regular) {
      num = 1;
      while (num * subangle < 2 * pi)
	num *= 2;
      if (angle >= 0.0)
	angle2 = 2 * pi / num * (long)(num * angle / (2 * pi) + 0.5);
      else
	angle2 = 2 * pi / num * (long)(num * angle / (2 * pi) - 0.5);
    } else
      angle2 = angle;
    if (uselengths)
      len = pp->back->oldlen;
    else
      len = 1.0;
    pp->back->r = sqrt(len * len + pr * pr + 2 * len * pr * cos(angle2 - ptheta));
    if (fabs(pr * cos(ptheta) + len * cos(angle2)) > epsilon)
      pp->back->theta = atan((pr * sin(ptheta) + len * sin(angle2)) /
			     (pr * cos(ptheta) + len * cos(angle2)));
    else if (pr * sin(ptheta) + len * sin(angle2) >= 0.0)
      pp->back->theta = pi / 2;
    else
      pp->back->theta = 1.5 * pi;
    if (pr * cos(ptheta) + len * cos(angle2) < -epsilon)
      pp->back->theta += pi;
    if (!pp->back->tip)
      plrtrans(pp->back, pp->back->theta,
		 angle - pp->back->width * subangle / 2.0,
		 angle + pp->back->width * subangle / 2.0);
    else
      pp->back->oldtheta = angle2;
    angle -= pp->back->width / 2.0 * subangle;
    pp = pp->next;
  } while (pp != qq);
}  /* plrtrans */

void coordtrav(p, xx, yy)
node *p;
double *xx, *yy;
{
  /* compute x and y coordinates */
  long i;
  node *pp;

  if (!p->tip) {
    pp = p->next;
    while (pp != p) {
      coordtrav(pp->back, xx,yy);
      pp = pp->next;
    }
  }
  if (p->tip) {
    i = 1;
    while (nodep[i - 1] != p)
      i++;
    textlength[i - 1] = (double)lengthtext(p->nayme, p->naymlength, font);
  }
  (*xx) = p->r * cos(p->theta);
  (*yy) = p->r * sin(p->theta);
  if ((*xx) > maxx)
    maxx = (*xx);
  if ((*xx) < minx)
    minx = (*xx);
  if ((*yy) > maxy)
    maxy = (*yy);
  if ((*yy) < miny)
    miny = (*yy);
  p->xcoord = (*xx);
  p->ycoord = (*yy);
}  /* coordtrav */


double angleof(x, y)
double x, y;
{
  /* compute the angle of a vector */
  double theta;

  if (fabs(x) > epsilon)
    theta = atan(y / x);
  else if (y >= 0.0)
    theta = pi / 2;
  else
    theta = 1.5 * pi;
  if (x < -epsilon)
    theta = pi + theta;
  if (theta > 2 * pi)
    theta -= 2 * pi;
  return theta;
}  /* angleof */

void polartrav(p,xx,yy,leftx,lefty,rightx,righty)
node    *p;
double  *xx,*yy,*leftx,*lefty,*rightx,*righty;
{
  /* go through subtree getting left and right vectors */
  double x, y;
  boolean lookatit;
  node *pp;

  lookatit = true;
  if (!p->tip)
    lookatit = (p->next->next->next != p || p->index != root->index);
  if (lookatit) {
    x = nodep[p->index - 1]->xcoord;
    y = nodep[p->index - 1]->ycoord;
    if ((y - (*yy)) * (*rightx) - (x - (*xx)) * (*righty) < 0.0) {
      (*rightx) = x - (*xx);
      (*righty) = y - (*yy);
    }
    if ((y - (*yy)) * (*leftx) - (x - (*xx)) * (*lefty) > 0.0 &&
	!(notfirst && fabs((*rightx) - x + (*xx)) +
	  fabs((*righty) - y + (*yy)) < epsilon)) {
      (*leftx) = x - (*xx);
      (*lefty) = y - (*yy);
      notfirst = true;
    }
  }
  if (p->tip)
    return;
  pp = p->next;
  while (pp != p) {
    if (pp != root)
      polartrav(pp->back,xx,yy,leftx,lefty,rightx,righty);
    pp = pp->next;
  }
}  /* polartrav */

void tilttrav(q,xx,yy,sinphi,cosphi)
node *q;
double *xx,*yy,*sinphi,*cosphi;
{
  /* traverse to move successive nodes */
  double x, y;
  node *pp;

  pp = nodep[q->index - 1];
  x = pp->xcoord;
  y = pp->ycoord;
  pp->xcoord = (*xx) + (x - (*xx)) * (*cosphi) + (y - (*yy)) * (*sinphi);
  pp->ycoord = (*yy) + ((*xx) - x) * (*sinphi) + (y - (*yy)) * (*cosphi);
  if (q->tip)
    return;
  pp = q->next;
  while (pp != q) {
    if (pp != root)
      tilttrav(pp->back,xx,yy,sinphi,cosphi);
    pp = pp->next;
  }
}  /* tilttrav */

void polarize(p,xx,yy)
node *p;
double *xx,*yy;
{
  double TEMP, TEMP1;

  if (fabs(p->xcoord - (*xx)) > epsilon)
    p->oldtheta = atan((p->ycoord - (*yy)) / (p->xcoord - (*xx)));
  else if (p->ycoord - (*yy) >= 0.0)
    p->oldtheta = pi / 2;
  else
    p->oldtheta = 1.5 * pi;
  if (p->xcoord - (*xx) < -epsilon)
    p->oldtheta += pi;
  if (fabs(p->xcoord - root->xcoord) > epsilon)
    p->theta = atan((p->ycoord - root->ycoord) / (p->xcoord - root->xcoord));
  else if (p->ycoord - root->ycoord >= 0.0)
    p->theta = pi / 2;
  else
    p->theta = 1.5 * pi;
  if (p->xcoord - root->xcoord < -epsilon)
    p->theta += pi;
  TEMP = p->xcoord - root->xcoord;
  TEMP1 = p->ycoord - root->ycoord;
  p->r = sqrt(TEMP * TEMP + TEMP1 * TEMP1);
}  /* polarize */

void improvtrav(p)
node *p;
{
  /* traverse tree trying different tiltings at each node */
  double xx, yy, leftx, lefty, rightx, righty, cosphi, sinphi;
  long n;
  double usedangle, langle, rangle, freeangle, meanangle, sumrot;
  node *pp, *qq;

  if (p->tip)
    return;
  xx = p->xcoord;
  yy = p->ycoord;
  if (p != root) {
    n = 0;
    usedangle = 0.0;
    pp = p->next;
    do {
      leftx = pp->back->xcoord - xx;
      lefty = pp->back->ycoord - yy;
      rightx = leftx;
      righty = lefty;
      notfirst = false;
      if (!pp->back->tip)
	polartrav(pp->back, &xx,&yy,&leftx,&lefty,&rightx,&righty);
      n++;
      langle = angleof(leftx, lefty);
      rangle = angleof(rightx, righty);
      if (rangle > langle)
	langle += 2 * pi;
      pp->lefttheta = langle;
      pp->righttheta = rangle;
      usedangle += langle - rangle;
      pp = pp->next;
    } while (pp != p->next);
    freeangle = 2 * pi - usedangle;
    meanangle = freeangle / n;
    sumrot = 0.0;
    qq = p;
    pp = p->next;
    while (pp != p) {
      langle = qq->righttheta;
      rangle = pp->lefttheta;
      if (rangle > langle)
	langle += 2 * pi;
      sumrot += enthusiasm * (meanangle - langle + rangle);
      cosphi = cos(sumrot);
      sinphi = sin(sumrot);
      if (pp != root)
	tilttrav(pp->back, &xx,&yy,&sinphi,&cosphi);
      qq = pp;
      pp = pp->next;
    }
  }
  pp = p->next;
  while (pp != p) {
    if (pp != root)
      polarize(pp->back, &xx,&yy);
    pp = pp->next;
  }
  pp = p->next;
  while (pp != p) {
    if (pp != root)
      improvtrav(pp->back);
    pp = pp->next;
  }
}  /* improvtrav */

void coordimprov(xx,yy)
double *xx,*yy;
{
  /* use angles calculation to improve node coordinate placement */
  long i;
  for (i = 1; i <= iterations; i++)
    improvtrav(root);
}  /* coordimprov */


void calculate()
{
  /* compute coordinates for tree */
  double xx, yy;
  long i;
  double nttot, fontheight, labangle, top, bot, rig, lef;

  for (i = 0; i < nextnode; i++)
    nodep[i]->width = 1.0;
  for (i = 0; i < nextnode; i++)
    nodep[i]->xcoord = 0.0;
  for (i = 0; i < nextnode; i++)
    nodep[i]->ycoord = 0.0;
  if (!uselengths) {
    for (i = 0; i < nextnode; i++)
      nodep[i]->length = 1.0;
  } else {
    for (i = 0; i < nextnode; i++)
      nodep[i]->length = nodep[i]->oldlen;
  }
  getwidth(root);
  nttot = root->width;
  for (i = 0; i < nextnode; i++)
    nodep[i]->width = nodep[i]->width * ntips / nttot;
  plrtrans(root, treeangle, treeangle - ark / 2.0, treeangle + ark / 2.0);
  maxx = 0.0;
  minx = 0.0;
  maxy = 0.0;
  miny = 0.0;
  coordtrav(root, &xx,&yy);
  if (improve) {
    coordimprov(&xx,&yy);
    coordtrav(root, &xx,&yy);
  }
  fontheight = font[2];
  if (labeldirec == fixed)
    labangle = pi * labelrotation / 180.0;
  for (i = 0; i < nextnode; i++) {
    if (nodep[i]->tip)
      textlength[i] /= fontheight;
  }
  if (ntips > 1)
    labelheight = charht * (maxx - minx) / (ntips - 1);
  else
    labelheight = charht * (maxx - minx);
  topoflabels = 0.0;
  bottomoflabels = 0.0;
  rightoflabels = 0.0;
  leftoflabels = 0.0;
  for (i = 0; i < nextnode; i++) {
    if (nodep[i]->tip) {
      if (labeldirec == radial)
	labangle = nodep[i]->theta;
      else if (labeldirec == along)
	labangle = nodep[i]->oldtheta;
      if (cos(labangle) < 0.0 && labeldirec != fixed)
	labangle -= pi;
      firstlet[i] = (double)lengthtext(nodep[i]->nayme,1L,font)/ fontheight;
      top = (nodep[i]->ycoord - maxy) / labelheight + sin(nodep[i]->oldtheta);
      rig = (nodep[i]->xcoord - maxx) / labelheight + cos(nodep[i]->oldtheta);
      bot = (miny - nodep[i]->ycoord) / labelheight - sin(nodep[i]->oldtheta);
      lef = (minx - nodep[i]->xcoord) / labelheight - cos(nodep[i]->oldtheta);
      if (cos(labangle) * cos(nodep[i]->oldtheta) +
	  sin(labangle) * sin(nodep[i]->oldtheta) > 0.0) {
	if (sin(labangle) > 0.0)
	  top += sin(labangle) * textlength[i];
	top += sin(labangle - 1.25 * pi) * gap * firstlet[i];
	if (sin(labangle) < 0.0)
	  bot -= sin(labangle) * textlength[i];
	bot -= sin(labangle - 0.75 * pi) * gap * firstlet[i];
	if (sin(labangle) > 0.0)
	  rig += cos(labangle - 0.75 * pi) * gap * firstlet[i];
	else
	  rig += cos(labangle - 1.25 * pi) * gap * firstlet[i];
	rig += cos(labangle) * textlength[i];
	if (sin(labangle) > 0.0)
	  lef -= cos(labangle - 1.25 * pi) * gap * firstlet[i];
	else
	  lef -= cos(labangle - 0.75 * pi) * gap * firstlet[i];
      } else {
	if (sin(labangle) < 0.0)
	  top -= sin(labangle) * textlength[i];
	top += sin(labangle + 0.25 * pi) * gap * firstlet[i];
	if (sin(labangle) > 0.0)
	  bot += sin(labangle) * textlength[i];
	bot -= sin(labangle - 0.25 * pi) * gap * firstlet[i];
	if (sin(labangle) > 0.0)
	  rig += cos(labangle - 0.25 * pi) * gap * firstlet[i];
	else
	  rig += cos(labangle + 0.25 * pi) * gap * firstlet[i];
	if (sin(labangle) < 0.0)
	  rig += cos(labangle) * textlength[i];
	if (sin(labangle) > 0.0)
	  lef -= cos(labangle + 0.25 * pi) * gap * firstlet[i];
	else
	  lef -= cos(labangle - 0.25 * pi) * gap * firstlet[i];
	lef += cos(labangle) * textlength[i];
      }
      if (top > topoflabels)
	topoflabels = top;
      if (bot > bottomoflabels)
	bottomoflabels = bot;
      if (rig > rightoflabels)
	rightoflabels = rig;
      if (lef > leftoflabels)
	leftoflabels = lef;
    }
  }
  topoflabels *= labelheight;
  bottomoflabels *= labelheight;
  leftoflabels *= labelheight;
  rightoflabels *= labelheight;
}  /* calculate */


void rescale()
{
  /* compute coordinates of tree for plot or preview device */
  long i;
  double treeheight, treewidth, extrax, extray, temp;

  treeheight = maxy - miny + topoflabels + bottomoflabels;
  treewidth = maxx - minx + rightoflabels + leftoflabels;
  if (grows == vertical) {
    if (!rescaled)
      expand = bscale;
    else {
      expand = (xsize - 2 * xmargin) / treewidth;
      if ((ysize - 2 * ymargin) / treeheight < expand)
	expand = (ysize - 2 * ymargin) / treeheight;
    }
    extrax = (xsize - 2 * xmargin - treewidth * expand) / 2.0;
    extray = (ysize - 2 * ymargin - treeheight * expand) / 2.0;
  } else {
    if (!rescaled)
      expand = bscale;
    else {
      expand = (ysize - 2 * ymargin) / treewidth;
      if ((xsize - 2 * xmargin) / treeheight < expand)
	expand = (xsize - 2 * xmargin) / treeheight;
    }
    extrax = (xsize - 2 * xmargin - treeheight * expand) / 2.0;
    extray = (ysize - 2 * ymargin - treewidth * expand) / 2.0;
  }
  for (i = 0; i < (nextnode); i++) {
    nodep[i]->xcoord = expand * (nodep[i]->xcoord - minx + leftoflabels);
    nodep[i]->ycoord = expand * (nodep[i]->ycoord - miny + bottomoflabels);
    if (grows == horizontal) {
      temp = nodep[i]->ycoord;
      nodep[i]->ycoord = expand * treewidth - nodep[i]->xcoord;
      nodep[i]->xcoord = temp;
    }
    nodep[i]->xcoord += xmargin + extrax;
    nodep[i]->ycoord += ymargin + extray;
  }
}  /* rescale */

void plottree(p, q)
node *p, *q;
{
  /* plot part or all of tree on the plotting device */
  double x1, y1, x2, y2;
  node *pp;

  x2 = xscale * (xoffset + p->xcoord);
  y2 = yscale * (yoffset + p->ycoord);
  if (p != root) {
    x1 = xscale * (xoffset + q->xcoord);
    y1 = yscale * (yoffset + q->ycoord);
    plot(penup, x1, y1);
    plot(pendown, x2, y2);
  }
  if (p->tip)
    return;
  pp = p->next;
  while (pp != p) {
    plottree(pp->back, p);
    pp = pp->next;
  }
}  /* plottree */


void plotlabels(fontname)
char *fontname;
{
  long i;
  double compr, dx, dy, labangle;
  boolean left, right;
  node *lp;

  compr = xunitspercm / yunitspercm;
  if (penchange == yes)
    changepen(labelpen);
  for (i = 0; i < (nextnode); i++) {
    if (nodep[i]->tip) {
      lp = nodep[i];
      labangle = labelrotation * pi / 180.0;
      if (labeldirec == radial)
	labangle = nodep[i]->theta;
      else if (labeldirec == along)
	labangle = nodep[i]->oldtheta;
      if (cos(labangle) < 0.0)
	labangle -= pi;
      right = (cos(labangle) * cos(nodep[i]->oldtheta) +
	       sin(labangle) * sin(nodep[i]->oldtheta) > 0.0);
      left = !right;
      dx = labelheight * expand * cos(nodep[i]->oldtheta);
      dy = labelheight * expand * sin(nodep[i]->oldtheta);
      if (right) {
	dx += labelheight * expand * 0.5 * firstlet[i] * cos(labangle - 0.75 * pi);
	dy += labelheight * expand * 0.5 * firstlet[i] * sin(labangle - 0.75 * pi);
      }
      if (left) {
	dx += labelheight * expand * 0.5 * firstlet[i] * cos(labangle - 0.25 * pi);
	dy += labelheight * expand * 0.5 * firstlet[i] * sin(labangle - 0.25 * pi);
	dx -= textlength[i] * labelheight * expand * cos(labangle);
	dy -= textlength[i] * labelheight * expand * sin(labangle);
      }
	
      plottext(lp->nayme, lp->naymlength,
	       labelheight * expand * xscale / compr, compr,
	       xscale * (lp->xcoord + dx + xoffset),
	       yscale * (lp->ycoord + dy + yoffset), -180 * labangle / pi,
	       font,fontname);
    }
  }
  if (penchange == yes)
    changepen(treepen);
}  /* plotlabels */



main(argc, argv)
     int argc;
     Char *argv[];
{
  long i,n,stripedepth;
  char filename1[100];
#ifdef MAC
  OSErr retcode;
  FInfo  fndrinfo;
  macsetup("Drawtree","Preview");
  argv[0] = "Drawtree";
#endif
#ifdef TURBOC
  if ((registerbgidriver(EGAVGA_driver) <0) ||
      (registerbgidriver(Herc_driver) <0)   ||
      (registerbgidriver(CGA_driver) <0)){
    fprintf(stderr,"Graphics error: %s ",grapherrormsg(graphresult()));
    exit(-1);}
#endif
  strcpy(fontname,"Hershey");

  openfile(&plotfile,PLOTFILE,"w",argv[0],pltfilename);
  openfile(&treefile,TREEFILE,"r",argv[0],NULL);

  printf("DRAWTREE from PHYLIP version %s\n",VERSION);
  printf("Reading tree ... \n");
  treeread();
  printf("Tree has been read.\n");
  printf("Loading the font ... \n");
  loadfont(font,argv[0]);
  printf("Font loaded.\n");
  previewing = false;
  initialparms();

  canbeplotted = false;
  while (!canbeplotted) {
    do {
      n=showparms();
      if ( n != -1)
	getparms(n);
    } while (n != -1);
    calculate();
    rescale();
    canbeplotted = true;
    if (preview)
      canbeplotted=plotpreview(fontname,&xoffset,&yoffset,&scale,ntips,root);
  }
  if (dotmatrix) {
     stripedepth = allocstripe(stripe,(strpwide/8),
			    ((long)(yunitspercm * ysize)));
     strpdeep = stripedepth;
     strpdiv  = stripedepth;
     }
  previewing = false;
  initplotter(ntips,fontname);
  numlines = dotmatrix ? ((long)floor(yunitspercm * ysize + 0.5)/strpdeep):1;   if (plotter != ibmpc)
    printf("Writing plot file ...\n");
  drawit(fontname,&xoffset,&yoffset,numlines,root);
  finishplotter();
  printf("Finished.\n");
  FClose(treefile);
  FClose(plotfile);
  printf("End of run.\n");
#ifdef MAC
  if (plotter == pict){
    strcpy(filename1,pltfilename);
    retcode=GetFInfo(CtoPstr(filename1),0,&fndrinfo);
    fndrinfo.fdType='PICT';
    fndrinfo.fdCreator='MDPL';
    strcpy(filename1,pltfilename);
    retcode=SetFInfo(CtoPstr(PLOTFILE),0,&fndrinfo);}
  if (plotter == lw){
    retcode=GetFInfo(CtoPstr(PLOTFILE),0,&fndrinfo);
    fndrinfo.fdType='TEXT';
    fndrinfo.fdCreator='gsVR';
    retcode=SetFInfo(CtoPstr(PLOTFILE),0,&fndrinfo);}
#endif
  exit(0);
}

int eof(f)
FILE *f;
{
    register long ch;

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
    register long ch;

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
