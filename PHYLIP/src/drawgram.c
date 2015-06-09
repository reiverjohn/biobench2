#include "drawgraphics.h"

/* Version 3.52c.  Copyright (c) 1986-1993 by Joseph Felsenstein and
  Christopher A. Meacham.  Additional code written by Hisashi Horino,
  Sean Lamont, Andrew Keefe, and Akiko Fuseki.
  Permission is granted to copy, distribute, and modify this
  program provided that (1) this copyright message is not removed
  and (2) no fee is charged for this program. */

FILE *treefile,  *plotfile;
char pltfilename[100];
long     ntips, nextnode,  strpwide, strpdeep,strpdiv,
        strptop, strpbottom, payge, numlines;
boolean  preview, previewing, dotmatrix,
         haslengths, uselengths, empty, rescaled;
double xmargin, ymargin, topoflabels, rightoflabels, leftoflabels,
       tipspacing,maxheight, scale, xscale, yscale, xoffset, yoffset,
       nodespace, stemlength, treedepth, xnow, ynow, xunitspercm, yunitspercm,
       xsize, ysize, xcorner, ycorner, labelheight,labelrotation,expand, rooty,
       bscale;
       striptype stripe;
       plottertype plotter, oldplotter, previewer;
       growth grows;
       treestyle style;
       node *root;
node *nodep[maxnodes];
       fonttype font;
       long filesize;
Char   ch;
char   fontname[64];

       enum {  yes, no} penchange,oldpenchange;
static enum {  weighted, intermediate, centered, inner, vshaped} nodeposition;


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
{  /* make ch upper-case */
   *ch = (islower(*ch) ?  toupper(*ch) : (*ch));
}  /* uppercase */

void getch(c)
Char *c;
{  /* get next nonblank character */
  do {
    *c=getc(treefile); }
  while ((*c == ' ')||(*c == '\n')||(*c == '\t'));
  }  /* getch */

Void processlength(p)
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
  while (((unsigned long)digit <= 9) | (ch == '.') || (ch == '-')){
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
  node *pfirst;
  long n;
  boolean notlast;

  nextnode++;
  *p = (node *)Malloc((long)sizeof(node));
  nodep[nextnode - 1] = *p;
  if (ch == '(') {
    (*p)->tip = false;
    (*p)->tipsabove = 0;
    pfirst = *p;
    notlast = true;
    while (notlast) {
      (*p)->next = (node *)Malloc((long)sizeof(node));
      *p = (*p)->next;
      (*p)->tip = false;
      getch(&ch);
      addelement(&(*p)->back, *p);
      pfirst->tipsabove += (*p)->back->tipsabove;
      if (ch == ')') {
        notlast = false;
        do {
          getch(&ch);
        } while (ch != ':' && ch != ',' && ch != ')' && ch != '[' && ch != ';');
      }
    }
    (*p)->next = pfirst;
    *p = pfirst;
  } else {
    (*p)->tip = true;
    (*p)->tipsabove = 1;
    ntips++;
    n = 1;
    do {
      if (ch == '_')
        ch = ' ';
      if (!ebcdic && (ch & 255) == 255)
        ch = '\'';
      if (!ebcdic && (ch & 255) > 175)
        ch -= 48;
      if (!ebcdic && (ch & (~127)) != 0)
        ch -= 64;
      if (n < maxnch)
        (*p)->nayme[n - 1] = ch;
      if (eoln(treefile)) {
        fscanf(treefile, "%*[^\n]");
        getc(treefile);
      }
      ch = getc(treefile);
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
}  /* addelement */


void treeread()
{
  /* read a tree from the treefile and set up nodes and pointers */
  haslengths = true;
  ntips = 0;
  nextnode = 0;
  getch(&ch);
  addelement(&root, NULL);
  fscanf(treefile, "%*[^\n]");
  getc(treefile);
  uselengths = haslengths;
}  /* treeread */



void initialparms()
{
  /* initialize parameters */
  getplotter();
  plotrparms();
  if (dotmatrix)
    numlines = (long)floor(yunitspercm * ysize + 0.5) / strpdeep;
  xmargin = 0.08 * xsize;
  ymargin = 0.08 * ysize;
  xscale = xunitspercm;
  yscale = yunitspercm;
  style =  cladogram;
  grows = vertical;
  labelrotation = 45.0;
  nodespace = 3.0;
  stemlength = 0.05;
  treedepth = 0.5 / 0.95;
  rescaled = true;
  bscale = 1.0;
  if (uselengths)
    nodeposition = intermediate;
  else
    nodeposition = vshaped;
}  /* initialparms */


long showparms()
{
  long i;
  long numtochange;
  char input[32];
  Char ch, ch2,trash;

  putchar('\n');
  if (previewer == tek)
    printf("%c\f", escape);
  else {
    for (i = 1; i <= 24; i++)
      putchar('\n');
  }
  printf("Here are the settings: \n\n");
  printf(" (1)               Tree grows:  ");
  printf((grows == vertical) ? "Vertically\n" : "Horizontally\n");
  printf(" (2)            Style of tree:  %s\n",
	 (style == cladogram) ? "Cladogram" :
         (style == phenogram)  ? "Phenogram" :
         (style == curvogram) ? "Curvogram" :
         (style == eurogram)  ? "Eurogram"  : "Swoopogram");


  printf(" (3)       Use branch lengths:  ");
  if (haslengths) {
    if (uselengths)
      printf("Yes\n");
    else
      printf("No\n");
  } else
    printf("(no branch lengths available)\n");
  printf(" (4)          Angle of labels:");
  if (labelrotation < 10.0)
    printf("%5.1f\n", labelrotation);
  else
    printf("%6.1f\n", labelrotation);
  if (plotter == ray) {
    printf(" (5)       Horizontal margins:%6.2f pixels\n", xmargin);
    printf(" (5)         Vertical margins:%6.2f pixels\n", ymargin);
  } else {
    printf(" (5)       Horizontal margins:%6.2f cm\n", xmargin);
    printf(" (5)         Vertical margins:%6.2f cm\n", ymargin);
  }
  printf(" (6)   Scale of branch length:");
  if (rescaled)
    printf("  Automatically rescaled\n");
  else
    printf("  Fixed:%6.2f cm per unit branch length\n", bscale);
  printf(" (7)    Depth/Breadth of tree:%6.2f\n", treedepth);
  printf(" (8)   Stem-length/tree-depth:%6.2f\n", stemlength);
  printf(" (9) Character ht / tip space:%8.4f\n", 1.0 / nodespace);
  printf("(10)          Ancestral nodes:  %s\n",
	 (nodeposition == weighted)     ? "Weighted"     :
	 (nodeposition == intermediate) ? "Intermediate" :
	 (nodeposition == centered)     ? "Centered"     :
	 (nodeposition == inner)        ? "Inner"        :
	 "So tree is V-shaped");
  if (plotter == lw)
    printf("(11)          Font           :  %s\n",fontname);

  printf("\n Do you want to accept these? (Yes or No)\n");
  for (;;) {
    printf(" Type Y or N or the number (1-%2ld) of the one to change:\n",
           ((plotter == lw) ? 11L : 10L));
    gets(input);
    uppercase(&input[0]);
    numtochange = atoi(input);
    ch = input[0];
    if ((ch == 'Y' || ch == 'N') || (numtochange >= 1 && numtochange <= 11))
      break;
  }
 return (ch == 'Y') ? -1 : numtochange;
}  /* showparms */


void getparms(numtochange)
long numtochange;
{
  /* get from user the relevant parameters for the plotter and diagram */
  Char ch,trash;
  boolean ok;

  if (numtochange == 0) {
    do {
      printf(" Type the number of one that you want to change (1-%2ld):\n",
             ((plotter == lw) ? 11L : 10L));
      scanf("%hd%*[^\n]", &numtochange);
      trash=getchar();
    } while (numtochange < 1 || numtochange > ((plotter == lw) ? 11 : 10));
  }
  switch (numtochange) {

  case 1:
    if (grows == vertical)
      grows = horizontal;
    else
      grows = vertical;
    break;

  case 2:
    printf("\nWhat style tree is this to be:\n");
    printf("   Cladogram, Phenogram, curVogram, Eurogram,");
    printf("  or Swoopogram\n");
    printf(" (C, P, V, E, or S)\n");
    do {
      printf(" Choose one: \n");
      scanf("%c%*[^\n]", &ch);
      trash=getchar();
      uppercase(&ch);
    } while (ch != 'C' && ch != 'P' && ch != 'V' && ch != 'E' && ch != 'S');
    switch (ch) {

    case 'C':
      style = cladogram;
      break;

    case 'P':
      style = phenogram;
      break;

    case 'E':
      style = eurogram;
      break;

    case 'S':
      style = swoopogram;
      break;

    case 'V':
      style = curvogram;
      break;
    }
    break;

  case 3:
    if (haslengths) {
      uselengths = !uselengths;
      if (!uselengths)
        nodeposition = vshaped;
      else
        nodeposition = intermediate;
    } else {
      printf("Cannot use lengths since not all of them exist\n");
      uselengths = false;
    }
    break;

  case 4:
    printf("\n(Considering the tree as if it \"grew\" vertically:)\n");
    printf("Are the labels to be plotted vertically (90),\n");
    printf(" horizontally (0), or at a 45-degree angle?\n");
    do {
      printf(" Choose an angle in degrees from 90 to 0:\n");
      scanf("%lf%*[^\n]", &labelrotation);
      trash=getchar();
      uppercase(&ch);
    } while (labelrotation < 0.0 && labelrotation > 90.0);
    break;

  case 5:
    printf("\nThe tree will be drawn to fit in a rectangle which has \n");
    printf(" margins in the horizontal and vertical directions of:\n");
    if (plotter == ray)
      printf("%6.2f pixels (horizontal margin) and%6.2f pixels (vertical margin)\n",
             xmargin, ymargin);
    else
      printf("%6.2f cm (horizontal margin) and%6.2f cm (vertical margin)\n",
             xmargin, ymargin);
    putchar('\n');
    do {
      if (plotter == ray)
        printf(" New value (in pixels) of horizontal margin?\n");
      else
        printf(" New value (in cm) of horizontal margin?\n");
      scanf("%lf%*[^\n]", &xmargin);
      trash=getchar();
      ok = ((unsigned)xmargin < xsize / 2.0);
      if (!ok)
        printf(" Impossible value.  Please retype it.\n");
    } while (!ok);
    do {
      if (plotter == ray)
        printf(" New value (in pixels) of vertical margin?\n");
      else
        printf(" New value (in cm) of vertical margin?\n");
      scanf("%lf%*[^\n]", &ymargin);
      trash=getchar();
      ok = ((unsigned)ymargin < ysize / 2.0);
      if (!ok)
        printf(" Impossible value.  Please retype it.\n");
    } while (!ok);
    break;

  case 6:
    rescaled = !rescaled;
    if (!rescaled) {
      printf("Centimeters per unit branch length?\n");
      scanf("%lf%*[^\n]", &bscale);
      trash=getchar();
    }
    break;

  case 7:
    printf("New value of depth of tree as fraction of its breadth?\n");
    scanf("%lf%*[^\n]", &treedepth);
    trash=getchar();
    break;

  case 8:
    do {
      printf("New value of stem length as fraction of tree depth?\n");
      scanf("%lf%*[^\n]", &stemlength);
      trash=getchar();
    } while ((unsigned)stemlength >= 0.9);
    break;

  case 9:
    printf("New value of character height as fraction of tip spacing?\n");
    scanf("%lf%*[^\n]", &nodespace);
    trash=getchar();
    nodespace = 1.0 / nodespace;
    break;

  case 10:
    printf("Should interior node positions:\n");
    printf(" be Intermediate between their immediate descendants,\n");
    printf("    Weighted average of tip positions\n");
    printf("    Centered among their ultimate descendants\n");
    printf("    iNnermost of immediate descendants\n");
    printf(" or so that tree is V-shaped\n");
    do {
      printf(" (type I, W, C, N or V):\n");
      scanf("%c%*[^\n]", &ch);
      trash=getchar();
      uppercase(&ch);
    } while (ch != 'I' && ch != 'W' && ch != 'C' && ch != 'N' && ch != 'V');
    switch (ch) {

    case 'W':
      nodeposition = weighted;
      break;

    case 'I':
      nodeposition = intermediate;
      break;

    case 'C':
      nodeposition = centered;
      break;

    case 'N':
      nodeposition = inner;
      break;

    case 'V':
      nodeposition = vshaped;
      break;
    }
    break;
  case 11:
    printf("Enter font name or \"Hershey\" for the default font\n");
    gets(fontname);
    break;
  }
}  /* getparms */



void calctraverse(p, lengthsum,tipx)
node *p;
double lengthsum;
double *tipx;
{
  /* traverse to establish initial node coordinates */
  double x1, y1, x2, y2, x3, w1, w2, sumwx, sumw, nodeheight, rr;
  node *pp, *plast;

  if (p == root)
    nodeheight = 0.0;
  else if (uselengths)
    nodeheight = lengthsum + p->oldlen;
  else
    nodeheight = 1.0;
  if (nodeheight > maxheight)
    maxheight = nodeheight;
  if (p->tip) {
    p->xcoord = *tipx;
    if (uselengths)
      p->ycoord = nodeheight;
    else
      p->ycoord = 1.0;
    *tipx += tipspacing;
    return;
  }
  sumwx = 0.0;
  sumw = 0.0;
  pp = p->next;
  x3 = 0.0;
  do {
    calctraverse(pp->back, nodeheight,tipx);
    sumw += pp->back->tipsabove;
    sumwx += pp->back->tipsabove * pp->back->xcoord;
    if (fabs(pp->back->xcoord - 0.5) < fabs(x3 - 0.5))
      x3 = pp->back->xcoord;
    plast = pp;
    pp = pp->next;
  } while (pp != p);
  x1 = p->next->back->xcoord;
  x2 = plast->back->xcoord;
  y1 = p->next->back->ycoord;
  y2 = plast->back->ycoord;
  rr = 2 * (1.0 - stemlength) * treedepth * maxheight;
  switch (nodeposition) {

  case weighted:
    w1 = y1 - nodeheight;
    w2 = y2 - nodeheight;
    if (w1 + w2 <= 0.0)
      p->xcoord = (x1 + x2) / 2.0;
    else
      p->xcoord = (w2 * x1 + w1 * x2) / (w1 + w2);
    break;

  case intermediate:
    p->xcoord = (x1 + x2) / 2.0;
    break;

  case centered:
    p->xcoord = sumwx / sumw;
    break;

  case inner:
    p->xcoord = x3;
    break;

  case vshaped:
    p->xcoord = (x1 + x2 + (y1 - y2) / rr) / 2.0;
    break;
  }
  if (uselengths) {
    p->ycoord = nodeheight;
    return;
  }
  if (nodeposition != inner) {
    p->ycoord = (y1 + y2 - sqrt((y1 + y2) * (y1 + y2) - 4 * (y1 * y2 -
                 rr * rr * (x2 - p->xcoord) * (p->xcoord - x1)))) / 2.0;

    return;
  }
  if (fabs(x1 - 0.5) > fabs(x2 - 0.5)) {
    p->ycoord = y1 + x1 - x2;
    w1 = y2 - p->ycoord;
  } else {
    p->ycoord = y2 + x1 - x2;
    w1 = y1 - p->ycoord;
  }
  if (w1 < epsilon)
    p->ycoord -= fabs(x1 - x2);
}  /* calctraverse */


void calculate()
{
  /* compute coordinates for tree */
  double tipx;
  double sum, maxtextlength, textlength, firstlet, fontheight, angle;
  long i;
  for (i = 0; i < nextnode; i++)
    nodep[i]->xcoord = 0.0;
  for (i = 0; i < nextnode; i++)
    nodep[i]->ycoord = 0.0;
  maxheight = 0.0;
  maxtextlength = 0.0;
  if (nodep[0]->naymlength > 0)
    firstlet = lengthtext(nodep[0]->nayme, 1L,font);
  else
    firstlet = 0.0;
  sum = 0.0;
  tipx = 0.0;
    for (i = 0; i < nextnode; i++) {
    if (nodep[i]->tip) {
      textlength = lengthtext(nodep[i]->nayme, nodep[i]->naymlength,font);
      if (textlength > maxtextlength)
        maxtextlength = textlength;
    }
  }
  fontheight = font[2];
  angle = pi * labelrotation / 180.0;
  maxtextlength /= fontheight;
  textlength /= fontheight;
  firstlet /= fontheight;
  if (ntips > 1)
    labelheight = 1.0 / (nodespace * (ntips - 1));
  else
    labelheight = 1.0 / nodespace;
  if (angle < pi / 6.0)
    tipspacing = (nodespace + cos(angle) * (maxtextlength - 0.5)) * labelheight;
  else if (ntips > 1)
    tipspacing = 1.0 / (ntips - 1.0);
  else
    tipspacing = 1.0;
  topoflabels = labelheight *
                (1.0 + sin(angle) * (maxtextlength - 0.5) + cos(angle) * 0.5);
  rightoflabels = labelheight *
                  (cos(angle) * (textlength - 0.5) + sin(angle) * 0.5);
  leftoflabels = labelheight * (cos(angle) * firstlet * 0.5 + sin(angle) * 0.5);
  calctraverse(root, sum, &tipx);
  rooty = root->ycoord;
  for (i = 0; i < nextnode; i++) {
    if (rescaled) {
      nodep[i]->xcoord *= 1.0 - stemlength;
      nodep[i]->ycoord = stemlength * treedepth + (1.0 - stemlength) *
            treedepth * (nodep[i]->ycoord - rooty) / (maxheight - rooty);
    } else {
      nodep[i]->xcoord = nodep[i]->xcoord * (maxheight - rooty) / treedepth;
      nodep[i]->ycoord = stemlength / (1 - stemlength) * (maxheight - rooty) +
                         nodep[i]->ycoord;
    }
  }
  rooty = 0.0;
}  /* calculate */


void rescale()
{
  /* compute coordinates of tree for plot or preview device */
  long i;
  double treeheight, treewidth, extrax, extray, temp;

  treeheight = 0.0;
  for (i = 0; i < nextnode; i++) {
    if (nodep[i]->ycoord > treeheight)
      treeheight = nodep[i]->ycoord;
  }
  treewidth = (ntips - 1) * tipspacing + rightoflabels + leftoflabels;
  if (rescaled) {
    leftoflabels *= 1.0 - stemlength;
    rightoflabels *= 1.0 - stemlength;
    treewidth *= 1.0 - stemlength;
  } else {
    if (uselengths) {
      labelheight = labelheight * (maxheight - rooty) / treedepth;
      topoflabels = topoflabels * (maxheight - rooty) / treedepth;
      leftoflabels = leftoflabels * (maxheight - rooty) / treedepth;
      rightoflabels = rightoflabels * (maxheight - rooty) / treedepth;
      treewidth = treewidth * (maxheight - rooty) / treedepth;
    }
  }
  treeheight += topoflabels;
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
  for (i = 0; i < nextnode; i++) {
    nodep[i]->xcoord = expand * (nodep[i]->xcoord + leftoflabels);
    nodep[i]->ycoord = expand * (nodep[i]->ycoord - rooty);
    if (grows == horizontal) {
      temp = nodep[i]->ycoord;
      nodep[i]->ycoord = expand * treewidth - nodep[i]->xcoord;
      nodep[i]->xcoord = temp;
    }
    nodep[i]->xcoord += xmargin + extrax;
    nodep[i]->ycoord += ymargin + extray;
  }
  if (grows == vertical)
    rooty = ymargin + extray;
  else
    rooty = xmargin + extrax;
}  /* rescale */



void plottree(p, q)
node *p, *q;
{
  /* plot part or all of tree on the plotting device */
  long i;
  double x1, y1, x2, y2, x3, y3, f, g, h, fract, minny, miny;
  node *pp;

  x2 = xscale * (xoffset + p->xcoord);
  y2 = yscale * (yoffset + p->ycoord);
  if (p != root) {
    x1 = xscale * (xoffset + q->xcoord);
    y1 = yscale * (yoffset + q->ycoord);
    plot(penup, x1, y1);
    switch (style) {

    case cladogram:
      plot(pendown, x2, y2);
      break;

    case phenogram:
      if (grows == vertical)
        plot(pendown, x2, y1);
      else
        plot(pendown, x1, y2);
      plot(pendown, x2, y2);
      break;

    case curvogram:
      for (i = 1; i <= segments; i++) {
        f = (double)i / segments;
        g = (double)i / segments;
        h = 1.0 - sqrt(1.0 - g * g);
        if (grows == vertical) {
          x3 = x1 * (1.0 - f) + x2 * f;
          y3 = y1 + (y2 - y1) * h;
        } else {
          x3 = x1 + (x2 - x1) * h;
          y3 = y1 * (1.0 - f) + y2 * f;
        }
        plot(pendown, x3, y3);
      }
      break;

    case eurogram:
      if (grows == vertical)
        plot(pendown, x2, (2 * y1 + y2) / 3);
      else
        plot(pendown, (2 * x1 + x2) / 3, y2);
      plot(pendown, x2, y2);
      break;

    case swoopogram:
      if ((grows == vertical && fabs(y1 - y2) >= epsilon) ||
          (grows == horizontal && fabs(x1 - x2) >= epsilon)) {
        if (grows == vertical)
          miny = p->ycoord;
        else
          miny = p->xcoord;
        pp = q->next;
        while (pp != q) {
          if (grows == vertical)
            minny = pp->back->ycoord;
          else
            minny = pp->back->xcoord;
          if (minny < miny)
            miny = minny;
          pp = pp->next;
        }
        if (grows == vertical)
          miny = yscale * (yoffset + miny);
        else
          miny = xscale * (xoffset + miny);
        if (grows == vertical)
          fract = 0.3333 * (miny - y1) / (y2 - y1);
        else
          fract = 0.3333 * (miny - x1) / (x2 - x1);
        for (i = 1; i <= segments; i++) {
          f = (double)i / segments;
          if (f < fract)
            g = f / fract;
          else
            g = (f - fract) / (1.0 - fract);
          if (f < fract)
            h = fract * sqrt(1.0 - (1.0 - g) * (1.0 - g));
          else
            h = fract + (1.0 - fract) * (1.000001 - sqrt(1.000001 - g * g));
          if (grows == vertical) {
            x3 = x1 * (1.0 - f) + x2 * f;
            y3 = y1 + (y2 - y1) * h;
          } else {
            x3 = x1 + (x2 - x1) * h;
            y3 = y1 * (1.0 - f) + y2 * f;
          }
          plot(pendown, x3, y3);
        }
      }
      break;
    }
  } else {
    if (grows == vertical) {
      x1 = xscale * (xoffset + p->xcoord);
      y1 = yscale * (yoffset + rooty);
    } else {
      x1 = xscale * (xoffset + rooty);
      y1 = yscale * (yoffset + p->ycoord);
    }
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
  double compr, dx, dy, angle;
  node *lp;

  compr = xunitspercm / yunitspercm;
  if (penchange == yes)
    changepen(labelpen);
  angle = labelrotation * pi / 180.0;
  for (i = 0; i < (nextnode); i++) {
    if (nodep[i]->tip) {
      lp = nodep[i];
      dx = labelheight * expand * -0.70710 * cos(angle + pi / 4.0);
      dy = labelheight * expand * (1.0 - 0.70710 * sin(angle + pi / 4.0));
      if (grows == vertical)
        plottext(lp->nayme, lp->naymlength,
                 labelheight * expand * xscale / compr, compr,
                 xscale * (lp->xcoord + dx + xoffset),
                 yscale * (lp->ycoord + dy + yoffset),
		 -labelrotation, font,fontname);
      else
        plottext(lp->nayme, lp->naymlength, labelheight * expand * yscale,
                 compr, xscale * (lp->xcoord + dy + xoffset),
                 yscale * (lp->ycoord - dx + yoffset), 90.0 - labelrotation,
                 font,fontname);
    }
  }
  if (penchange == yes)
    changepen(treepen);
}  /* plotlabels */



main(argc, argv)
long argc;
Char *argv[];
{
  long i,n,stripedepth;
  boolean canbeplotted;
  char filename1[100];

#ifdef MAC
  OSErr retcode;
  FInfo  fndrinfo;
  macsetup("Drawgram","Preview");
  argv[0] = "Drawgram";
#endif
#ifdef TURBOC
  if ((registerbgidriver(EGAVGA_driver) <0) ||
      (registerbgidriver(Herc_driver) <0)   ||
      (registerbgidriver(CGA_driver) <0)){
    printf("Graphics error: %s ",grapherrormsg(graphresult()));
    exit(-1);}
#endif
  strcpy(fontname,"Hershey");

  openfile(&plotfile,PLOTFILE,"w",argv[0],pltfilename);
  openfile(&treefile,TREEFILE,"r",argv[0],NULL);

  printf("DRAWGRAM from PHYLIP version %s\n",VERSION);
  printf("Reading tree ... \n");
  treeread();

  printf("Tree has been read.\nLoading the font .... \n");
  loadfont(font,argv[0]);
  printf("Font loaded.\n");
  previewing = false;
  initialparms();
  canbeplotted = false;
  while (!canbeplotted) {
    do {
      n=showparms();
      if (n  != -1)
        getparms(n);
    } while (n  != -1);
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
  numlines = dotmatrix ? ((long)floor(yunitspercm * ysize + 0.5)/strpdeep) : 1;
  if (plotter != ibmpc)
    printf("Writing plot file ...\n");
  drawit(fontname,&xoffset,&yoffset,numlines,root);
  finishplotter();
  FClose(plotfile);
  FClose(treefile);
  printf("Finished.\nEnd of run.\n");
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
mem = (MALLOCRETURN *)malloc((size_t)x);
if (!mem)
     memerror();
else
     return (MALLOCRETURN *)mem;

}

