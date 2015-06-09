#include "phylip.h"

/* version 3.52c. (c) Copyright 1993 by Joseph Felsenstein.
   Written by Joseph Felsenstein.  Permission is granted to copy and use
   this program provided no fee is charged for it and provided that this
   copyright notice is not removed. */

#define maxsp           500   /* maximum number of species               */
#define maxsz           999   /* size of pointer array.  >= 2*maxsp - 1  */
                              /* this can be large without eating memory */
#define nmlngth         30    /* number of characters in species name    */

#define ibmpc0          false
#define ansi0           true
#define vt520           false
#define downn           2
#define overr           4

  typedef enum {    valid, remoov, quit    } reslttype;


typedef enum {
  horiz, vert, up, over, upcorner, downcorner, aa, cc, gg, tt, quest
  } chartype;
/* nodes will form a binary tree */

typedef struct node {       /* describes a tip species or an ancestor        */
  struct node *next, *back; /* pointers to nodes                             */
  long index;               /* number of the node                            */
  boolean tip;              /* present species are tips of tree              */
  boolean hasname;          /* true if tip has a name                        */
  char    name[nmlngth];    /* if it does, its name                          */
  boolean haslength;        /* true if node has a valid length               */
  double length;            /* length of this branch                         */
  double lbeyond;           /* distance beyond this node to most distant tip */
                            /* on left side                                  */
  double rbeyond;           /* distance beyond this node to most distant tip */
                            /* on right side                                 */
  long xcoord, ycoord;      /* used by printree                              */
  long ymin, ymax;
} node;

typedef node *pointarray[maxsz];
typedef enum {  arb, use, spec } howtree;
typedef enum {  left, down, upp, right } adjwindow;

typedef enum {
  rearr, transp, flipp, reroott, mpoint, lengthh, none
  } rearrtype;


typedef struct lenstore {
  long index;
  double length;
  boolean haslength;
} lenstore;

typedef struct midstoret {
  node *leftp, *rightp;
  double leftlength, rightlength;
} midstoret;


node *root, *garbage;
FILE *intree, *outtree;
long spp, spp2, nonodes, outgrno, screenwidth, vscreenwidth,
     screenlines, col, treenumber, leftedge, topedge, treelines,
     hscroll, vscroll, scrollinc;
/* spp    = number of species
  nonodes = number of nodes in tree
  outgrno indicates outgroup */

double     backlength, backlength2, rlen1, llen1, rlen2, llen2, trweight;
boolean    waswritten, ibmpc, vt52, ansi, hasmult, haslengths,
           nolengths, nexus;
pointarray treenode;                /* pointers to all nodes in tree */

boolean    reversed[11];
boolean    graphic[11];
Char       cch[11];
howtree    how;
lenstore   oldls, templs;
midstoret  midstore;
char       intreename[100],outtreename[100];

long       atwhat, what, fromwhere, towhere, oldoutgrno;
boolean    subtree, written, oldwritten, restoring, wasleft, oldleft, readnext;
double     rtlength;
node      *nuroot;
rearrtype  lastop;
node      *washere;
Char      ch;


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

long readlong(prompt)
char *prompt;
{
int err;
long res;
char string[100];
do {
  printf("%s",prompt);
  gets(string);
  if (sscanf(string,"%ld",&res) == 1)
    break;
 } while (1);
return res;
}

void inpnum(n, success)
     long *n;
     boolean *success;
{
  int fields;
  char line[100];
  gets(line);
  *n = atof(line);
  fields = sscanf(line,"%ld",n);
  *success = (fields == 1);

}  /* inpnum */

void gnu(p)
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


void chuck(p)
     node *p;
{
  /* collect garbage on p -- put it on front of garbage list */
  p->next = garbage;
  garbage = p;
}  /* chuck */


void gdispose(p)
     node *p;
{
  /* go through tree throwing away nodes */
  node *q, *r;

  if (p->tip)
    return;
  q = p->next;
  while (q != p) {
    gdispose(q->back);
    q->tip = false;
    q->hasname = false;
    q->haslength = false;
    r = q;
    q = q->next;
    chuck(r);
  }
  q->tip = false;
  q->hasname = false;
  q->haslength = false;
  chuck(q);
}  /*  gdispose */


void maketriad(p, index)
     node **p;
     long index;
{
  long i;
  node *q;
  q = NULL;
  for (i = 1; i <= 3; i++) {
    gnu(p);
    (*p)->index = index;
    (*p)->hasname = false;
    (*p)->haslength = false;
    (*p)->next = q;
    q = *p;
  }
  (*p)->next->next->next = *p;
  q = (*p)->next;
  while (*p != q) {
    (*p)->back = NULL;
    (*p)->tip = false;
    *p = (*p)->next;
  }
  treenode[index - 1] = *p;
}  /* maketriad */


void maketip(p, index)
     node **p;
     long index;
{
  gnu(p);
  (*p)->index = index;
  (*p)->tip = true;
  (*p)->hasname = false;
  (*p)->haslength = false;
  treenode[index - 1] = *p;
}  /* maketip */


void getoptions()
{
  /* interactively set options */
  Char ch;
  boolean done, gotopt;

  how = use;
  outgrno = 1;
  do {
    if (ansi || ibmpc)
      printf("\033[2J\033[H");
    else if (vt52)
      printf("\033E\033H");
    else
      putchar('\n');
    printf("\nTree Rearrangement, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  U      Initial tree (arbitrary, user, specify)?");
    if (how == arb)
      printf("  Arbitrary\n");
    else if (how == use)
      printf("  User tree from tree file\n");
    else
      printf("  Tree you specify\n");
    printf("  N      Use the Nexus format to write out trees?");
    if (nexus)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  0           Graphics type (IBM PC, VT52, ANSI)?");
    if (ibmpc)
      printf("  IBM PC\n");
    if (ansi )
      printf("  ANSI\n");
    if (vt52)
      printf("  VT52\n");
    if (!(ibmpc || vt52 || ansi))
      printf("  (none)\n");
    printf("  W   Width of terminal screen, of plotting area?");
    printf("%4ld, %2ld\n", screenwidth, vscreenwidth);
    printf("  L                    Number of lines on screen?");
    printf("%4ld\n", screenlines);
    do {
      printf("\nAre these settings correct?");
      printf(" (type Y or the letter for one to change)\n");
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
	ch = ' ';
      ch = (isupper(ch)) ? ch : toupper(ch);
      done = (ch == 'Y');
      gotopt = (ch == 'U' || ch == 'N' || ch == '0' || ch == 'L' || ch == 'W');
      if (gotopt) {
	switch (ch) {
	
	case 'U':
	  if (how == arb)
	    how = use;
	  else if (how == use)
	    how = spec;
	  else
	    how = arb;
	  break;

	case 'N':
	  nexus = !nexus;
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

	case 'L':
	  do {
	    screenlines = readlong("Number of lines on screen?\n");
	  } while (screenlines <= 12);
	  break;

	case 'W':
	  screenwidth= readlong("Width of terminal screen (in characters)?\n");
	  vscreenwidth=readlong("Width of plotting area (in characters)?\n");
	  break;
	}
      }
      if (!(gotopt || done))
	printf("Not a possible option!\n");
    } while (!(gotopt || done));
  } while (!done);
  if (scrollinc < screenwidth / 2.0)
    hscroll = scrollinc;
  else
    hscroll = screenwidth / 2;
  if (scrollinc < screenlines / 2.0)
    vscroll = scrollinc;
  else
    vscroll = screenlines / 2;
}  /* getoptions */




void configure()
{
  /* configure to machine -- set up special characters */
  chartype a;

  for (a = horiz; (long)a <= (long)quest; a = (chartype)((long)a + 1))
    reversed[(long)a] = false;
  for (a = horiz; (long)a <= (long)quest; a = (chartype)((long)a + 1))
    graphic[(long)a] = false;
  if (ibmpc) {
    cch[(long)horiz] = '>';
    cch[(long)vert] = 186;
    graphic[(long)vert] = true;
    cch[(long)up] = 186;
    graphic[(long)up] = true;
    cch[(long)over] = 205;
    graphic[(long)over] = true;
    cch[(long)upcorner] = 200;
    graphic[(long)upcorner] = true;
    cch[(long)downcorner] = 201;
    graphic[(long)downcorner] = true;
    return;
  }
  if (vt52) {
    cch[(long)horiz] = '>';
    cch[(long)vert] = cch[(long)horiz];
    reversed[(long)vert] = true;
    cch[(long)up] = '`';
    graphic[(long)up] = true;
    cch[(long)over] = 'a';
    graphic[(long)over] = true;
    cch[(long)upcorner] = 'e';
    graphic[(long)upcorner] = true;
    cch[(long)downcorner] = 'f';
    graphic[(long)downcorner] = true;
    return;
  }
  if (ansi) {
    cch[(long)horiz] = '>';
    cch[(long)vert] = cch[(long)horiz];
    reversed[(long)vert] = true;
    cch[(long)up] = 'x';
    graphic[(long)up] = true;
    cch[(long)over] = 'q';
    graphic[(long)over] = true;
    cch[(long)upcorner] = 'm';
    graphic[(long)upcorner] = true;
    cch[(long)downcorner] = 'l';
    graphic[(long)downcorner] = true;
    return;
  }
  cch[(long)horiz] = '>';
  cch[(long)vert] = ' ';
  cch[(long)up] = '!';
  cch[(long)upcorner] = '`';
  cch[(long)downcorner] = ',';
  cch[(long)over] = '-';
}  /* configure */


void pregraph()
{
  /* turn on graphic characters */
  if (vt52)
    printf("\033F");
  if (ansi) {
    printf("\033(0");
    printf("\033[10m");
  }
}  /* pregraph */

void prereverse()
{
  /* turn on reverse video */
  if (vt52)
    printf("\033p");
  if (ansi)
    printf("\033[7m");
}  /* prereverse */


void prefix(a)
chartype a;
{
  /* give prefix appropriate for this character */
  if (reversed[(long)a])
    prereverse();
  if (graphic[(long)a])
    pregraph();
}  /* prefix */


void postgraph()
{
  /* turn off graphic characters */
  if (vt52)
    printf("\033G");
  if (ansi) {
    printf("\033[11m");
    printf("\033(B");
  }
}  /* postgraph */

void postreverse()
{
  /* turn off reverse video */
  if (vt52)
    printf("\033q");
  if (ansi)
    printf("\033[0m");
}  /* postreverse */


void postfix(a)
chartype a;
{
  /* give postfix appropriate for this character */
  if (reversed[(long)a])
    postreverse();
  if (graphic[(long)a])
    postgraph();
}  /* postfix */



void ltrav(p, localhl)
node *p;
boolean *localhl;
{
  node *q;

  if (p->tip) {
    (*localhl) = ((*localhl) && p->haslength);
    return;
  }
  q = p->next;
  do {
    (*localhl) = ((*localhl) && q->haslength);
    if ((*localhl))
      ltrav(q->back, localhl);
    q = q->next;
  } while (p != q);
}  /* ltrav */

boolean ifhaslengths()
{
  /* return true if every branch in tree has a length */
  /* Local variables for ifhaslengths: */
  boolean localhl;
  localhl = true;
  ltrav(root, &localhl);
  return localhl;
}  /* ifhaslengths */

void add(below, newtip, newfork, restoring, wasleft)
node *below, *newtip, *newfork;
boolean *restoring,*wasleft;
{
  /* inserts the nodes newfork and its left descendant, newtip,
    to the tree.  below becomes newfork's right descendant */
  boolean putleft;
  node *leftdesc, *rtdesc;
  double length;

  if (below != treenode[below->index - 1])
    below = treenode[below->index - 1];
  if (below->back != NULL)
    below->back->back = newfork;
  newfork->back = below->back;
  putleft = true;
  if (*restoring)
    putleft = (*wasleft);
  if (putleft) {
    leftdesc = newtip;
    rtdesc = below;
  } else {
    leftdesc = below;
    rtdesc = newtip;
  }
  rtdesc->back = newfork->next->next;
  newfork->next->next->back = rtdesc;
  newfork->next->back = leftdesc;
  leftdesc->back = newfork->next;
  if (root == below)
    root = newfork;
  root->back = NULL;
  if (!haslengths)
    return;
  if (*restoring) {
    if (newfork->back != NULL) {
      length = below->length - backlength;
      below->length = length;
      below->back->length = length;
      newfork->length = backlength;
      newfork->back->length = backlength;
    } else {
      newfork->next->length = llen1;
      newfork->next->back->length = llen1;
      newfork->next->next->length = rlen1;
      newfork->next->next->back->length = rlen1;
      rlen1 = rlen2;
      llen1 = llen2;
    }
    backlength = backlength2;
  } else {
    if (newfork->back != NULL) {
      length = newfork->back->length / 2.0;
      newfork->length = length;
      newfork->back->length = length;
      below->length = length;
      below->back->length = length;
    } else {
      length = newtip->length / 2.0;
      newtip->length = length;
      newtip->back->length = length;
      below->length = length;
      below->back->length = length;
      below->haslength = true;
    }
  }
  newtip->back->length = newtip->length;
}  /* add */

void re_move(item, fork, restoring, wasleft)
node **item, **fork;
boolean *restoring,*wasleft;
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
  if (*item == (*fork)->next->back) {
    if (root == *fork)
      root = (*fork)->next->next->back;
    (*wasleft) = true;
  } else {
    if (root == *fork)
      root = (*fork)->next->back;
    (*wasleft) = false;
  }
  if (haslengths) {
    if (*restoring) {
      rlen2 = (*fork)->next->next->length;
      llen2 = (*fork)->next->length;
    } else {
      rlen1 = (*fork)->next->next->length;
      llen1 = (*fork)->next->length;
    }
  }
  p = (*item)->back->next->back;
  q = (*item)->back->next->next->back;
  if (p != NULL)
    p->back = q;
  if (q != NULL)
    q->back = p;
  if (haslengths) {
    if (*restoring) {
      if ((*fork)->back == NULL) {
	(*item)->length *= 2.0;
	(*item)->back->length = (*item)->length;
      }
      backlength2 = (*fork)->length;
    } else
      backlength = (*fork)->length;
    if (p != NULL && q != NULL) {
      p->length += q->length;
      q->length = p->length;
    } else
      (*item)->length = (*fork)->next->length + (*fork)->next->next->length;
  }
  (*fork)->back = NULL;
  p = (*fork)->next;
  while (p != *fork) {
    p->back = NULL;
    p = p->next;
  }
  (*item)->back = NULL;
}  /* re_move */

void reroot(outgroup)
node *outgroup;
{
  /* reorients tree, putting outgroup in desired position. */
  node *p, *q, *newbottom, *oldbottom, *tempwashere;
  double outgrlength, temprtlength;

  if (!(outgroup->back->index != root->index || restoring))
    return;
  oldoutgrno = root->next->back->index;
  if (haslengths && !restoring)
    outgrlength = outgroup->length;
  if (restoring)
    outgroup = washere;
  newbottom = outgroup->back;
  p = treenode[newbottom->index - 1]->back;
  while (p->index != root->index) {
    oldbottom = treenode[p->index - 1];
    treenode[p->index - 1] = p;
    p = oldbottom->back;
  }
  p = root->next;
  q = root->next->next;
  tempwashere = root->next->next->back;
  if (haslengths)
    temprtlength = root->next->next->length;
  p->back->back = q->back;
  q->back->back = p->back;
  if (haslengths) {
    p->back->length = p->length + q->length;
    q->back->length = p->back->length;
  }
  if (restoring) {
    p->back = washere->back;
    p->back->back = p;
    q->back = washere;
    q->back->back = q;
    if (haslengths) {
      p->back->length -= rtlength;
      p->length = p->back->length;
      q->back->length = rtlength;
      q->length = rtlength;
      rtlength = temprtlength;
    }
  } else {
    p->back = outgroup;
    q->back = outgroup->back;
    outgroup->back->back = root->next->next;
    outgroup->back = root->next;
    if (haslengths) {
      root->next->length = outgrlength / 2.0;
      root->next->back->length = outgrlength / 2.0;
      root->next->next->length = outgrlength / 2.0;
      root->next->next->back->length = outgrlength / 2.0;
    }
  }
  washere = tempwashere;
  if (haslengths)
    rtlength = temprtlength;
  treenode[newbottom->index - 1] = newbottom;
}  /* reroot */


void ltrav_(p, lengthsum, lmin, tipmax, across,maxchar)
node *p;
double  lengthsum, lmin;
double *tipmax;
long   *across;
long   *maxchar;
{
  node *q;
  long rchar, nl;
  double sublength;

  if (p->tip) {
    if (lengthsum > (*tipmax))
      (*tipmax) = lengthsum;
    if (lmin == 0.0)
      return;
    rchar = (long)(lengthsum / (*tipmax) * (*across) + 0.5);
    nl = strlen(treenode[p->index - 1]->name);
    if (rchar + nl > (*maxchar))
      (*across) = (*maxchar) - (long)(nl * (*tipmax) / lengthsum + 0.5);
    return;
  }
  q = p->next;
  do {
    if (q->length >= lmin)
      sublength = q->length;
    else
      sublength = lmin;
    ltrav_(q->back, lengthsum + sublength, lmin,tipmax,across,maxchar);
    q = q->next;
  } while (p != q);
}


void precoord(nuroot, subtree,tipmax,across)
node    *nuroot;
boolean *subtree;
double  *tipmax;
long    *across;
{
  /*                                                              *
   *  set tipmax and across so that tree is scaled to screenwidth */

  double oldtipmax, minimum;
  long i;
  long maxchar;

  (*tipmax) = 0.0;
  if ((*subtree))
    maxchar = vscreenwidth - 13;
  else
    maxchar = vscreenwidth - 5;
  (*across) = maxchar;
  ltrav_(nuroot, 0.0, 0.0, tipmax,across,&maxchar);
  i = 0;
  do {
    oldtipmax = (*tipmax);
    minimum = 3.0 / (*across) * (*tipmax);
    ltrav_(nuroot, 0.0, minimum,tipmax,across,&maxchar);
    i++;
  } while (fabs((*tipmax) - oldtipmax) > 0.01 * oldtipmax && i <= 40);
}  /* precoord */

void coordinates(p, lengthsum, across, tipy, tipmax)
node   *p;
double lengthsum;
long   *across,*tipy;
double *tipmax;
{
  /* establishes coordinates of nodes for display with lengths */
  node *q, *first, *last;

  if (p->tip) {
    p->xcoord = (long)((*across) * lengthsum / (*tipmax) + 0.5);
    p->ycoord = (*tipy);
    p->ymin   = (*tipy);
    p->ymax   = (*tipy);
    (*tipy)  += downn;
    return;
  }
  q = p->next;
  do {
   coordinates(q->back, lengthsum + q->length, across, tipy, tipmax);
    q = q->next;
  } while (p != q);
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (long)((*across) * lengthsum / (*tipmax) + 0.5);
  p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */

void flatcoordinates(p, tipy)
node *p;
long *tipy;
{
  /* establishes coordinates of nodes for display without lengths */
  node *q, *first, *last;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = (*tipy);
    p->ymin   = (*tipy);
    p->ymax   = (*tipy);
    (*tipy) += downn;
    return;
  }
  q = p->next;
  do {
    flatcoordinates(q->back, tipy);
    q = q->next;
  } while (p != q);
  first = p->next->back;
  q = p->next;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (last->ymax - first->ymin) * 3 / 2;
  p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* flatcoordinates */


void grwrite(c, num, pos)
chartype c;
long num;
long *pos;
{
  long i;

  prefix(c);
  for (i = 1; i <= num; i++) {
    if ((*pos) >= leftedge && (*pos) - leftedge + 1 < screenwidth)
      putchar(cch[(long)c]);
    (*pos)++;
  }
  postfix(c);
}

void chwrite(ch, num, pos)
Char ch;
long num;
long *pos;
{
  long i;

  for (i = 1; i <= num; i++) {
    if ((*pos) >= leftedge && (*pos) - leftedge + 1 < screenwidth)
      putchar(ch);
    (*pos)++;
  }
}

void nnwrite(nodenum, num, pos)
long nodenum, num;
long *pos;
{
  int i;
  long leftx;

  leftx = leftedge - (*pos);
  if ((*pos) >= leftedge && (*pos) - leftedge + num < screenwidth)
    printf("%*ld", (int)num, nodenum);
  else if (leftx > 0 && leftx < 3)
    for(i=0;i<(int)(num-leftx);i++)
      printf(" ");
  (*pos) += num;
}

void stwrite(s, length, pos)
Char *s;
long length;
long *pos;
{
  if ((*pos) >= leftedge && (*pos) - leftedge + 1 < screenwidth)
    printf("%*s", (int)length, s);
  (*pos) += length;
}

void drawline(i, nuroot, subtree)
long i;
node *nuroot;
boolean *subtree;
{
  /* draws one row of the tree diagram by moving up tree */
  long pos;
  node *p, *q, *r, *first, *last;
  long n, j;
  boolean extra, done;
  chartype c, d;
  pos = 1;
  p = nuroot;
  q = nuroot;

  extra = false;
  if (i == p->ycoord && (p == root || (*subtree))) {
    c = over;
    if ((*subtree))
      stwrite("Subtree:", 8L, &pos);
    if (p->index >= 100)
      nnwrite(p->index, 3L, &pos);
    else if (p->index >= 10) {
      grwrite(c, 1L, &pos);
      nnwrite(p->index, 2L, &pos);
    } else {
      grwrite(c, 2L, &pos);
      nnwrite(p->index, 1L, &pos);
    }
    extra = true;
  } else {
    if ((*subtree))
      stwrite("          ", 10L, &pos);
    else
      stwrite("  ", 2L, &pos);
  }
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
    if (haslengths && !nolengths)
      n = q->xcoord - p->xcoord;
    else
      n = p->xcoord - q->xcoord;
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if (q->ycoord == i && !done) {
      if (q->ycoord > p->ycoord)
	d = upcorner;
      else
	d = downcorner;
      c = over;
      if (!haslengths && !q->haslength)
	c = horiz;
      if (n > 1 || q->tip) {
	grwrite(d, 1L, &pos);
	grwrite(c, n - 3, &pos);
      }
      if (q->index >= 100)
	nnwrite(q->index, 3L, &pos);
      else if (q->index >= 10) {
	grwrite(c, 1L, &pos);
	nnwrite(q->index, 2L, &pos);
      } else {
	grwrite(c, 2L, &pos);
	nnwrite(q->index, 1L, &pos);
      }
      extra = true;
    } else if (!q->tip) {
      if (last->ycoord > i && first->ycoord < i && i != p->ycoord) {
	c = up;
	grwrite(c, 1L, &pos);
	chwrite(' ', n - 1, &pos);
      } else
	chwrite(' ', n, &pos);
    } else
      chwrite(' ', n, &pos);
    if (p != q)
      p = q;
  } while (!done);
  if (p->ycoord == i && p->tip) {
    if (p->hasname) {
      n = 0;
      for (j = 1; j <= nmlngth; j++) {
	if (treenode[p->index - 1]->name[j - 1] != '\0')
	  n = j;
      }
      chwrite(':', 1L, &pos);
      for (j = 0; j < n; j++)
	chwrite(treenode[p->index - 1]->name[j], 1L, &pos);
    }
  }
  putchar('\n');
}  /* drawline */

void printree()
{
  /* prints out diagram of the tree */
  /* Local variables for printree:  */
  long across, tipy;
  double tipmax;
  long i, rover, dow, vmargin;

  haslengths = ifhaslengths();
  if (!subtree)
    nuroot = root;
  printf((ansi || ibmpc) ? "\033[2J\033[H" :
         (vt52)          ? "\033E\033H"    : "\n");
  tipy = 1;
  rover = 100 / spp;
  if (rover > overr)
    rover = overr;
  dow = downn;
  if (spp * dow > screenlines && !subtree) {
    dow--;
    rover--;
  }
  if (haslengths && !nolengths) {
    precoord(nuroot, &subtree,&tipmax,&across);
    coordinates(nuroot, 0.0, &across,&tipy,&tipmax);
  } else
    flatcoordinates(nuroot, &tipy);

  vmargin = 2;
  treelines = tipy - dow;
  if (topedge != 1) {
    printf("** %ld lines above screen **\n", topedge - 1);
    vmargin++;
  }
  if ((treelines - topedge + 1) > (screenlines - vmargin))
    vmargin++;
  for (i = 1; i <= treelines; i++) {
    if (i >= topedge && i < topedge + screenlines - vmargin)
      drawline(i, nuroot,&subtree);
  }
  if ((treelines - topedge + 1) > (screenlines - vmargin)) {
    printf("** %ld", treelines - (topedge - 1 + screenlines - vmargin));
    printf(" lines below screen **\n");
  }
  if (treelines - topedge + vmargin + 1 < screenlines)
    putchar('\n');
}  /* printree */

#undef downn
#undef overr

void togglelengths()
{
  nolengths = !nolengths;
  printree();
}



void arbitree()
{
  long i;
  node *newtip, *newfork;
  spp = readlong("How many species?\n");
  nonodes = spp * 2 - 1;
  maketip(&root, 1L);
  maketip(&newtip, 2L);
  maketriad(&newfork, spp + 1);
  add(root, newtip, newfork, &restoring,&wasleft);
  for (i = 3; i <= spp; i++) {
    maketip(&newtip, i);
    maketriad(&newfork, spp + i - 1);
    add(treenode[spp + i - 3], newtip, newfork,&restoring, &wasleft);
  }
}  /* arbitree */

void yourtree()
{
  long i, j, k;
  boolean ok, done;
  node *newtip, *newfork;

  spp = 2;
  nonodes = spp * 2 - 1;
  maketip(&root, 1L);
  maketip(&newtip, 2L);
  maketriad(&newfork, spp + 2);
  add(root, newtip, newfork,&restoring,&wasleft);
  i = 2;
  do {
    i++;
    printree();
    printf("Enter 0 to stop building tree.\n");
    printf("Add species%3ld", i);
    do {
      printf("\n before node (type number): ");
      inpnum(&j, &ok);
      ok = (ok && ((unsigned long)j < i || (j > spp + 1 && j < spp + i)));
      if (!ok)
	printf("Impossible number. Please try again:\n");
    } while (!ok);
    done = (j == 0);
    if (!done) {
      k = spp * 2 + 1;
      do {
	treenode[k - 1] = treenode[k - 2];
	treenode[k - 1]->index = k;
	treenode[k - 1]->next->index = k;
	treenode[k - 1]->next->next->index = k;
	k--;
      } while (k != spp + 2);
      if (j > spp + 1)
	j++;
      spp++;
      nonodes = spp * 2 - 1;
      maketip(&newtip, i);
      maketriad(&newfork, spp + i);
      add(treenode[j - 1], newtip, newfork, &restoring, &wasleft);
    }
  } while (!done);
  for (i = spp + 1; i < (spp * 2); i++) {
    treenode[i - 1] = treenode[i];
    treenode[i - 1]->index = i;
    treenode[i - 1]->next->index = i;
    treenode[i - 1]->next->next->index = i;
  }
  nonodes = spp * 2 - 1;
}  /* yourtree */

void getch(c)
Char *c;
{
  /* get next nonblank character */
  do {
    if (eoln(intree)) {
      fscanf(intree, "%*[^\n]");
      getc(intree);
    }
    *c = getc(intree);
    if (*c == '\n')
      *c = ' ';
  } while (*c == ' ');
}  /* getch */

void findch(c)
Char c;
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
	if (ch == ')')
	  done = true;
      }
    } else if (c == ';') {
      if (ch != ';') {
	printf("\nERROR IN USER TREE: UNMATCHED PARENTHESIS OR MISSING SEMICOLON\n");
	exit(-1);
      } else
	done = true;
    }
    if (ch != ')')
      continue;
    if (eoln(intree)) {
      fscanf(intree, "%*[^\n]");
      getc(intree);
    }
    ch = getc(intree);
    if (ch == '\n')
      ch = ' ';
  }
}  /* findch */


void processlength(p)
node *p;
{
int i=0;
char tmp[100];

do {(getch(&tmp[i]));}
  while (strchr("0123456789-.+-",tmp[i++]) != NULL);
ch = tmp[i-1];
tmp[--i]=0;
#ifdef VMS
sscanf(tmp,"%lf",&(p->length));
#else
p->length    = atof(tmp);
#endif
p->haslength = true;
}


void addelement(p, endoffile,nextnode,lparens,n)
node **p;
boolean *endoffile;
long *nextnode,*lparens,*n;
{
  /* recursive procedure adds nodes to user-defined tree */
  node *q, *x, *temp;
  long i;
  char name[100];

  do {
    (*endoffile) = eof(intree);
    if (!(*endoffile)) {
      if (eoln(intree)) {
	fscanf(intree, "%*[^\n]");
	getc(intree);
      }
    }
    (*endoffile) = eof(intree);
    if (!(*endoffile)) {
      ch = getc(intree);
      if (ch == '\n')
	ch = ' ';
    }
  } while (!(ch != ' ' || (*endoffile)));
  if (*endoffile)
    exit(-1);
  if (*p == root && (*endoffile))
    printf("ERROR IN USER TREE: INPUT FILE TRUNCATED\n");
  if (ch == '(' ) {
    if ((*lparens) >= spp - 1) {
      printf("\nERROR IN USER TREE: TOO MANY LEFT PARENTHESES\n");
      exit(-1);
    } else {
      (*nextnode)++;
      (*lparens)++;
      maketriad(&q, (*nextnode));
      addelement(&q->next->back,endoffile,nextnode,lparens,n);
      q->next->back->back = q->next;
      q->next->haslength = q->next->back->haslength;
      if (q->next->haslength)
	q->next->length = q->next->back->length;
      findch(',');
      addelement(&q->next->next->back,endoffile,nextnode,lparens,n);
      q->next->next->back->back = q->next->next;
      q->next->next->haslength = q->next->next->back->haslength;
      if (q->next->next->haslength)
	q->next->next->length = q->next->next->back->length;
      if (ch == ',') {
	if (*p == root) {
	  (*nextnode)++;
	  maketriad(&x, (*nextnode));
	  temp = q->next->next->back;
	  x->back = q->next->next;
	  q->next->next->back = x;
	  temp->back = x->next;
	  x->next->back = temp;
	  x->next->length = x->next->back->length;
	  x->next->haslength = x->next->back->haslength;
	  x->haslength = true;
	  x->back->haslength = true;
	  x->length = 0.0;
	  x->back->length = 0.0;
	  addelement(&x->next->next->back,endoffile,nextnode,lparens,n);
	  x->next->next->back->back = x->next->next;
	  x->next->next->haslength = x->next->next->back->haslength;
	  if (x->next->next->haslength)
	    x->next->next->length = x->next->next->back->length;
	} else {
	  printf("ERROR IN USER TREE: TRIFURCATION NOT ALLOWED ");
	  printf("EXCEPT AT ROOT\n");
	  exit(-1);
	}
      }
      findch(')');
      *p = q;
    }
  } else {
    spp2++;
     for (i=0;i<nmlngth;++i)
       name[i] = ' ';
    (*n) = 1;
    while (ch != ',' && ch != ')' && ch != ':' && (*n) <= nmlngth) {
      if (ch == '_')
	ch = ' ';
      name[(*n) - 1] = ch;
      name[(*n)] = '\0';
      if (eoln(intree)) {
	fscanf(intree, "%*[^\n]");
	getc(intree);
      }
      ch = getc(intree);
      if (ch == '\n')
	ch = ' ';
      (*n)++;
    }
    maketip(p, spp2);
    if ((*n) > 1){
      (*p)->hasname = true;
      strncpy(treenode[spp2 - 1]->name,name,nmlngth);
    }
  }
  if (ch == ':')
    processlength(*p);
}  /* addelement */


void treeread()
{
  /* read in user-defined tree and set it up */
  long n, nextnode, lparens;
  boolean endoffile;
  long i;

  hasmult = false;
  if (!readnext)
    openfile(&intree,INTREE,"r","retree",intreename);
  spp2 = 0;
  spp = maxsp;
  nonodes = spp * 2 - 1;
  maketriad(&root, spp + 1);
  nextnode = spp;
  root->back = NULL;
  lparens = 0;
  addelement(&root, &endoffile,&nextnode,&lparens,&n);
  if (ch == '[' ) {
    if (!eoln(intree)) {
      hasmult = true;
      fscanf(intree, "%lf", &trweight);
     getch(&ch);
      if (ch != ']') {
	printf("ERROR: MISSING RIGHT SQUARE BRACKET\n");
	exit(-1);
      } else {
	getch(&ch);
	if (ch != ';') {
	  printf("ERROR: MISSING SEMICOLON AFTER SQUARE BRACKETS\n");
	  exit(-1);
	}
      }
    }
  }
  findch(';');
  for (i = 1; i < spp2; i++) {
    treenode[spp2 + i - 1] = treenode[spp + i - 1];
    treenode[spp2 + i - 1]->index = spp2 + i;
    treenode[spp2 + i - 1]->next->index = spp2 + i;
    treenode[spp2 + i - 1]->next->next->index = spp2 + i;
  }
  spp = spp2;
  nonodes = spp * 2 - 1;
  root = treenode[spp];
  printf("\n\n");
}  /* treeread */




void buildtree()
{
  long i;

  switch (how) {

  case arb:
    arbitree();
    break;

  case use:
    printf("\nReading tree file ...\n\n");
    treeread();
    break;

  case spec:
    yourtree();
    break;
  }
    outgrno = root->next->back->index;
    reroot(treenode[outgrno - 1]);
}  /* buildtree */

void unbuildtree()
{
  /* throw all nodes of the tree onto the garbage heap */
  long i;

  gdispose(root);
  for (i = 0; i < nonodes; i++)
    treenode[i] = NULL;
}


void help()
{
  /* display help information */
  char tmp[100];
  printf("\n\nR Rearrange a tree by moving a node or group\n");
  printf(". Redisplay the same tree again\n");
  if (haslengths) {
    printf("= Redisplay the same tree with");
    if (!nolengths)
      printf("out");
    printf(" lengths\n");
  }
  printf("U Undo the most recent rearrangement\n");
  printf("W Write tree to a file\n");
  printf("O select an Outgroup for the tree\n");
  if (haslengths)
    printf("M Midpoint root the tree\n");
  printf("T Transpose immediate branches at a node\n");
  printf("F Flip (rotate) subtree at a node\n");
  printf("B Change or specify the length of a branch\n");
  printf("N Change or specify the name(s) of tip(s)\n");
  printf("H Move viewing window to the left\n");
  printf("J Move viewing window downward\n");
  printf("K Move viewing window upward\n");
  printf("L Move viewing window to the right\n");
  printf("C show only one Clade (subtree) (might be useful if tree is too big)\n");
  printf("+ Read next tree from file\n");
  printf("? Help (this screen)\n");
  printf("Q (Quit) Exit from program\n");
  printf("X Exit from program\n\n\n");
  printf("TO CONTINUE, PRESS ON THE Return OR Enter KEY");
  gets(tmp);
  printree();
}  /* help */

void rearrange()
{
  long i, j;
  boolean ok1, ok2;
  node *p, *q;

  printf("Remove everything to the right of which node? ");
  inpnum(&i, &ok1);
  ok1 = (ok1 && i >= 1 && i < spp * 2 && i != root->index);
  if (ok1) {
    printf("Add before which node? ");
    inpnum(&j, &ok2);
    ok2 = (ok2 && j >= 1 && j < spp * 2);
    if (ok2) {
      ok2 = (treenode[j - 1] != treenode[treenode[i - 1]->back->index - 1]);
      p = treenode[j - 1];
      while (p != root) {
	ok2 = (ok2 && p != treenode[i - 1]);
	p = treenode[p->back->index - 1];
      }
      if (ok1 && ok2) {
	what = i;
	q = treenode[treenode[i - 1]->back->index - 1];
	if (q->next->back->index == i)
	  fromwhere = q->next->next->back->index;
	else
	  fromwhere = q->next->back->index;
	towhere = j;
	re_move(&treenode[i - 1], &q, &restoring,&wasleft);
	add(treenode[j - 1], treenode[i - 1], q, &restoring, &wasleft);
	lastop = rearr;
      }
    }
  }
  printree();
  if (!(ok1 && ok2))
    printf("Not a possible rearrangement.  Try again: \n");
  else {
    oldwritten = written;
    written = false;
  }
}  /* rearrange */

void transpose(atnode)
long atnode;
{
  /* transpose at a node left-right */
  long i;
  boolean ok;
  node *p;
  double templen;

  if (atnode == 0) {
    printf("Transpose branches at which node? ");
    inpnum(&i, &ok);
    ok = (ok && i > spp && i <= nonodes);
  } else {
    i = atnode;
    ok = true;
  }
  if (ok) {
    p = treenode[i - 1]->next->back;
    treenode[i - 1]->next->back = treenode[i - 1]->next->next->back;
    treenode[i - 1]->next->next->back = p;
    treenode[i - 1]->next->back->back = treenode[i - 1]->next;
    treenode[i - 1]->next->next->back->back = treenode[i - 1]->next->next;
    if (haslengths) {
      templen = treenode[i - 1]->next->length;
      treenode[i - 1]->next->length = treenode[i - 1]->next->next->length;
      treenode[i - 1]->next->back->length = treenode[i - 1]->next->next->length;
      treenode[i - 1]->next->next->length = templen;
      treenode[i - 1]->next->next->back->length = templen;
    }
    atwhat = i;
    lastop = transp;
  }
  if (atnode == 0)
    printree();
  if (ok) {
    oldwritten = written;
    written = false;
    return;
  }
  if (i >= 1 && i <= spp)
    printf("Can't transpose there. ");
  else
    printf("No such node. ");
}  /* transpose */

void ltrav__(p)
node *p;
{
  node *temp;
  double templen;

  if (p->tip)
    return;
  temp = p->next->back;
  p->next->back = p->next->next->back;
  p->next->next->back = temp;
  p->next->back->back = p->next;
  p->next->next->back->back = p->next->next;
  if (haslengths) {
    templen = p->next->length;
    p->next->length = p->next->next->length;
    p->next->back->length = p->next->next->length;
    p->next->next->length = templen;
    p->next->next->back->length = templen;
  }
  ltrav__(p->next->back);
  ltrav__(p->next->next->back);
}  /* ltrav */

void flip(atnode)
long atnode;
{
  /* flip at a node left-right */
  long i;
  boolean ok;

  if (atnode == 0) {
    printf("Flip branches at which node? ");
    inpnum(&i, &ok);
    ok = (ok && i > spp && i <= nonodes);
  } else {
    i = atnode;
    ok = true;
  }
  if (ok) {
    ltrav__(treenode[i - 1]);
    atwhat = i;
    lastop = flipp;
  }
  if (atnode == 0)
    printree();
  if (ok) {
    oldwritten = written;
    written = false;
    return;
  }
  if (i >= 1 && i <= spp)
    printf("Can't flip there. ");
  else
    printf("No such node. ");
}  /* flip */

void unmidpoint()
{
  midstoret tempmid;

  tempmid.leftp = root->next->back;
  tempmid.rightp = root->next->next->back;
  tempmid.leftlength = root->next->length;
  tempmid.rightlength = root->next->next->length;
  if (root->next->back != midstore.leftp &&
      root->next->next->back != midstore.leftp)
    reroot(midstore.leftp);
  else if (root->next->back != midstore.rightp &&
	   root->next->next->back != midstore.rightp) {
    reroot(midstore.rightp);
    transpose(root->index);
  }
  root->next->length = midstore.leftlength;
  root->next->back->length = midstore.leftlength;
  root->next->next->length = midstore.rightlength;
  root->next->next->back->length = midstore.rightlength;
  midstore = tempmid;
  lastop = mpoint;
}  /* unmidpoint */

double ltrav___(p)
node *p;
{
  if (p->tip) {
    p->lbeyond = 0.0;
    p->rbeyond = 0.0;
    return 0.0;
  } else {
    p->lbeyond = p->next->length + ltrav___(p->next->back);
    p->rbeyond = p->next->next->length + ltrav___(p->next->next->back);
    if (p->lbeyond > p->rbeyond)
      return (p->lbeyond);
    else
      return (p->rbeyond);
  }
}  /* ltrav */

void outlength()
{
  /* compute the distance to the farthest tip in both directions */
  /* from each node */
  double dummy;

  dummy = ltrav___(root);
}  /* outlength */

void midpoint()
{
  /* Midpoint root the tree */
  double balance;
  node *p, *q;

  midstore.leftp = root->next->back;
  midstore.rightp = root->next->next->back;
  midstore.leftlength = root->next->length;
  midstore.rightlength = root->next->next->length;
  p = root;
  outlength();
  balance = p->lbeyond - (p->lbeyond + p->rbeyond) / 2.0;
  do {
    if (balance > p->next->length) {
      q = p->next->back;
      if (q->lbeyond > q->rbeyond)
	reroot(q->next->back);
      else {
	reroot(q->next->next->back);
	transpose(root->index);
      }
    }
    if (-balance > p->next->next->length) {
      q = p->next->next->back;
      if (q->lbeyond > q->rbeyond)
	reroot(q->next->back);
      else {
	reroot(q->next->next->back);
	transpose(root->index);
      }
    }
    outlength();
    p = root;
    balance = p->lbeyond - (p->lbeyond + p->rbeyond) / 2.0;
  } while (balance > p->next->length || -balance > p->next->next->length);
  p->next->length -= balance;
  p->next->back->length = p->next->length;
  p->next->next->length += balance;
  p->next->next->back->length = p->next->next->length;
  lastop = mpoint;
  printree();
}  /* midpoint */

void undo()
{
  /* restore to tree before last rearrangement */
  long temp;
  boolean btemp;
  node *q;

  switch (lastop) {

  case rearr:
    restoring = true;
    oldleft   = wasleft;
    re_move(&treenode[what - 1], &q, &restoring,&wasleft);
    btemp = wasleft;
    wasleft = oldleft;
    add(treenode[fromwhere - 1], treenode[what - 1],q, &restoring,&wasleft);
    wasleft = btemp;
    restoring = false;
    temp = fromwhere;
    fromwhere = towhere;
    towhere = temp;
    break;

  case flipp:
    flip(atwhat);
    break;

  case transp:
    transpose(atwhat);
    break;

  case reroott:
    restoring = true;
    temp = oldoutgrno;
    oldoutgrno = outgrno;
    outgrno = temp;
    reroot(treenode[outgrno - 1]);
    restoring = false;
    break;

  case mpoint:
    unmidpoint();
    break;

  case lengthh:
    temp = oldls.index;
    templs = oldls;
    oldls.haslength = treenode[temp - 1]->haslength;
    oldls.length = treenode[temp - 1]->length;
    treenode[temp - 1]->length = templs.length;
    treenode[temp - 1]->haslength = templs.haslength;
    treenode[temp - 1]->back->length = templs.length;
    treenode[temp - 1]->back->haslength = templs.haslength;
    break;

  case none:
    /* blank case */
    break;
  }
  printree();
  if (lastop == none) {
    printf("No operation to undo! ");
    return;
  }
  btemp = oldwritten;
  oldwritten = written;
  written = btemp;
}  /* undo */


void treeout(p, writeparens, addlength)
node *p;
boolean writeparens;
double addlength;
{
  /* write out file with representation of final tree */
  long i, n, w;
  Char c;
  double x;

  if (p->tip) {
    if (p->hasname) {
      n = 0;
      for (i = 1; i <= nmlngth; i++) {
	if (treenode[p->index - 1]->name[i - 1] != '\0')
	  n = i;
      }
      for (i = 0; i < n; i++) {
	c = treenode[p->index - 1]->name[i];
	if (c == ' ')
	  c = '_';
	putc(c, outtree);
      }
      col += n;
    }
  } else {
    if (writeparens)
      putc('(', outtree);
    col++;
    treeout(p->next->back, true, 0.0);
    putc(',', outtree);
    col++;
    if (col > 65) {
      putc('\n', outtree);
      col = 0;
    }
    treeout(p->next->next->back, true, 0.0);
    if (writeparens)
      putc(')', outtree);
    col++;
  }
  if (!(p->haslength && writeparens))
    return;
  x = p->length + addlength;
  if (x > 0.0)
    w = (long)(0.43429448222 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.43429448222 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  fprintf(outtree, ":%*.5f", (int)(w + 7), x);
  col += w + 8;
}  /* treeout */

void roottreeout(p,rooted)
node *p;
boolean *rooted;
{
  /* write out file with representation of final tree */
  long trnum, trnumwide;
  double rtlength, leftlength;

  if (nexus) {
    if (!(*rooted))
      putc('U', outtree);
    trnum = treenumber;
    trnumwide = 1;
    while (trnum >= 10) {
      trnum /= 10;
      trnumwide++;
    }
    fprintf(outtree, "TREE PHYLIP_%*ld ", (int)trnumwide, treenumber);
    col += 15;
  }
  if (p->next->back->haslength)
    leftlength = p->next->back->length;
  else
    leftlength = 0.0;
  if (p->next->next->back->haslength)
    rtlength = p->next->next->back->length;
  else
    rtlength = 0.0;
  putc('(', outtree);
  col++;
  if ((*rooted) || (p->next->back->tip && p->next->next->back->tip))
    treeout(p->next->back, true, 0.0);
  else if (!p->next->next->back->tip)
    treeout(p->next->back, true, rtlength);
  else
    treeout(p->next->back, false, 0.0);
  putc(',', outtree);
  col++;
  if (col > 65) {
    putc('\n', outtree);
    col = 0;
  }
  if ((*rooted) || (p->next->back->tip && p->next->next->back->tip))
    treeout(p->next->next->back, true, 0.0);
  else if (!p->next->next->back->tip)
    treeout(p->next->next->back, false, 0.0);
  else
    treeout(p->next->next->back, true, leftlength);
  if (hasmult)
    fprintf(outtree, ")[%6.4f];\n", trweight);
  else
    fprintf(outtree, ");\n");
}  /* roottreeout */

void treewrite(done)
boolean *done;
{
  /* write out tree to a file */
  /* Local variables for treewrite: */
  boolean rooted;

  if (waswritten) {
    printf("\nTree file already was open.\n");
    printf("   A   Add to this tree to tree file\n");
    printf("   R   Replace tree file contents by this tree\n");
    printf("   N   Do Not write out this tree\n");
    do {
      printf("Which should we do? ");
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
	ch = ' ';
      ch = isupper(ch) ? ch : toupper(ch);
    } while (ch != 'A' && ch != 'R' && ch != 'N');
  } else {
    openfile(&outtree,OUTTREE,"w","retree",outtreename);
    if (nexus)
      fprintf(outtree, "BEGIN TREES\n");
    treenumber = 1;
  }
  if (waswritten && ch == 'R') {
    openfile(&outtree,OUTTREE,"w","retree",outtreename);
    if (nexus)
      fprintf(outtree, "BEGIN TREES\n");
    treenumber = 1;
  }
  if (waswritten && ch == 'A')
    openfile(&outtree,OUTTREE,"a","retree",outtreename);
  if (!waswritten || ch == 'A' || ch == 'R') {
    do {
      printf("Enter R if the tree is to be rooted\n");
      printf("OR enter U if the tree is to be unrooted: ");
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
	ch = ' ';
      ch = (isupper(ch)) ? ch : toupper(ch);
    } while (ch != 'R' && ch != 'U');
    col = 0;
    rooted = (ch == 'R');
    roottreeout(root, &rooted);
    treenumber++;
    printf("\nTree written to file\n\n");
    waswritten = true;
    written = true;
  }
  if (!(*done))
    printree();
  fclose(outtree);
}  /* treewrite */

void window(action)
adjwindow action;
{
  /* move viewing window of tree */
  switch (action) {

  case left:
    if (leftedge != 1)
      leftedge -= hscroll;
    break;

  case down:
    /* The 'topedge + 3' is needed to allow downward scrolling
       when part of the tree is above the screen and only 1 or 2 lines
       are below it. */
    if (treelines - topedge + 3 >= screenlines)
      topedge += vscroll;
    break;

  case upp:
    if (topedge != 1)
      topedge -= vscroll;
    break;

  case right:
    if (leftedge < vscreenwidth)
      leftedge += hscroll;
    break;
  }
  printree();
}  /* window */


void getlength(length, reslt, hslngth)
double *length;
reslttype *reslt;
boolean *hslngth;
{
  long digit, ordzero;
  double valyew, divisor;
  Char ch;
  boolean done, pointread;
  char tmp[100];
  valyew = 0.0;
  do {
    printf("\nEnter the new branch length\n");
    printf("OR enter U to leave the length unchanged\n");
    if (*hslngth)
      printf("OR enter R to remove the length from this branch: \n");
    gets(tmp);

    if (tmp[0] == 'u' || tmp[0] == 'U'){
      *reslt = quit;
      break;
    }
    else if (tmp[0] == 'r' || tmp[0] == 'R') {
      (*reslt) = remoov;
      break;}
    else if (sscanf(tmp,"%lf",&valyew) == 1){
      (*reslt) = valid;
      break;}
  } while (1);
  (*length) = valyew;
}  /* getlength */


void changelength()
{
  /* change or specify the length of a tip */
  /* Local variables for changelength: */
  boolean hslngth;
  boolean ok;
  long i, w;
  double length, x;
  Char ch;
  reslttype reslt;
  node *p;

  do {
    printf("Specify length of which branch (0 = all branches)? ");
    inpnum(&i, &ok);
    ok = (ok && (unsigned long)i <= nonodes);
    if (i == 0)
      ok = (treenode[i - 1] != root);
  } while (!ok);
  if (i != 0) {
    p = treenode[i - 1];
    putchar('\n');
    if (p->haslength) {
      x = p->length;
      if (x > 0.0)
	w = (long)(0.43429448222 * log(x));
      else if (x == 0.0)
	w = 0;
      else
	w = (long)(0.43429448222 * log(-x)) + 1;
      if (w < 0)
	w = 0;
      printf("The current length of this branch is %*.5f\n", (int)(w + 7), x);
    } else
      printf("This branch does not have a length\n");
    hslngth = p->haslength;
    getlength(&length, &reslt, &hslngth);
    if (reslt != quit) {
      oldls.length = p->length;
      oldls.haslength = p->haslength;
      oldls.index = p->index;
      lastop = lengthh;
    }
    switch (reslt) {

    case valid:
      p->length = length;
      p->haslength = true;
      if (p->back != NULL) {
	p->back->length = length;
	p->back->haslength = true;
      }
      break;

    case remoov:
      p->haslength = false;
      if (p->back != NULL)
	p->back->haslength = false;
      break;

    case quit:
      /* blank case */
      break;
    }
  } else {
    printf("\n (this operation cannot be undone)\n");
    do {
      printf("\n   enter U to leave the lengths unchanged\n");
      printf("OR enter R to remove the lengths from all branches: \n");
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
	ch = ' ';
    } while (ch != 'U' && ch != 'u' && ch != 'R' && ch != 'r');
    if (ch == 'R' || ch == 'r') {
      for (i = 0; i < spp; i++)
	treenode[i]->haslength = false;
      for (i = spp; i < nonodes; i++) {
	treenode[i]->haslength = false;
	treenode[i]->next->haslength = false;
	treenode[i]->next->next->haslength = false;
      }
    }
  }
  printree();
}  /* changelength */

void changename()
{
  /* change or specify the name of a tip */
  boolean ok;
  long i, n, tipno;
  char tipname[100];

  for(;;) {
   for(;;) {
    printf("Specify name of which tip? (enter its number or 0 to quit): ");
    inpnum(&i, &ok);
   if (((unsigned long)i <= spp) && ok)
       {tipno = i;
        break;};
   }
   if (tipno == 0)
       break;
    if (treenode[tipno - 1]->hasname) {
      n = 0;
      for (i = 1; i <= nmlngth; i++) {
	if (treenode[tipno - 1]->name[i - 1] != '\0')
	  n = i;
      }
      printf("The current name of tip %ld is \"", tipno);
      for (i = 0; i < n; i++)
	putchar(treenode[tipno - 1]->name[i]);
      printf("\"\n");
    }
    for (i = 0; i < nmlngth; i++)
      treenode[tipno - 1]->name[i] = ' ';
    printf("Enter new tip name: ");
    i = 1;
    gets(tipname);
    strncpy(treenode[tipno-1]->name,tipname,nmlngth);
    treenode[tipno - 1]->hasname = true;
    printree();

  }
  printree();
}  /* changename */

void clade()
{
  /* pick a subtree and show only that on screen */
  long i;
  boolean ok;

  printf("Select subtree rooted at which node (0 for whole tree)? ");
  inpnum(&i, &ok);
  ok = (ok && (unsigned long)i <= nonodes);
  if (ok) {
    subtree = (i > 0);
    if (subtree)
      nuroot = treenode[i - 1];
    else
      nuroot = root;
  }
  printree();
  if (!ok)
    printf("Not possible to use this node. ");
}  /* clade */

void changeoutgroup()
{
  long i;
  boolean ok;

  do {
    printf("Which node should be the new outgroup? ");
    inpnum(&i, &ok);
    ok = (ok && i >= 1 && i <= nonodes && i != root->index);
    if (ok)
      outgrno = i;
  } while (!ok);
  reroot(treenode[outgrno - 1]);
  lastop = reroott;
  printree();
  oldwritten = written;
  written = false;
}  /* changeoutgroup */

void redisplay()
{
  boolean done;
  char    ch;

  done = false;
  do {
    printf("\nNEXT? (Options: R . ");
    if (haslengths)
      printf("= ");
    printf("U W O ");
    if (haslengths)
      printf("M ");
    printf("T F B N H J K L C + ? X Q) (? for Help) ");
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    ch = isupper(ch) ? ch : toupper(ch);
    if (ch == 'C' || ch == 'F' || ch == 'O' || ch == 'R' ||
	ch == 'U' || ch == 'X' || ch == 'Q' || ch == '.' ||
	ch == 'W' || ch == 'B' || ch == 'N' || ch == '?' ||
	ch == 'H' || ch == 'J' || ch == 'K' || ch == 'L' ||
	ch == '+' || ch == 'T' || (haslengths && ch == 'M') ||
	(haslengths && ch == '=')) {

      switch (ch) {

      case 'R':
	rearrange();
	break;

      case '.':
	printree();
	break;

      case '=':
	togglelengths();
	break;

      case 'U':
	undo();
	break;

      case 'W':
	treewrite(&done);
	break;

      case 'O':
	changeoutgroup();
	break;

      case 'M':
	midpoint();
	break;

      case 'T':
	transpose(0L);
	break;

      case 'F':
	flip(0L);
	break;

      case 'C':
	clade();
	break;

      case 'B':
	changelength();
	break;

      case 'N':
	changename();
	break;

      case 'H':
	window(left);
	break;

      case 'J':
	window(down);
	break;

      case 'K':
	window(upp);
	break;

      case 'L':
	window(right);
	break;

      case '?':
	help();
	break;

      case '+':
	readnext = true;
	done = true;
	break;

      case 'X':
	done = true;
	break;

      case 'Q':
	done = true;
	break;
      }
    }
  } while (!done);
  if (!written) {
    do {
      printf("Do you want to write out the tree to a file? (Y or N) ");
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
	ch = ' ';
      if (ch == 'Y' || ch == 'y')
	treewrite(&done);
    } while (ch != 'Y' && ch != 'y' && ch != 'N' && ch != 'n');
  }
  if (!readnext && nexus)
    fprintf(outtree, "ENDBLOCK;\n");
}  /* redisplay */


void treeconstruct()
{
  /* constructs a binary tree from the pointers in treenode. */

  waswritten = false;
  readnext = false;
  do {
    restoring = false;
    subtree = false;
    topedge = 1;
    leftedge = 1;
    buildtree();
    readnext = false;
    written = false;
    printree();
    lastop = none;
    redisplay();
    if (readnext)
      unbuildtree();
  } while (readnext);
}  /* treeconstruct */


main(argc, argv)
     int argc;
     Char *argv[];
{
  /* reads in spp. Then calls treeconstruct to                               *
     construct the tree and query the user
                             */

#ifdef MAC
  macsetup("Retree","");
  argv[0] = "Retree";
#endif
  nexus     = false;
  nolengths = false;
  scrollinc = 20;
  screenlines = 24;
  screenwidth = 80;
  vscreenwidth = 80;
  ibmpc = ibmpc0;
  ansi  = ansi0;
  vt52  = vt520;
  getoptions();
  configure();
  treeconstruct();
  if (waswritten)
    FClose(outtree);
  FClose(intree);
  FClose(outtree);
#ifdef MAC
  fixmacfile(intreename);
  fixmacfile(outtreename);
#endif
  exit(0);
}  /* Retree */

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
  mem = (MALLOCRETURN *)calloc(1,x);
  if (!mem)
    memerror();
  else
    return (MALLOCRETURN *)mem;
}


