#include "phylip.h"

#define maxnodes        1200
#define maxnch          30
#define point           '.'
#define minus           '-'
#define stripewidth     3000L
#define maxstripedepth  3500
#define fontsize        3800
#define pi              3.141592653
#define epsilon         0.00001
#define ebcdic          EBCDIC
#define segments        40
#define xstart          10
#define ystart          35
#define LF              10
#define CR              13
#define escape  (ebcdic ?  '\'' :  '\033')
#define null  '\000'


typedef enum {  treepen, labelpen} pentype;
typedef enum { lw,hp,tek,ibmpc,mac,houston, decregis,epson, oki,fig,
                 citoh,toshiba,pcx,pcl,pict,ray,xbm,other} plottertype;

typedef enum {  vertical, horizontal} growth;
typedef enum {cladogram,phenogram,curvogram,eurogram,swoopogram} treestyle;
typedef enum { penup,pendown} pensttstype;

typedef Char plotstring[maxnch];
typedef short fonttype[fontsize];
typedef Char *striparray;
typedef striparray striptype[maxstripedepth];

typedef struct node {
  struct node *next, *back;
  boolean tip;
  plotstring nayme;
  long naymlength, tipsabove, index;
  double xcoord,ycoord,oldlen,length,
         r,theta,oldtheta,width,depth,tipdist,lefttheta,righttheta;
} node;

struct LOC_plottext {              /* Local variables for plottext: */
  double height, compress;
  short *font;
  short coord;
  double heightfont, xfactor, yfactor, xfont, yfont, xplot, yplot, sinslope,
         cosslope, xx, yy;
  pensttstype penstatus;
} ;

typedef struct colortype {
  Char *name;
  double red, green, blue;
} colortype;

double lengthtext();

