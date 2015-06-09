#include "drawgraphics.h"

#ifdef QUICKC
struct videoconfig myscreen;
void   setupgraphics();
#endif

Static colortype colors[7] = {
  {"White    ",0.9,0.9,0.9},
  {"Red      ",1.0,0.3,0.3},
  {"Orange   ",1.0,0.6,0.6},
  {"Yellow   ",1.0,0.9,0.4},
  {"Green    ",0.3,0.8,0.3},
  {"Blue     ",0.5,0.5,1.0},
  {"Violet   ",0.6,0.4,0.8},
};

/*
 * Used ONLY here: */

static long eb[]={
  0 , 1 ,2 ,3 ,55,45,46,47,22,5,37,11,12,13,14,15,16,17,18,19,60,61,50,38,
  24, 25,63,39,28,29,30,31,64,90,127,123,91,108,80,125,77,93,92,78,107,96,
  75,97,240,241,242,243,244,245,246,247,248,249,122,94,76,126,110,111, 124,
  193,194,195,196,197,198,199,200,201,209,210,211, 212,213,214,215,216,217,
  226,227,228,229,230,231,232,233,173,224,189, 95,109,121,129,130,131,132,
  133,134,135,136,137,145,146,147,148,149,150, 151, 152,153,162,163,164,165,
  166,167,168,169,192,79,208,161,7};

long   hpresolution,nmoves,oldpictint,oldx,oldy,bytewrite;
double  labelline,linewidth,oldxhigh,oldxlow,oldyhigh,oldylow,oldxreal,
        oldyreal,raylinewidth,treeline,oldxsize,oldysize,oldxunitspercm,
        oldyunitspercm,oldxcorner,oldycorner;
long rootmatrix[51][51];
long  HiMode,GraphDriver,GraphMode,LoMode;
boolean didenter,didexit,didcompute;

/* externals will move to .h file later. */
extern long         strpbottom,strptop,strpwide,strpdeep,strpdiv;
extern boolean       dotmatrix,empty,preview,previewing;
extern double        expand,xcorner,xnow,xsize,xscale,xunitspercm,
                                      ycorner,ynow,ysize,yscale,yunitspercm,
                                      labelheight,ymargin;
extern long          filesize,bytewrite;
extern growth        grows;
extern enum {yes,no} penchange,oldpenchange;
extern FILE          *plotfile;
extern plottertype   plotter,oldplotter,previewer;
extern striptype     stripe;

extern char pltfilename[100];

void pout(n)
long n;
{
#ifdef MAC
  if (previewing)
    printf("%*ld", (int)((long)(0.434295 * log((double)n) + 0.0001)), n);
  else
    fprintf(plotfile, "%*ld",
	    (int)((long)(0.434295 * log((double)n) + 0.0001)), n);
#else
  if (previewing)
    printf("%*d", (int)((long)(0.434295 * log((double)n) + 0.0001)), n);
  else
    fprintf(plotfile, "%*d",
            (int)((long)(0.434295 * log((double)n) + 0.0001)), n);
#endif
}  /* pout */


long upbyte(num)
long num;
{
  /* get upper nibble of byte */
  long Result, i, j, bytenum, nibcount;
  boolean done;

  bytenum = 0;
  done = false;
  nibcount = 0;
  i = num / 16;
  i /= 16;
  j = 1;
  while (!done) {
    bytenum += (i & 15) * j;
    nibcount++;
    if (nibcount == 2) {
      Result = bytenum;
      done = true;
    } else {
      j *= 16;
      i /= 16;
    }
  }
  return Result;
}  /* upbyte */



Local long lobyte(num)
long num;
{
  /* get low order nibble of byte */
  long Result, i, j, bytenum, nibcount;
  boolean done;

  bytenum = 0;
  done = false;
  nibcount = 0;
  i = num;
  j = 1;
  while (!done) {
    bytenum += (i & 15) * j;
    nibcount++;
    if (nibcount == 2) {
      Result = bytenum;
      done = true;
    } else {
      j *= 16;
      i /= 16;
    }
  }
  return Result;
}  /* lobyte */


void plotdot(ix, iy)
long ix, iy;
{
  /* plot one dot at ix, iy */
  long ix0, iy0, iy1, iy2;

  iy0 = strptop - iy;
  if ((unsigned)iy0 > strpdeep || ix <= 0 || ix > strpwide)
    return;
  empty = false;
  ix0 = ix;
  switch (plotter) {

  case citoh:
    iy1 = 1;
    iy2 = iy0;
    break;

  case epson:
    iy1 = 1;
    iy2 = 7 - iy0;
    break;

  case oki:
    iy1 = 1;
    iy2 = 7 - iy0;
    break;

  case toshiba:
    iy1 = iy0 / 6 + 1;
    iy2 = 5 - iy0 % 6;
    break;

  case pcx:
    iy1 = iy0 + 1;
    ix0 = (ix - 1) / 8 + 1;
    iy2 = 7 - ((ix - 1) & 7);
    break;

  case pcl:
    iy1 = iy0 + 1;
    ix0 = (ix - 1) / 8 + 1;
    iy2 = 7 - ((ix - 1) & 7);
    break;

  case xbm:
    iy1 = iy0 + 1;
    ix0 = (ix - 1) / 8 + 1;
    iy2 = (ix - 1) & 7;
    break;

  case other:
    break;
    /* code for making dot array for a new printer
      goes here */
  }
  stripe[iy1 - 1][ix0 - 1] |= (unsigned char)1<<iy2;
}  /* plotdot */



void drawpen(x, y, width)
long x, y, width;
{
  long i, j, radius, low, hi;

  radius = (long)floor(width / 2.0 + 0.5);
  low = (long)floor(0.5 - width / 2.0);
  hi = width - low;
  if (y + hi < strpbottom || y + low > strptop) {
    if (didenter)
      didexit = true;
    return;
  }
  if (!didenter)
    didenter = true;
  for (i = low; i <= hi; i++) {
    for (j = low; j <= hi; j++) {
      if (rootmatrix[abs(i)][abs(j)] <= radius){
	plotdot(x + i, y + j);
	plotdot(x + i, y - j);
	plotdot(x - i, y + j);
	plotdot(x - 1, y - j);}
    }
  }
}  /* drawpen */




void drawfatline(ixabs, iyabs, ixnow, iynow, penwide)
long ixabs, iyabs, ixnow, iynow, penwide;
{
  long temp, xdiff, ydiff, err, x1, y1;

  didenter = false;
  didexit = false;

  if (ixabs < ixnow) {
    temp = ixnow;
    ixnow = ixabs;
    ixabs = temp;
    temp = iynow;
    iynow = iyabs;
    iyabs = temp;
  }
  xdiff = ixabs - ixnow;
  ydiff = iyabs - iynow;
  if (ydiff >= 0) {
    if (xdiff >= ydiff) {
      err = -(xdiff / 2);
      x1 = ixnow;
      while (x1 <= ixabs && !(didenter && didexit)) {
	drawpen(x1, iynow, penwide);
	err += ydiff;
	if (err > 0) {
	  iynow++;
	  err -= xdiff;
	}
	x1++;
      }
      return;
    }
    err = -(ydiff / 2);
    y1 = iynow;
    while (y1 < iyabs && !(didenter && didexit)) {
      drawpen(ixnow, y1, penwide);
      err += xdiff;
      if (err > 0) {
	ixnow++;
	err -= ydiff;
      }
      y1++;
    }
    return;
  }
  if (xdiff < -ydiff) {
    err = ydiff / 2;
    y1 = iynow;
    while (y1 >= iyabs && !(didenter && didexit)) {
      drawpen(ixnow, y1, penwide);
      err += xdiff;
      if (err > 0) {
	ixnow++;
	err += ydiff;
      }
      y1--;
    }
    return;
  }
  err = -(xdiff / 2);
  x1 = ixnow;
  while (x1 <= ixabs && !(didenter && didexit)) {
    drawpen(x1, iynow, penwide);
    err -= ydiff;
    if (err > 0) {
      iynow--;
      err -= xdiff;
    }
    x1++;
  }
}  /* drawfatline */


void plot(pen, xabs, yabs)
pensttstype pen;
double xabs, yabs;
{
  long xhigh, yhigh, xlow, ylow, newx, newy, ixnow, iynow, ixabs, iyabs,
       ixleft, iyleft, ixbot, iybot, ixtop, iytop, cdx, cdy, temp;
  long pictint;
  double c, dx, dy, lscale, dxreal, dyreal;
  Char picthi, pictlo;
#ifdef MAC
  queryevent();
#endif
  if (!dotmatrix || previewing) {
    switch (plotter) {

    case tek:
      if (pen == penup) {
        if (previewing)
          putchar('\035');
        else
          putc('\035', plotfile);
      }
      ixnow = (long)floor(xabs + 0.5);
      iynow = (long)floor(yabs + 0.5);
      xhigh = ixnow / 32;
      yhigh = iynow / 32;
      xlow = ixnow & 31;
      ylow = iynow & 31;
      if (!ebcdic) {
        if (yhigh != oldyhigh) {
          if (previewing)
            putchar(yhigh + 32);
          else
            putc(yhigh + 32, plotfile);
        }
        if (ylow != oldylow || xhigh != oldxhigh) {
          if (previewing)
            putchar(ylow + 96);
          else
            putc(ylow + 96, plotfile);
        }
        if (xhigh != oldxhigh) {
          if (previewing)
            putchar(xhigh + 32);
          else
            putc(xhigh + 32, plotfile);
        }
        if (previewing)
          putchar(xlow + 64);
        else
          putc(xlow + 64, plotfile);
      } else {  /* DLS/JMH -- for systems that use EBCDIC coding */
        if (yhigh != oldyhigh) {
          if (previewing)
            putchar(eb[yhigh + 32]);
          else
            putc(eb[yhigh + 32], plotfile);
        }
        if (ylow != oldylow || xhigh != oldxhigh) {
          if (previewing)
            putchar(eb[ylow + 96]);
          else
            putc(eb[ylow + 96], plotfile);
        }
        if (xhigh != oldxhigh) {
          if (previewing)
            putchar(eb[xhigh + 32]);
          else
            putc(eb[xhigh + 32], plotfile);
        }
        if (previewing)
          putchar(eb[xlow + 64]);
        else
          putc(eb[xlow + 64], plotfile);
      }

      oldxhigh = xhigh;
      oldxlow = xlow;
      oldyhigh = yhigh;
      oldylow = ylow;
      break;

    case hp:
      if (pen == pendown)
        fprintf(plotfile, "PD");
      else
        fprintf(plotfile, "PU");
      pout((long)floor(xabs + 0.5));
      putc(',', plotfile);
      pout((long)floor(yabs + 0.5));
      fprintf(plotfile, ";\n");
      break;

    case pict:
      newx = (long)floor(xabs + 0.5);
      newy = (long)floor(ysize * yunitspercm - yabs + 0.5);
      if (pen == pendown) {
        if (linewidth > 5) {
          dxreal = xabs - oldxreal;
          dyreal = yabs - oldyreal;
          lscale = sqrt(dxreal * dxreal + dyreal * dyreal) /
            (fabs(dxreal) + fabs(dyreal));
          pictint = (long)(lscale * linewidth + 0.5);

          if (pictint == 0)
            pictint = 1;
          if (pictint != oldpictint) {
            picthi = (Char)(pictint / 256);
            pictlo = (Char)(pictint & 255);
            fprintf(plotfile, "\007%c%c%c%c", picthi, pictlo, picthi, pictlo);
          }
          oldpictint = pictint;
        }
        fprintf(plotfile, " %c%c%c%c",
                (Char)(oldy / 256), (Char)(oldy & 255), (Char)(oldx / 256),
                (Char)(oldx & 255));
        fprintf(plotfile, "%c%c%c%c",
                (Char)(newy / 256), (Char)(newy & 255), (Char)(newx / 256),
                (Char)(newx & 255));
        }
      oldxreal = xabs;
      oldyreal = yabs;
      oldx = newx;
      oldy = newy;

      break;
#ifndef MAC
    case ray:
      if (pen == pendown) {
        if (linewidth != treeline) {
          if (raylinewidth > labelline) {
            raylinewidth = labelline;
            fprintf(plotfile, "end\n\n");
            fprintf(plotfile, "name species_names\n");
            fprintf(plotfile, "grid 22 22 22\n");
          }
        }

        if (oldxreal != xabs || oldyreal != yabs) {
          raylinewidth *= 0.99999;
          fprintf(plotfile, "cylinder %8.7f %6.3f 0 %6.3f %6.3f 0 %6.3f\n",
                  raylinewidth, oldxreal, oldyreal, xabs, yabs);
          fprintf(plotfile, "sphere %8.7f %6.3f 0 %6.3f\n",
                  raylinewidth, xabs, yabs);
        }
      }
      oldxreal = xabs;
      oldyreal = yabs;
      break;
#endif /* ifndef MAC */
    case lw:
      if (pen == pendown)
        fprintf(plotfile, "%8.2f%8.2f lineto\n", xabs, yabs);
      else
        fprintf(plotfile, "stroke %8.2f%8.2f moveto\n", xabs, yabs);
      break;

    case ibmpc:
#ifdef TURBOC
    newx =  (long)(floor(xabs + 0.5));
    newy = abs((long)(floor(yabs)) - getmaxy());
    if (pen == pendown)
        line(oldx,oldy,newx,newy);
    oldx=newx;
    oldy=newy;
#endif
#ifdef QUICKC
    newx =  (long)(floor(xabs + 0.5));
    newy = abs((long)(floor(yabs)) - myscreen.numypixels);

    if (pen == pendown)
        _lineto((long)newx,(long)newy);
    else
        _moveto((long)newx,(long)newy);
    oldx=newx;
    oldy=newy;

#endif
    break;
     case mac:
#ifdef MAC
      if (pen == pendown){
	LineTo((int)floor((double)xabs + 0.5),
	       342 - (long)floor((double)yabs + 0.5));}
      else{
	MoveTo((int)floor((double)xabs + 0.5),
	       342 - (long)floor((double)yabs + 0.5));}
#endif

      break;

    case houston:
      if (pen == pendown)
        fprintf(plotfile, "D ");
      else
        fprintf(plotfile, "U ");
      pout((long)((long)floor(xabs + 0.5)));
      putc(',', plotfile);
      pout((long)((long)floor(yabs + 0.5)));
      putc('\n', plotfile);
      break;

    case decregis:
      newx = (long)floor(xabs + 0.5);
      newy = (long)abs((long)floor(yabs + 0.5) - 479);
      if (pen == pendown) {
        if (previewing) {
          printf("P[");
          pout(oldx);
          putchar(',');
          pout(oldy);
          printf("]V[");
          pout(newx);
          putchar(',');
          pout(newy);
          putchar(']');
        } else {
          fprintf(plotfile, "P[");
          pout(oldx);
          putc(',', plotfile);
          pout(oldy);
          fprintf(plotfile, "]V[");
          pout(newx);
          putc(',', plotfile);
          pout(newy);
          putc(']', plotfile);
        }
        nmoves++;
        if (nmoves == 3) {
          nmoves = 0;
          if (previewing)
            putchar('\n');
          else
            putc('\n', plotfile);
        }
      }
      oldx = newx;
      oldy = newy;
      break;

    case fig:
      newx = (long)floor(xabs + 0.5);
      newy = (long)floor(yabs + 0.5);
      if (pen == pendown) {
	fprintf(plotfile, "2 1 0 %5ld 0 0 0 0 0.000 0 0\n",
		(long)floor(linewidth + 0.5) + 1);
	fprintf(plotfile, "%5ld%5ld%5ld%5ld 9999 9999\n",
		oldx, 606 - oldy, newx, 606 - newy);

	fprintf(plotfile,
	  "1 3 0  1 0 0 0 21 0.00 1 0.0 %5ld%5ld%5ld %5ld %5ld%5ld%5ld 349\n",
	  oldx, 606 - oldy, (long)floor(linewidth / 2 + 0.5),
	  (long)floor(linewidth / 2 + 0.5), oldx, 606 - oldy, 606 - oldy);
	fprintf(plotfile,
	  "1 3 0  1 0 0 0 21 0.00 1 0.0 %5ld%5ld%5ld %5ld %5ld%5ld%5ld 349\n",
	  newx, 606 - newy, (long)floor(linewidth / 2 + 0.5),
	  (long)floor(linewidth / 2 + 0.5), newx, 606 - newy, 606 - newy);

      }
      oldx = newx;
      oldy = newy;
      break;
    case other:
      break;
      /* code for a pen move on a new plotter goes here */
    }
    return;
  }
  if (pen == pendown) {
    ixabs = (long)floor(xabs + 0.5);
    iyabs = (long)floor(yabs + 0.5);
    ixnow = (long)floor(xnow + 0.5);
    iynow = (long)floor(ynow + 0.5);
    if (ixnow > ixabs) {
      temp = ixnow;
      ixnow = ixabs;
      ixabs = temp;
      temp = iynow;
      iynow = iyabs;
      iyabs = temp;
    }
    dx = ixabs - ixnow;
    dy = iyabs - iynow;
   /* if (dx + fabs(dy) <= 0.0)
      c = 0.0;
    else
      c = 0.5 * linewidth / sqrt(dx * dx + dy * dy); */
    cdx = (long)floor(linewidth + 0.5);
    cdy = (long)floor(linewidth + 0.5);
    if ((iyabs + cdx >= strpbottom || iynow + cdx >= strpbottom) &&
        (iyabs - cdx <= strptop || iynow - cdx <= strptop)) {
      drawfatline(ixnow,iynow,ixabs,iyabs,(long)floor(linewidth+0.5));}
  }

  xnow = xabs;
  ynow = yabs;

  /* Bitmap Code to plot (xnow,ynow) to (xabs,yabs)                 */
} /* plot                                                           */


void pictoutint(file,pictint)
FILE *file;
long  pictint;
{
char picthi, pictlo;

picthi = (char)(pictint / 256);
pictlo = (char)(pictint % 256);
fprintf(file, "%c%c", picthi, pictlo);
}

void initplotter(ntips,fontname)
long ntips;
char *fontname;
{
  long i,j, hres, vres;
  Char picthi, pictlo;
  long pictint;

  treeline = 0.18 * labelheight * yscale * expand;
  labelline = 0.06 * labelheight * yscale * expand;
  linewidth = treeline;
  if (dotmatrix ) {
    for (i = 0; i <= 50; i++) {   /* for fast circle calculations */
     for (j = 0; j <= 50; j++){
       rootmatrix[i][j] =
	   (long)floor(sqrt((double)(i * i + j * j)) + 0.5);}
   }
  }
 switch (plotter) {

  case tek:
    oldxhigh = -1.0;
    oldxlow = -1.0;
    oldyhigh = -1.0;
    oldylow = -1.0;
    nmoves = 0;       /* DLS/JMH -- See function  PLOT                  */
    if (previewing)   /* DLS/JMH                                        */
      printf("%c\f", escape);   /* DLS/JMH */
    else
      fprintf(plotfile, "%c\f", escape);
    break;

  case hp:
    fprintf(plotfile, "IN;SP1;VS10.0;\n");
    break;
#ifndef MAC
  case ray:
    treeline = 0.27 * labelheight * yscale * expand;
    linewidth = treeline;
    raylinewidth = treeline;
    if (grows == vertical)
      fprintf(plotfile, "plane backcolor 0 0 %2.4f 0 0 1\n", ymargin);
    else
      fprintf(plotfile, "plane backcolor 0 0 %2.4f 0 0 1\n",
              ymargin - ysize / (ntips - 1));

    fprintf(plotfile, "\nname tree\n");
    fprintf(plotfile, "grid 22 22 22\n");
    break;
#endif /* ifndef mac */
   case pict:
    plotfile = freopen(pltfilename,"wb",plotfile);
    for (i=0;i<512;++i)
      putc('\000',plotfile);
    pictoutint(plotfile,1000); /* size...replaced later with seek */
    pictoutint(plotfile,1);    /* bbx0   */
    pictoutint(plotfile,1);    /* bby0   */
    pictoutint(plotfile,612);  /* bbx1   */
    pictoutint(plotfile,792);  /* bby1   */
    fprintf(plotfile,"%c%c",0x11,0x01); /* version "1" (B&W) pict */
    fprintf(plotfile,"%c%c%c",0xa0,0x00,0x82);
    fprintf(plotfile,"%c",1);    /* clip rect */
    pictoutint(plotfile,10);  /* region size, bytes. */
    pictoutint(plotfile,1);   /* clip x0             */
    pictoutint(plotfile,1);   /* clip y0             */
    pictoutint(plotfile,612); /* clip x1             */
    pictoutint(plotfile,792); /* clip y1             */
    
    bytewrite=543;

    oldpictint = 0;
    pictint = (long)(linewidth + 0.5);
    if (pictint == 0)
      pictint = 1;
    picthi = (Char)(pictint / 256);
    pictlo = (Char)(pictint % 256);
    fprintf(plotfile, "\007%c%c%c%c", picthi, pictlo, picthi, pictlo);
    /* Set pen size for drawing tree. */
    break;

  case xbm:  /* what a completely verbose data representation format.*/

    fprintf(plotfile, "#define drawgram_width %5ld\n",
            (long)(xunitspercm * xsize));
    fprintf(plotfile, "#define drawgram_height %5ld\n",
            (long)(yunitspercm * ysize));
    fprintf(plotfile, "static char drawgram_bits[] = {\n");
    /*filesize := 53;  */
    break;
  case lw:     /* write conforming encapsulated postscript */
    fprintf(plotfile, "%%!PS-Adobe-2.0 EPSF-2.0\n");
    fprintf(plotfile,"%%%%Creator: Phylip\n");
    fprintf(plotfile,"%%%%Title:  phylip.ps\n%%%%Pages: 1\n");
    fprintf(plotfile,"%%%%BoundingBox: 0 0 612 792\n");
    fprintf(plotfile,"%%%%EndComments\n%%%%EndProlog\n%%%%Page: 1 1\n\n");
    fprintf(plotfile," 1 setlinecap \n 1 setlinejoin  \n");
    fprintf(plotfile, "%8.2f setlinewidth newpath \n", treeline);
    break;

  case ibmpc:
#ifdef TURBOC
    initgraph(&GraphDriver,&HiMode,"");
#endif
#ifdef QUICKC
    setupgraphics();
#endif
    break;

  case mac:
#ifdef MAC
    gfxmode();
    pictint=(long)(linewidth + 0.5);
    if (pictint == 0)
         pictint=1;
    PenSize((int)pictint,(int)pictint);
#endif
    break;

  case houston:
    break;

  case decregis:
    oldx = 300;
    oldy = 1;
    nmoves = 0;
    if (previewing)
	  printf("%c[2J%cPpW(I3);SA[0,0][799,479];S(I(W))S(E);S(C00;W(I(D))\n",
		 escape,escape);
     else
       fprintf(plotfile,
	       "%c[2J%cPpW(I3);S(A[0,0][799,479]);S(I(W))S(E);S(C0);W(I(D))\n",
	       escape,escape);
    break;

  case epson:
    plotfile = freopen(pltfilename,"wb",plotfile);
    fprintf(plotfile, "\0333\030");
    break;

  case oki:
    plotfile = freopen(pltfilename,"wb",plotfile);
    fprintf(plotfile, "\033%%9\020");
    break;

  case citoh:
    plotfile = freopen(pltfilename,"wb",plotfile);
    fprintf(plotfile, "\033T16");
    break;

  case toshiba: /* reopen in binary since we always need \n\r on the file */
                /* and dos in text mode puts it, but unix does not        */
    plotfile = freopen(pltfilename,"wb",plotfile);
    fprintf(plotfile, "\033\032I\n\r\n\r");
    fprintf(plotfile, "\033L06\n\r");
    break;

  case pcl:
    plotfile = freopen(pltfilename,"wb",plotfile);
    /* fprintf(plotfile, "\033&f0S");   Push current cursor */
    if (hpresolution == 150 || hpresolution == 300)
      fprintf(plotfile, "\033*t%3ldR", hpresolution);
    else if (hpresolution == 75)
      fprintf(plotfile, "\033*t75R");
    break;

  case pcx:
    plotfile = freopen(pltfilename,"wb",plotfile);
    fprintf(plotfile,"\012\003\001\001%c%c%c%c",0,0,0,0);
  /* Manufacturer version (1 byte) version (1 byte), encoding (1 byte),
     bits per pixel (1 byte), xmin (2 bytes) ymin (2 bytes),
     Version */
    hres = strpwide;
    vres = (long)floor(yunitspercm * ysize + 0.5);
    fprintf(plotfile, "%c%c", (unsigned char)lobyte(hres - 1),
	    (unsigned char)upbyte(hres - 1)); /* Xmax */
    fprintf(plotfile, "%c%c", (unsigned char)lobyte(vres - 1),
	    (unsigned char)upbyte(vres - 1)); /* Ymax */
    fprintf(plotfile, "%c%c", (unsigned char)lobyte(hres),
	    (unsigned char)upbyte(hres));
    /* Horizontal resolution */
    fprintf(plotfile, "%c%c", (unsigned char)lobyte(vres),
	    (unsigned char)upbyte(vres));
    /* Vertical resolution */
    for (i = 1; i <= 48; i++)  /* fill color map with 0 */
      putc('\000', plotfile);
    putc('\000', plotfile);
    putc('\001', plotfile);   /* Num Planes */
    putc(hres / 8, plotfile);   /* Bytes per line */
    putc('\000',plotfile);
    for (i = 1; i <= 60; i++)   /* Filler */
      putc('\000',plotfile);
    break;
  case fig:
    fprintf(plotfile, "#FIG 2.0\n");
    fprintf(plotfile, "80 2\n");
    break;
  case other:
    break;
    /* initialization code for a new plotter goes here */
  }
}  /* initplotter */




void finishplotter()
{
char trash;
  switch (plotter) {

  case tek:
    if (previewing) {
      scanf("%c%*[^\n]", &trash);
      trash=getchar();
      printf("%c\f", escape);
    } else {
      putc('\n', plotfile);
      plot(penup, 1.0, 1.0);
    }
    break;

  case hp:
    plot(penup, 1.0, 1.0);
    fprintf(plotfile, "SP;\n");
    break;

#ifndef MAC
  case ray:
    fprintf(plotfile,"end\n\nobject treecolor tree\n");
    fprintf(plotfile,"object namecolor species_names\n");
    break;
#endif

  case pict:
    fprintf(plotfile,"%c%c%c%c%c",0xa0,0x00,0x82,0xff,0x00);
    bytewrite+=5;
#ifndef SEEK_SET
#define SEEK_SET 0
#endif
    fseek(plotfile,512L,SEEK_SET);
    pictoutint(plotfile,bytewrite);
    break;


  case lw:
    fprintf(plotfile, "stroke showpage \n\n");
    break;

  case ibmpc:
#ifdef TURBOC
    trash=getchar();
    restorecrtmode();
#endif
#ifdef QUICKC
    trash=getchar();
    _clearscreen(_GCLEARSCREEN);
    _setvideomode(_DEFAULTMODE);
#endif
    break;

  case mac:
#ifdef MAC
    if (previewing) {
      scanf("%c%*[^\n]", &trash);
      trash=getchar();
      textmode();}
#endif
    break;

  case houston:
    break;

  case decregis:
    plot(penup, 1.0, 1.0);
    if (previewing)
      printf("%c\\", escape);
    else
      fprintf(plotfile, "%c\\", escape);
    if (previewing) {
      trash = getchar();
      printf("%c[2J",escape);
    }
    
    break;

  case epson:
    fprintf(plotfile, "\0333$");
    break;

  case oki:
    /* blank case */
    break;

  case citoh:
    fprintf(plotfile, "\033A");
    break;

  case toshiba:
    fprintf(plotfile, "\033\032I\n\r");
    break;

  case pcl:
    /* fprintf(plotfile, "\033&f1S");   pop cursor         */
    fprintf(plotfile, "\033*rB");    /* Exit graphics mode */
    putc('\f', plotfile);            /* just to make sure? */
    break;

  case pcx:
    /* blank case */
    break;

  case xbm:
    fprintf(plotfile, "}\n");
    break;

  case fig:
    /* blank case */
    break;
  case other:
    break;
    /* termination code for a new plotter goes here */
  }
}  /* finishplotter */


Local long SFactor()
{
  /* the dot-skip is resolution-independent. */
  /* this makes all the point-skip instructions skip the same # of dots. */
  long Result;

  if (hpresolution == 150)
    Result = 2;
  if (hpresolution == 300)
    Result = 1;
  if (hpresolution == 75)
    return 4;
  return Result;
}  /* SFactor */

long DigitsInt(x)
long x;
{
  if (x < 10)
    return 1;
  else if (x >= 10 && x < 100)
    return 2;
  else
    return 3;
}  /* DigistInt */

Local boolean IsColumnEmpty(mystripe, pos,deep)
striparray *mystripe;
 long pos,deep;
{
  long j;
  boolean ok;

  ok = true;
  j = 1;
  while (ok && j <= deep) {
    ok = (ok && mystripe[j - 1][pos - 1] == null);
    j++;
  }
  return ok;
}  /* IsColumnEmpty */

void Skip(Amount)
 long Amount;
{
  /* assume we're not in gfx mode. */
  fprintf(plotfile, "\033&f1S");   /* Pop the graphics cursor    */
#ifdef MAC
  fprintf(plotfile, "\033*p+%*ldX",
	  (int)DigitsInt(Amount * SFactor()), Amount * SFactor());
#else
  fprintf(plotfile, "\033*p+%*dX",
	  (int)DigitsInt(Amount * SFactor()), Amount * SFactor());
#endif
  fprintf(plotfile, "\033&f0S");   /* Push the cursor to new location */
  filesize += 15 + DigitsInt(Amount * SFactor());
}  /* Skip */

Local long FirstBlack(mystripe, startpos,deep)
striparray *mystripe;
 long startpos,deep;
{
  /* returns, given a strip and a position, next x with some y's nonzero */
  long i;
  boolean columnempty;

  i = startpos;
  columnempty = true;
  while (columnempty && i < strpwide / 8) {
    columnempty = (columnempty && IsColumnEmpty(mystripe, i,deep));
    if (columnempty)
      i++;
  }
  return i;
}  /* FirstBlack */

Local long FirstWhite(mystripe, startpos,deep)
striparray *mystripe;
 long startpos,deep;
{
  /* returns, given a strip and a position, the next x with all y's zero */
  long i;
  boolean columnempty;

  i = startpos;
  columnempty = false;
  while (!columnempty && i < strpwide / 8) {
    columnempty = IsColumnEmpty(mystripe, i,deep);
    if (!columnempty)
      i++;
  }
  return i;
}  /* FirstWhite */

Local boolean IsBlankStrip(mystripe,deep)
striparray *mystripe;
 long deep;
{
  long i, j;
  boolean ok;

  ok = true;
  i = 1;
  while (ok && i <= strpwide / 8) {
    for (j = 0; j < (deep); j++)
      ok = (ok && mystripe[j][i - 1] == '\0');
    i++;
  }
  return ok;
}  /* IsBlankStrip */


void striprint(div,deep)
 long div,deep;
{
  long i, j, t, x, theend, width;
  Char counter;
  boolean done;
  done = false;
  width = strpwide;
  if (plotter != pcx && plotter != pcl) {
    while (!done) {
      for (i = 0; i < div; i++)
        done = (done || (stripe[i] && stripe[i][width - 1] != null));
      if (!done)
        width--;
      done = (done || width == 0);
    }
  }
  switch (plotter) {

  case epson:
    if (!empty) {
      fprintf(plotfile, "\033L%c%c", width & 255, width / 256);
      for (i = 0; i < width; i++)
        putc(stripe[0][i], plotfile);
      filesize += width + 4;
    }
    putc('\n', plotfile);
    putc('\r', plotfile);
    break;

  case oki:
    if (!empty) {
      fprintf(plotfile, "\033%%1%c%c", width / 128, width & 127);
      for (i = 0; i < width; i++)
        putc(stripe[0][i], plotfile);
      filesize += width + 5;
    }
    putc('\n', plotfile);
    putc('\r', plotfile);
    break;

  case citoh:
    if (!empty) {
      fprintf(plotfile, "\033S%04ld",width);
      for (i = 0; i < width; i++)
        putc(stripe[0][i], plotfile);
      filesize += width + 6;
    }

    putc('\n', plotfile);
    putc('\r', plotfile);
    break;

  case toshiba:
    if (!empty) {
      for (i = 0; i < width; i++) {
        for (j = 0; j <= 3; j++)
          stripe[j][i] += 64;
      }
      fprintf(plotfile, "\033;%04ld",width);

      for (i = 0; i < width; i++)
        fprintf(plotfile, "%c%c%c%c",
                stripe[0][i], stripe[1][i], stripe[2][i], stripe[3][i]);
      filesize += width * 4 + 6;
    }
    putc('\n', plotfile);
    putc('\r', plotfile);
    break;

  case pcx:
    width = strpwide / 8;
    for (j = 0; j < div; j++) {
      t = 1;
      while (1) {
        i = 0; /* i == RLE count ???? */
        while ((stripe[j][t + i - 1]) == (stripe[j][t + i])
               && t + i < width && i < 63)
          i++;
        if (i > 0) {
          counter = 192;
          counter += i;
          putc(counter, plotfile);
          putc(255 - stripe[j][t - 1], plotfile);
          t += i;
          filesize += 2;
        } else {
          if (255 - (stripe[j][t - 1] & 255) >= 192) {
            putc(193, plotfile);
            filesize++;
          }
          putc(255 - stripe[j][t - 1], plotfile);
          t++;
          filesize++;

        }
        if (t >width) break;
      }
    }
    break;

  case pcl:
    width = strpwide / 8;
    if (IsBlankStrip(stripe,deep)) {
#ifdef MAC
      fprintf(plotfile, "\033&f1S\033*p0X\033*p+%*ldY\033&f0S",
	      (int)DigitsInt(deep * SFactor()), deep * SFactor());
#else
      fprintf(plotfile, "\033&f1S\033*p0X\033*p+%*dY\033&f0S",
	      (int)DigitsInt(deep * SFactor()), deep * SFactor());
#endif
      filesize += 20 + DigitsInt(deep * SFactor());
    } else {  /* plotting the actual strip as bitmap data */
      x = 1;
      theend = 1;
      while (x < width) {
        x = FirstBlack(stripe, x,deep);    /* all-black strip is now    */
        Skip((x - theend - 1) * 8);        /* x..theend                 */
        theend = FirstWhite(stripe, x,deep) - 1;/* like lastblack            */
        fprintf(plotfile, "\033*r1A");     /* enter gfx mode            */
        for (j = 0; j < div; j++) {
#ifdef MAC
	  fprintf(plotfile, "\033*b%*ldW",
		  (int)DigitsInt(theend - x + 1), theend - x + 1);
#else
	  fprintf(plotfile, "\033*b%*dW",
                  (int)DigitsInt(theend - x + 1), theend - x + 1);
#endif
              /* dump theend-x+1 bytes */
          for (t = x - 1; t < theend; t++)
            putc(stripe[j][t], plotfile);
          filesize += theend - x + DigitsInt(theend - x + 1) + 5;
        }
        fprintf(plotfile, "\033*rB");   /* end gfx mode */
        Skip((theend - x + 1) * 8);
        filesize += 9;
        x = theend + 1;
      }
      fprintf(plotfile, "\033&f1S");   /* Pop cursor  */
#ifdef MAC
      fprintf(plotfile, "\033*p0X\033*p+%*ldY",
	      (int)DigitsInt(deep * SFactor()), deep * SFactor());
#else
      fprintf(plotfile, "\033*p0X\033*p+%*dY",
	      (int)DigitsInt(deep * SFactor()), deep * SFactor());
#endif
      filesize += 20 + DigitsInt(deep * SFactor());
      fprintf(plotfile, "\033&f0S");   /* Push cursor  */
    }
    break;
    /* case for hpcl code */
  case xbm:
    x = 0;   /* count up # of bytes so we can put returns. */
    width = ((strpwide -1) / 8) +1;
    for (j = 0; j <  div; j++) {
      for (i = 0; i < width; i++) {
        fprintf(plotfile, "0x%02x,",(unsigned char)stripe[j][i]);
        filesize += 5;
        x++;
        if ((x % 15) == 0) {
          putc('\n', plotfile);
          filesize++;
        }
      }
    }
   putc('\n',plotfile);
   break;

  case other:
    break;
    /* graphics print code for a new printer goes here */
  }
 }  /* striprint */

#ifdef QUICKC
void setupgraphics()
{
_getvideoconfig(&myscreen);
#ifndef WATCOM
switch(myscreen.adapter){
  case _CGA:
  case _OCGA:
   _setvideomode(_HRESBW);
    break;
  case _EGA:
  case _OEGA:
    _setvideomode(_ERESNOCOLOR);
  case _VGA:
  case _OVGA:
  case _MCGA:
    _setvideomode(_VRES2COLOR);
     break;
  case _HGC:
    _setvideomode(_HERCMONO);
     break;
  default:
     printf("Your display hardware is unsupported by this program.\n");
      break;
}
#else
switch(myscreen.adapter){
  case _VGA:
  case _SVGA:
      _setvideomode(_VRES16COLOR);
      break;
  case _MCGA:
      _setvideomode(_MRES256COLOR);
      break;
  case _EGA:
     _setvideomode(_ERESNOCOLOR);
     break;
  case _CGA:
     _setvideomode(_MRES4COLOR);
     break;
  case _HERCULES:
     _setvideomode(_HERCMONO);
     break;
  default:
     printf("Your display hardware is unsupported by this program.\n");
     exit(-1);
     break;
   }
#endif
_getvideoconfig(&myscreen);
_setlinestyle(0xffff);
xunitspercm=myscreen.numxpixels / 25;
yunitspercm=myscreen.numypixels / 17.5;
xsize = 25.0;
ysize = 17.5;
}
#endif

void loadfont(font,application)
short *font;
char *application;
{

FILE *fontfile;
 long i, charstart, dummy;
Char trash,ch = 'A';
  i=0;
  openfile(&fontfile,FONTFILE,"r",application,NULL);

  while (!(eof(fontfile) || ch == ' ')) {
    charstart = i + 1;
    fscanf(fontfile, "%c%c%hd%hd%hd", &ch, &ch, &dummy, &font[charstart + 1],
           &font[charstart + 2]);
    font[charstart] = ch;
    i = charstart + 3;
    do {
      if ((i - charstart - 3) % 10 == 0) {
        fscanf(fontfile, "%*[^\n]");
        getc(fontfile);
      }
      i++;
      fscanf(fontfile, "%hd", &font[i - 1]);
    } while (abs(font[i - 1]) < 10000);
    fscanf(fontfile, "%*[^\n]");
#ifdef MAC
    queryevent();
#endif
    getc(fontfile);
    font[charstart - 1] = i + 1;
  }
  font[charstart - 1] = 0;
 FClose(fontfile);
}  /* loadfont */

#ifndef MAC

 long  showrayparms(treecolor,namecolor,backcolor,rx,ry)
 long treecolor,namecolor,backcolor,rx,ry;
{
  long i;
  Char ch,input[32];
  long numtochange;

  if (previewer == tek)
    printf("%c\f", escape);
  else {
    for (i = 1; i <= 24; i++)
      putchar('\n');
  }
  printf("Settings for Rayshade file: \n\n");
  printf(" (1)               Tree color:  %.10s\n",colors[treecolor-1].name);
  printf(" (2)      Species names color:  %.10s\n",colors[namecolor-1].name);
  printf(" (3)         Background color:  %.10s\n",colors[backcolor-1].name);
  printf(" (4)               Resolution:  %2ld X %2ld\n\n",rx,ry);

  printf(" Do you want to accept these? (Yes or No)\n");
  for (;;) {
    printf(" Type Y or N or the number (1-4) of the one to change: \n");
    gets(input);
    numtochange=atoi(input);
    uppercase(&input[0]);
    ch=input[0];
    if (ch == 'Y' || ch == 'N' || (numtochange >= 1 && numtochange <= 4))
      break;
  }
 return (ch == 'Y') ? -1 : numtochange;
}  /* showrayparms */


void getrayparms(treecolor,namecolor,backcolor,rx,ry,numtochange)
 long *treecolor,*namecolor,*backcolor,*rx,*ry;
 long numtochange;
{
  Char ch;
  long i;

  if (numtochange == 0) {
    do {
      printf(" Type the number of one that you want to change (1-4):\n");
      scanf("%ld%*[^\n]", &numtochange);
      getchar();
    } while (numtochange < 1 || numtochange > 10);
  }
  switch (numtochange) {

  case 1:
    printf("\nWhich of these colors will the tree be?:\n");
    printf("   White, Red, Orange, Yellow, Green, Blue, or Violet\n");
    printf(" (W, R, O, Y, G, B, or V)\n");
    do {
      printf(" Choose one: \n");
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
        ch = ' ';
      uppercase(&ch);
      (*treecolor) = 0;
      for (i = 1; i <= 7; i++) {
        if (ch == colors[i - 1].name[0]) {
          (*treecolor) = i;
          return;
        }
      }
    } while ((*treecolor) == 0);
    break;

  case 2:
    printf("\nWhich of these colors will the species names be?:\n");
    printf("   White, Red, Orange, Yellow, Green, Blue, or Violet\n");
    printf(" (W, R, O, Y, G, B, or V)\n");
    do {
      printf(" Choose one: \n");
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
        ch = ' ';
      uppercase(&ch);
      (*namecolor) = 0;
      for (i = 1; i <= 7; i++) {
        if (ch == colors[i - 1].name[0]) {
          (*namecolor) = i;
          return;
        }
      }
    } while ((*namecolor) == 0);
    break;

  case 3:
    printf("\nWhich of these colors will the background be?:\n");
    printf("   White, Red, Orange, Yellow, Green, Blue, or Violet\n");
    printf(" (W, R, O, Y, G, B, or V)\n");
    do {
      printf(" Choose one: \n");
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
        ch = ' ';
      uppercase(&ch);
      (*backcolor) = 0;
      for (i = 1; i <= 7; i++) {
        if (ch == colors[i - 1].name[0]) {
          (*backcolor) = i;
          return;
        }
      }
    } while ((*backcolor) == 0);
    break;

  case 4:
    printf("\nEnter the X resolution:\n");
    scanf("%ld%*[^\n]", rx);
    getchar();
    printf("Enter the Y resolution:\n");
    scanf("%ld%*[^\n]",ry);
    getchar();
    break;
  }
}  /* getrayparms */

#endif

void plotrparms()
{
  /* set up initial characteristics of plotter or printer */
  Char trash,ch;                       /* colors is declared globally */
  long treecolor, namecolor, backcolor, i, rayresx, rayresy;
  double viewangle;
  long n;

  penchange = no;
  xcorner = 0.0;
  ycorner = 0.0;
  if (dotmatrix && (!previewing))
    strpdiv = 1;
  switch (plotter) {
#ifndef MAC
  case ray:
    penchange = yes;
    xunitspercm = 1.0;
    yunitspercm = 1.0;
    xsize = 10.0;
    ysize = 10.0;
    rayresx = 512;
    rayresy = 512;
    treecolor = 6;
    namecolor = 4;
    backcolor = 1;
    do {
      n=showrayparms(treecolor,namecolor,backcolor,rayresx,rayresy);
      if (n != -1)
        getrayparms(&treecolor,&namecolor,&backcolor,&rayresx,&rayresy,n);
    } while (n != -1);
    xsize = rayresx;
    ysize = rayresy;

    fprintf(plotfile, "report verbose\n");

    fprintf(plotfile, "screen %ld %ld\n", rayresx, rayresy);
    if (ysize >= xsize) {
      viewangle = 2 * atan(ysize / (2 * 1.21 * xsize)) * 180 / pi;
      fprintf(plotfile, "fov 45 %3.1f\n", viewangle);
      fprintf(plotfile, "light 1 point 0 %6.2f %6.2f\n",
              -xsize * 1.8, xsize * 1.5);
      fprintf(plotfile, "eyep %6.2f %6.2f %6.2f\n",
              xsize * 0.5, -xsize * 1.21, ysize * 0.55);
    } else {
      viewangle = 2 * atan(xsize / (2 * 1.21 * ysize)) * 180 / pi;
      fprintf(plotfile, "fov %3.1f 45\n", viewangle);
      fprintf(plotfile, "light 1 point 0 %6.2f %6.2f\n",
              -ysize * 1.8, ysize * 1.5);
      fprintf(plotfile, "eyep %6.2f %6.2f %6.2f\n",
              xsize * 0.5, -ysize * 1.21, ysize * 0.55);
    }

    fprintf(plotfile, "lookp %6.2f 0 %6.2f\n", xsize * 0.5, ysize * 0.5);
    fprintf(plotfile, "/* %.10s */\n", colors[treecolor - 1].name);
    fprintf(plotfile,
            "surface treecolor diffuse %5.2f%5.2f%5.2f specular 1 1 1 specpow 30\n",
            colors[treecolor - 1].red, colors[treecolor - 1].green,
            colors[treecolor - 1].blue);
    fprintf(plotfile, "/* %.10s */\n", colors[namecolor - 1].name);
    fprintf(plotfile,
            "surface namecolor diffuse %5.2f%5.2f%5.2f specular 1 1 1 specpow 30\n",
            colors[namecolor - 1].red, colors[namecolor - 1].green,
            colors[namecolor - 1].blue);
    fprintf(plotfile, "/* %.10s */\n", colors[backcolor - 1].name);
    fprintf(plotfile, "surface backcolor diffuse %5.2f%5.2f%5.2f\n\n",
            colors[backcolor - 1].red, colors[backcolor - 1].green,
            colors[backcolor - 1].blue);
    break;
#endif /* ifndef mac */
  case pict:
    penchange = yes;
    xunitspercm = 28.346456693;
    yunitspercm = 28.346456693;
    /*7.5 x 10 inch default PICT page size*/
    xsize = 19.05;
    ysize = 25.40;
    break;

  case lw:
    penchange = yes;
    xunitspercm = 28.346456693;
    yunitspercm = 28.346456693;
    xsize = 21.59;
    ysize = 27.94;
    break;

  case hp:
    penchange = yes;
    xunitspercm = 400.0;
    yunitspercm = 400.0;
    xsize = 24.0;
    ysize = 18.0;
    break;

  case tek:
    xunitspercm = 50.0;
    yunitspercm = 50.0;
    xsize = 20.46;
    ysize = 15.6;
    break;

  case ibmpc:
#ifdef TURBOC
  GraphDriver = 0;
  detectgraph(&GraphDriver,&GraphMode);
  getmoderange(GraphDriver,&LoMode,&HiMode);
  initgraph(&GraphDriver,&HiMode,"");
  xunitspercm = getmaxx()/25;
  yunitspercm = getmaxy() / 17.5;
  restorecrtmode();
  xsize = 25.0;
  ysize = 17.5;
#endif
#ifdef QUICKC
setupgraphics();

#endif
  break;

  case mac:
    penchange = yes;
    xunitspercm = 33.5958;
    yunitspercm = 33.6624;
    xsize = 15.24;
    ysize = 10.00;
    break;

  case houston:
    penchange = yes;
    xunitspercm = 100.0;
    yunitspercm = 100.0;
    xsize = 24.5;
    ysize = 17.5;
    break;

  case decregis:
    xunitspercm = 30.0;
    yunitspercm = 30.0;
    xsize = 25.0;
    ysize = 15.0;
    break;

  case epson:
    penchange = yes;
    xunitspercm = 47.244;
    yunitspercm = 28.346;
    xsize = 18.70;
    ysize = 22.0;
    strpwide = 960;
    strpdeep = 8;
    strpdiv = 1;
    break;

  case oki:
    penchange = yes;
    xunitspercm = 56.692;
    yunitspercm = 28.346;
    xsize = 19.0;
    ysize = 22.0;
    strpwide = 1100;
    strpdeep = 8;
    strpdiv = 1;
    break;

  case citoh:
    penchange = yes;
    xunitspercm = 28.346;
    yunitspercm = 28.346;
    xsize = 22.3;
    ysize = 26.0;
    strpwide = 640;
    strpdeep = 8;
    strpdiv = 1;
    break;

  case toshiba:
    penchange = yes;
    xunitspercm = 70.866;
    yunitspercm = 70.866;
    xsize = 19.0;
    ysize = 25.0;
    strpwide = 1350;
    strpdeep = 24;
    strpdiv = 4;
    break;

  case pcl:
    penchange = yes;
    xsize = 21.59;
    ysize = 27.94;
    xunitspercm = 118.11023622;   /* 300 DPI = 118.1 DPC                    */
    yunitspercm = 118.11023622;
    strpwide = 2550;   /* 8.5 * 300 DPI                                     */
    strpdeep = 20;     /* height of the strip                               */
    strpdiv = 20;      /* in this case == strpdeep                          */
                       /* this is information for 300 DPI resolution        */
    printf("Please select Laserjet resolution\n\n");
    printf("1:  75 DPI\n2:  150 DPI\n3:  300 DPI\n\n");
    do {
      scanf("%c%*[^\n]", &ch);
      trash=getchar();
      uppercase(&ch);
    } while (ch != '1' && ch != '2' && ch != '3');
    switch (ch) {

    case '1':
      strpwide /= 4;
      xunitspercm /= 4.0;
      yunitspercm /= 4.0;
      hpresolution = 75;
      break;

    case '2':
      strpwide /= 2;
      xunitspercm /= 2.0;
      yunitspercm /= 2.0;
      hpresolution = 150;
      break;

    case '3':
      hpresolution = 300;
      break;
    }
    break;
  case xbm:            /* since it's resolution dependent, make 1x1 pixels  */
    penchange = yes;   /* per square cm for easier math.                    */
    xunitspercm = 1.0;
    yunitspercm = 1.0;
    strpdeep = 10;
    strpdiv = 10;
    printf("Please select the X-bitmap file resolution\n");
    printf("X resolution?\n");
    scanf("%lf%*[^\n]", &xsize);
    getchar();
    printf("Y resolution?\n");
    scanf("%lf%*[^\n]", &ysize);
    getchar();
    strpwide = (long)xsize;
    xsize /= xunitspercm;
    ysize /= yunitspercm;

    break;

  case pcx:
    penchange = yes;
    xsize = 21.16;
    ysize = 15.88;
    strpdeep = 10;
    strpdiv = 10;
    printf("Please select the PCX file resolution\n\n");
    printf("1: EGA 640  X 350\n");
    printf("2: VGA 800  X 600\n");
    printf("3: VGA 1024 X 768\n\n");
    do {
      scanf("%c%*[^\n]", &ch);
      trash=getchar();
      uppercase(&ch);
    } while (ch != '1' && ch != '2' && ch != '3');
    switch (ch) {

    case '1':
      strpwide = 640;
      yunitspercm = 350 / ysize;
      break;

    case '2':
      strpwide = 800;
      yunitspercm = 600 / ysize;
      break;

    case '3':
      strpwide = 1024;
      yunitspercm = 768 / ysize;
      break;
    }
    xunitspercm = strpwide / xsize;
    break;
  case fig:
    penchange = yes;
    xunitspercm = 31.011;
    yunitspercm = 29.78;
    xsize = 25.4;
    ysize = 20.32;
    break;
  case other:
    break;
    /* initial parameter settings for a new plotter go here */
  }
  if (previewing)
    return;
}  /* plotrparms */


void getplotter()
{
  Char ch,trash;

  printf("\nWhich plotter or printer will the tree be drawn on?\n");
  printf("(many other brands or models are compatible with these)\n\n");
  printf("   type:       to choose one compatible with:\n\n");
  printf("        L         Apple Laserwriter (with Postscript)\n");
  printf("        M         MacDraw PICT format\n");
#ifndef MAC
  printf("        R         Rayshade 3D rendering program file\n");
#endif
  printf("        J         Hewlett-Packard Laserjet\n");
  printf("        K         TeKtronix 4010 graphics terminal\n");
  printf("        H         Hewlett-Packard 7470 plotter\n");
#ifdef DOS
  printf("        I         IBM PC graphics screens\n");
#endif
  printf("        D         DEC ReGIS graphics (VT240 terminal)\n");
  printf("        B         Houston Instruments plotter\n");
  printf("        E         Epson MX-80 dot-matrix printer\n");
  printf("        C         Prowriter/Imagewriter dot-matrix printer\n");
  printf("        O         Okidata dot-matrix printer\n");
  printf("        T         Toshiba 24-pin dot-matrix printer\n");
  printf("        P         PC Paintbrush monochrome PCX file format\n");
  printf("        X         X Bitmap format                         \n");
  printf("        F         FIG 2.0 format                          \n");
  printf("        U         other: one you have inserted code for\n");
  do {
    printf(" Choose one: \n");
    scanf("%c%*[^\n]", &ch);
    trash=getchar();
    uppercase(&ch);
  }
#ifdef DOS
while (strchr("LJKHIDBECOTUPXRMF",ch) == NULL);
#else
while (strchr("LJKHDBECOTUPXRMF",ch) == NULL);
#endif
  switch (ch) {

  case 'L':
    plotter = lw;
    break;

  case 'M':
    plotter = pict;
    break;

 case 'R':
    plotter = ray;
    break;

  case 'J':
    plotter = pcl;
    break;

  case 'K':
    plotter = tek;
    break;

  case 'H':
    plotter = hp;
    break;

  case 'I':
    plotter = ibmpc;
    break;

  case 'D':
    plotter = decregis;
    break;

  case 'B':
    plotter = houston;
    break;

  case 'E':
    plotter = epson;
    break;

  case 'C':
    plotter = citoh;
    break;

  case 'O':
    plotter = oki;
    break;

  case 'T':
    plotter = toshiba;
    break;

  case 'P':
    plotter = pcx;
    break;

  case 'X':
    plotter = xbm;
    break;

  case 'F':
    plotter = fig;
    break;

  case 'U':
    plotter = other;
    break;
  }
  dotmatrix = (plotter == epson || plotter == oki || plotter == citoh ||
               plotter == toshiba || plotter == pcx || plotter == pcl ||
               plotter == xbm);
  printf("\nWhich type of screen will it be previewed on?\n\n");
  printf("   type:       to choose one compatible with:\n\n");
  printf("        N         will not be previewed\n");
#ifdef DOS
  printf("        I         IBM PC graphics screens\n");
#else
# ifdef MAC
  printf("        M         Macintosh screens\n");
# else
  printf("        K         TeKtronix 4010 graphics terminal\n");
  printf("        D         DEC ReGIS graphics (VT240 terminal)\n");
  printf("        U         other: one you have inserted code for\n");
# endif
#endif
  do {
    printf(" Choose one: \n");
    scanf("%c%*[^\n]", &ch);
    trash=getchar();
    uppercase(&ch);
  }
#ifdef DOS
  while (strchr("NIKDU",ch) == NULL);
#else
# ifdef MAC
  while (strchr("NMKDU",ch) == NULL);
#  else
  while (strchr("NKDU",ch) == NULL);
#  endif
#endif
  preview = true;
  switch (ch) {

  case 'N':
    preview = false;
    break;

  case 'I':
    previewer = ibmpc;
    break;

  case 'M':
    previewer = mac;
    break;

  case 'K':
    previewer = tek;
    break;

  case 'D':
    previewer = decregis;
    break;

  case 'U':
    previewer = other;
    break;
  }
  printf("\n\n\n");
}  /* getplotter */


void changepen(pen)
pentype pen;
{
  Char picthi, pictlo;
  long  pictint;

 switch (pen) {

  case treepen:
    linewidth = treeline;
    if (plotter == hp)
      fprintf(plotfile, "SP1;\n");
    if (plotter == lw) {
      fprintf(plotfile, "stroke %8.2f setlinewidth \n", treeline);
      fprintf(plotfile, " 1 setlinecap 1 setlinejoin \n");
    }
    break;

  case labelpen:
    linewidth = labelline;
    if (plotter == hp)
      fprintf(plotfile, "SP2;\n");
    if (plotter == lw) {
      fprintf(plotfile, " stroke%8.2f setlinewidth \n", labelline);
      fprintf(plotfile, "1 setlinecap 1 setlinejoin \n");
    }
    break;
  }
#ifdef MAC
if (plotter == mac){
      pictint = ( long)(linewidth + 0.5);
      if (pictint ==0)
           pictint = 1;
      PenSize((int)pictint,(int)pictint);}
#endif

  if (plotter != pict)
    return;
  pictint = ( long)(linewidth + 0.5);
  if (pictint == 0)
    pictint = 1;
  picthi = (Char)(pictint / 256);
  pictlo = (Char)(pictint & 255);
  fprintf(plotfile, "\007%c%c%c%c", picthi, pictlo, picthi, pictlo);
}  /* changepen */


double lengthtext(pstring, nchars,font)
Char *pstring;
 long nchars;
fonttype font;
{  /* lengthext */
  long i, j, code;
  static double sumlength;
  sumlength = 0.0;
  for (i = 0; i < nchars; i++) {
    code = pstring[i];
    j = 1;
    while (font[j] != code && font[j - 1] != 0)
      j = font[j - 1];
    if (font[j] == code)
      sumlength += font[j + 2];
  }
  return sumlength;
}  /* lengthtext */


void plotchar(place,text)
 long *place;
struct LOC_plottext *text;         /* variables passed from plottext */
{
  text->heightfont = text->font[*place + 1];
  text->yfactor = text->height / text->heightfont;
  text->xfactor = text->yfactor;
  *place += 3;
  do {
    (*place)++;
    text->coord = text->font[*place - 1];
    if (text->coord > 0)
      text->penstatus = pendown;
    else
      text->penstatus = penup;
    text->coord = abs(text->coord);
    text->coord %= 10000;
    text->xfont = (text->coord / 100 - xstart) * text->xfactor;
    text->yfont = (text->coord % 100 - ystart) * text->yfactor;
    text->xplot = text->xx + (text->xfont * text->cosslope +
			      text->yfont * text->sinslope) * text->compress;
    text->yplot = text->yy - text->xfont * text->sinslope +
      text->yfont * text->cosslope;
    plot(text->penstatus, text->xplot, text->yplot);
  } while (abs(text->font[*place - 1]) < 10000);
  text->xx = text->xplot;
  text->yy = text->yplot;
}  /* plotchar */

swap(one,two)
char **one,**two;
{
char *tmp = (*one);
(*one)= (*two);
(*two) = tmp;
return;
}

void drawit(fontname,xoffset,yoffset,numlines,root)
char *fontname;
double *xoffset,*yoffset;
 long  numlines;
node *root;
{
  long i, j, line;
  long deep,iterations;
  (*xoffset) = 0.0;
  (*yoffset) = 0.0;
  if (dotmatrix){
    strptop    = ( long)(ysize * yunitspercm);
    strpbottom = numlines*strpdeep + 1;
  }
  else {
    plottree(root,root);
    plotlabels(fontname);
  }
  if (dotmatrix){
    striprint(( long)((ysize * yunitspercm)- (numlines * strpdeep)),
	      ( long)((ysize * yunitspercm)- (numlines * strpdeep)));
    strptop = numlines * strpdeep;
    strpbottom = strptop - strpdeep + 1;
    printf(" writing%3ld lines ...\n", numlines);
    printf("  Line     Output file size\n");
    printf("  ----     ------ ---- ----\n");
    for (line = 1; line <= numlines ; line++) {
      for (i = 0; i <= strpdeep ; i++){
        for (j=0; j<=(strpwide/8);++j)
	  stripe[i][j] = 0;}
      empty = true;
      xnow = strpwide / 2.0;
      ynow = 0.0;
      plottree(root, root);
      plotlabels(fontname);
      strptop = strpbottom - 1;
      strpbottom -= strpdeep;
      deep=20;
      if (strpdeep > deep){              /* large stripe, do in 20-line     */
        for (i=0;i<strpdeep;++i){
	  swap(&stripe[i%deep],&stripe[i]);
	  if ((i%deep) == (deep -1)){
	    striprint(deep,deep);}
        }
        striprint(strpdeep%deep,strpdeep%deep);
      }
      else{                          /* small stripe, do it all now.     */
        striprint(strpdiv,strpdeep);
        if (line % 5 == 0)
          printf("%5ld%16ld\n", line, filesize);
      }
    }
  }
}  /* drawit */


void plottext(pstring, nchars, height_, cmpress2, x, y, slope, font_,fontname)
Char *pstring;
 long nchars;
double height_, cmpress2, x, y, slope;
short *font_;
char *fontname;
{
  struct LOC_plottext text;
  long i, j, code;
  double pointsize;
  text.heightfont = font_[2];

  pointsize = (1.2*(height_/text.heightfont)*cmpress2 / 2.54) * 72.0;
  text.height = height_;
  text.compress = cmpress2;
  text.font = font_;
  text.xx = x;
  text.yy = y;
  text.sinslope = sin(pi * slope / 180.0);
  text.cosslope = cos(pi * slope / 180.0);
  if (previewing || (strcmp(fontname,"Hershey") == 0)){
    for (i = 0; i < nchars; i++) {
      code = pstring[i];
      j = 1;
      while (text.font[j] != code && text.font[j - 1] != 0)
	j = text.font[j - 1];
      plotchar(&j,  &text);
    }
  }

  else if (plotter == lw)  {
 /* print alternate font.  Right now, only postscript. */
    fprintf(plotfile,"gsave\n");
    fprintf(plotfile,"/%s findfont %f scalefont setfont\n",fontname,
	    pointsize);
    fprintf(plotfile,"%f %f translate %f rotate\n",x,y,-slope);
    fprintf(plotfile,"0 0 moveto\n");
    fprintf(plotfile,"(%s) show\n",pstring);
    fprintf(plotfile,"grestore\n");
  }
}  /* plottext */


void makebox(fn,xo,yo,scale,ntips)
char *fn;                          /* fontname                       */
double *xo,*yo;                    /* x and y offsets                 */
double *scale;
 long ntips;
{
  /* draw the box on screen which represents plotting area.        */
  char ch;

  printf("\nWe now will preview the tree.  The box that will be\n");
  printf("plotted on the screen represents the boundary of the\n");
  printf("final plotting surface.  To see the preview, press on\n");
  printf("the ENTER or RETURN key (you may need to do it twice).\n");
  printf("When finished viewing it, press on that key again.\n");
  oldpenchange   = penchange;
  oldxsize       = xsize;
  oldysize       = ysize;
  oldxunitspercm = xunitspercm;
  oldyunitspercm = yunitspercm;
  oldxcorner     = xcorner;
  oldycorner     = ycorner;
  oldplotter     = plotter;
  plotter        = previewer;
  scanf("%c%*[^\n]", &ch);
  (void)getchar();
  if (ch == '\n')
    ch = ' ';
  plotrparms();
  xcorner += 0.05 * xsize;
  ycorner += 0.05 * ysize;
  xsize *= 0.9;
  ysize *= 0.9;
  (*scale) = ysize / oldysize;
  if (xsize / oldxsize < (*scale))
    (*scale) = xsize / oldxsize;
  (*xo) = (xcorner + (xsize - oldxsize * (*scale)) / 2.0) / (*scale);
  (*yo) = (ycorner   + (ysize - oldysize * (*scale)) / 2.0) / (*scale);
  xscale = (*scale) * xunitspercm;
  yscale = (*scale) * yunitspercm;
  initplotter(ntips,fn);
  plot(penup, xscale * (*xo), yscale * (*yo));
  plot(pendown, xscale * (*xo), yscale * ((*yo) + oldysize));
  plot(pendown, xscale * ((*xo) + oldxsize), yscale * ((*yo) + oldysize));
  plot(pendown, xscale * ((*xo) + oldxsize), yscale * (*yo));
  plot(pendown, xscale * (*xo), yscale * (*yo));
}  /* makebox */


boolean plotpreview(fn,xo,yo,scale,nt,root)
char *fn;        /* font name                                              */
double *xo,*yo;
double *scale;
 long nt;        /* ntips                                                  */
node *root;
{
  boolean canbeplotted;
  Char ch;

  previewing = true;
  makebox(fn,xo,yo,scale,nt);
  plottree(root, root);
  plotlabels(fn);
  finishplotter();
  penchange = oldpenchange;
  xsize = oldxsize;
  ysize = oldysize;
  xunitspercm = oldxunitspercm;
  yunitspercm = oldyunitspercm;
  xscale = xunitspercm;
  yscale = yunitspercm;
  plotter = oldplotter;
  xcorner = oldxcorner;
  ycorner = oldycorner;
  printf(" Is the tree ready to be plotted? (Answer Y or N)\n");
  do {
    printf("Type Y or N:\n");
    scanf("%c%*[^\n]", &ch);
    (void)getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
  } while (ch != 'Y' && ch != 'N');
  canbeplotted = (ch == 'Y');
  return canbeplotted;
}  /* plotpreview */


long allocstripe(stripe,x,y)
striptype stripe;
long x,y;
{
  long i,stripedepth;
  for (i=0 ; i<=y;++i){
    stripe[i] = (MALLOCRETURN *)malloc((x+1)*sizeof(Char));
    if (!stripe[i])
      break;
   }
return i-1;
}

