/* stderr */

/* Interface
   By Sean T. Lamont
   For use with the Macntosh version of the Phylogeny Inference Package,
   version 3.5c. (c) Copyright 1992 by Joseph Felsenstein.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed.

   This file defines a 2-window environment which replicates many
   of the standard c I/O functions (puts, printf,gets,scanf)
   in a maclike environment, with a scrollback buffer.  This
   is necessary because of the very weak implementation of IO under
   straight C, and the difficulty of using separate windows. It also
   adds the ability to "quit" out of an application while it is running,
   and some very basic point-and-click interface type of stuff on the
   screen.

   Functions you need to know how to use:
   macsetup(char *name, char *name):  initializes the interface, brings up a
                                  windows of the name of the argument, for I/O.
   macprintf(char *string,arg,arg):  like printf but into the window
   macgets(char *string);            like puts   but into the window
   macscanf(char *string, arg,arg)   like scanf  but from the window
   macgets (char *string,arg,arg)    like gets but from the window.

   It is recommended that you include "interface.h" file, which sets
   up these #define's for you (IE all calls to scanf will call
   macscanf, etc.)

  textmode(); makes the current drawing window the text window
   gfxmode(); makes the current drawing window the gfx window,
              unhides the graphics window.

  eventloop():  process mouse events, menus, etc., and wait for
                "go away" event.  You are implicitly doing this\
                 when you run macgets/macscanf.  If you want to query
                 the event handler during non-IO periods of time this
                 would be what to call (or alternatively, printf("");
 */


#include <stdarg.h>
#include <stdio.h>
#include "interface.h"

#define SCROLLLINES 250

#define SCROLLCOLUMNS 80
#define SCROLLMEM    22000
#define LINEDEPTH  (int)12
#define PAGELINES (int)(260 / 12)
#define TEXT 0
#define GFX 1

/* Global variables, most for use with input/output */
int disablescrollback = 0;
int mode = TEXT;
WindowPtr gfxWindow;
WindowPtr textWindow;
MenuHandle appleMenu,fileMenu;
char *lines[SCROLLLINES]; /* the scrollback buffer */
int cursorx=0; /* the current x position that the text cursor is at */
int cursory=0; /* the current y position that the text cursor is at */
int numlines_=0;/* the next place we put a new line.                 */
int memptr=0;  /* the next place in our space we get memory from.   */
int linectr=0;
char memory[SCROLLMEM];
char line[256];/* used to accumulate output before newlines.        */
Rect textBounds = {40,5,300,500}; /* position of window 1          */
Rect gfxBounds  = {0,0,350,500}; /* position of window 2          */
Rect textrect   = {0,0,260,479};   /* region to clear within window */
Rect barBounds  = {0,480,260,495}; /* position of the scroll bar    */
RgnHandle rgn;                     /* used by scrollrect            */
ControlHandle bar;                 /* this is the scrollbar control */
int lasttop=0;                     /* used by update / scrollrect   */
char inputs[256];                  /* used by the input routines    */
int  inputcount = 0;               /* offset into the string.       */
int  collect    = false;           /* boolean: collect chars?       */

void macsetup(tname,gname)
char *tname,*gname;
{
/*
static char buf1[128];
static char buf2[128];
*/
Str255 buf1,buf2; /*PB feb95*/
Str255 title="";  /*PB feb95*/

strcpy(buf1+1,tname);
strcpy(buf2+1,gname);
buf1[0]=strlen(tname);
buf2[0]=strlen(gname);
MaxApplZone();
#ifdef METRO
InitGraf(&qd.thePort);  /*PB feb95, for metrowerks codewarrior*/
#else
InitGraf(&thePort)
#endif
InitFonts();
FlushEvents(everyEvent,0);
InitWindows();
InitMenus();
TEInit();
InitDialogs(0L);
InitCursor();
textWindow = NewWindow(0L,&textBounds,buf1,true,documentProc,
            (WindowPtr) -1L, true,0);
gfxWindow = NewWindow(0L,&gfxBounds,buf2,false,noGrowDocProc,
            (WindowPtr) -1L, true,0);
rgn=NewRgn();
SetPort(textWindow);
bar=NewControl(textWindow,&barBounds,title,true,0,0,
	           0,scrollBarProc,0);
InsertMenu(appleMenu=NewMenu(1,"\p\024"),0);     /* add apple menu  */
InsertMenu(fileMenu=NewMenu(2,"\pFile"),0);      /* add file menu   */
TextFont(courier);
TextSize(10);
DrawMenuBar();
AppendMenu(fileMenu,"\pQuit/Q");
AddResMenu(appleMenu, 'DRVR');
macprintf("\n");
}


void queryevent()
{
int status;
status=handleevent();
if (status <= 0)
     process_window_closure(status);
}

void eventloop()
{
int status;

while (1){
	status=handleevent();
if (status <= 0)
		process_window_closure(status);  }
}


process_window_closure(status)
int status;
{
if (status == -1){
 		CloseWindow(gfxWindow);   /* "Close main window", so run all the */
		CloseWindow(textWindow);  /* cleanup stuff.                      */
		DisposeRgn(rgn);
#undef exit
		exit(0);
#define exit(status) eventloop()
}
else if (status == 0)
   HideWindow(gfxWindow);}



int handleevent()
{
pascal void scroll();
char pstring[256];
OSErr res;
EventRecord ev;
WindowPtr win;
ControlHandle ctrl;
short menuid,menuitem;
Point MouseLoc;
long menu,i;
int cx,cy;
int PathRefNum;
long  count=80;
SFReply fileinfo;
char c,pasteword[64];
Str255 name;
GrafPtr savePort;
FILE *fp;
Rect drect;
Rect      rect  = {0,0,1000,1000};   /* the limit for dragging windows */
int ok=GetNextEvent(everyEvent,&ev);
int where=FindWindow(ev.where,&win);
if (ev.what == keyDown && collect && mode != GFX){
     SetCtlValue(bar,GetCtlMax(bar));
     redraw(numlines_ -  (GetCtlMax(bar) - GetCtlValue(bar)) - PAGELINES);
     }
if ((ev.what == keyDown) &&
    (ev.modifiers &  cmdKey) &&
    (toupper((char)(ev.message & charCodeMask)) == 'Q')  )
         return -1;
if (( ev.what == keyDown ) && collect && mode == GFX){
     textmode();
     process_char((char)(ev.message & charCodeMask));}
if (ev.what == mouseDown && !disablescrollback && where ==inContent & mode == TEXT){
        SelectWindow(win);
        GlobalToLocal(&ev.where);
         switch (FindControl(ev.where,win,&ctrl)){
       	case inThumb:
    	TrackControl(ctrl,ev.where,nil);
       	 break;
        case inUpButton:
         /*TrackControl(ctrl,ev.where, scroll); PB feb95*/
         TrackControl(ctrl,ev.where,(ControlActionUPP) scroll);
         break;
        case inPageUp:
           SetCtlValue(bar,GetCtlValue(bar) - 10);
         break;
        case inDownButton:
          /*TrackControl(ctrl,ev.where, scroll); PB feb95*/
         TrackControl(ctrl,ev.where,(ControlActionUPP) scroll);
         break;
        case inPageDown:
            SetCtlValue(bar,GetCtlValue(bar) + 10);
         break;
        default:
         if (collect && mode == TEXT){
            GetMouse(&MouseLoc);
            cy = (int)((double)MouseLoc.v / (double)(LINEDEPTH)) +
                 (numlines_ - (GetCtlMax(bar) - GetCtlValue(bar))
                   - PAGELINES + 1);
            for (i=0;i<strlen(lines[cy%SCROLLLINES]);++i)
                 if ((lines[cy%SCROLLLINES])[i] == '(' ||
                     (lines[cy%SCROLLLINES])[i] == '('   )
                      (lines[cy%SCROLLLINES])[i]=' ';
            sscanf(lines[cy%SCROLLLINES]," %[0-9a-zA-Z]",pasteword);
              SetCtlValue(bar,GetCtlMax(bar));
              redraw(numlines_ -  (GetCtlMax(bar) - GetCtlValue(bar)) - PAGELINES);
            for (i=0;i<strlen(pasteword);++i)
                 process_char(pasteword[i]);
            process_char(0x0d);
            }

            break;
           }
	redraw(numlines_ -  (GetCtlMax(bar) - GetCtlValue(bar)) - PAGELINES);
    if (GetCtlMax(bar) == GetCtlValue(bar)){
    	macflush();
   		numlines_--;}
    }

else if (ev.what == activateEvt)
    InvalRect(&win->portRect);
  else if (ev.what == updateEvt && mode == TEXT){
 	BeginUpdate(textWindow);
	lasttop=100000;
    redraw(numlines_ -  (GetCtlMax(bar) - GetCtlValue(bar)) - PAGELINES);
	DrawControls(textWindow);
 	EndUpdate(textWindow);
	}

else if (ev.what == mouseDown && where == inSysWindow)
    SystemClick(&ev,win);
else if (ev.what == mouseDown && where == inDrag) {
	DragWindow(win,ev.where,&rect);}
	
else if (ev.what == mouseDown  && where == inGoAway)
     return (win == gfxWindow ? 1: -1);

if (ev.what == mouseDown && where == inMenuBar){
	menu=MenuSelect(ev.where);
	menuitem = LoWord(menu);
	menuid   = HiWord(menu);
if (menuid == 2 && menuitem == 1)
		return -1;
if (menuid == 1){
	GetPort(&savePort);
	GetItem(appleMenu, menuitem, name);
    OpenDeskAcc(name);
	SetPort(savePort);}

	}
else if (collect && mode == TEXT && (ev.what == keyDown))
	   process_char((char)(ev.message & charCodeMask));
	
return 1;
}


void macgets(s)
char *s;
{
int status;
collect=true;
if (mode == GFX){
  do {status=handleevent();} /* loop until this is false, or hit cr  */
     while (collect && status);

if (status<= 0) process_window_closure(status);

 }

else {
    macflush();               /* flush any waiting output (prompt?)   */
	numlines_--;
    inputcount=0;
    collect=true;              /* tell the event loop to colect chars  */
    do {status=handleevent();} /* loop until this is false, or hit cr  */
       while (collect&&(status>0));
    if (status<= 0) process_window_closure(status);
    inputs[inputcount]=0;
    strcpy(s,inputs);
    macprintf("\n");
    }
}

int macscanf(char *s,...)
{
int i;
char buf[256];
void *p[10];
va_list args;
gets(buf);

va_start(args,s);
for (i=0;i<10;++i)
	p[i]=va_arg(args,void *);
va_end(args);
return sscanf(buf,s,p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9]);
}

void macputs(s)
char *s;
{
macprintf("%s\n",s);
}

macputchar(c)
char c;
{
macprintf("%c",c);
}

void macprintf(char *s,...)
{
char buf[256];

va_list args;
int i;
if (mode == GFX)
     return;
if (strcmp(s,"\033[2J\033[H") ==0)
     return;
disablescrollback=1;
va_start(args,s);
vsprintf(buf,s,args);
va_end(args);
queryevent();
SetPort(textWindow);
MoveTo(cursorx,cursory);
for (i=0;i<strlen(buf);++i){
	if (isprint(buf[i]))
		line[linectr++]=buf[i];	
	else if (buf[i] == '\n'){
		macflush();
		linectr=0;
		macnewline();
		}
	}
	disablescrollback=0;
}

macflush()
{
line[linectr]=0;
if ((memptr+strlen(line) +1) > SCROLLMEM)
     memptr=0;

lines[(numlines_)%SCROLLLINES]=&memory[(numlines_%SCROLLLINES)*SCROLLCOLUMNS];
strcpy(lines[(numlines_++)%SCROLLLINES],line);
MoveTo(0,cursory);
putstring(line);
}

macclear()
{
cursorx=0;
cursory=0;
EraseRect(&textrect);
}

macnewline()
{
int i;
cursorx=0;
cursory+=LINEDEPTH;
if (cursory >  260  ){
	cursory-=LINEDEPTH;}
MoveTo(cursorx,cursory);
 SetCtlMax(bar,((numlines_ > SCROLLLINES ) ? SCROLLLINES  :
               (numlines_ )) - PAGELINES);
SetCtlValue(bar,GetCtlMax(bar));
 redraw(numlines_ -  (GetCtlMax(bar) - GetCtlValue(bar)) - PAGELINES);

}


redraw(pos)
int pos;
{
int i,j,y=0;
int lastbot = lasttop+pos;
int delta   = pos-lasttop;
int delta2  = lasttop-pos;
int lbound  = pos;
int ubound  = pos+PAGELINES;
int tmpx    = cursorx;
int tmpy    = cursory;
if (numlines_ < PAGELINES)
	delta=delta2=10000;
textmode();
pos = (pos < 0 ? 0 : pos);
if (delta >= 0 && delta < SCROLLLINES ){
	ScrollRect(&textrect,0,-delta*LINEDEPTH,rgn);
	lbound = pos+PAGELINES-delta;
	}
else if (delta2 > 0 && delta2 < SCROLLLINES ){
	ScrollRect(&textrect,0,delta2*LINEDEPTH,rgn);
	ubound=delta2+pos;

}

else if (numlines_ >= PAGELINES)
     EraseRect(&textrect);
if (numlines_ == ubound)
     ubound--;
for (i=pos;i<pos-1+((PAGELINES < numlines_) ? numlines_ : PAGELINES);++i){
	MoveTo(0,y);
	if (i>= lbound && i<= ubound)
		putstring(lines[i%SCROLLLINES]);
	y+=LINEDEPTH;
	
 }
 MoveTo(0,y);
 y+=LINEDEPTH;
lasttop = pos;
MoveTo(tmpx,tmpy);
}

putstring(string)
char *string;
{
/*char buf[256]; PB feb 95*/
Str255 buf;
strcpy(buf+1,string);
buf[0]=strlen(string);
DrawString(buf);
cursorx+=StringWidth(buf);
}

textmode()
{
SetPort(textWindow);
SelectWindow(textWindow);
ShowWindow(textWindow);
HideWindow(gfxWindow);
mode = TEXT;
}

gfxmode()
{
int status,ok;
EventRecord ev;
char c;

SetPort(gfxWindow);
ShowWindow(gfxWindow);
HideWindow(textWindow);
SelectWindow(gfxWindow);
mode = GFX;
}

pascal void scroll(c,p)
ControlHandle c;
int p;
{
int direction=((p == inDownButton) ? 1 : (p ==  0)  ? 0 : -1);
SetCtlValue(bar,GetCtlValue(bar) + direction);
/* redraw(numlines_ -  (GetCtlMax(bar) - GetCtlValue(bar)) - PAGELINES);
*/
}

 process_char(c)
char c;
{
Rect drect;

if (isprint(c)){
		inputs[inputcount++]=c;
		line[linectr++]=c;
		if (GetCtlValue(bar) == GetCtlMax(bar) && mode==TEXT) {
		macflush();
		numlines_--;}
		}
	else{
	switch (c){
		case 0x03:     /* if it's the enter key */
		case 0x0d:     /* or the return key     */
			collect=false; /* stop collecting chars */
			break;
		case 0x08: /* delete */
		case 0x1c: /* or back space */
			if (inputcount > 0 && mode == TEXT) {
				cursorx-=CharWidth(inputs[--inputcount]);
				MoveTo(cursorx,cursory);
				inputs[inputcount]=0;
	     		drect.top=cursory-LINEDEPTH;
	     		drect.left=0;
	     		drect.bottom=cursory+3;
	     		drect.right=344;
	     		EraseRect(&drect);
	     		line[--linectr]=0;
				if (GetCtlValue(bar) == GetCtlMax(bar)){
				macflush();
				numlines_--;}
				}
			break;
		default:
			break;}
	  }
}

void fixmacfile(filename)
char *filename;
{
OSErr retcode;
FInfo  fndrinfo;
char filename1[100];
char filename2[100];
strcpy(filename1,filename);
strcpy(filename2,filename);
retcode=GetFInfo(CtoPstr(filename1),0,&fndrinfo);
fndrinfo.fdType='TEXT';
fndrinfo.fdCreator='MSWD';
retcode=SetFInfo(CtoPstr(filename2),0,&fndrinfo);
}

