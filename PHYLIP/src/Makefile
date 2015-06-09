LIBS    = -lm
BINDIR  = ../bin
CFLAGS  = -O4
CC        = gcc $(CFLAGS)
PROGS   =          clique consense contml contrast dnacomp dnadist \
                   dnainvar dnaml dnamlk dnamove dnapars dnapenny \
                   dolmove dollop dolpenny drawgram drawtree factor \
                   fitch gendist kitsch mix move neighbor penny \
                   protdist protpars restml retree seqboot

all:        introduce $(PROGS)
	@echo done.

install: all
	@echo "Installing PHYLIP v3.5c binaries in $(BINDIR)"
	@cp $(PROGS) $(BINDIR)

clique:       clique.c phylip.h
	$(CC) clique.c $(LIBS) -o clique

consense:     consense.c phylip.h
	$(CC) consense.c $(LIBS) -o consense

contml:       contml.c phylip.h
	$(CC) contml.c $(LIBS) -o contml

contrast:     contrast.c phylip.h
	$(CC) contrast.c $(LIBS) -o contrast

dnacomp:      dnacomp.c phylip.h
	$(CC) dnacomp.c $(LIBS) -o dnacomp

dnadist:      dnadist.c phylip.h
	$(CC) dnadist.c $(LIBS) -o dnadist

dnainvar:     dnainvar.c phylip.h
	$(CC) dnainvar.c $(LIBS) -o dnainvar

dnaml:        dnaml.o dnaml2.o
	$(CC) dnaml.o dnaml2.o $(LIBS) -o dnaml

dnaml.o:        dnaml.c phylip.h
	$(CC) -c dnaml.c

dnaml2.o:        dnaml2.c phylip.h
	$(CC) -c dnaml2.c

dnamlk:       dnamlk.o dnamlk2.o
	$(CC) dnamlk.o dnamlk2.o $(LIBS) -o dnamlk

dnamlk.o: dnamlk.c phylip.h
	$(CC) -c dnamlk.c

dnamlk2.o: dnamlk2.c phylip.h
	$(CC) -c dnamlk2.c

dnamove:      dnamove.c phylip.h
	$(CC) dnamove.c $(LIBS) -o dnamove

dnapenny: dnapenny.c phylip.h
	$(CC) dnapenny.c $(LIBS) -o dnapenny

dnapars:      dnapars.c phylip.h
	$(CC) dnapars.c $(LIBS) -o dnapars

dolmove:      dolmove.c phylip.h
	$(CC) dolmove.c $(LIBS) -o dolmove

dollop:       dollop.c phylip.h
	$(CC) dollop.c $(LIBS) -o dollop

dolpenny:     dolpenny.c phylip.h
	$(CC) dolpenny.c $(LIBS) -o dolpenny

drawgraphics.o:   drawgraphics.c phylip.h drawgraphics.h
	$(CC) -c drawgraphics.c

drawtree.o:     drawtree.c
	$(CC) -c drawtree.c

drawgram.o:     drawgram.c
	$(CC) -c drawgram.c

drawgram:     drawgram.o phylip.h drawgraphics.o
	$(CC) drawgraphics.o drawgram.o $(LIBS) -o drawgram

drawtree:     drawtree.o phylip.h drawgraphics.o
	$(CC) drawgraphics.o drawtree.o $(LIBS) -o drawtree

factor:       factor.c phylip.h
	$(CC) factor.c $(LIBS) -o factor

fitch:        fitch.c phylip.h
	$(CC) fitch.c $(LIBS) -o fitch

gendist:      gendist.c phylip.h
	$(CC) gendist.c $(LIBS) -o gendist

kitsch:       kitsch.c phylip.h
	$(CC) kitsch.c $(LIBS) -o kitsch

mix:        mix.o mix2.o
	$(CC) mix.o mix2.o $(LIBS) -o mix

mix.o:        mix.c phylip.h
	$(CC) -c mix.c

mix2.o:        mix2.c phylip.h
	$(CC) -c mix2.c

move:        move.c phylip.h
	$(CC) move.c $(LIBS) -o move

neighbor:     neighbor.c phylip.h
	$(CC) neighbor.c $(LIBS) -o neighbor

penny:  penny.c phylip.h
	$(CC) penny.c $(LIBS) -o penny

protdist:      protdist.c phylip.h
	$(CC) protdist.c $(LIBS) -o protdist

protpars: protpars.c phylip.h
	$(CC) protpars.c $(LIBS) -o protpars

restml:       restml.o restml2.o phylip.h
	$(CC) restml.o restml2.o $(LIBS) -o restml

restml.o: restml.c
	$(CC) -c restml.c

restml2.o: restml2.c
	$(CC) -c restml2.c

retree:       retree.c phylip.h
	$(CC) retree.c $(LIBS) -o retree

seqboot:      seqboot.c phylip.h
	$(CC) seqboot.c $(LIBS) -o seqboot

introduce:
	@echo "Building PHYLIP version 3.5c "

clean:
	@rm $(PROGS) *.o
