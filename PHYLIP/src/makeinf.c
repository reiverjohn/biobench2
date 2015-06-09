# include <stdio.h>

char aliname[80], nucname[80], destname[80], taxname[11][255];
FILE *sourceali, *sourcenuc, *destali;
        /* files: their names and pointers
        */
int totnum, seqnum, blocknum, amslr, dpvou,
             nucoraa, taxnum, globcount, position;
        /* totnum is the total number of sequences in the alignment
           seqnum is the number of sequences that will be used
           counters that get incremented or decremented:
            amslr, bntms, count, dpvou, etc
           nucoraa: 1 for nuc alignments, 3 for amino acid alignments
           taxnum is the ID-number of the current taxon
           globcount is the counter for the main loop
                position is the number for the codon position to be analyzed
        */
int oldcheck, newcheck;   /* for checking the number of characters */
int nucmito, yorm, namebool;
        /* for nuclear genetic code: nucmito == 0
           for mitochondrial code : nucmito == 1
           if yorm == 0: no conversion of leu and arg first positions;
           if yorm == 1: conversion occurs
        */
int putnucbool;
        /* for exclusion of regions where homology is uncertain */
long begin, remember;
        /* current pointer position within sourceali and destali */
double  actualcount, count[9];


main()
{
long temp;                  /* holds a temporary position of pointer
                               within destali */

/* INITIALIZING GLOBAL VARIABLES */

position = 5;			/* initialize to analyze all positions */
putnucbool = 0;			/* initialize exclude-boolean */
actualcount = 0;		/* initialize base counter */
yorm = 2;			/* initialize 'Y' or 'M' option */
nucmito = 2;			/* initialize type of genetic code option */
for (amslr = 0; amslr < 9; amslr++) {
	count[amslr] = 0;
}				/* initialize base counters */

/* FUNCTION CALLS START HERE */

parstuff();			/* gets info from user */
findalignment();		/* finds the string ALIGNMENT */
findbeginning();		/* finds the start of the sequence alignment */
blocknum = countblocks();	/* counts the number of blocks in the
				 * alignment */
goback();			/* pointer back to last \n before alignment */

/* MAIN LOOP. ONE EXECUTION OF IT PROCESSES ONE SEQUENCE */

for (globcount = 0; globcount < seqnum; globcount++) {
	if (globcount > 0) {	/* if it's not the first sequence */
		jump(1);	/* puts pointer at beginning of sequence line */
		begin = ftell(sourceali);
	}
	findname();		/* finds the number of this sequence in
				 * sourceali */
	goback();		/* goes back to beginning of line */
	findnucseq();		/* finds the corresponding sequence in
				 * sourcenuc */

	for (dpvou = 0; dpvou < blocknum; dpvou++)
		/*
		 * loops through the alignment as many times as there are
		 * blocks
		 */
	{
		putnucs();	/* grabs nucleotides, puts them into destali */
		if (dpvou != (blocknum - 1))	/* if it's not the last
						 * block.. */
			jump(seqnum - 1);	/* ..jump */
		else if (dpvou == (blocknum - 1)) {
			tailcheck(globcount);	/* checks if there are bases
					         * remaining beyond the end
					         * of the sequence */
			goback();	/* go back to beginning of previous
					 * sequence in alignment */
		}
	}			/* end of single sequence processing */

	if (globcount == 0) {
		oldcheck = newcheck;	/* length of first sequence is
					 * standard */
		temp = ftell(destali);
		fseek(destali, remember, 0);
		fprintf(destali, "%10d", oldcheck);
		fseek(destali, temp, 0);
	} else if (oldcheck != newcheck) {
		printf("CRASH! Lengths of sequences do not agree. Check %s\n", aliname);
		exit(0);
	}
	if ((newcheck % 60) != 0)
		fputc('\n', destali);
	newcheck = 0;		/* resets counter */
}

/* END OF MAIN LOOP */

/* CALCULATE BASE COMPOSITION IF FIRST POSITIONS OF CODONS ARE INCLUDED */
if ((nucoraa == 3) && (position != 2) &&
    (position != 3) && (yorm == 1))
	summer();

fclose ( destali );
fclose ( sourcenuc );
fclose ( sourceali );

return 0;
}        /* end of main */


parstuff()
{				/* gets info from user */
	int             bool, dumnum;
	char            dummy, nucchar;

	do {
		printf("Alignment file to be read:      ");
		scanf("%s", aliname);
		if (NULL == (sourceali = fopen(aliname, "r"))) {
			bool = 0;
			printf("Cannot open %s. Check path or spelling.\n", aliname);
		} else
			bool = 1;
	}
	while (bool == 0);

	do {
		printf("Nucleotide file to be read:     ");
		scanf("%s", nucname);
		if (NULL == (sourcenuc = fopen(nucname, "r"))) {
			bool = 0;
			printf("Cannot open %s. Check path or spelling.\n", nucname);
		} else
			bool = 1;
	}
	while (bool == 0);

	do {
		printf("Destination file to be written: ");
		scanf("%s", destname);
		if (NULL == (destali = fopen(destname, "w"))) {
			bool = 0;
			printf("Cannot open %s. Check path or spelling.\n", destname);
		} else
			bool = 1;
	}
	while (bool == 0);

	do {
		printf("\nTotal number of sequences in alignment: ");
		scanf("%d", &totnum);
		if (totnum < 2) {
			printf("%d; CRASH! That won't work.\n", &totnum);
			exit(0);
		} else
			bool = 1;
	}
	while (bool == 0);
	scanf("%c", &dummy);

	do {
		printf("Number of sequences to be used:         ");
		scanf("%d", &seqnum);
		if (seqnum < 2) {
			printf("%d; CRASH! That won't work.\n", &seqnum);
			exit(0);
		} else
			bool = 1;
	}
	while ((bool == 0) || (seqnum > totnum));
	scanf("%c", &dummy);

	fprintf(destali, "  %d", seqnum);	/* Prints the number of
						 * sequences.. */
	remember = ftell(destali);	/* ..remembers where the pointer is.. */
	fprintf(destali, "    LENGTH\n");	/* ..and keeps ten spaces for
						 * the length to be filled in
						 * later  */
	do {
		printf("Nucleic acid or Protein coding sequence?          (n/p): ");
		scanf("%c", &nucchar);
		if ((nucchar == 'n') || (nucchar == 'N')) {
			nucoraa = 1;
			bool = 1;
		} else if ((nucchar == 'p') || (nucchar == 'P')) {
			nucoraa = 3;
			bool = 1;
		} else {
			printf("Try again:\n");
			scanf("%c", &dummy);
			bool = 0;
		}
	}
	while (bool == 0);
	scanf("%c", &dummy);

	if (nucoraa == 3) {
		do {
			printf("Nuclear or mitochondrial genetic code?            (n/m): ");
			scanf("%c", &nucchar);
			if ((nucchar == 'n') || (nucchar == 'N')) {
				nucmito = 0;
				bool = 1;
			} else if ((nucchar == 'm') || (nucchar == 'M')) {
				nucmito = 1;
				bool = 1;
			} else {
				printf("Try again:\n");
				scanf("%c", &dummy);
				bool = 0;
			}
		}
		while (bool == 0);
		scanf("%c", &dummy);

		do {
			printf("\nEnter a number between 1 and 5, for the\n");
			printf("   codon position you wish to analyze:\n");
			printf("   1 for first, 2 for second, 3 for third\n");
			printf("   4 for first plus second, 5 for all.            (1-5): ");
			scanf("%d", &position);
			if ((position < 1) || (5 < position)) {
				printf("Bull....  Try again\n");
			}
		}
		while ((position < 1) || (5 < position));
		scanf("%c", &dummy);

		if ((position != 2) && (position != 3)) {
			do {
				printf("Conversion of first positions to degenerate base? (y/n): ");
				scanf("%c", &nucchar);
				if ((nucchar == 'y') || (nucchar == 'Y')) {
					yorm = 1;
					bool = 1;
				} else if ((nucchar == 'n') || (nucchar == 'N')) {
					yorm = 0;
					bool = 1;
				} else {
					printf("Try again:\n");
					scanf("%c", &dummy);
					bool = 0;
				}
			}
			while (bool == 0);
			scanf("%c", &dummy);
		}
	}
	do {
		printf("Use nAmes or nUmbers as identifiers?              (a/u): ");
		scanf("%c", &nucchar);
		if ((nucchar == 'a') || (nucchar == 'A')) {
			namebool = 0;
			bool = 1;
		} else if ((nucchar == 'u') || (nucchar == 'U')) {
			namebool = 1;
			bool = 1;
		} else {
			printf("Try again:\n");
			scanf("%c", &dummy);
			bool = 0;
		}
	}
	while (bool == 0);
	scanf("%c", &dummy);

	printf("\nNucleotide sequences in:           %12s\n", nucname);
	if (nucoraa == 1) {
		printf("Nucleotide alignment source:       %12s\n", aliname);
		printf("Nucleotide alignment destination:  %12s\n", destname);
	} else if (nucoraa == 3) {
		printf("Amino acid alignment source:       %12s\n", aliname);
		printf("Nucleotide alignment destination:  %12s\n", destname);
		if (position == 1) {
			printf("First ");
		} else if (position == 2) {
			printf("Second ");
		} else if (position == 3) {
			printf("Third ");
		} else if (position == 4) {
			printf("First plus second ");
		} else if (position == 5) {
			printf("All ");
		}
		printf("codon positions will be used.\n");

		if ((yorm == 0) && (position != 2) && (position != 3)) {
			if (nucmito == 0) {
				printf("No conversion of L and R 1st positions.\n");
			} else if (nucmito == 1) {
				printf("No conversion of L 1st positions.\n");
			}
		} else if ((yorm == 1) && (position != 2) && (position != 3)) {
			if (nucmito == 0) {
				printf("L and R 1st positions will be converted to Y and M.\n");
			} else if (nucmito == 1) {
				printf("L 1st positions will be converted to Y.\n");
			}
		}		/* end of if yorm == 0 */
	}			/* end of nucoraa == 3 */
	if (namebool == 0) {
		printf("Names ");
	} else if (namebool == 1) {
		printf("Numbers ");
	}
	printf("will be used to identify sequences.\n\n");

	return 0;
}


findalignment()
{				/* finds the string "ALIGNMENT" in the
				 * alignment */
	char            testchar[9], dummy, seqname[200];
	int             bntms, count, eqwpv;

	for (count = 0; count < totnum; count++) {
		fscanf(sourceali, "%d", &taxnum);
		dummy = fgetc(sourceali);
		fgets(seqname, 200, sourceali);
		bntms = 0;
		while ((seqname[bntms] != '\n') && (bntms < 10)) {
			if (seqname[bntms] != '\n') {
				taxname[bntms][taxnum] = seqname[bntms];
			}
			bntms++;
		}
		if (bntms < 10) {
			for (eqwpv = bntms; eqwpv < 10; eqwpv++)
				taxname[eqwpv][taxnum] = ' ';
		}
	}

	testchar[0] = '\n';
	while (strncmp(testchar, "ALIGNMENT", 9) != 0) {
		for (count = 0; count < 9; count++) {
			testchar[count] = fgetc(sourceali);
			if (testchar[count] == '\n') {
				break;
			}
		}
		dummy = testchar[count];
		if (dummy == '\n') {
			continue;
		}
		while ((dummy = fgetc(sourceali)) != '\n');
	}
	return 0;
}

findbeginning()
{				/* finds the first line of sequence alignment */
	char            testchar;

	while ((testchar = fgetc(sourceali)) != '\n');
	while ((testchar != 45 && testchar < 65) ||
	       (testchar > 91 && testchar != '{'))
		testchar = fgetc(sourceali);
	begin = ftell(sourceali) - 1;
	/* pointer will be put to last \n before first amino acid */
	return 0;
}

countblocks()
{				/* counts the number of seqeunce blocks in
				 * the alignment */
	char            dummy;
	int             dumnuc;

	dumnuc = -1;
	do {
		dummy = fgetc(sourceali);
		while ((dummy != '*') && (dummy != EOF))
			dummy = fgetc(sourceali);
		dumnuc++;
		while ((dummy != '\n') && (dummy != EOF))
			dummy = fgetc(sourceali);
	}
	while (dummy != EOF);
	return dumnuc;
}

jump(howmany)
{				/* jumps howmany lines in sourceali */
	char            linedump, aminocheck;
	int             eqwpv;
	long            lastpos;
	eqwpv = 0;

	do {
		linedump = fgetc(sourceali);
		/*
		 * checks if the first character on the line is an amino acid
		 * or a [
		 */
		if ((linedump < 65 || 91 < linedump) && (linedump != 45))
			/* if it's not an amino acid.. */
		{
			--eqwpv;/* ..do not increase counter */
			if (linedump == '{')	/* if it's an ignore
						 * character.. */
				while ((linedump = fgetc(sourceali)) != '}');
			/* ..go to resume char */
		}
		if (linedump != '\n')	/* if not at end of line.. */
			while ((linedump = fgetc(sourceali)) != '\n');
		/* ..goes to end of line */
		eqwpv++;
	}
	while (eqwpv < howmany);

	do {			/* this loop ties up the end by ensuring that
		 * the next line starts with an amino acid or a [ */
		lastpos = ftell(sourceali);	/* remembers the position in
						 * the sequence */
		aminocheck = fgetc(sourceali);	/* amino acid, hyphen, or '['
						 * ? */
		if ((aminocheck < 65 || 91 < aminocheck) && aminocheck != 45) {
			if (aminocheck == '{')	/* if it's an ignore
						 * character.. */
				while ((aminocheck = fgetc(sourceali)) != '}');
			/* ..goes to resume char */
			if (aminocheck != '\n')	/* if not at end of line ... */
				while ((aminocheck = fgetc(sourceali)) != '\n');
			/* ..goes to end of line */
		}
	}
	while ((aminocheck < 65 || 91 < aminocheck) && (aminocheck != 45));

	fseek(sourceali, lastpos, 0);	/* loop is ended by amino acid; now
					 * go back to the beginning of the
					 * current line */
	return 0;
}				/* end of jump */


findname()
{				/* finds the number of the sequence line that
				 * it's looking at and prints it or the
				 * corresponding name at the beginning of the
				 * destfile's corresponding sequence */
	char            dummy;
	int             count, totnumloop, testcount, namefield;
	long            currpos;

	amslr = 0;
	dummy = fgetc(sourceali);
	if (dummy == '{') {	/* if it's an ignore char.. */
		do
			dummy = fgetc(sourceali);
		while (dummy != '}');
		do
			dummy = fgetc(sourceali);
		while (dummy != '\n');
		begin = ftell(sourceali);	/* for putnucs (\n) */
		dummy = fgetc(sourceali);	/* pointer on first aa */
	}			/* ..read until on the first amino acid after
				 * the resume char */
	while (dummy == 32 || dummy == 45 ||
	       (65 <= dummy && dummy <= 91) || dummy == 93) {
		dummy = fgetc(sourceali);
		amslr++;
	}

/* USER: IF YOUR ALIGNMENT ONLY HAS THE TAXON NUMBER AT THE END OF EACH
 * LINE, AND NOT ALSO THE POSITION NUMBER, DELETE ONE OF THE FOLLOWING
 * TWO LINES BEFORE COMPILING! */

	fscanf(sourceali, "%d", &taxnum);	/* dumps the length number */
	fscanf(sourceali, "%d", &taxnum);	/* correct number of taxon */

	if (namebool == 1) {	/* if use numbers */
		fprintf(destali, "%d", taxnum);
		if (taxnum < 10) {
			testcount = 1;
		}
		 /* writes the number ... */
		else if ((taxnum < 100) && (taxnum > 9)) {
			testcount = 2;
		} else if ((taxnum < 512) && (taxnum > 99)) {
			testcount = 3;
		} else {
			printf("Number of sequence out of range.\n");
			exit(0);
		}

		for (namefield = testcount; namefield < 10; namefield++) {
			fputc(' ', destali);
		}		/* ..and fills up the namefield to 10 chars */
	} else if (namebool == 0) {	/* if use names, not numbers.. */
		for (count = 0; count < 10; count++) {
			fprintf(destali, "%c", taxname[count][taxnum]);
		}		/* ..writes the name */
	}
	while (dummy != '\n') {
		dummy = fgetc(sourceali);
	}

	return 0;
}

findnucseq()
{				/* first puts file-pointer at > of sequence
				 * and then puts it at first nucleotide */
	char            dummy;
	int             countseq;

	fseek(sourcenuc, 0, 0);

	countseq = -1;		/* sequence numbers start with 0 */
	while (countseq < taxnum) {
		countseq++;
		while ((dummy = fgetc(sourcenuc)) != '>');	/* marker */
		while ((dummy = fgetc(sourcenuc)) != '\n');	/* marker line end */
		while ((dummy = fgetc(sourcenuc)) != '\n');	/* name line end */
	}
	return 0;
}


goback()
{
	fseek(sourceali, begin, 0);
	return 0;
}


countnuc(nuc)
{
	actualcount++;
	switch (nuc) {
	case 'G':
		count[0]++;
		break;
	case 'A':
		count[1]++;
		break;
	case 'T':
	case 'U':
		count[2]++;
		break;
	case 'C':
		count[3]++;
		break;
	case 'Y':
		count[4]++;
		break;
	case 'R':
		count[5]++;
		break;
	case 'M':
		count[6]++;
		break;
	case 'N':
		count[7]++;
		break;
	default:
		count[8]++;
		break;
	}
	return 0;
}


putnucs()
{
char            base, dummy;
int             zahl;

while ((dummy = fgetc(sourceali)) != '\n') {
 if (dummy == 91) {	/* if it's a start-exclude character ([) */
  putnucbool = 1;
 } else if (dummy == 93) { /* if it's a stop-exclude character (]) */
  putnucbool = 0;
 } else if (65 <= dummy && dummy <= 90) { /* if it's an authentic
                                           * sequence character */
  for (zahl = 0; zahl < nucoraa; zahl++) {
   base = fgetc(sourcenuc);
   if ((65 <= base && base <= 90) || (97 <= base && base <= 122)) {
    if (putnucbool == 0) {
     if (nucoraa == 3) { /* if it's an amino acid sequence */
      if (yorm == 1) { /* if conversion is desired */
       if ((dummy == 'L') && (zahl == 0)) {
         base = 'Y'; /* if amino acid is leucine */
       }
       if (nucmito == 0) { /* if code is nuclear */
        if ((dummy == 'R') && (zahl == 0)) {
         base = 'M'; /* if amino acid is arginine */
        }
       }
      }
      if (97 <= base && base <= 122) {
       base -= 32;
      }
     }
     if ((position == 5) || ((position == 1) && (zahl == 0)) ||
        ((position == 2) && (zahl == 1)) ||
        ((position == 3) && (zahl == 2)) ||
        ((position == 4) && ((zahl == 0) || (zahl == 1)))) {
      countnuc(base);
      newcheck++;
      fputc(base, destali);
      if (newcheck % 60 == 0)
       fputc('\n', destali);
     }
    } /* end of if conditional on putnucbool */
   } /* end of if conditional on authentic nucleotide */
   else if (base == '*') {
    printf("CRASH! not enough nucleotides.\n");
    exit(0);
   } else
     --zahl;
  } /* end of loop counted by nucoraa and zahl */
 } /* end of if conditional on authentic sequence character in alignment */
 else if (dummy == 45) { /* if it is a hyphen */
  if (putnucbool == 0) {
   for (zahl = 0; zahl < nucoraa; zahl++) {
    if ((position == 5) || ((position == 1) && (zahl == 0)) ||
       ((position == 2) && (zahl == 1)) ||
       ((position == 3) && (zahl == 2)) ||
       ((position == 4) && ((zahl == 0) || (zahl == 1)))) {
     newcheck++;
     fputc('-', destali);
     if (newcheck % 60 == 0)
      fputc('\n', destali);
    }
   }
  }
 }
}  /* end of while */
return 0;
}  /* end of putnucs */


tailcheck(seq)
{
	char            checkbool, pacman;
	checkbool = 0;
	do {
		pacman = fgetc(sourcenuc);
		if ((65 <= pacman && pacman <= 90) ||
		    (97 <= pacman && pacman <= 122)) {
			checkbool = 1;
			if (checkbool == 1)
				printf("%c", pacman);
		}
	}
	while (pacman != '*');
	if (checkbool == 1) {
		printf("\nWARNING. Sequence # %d in %s contains the above\n", seq, nucname);
		printf(" bases beyond the number of bases it should contain\n");
		printf(" according to the alignment in %s.\n", aliname );
		printf(" Please compare the two files.\n\n");
	}
}

summer()
{
	int             cou;

	count[1] += (2 * count[6] / 3);
	count[2] += (count[4] / 3);
	count[3] += ((2 * count[4] + count[6]) / 3);
	printf("Frequencies of A, C, G, T:\n");
	printf("%10.5lf", count[1] / actualcount);
	printf("%10.5lf", count[3] / actualcount);
	printf("%10.5lf", count[0] / actualcount);
	printf("%10.5lf\n", count[2] / actualcount);
	return 0;
}
/* end of program */

