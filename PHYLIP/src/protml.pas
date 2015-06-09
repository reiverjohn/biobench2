PROGRAM PROTML (* Maximum Likelihood Inference of Protein Phylogeny *)
   (seqfile, lklfile, tpmfile, output);
(* (seqfile, output); *)

(*            PROTML  Ver 1.01    January 6, 1993          *)
(*               J.Adachi  and  M.Hasegawa                 *)
(*        The Institute of Statistical Mathematics         *)
(*     4-6-7 Minami-Azabu, Minato-ku, Tokyo 106, Japan     *)

 CONST
   maxsp     =      20; (* maximum number of species                 *)
   maxnode   =      37; (* maxsp * 2 - 3                             *)
   maxpair   =     190; (* maxsp * (maxsp-1) / 2                     *)
   maxsite   =     500; (* maximum number of sites                   *)
   maxptrn   =     500; (* maximum number of different site patterns *)
   maxtree   =      15; (* maximum number of user trees              *)
   maxsmooth =      30; (* number of smoothing algorithm             *)
   maxiterat =      10; (* number of iterates of Newton method       *)
   epsilon   =  1.0e-5; (* stopping value of error                   *)
   minarc    =  1.0e-3; (* lower limit on number of subsutitutions   *)
   maxarc    =  3.0e+2; (* upper limit on number of subsutitutions   *)
   prprtn    =       1; (* proportion of branch length               *)
   maxboot   =   10000; (* number of bootstrap resamplings           *)
   maxexe    =       1; (* number of jobs                            *)
   maxline   =      60; (* length of sequence output per line        *)
   maxname   =      10; (* max. number of characters in species name *)
   maxami    =      20; (* number of amino acids                     *)
   minreal   = 1.0e-55; (* if job is in underflow error,
                           then increase this value                  *)
   seqfname  = 'seqfile';         (* input file of sequences data    *)
   tpmfname  = 'default_not_use'; (* input file of trans probability *)
   lklfname  = 'default_not_use'; (* output file of ln likelihood    *)
(* maxnode   = maxsp*2 -3;            *)
(* maxpair   = maxsp*(maxsp-1) DIV 2; *)

 TYPE
   nacidty      = (ami);
   amity        = (AA,RR,NN,DD,CC,QQ,EE,GG,HH,II,
                   LL,KK,MM,FF,PP,SS,TT,WW,YY,VV);
   aryamity     = ARRAY[AA..VV] OF REAL;
   iaryamity    = ARRAY[AA..VV] OF INTEGER;
   daryamity    = ARRAY[AA..VV] OF aryamity;
   taryamity    = ARRAY[AA..VV] OF daryamity;
   niaryamity   = ARRAY[1..maxami] OF INTEGER;
   naryamity    = ARRAY[1..maxami] OF REAL;
   ndaryamity   = ARRAY[1..maxami] OF naryamity;
   probamity    = ARRAY[1..maxptrn] OF aryamity;
   charseqty    = ARRAY[0..maxsp] OF ARRAY [1..maxsite] OF CHAR;
   namety       = ARRAY[1..maxname] OF CHAR;
   nodepty      = ^node;
   arynodepty   = ARRAY[0..maxnode] OF nodepty;
   lengthty     = ARRAY[1..maxnode] OF REAL;
   ilengthty    = ARRAY[1..maxnode] OF INTEGER;
   dlengthty    = ARRAY[1..maxnode] OF lengthty;
   sitety       = ARRAY[1..maxsite] OF INTEGER;
   itreety      = ARRAY[0..maxtree] OF INTEGER;
   rtreety      = ARRAY[0..maxtree] OF REAL;
   lptrnty      = ARRAY[1..maxptrn] OF REAL;
   ltrptty      = ARRAY[0..maxtree] OF lptrnty;
   spity        = ARRAY[1..maxsp] OF INTEGER;
   rpairty      = ARRAY[1..maxpair] OF REAL;
   rnodety      = ARRAY[1..maxnode] OF REAL;
   dpanoty      = ARRAY[1..maxpair] OF rnodety;
   dnopaty      = ARRAY[1..maxnode] OF rpairty;
   dnonoty      = ARRAY[1..maxnode] OF rnodety;
   umbty        = ARRAY[0..maxsp] OF BOOLEAN;
   lspty        = ARRAY[0..maxsp] OF INTEGER;
   longintty    = ARRAY[0..5] OF INTEGER;

   tree = RECORD             (* tree topology data *)
      brnchp   : arynodepty;    (* point to node                *)
      lklhd    : REAL;          (* ln likelihood of this tree   *)
      vrilkl   : REAL;          (* variance of likelihood       *)
      vrilnga  : lengthty;      (* variance of length(branch)   *)
      vrilnga2 : lengthty;      (* variance2 of length(branch)  *)
      blklhd   : lengthty;      (* ln likelihood(branch)        *)
      startp   : nodepty;       (* point to a basal node        *)
      aic      : REAL;          (* Akaike Information Criterion *)
      cnvrgnc  : BOOLEAN        (* convergence of branch length *)
   END;

   node = RECORD             (* a node of tree *)
      isop,                     (* pointer to next subnode           *)
      kinp   : nodepty;         (* pointer to an ascendant node      *)
      diverg : INTEGER;         (* number of divergences             *)
      number : INTEGER;         (* node number                       *)
      namesp : namety;          (* species name. if this node is tip *)
      descen : spity;           (* number of descenant nodes         *)
      CASE nadata : nacidty OF
       ami :
        ( prba : probamity;     (* a partial likelihood *)
          lnga : REAL )         (* branch length        *)
   END;

 VAR
   seqfile,                  (* data file of amino acid sequences    *)
   lklfile,                  (* output file of ln likelihood of site *)
   tpmfile                   (* transition probability matrix data   *)
    : TEXT;
   numsp,                    (* number of species                   *)
   ibrnch1,                  (* first numbering of internal branch  *)
   ibrnch2,                  (* last numbering of internal branch   *)
   numpair,                  (* number of species pairs             *)
   endsite,                  (* number of sites                     *)
   endptrn,                  (* max number of site patterns         *)
   numnw,                    (* curent numbering of Newton method   *)
   numsm,                    (* curent numbering of smoothing       *)
   numexe,                   (* curent numbering of executes        *)
   numtrees,                 (* total number of tree topologies     *)
   notree,                   (* current numbering of tree           *)
   maxltr,                   (* numbering of max ln likelihood tree *)
   minatr,                   (* numbering of min AIC tree           *)
   stage                     (* numbering of decomposable stage     *)
    : INTEGER;
   maxlkl,                   (* max ln likelihood of trees *)
   minaic                    (* min AIC of trees           *)
    : REAL;
   chs                       (* current character of seqfile *)
    : CHAR;
   normal,                   (* break out error                   *)
   usertree,                 (* U option  designate user trees    *)
   semiaoptn,                (* S option  semi-auto decomposition *)
   bootsoptn,                (* B option  bootstrap probability   *)
   writeoptn,                (* W option  print output data       *)
   debugoptn,                (* D option  print debug data        *)
   putlkoptn,                (* P option *)
   firstoptn,                (* F option *)
   lastoptn,                 (* L option *)
   readtoptn,                (* R option *)
   smoothed                  (*          *)
    : BOOLEAN;

   ctree       : tree;       (* current tree                         *)
   freqa       : aryamity;   (* amino acid frequency(acid type)      *)
   chsequen    : charseqty;  (* acid sequence data(species,site)     *)
   weight,                   (* total number of a new site(new site) *)
   alias,                    (* number of old site(new site)         *)
   weightw     : sitety;     (* work area of weight                  *)
   freqdyhf,
   eval,
   eval2       : aryamity;
   ev,
   iev         : daryamity;
   coefp       : taryamity;
   freenode    : nodepty;
   paratree    : itreety;   (* number of parameters         *)
   lklhdtree,               (* ln likelihood                *)
   lklboottree,             (* ln likelihood (bootstrap)    *)
   aictree,                 (* Akaike Information Criterion *)
   boottree    : rtreety;   (* bootstrap probability        *)
   lklhdtrpt   : ltrptty;   (* ln likelihood                *)

  (*** print line ***)
   PROCEDURE printline (* print '-' line to standard output *)
    ( num : INTEGER );
   VAR i : INTEGER;
   BEGIN (* printline *)
      write(output,' ');
      FOR i := 1 TO num DO write(output,'-');
      writeln(output);
   END;  (* printline *)

  (*******************************************************)
  (*****    READ DATA, INITIALIZATION AND SET UP     *****)
  (*******************************************************)

     (*** GET NUMBERS OF SPECIES AND SITES ***)
      PROCEDURE getnums; (* get species-number,sites-number from file*)
      (* input number of species, number of sites *)
      BEGIN (* getnums *)
         read(seqfile, numsp, endsite);
         IF usertree THEN
         BEGIN
            ibrnch1 := numsp + 1;
            ibrnch2 := numsp*2 -3;
         END
         ELSE
         BEGIN
            ibrnch1 := numsp + 1;
            ibrnch2 := numsp;
         END;
      (* numpair  := numsp*(numsp-1) DIV 2; *)
         numpair  := trunc ( (numsp*numsp-numsp) / 2.0 );
         IF numsp > maxsp THEN
            writeln(output, 'TOO MANY SPECIES: adjust CONSTants');
         IF endsite > maxsite THEN
            writeln(output, 'TOO MANY SITES:   adjust CONSTants');
         normal := (numsp <= maxsp) AND (endsite <= maxsite);
      END; (* getnums *)

     (*** CONVERT CHARACTER TO UPPER CASE ***)
      PROCEDURE uppercase (* convert character to upper case *)
       (VAR chru : CHAR);
      (* convert chru to upper case -- either ASCII or EBCDIC *)
         FUNCTION letter (chru : CHAR) : BOOLEAN;
         (* tests whether chru is a letter, (ASCII and EBCDIC) *)
         BEGIN (* letter *)
            letter := ((chru >= 'A') AND (chru <= 'I'))
                   OR ((chru >= 'J') AND (chru <= 'R'))
                   OR ((chru >= 'S') AND (chru <= 'Z'))
                   OR ((chru >= 'a') AND (chru <= 'i'))
                   OR ((chru >= 'j') AND (chru <= 'r'))
                   OR ((chru >= 's') AND (chru <= 'z'));
         END; (* letter *)
      BEGIN (* uppercase *)
         IF letter(chru) AND
            (('a' > 'A') AND (chru >= 'a') OR
             ('a' < 'A') AND (chru <  'A')) THEN
            chru := CHR( ORD(chru) +ORD('A') -ORD('a') );
      END; (* uppercase *)

     (*** GET OPERATIONAL OPTIONS ***)
      PROCEDURE getoptions; (* get option from sequences file *)
         VAR chrop : CHAR;
      BEGIN (* getoptions *)
         usertree  := FALSE;
         semiaoptn := FALSE;
         firstoptn := FALSE;
         lastoptn  := FALSE;
         putlkoptn := FALSE;
         bootsoptn := FALSE;
         writeoptn := FALSE;
         debugoptn := FALSE;
         readtoptn := FALSE;
         WHILE NOT EOLN(seqfile) DO
         BEGIN
            read(seqfile, chrop);
            uppercase (chrop);
            IF (chrop = 'U') OR (chrop = 'S') OR
               (chrop = 'F') OR (chrop = 'L') OR
               (chrop = 'B') OR (chrop = 'P') OR
               (chrop = 'W') OR (chrop = 'D') OR
               (chrop = 'R') THEN
               CASE chrop OF
                  'U' : usertree  := TRUE;
                  'S' : semiaoptn := TRUE;
                  'F' : firstoptn := TRUE;
                  'L' : lastoptn  := TRUE;
                  'P' : putlkoptn := TRUE;
                  'B' : bootsoptn := TRUE;
                  'W' : writeoptn := TRUE;
                  'D' : debugoptn := TRUE;
                  'R' : readtoptn := TRUE;
               END
            ELSE
               IF chrop <> ' ' THEN
               BEGIN
                  writeln(output, ' BAD OPTION CHARACTER: ', chrop);
                  normal := FALSE;
               END;
         END;
         IF (numexe = 1) AND writeoptn THEN
         BEGIN
         writeln(output);
         write  ('Users Option':15);
         write  ('Usertree  : ':14,usertree:5);
         write  ('Semiaoptn : ':14,semiaoptn:5);
         writeln('Bootsoptn : ':14,bootsoptn:5);
         write  (' ':15);
         write  ('Putlkoptn : ':14,putlkoptn:5);
         write  ('Writeoptn : ':14,writeoptn:5);
         writeln('Debugoptn : ':14,debugoptn:5);
         write  (' ':15);
         write  ('Firstoptn : ':14,firstoptn:5);
         writeln('Lastoptn  : ':14,lastoptn:5);
         writeln;
         END;
      END; (* getoptions *)

     (*** GET AMINO ACID FREQENCY FROM DATASET ***)
      PROCEDURE readbasefreqs; (* get freqency from sequences file *)
         VAR
         ba : amity;
      BEGIN (* readbasefreqs *)
         readln(seqfile);
         FOR ba := AA TO VV DO read(seqfile, freqa[ba]);
      END; (* readbasefreqs *)

     (*** SETUP DATASTRUCTURE OF TREE ***)
      PROCEDURE setuptree (* setup datastructure of tree *)
       (VAR tr : tree);
         VAR
         i, k : INTEGER;
         p : nodepty;
      BEGIN (* setuptree *)
         FOR i := 1 TO numsp DO
         BEGIN
            NEW(p);
            p^.number := i;
            p^.lnga   := 0.0;
            p^.isop   := NIL;
            p^.diverg := 0;
            FOR k := 1 TO numsp DO p^.descen[k] := 0;
            p^.descen[i] := 1;
            tr.brnchp[i] := p;
            p^.nadata := ami;
         END;
         tr.lklhd   := -999999.0;
         tr.aic     := 0.0;
         tr.startp  := tr.brnchp[1];
         tr.cnvrgnc := false;
      END; (* setuptree *)

     (*** GET SEQUENCES AND PRINTOUT SEQUENCES ***)
      PROCEDURE getsequence; (* get and print sequences from file *)
         VAR
         is, js, ks, ls, ms, ns,
         ichr, numchr, maxchr  : INTEGER;
         chrs           : CHAR;
         chrcnt  : array[1..30] OF CHAR;
         ichrcnt : array[1..30] OF INTEGER;
         namechar, existenced : BOOLEAN;
      BEGIN (* getsequence *)
         readln(seqfile);
         IF debugoptn THEN writeln(output);
         FOR is := 1 TO numsp DO
         BEGIN
         (* readln(seqfile); *)
            REPEAT
               IF EOLN(seqfile) THEN readln(seqfile);
               read(seqfile, chrs);
               IF debugoptn THEN
                  IF chrs = ' ' THEN write(output, chrs:1);
            UNTIL (chrs <> ' ');
            IF debugoptn THEN
            BEGIN
               writeln(output,'<SPACE>');
               write(output, chrs:1);
            END;
            namechar := TRUE;
            ns := 0;
            REPEAT
               ns := ns +1;
               IF ns <= maxname THEN
                  ctree.brnchp[is]^.namesp[ns] := chrs;
               IF EOLN(seqfile) THEN
                  namechar := FALSE
               ELSE
               BEGIN
                  read(seqfile, chrs);
                  IF chrs = ' ' THEN namechar := FALSE
               END;
               IF debugoptn THEN
                  IF namechar THEN write(output, chrs:1);
            UNTIL (NOT namechar);
            IF ns < maxname THEN
            BEGIN
               FOR js := ns+1 TO maxname DO
                  ctree.brnchp[is]^.namesp[js] := ' ';
            END;
        (*  FOR js := 1 TO maxname DO
            BEGIN
               read(seqfile, chrs);
               ctree.brnchp[is]^.namesp[js] := chrs;
            END;  *)
            IF debugoptn THEN
            BEGIN
               writeln(output, '<NAME>');
               write(output, chrs:1);
            END;
            FOR js := 1 TO endsite DO IF normal THEN
            BEGIN
               REPEAT
                  IF EOLN(seqfile) THEN readln(seqfile);
                  read(seqfile, chrs);
               UNTIL (chrs <> ' ') AND ((chrs < '0') OR (chrs > '9'));
               uppercase (chrs);
               IF debugoptn THEN write(output, chrs);
               BEGIN
                  IF NOT  ((chrs = 'A') OR (chrs = 'C') OR (chrs = 'D')
                        OR (chrs = 'E') OR (chrs = 'F') OR (chrs = 'G')
                        OR (chrs = 'H') OR (chrs = 'I') OR (chrs = 'K')
                        OR (chrs = 'L') OR (chrs = 'M') OR (chrs = 'N')
                        OR (chrs = 'P') OR (chrs = 'Q') OR (chrs = 'R')
                        OR (chrs = 'S') OR (chrs = 'T') OR (chrs = 'V')
                        OR (chrs = 'W') OR (chrs = 'Y')
                        OR (chrs = 'B') OR (chrs = 'Z') OR (chrs = 'V')
                        OR (chrs = 'X') OR (chrs = '*') OR (chrs = '-'))
                         THEN
                  BEGIN
                     writeln(output);
                     writeln(output, ' WARNING -- BAD AMINO ACID "',
                     chrs, '" AT POSITION', js:5, ' OF SPECIES', is:3);
                     normal := FALSE;
                  END;
               END;
               chsequen[is][js] := chrs;
            END;
            IF debugoptn THEN writeln(output);
         END;
         readln(seqfile);

         FOR ks := 1 TO endsite DO
         BEGIN
            ichr := 0;
            FOR js := 1 TO numsp DO ichrcnt[js] := 0;
            FOR js := 1 TO numsp DO chrcnt[js] := ' ';
            FOR js := 1 TO numsp DO
            BEGIN
              chrs := chsequen[js][ks];
              existenced := FALSE;
              FOR ms := 1 TO ichr DO
              BEGIN
                IF chrs = chrcnt[ms] THEN
                BEGIN
                  ichrcnt[ms] := ichrcnt[ms] +1;
                  existenced := TRUE;
                END;
              END;
              IF NOT existenced THEN
              BEGIN
                ichr := ichr+1;
                chrcnt[ichr] := chrs;
                ichrcnt[ichr] := 1;
              END;
            END;
            maxchr := 0;
            numchr := 0;
            FOR ms := 1 TO ichr DO
            BEGIN
              IF ichrcnt[ms] > maxchr THEN
              BEGIN
                maxchr := ichrcnt[ms];
                numchr := ms;
              END;
            END;
            IF numchr <> 0 THEN chsequen[0][ks] := chrcnt[numchr]
            ELSE                chsequen[0][ks] := '?';
         END;

       IF numexe = 1 THEN
       BEGIN
       writeln(output);
       write  (output, 'Sequences Data':15);
       writeln(output, numsp:5,'Species,':9, endsite:5,'Sites':6);
       writeln(output);
       writeln(output, ' Species',' ':5,'Amino Acid Sequences');
       writeln(output, ' -------',' ':5,'--------------------');
       writeln(output);
       FOR is := 1 TO ((endsite-1) DIV maxline)+1 DO
       BEGIN
          FOR js := 0 TO numsp DO
          BEGIN
             IF js = 0 THEN
             BEGIN
                write(output, ' Consensus'); (* Consensus Consent *)
                FOR ks := 1 TO maxname-7 DO write(output, ' ');

             END
             ELSE
             BEGIN
                write(output, ' ');
                FOR ks := 1 TO maxname DO
                   write(output, ctree.brnchp[js]^.namesp[ks]);
                write(output, '  ');
             END;
             ls := maxline*is;
             IF ls > endsite THEN ls := endsite;
             FOR ks := maxline*(is-1)+1 TO ls DO
                IF normal THEN
                BEGIN
                   IF (js>0)AND(chsequen[js][ks]=chsequen[0][ks]) THEN
                      write(output, '.')
                   ELSE
                      write(output, chsequen[js][ks]);
                   IF ((ks MOD 10) = 0)AND((ks MOD maxline) <> 0) THEN
                      write(output, ' ');
                END;
             writeln(output);
          END;
          writeln(output);
       END;
       END;
      END; (* getsequence *)

     (*************************************)
     (***      READ SEQUENCES DATA      ***)
     (*************************************)

      PROCEDURE getdata; (* read data from sequences file *)
      BEGIN (* getdata *)
         IF numexe = 1 THEN
         BEGIN
         printline( 57 );
         writeln(output,'              PROTML Ver.1.00b',
                         '  in MOLPHY');
         writeln(output,'    Maximum Likelihood Inference',
                         ' of Protein Phylogeny');
         writeln(output,'                based on Dayhoff model');
         printline( 57 );
         END;
         getnums;
         IF normal THEN getoptions;
      (* IF NOT freqsfrom THEN readbasefreqs; *)
         IF numexe = 1 THEN
         BEGIN
            IF normal THEN setuptree (ctree);
         END;
         IF normal THEN getsequence;
      END; (* getdata *)

        (*** SORT OF SEQUENCES ***)
         PROCEDURE sitesort; (* sort sequences *)
            (* Shell sort keeping sites, weights in same order *)
            VAR
            gap, i, j, jj, jg, k, itemp : INTEGER;
            flip, tied : BOOLEAN;
         BEGIN (* sitesort *)
            FOR i := 1 TO endsite DO
            BEGIN
               alias[i]   := i;
               weightw[i] := 1;
            END;
            gap := endsite DIV 2;
            WHILE (gap > 0) DO
            BEGIN
               FOR i := gap TO endsite-1 DO
               BEGIN
                  j    := i - gap + 1;
                  flip := TRUE;
                  WHILE ((j>0) AND flip) DO
                  BEGIN
                     jj   := alias[j];
                     jg   := alias[j+gap];
                     flip := FALSE; (* fel *)
                     tied := TRUE;  (* fel *)
                     k    := 1;
                     WHILE (k <= numsp) AND tied DO
                     BEGIN
                        flip := chsequen[k][jj] > chsequen[k][jg];
                        tied := tied
                            AND (chsequen[k][jj]=chsequen[k][jg]);
                        k    := k + 1;
                     END;
                     IF flip THEN
                     BEGIN
                        itemp        := alias[j];
                        alias[j]     := alias[j+gap];
                        alias[j+gap] := itemp;
                        itemp          := weightw[j];
                        weightw[j]     := weightw[j+gap];
                        weightw[j+gap] := itemp;
                        j := j - gap;
                     END;
                  END;
               END;
               gap := gap DIV 2;
            END;
         END; (* sitesort *)

        (*** COMBINATION OF SITE ***)
         PROCEDURE sitecombine; (* combine sites *)
            (* combine sites that have identical patterns *)
            VAR
            i, j, k : INTEGER;
            tied    : BOOLEAN;
         BEGIN (* sitecombine *)
            i := 1;
            WHILE i < endsite DO
            BEGIN
               j    := i+1;
               tied := TRUE;
               WHILE (j <= endsite) AND tied DO
               BEGIN
                  k    := 1; (* fel *)
                  WHILE (k <= numsp) AND tied DO
                  BEGIN
                     tied := tied AND
                        (chsequen[k][alias[i]]=chsequen[k][alias[j]]);
                     k    := k + 1;
                  END;
                  IF tied AND (weightw[j] > 0) THEN
                  BEGIN
                     weightw[i] := weightw[i] +weightw[j];
                     weightw[j] := 0;
                     alias[j]   := alias[i];
                  END;
                  j := j + 1;
               END;
              i := j - 1;
            END;
         END; (* sitecombine *)

        (*** SCRUNCH OF SITE ***)
         PROCEDURE sitescrunch;
            (* move so positively weighted sites come first *)
            VAR
            i, j, itemp : INTEGER;
            done, found : BOOLEAN;
         BEGIN (* sitescrunch *)
            done := FALSE;
            i    := 1;
            j    := 2;
            WHILE NOT done DO
            BEGIN
               found := FALSE;
               IF weightw[i] > 0 THEN
                  i := i + 1
               ELSE
               BEGIN
                  IF j <= i THEN j := i + 1;
                  IF j <= endsite THEN
                  BEGIN
                     found := FALSE;
                     REPEAT
                        found := weightw[j] > 0;
                        j     := j + 1;
                     UNTIL found OR (j > endsite);
                     IF found THEN
                     BEGIN
                        j := j - 1;
                        itemp    := alias[i];
                        alias[i] := alias[j];
                        alias[j] := itemp;
                        itemp      := weightw[i];
                        weightw[i] := weightw[j];
                        weightw[j] := itemp;
                     END
                     ELSE
                        done := TRUE;
                  END
                  ELSE
                     done := TRUE;
               END;
               done := done OR (i >= endsite);
            END;
         END; (* sitescrunch *)

      (*** PRINT PATTERN ***)
       PROCEDURE printpattern (* print patterned sequences *)
        ( maxorder : INTEGER );
       (* print pattern *)
       VAR i, j, k, l, n, m, big, sml, kweight : INTEGER;
       BEGIN (* printpattern *)
(* *)
          IF debugoptn THEN
          BEGIN
             writeln(output);
             writeln(output,'alias  :' );
             FOR i := 1 TO endptrn DO write (alias[i]:4);
             writeln(output); writeln(output);
             writeln(output,'weight :' );
             FOR i := 1 TO endptrn DO write (weight[i]:4);
             writeln(output); writeln(output);
             writeln(output,'endptrn =',endptrn:5);
             writeln(output); writeln(output);
          END;
          IF debugoptn THEN
          BEGIN
            writeln(output,'  num alias weight  pattern');
            FOR i:=1 TO endptrn DO
            BEGIN
              write(output,i:4,alias[i]:6,weight[i]:7,'   ');
              FOR j := 1 TO numsp DO
                 write(output,chsequen[j][alias[i]]:2);
              writeln(output);
            END;
            writeln(output);
          END;
(* *)
          IF writeoptn THEN
          BEGIN
          writeln(output);
          writeln(output, ' Species',' ':5,'Patternized Sequences');
          writeln(output, ' -------',' ':5,'---------------------');
          writeln(output);
          FOR i := 1 TO ((endptrn-1) DIV maxline)+1 DO
          BEGIN
             l := maxline*i;
             IF l > endptrn THEN l := endptrn;
             FOR j := 1 TO numsp DO
             BEGIN
                write(output, ' ');
                FOR k := 1 TO maxname DO
                   write(output, ctree.brnchp[j]^.namesp[k]);
                write(output, '  ');
                FOR k := maxline*(i-1)+1 TO l DO
                IF normal THEN
                BEGIN
                   IF (j>1) AND
                      (chsequen[j][alias[k]]=chsequen[1][alias[k]]) THEN
                      write(output, '.')
                   ELSE
                      write(output, chsequen[j][alias[k]]);
                   IF ((k MOD 10) = 0) AND ((k MOD maxline) <> 0) THEN
                      write(output, ' ');
                END;  writeln(output);
             END;  writeln(output);

             FOR n := maxorder DOWNTO 1 DO
             BEGIN
                write(output,' ':3);
                FOR k := 1 TO maxname DO write(output,' ');
                big := 1;
                FOR m := 1 TO n DO big := big * 10;
                sml := big DIV 10;
                FOR k := maxline*(i-1)+1 TO l DO
                IF normal THEN
                BEGIN
                   kweight := ( weight[k] MOD big ) DIV sml;
                   IF (kweight > 0) AND (kweight < 10) THEN
                      write(output, kweight:1)
                   ELSE IF kweight = 0 THEN
                   BEGIN
                      IF (weight[k] MOD big) = weight[k] THEN
                         write(output, ' ':1)
                      ELSE IF (weight[k] MOD big) < weight[k] THEN
                         write(output, '0':1)
                      ELSE
                         write(output, '*':1);
                   END
                   ELSE
                      write(output, '?':1);
                   IF ((k MOD 10) = 0) AND ((k MOD maxline) <> 0) THEN
                      write(output,' ');
                END;  writeln(output);
             END;  writeln(output);
          END;
          END;

       END; (* printpattern *)

     (***********************************)
     (***  ARRANGE SITES OF SEQUENCES  ***)
     (***********************************)

      PROCEDURE makeweights; (* condense same site-pattern *)
      (* make up weights vector to avoid duplicate computations *)
      VAR iw, jw, kw, maxorder, maxweight, nweight, nw : INTEGER;
      BEGIN (* makeweights *)
         sitesort;
         sitecombine;
         sitescrunch;
         maxorder  := 0;
         maxweight := 0;
         FOR iw := 1 TO endsite DO
         BEGIN
            weight[iw] := weightw[iw];
            IF weight[iw] > 0 THEN endptrn := iw;
            IF weight[iw] > maxweight THEN
            BEGIN
               maxweight := weight[iw];
               nweight   := weight[iw];
               nw        := 0;
               REPEAT
                  nweight := nweight DIV 10;
                  nw      := nw + 1;
               UNTIL nweight = 0;
               IF nw > maxorder THEN maxorder := nw;
            END;
         END;
         IF endptrn > maxptrn THEN
         BEGIN
            writeln(output, ' TOO MANY PATTERNS: increase',
            ' CONSTant maxptrn to at least ', endptrn:5);
            normal := FALSE;
         END;

         kw := 0;
         FOR iw := 1 TO endptrn DO
         BEGIN
            FOR jw := 1 TO weight[iw] DO
            BEGIN
               kw := kw +1;
               weightw[kw] := iw;
            END;
         END;

         printpattern ( maxorder );

      END; (* makeweights *)

     (*****************************************)
     (***  SET PARTIAL LIKELIHOODS AT TIPS  ***)
     (*****************************************)

      PROCEDURE makevalues; (* set up fractional likelihoods *)
         (* set up fractional likelihoods at tips *)
         VAR
         i, j, k : INTEGER;
         ba      : amity;
      BEGIN (* makevalues *)
         FOR k := 1 TO endptrn DO
         BEGIN
            j := alias[k];
            FOR i := 1 TO numsp DO
            BEGIN
               FOR ba := AA TO VV DO
                  ctree.brnchp[i]^.prba[k][ba] := 0.0;
               CASE chsequen[i][j] OF
                  'A' : ctree.brnchp[i]^.prba[k][AA] := 1.0;
                  'R' : ctree.brnchp[i]^.prba[k][RR] := 1.0;
                  'N' : ctree.brnchp[i]^.prba[k][NN] := 1.0;
                  'D' : ctree.brnchp[i]^.prba[k][DD] := 1.0;
                  'C' : ctree.brnchp[i]^.prba[k][CC] := 1.0;
                  'Q' : ctree.brnchp[i]^.prba[k][QQ] := 1.0;
                  'E' : ctree.brnchp[i]^.prba[k][EE] := 1.0;
                  'G' : ctree.brnchp[i]^.prba[k][GG] := 1.0;
                  'H' : ctree.brnchp[i]^.prba[k][HH] := 1.0;
                  'I' : ctree.brnchp[i]^.prba[k][II] := 1.0;
                  'L' : ctree.brnchp[i]^.prba[k][LL] := 1.0;
                  'K' : ctree.brnchp[i]^.prba[k][KK] := 1.0;
                  'M' : ctree.brnchp[i]^.prba[k][MM] := 1.0;
                  'F' : ctree.brnchp[i]^.prba[k][FF] := 1.0;
                  'P' : ctree.brnchp[i]^.prba[k][PP] := 1.0;
                  'S' : ctree.brnchp[i]^.prba[k][SS] := 1.0;
                  'T' : ctree.brnchp[i]^.prba[k][TT] := 1.0;
                  'W' : ctree.brnchp[i]^.prba[k][WW] := 1.0;
                  'Y' : ctree.brnchp[i]^.prba[k][YY] := 1.0;
                  'V' : ctree.brnchp[i]^.prba[k][VV] := 1.0;
                  'B' : BEGIN
                          ctree.brnchp[i]^.prba[k][DD] := 0.5;
                          ctree.brnchp[i]^.prba[k][NN] := 0.5;
                        END;
                  'Z' : BEGIN
                          ctree.brnchp[i]^.prba[k][EE] := 0.5;
                          ctree.brnchp[i]^.prba[k][QQ] := 0.5;
                        END;
                  'X' : FOR ba := AA TO VV DO
                        ctree.brnchp[i]^.prba[k][ba] := 1.0;
                  '?' : FOR ba := AA TO VV DO
                        ctree.brnchp[i]^.prba[k][ba] := 1.0;
                  '-' : FOR ba := AA TO VV DO
                        ctree.brnchp[i]^.prba[k][ba] := 1.0;
               END;
            END;
         END;
      END; (* makevalues *)

      (*** EMPIRICAL FREQENCIES OF AMINO ACIDS ***)
       PROCEDURE empirifreqsA; (* calculate empirical frequencies *)
          (* Get empirical frequencies of amino acids from the data *)
          VAR
          ia, ja : INTEGER;
          ba     : amity;
          sfreqa : aryamity;
          sum    : REAL;
       BEGIN (* empirifreqsA *)
          FOR ba := AA TO VV DO sfreqa[ba] := 0.0;
          sum := 0.0;
          FOR ia := 1 TO numsp DO
            FOR ja := 1 TO endptrn DO
               FOR ba := AA TO VV DO
                  sfreqa[ba] := sfreqa[ba]
                     +weight[ja]*ctree.brnchp[ia]^.prba[ja][ba];
          FOR ba := AA TO VV DO sum       := sum +sfreqa[ba];
          FOR ba := AA TO VV DO freqa[ba] := sfreqa[ba] /sum;
          IF (numexe = 1) AND writeoptn THEN
          BEGIN
          writeln(output);
          writeln(output,' Total acid :':13,sum:7:1);
          write  (output,'  A,R,N,D,C :':13);
          FOR ba := AA TO CC DO write(output,sfreqa[ba]:7:1);
          writeln(output);
          write  (output,'  Q,E,G,H,I :':13);
          FOR ba := QQ TO II DO write(output,sfreqa[ba]:7:1);
          writeln(output);
          write  (output,'  L,K,M,F,P :':13);
          FOR ba := LL TO PP DO write(output,sfreqa[ba]:7:1);
          writeln(output);
          write  (output,'  S,T,W,Y,V :':13);
          FOR ba := SS TO VV DO write(output,sfreqa[ba]:7:1);
          writeln(output);
          END;
       END; (* empirifreqsA *)

     (*************************************)
     (***      GET FREQUENCY            ***)
     (*************************************)

      PROCEDURE getbasefreqs; (* print frequencies *)
         VAR ba : amity;
      BEGIN (* getbasefreqs *)
         BEGIN
            empirifreqsA;
            IF (numexe = 1) AND writeoptn THEN
            BEGIN
            writeln(output);
            write(output, ' Empirical');      (* IF freqsfrom THEN *)
            writeln(output, ' Amino Acid Frequencies:');
            write  (output,'A,R,N,D,C :':13);
            FOR ba := AA TO CC DO write(output, freqa[ba]:7:3);
            writeln(output);
            write  (output,'Q,E,G,H,I :':13);
            FOR ba := QQ TO II DO write(output, freqa[ba]:7:3);
            writeln(output);
            write  (output,'L,K,M,F,P :':13);
            FOR ba := LL TO PP DO write(output, freqa[ba]:7:3);
            writeln(output);
            write  (output,'S,T,W,Y,V :':13);
            FOR ba := SS TO VV DO write(output, freqa[ba]:7:3);
            writeln(output);
            END;
         END;
      END; (* getbasefreqs *)

   PROCEDURE getinput; (* read data and set up data *)
   (* read input file *)
   BEGIN (* getinput *)
      normal := TRUE;
      IF normal THEN getdata;
      IF normal THEN makeweights;
      IF normal THEN makevalues;
      IF normal THEN getbasefreqs
   END; (* getinput *)

      PROCEDURE printnmtrx ( nmtrx : ndaryamity );
      VAR i, j : INTEGER;
      BEGIN (* printnmtrx *)
         writeln(output,' MATRIX(numtype)');
         FOR i := 1 TO 20 DO
         BEGIN
            FOR j := 1 TO 20 DO
            BEGIN
               write(output, nmtrx[i][j]:7:3);
               IF (j=10)OR(j=20) THEN
                  writeln(output);
            END;  writeln(output);
         END;
      END;  (* printnmtrx *)

   PROCEDURE luinverse (* INVERSION OF MATRIX ON LU DECOMPOSITION *)
   ( VAR omtrx, imtrx : ndaryamity; nmtrx : INTEGER );
   CONST
      eps = 1.0e-20;
   VAR
      i, j, k, l, maxi, idx, ix, jx : INTEGER;
      sum, tmp, maxb, aw : REAL;
      index : niaryamity;
      wk    : ^naryamity;
   BEGIN (* luinverse *)
      NEW(wk);
      aw := 1.0;
      FOR i := 1 TO nmtrx DO
      BEGIN
         maxb := 0.0;
         FOR j := 1 TO nmtrx DO
            IF ABS(omtrx[i,j]) > maxb THEN maxb := ABS(omtrx[i,j]);
         IF maxb = 0.0 THEN
         BEGIN
            writeln('PROC. LUINVERSE:  singular');
         END;
         wk^[i] := 1.0/maxb
      END;
      FOR j := 1 TO nmtrx DO
      BEGIN
         FOR i := 1 TO j-1 DO
         BEGIN
            sum := omtrx[i,j];
            FOR k := 1 TO i-1 DO
               sum := sum -omtrx[i,k]*omtrx[k,j];
            omtrx[i,j] := sum;
         END;
         maxb := 0.0;
         FOR i := j TO nmtrx DO
         BEGIN
            sum := omtrx[i,j];
            FOR k := 1 TO j-1 DO
               sum := sum -omtrx[i,k]*omtrx[k,j];
            omtrx[i,j] := sum;
            tmp := wk^[i] *ABS(sum);
            IF tmp >= maxb THEN
            BEGIN
               maxb := tmp;
               maxi := i
            END
         END;
         IF j <> maxi THEN
         BEGIN
            FOR k := 1 TO nmtrx DO
            BEGIN
               tmp           := omtrx[maxi,k];
               omtrx[maxi,k] := omtrx[j,k];
               omtrx[j,k]    := tmp
            END;
            aw        := -aw;
            wk^[maxi] := wk^[j]
         END;
         index[j] := maxi;
         IF omtrx[j,j] = 0.0 THEN omtrx[j,j] := eps;
         IF j <> nmtrx THEN
         BEGIN
            tmp := 1.0/omtrx[j,j];
            FOR i := j+1 TO nmtrx DO
               omtrx[i,j] := omtrx[i,j]*tmp
         END
      END;
      FOR jx := 1 TO nmtrx DO
      BEGIN
         FOR ix := 1 TO nmtrx DO wk^[ix] := 0.0;
         wk^[jx] := 1.0;
         l := 0;
         FOR i := 1 TO nmtrx DO
         BEGIN
            idx := index[i];
            sum := wk^[idx];
            wk^[idx] := wk^[i];
            IF l <> 0 THEN
               FOR j := l TO i-1 DO
                  sum := sum -omtrx[i,j]*wk^[j]
            ELSE IF sum <> 0.0 THEN
               l := i;
            wk^[i] := sum
         END;
         FOR i := nmtrx DOWNTO 1 DO
         BEGIN
            sum := wk^[i];
            FOR j := i+1 TO nmtrx DO
               sum := sum -omtrx[i,j]*wk^[j];
               wk^[i] := sum/omtrx[i,i]
         END;
         FOR ix := 1 TO nmtrx DO imtrx[ix,jx] := wk^[ix]
      END;
      DISPOSE(wk)
   END;  (* luinverse *)

      PROCEDURE mproduct ( am, bm : ndaryamity;
                           VAR cm : ndaryamity;
                       na, nb, nc : INTEGER );
      VAR
         ia, ib, ic : INTEGER;
         sum        : REAL;
      BEGIN (* mproduct *)
         FOR ia := 1 TO na DO
         FOR ic := 1 TO nc DO
         BEGIN
            sum := 0.0;
            FOR ib := 1 TO nb DO
               sum := sum +am[ia][ib]*bm[ib][ic];
            cm[ia][ic] := sum;
         END;
      END;  (* mproduct *)

      PROCEDURE convermtrx ( amtrx : daryamity;
                         VAR nmtrx : ndaryamity );
      VAR ba1, ba2 : amity;
      BEGIN (* convermtrx *)
         FOR ba1 := AA TO VV DO
            FOR ba2 := AA TO VV DO
               nmtrx[ORD(ba1)+1][ORD(ba2)+1] := amtrx[ba1][ba2];
      END;  (* convermtrx *)

      PROCEDURE reconvermtrx( nmtrx : ndaryamity;
                          VAR amtrx : daryamity );
      VAR ba1, ba2 : amity;
      BEGIN (* reconvermtrx *)
         FOR ba1 := AA TO VV DO
            FOR ba2 := AA TO VV DO
               amtrx[ba1][ba2] := nmtrx[ORD(ba1)+1][ORD(ba2)+1];
      END;  (* reconvermtrx *)

      PROCEDURE printeigen;
      VAR ba1, ba2 : amity;
      BEGIN (* printeigen *)
         writeln(output,' EIGEN VECTOR');
         FOR ba1 := AA TO VV DO
         BEGIN
            FOR ba2 := AA TO VV DO
            BEGIN
               write(output, ev[ba1][ba2]:7:3);
               IF (ba2 = II)OR(ba2 = VV) THEN
                  writeln(output);
            END;  writeln(output);
         END;
         writeln(output,' EIGEN VALUE');
         FOR ba1 := AA TO VV DO
         BEGIN
            write(output, eval[ba1]:7:3);
            IF (ba1 = II)OR(ba1 = VV) THEN
               writeln(output);
         END;
      END;  (* printeigen *)

      PROCEDURE checkevector ( imtrx : ndaryamity; nn :INTEGER );
      VAR i, j : INTEGER;
      BEGIN (* checkevector *)
         FOR i := 1 TO nn DO
         BEGIN
            FOR j := 1 TO nn DO
            BEGIN
               IF i = j THEN
               BEGIN
                  IF ABS(imtrx[i][j] -1.0) > 1.0e-5 THEN
                     writeln(output,
                     ' error1 eigen vector (checkevector)');
               END
               ELSE
               BEGIN
                  IF ABS(imtrx[i][j]) > 1.0e-5 THEN
                     writeln(output,
                     ' error2 eigen vector (checkevector)');
               END;
            END;
         END;
      END;  (* checkevector *)

      PROCEDURE printamtrx ( amtrx : daryamity );
      VAR ba1, ba2 : amity;
      BEGIN (* printamtrx *)
         writeln(output,' MATRIX(amitype)');
         FOR ba1 := AA TO VV DO
         BEGIN
            FOR ba2 := AA TO VV DO
            BEGIN
               IF (amtrx[ba1][ba2] > 0.1e-5) THEN
                  write(output, amtrx[ba1][ba2]:14:8)
               ELSE IF (amtrx[ba1][ba2] <= 0.0) THEN
                  write(output, ' ':6,amtrx[ba1][ba2]:8)
               ELSE
                  write(output, ' ':3,amtrx[ba1][ba2]:11);
               IF (ba2=II)OR(ba2=VV)OR
                  (ba2=CC)OR(ba2=PP) THEN
                  writeln(output);
            END;  writeln(output);
         END;
      END;  (* printamtrx *)

      PROCEDURE printfreqdyhf;
      VAR ba : amity;
      BEGIN (* printfreqdyhf *)
         writeln(output);
         writeln(output, ' Dayhoff''s Amino Acid Frequencies:');
         write  (output,'A,R,N,D,C :':13);
         FOR ba := AA TO CC DO write(output, freqdyhf[ba]:7:3);
         writeln(output);
         write  (output,'Q,E,G,H,I :':13);
         FOR ba := QQ TO II DO write(output, freqdyhf[ba]:7:3);
         writeln(output);
         write  (output,'L,K,M,F,P :':13);
         FOR ba := LL TO PP DO write(output, freqdyhf[ba]:7:3);
         writeln(output);
         write  (output,'S,T,W,Y,V :':13);
         FOR ba := SS TO VV DO write(output, freqdyhf[ba]:7:3);
         writeln(output);
      END;  (* printfreqdyhf *)

     (*************************************)
     (*** TRANSITION PROBABILITY MATRIX ***)
     (*************************************)
     (*** READ MATRIX ( EIGEN VALUE,VECTOR AND FREQUENCY ) ***)
      PROCEDURE readeigenv;
      VAR
         ba1, ba2 : amity;
      BEGIN (* readeigenv *)
         FOR ba2 := AA TO VV DO
         BEGIN
            FOR ba1 := AA TO VV DO
            BEGIN
               IF EOLN(tpmfile) THEN readln(tpmfile);
               read(tpmfile, ev[ba1][ba2]);
            END;
         END;  readln(tpmfile);
         FOR ba1 := AA TO VV DO
         BEGIN
            IF EOLN(tpmfile) THEN readln(tpmfile);
            read(tpmfile, eval[ba1]);
         END;  readln(tpmfile);
         FOR ba1 := AA TO VV DO
         BEGIN
            IF EOLN(tpmfile) THEN readln(tpmfile);
            readln(tpmfile, freqdyhf[ba1]);
         END;
      END;  (* readeigenv *)

     (*** DATA MATRIX ( EIGEN VALUE,VECTOR AND FREQUENCY ) ***)
      PROCEDURE dataeigenv;
      BEGIN (* dataeigenv *)
      ev[AA,AA] := -2.18920E-01;  ev[RR,AA] := -4.01968E-02;
      ev[NN,AA] := -2.27758E-01;  ev[DD,AA] := -4.04699E-01;
      ev[CC,AA] := -3.49083E-02;  ev[QQ,AA] := -2.38091E-01;
      ev[EE,AA] :=  5.24114E-01;  ev[GG,AA] := -2.34117E-02;
      ev[HH,AA] :=  9.97954E-02;  ev[II,AA] := -8.45249E-03;
      ev[LL,AA] :=  5.28066E-03;  ev[KK,AA] :=  2.05284E-02;
      ev[MM,AA] := -1.23853E-02;  ev[FF,AA] := -9.50275E-03;
      ev[PP,AA] := -3.24017E-02;  ev[SS,AA] :=  5.89404E-01;
      ev[TT,AA] := -1.99416E-01;  ev[WW,AA] := -1.52404E-02;
      ev[YY,AA] := -1.81586E-03;  ev[VV,AA] :=  4.98109E-02;

      ev[AA,RR] :=  4.00877E-02;  ev[RR,RR] :=  3.76303E-02;
      ev[NN,RR] :=  7.25090E-01;  ev[DD,RR] := -5.28042E-01;
      ev[CC,RR] :=  1.46525E-02;  ev[QQ,RR] := -8.66164E-02;
      ev[EE,RR] :=  3.20299E-01;  ev[GG,RR] :=  7.74478E-03;
      ev[HH,RR] := -8.90059E-02;  ev[II,RR] := -2.94900E-02;
      ev[LL,RR] := -3.50400E-03;  ev[KK,RR] := -5.22374E-02;
      ev[MM,RR] :=  1.94473E-02;  ev[FF,RR] :=  5.80870E-03;
      ev[PP,RR] :=  1.62922E-02;  ev[SS,RR] := -2.60397E-01;
      ev[TT,RR] :=  4.22891E-02;  ev[WW,RR] :=  2.42879E-03;
      ev[YY,RR] := -1.46244E-02;  ev[VV,RR] := -5.05405E-04;

      ev[AA,NN] := -4.62992E-01;  ev[RR,NN] :=  2.14018E-02;
      ev[NN,NN] :=  1.71750E-01;  ev[DD,NN] :=  1.02440E-01;
      ev[CC,NN] :=  3.17009E-03;  ev[QQ,NN] :=  1.50618E-01;
      ev[EE,NN] := -1.07848E-01;  ev[GG,NN] :=  5.62329E-02;
      ev[HH,NN] := -1.09388E-01;  ev[II,NN] := -6.52572E-01;
      ev[LL,NN] :=  1.46263E-02;  ev[KK,NN] := -5.22977E-02;
      ev[MM,NN] :=  4.61043E-02;  ev[FF,NN] :=  4.16223E-02;
      ev[PP,NN] :=  6.60296E-02;  ev[SS,NN] :=  4.86194E-02;
      ev[TT,NN] :=  3.28837E-01;  ev[WW,NN] := -4.61826E-03;
      ev[YY,NN] := -8.60905E-03;  ev[VV,NN] :=  3.84867E-01;

      ev[AA,DD] :=  2.47117E-02;  ev[RR,DD] := -3.76030E-02;
      ev[NN,DD] :=  5.56820E-01;  ev[DD,DD] :=  5.60598E-02;
      ev[CC,DD] := -2.51572E-02;  ev[QQ,DD] :=  3.40304E-01;
      ev[EE,DD] := -4.10919E-01;  ev[GG,DD] := -7.16450E-02;
      ev[HH,DD] := -2.14248E-01;  ev[II,DD] :=  5.25788E-02;
      ev[LL,DD] := -1.22652E-02;  ev[KK,DD] := -5.68214E-02;
      ev[MM,DD] :=  1.75995E-02;  ev[FF,DD] := -7.65321E-03;
      ev[PP,DD] := -5.79082E-02;  ev[SS,DD] :=  3.84128E-01;
      ev[TT,DD] := -4.36123E-01;  ev[WW,DD] := -1.26346E-02;
      ev[YY,DD] := -3.19921E-03;  ev[VV,DD] :=  2.57284E-02;

      ev[AA,CC] := -2.23607E-01;  ev[RR,CC] := -2.23607E-01;
      ev[NN,CC] := -2.23607E-01;  ev[DD,CC] := -2.23607E-01;
      ev[CC,CC] := -2.23607E-01;  ev[QQ,CC] := -2.23607E-01;
      ev[EE,CC] := -2.23607E-01;  ev[GG,CC] := -2.23607E-01;
      ev[HH,CC] := -2.23607E-01;  ev[II,CC] := -2.23607E-01;
      ev[LL,CC] := -2.23607E-01;  ev[KK,CC] := -2.23607E-01;
      ev[MM,CC] := -2.23607E-01;  ev[FF,CC] := -2.23607E-01;
      ev[PP,CC] := -2.23607E-01;  ev[SS,CC] := -2.23607E-01;
      ev[TT,CC] := -2.23607E-01;  ev[WW,CC] := -2.23607E-01;
      ev[YY,CC] := -2.23607E-01;  ev[VV,CC] := -2.23607E-01;

      ev[AA,QQ] :=  4.04873E-01;  ev[RR,QQ] :=  6.49103E-03;
      ev[NN,QQ] := -1.19891E-01;  ev[DD,QQ] := -3.67766E-03;
      ev[CC,QQ] := -9.46397E-04;  ev[QQ,QQ] := -1.89655E-01;
      ev[EE,QQ] :=  5.91497E-02;  ev[GG,QQ] := -8.40994E-02;
      ev[HH,QQ] :=  8.85589E-02;  ev[II,QQ] := -7.24798E-01;
      ev[LL,QQ] :=  2.13078E-02;  ev[KK,QQ] :=  6.70215E-02;
      ev[MM,QQ] :=  4.33730E-02;  ev[FF,QQ] :=  4.68020E-02;
      ev[PP,QQ] := -7.75734E-02;  ev[SS,QQ] := -6.04936E-02;
      ev[TT,QQ] := -3.20011E-01;  ev[WW,QQ] :=  6.03313E-04;
      ev[YY,QQ] := -9.61127E-03;  ev[VV,QQ] :=  3.47473E-01;

      ev[AA,EE] := -4.99671E-02;  ev[RR,EE] :=  4.11185E-02;
      ev[NN,EE] :=  1.47778E-01;  ev[DD,EE] :=  1.64073E-01;
      ev[CC,EE] := -1.66608E-03;  ev[QQ,EE] := -3.40600E-01;
      ev[EE,EE] := -2.75784E-02;  ev[GG,EE] := -6.14738E-03;
      ev[HH,EE] :=  7.45280E-02;  ev[II,EE] :=  6.19257E-02;
      ev[LL,EE] :=  8.66441E-02;  ev[KK,EE] :=  3.88639E-02;
      ev[MM,EE] := -8.97628E-01;  ev[FF,EE] := -3.67076E-03;
      ev[PP,EE] :=  4.14761E-02;  ev[SS,EE] :=  1.90444E-03;
      ev[TT,EE] := -4.52208E-02;  ev[WW,EE] := -8.08550E-03;
      ev[YY,EE] := -1.28114E-02;  ev[VV,EE] :=  4.57149E-02;

      ev[AA,GG] := -6.76985E-02;  ev[RR,GG] :=  1.07693E-01;
      ev[NN,GG] :=  1.83827E-01;  ev[DD,GG] :=  2.57799E-01;
      ev[CC,GG] := -8.80400E-04;  ev[QQ,GG] := -5.07282E-01;
      ev[EE,GG] := -2.08116E-02;  ev[GG,GG] := -9.91077E-03;
      ev[HH,GG] :=  1.32778E-01;  ev[II,GG] :=  3.30619E-02;
      ev[LL,GG] := -5.54361E-02;  ev[KK,GG] := -7.58113E-02;
      ev[MM,GG] :=  7.67118E-01;  ev[FF,GG] := -8.15430E-03;
      ev[PP,GG] :=  5.85005E-02;  ev[SS,GG] :=  1.59381E-02;
      ev[TT,GG] := -6.67189E-02;  ev[WW,GG] := -9.06354E-03;
      ev[YY,GG] := -6.80482E-03;  ev[VV,GG] := -3.69299E-02;

      ev[AA,HH] :=  2.37414E-02;  ev[RR,HH] := -1.31704E-02;
      ev[NN,HH] :=  1.93358E-02;  ev[DD,HH] :=  2.78387E-02;
      ev[CC,HH] :=  6.35183E-02;  ev[QQ,HH] :=  1.94232E-02;
      ev[EE,HH] :=  2.72120E-02;  ev[GG,HH] :=  3.14653E-02;
      ev[HH,HH] :=  8.67994E-03;  ev[II,HH] :=  3.64294E-03;
      ev[LL,HH] := -1.65127E-02;  ev[KK,HH] :=  1.45011E-02;
      ev[MM,HH] :=  1.53344E-04;  ev[FF,HH] := -7.14500E-02;
      ev[PP,HH] :=  2.47183E-02;  ev[SS,HH] :=  1.83958E-02;
      ev[TT,HH] :=  2.03333E-02;  ev[WW,HH] := -9.90001E-01;
      ev[YY,HH] := -6.85327E-02;  ev[VV,HH] :=  1.15883E-02;

      ev[AA,II] :=  3.86884E-02;  ev[RR,II] :=  6.83612E-02;
      ev[NN,II] :=  5.94583E-02;  ev[DD,II] :=  7.89097E-02;
      ev[CC,II] := -9.49564E-01;  ev[QQ,II] :=  7.75597E-02;
      ev[EE,II] :=  7.93973E-02;  ev[GG,II] :=  5.96830E-02;
      ev[HH,II] :=  5.15293E-02;  ev[II,II] :=  7.67383E-03;
      ev[LL,II] :=  2.21331E-02;  ev[KK,II] :=  8.30887E-02;
      ev[MM,II] :=  3.95955E-02;  ev[FF,II] := -1.11612E-01;
      ev[PP,II] :=  5.03976E-02;  ev[SS,II] :=  1.59569E-02;
      ev[TT,II] :=  3.72414E-02;  ev[WW,II] := -5.52876E-02;
      ev[YY,II] := -1.86929E-01;  ev[VV,II] :=  1.41872E-02;

      ev[AA,LL] :=  5.80295E-02;  ev[RR,LL] :=  1.02103E-01;
      ev[NN,LL] :=  7.00404E-02;  ev[DD,LL] :=  9.76903E-02;
      ev[CC,LL] :=  2.47574E-01;  ev[QQ,LL] :=  7.38092E-02;
      ev[EE,LL] :=  9.10580E-02;  ev[GG,LL] :=  1.00958E-01;
      ev[HH,LL] :=  3.49872E-02;  ev[II,LL] := -1.13481E-01;
      ev[LL,LL] := -2.12223E-01;  ev[KK,LL] :=  9.20232E-02;
      ev[MM,LL] := -1.06117E-01;  ev[FF,LL] := -5.51278E-01;
      ev[PP,LL] :=  8.17221E-02;  ev[SS,LL] :=  7.48437E-02;
      ev[TT,LL] :=  4.59234E-02;  ev[WW,LL] :=  4.34903E-01;
      ev[YY,LL] := -5.43823E-01;  ev[VV,LL] := -6.69564E-02;

      ev[AA,KK] :=  1.68796E-01;  ev[RR,KK] :=  6.98733E-01;
      ev[NN,KK] := -3.64477E-02;  ev[DD,KK] :=  4.79840E-02;
      ev[CC,KK] := -2.83352E-02;  ev[QQ,KK] :=  2.07546E-03;
      ev[EE,KK] :=  6.44875E-02;  ev[GG,KK] := -1.24957E-01;
      ev[HH,KK] := -2.31511E-01;  ev[II,KK] := -1.08720E-01;
      ev[LL,KK] :=  6.27788E-02;  ev[KK,KK] := -4.17790E-01;
      ev[MM,KK] := -1.98195E-01;  ev[FF,KK] := -1.01464E-02;
      ev[PP,KK] := -2.34527E-01;  ev[SS,KK] :=  1.80480E-01;
      ev[TT,KK] :=  2.62775E-01;  ev[WW,KK] := -7.67023E-02;
      ev[YY,KK] :=  8.53465E-03;  ev[VV,KK] := -1.14869E-01;

      ev[AA,MM] :=  1.14030E-02;  ev[RR,MM] :=  1.03453E-01;
      ev[NN,MM] :=  1.03058E-01;  ev[DD,MM] :=  1.25170E-01;
      ev[CC,MM] := -7.83690E-02;  ev[QQ,MM] :=  8.15278E-02;
      ev[EE,MM] :=  1.09160E-01;  ev[GG,MM] :=  9.53752E-02;
      ev[HH,MM] :=  1.40899E-01;  ev[II,MM] := -2.81849E-01;
      ev[LL,MM] := -4.92274E-01;  ev[KK,MM] :=  9.83942E-02;
      ev[MM,MM] := -3.14326E-01;  ev[FF,MM] :=  2.70954E-01;
      ev[PP,MM] :=  4.05413E-02;  ev[SS,MM] :=  5.49397E-02;
      ev[TT,MM] := -2.30001E-04;  ev[WW,MM] := -6.60170E-02;
      ev[YY,MM] :=  5.64557E-01;  ev[VV,MM] := -2.78943E-01;

      ev[AA,FF] :=  1.68903E-01;  ev[RR,FF] := -4.16455E-01;
      ev[NN,FF] :=  1.75235E-01;  ev[DD,FF] := -1.57281E-01;
      ev[CC,FF] := -3.02415E-02;  ev[QQ,FF] := -1.99945E-01;
      ev[EE,FF] := -2.80839E-01;  ev[GG,FF] := -1.65168E-01;
      ev[HH,FF] :=  3.84179E-01;  ev[II,FF] := -1.72794E-01;
      ev[LL,FF] :=  4.08470E-02;  ev[KK,FF] :=  1.09632E-01;
      ev[MM,FF] :=  3.10366E-02;  ev[FF,FF] :=  1.78854E-02;
      ev[PP,FF] := -2.43651E-01;  ev[SS,FF] :=  2.58272E-01;
      ev[TT,FF] :=  4.85439E-01;  ev[WW,FF] :=  1.87008E-02;
      ev[YY,FF] := -8.19317E-02;  ev[VV,FF] := -1.85319E-01;

      ev[AA,PP] :=  1.89407E-01;  ev[RR,PP] := -5.67638E-01;
      ev[NN,PP] := -2.92067E-02;  ev[DD,PP] :=  4.63283E-02;
      ev[CC,PP] := -7.50547E-02;  ev[QQ,PP] := -1.77709E-01;
      ev[EE,PP] :=  2.11012E-02;  ev[GG,PP] :=  4.57586E-01;
      ev[HH,PP] := -2.73917E-01;  ev[II,PP] :=  3.02761E-02;
      ev[LL,PP] := -8.55620E-02;  ev[KK,PP] := -4.68228E-01;
      ev[MM,PP] := -1.49034E-01;  ev[FF,PP] :=  4.32709E-02;
      ev[PP,PP] :=  1.13488E-01;  ev[SS,PP] :=  1.12224E-01;
      ev[TT,PP] :=  8.54433E-02;  ev[WW,PP] :=  1.49702E-01;
      ev[YY,PP] :=  1.44889E-02;  ev[VV,PP] :=  9.94179E-02;

      ev[AA,SS] := -4.51742E-02;  ev[RR,SS] :=  2.72843E-01;
      ev[NN,SS] := -6.34739E-02;  ev[DD,SS] := -2.99613E-01;
      ev[CC,SS] := -1.26065E-02;  ev[QQ,SS] :=  6.62398E-02;
      ev[EE,SS] := -2.93927E-01;  ev[GG,SS] :=  2.07444E-01;
      ev[HH,SS] :=  7.61271E-01;  ev[II,SS] :=  1.22763E-01;
      ev[LL,SS] := -9.12377E-02;  ev[KK,SS] := -1.67370E-01;
      ev[MM,SS] := -7.69871E-02;  ev[FF,SS] :=  3.25168E-02;
      ev[PP,SS] := -8.69178E-02;  ev[SS,SS] := -4.91352E-02;
      ev[TT,SS] := -9.28875E-02;  ev[WW,SS] := -3.40158E-02;
      ev[YY,SS] := -9.70949E-02;  ev[VV,SS] :=  1.69233E-01;

      ev[AA,TT] := -6.57152E-02;  ev[RR,TT] := -2.96127E-01;
      ev[NN,TT] :=  1.40073E-01;  ev[DD,TT] :=  3.67827E-01;
      ev[CC,TT] :=  3.88665E-02;  ev[QQ,TT] :=  3.54579E-01;
      ev[EE,TT] :=  3.95304E-01;  ev[GG,TT] := -1.30872E-01;
      ev[HH,TT] :=  5.22294E-01;  ev[II,TT] := -9.20181E-02;
      ev[LL,TT] :=  1.35098E-01;  ev[KK,TT] := -2.61406E-01;
      ev[MM,TT] := -5.72654E-02;  ev[FF,TT] := -1.06748E-01;
      ev[PP,TT] := -1.86932E-01;  ev[SS,TT] := -8.60432E-02;
      ev[TT,TT] := -1.17435E-01;  ev[WW,TT] :=  4.21809E-02;
      ev[YY,TT] :=  3.66170E-02;  ev[VV,TT] := -1.03287E-01;

      ev[AA,WW] :=  8.12248E-02;  ev[RR,WW] := -5.52327E-02;
      ev[NN,WW] := -5.60451E-02;  ev[DD,WW] := -1.10967E-01;
      ev[CC,WW] := -4.36007E-02;  ev[QQ,WW] :=  7.49285E-02;
      ev[EE,WW] := -6.08685E-02;  ev[GG,WW] := -3.69332E-01;
      ev[HH,WW] :=  1.94411E-01;  ev[II,WW] :=  3.55141E-02;
      ev[LL,WW] := -8.63681E-02;  ev[KK,WW] := -2.09830E-01;
      ev[MM,WW] := -9.78879E-02;  ev[FF,WW] := -7.98604E-02;
      ev[PP,WW] :=  8.35838E-01;  ev[SS,WW] :=  4.95931E-02;
      ev[TT,WW] :=  7.86337E-02;  ev[WW,WW] :=  1.01576E-02;
      ev[YY,WW] :=  8.79878E-02;  ev[VV,WW] :=  7.51898E-02;

      ev[AA,YY] := -3.45102E-02;  ev[RR,YY] :=  6.45521E-02;
      ev[NN,YY] := -2.22780E-02;  ev[DD,YY] := -3.19792E-02;
      ev[CC,YY] :=  5.97900E-02;  ev[QQ,YY] :=  4.40846E-02;
      ev[EE,YY] := -2.90400E-02;  ev[GG,YY] :=  1.32178E-01;
      ev[HH,YY] :=  4.83456E-02;  ev[II,YY] := -3.33680E-01;
      ev[LL,YY] :=  1.96787E-01;  ev[KK,YY] := -1.27120E-02;
      ev[MM,YY] := -1.01360E-02;  ev[FF,YY] :=  5.00463E-01;
      ev[PP,YY] :=  2.65762E-01;  ev[SS,YY] :=  8.18628E-03;
      ev[TT,YY] := -1.15726E-01;  ev[WW,YY] := -3.46187E-02;
      ev[YY,YY] := -5.64673E-01;  ev[VV,YY] := -4.02511E-01;

      ev[AA,VV] := -2.37427E-02;  ev[RR,VV] :=  6.64165E-02;
      ev[NN,VV] := -3.07931E-02;  ev[DD,VV] := -1.35114E-01;
      ev[CC,VV] := -1.26575E-02;  ev[QQ,VV] := -4.34887E-02;
      ev[EE,VV] := -1.42949E-01;  ev[GG,VV] :=  1.85463E-01;
      ev[HH,VV] :=  4.82054E-02;  ev[II,VV] := -2.88402E-01;
      ev[LL,VV] :=  2.93123E-01;  ev[KK,VV] :=  1.32034E-02;
      ev[MM,VV] :=  6.77690E-02;  ev[FF,VV] := -5.43584E-01;
      ev[PP,VV] :=  1.46520E-01;  ev[SS,VV] :=  3.28990E-03;
      ev[TT,VV] := -7.67072E-02;  ev[WW,VV] := -2.03106E-02;
      ev[YY,VV] :=  5.89467E-01;  ev[VV,VV] := -2.68367E-01;


      eval[AA] :=  -2.03036E-02;  eval[RR] :=  -2.33761E-02;
      eval[NN] :=  -1.71812E-02;  eval[DD] :=  -1.82705E-02;
      eval[CC] :=  -1.55431E-15;  eval[QQ] :=  -1.60581E-02;
      eval[EE] :=  -1.34008E-02;  eval[GG] :=  -1.38363E-02;
      eval[HH] :=  -2.36636E-03;  eval[II] :=  -2.68394E-03;
      eval[LL] :=  -2.91392E-03;  eval[KK] :=  -1.13524E-02;
      eval[MM] :=  -4.34547E-03;  eval[FF] :=  -1.06999E-02;
      eval[PP] :=  -5.48466E-03;  eval[SS] :=  -8.98371E-03;
      eval[TT] :=  -7.23244E-03;  eval[WW] :=  -7.37540E-03;
      eval[YY] :=  -7.91291E-03;  eval[VV] :=  -8.24441E-03;


      freqdyhf[AA] := 0.8712669E-01;  freqdyhf[RR] := 0.4090396E-01;
      freqdyhf[NN] := 0.4043229E-01;  freqdyhf[DD] := 0.4687196E-01;
      freqdyhf[CC] := 0.3347348E-01;  freqdyhf[QQ] := 0.3825543E-01;
      freqdyhf[EE] := 0.4953036E-01;  freqdyhf[GG] := 0.8861207E-01;
      freqdyhf[HH] := 0.3361838E-01;  freqdyhf[II] := 0.3688570E-01;
      freqdyhf[LL] := 0.8535736E-01;  freqdyhf[KK] := 0.8048168E-01;
      freqdyhf[MM] := 0.1475275E-01;  freqdyhf[FF] := 0.3977166E-01;
      freqdyhf[PP] := 0.5067984E-01;  freqdyhf[SS] := 0.6957710E-01;
      freqdyhf[TT] := 0.5854172E-01;  freqdyhf[WW] := 0.1049367E-01;
      freqdyhf[YY] := 0.2991627E-01;  freqdyhf[VV] := 0.6471762E-01;
      END;  (* dataeigenv *)

     (*** TRANSITION PROBABILITY MATRIX ***)
      PROCEDURE tprobmtrx ( arc : REAL;
                        VAR tpr : daryamity );
         (* transition probability matrix *)
         VAR
         sum : REAL;
         iba, jba, kba : amity;
         exparc : aryamity;
      (* negative : BOOLEAN; *)
      BEGIN (* tprobmtrx *)
      (* negative := FALSE; *)
         FOR kba := AA TO VV DO
            exparc[kba] := EXP(arc*eval[kba]);
         FOR iba := AA TO VV DO
         BEGIN
            FOR jba := AA TO VV DO
            BEGIN
               sum := coefp[iba][jba][AA]*exparc[AA]
                    + coefp[iba][jba][RR]*exparc[RR]
                    + coefp[iba][jba][NN]*exparc[NN]
                    + coefp[iba][jba][DD]*exparc[DD]
                    + coefp[iba][jba][CC]*exparc[CC]
                    + coefp[iba][jba][QQ]*exparc[QQ]
                    + coefp[iba][jba][EE]*exparc[EE]
                    + coefp[iba][jba][GG]*exparc[GG]
                    + coefp[iba][jba][HH]*exparc[HH]
                    + coefp[iba][jba][II]*exparc[II]
                    + coefp[iba][jba][LL]*exparc[LL]
                    + coefp[iba][jba][KK]*exparc[KK]
                    + coefp[iba][jba][MM]*exparc[MM]
                    + coefp[iba][jba][FF]*exparc[FF]
                    + coefp[iba][jba][PP]*exparc[PP]
                    + coefp[iba][jba][SS]*exparc[SS]
                    + coefp[iba][jba][TT]*exparc[TT]
                    + coefp[iba][jba][WW]*exparc[WW]
                    + coefp[iba][jba][YY]*exparc[YY]
                    + coefp[iba][jba][VV]*exparc[VV];
               IF sum < minreal THEN         (* negative := TRUE; *)
                 tpr[iba][jba] := 0.0                   (* attention *)
               ELSE
                 tpr[iba][jba] := sum;
            END;
         END;
      (* IF negative THEN
         IF debugoptn THEN
         BEGIN
            write(output,' Trans. prob. 1 is negative !');
            writeln(output,'  arc =',arc:10:5);
         END; *)
      END;  (* tprobmtrx *)

     (*** TRANSITION PROBABILITY AND DIFFRENCE MATRIX ***)
      PROCEDURE tdiffmtrx ( arc : REAL;
                        VAR tpr, td1, td2 : daryamity );
      (* transition probability matrix *)
      VAR
         sum, sumd1, sumd2,
         aaa, rrr, nnn, ddd, ccc, qqq, eee, ggg, hhh, iii,
         lll, kkk, mmm, fff, ppp, sss, ttt, www, yyy, vvv : REAL;
         iba, jba, kba : amity;
         exparc : aryamity;
      BEGIN (* tdiffmtrx *)
         FOR kba := AA TO VV DO
            exparc[kba] := EXP(arc*eval[kba]);
         FOR iba := AA TO VV DO
         BEGIN
           FOR jba := AA TO VV DO
           BEGIN
             aaa := coefp[iba][jba][AA]*exparc[AA];
             rrr := coefp[iba][jba][RR]*exparc[RR];
             nnn := coefp[iba][jba][NN]*exparc[NN];
             ddd := coefp[iba][jba][DD]*exparc[DD];
             ccc := coefp[iba][jba][CC]*exparc[CC];
             qqq := coefp[iba][jba][QQ]*exparc[QQ];
             eee := coefp[iba][jba][EE]*exparc[EE];
             ggg := coefp[iba][jba][GG]*exparc[GG];
             hhh := coefp[iba][jba][HH]*exparc[HH];
             iii := coefp[iba][jba][II]*exparc[II];
             lll := coefp[iba][jba][LL]*exparc[LL];
             kkk := coefp[iba][jba][KK]*exparc[KK];
             mmm := coefp[iba][jba][MM]*exparc[MM];
             fff := coefp[iba][jba][FF]*exparc[FF];
             ppp := coefp[iba][jba][PP]*exparc[PP];
             sss := coefp[iba][jba][SS]*exparc[SS];
             ttt := coefp[iba][jba][TT]*exparc[TT];
             www := coefp[iba][jba][WW]*exparc[WW];
             yyy := coefp[iba][jba][YY]*exparc[YY];
             vvv := coefp[iba][jba][VV]*exparc[VV];

             sum := aaa +rrr +nnn +ddd +ccc
                   +qqq +eee +ggg +hhh +iii
                   +lll +kkk +mmm +fff +ppp
                   +sss +ttt +www +yyy +vvv;
             sumd1 :=
                aaa*eval[AA] +rrr*eval[RR] +nnn*eval[NN]
               +ddd*eval[DD] +ccc*eval[CC] +qqq*eval[QQ]
               +eee*eval[EE] +ggg*eval[GG] +hhh*eval[HH]
               +iii*eval[II] +lll*eval[LL] +kkk*eval[KK]
               +mmm*eval[MM] +fff*eval[FF] +ppp*eval[PP]
               +sss*eval[SS] +ttt*eval[TT] +www*eval[WW]
               +yyy*eval[YY] +vvv*eval[VV];
             sumd2 :=
                aaa*eval2[AA] +rrr*eval2[RR] +nnn*eval2[NN]
               +ddd*eval2[DD] +ccc*eval2[CC] +qqq*eval2[QQ]
               +eee*eval2[EE] +ggg*eval2[GG] +hhh*eval2[HH]
               +iii*eval2[II] +lll*eval2[LL] +kkk*eval2[KK]
               +mmm*eval2[MM] +fff*eval2[FF] +ppp*eval2[PP]
               +sss*eval2[SS] +ttt*eval2[TT] +www*eval2[WW]
               +yyy*eval2[YY] +vvv*eval2[VV];
             IF sum < minreal THEN                      (* attention *)
             BEGIN
                tpr[iba][jba] := 0.0;
                td1[iba][jba] := 0.0;
                td2[iba][jba] := 0.0;
             END
             ELSE
             BEGIN
                tpr[iba][jba] := sum;
                td1[iba][jba] := sumd1;
                td2[iba][jba] := sumd2;
             END;
           END;
         END;
      (* IF negative THEN
         IF debugoptn THEN
         BEGIN
            write(output,' Trans. prob. 2 is negative !');
            writeln(output,'  arc =',arc:10:5);
         END; *)
      END;  (* tdiffmtrx *)

 (*** MAKE TRANSITION PROBABILITY MATRIX ***)
  PROCEDURE maketransprob;
     (* make transition probability matrix *)
     VAR
     amatrix, bmatrix, cmatrix : ndaryamity;
     iba, jba, kba : amity;
  BEGIN  (* maketransprob *)
     IF readtoptn THEN readeigenv
     ELSE              dataeigenv;
     IF writeoptn THEN printfreqdyhf;

     IF readtoptn THEN writeln(output,'read trans. matrix');
     IF debugoptn THEN printeigen;

     convermtrx(ev,amatrix);
     cmatrix := amatrix;
     luinverse ( cmatrix, bmatrix, 20 );
   (* IF debugoptn printnmtrx (bmatrix); *)
     mproduct ( amatrix, bmatrix, cmatrix, 20, 20, 20 );
     checkevector ( cmatrix, 20 );
   (* IF debugoptn printnmtrx (cmatrix); *)
     reconvermtrx(bmatrix,iev);
     FOR iba := AA TO VV DO
       FOR jba := AA TO VV DO
         FOR kba := AA TO VV DO
         BEGIN
           coefp[iba][jba][kba] := ev[iba][kba]
                                 *iev[kba][jba];
         END;
     FOR kba := AA TO VV DO
       eval2[kba] := eval[kba]*eval[kba];
  (*
     tprobmtrx ( 100.0, tprobt );
     printamtrx ( tprobt );
  *)
  END;  (* maketransprob *)

  FUNCTION randum (VAR seed : longintty) : REAL;
    (* random number generator -- slow but machine independent *)
  VAR
     i, j, k, sum : INTEGER;
     mult, newseed : longintty;
     x : REAL;
  BEGIN (* randum *)
    mult[0] := 13;
    mult[1] := 24;
    mult[2] := 22;
    mult[3] :=  6;
    FOR i := 0 TO 5 DO newseed[i] := 0;
    FOR i := 0 TO 5 DO
    BEGIN
      sum := newseed[i];
      k := i;
      IF i > 3 THEN k := 3;
      FOR j := 0 TO k DO sum := sum +mult[j]*seed[i-j];
      newseed[i] := sum;
      FOR j := i TO 4 DO
      BEGIN
        newseed[j+1] := newseed[j+1] + newseed[j] DIV 64;
        newseed[j] := newseed[j] MOD 64;
      END;
    END;
    seed := newseed;
    seed[5] := seed[5] MOD 4;
    x := 0.0;
    FOR i := 0 TO 5 DO x := x/64.0 + seed[i];
    x := x / 4.0;
    randum := x;
  END; (* randum *)

  (*******************************************************)
  (*****          ESTIMATE TREE TOPOLOGY             *****)
  (*******************************************************)

     (**************************************)
     (***       READ TREE STRUCTURE      ***)
     (**************************************)

       (*** SERCH OF END CHARACTER ***)
        PROCEDURE serchchend;
        BEGIN (* serchchend *)
          IF (chs <> ',') AND (chs <> ')') THEN
            REPEAT
              IF EOLN(seqfile) THEN readln(seqfile);
              read(seqfile, chs);
            UNTIL (chs = ',') OR (chs = ')');
        END;  (* serchchend *)

       (*** NEXT CHARACTER ***)
        PROCEDURE nextchar ( VAR ch : CHAR );
        BEGIN (* nextchar *)
          REPEAT
            IF EOLN(seqfile) THEN readln(seqfile);
            read(seqfile, ch);
          UNTIL ch <> ' ';
        END;  (* nextchar *)

       (*** CREATE DATASTRUCTURE OF NODE ***)
        PROCEDURE clearnode ( p : nodepty);
        VAR i : INTEGER;
        BEGIN (* clearnode *)
           with p^ DO
           BEGIN
              isop   := nil;
              kinp   := nil;
              diverg := 1;
              number := 0;
              FOR i := 1 TO maxname DO namesp[i] := ' ';
              FOR i := 1 TO numsp DO descen[i] := 0;
              nadata := ami;
              lnga   := 0.0;
           END;
        END;  (* clearnode *)

       (*** JUDGMENT SPEICIES NAME ***)
        PROCEDURE judgspname ( VAR tr:tree; VAR num:INTEGER );
        VAR
          ie      : INTEGER; (* number of name strings *)
          namestr : namety; (* current species name of seqfile *)
          found   : BOOLEAN;
        BEGIN (* judgspname *)
          FOR ie := 1 TO maxname DO namestr[ie] := ' ';
          ie := 1;
          REPEAT
            IF (chs = '_') THEN chs := ' ';
            namestr[ie] := chs;
            IF EOLN(seqfile) THEN readln(seqfile);
            read(seqfile, chs);
            ie := ie + 1;
          UNTIL ((chs=':')OR(chs=',')OR(chs=')')OR(ie>maxname));
          num := 1;
          REPEAT
            found := TRUE;
            FOR ie := 1 TO maxname DO
              found := found AND
                 (namestr[ie]=tr.brnchp[num]^.namesp[ie]);
            IF NOT found THEN num := num + 1;
          UNTIL (num > numsp) OR found;
          IF num > numsp THEN
          BEGIN
            write(output, ' Cannot find species: ');
            FOR ie := 1 TO maxname DO write(output, namestr[ie]);
            writeln(output);
          END;
        END;  (* judgspname *)

      (*** ADD EXTERNAL NODE ***)
        PROCEDURE externalnode ( VAR tr:tree; up:nodepty );
        VAR
          num : INTEGER; (* number of external nodes *)
        BEGIN (* externalnode *)
          judgspname ( tr,num );
          tr.brnchp[num]^.kinp := up;
          up^.kinp := tr.brnchp[num];
          IF tr.startp^.number > num THEN
             tr.startp := tr.brnchp[num];
          tr.brnchp[num]^.lnga := 0.0;
          up^.lnga := 0.0;
        END;  (* externalnode *)

       (*** ADD INTERNAL NODE ***)
        PROCEDURE internalnode ( VAR tr:tree; up:nodepty;
                                 VAR ninode:INTEGER );
        VAR
           np,            (* new pointer to internal node *)
           cp  : nodepty; (* current pointer to internal node *)
           i, nd,
           dvg : INTEGER;
        BEGIN (* internalnode *)
           nextchar(chs);
           IF chs = '(' THEN
           BEGIN
             ninode := ninode + 1;
             nd := ninode;
             IF freenode = NIL THEN
                NEW ( np )
             ELSE
             BEGIN
                np := freenode;
                freenode := np^.isop;
                np^.isop := NIL;
             END;
             clearnode ( np );
             np^.isop := np;
             cp := np;
             np := NIL;
             cp^.number := nd;
             tr.brnchp[ninode] := cp;
             up^.kinp := cp;
             cp^.kinp := up;
             dvg := 0;
             WHILE chs <> ')' DO
             BEGIN
                dvg := dvg +1;
                IF freenode = NIL THEN
                   NEW ( np )
                ELSE
                BEGIN
                   np := freenode;
                   freenode := np^.isop;
                   np^.isop := NIL;
                END;
                clearnode ( np );
                np^.isop := cp^.isop;
                cp^.isop := np;
                cp := np;
                np := NIL;
                cp^.number := nd;
                internalnode ( tr,cp,ninode );
                serchchend;
             END;
             FOR i := 0 TO dvg DO
             BEGIN
                cp^.diverg := dvg;
                cp := cp^.isop;
             END;
             nextchar(chs);          (* *)
           END
           ELSE
           BEGIN
             externalnode ( tr,up );
           END;
        END; (* internalnode *)

     (*** MAKE STRUCTURE OF TREE ***)
      PROCEDURE treestructure ( VAR tr:tree );
      VAR
        i,
        dvg,
        ninode : INTEGER;  (* number of internal nodes *)
        np,                (* new pointer *)
        cp     : nodepty;  (* current pointer to zero node(root) *)
      BEGIN (* treestructure *)
        ninode := numsp;
        nextchar ( chs );
        IF chs = '(' THEN
        BEGIN
           dvg := -1;
           WHILE chs <> ')' DO
           BEGIN
              IF freenode = NIL THEN
                 NEW ( np )
              ELSE
              BEGIN
                 np := freenode;
                 freenode := np^.isop;
                 np^.isop := NIL;
              END;
              clearnode ( np );
              internalnode ( tr,np,ninode );
              dvg := dvg +1;
              IF dvg = 0 THEN
                 np^.isop := np
              ELSE
              BEGIN
                 np^.isop := cp^.isop;
                 cp^.isop := np;
              END;
              cp := np;
              np := NIL;
              serchchend;
           END;
           tr.brnchp[0] := cp^.isop;
           tr.startp := tr.brnchp[numsp];
           FOR i := 0 TO dvg DO
           BEGIN
              cp^.diverg := dvg;
              cp := cp^.isop;
           END;
        END;
        readln(seqfile);
        ibrnch2 := ninode;
      END; (* treestructure *)

     (*** MAKE STAR STRUCTURE OF TREE ***)
      PROCEDURE starstructure ( VAR tr:tree );
      VAR
        i,
        dvg                (* *)
               : INTEGER;
        np,                (* new pointer *)
        cp     : nodepty;  (* current pointer to zero node(root) *)
      BEGIN (* starstructure *)

        dvg := -1;
        FOR i := 1 TO numsp DO
        BEGIN
             IF freenode = NIL THEN
                 NEW ( np )
              ELSE
              BEGIN
                 np := freenode;
                 freenode := np^.isop;
                 np^.isop := NIL;
              END;
              clearnode ( np );

              tr.brnchp[i]^.kinp := np;
              np^.kinp := tr.brnchp[i];

              dvg := dvg +1;
              IF dvg = 0 THEN
                 np^.isop := np
              ELSE
              BEGIN
                 np^.isop := cp^.isop;
                 cp^.isop := np;
              END;
              cp := np;
              np := NIL;
        END;
        tr.brnchp[0] := cp^.isop;
        tr.startp := tr.brnchp[numsp];
        FOR i := 0 TO dvg DO
        BEGIN
           cp^.diverg := dvg;
           cp := cp^.isop;
        END;
        ibrnch2 := numsp;
      END; (* starstructure *)

     (* INSERT BRANCH IN TREE *)
      PROCEDURE insertbranch ( VAR tr:tree;
                               cp1,cp2,bp1,bp2,np1,np2 :nodepty );
      VAR
         i, num,
         dvg : INTEGER;
         ap  : nodepty;
      BEGIN (* insertbranch *)
         i := cp1^.number;
         IF cp1 = tr.brnchp[i] THEN
         BEGIN
            IF i <> 0 THEN
            BEGIN
               ap  := np1;
               np1 := np2;
               np2 := ap;
            END
            ELSE
               tr.brnchp[i] := np2;
         END;
     (*  np2^.number := bp1^.number;
         cp1^.number := np1^.number;
         cp2^.number := np1^.number;  *)
         bp1^.isop := np2;
         IF cp1 = bp2 THEN
            np2^.isop := cp2^.isop
         ELSE
            np2^.isop := cp1^.isop;
         IF cp1 <> bp2 THEN
            bp2^.isop := cp2^.isop;
         np1^.isop := cp1;
         IF cp1 <> bp2 THEN
            cp1^.isop := cp2;
         cp2^.isop := np1;
         tr.brnchp[ibrnch2]^.kinp^.number :=
            tr.brnchp[ibrnch2]^.kinp^.isop^.number;
         ap  := tr.brnchp[ibrnch2];
         num := tr.brnchp[ibrnch2]^.number;
         REPEAT
            ap := ap^.isop;
            ap^.number := num;
         UNTIL ap^.isop = tr.brnchp[ibrnch2];
         dvg := 0;
         ap := np1;
         REPEAT
            ap := ap^.isop;
            dvg := dvg +1;
         UNTIL ap^.isop = np1;
         ap := np1;
         REPEAT
            ap := ap^.isop;
            ap^.diverg := dvg;
         UNTIL ap = np1;
         dvg := 0;
         ap := np2;
         REPEAT
            ap := ap^.isop;
            dvg := dvg +1;
         UNTIL ap^.isop = np2;
         ap := np2;
         REPEAT
            ap := ap^.isop;
            ap^.diverg := dvg;
         UNTIL ap = np2;
      END;  (* insertbranch *)

     (* DELETE BRANCH IN TREE *)
      PROCEDURE deletebranch ( VAR tr:tree; cnode:INTEGER;
                               cp1,cp2,bp1,bp2,np1,np2 :nodepty );
      VAR
      (* i, *)
         dvg : INTEGER;
         ap  : nodepty;
      BEGIN (* deletebranch *)
         IF cnode <> 0 THEN
         BEGIN
            IF tr.brnchp[cnode] = cp1 THEN
            BEGIN
               ap  := np1;
               np1 := np2;
               np2 := ap;
            END
         END
         ELSE
            IF tr.brnchp[cnode] = np2 THEN tr.brnchp[cnode] := cp1;
      (* i := np2^.number;
         IF i = 0 THEN
            IF np2 = tr.brnchp[i] THEN tr.brnchp[i] := cp1
         ELSE
            IF cp1 = tr.brnchp[i] THEN
            BEGIN
               ap  := np1;
               np1 := np2;
               np2 := ap;
            END; *)
      (* cp1^.number := np2^.number;
         cp2^.number := np2^.number; *)
         IF cp1 = bp2 THEN
            cp2^.isop := np2^.isop
         ELSE
            cp2^.isop := bp2^.isop;
         IF cp1 <> bp2 THEN
            cp1^.isop := np2^.isop;
         np1^.isop := NIL;
         IF cp1 <> bp2 THEN
            bp2^.isop := cp2;
         np2^.isop := NIL;
         bp1^.isop := cp1;
         dvg := 0;
         ap := cp1;
         REPEAT
            ap := ap^.isop;
            dvg := dvg +1;
         UNTIL ap^.isop = cp1;
         ap := cp1;
         REPEAT
            ap := ap^.isop;
            ap^.diverg := dvg;
            ap^.number := cnode;
         UNTIL ap = cp1;
         np1^.diverg := 0;
         np2^.diverg := 0;
      END;  (* deletebranch *)

     (*** PRINT STRUCTURE OF TREE ***)
      PROCEDURE printcurtree ( VAR tr:tree );
      VAR
        ap : nodepty;
        i, j, k, num : INTEGER;
      BEGIN (* printcurtree *)
        writeln(output);
        writeln(output,'Structure of Tree');
        writeln(output,'number':7,'kinp':5,'isop':5,'diverg':8,
           'namesp':9,'length':11,'descendant':13);
        FOR i := 1 TO numsp DO
        BEGIN
          ap := tr.brnchp[i];
          write (ap^.number:5);
          write (ap^.kinp^.number:6);
          IF ap^.isop = NIL THEN write ('nil':6)
          ELSE                   write (ap^.isop^.number:6);
          write (ap^.diverg:7);
          write (' ':3);
          FOR j := 1 TO maxname DO write (ap^.namesp[j]:1);
          write (ap^.lnga:8:3);
          write (' ':3);
          FOR j := 1 TO numsp DO write (ap^.descen[j]:2);
          writeln(output);
        END;
        FOR i:=ibrnch1 TO ibrnch2+1 DO
        BEGIN
          IF i = ibrnch2+1 THEN num := 0
          ELSE                  num := i;
          k := 0;
          ap := tr.brnchp[num];
          REPEAT
            k := k +1;
            IF (ap^.number = tr.brnchp[num]^.number)
               AND (k > 1) THEN
               write ('.':5)
            ELSE
               write (ap^.number:5);
            write (ap^.kinp^.number:6);
            IF (ap^.isop^.number = tr.brnchp[num]^.number)
               AND (k > 0) THEN                           (* 1 *)
               write ('.':6)
            ELSE
               write (ap^.isop^.number:6);
            write (ap^.diverg:7);
            write (' ':3);
            FOR j:=1 TO maxname DO write (' ':1);
            write (ap^.lnga:8:3);
            write (' ':3);
            FOR j := 1 TO numsp DO write (ap^.descen[j]:2);
            writeln(output);
            ap := ap^.isop;
          UNTIL ap = tr.brnchp[num];
        END;
      END; (* printcurtree *)

     (**************************************)
     (***  ESTIMATE LIKELIHOOD OF TREE   ***)
     (**************************************)

      PROCEDURE initdescen ( VAR tr:tree );
      VAR
        ap : nodepty;
        i, j, k, num : INTEGER;
      BEGIN (* initdescen *)
        FOR i := 1 TO numsp DO
        BEGIN
          ap := tr.brnchp[i];
          FOR j := 1 TO numsp DO ap^.descen[j] := 0;
          ap^.descen[i] := 1;
        END;
        FOR i:=ibrnch1 TO ibrnch2+1 DO
        BEGIN
          IF i = ibrnch2+1 THEN num := 0
          ELSE                  num := i;
          ap := tr.brnchp[num];
          FOR k:=1 TO tr.brnchp[num]^.diverg+1 DO
          BEGIN
            FOR j := 1 TO numsp DO ap^.descen[j] := 0;
            ap := ap^.isop;
          END;
        END;
      END; (* initdescen *)

      PROCEDURE initnodeA (p : nodepty);
         VAR
         cp     : nodepty;
         i, n   : INTEGER;
         sumprb : REAL;
         ba     : amity;
      BEGIN (* initnodeA *)
         IF p^.isop <> NIL THEN        (* TIP *)
         BEGIN
            cp := p^.isop;
            FOR n := 1 TO p^.diverg DO
            BEGIN
               initnodeA(cp^.kinp);
               cp := cp^.isop;
            END;
            FOR i  := 1 TO endptrn DO
            FOR ba := AA TO VV DO
            BEGIN
               sumprb := 0.0;
               cp := p^.isop;
               FOR n := 1 TO p^.diverg DO
               BEGIN
                  sumprb := sumprb +cp^.kinp^.prba[i][ba];
                  cp := cp^.isop;
               END;
               IF sumprb <> 0.0 THEN
                  p^.prba[i][ba] := sumprb / p^.diverg
               ELSE
                  p^.prba[i][ba] := 0.0;
            END;
            p^.lnga       := 10.0;                      (* attention *)
            p^.kinp^.lnga := 10.0;                      (* attention *)
            FOR i := 1 TO numsp DO
            BEGIN
               cp := p^.isop;
               FOR n := 1 TO p^.diverg DO
               BEGIN
                  IF cp^.kinp^.descen[i] > 0 THEN p^.descen[i] := 1;
                  cp := cp^.isop;
               END;
            END;
        (*  IF debugoptn THEN
            BEGIN
               write(output,p^.kinp^.number:5,'-':2,p^.number:3,' ':3);
               '(':2,q^.number:3,r^.number:3,')':2,
               FOR i := 1 TO numsp DO  write(output,p^.descen[i]:2);
               writeln(output);
            END;  *)
         END;
      END; (* initnodeA *)

      PROCEDURE initbranch ( p : nodepty );
      VAR
         n  : INTEGER;
         cp : nodepty;
      BEGIN (* initbranch *)
         IF p^.isop = NIL THEN        (* TIP *)
         BEGIN
            initnodeA(p^.kinp);
         END
         ELSE
         BEGIN
            cp := p^.isop;
            FOR n := 1 TO p^.diverg DO
            BEGIN
               initbranch (cp^.kinp);
               cp := cp^.isop;
            END;
         END;
      END; (* initbranch *)

      PROCEDURE evaluateA (VAR tr:tree );
         VAR
         arc, lkl, sum, prod, lnlkl : REAL;
         i        : INTEGER;
         p, q     : nodepty;
         xa1, xa2 : aryamity;
         ai, aj   : amity;
         tprobe   : daryamity;
      BEGIN (* evaluateA *)
         p   := tr.startp;
         q   := p^.kinp;
         arc := p^.lnga;
         IF arc < minarc THEN arc := minarc;
         tprobmtrx ( arc, tprobe );
         lkl := 0.0;
         FOR i := 1 TO endptrn DO
         BEGIN
            xa1 := p^.prba[i];
            xa2 := q^.prba[i];
            sum := 0.0;
            FOR ai := AA TO VV DO
            BEGIN
               prod := freqdyhf[ai]*xa1[ai];
               FOR aj := AA TO VV DO
                  sum := sum +prod *tprobe[ai][aj] *xa2[aj];
            END;
            IF sum > 0.0 THEN lnlkl := LN(sum)
            ELSE              lnlkl := 0.0;
            lkl := lkl +lnlkl*weight[i];
            lklhdtrpt[notree][i] := lnlkl;
         END;
         tr.lklhd := lkl;
         tr.aic   := -2.0*lkl +ibrnch2*2;
      END; (* evaluateA *)

      PROCEDURE sublklhdA (p : nodepty);
         VAR
         i, n      : INTEGER;
         arc       : REAL;
         cp, sp    : nodepty;
         ai        : amity;
         tprob     : daryamity;
      BEGIN (* sublklhdA *)
         cp := p^.isop;
         FOR n := 1 TO p^.diverg DO
         BEGIN
            sp := cp^.kinp;
            arc := sp^.lnga;
            tprobmtrx ( arc, tprob );
            FOR i := 1 TO endptrn DO
            BEGIN
               FOR ai := AA TO VV DO
               BEGIN
                  IF n = 1 THEN p^.prba[i][ai] := 1.0;
                  IF p^.prba[i][ai] < minreal THEN
                     p^.prba[i][ai] := 0.0
                  ELSE                                  (* attention *)
                  p^.prba[i][ai] := p^.prba[i][ai]*
                   (  tprob[ai][AA]*sp^.prba[i][AA]
                    + tprob[ai][RR]*sp^.prba[i][RR]
                    + tprob[ai][NN]*sp^.prba[i][NN]
                    + tprob[ai][DD]*sp^.prba[i][DD]
                    + tprob[ai][CC]*sp^.prba[i][CC]
                    + tprob[ai][QQ]*sp^.prba[i][QQ]
                    + tprob[ai][EE]*sp^.prba[i][EE]
                    + tprob[ai][GG]*sp^.prba[i][GG]
                    + tprob[ai][HH]*sp^.prba[i][HH]
                    + tprob[ai][II]*sp^.prba[i][II]
                    + tprob[ai][LL]*sp^.prba[i][LL]
                    + tprob[ai][KK]*sp^.prba[i][KK]
                    + tprob[ai][MM]*sp^.prba[i][MM]
                    + tprob[ai][FF]*sp^.prba[i][FF]
                    + tprob[ai][PP]*sp^.prba[i][PP]
                    + tprob[ai][SS]*sp^.prba[i][SS]
                    + tprob[ai][TT]*sp^.prba[i][TT]
                    + tprob[ai][WW]*sp^.prba[i][WW]
                    + tprob[ai][YY]*sp^.prba[i][YY]
                    + tprob[ai][VV]*sp^.prba[i][VV] );
               END;
            END;
            cp := cp^.isop;
         END;
      END; (* sublklhdA *)

      PROCEDURE branchlengthA ( p:nodepty; VAR it:INTEGER );
         VAR
         i, numloop : INTEGER;
         arc, arcold,
         sum, sumd1, sumd2,
         lkl, lkld1, lkld2,
         prod1, prod2 : REAL;
         done     : BOOLEAN;
         q        : nodepty;
         ai, aj   : amity;
         tprobl, tdiff1, tdiff2 : daryamity;
      BEGIN (* branchlengthA *)
         q    := p^.kinp;
         arc  := p^.lnga;
         done := FALSE;
         it   := 0;
         IF numsm < 3 THEN numloop := 3 ELSE numloop := maxiterat;

         WHILE (it < numloop) AND (NOT done) DO
         BEGIN
            it   := it + 1;
            IF arc < minarc THEN arc := minarc;
            IF arc > maxarc THEN arc := minarc;         (* attention *)
            tdiffmtrx ( arc, tprobl, tdiff1, tdiff2 );
            lkl := 0.0; lkld1 := 0.0; lkld2 := 0.0;
            FOR i := 1 TO endptrn DO
            BEGIN
               sum := 0.0; sumd1 := 0.0; sumd2 := 0.0;
               FOR ai := AA TO VV DO
               BEGIN
                  prod1 := freqdyhf[ai]*p^.prba[i][ai];
                  IF prod1 < minreal THEN
                     prod1 := 0.0;                      (* attention *)
                  FOR aj := AA TO VV DO
                  BEGIN
                     prod2 := prod1*q^.prba[i][aj];
                     IF prod2 < minreal THEN
                        prod2 := 0.0;                   (* attention *)
                     sum   := sum   +prod2*tprobl[ai][aj];
                     sumd1 := sumd1 +prod2*tdiff1[ai][aj];
                     sumd2 := sumd2 +prod2*tdiff2[ai][aj];
                  END;
               END;
               IF sum > minreal THEN
               BEGIN                                   (* attention *)
                  lkl   := lkl   +LN(sum)*weight[i];
                  lkld1 := lkld1 +(sumd1/sum)*weight[i];
                  IF (sum*sum) > minreal THEN
                  lkld2 := lkld2 +((sumd2*sum-sumd1*sumd1) /sum/sum)
                                 *weight[i];
               END
               ELSE
               BEGIN
                  IF debugoptn THEN writeln(output
                  ,' *check branchlength1*',p^.number:4,i:4 );
               END;                                    (* attention *)
            END;
            arcold := arc;
            IF (lkld1 <> lkld2) THEN
               arc  := arc -(lkld1/lkld2)
            ELSE
            BEGIN
               arcold := arc + epsilon*0.1;
               IF debugoptn THEN
                  write(output,' *check branchlength2*');
            END;
            IF arc > maxarc THEN arc := minarc;         (* attention *)
            done := ABS(arcold-arc) < epsilon;
         (* IF debugoptn THEN writeln(output,' ':10,arc:10:5,
            arcold:10:5, -(lkld1/lkld2):10:5); *)
         END;
         smoothed := smoothed AND done;
         IF arc < minarc THEN arc := minarc;
         IF arc > maxarc THEN arc := minarc;            (* attention *)
         p^.lnga := arc;
         q^.lnga := arc;
      END; (* branchlengthA *)

      PROCEDURE printupdate ( VAR tr:tree; p:nodepty;
                              vold:REAL; it:INTEGER );
      BEGIN (* printupdate *)
         IF p = tr.startp^.kinp THEN write (numsm:4)
         ELSE                             write (' ':4);
         IF p^.isop = NIL THEN        (* TIP *)
            write  (output, ' ':12)
         ELSE
         BEGIN
            write  (output, '(':4, p^.isop^.kinp^.number:3);
            write  (output,  p^.isop^.isop^.kinp^.number:3, ' )');
         END;
         write  (output, p^.number:3, ' -',p^.kinp^.number:3);
         IF p^.kinp^.isop = NIL THEN        (* TIP *)
            write  (output, ' ':10)
         ELSE
         BEGIN
            write  (output, ' (', p^.kinp^.isop^.kinp^.number:3);
            write  (output, p^.kinp^.isop^.isop^.kinp^.number:3,' )');
         END;
         write  (output, ' ':2, p^.lnga:10:5);
         write  (output, ' (',p^.lnga-vold:10:5,' )');
         writeln(output, it:4);
      END;  (* printupdate *)


      PROCEDURE updateA ( VAR tr:tree; p:nodepty; VAR it:INTEGER );
      VAR
         vold : REAL;
      BEGIN (* updateA *)
         vold := p^.lnga;
         IF p^.isop <> NIL       THEN sublklhdA(p);        (* TIP *)
         IF p^.kinp^.isop <> NIL THEN sublklhdA(p^.kinp);  (* TIP *)
         branchlengthA( p, it );
      (* IF debugoptn AND ((numsm = 1) OR (numsm = maxsmooth)) THEN *)
         IF debugoptn AND (numsm = 1) THEN
            printupdate( tr, p, vold, it );
      END;  (* updateA *)

      PROCEDURE smooth2 ( VAR tr:tree; VAR numit:INTEGER );
      VAR
         n, it : INTEGER;

      BEGIN (* smooth *)
         FOR n := 1 TO numsp DO
         BEGIN
            updateA(tr,tr.brnchp[n],it);
            numit := numit +it;
         END;
         FOR n := ibrnch2 downto ibrnch1 DO
         BEGIN
            updateA(tr,tr.brnchp[n],it);
            numit := numit +it;
         END;
      END;  (* smooth *)


      PROCEDURE smooth ( VAR tr:tree; p:nodepty; VAR numit:INTEGER );
      VAR
         n, it : INTEGER;
         vold  : REAL;
         cp    : nodepty;
      BEGIN (* smooth *)
         BEGIN
            vold := p^.lnga;
            IF p^.isop <> NIL       THEN sublklhdA(p);        (* TIP *)
            IF p^.kinp^.isop <> NIL THEN sublklhdA(p^.kinp);  (* TIP *)
            branchlengthA( p, it );
            numit := numit +it;
            IF debugoptn AND ((numsm = 1) OR (numsm = maxsmooth)) THEN
               printupdate( tr, p, vold, it );
         END;
         IF p^.isop <> NIL THEN        (* TIP *)
         BEGIN
         (* smooth ( tr,p^.isop^.kinp,numit );
            smooth ( tr,p^.isop^.isop^.kinp,numit ); *)
            cp := p^.isop;
            FOR n := 1 TO p^.diverg DO
            BEGIN
               smooth ( tr,cp^.kinp,numit );
               cp := cp^.isop;
            END;
         END;
      END;  (* smooth *)

      PROCEDURE printsmooth( VAR tr:tree; numit:INTEGER );
         VAR i : INTEGER;
      BEGIN (* printsmooth *)
         IF NOT debugoptn THEN
         BEGIN
            write  (output,' ',numsm:3);
            write  (output,numit:4);
            FOR i:=1 TO ibrnch2 DO
               write(output,tr.brnchp[i]^.lnga:5:1);
            writeln(output);
         END;
      END; (* printsmooth *)

      PROCEDURE leastsquares( am : dnonoty;
                          VAR ym : rnodety );
         VAR
         i, j, k : INTEGER;
         pivot, element : REAL;
         im : dnonoty;
      BEGIN (* leastsquares *)
         FOR i := 1 TO ibrnch2 DO
         FOR j := 1 TO ibrnch2 DO
         BEGIN
            IF i = j THEN im[i][j] := 1.0
            ELSE          im[i][j] := 0.0;
         END;
         FOR k := 1 TO ibrnch2 DO
         BEGIN
            pivot := am[k][k];
            ym[k] := ym[k]/pivot;
            FOR j := 1 TO ibrnch2 DO
            BEGIN
               am[k][j] := am[k][j]/pivot;
               im[k][j] := im[k][j]/pivot;
            END;
            FOR i := 1 TO ibrnch2 DO
            BEGIN
               IF k <> i THEN
               BEGIN
                  element := am[i][k];
                  ym[i] := ym[i] -element*ym[k];
                  FOR j := 1 TO ibrnch2 DO
                  BEGIN
                     am[i][j] := am[i][j] -element*am[k][j];
                     im[i][j] := im[i][j] -element*im[k][j];
                  END;
               END;
            END;
         END;
      END;  (* leastsquares *)

      PROCEDURE initlength ( VAR tr:tree );
         VAR
         ia, j1, j2, k, na, n1, n2 : INTEGER;
         suma, sumy : REAL;
         des   : spity;
         dfpair,
         ymt   : rpairty;
         atymt : rnodety;
         amt   : dpanoty;
         atmt  : dnopaty;
         atamt : dnonoty;
      BEGIN (* initlength *)
         FOR ia := 1 TO ibrnch2 DO
         BEGIN
            des := tr.brnchp[ia]^.descen;
            na := 0;
            FOR j1 := 1 TO numsp-1 DO
            BEGIN
               FOR j2 := j1+1 TO numsp DO
               BEGIN
                  na := succ(na);
               (* writeln(output,' ',ia:3,j1:3,j2:3,na:3); *)
                  IF des[j1] <> des[j2] THEN amt[na][ia] := 1.0
                  ELSE                       amt[na][ia] := 0.0;
                  IF ia = 1 THEN
                  BEGIN
                     dfpair[na] := 0;
                     FOR k := 1 TO endsite DO
                     BEGIN
                        IF chsequen[j1][k] <> chsequen[j2][k] THEN
                           dfpair[na] := dfpair[na] + 1.0;
                     END;
                  END;
                  IF dfpair[na] > 0.0 THEN
                     ymt[na] := - LN( 1.0 -dfpair[na]/endsite )
                  ELSE
                     ymt[na] := - LN( 1.0 );
               END;
            END;
         END;

         IF debugoptn THEN
         BEGIN
         writeln(output);
         FOR na := 1 TO numpair DO
         BEGIN
            write(output,' ':1,na:3,' ':1);
            FOR ia := 1 TO ibrnch2 DO
               write(output,trunc(amt[na][ia]):3);
            writeln(output,ymt[na]:8:3,trunc(dfpair[na]):5);
         END;
         END;

         FOR ia := 1 TO numpair  DO
         FOR na := 1 TO ibrnch2 DO
         BEGIN
            atmt[na][ia] := amt[ia][na];
         END;
         FOR n1 := 1 TO ibrnch2 DO
         FOR n2 := 1 TO ibrnch2 DO
         BEGIN
            suma := 0.0;
            FOR ia := 1 TO numpair DO
               suma := suma +atmt[n1][ia]*amt[ia][n2];
            atamt[n1][n2] := suma;
         END;
         FOR n1 := 1 TO ibrnch2 DO
         BEGIN
            sumy := 0.0;
            FOR ia := 1 TO numpair DO
               sumy := sumy +atmt[n1][ia]*ymt[ia];
            atymt[n1] := sumy;
         END;

         IF debugoptn THEN
         BEGIN
         writeln(output);
         FOR n1 := 1 TO ibrnch2 DO
         BEGIN
            write(output,' ':1,n1:3,' ':1);
            FOR n2 := 1 TO ibrnch2 DO
               write(output,trunc(atamt[n1][n2]):3);
            writeln(output,atymt[n1]:8:3);
         END;
         END;

         leastsquares( atamt, atymt );
         FOR na := 1 TO ibrnch2 DO atymt[na] := 100.0*atymt[na];

         IF NOT debugoptn THEN
         IF writeoptn THEN
         BEGIN
         writeln(output);
         write  (output,' arc':5);
         write  (output,'  it':4);
         FOR na := 1 TO ibrnch2 DO write(output,na:3,' ':2);
         writeln(output);
         write  (output,'0':4);
         write  (output,'0':4);
         FOR na := 1 TO ibrnch2 DO write(output,atymt[na]:5:1);
         writeln(output);
         END;

         FOR na := 1 TO ibrnch2 DO
         BEGIN
            tr.brnchp[na]^.lnga       := atymt[na];
            tr.brnchp[na]^.kinp^.lnga := atymt[na];
         END;

      END;  (* initlength *)

      PROCEDURE judgconverg ( oldarcs, newarcs : lengthty;
                                     VAR cnvrg : BOOLEAN );
      (* judgment of convergence *)
         VAR
         i : INTEGER;
         same : BOOLEAN;
      BEGIN (* judgconverg *)
         same := TRUE;
         FOR i := 1 TO ibrnch2 DO
         BEGIN
            IF ABS(newarcs[i]-oldarcs[i]) > epsilon THEN
               same := FALSE;
         END;
         cnvrg := same;
      END;  (* judgconverg *)

      PROCEDURE variance ( VAR tr:tree );
         VAR
         k, m : INTEGER;
         xa1, xa2 : aryamity;
         ai,  aj  : amity;
         arcm, alkl, lkl, lkl2,
         sarc, sarc2, sblklhd,
         ld1, ld2, prod, vlkl : REAL;
         varc, varc2, blklhdm : lengthty;
         tprobm, tdifm1, tdifm2 : daryamity;
      BEGIN (* variance *)
         alkl := tr.lklhd / endsite;
         vlkl := 0.0;
         FOR k := 1 TO endptrn DO
            vlkl := vlkl +SQR(lklhdtrpt[notree][k]-alkl)*weight[k];
         tr.vrilkl  := vlkl;

         FOR m := 1 TO ibrnch2 DO
         BEGIN
            arcm := tr.brnchp[m]^.lnga;
            tdiffmtrx ( arcm, tprobm, tdifm1, tdifm2 );
            sarc  := 0.0;
            sarc2 := 0.0;
            sblklhd := 0.0;
            FOR k := 1 TO endptrn DO
            BEGIN
               xa1 := tr.brnchp[m]^.prba[k];
               xa2 := tr.brnchp[m]^.kinp^.prba[k];
               lkl2 := 0.0;
               ld1  := 0.0;
               ld2  := 0.0;
               FOR ai := AA TO VV DO
               BEGIN
                  prod := freqdyhf[ai]*xa1[ai];
                  IF prod < minreal THEN prod := 0.0;   (* attention *)
                  FOR aj := AA TO VV DO
                  BEGIN
                     IF xa2[aj] < minreal THEN
                        xa2[aj] := 0.0;                 (* attention *)
                     lkl2 := lkl2 +prod*tprobm[ai][aj]*xa2[aj];
                     ld1  := ld1  +prod*tdifm1[ai][aj]*xa2[aj];
                     ld2  := ld2  +prod*tdifm2[ai][aj]*xa2[aj];
                  END;
               END;
               lkl := EXP(lklhdtrpt[notree][k]);
               IF (lkl*lkl) > minreal THEN              (* attention *)
               sarc  := sarc
                  +((ld2*lkl -ld1*ld1) /lkl/lkl )*weight[k];
               IF (lkl2*lkl2) > minreal THEN            (* attention *)
               sarc2 := sarc2
                  +((ld2*lkl2-ld1*ld1) /lkl2/lkl2 )*weight[k];
               IF lkl2 > minreal THEN                   (* attention *)
               sblklhd := sblklhd +LN(lkl2) *weight[k];
            END;
             varc[m]  := sarc;
             varc2[m] := sarc2;
             blklhdm[m] := sblklhd;

(*          IF debugoptn THEN BEGIN
               write  (output,m:4);
               FOR m := 1 TO ibrnch2 DO write(output,sarc[m]:9:6);
               writeln(output);
            END; *)

         END;
         FOR m := 1 TO ibrnch2 DO varc[m]  := 1.0/varc[m];
         FOR m := 1 TO ibrnch2 DO varc2[m] := 1.0/varc2[m];
         tr.vrilnga := varc;
         tr.vrilnga2 := varc2;
         tr.blklhd   := blklhdm;
      END; (* variance *)

     PROCEDURE manuvaluate ( VAR tr:tree );
         VAR
         n, it : INTEGER;
         newarcs, oldarcs : lengthty;
         cnvrg : BOOLEAN;
      BEGIN (* manuvaluate *)
         IF writeoptn THEN writeln(output);
         IF debugoptn THEN writeln(output,' MANUALVALUATE');
     (*  IF debugoptn THEN printcurtree ( tr );  *)
         IF debugoptn THEN writeln(output);
         initdescen( tr );
         initbranch (tr.startp);
         initbranch (tr.startp^.kinp);
     (*  IF debugoptn THEN printcurtree ( tr );  *)
         initlength ( tr );
         IF debugoptn THEN printcurtree ( tr );

         IF debugoptn THEN writeln(output,' SMOOTH');
         cnvrg := FALSE;
         FOR n := 1 TO ibrnch2 DO
            newarcs[n] := tr.brnchp[n]^.lnga;
         numnw := 0;
         numsm := 0;
         REPEAT
            numsm := SUCC(numsm);
            it := 0;
            smooth ( tr,tr.startp^.kinp,it );
            oldarcs := newarcs;
            FOR n := 1 TO ibrnch2 DO
               newarcs[n] := tr.brnchp[n]^.lnga;
            judgconverg(oldarcs,newarcs,cnvrg);
            numnw := numnw +it;
            IF writeoptn THEN printsmooth( tr,it );
         UNTIL ( (numsm >= maxsmooth) OR cnvrg );
         tr.cnvrgnc := cnvrg;

      (* IF debugoptn THEN printcurtree ( tr ); *)
         evaluateA ( tr );
         variance ( tr );
         IF usertree THEN
         BEGIN
            lklhdtree[notree] := tr.lklhd;
            aictree  [notree] := tr.aic;
            paratree [notree] := ibrnch2;
            IF notree = 1 THEN
            BEGIN
               maxltr := 1;
               maxlkl := tr.lklhd;
               minatr := 1;
               minaic := tr.aic;
            END
            ELSE
            BEGIN
               IF tr.lklhd > maxlkl THEN
               BEGIN
                  maxltr := notree;
                  maxlkl := tr.lklhd;
               END;
               IF tr.aic < minaic THEN
               BEGIN
                  minatr := notree;
                  minaic := tr.aic;
               END;
            END;
         END;
      END; (* manuvaluate *)

     PROCEDURE autovaluate ( VAR tr:tree );
         VAR
         n, it : INTEGER;
         newarcs, oldarcs : lengthty;
         cnvrg : BOOLEAN;
      BEGIN (* autovaluate *)
         IF writeoptn THEN writeln(output);
         IF debugoptn THEN writeln(output,' AUTOVALUATE');
      (*
         IF debugoptn THEN printcurtree( tr );
         IF debugoptn THEN writeln(output);
         initdescen( tr );
         initbranch (tr.startp);
         initbranch (tr.startp^.kinp);
         IF debugoptn THEN printcurtree ( tr );
         initlength ( tr );
         IF debugoptn THEN printcurtree ( tr );
      *)
         it := 0;
         FOR n := 1 TO 3 DO
         BEGIN
            sublklhdA(tr.brnchp[ibrnch2]);
            sublklhdA(tr.brnchp[ibrnch2]^.kinp);
            branchlengthA( tr.brnchp[ibrnch2], it );
         END;
         IF debugoptn THEN printcurtree ( tr );

         IF debugoptn THEN writeln(output,' SMOOTH');
         cnvrg := FALSE;
         FOR n := 1 TO ibrnch2 DO
            newarcs[n] := tr.brnchp[n]^.lnga;
         numnw := 0;
         numsm := 0;
         REPEAT
            numsm := SUCC(numsm);
            it := 0;
            smooth ( tr,tr.startp^.kinp,it );
            oldarcs := newarcs;
            FOR n := 1 TO ibrnch2 DO
               newarcs[n] := tr.brnchp[n]^.lnga;
            judgconverg(oldarcs,newarcs,cnvrg);
            numnw := numnw +it;
            IF writeoptn THEN printsmooth( tr,it );
         UNTIL ( (numsm >= maxsmooth) OR cnvrg );
         tr.cnvrgnc := cnvrg;

      (* IF debugoptn THEN printcurtree ( tr ); *)
         evaluateA ( tr );
         variance ( tr );
         IF NOT usertree THEN
         BEGIN
            lklhdtree[notree] := tr.lklhd;
            aictree  [notree] := tr.aic;
            paratree [notree] := ibrnch2;
            IF notree = 1 THEN
            BEGIN
               maxltr := 1;
               maxlkl := tr.lklhd;
               minatr := 1;
               minaic := tr.aic;
            END
            ELSE
            BEGIN
               IF tr.lklhd > maxlkl THEN
               BEGIN
                  maxltr := notree;
                  maxlkl := tr.lklhd;
               END;
               IF tr.aic < minaic THEN
               BEGIN
                  minatr := notree;
                  minaic := tr.aic;
               END;
            END;
         END;
      END; (* autovaluate *)

     (**************************************)
     (***    OUTPUT OF TREE TOPOLOGY     ***)
     (**************************************)

      PROCEDURE prbranch
      ( up:nodepty; depth, m ,maxm:INTEGER;
        VAR umbrella:umbty; VAR length:lspty );
      CONST
         maxover = 50;
         maxleng = 30;
      VAR
         i, j, n, d, maxn : INTEGER;
         cp   : nodepty;
         over : BOOLEAN;
      BEGIN (* prbranch *)
         over := false;
         d  := depth +1;
         IF trunc(up^.lnga*prprtn) >= maxover THEN
         BEGIN
            over := true;
            length[d] := maxleng;
         END
         ELSE
            length[d] := trunc(up^.lnga*prprtn) +3;   (* +4 *)
         IF up^.isop = NIL THEN        (* TIP *)
         BEGIN
            IF m = 1 THEN umbrella[d-1] := true;
            FOR j := 0 TO d-1 DO
            BEGIN
               IF umbrella[j] THEN write(output,':':length[j])
               ELSE                write(output,' ':length[j]);
            END;
            IF m = maxm THEN umbrella[d-1] := false;
            IF over THEN
               FOR i := 1 TO length[d]-3 DO write(output,'+')
            ELSE
               FOR i := 1 TO length[d]-3 DO write(output,'-');
            IF up^.number < 10 THEN
               write(output,'--',up^.number:1,'.')
            ELSE IF up^.number < 100 THEN
               write(output,'-',up^.number:2,'.')
            ELSE
               write(output,up^.number:3,'.');
            FOR i := 1 TO maxname DO write(output, up^.namesp[i]);
         (* write(output,d:36); *)
            writeln(output);
         END
         ELSE
         BEGIN
            cp := up^.isop;
            maxn := up^.diverg;
            FOR n := 1 TO maxn DO
            BEGIN
               prbranch ( cp^.kinp,d,n,maxn,umbrella,length );
               cp := cp^.isop;
               IF n = (maxn DIV 2) THEN
               BEGIN
                  IF m = 1 THEN umbrella[d-1] := true;
                  IF n = 1 THEN umbrella[d-1] := true;
                  FOR j := 0 TO d-1 DO
                  BEGIN
                     IF umbrella[j] THEN write(output,':':length[j])
                     ELSE                write(output,' ':length[j]);
                  END;
                  IF n = maxn THEN umbrella[d-1] := false;
                  IF m = maxm THEN umbrella[d-1] := false;
                  IF over THEN
                     FOR i := 1 TO length[d]-3 DO write(output,'+')
                  ELSE
                     FOR i := 1 TO length[d]-3 DO write(output,'-');
                  IF up^.number < 10 THEN
                     write(output,'--',up^.number:1 (* ,':' *) )
                  ELSE IF up^.number < 100 THEN
                     write(output,'-',up^.number:2 (* ,':' *) )
                  ELSE
                     write(output,up^.number:3 (* ,':' *) );
               (* write(output,d:50); *)
                  writeln(output);
               END
               ELSE IF n <> maxn THEN
               BEGIN
                  FOR j := 0 TO d DO
                  BEGIN
                     IF umbrella[j] THEN write(output,':':length[j])
                     ELSE                write(output,' ':length[j]);
                  END;
                  writeln(output);
               END;
            END;
         END;
      END;  (* prbranch *)

      PROCEDURE prtopology ( VAR tr:tree );
      VAR
         n, maxn,
         depth    : INTEGER;
         cp       : nodepty;
         umbrella : umbty;
         length   : lspty;
      BEGIN (* prtopology *)
         FOR n := 0 TO maxsp DO
         BEGIN
            umbrella[n] := false;
            IF n = 0 THEN length[n] := 3
            ELSE          length[n] := 6;
         END;
         cp := tr.brnchp[0];
         maxn := cp^.diverg+1;
         writeln(output);
         FOR n := 1 TO maxn DO
         BEGIN
            depth := 0;
            prbranch ( cp^.kinp,depth,n,maxn,umbrella,length );
            cp := cp^.isop;
            IF n = (maxn DIV 2) THEN
               writeln(output,'0:':length[0])
            ELSE IF n <> maxn THEN
               writeln(output,':':length[0]);
         END;
      END;  (* prtopology *)

     (*********************************)
     (***  OUTPUT OF TREE TOPOLOGY  ***)
     (*********************************)

      PROCEDURE charasubtoporogy( up:nodepty; VAR nsp,nch:INTEGER );
      VAR
         n, maxn : INTEGER;
         cp   : nodepty;
      BEGIN (* charasubtoporogy *)
         IF up^.isop = NIL THEN        (* TIP *)
         BEGIN
            FOR n := 1 TO maxname DO
               IF up^.namesp[n] <> ' ' THEN
               BEGIN
                  write(output,up^.namesp[n]);
                  nch := nch +1;
               END;
            nsp := nsp +1;
         END
         ELSE
         BEGIN
            cp := up^.isop;
            maxn := up^.diverg;
            write  (output, '(':1);
            nch := nch +1;
            FOR n := 1 TO maxn DO
            BEGIN
               charasubtoporogy ( cp^.kinp,nsp,nch );
               cp := cp^.isop;
               IF n <> maxn THEN
               BEGIN
                  write(output,',':1);
                  nch := nch +1;
                  IF nch > maxline-10 THEN
                  BEGIN
                     writeln(output);
                     write  (output,' ':3);
                     nch := 3;
                  END;
               END;
            END;
            write  (output, ')':1);
            nch := nch +1;
         END;
      END;  (* charasubtoporogy *)

      PROCEDURE charatopology ( VAR tr:tree );
      VAR
         n, maxn,
         nsp, nch : INTEGER;
         cp       : nodepty;
      BEGIN (* charatopology *)
         nsp := 0;
         nch := 3;
         cp := tr.brnchp[0];
         maxn := cp^.diverg+1;
         writeln(output);
         write  (output, ' ( ':3);
         FOR n := 1 TO maxn DO
         BEGIN
            charasubtoporogy ( cp^.kinp,nsp,nch );
            cp := cp^.isop;
            IF n <> maxn THEN
            BEGIN
               write(output,', ':2);
               nch := nch +2;
               IF nch > maxline-15 THEN
               BEGIN
                  writeln(output);
                  write  (output,' ':3);
                  nch := 3;
               END;
            END;
         END;
         writeln(output, ' )');
      END;  (* charatopology *)

     (**************************************)
     (***   OUTPUT OF THE FINAL RESULT   ***)
     (**************************************)

      PROCEDURE printbranch( VAR tr:tree );
         VAR
         i, j : INTEGER;
         p, q : nodepty;
      BEGIN (* printbranch *)
         FOR i := 1 TO ibrnch2 DO
         BEGIN
            p := tr.brnchp[i];
            q := p^.kinp;
            write(output,' ':5);
            IF p^.isop = NIL THEN        (* TIP *)
               FOR j := 1 TO maxname DO write(output, p^.namesp[j])
            ELSE
               FOR j := 1 TO maxname DO write(output, ' ':1);
            write(output, p^.number:5);
            BEGIN
               IF p^.lnga >= maxarc THEN
                  write(output, ' infinity':12)
               ELSE IF p^.lnga <= minarc THEN
                  write(output, ' zero    ':12)
               ELSE
                  write(output, p^.lnga:12:5);
               write(output, ' (':3,
                     SQRT(ABS(tr.vrilnga[i])):9:5,' )':2);
               IF debugoptn THEN
               BEGIN
                  write(output, SQRT(ABS(tr.vrilnga2[i])):12:5);
                  write(output, tr.blklhd[i]:13:5);
               END;
               writeln(output);
            END;
         END;

      END;  (* printbranch *)

      PROCEDURE summarize ( VAR tr:tree );
      BEGIN (* summarize *)
         writeln(output);
      (* printline( 46 ); *)
         IF usertree THEN
            write(output, ' No.':4,notree:3,' ':8)
         ELSE
            write(output, ' No.':4,stage:3,' -':2,notree:3,' ':3);
         write  (output, 'number':7,'Length':9,'S.E.':11);
         IF tr.cnvrgnc THEN
         (* write  (output, 'convergence    ':21) *)
         ELSE
            write  (output, 'non convergence':21);
         IF writeoptn THEN
            write  (output, ':':2,numsm:3,',':2,numnw:4);
         writeln(output);
         printline( 46 );
         printbranch( tr );
         printline( 46 );
         write  (output, 'ln L :':8);
         write  (output, tr.lklhd:10:3);
         write  (output, '(':2,SQRT(ABS(tr.vrilkl)):8:3,' )':2);
         writeln(output, 'AIC :':7,tr.aic:9:3);
         printline( 46 );
(*       writeln(output); *)
      END; (* summarize *)

      PROCEDURE cleartree( VAR tr:tree );
      VAR
         i, n       : INTEGER;
         cp, sp, dp : nodepty;
      BEGIN (* cleartree *)
         FOR i := 1 TO ibrnch2 DO
         BEGIN
            tr.brnchp[i]^.kinp^.kinp := NIL;
            tr.brnchp[i]^.kinp       := NIL;
         END;
         FOR i := ibrnch1 TO ibrnch2+1 DO
         BEGIN
            IF i = ibrnch2+1 THEN n := 0 ELSE n := i;
            sp := tr.brnchp[n];
            cp := sp^.isop;
            sp^.isop := NIL;
            WHILE cp <> sp DO
            BEGIN
               dp := cp;
               cp := cp^.isop;
               dp^.isop := freenode;
               freenode := dp;
               dp := NIL;
            (* dp^.isop := NIL;
               DISPOSE( dp ); *)
            END;
            sp := NIL;
            tr.brnchp[n] := NIL;
            cp^.isop := freenode;
            freenode := cp;
            cp := NIL;
         (* DISPOSE( cp ); *)
         END;
      END;  (* cleartree *)

      PROCEDURE bootstrap;
      VAR
         ib, jb, nb,
         bsite, maxboottree : INTEGER;
         inseed : INTEGER;
         maxlklboot  : REAL;
         seed : longintty;
      BEGIN (* bootstrap *)
         inseed := 12345;
         FOR ib := 0 TO 5 DO seed[ib] := 0;
         ib := 0;
         REPEAT
           seed[ib] := inseed MOD 64;
           inseed   := inseed DIV 64;
           ib       := ib + 1;
         UNTIL inseed = 0;

         FOR jb := 1 TO numtrees DO boottree[jb] := 0.0;
         FOR nb := 1 TO maxboot DO
         BEGIN
            FOR jb := 1 TO numtrees DO
            BEGIN
               lklboottree[jb] := 0.0;
            END;
           FOR ib := 1 TO endsite DO
           BEGIN
              bsite := weightw[ TRUNC(randum(seed)*endsite+1) ];
              FOR jb := 1 TO numtrees DO
              BEGIN
                 lklboottree[jb] := lklboottree[jb]
                                   +lklhdtrpt[jb][bsite];
              END;
           END;

           maxlklboot  := lklboottree[1];
           maxboottree := 1;
           IF debugoptn THEN write(output,nb:5);
           FOR jb := 1 TO numtrees DO
           BEGIN
              IF lklboottree[jb] > maxlklboot THEN
              BEGIN
                 maxlklboot  := lklboottree[jb];
                 maxboottree := jb;
              END;
              IF debugoptn THEN write(output,lklboottree[jb]:8:1);
           END;
           IF debugoptn THEN writeln(output,maxboottree:4);
           boottree[maxboottree] := boottree[maxboottree] +1.0;
         END;
         IF maxboot = 0 THEN
         FOR jb := 1 TO numtrees DO
         BEGIN
            boottree[jb] := boottree[jb] / maxboot;
         END;
      END;  (* bootstrap *)

      PROCEDURE printlklhd;
         VAR
         i, j, cul : INTEGER;
         suml,            (* difference of ln lklhd *)
         suml2,           (* absolute value of difference *)
         suma,            (* difference of AIC *)
         suma2,           (* absolute value of AIC *)
         lklsite,         (* likelihood of site *)
         aicsite,         (* AIC of site *)
         sdl,             (* standard error ln lklhd *)
         sda  : REAL;     (* standard error of AIC *)
      BEGIN (* printlklhd *)
         cul := 71;
         writeln(output);
         IF usertree THEN
            writeln(output, ' ':10,numtrees:4, ' user trees')
         ELSE
            writeln(output, ' NO.':10,stage:3,'stage':6,
            notree:5, 'trees':6);
         IF (numtrees > 0) AND (endsite > 1) THEN
         BEGIN
            printline ( cul );
            writeln(output,' Tree':6,'ln L':6,'Diff ln L':12,
               '#Para':12,'AIC':6,'Diff AIC':12,'Boot P':17);
            printline ( cul );
            FOR i := 1 TO numtrees DO
            BEGIN
               write(output, i:4, lklhdtree[i]:9:2);
               IF maxltr = i THEN
                  write(output, ' <-- best',' ':9)
               ELSE
               BEGIN
                  suml  := 0.0;
                  suml2 := 0.0;
                  FOR j := 1 TO endptrn DO
                  BEGIN
                     lklsite := lklhdtrpt[maxltr][j]-lklhdtrpt[i][j];
                     suml  := suml  +weight[j] *lklsite;
                     suml2 := suml2 +weight[j] *SQR(lklsite);
                  END;
                  sdl := SQRT( (endsite/(endsite-1))
                             *(suml2-SQR(suml/endsite)) );
                  write(output, lklhdtree[i]-maxlkl:9:2);
                  write(output, '+-':3,sdl:6:2);
               END;
               write(output, paratree[i]:4);
               write(output, aictree[i]:9:2);
               IF minatr = i THEN
                  write(output, ' <-- best',' ':9)
               ELSE
               BEGIN
                  suma  := 0.0;
                  suma2 := 0.0;
                  FOR j := 1 TO endptrn DO
                  BEGIN
                     aicsite := -2.0
                        *(lklhdtrpt[maxltr][j]-lklhdtrpt[i][j]);
                     suma  := suma  +weight[j] *aicsite;
                     suma2 := suma2 +weight[j] *SQR(aicsite);
                  END;
                  sda := SQRT( (endsite/(endsite-1))
                             *(suma2-SQR(suma/endsite)) );
                  write(output, aictree[i]-minaic:9:2);
                  write(output, '+-':3,sda:6:2);
               END;
               IF usertree THEN write(output, boottree[i]:9:4);
               writeln(output);
            END;
            printline ( cul );
         END;
      END; (* printlklhd *)

      PROCEDURE outlklhd;
      VAR i, j, k, nt : INTEGER;
      BEGIN (* outlklhd *)
         writeln(lklfile, numsp:5, endsite:5 );
         i := 0;
         FOR nt := 1 TO numtrees DO
         BEGIN
            writeln(lklfile, nt:5 );
            FOR j := 1 TO endptrn DO
            BEGIN
               FOR k := 1 TO weight[j] DO
               BEGIN
                  write(lklfile, lklhdtrpt[nt][j]:12);
                  i := i +1;
                  IF ( i = 6 ) THEN
                  BEGIN
                     writeln(lklfile);
                     i := 0;
                  END;
               END;
            END;
            IF ( i <> 0 ) THEN writeln(lklfile);
            i := 0;
         END;
      END; (* outlklhd *)

   PROCEDURE manutree;  VAR ntree : INTEGER;
   BEGIN (* manutree *)
      (*PAGE(output);*)
      readln (seqfile, numtrees);
      writeln(output);
      writeln(output,' ',numtrees:5,' user trees');
      writeln(output);
      FOR ntree := 1 TO numtrees DO
      BEGIN  notree := ntree;
         treestructure( ctree );
         manuvaluate ( ctree );
         prtopology ( ctree );
         summarize ( ctree );
         cleartree ( ctree );
         (*PAGE(output);*)
      END;
      IF bootsoptn THEN bootstrap;
      printlklhd;
      IF putlkoptn THEN outlklhd;
   END; (* manutree *)

    PROCEDURE newbranch ( nbranch:INTEGER; VAR np1,np2:nodepty );
    VAR
       j : INTEGER;
      np : nodepty;
    BEGIN (* newbranch *)
       FOR j := 1 TO 2 DO
       BEGIN
          IF freenode = NIL THEN
             NEW ( np )
          ELSE
          BEGIN
             np := freenode;
             freenode := np^.isop;
             np^.isop := NIL;
          END;
          clearnode ( np );
          IF j = 1 THEN np1 := np
          ELSE          np2 := np;
          np := NIL;
       END;
       np1^.kinp := np2;
       np2^.kinp := np1;
       ctree.brnchp[nbranch] := np1;
       ctree.brnchp[nbranch]^.number := ibrnch2;
    END;  (* newbranch *)

    PROCEDURE decomposition ( cnode : INTEGER );
    VAR
       i1, i2,
       maxdvg : INTEGER;
       cp1, cp2,
       bp1, bp2,
       lp1, lp2,
       np1, np2 : nodepty;
       maxlkls  : REAL;
    BEGIN (* decomposition *)
      IF ctree.brnchp[cnode]^.diverg > 2 THEN
      BEGIN
      (*
         PAGE(output);
      *)
         stage := stage +1;
         notree := 0;
         ibrnch2 := ibrnch2 +1;
         newbranch( ibrnch2, np1, np2 );
         cp1 := ctree.brnchp[cnode];
         cp2 := cp1;
         bp1 := cp1;
         REPEAT bp1 := bp1^.isop UNTIL bp1^.isop = cp1;
         bp2 := bp1;
         maxlkls := -999999.0;
         maxdvg := ctree.brnchp[cnode]^.diverg;
         IF maxdvg = 3 THEN
         BEGIN
            maxdvg := maxdvg -1;
            cp1 := cp1^.isop;
            cp2 := cp1;
            bp1 := bp1^.isop;
            bp2 := bp1;
         END;
         FOR i1 := 1 TO maxdvg DO
         BEGIN
            FOR i2 := i1+1 TO maxdvg+1 DO
            BEGIN
               IF debugoptn THEN
               BEGIN
                  ibrnch2 := ibrnch2 -1;
                  printcurtree ( ctree );
                  ibrnch2 := ibrnch2 +1;
               END;
               notree := notree +1;
               bp2 := cp2;
               cp2 := cp2^.isop;
               IF debugoptn THEN
                  writeln(output,' AUTO-INS',i1:3,i2:3,
                  'c':4,cp1^.kinp^.number:3,cp2^.kinp^.number:3,
                  'b':4,bp1^.kinp^.number:3,bp2^.kinp^.number:3,
                  'n':4,np1^.number:3,np2^.number:3);
               insertbranch ( ctree,cp1,cp2,bp1,bp2,np1,np2 );
               autovaluate ( ctree );
               prtopology ( ctree );
               charatopology ( ctree );
               summarize ( ctree );
               IF ctree.lklhd > maxlkls THEN
               BEGIN
                  maxlkls := ctree.lklhd;
                  lp1 := bp1;
                  lp2 := bp2;
               END;
               IF debugoptn THEN
                  writeln(output,' AUTO-DEL',i1:3,i2:3,
                  'c':4,cp1^.kinp^.number:3,cp2^.kinp^.number:3,
                  'b':4,bp1^.kinp^.number:3,bp2^.kinp^.number:3,
                  'n':4,np1^.number:3,np2^.number:3);
               deletebranch ( ctree,cnode, cp1,cp2,bp1,bp2,np1,np2 );
            END;
            bp1 := cp1;
            cp1 := bp1^.isop;
            bp2 := bp1;
            cp2 := bp1^.isop;
         END; (* FOR *)
      (*
         PAGE(output);
      *)
         numtrees := notree;
         printlklhd;
         notree := 0;
         cp1 := lp1^.isop;
         cp2 := lp2^.isop;
         bp1 := lp1;
         bp2 := lp2;
         IF debugoptn THEN
            writeln(output,' AUTO-MAX',
            'c':4,cp1^.kinp^.number:3,cp2^.kinp^.number:3,
            'b':4,bp1^.kinp^.number:3,bp2^.kinp^.number:3,
            'n':4,np1^.number:3,np2^.number:3);
         insertbranch ( ctree,cp1,cp2,bp1,bp2,np1,np2 );
         autovaluate ( ctree );
         prtopology ( ctree );
         summarize ( ctree );
      END; (* IF *)
    END;  (* decomposition *)

   PROCEDURE autotree;
   VAR
      i, j,
      dvg,
      cnode : INTEGER;
   BEGIN (* autotree *)
      stage  := 0;
      notree := 1;
      IF semiaoptn THEN treestructure( ctree )
      ELSE              starstructure ( ctree );
      manuvaluate ( ctree );
      prtopology ( ctree );
      summarize ( ctree );

      IF semiaoptn THEN
      BEGIN
        IF firstoptn THEN
        BEGIN
           REPEAT
              decomposition ( 0 );
           UNTIL ctree.brnchp[0]^.diverg < 3;
           REPEAT
              dvg := 2;            (* ctree.brnchp[0]^.diverg *)
              cnode := 0;
              FOR i := ibrnch1 TO ibrnch2 DO
              BEGIN
                 IF ctree.brnchp[i]^.diverg > dvg THEN
                 BEGIN
                    dvg := ctree.brnchp[i]^.diverg;
                    cnode := i;
                 END;
              END;
              decomposition ( cnode );
           UNTIL dvg < 3;
        END
        ELSE IF lastoptn THEN
        BEGIN
           REPEAT
              dvg := 2;            (* ctree.brnchp[0]^.diverg *)
              cnode := 0;
              FOR i := ibrnch1 TO ibrnch2 DO
              BEGIN
                 IF ctree.brnchp[i]^.diverg > dvg THEN
                 BEGIN
                    dvg := ctree.brnchp[i]^.diverg;
                    cnode := i;
                 END;
              END;
              decomposition ( cnode );
           UNTIL dvg < 3;
           REPEAT
              decomposition ( 0 );
           UNTIL ctree.brnchp[0]^.diverg < 3;
        END
        ELSE
        BEGIN
           REPEAT
              dvg := 2;            (* ctree.brnchp[0]^.diverg *)
              cnode := 0;
              FOR i := ibrnch1 TO ibrnch2+1 DO
              BEGIN
                 IF i = ibrnch2+1 THEN j := 0
                 ELSE                  j := i;
                 IF ctree.brnchp[j]^.diverg > dvg THEN
                 BEGIN
                    dvg := ctree.brnchp[j]^.diverg;
                    cnode := j;
                 END;
              END;
              decomposition ( cnode );
           UNTIL dvg < 3;
        END;
      END
      ELSE (* auto mode *)
      BEGIN
        REPEAT
           decomposition ( 0 );
        UNTIL ctree.brnchp[0]^.diverg < 3;
      END;

   (* printlklhd; *)
   (* IF putlkoptn THEN outlklhd; *)
   END; (* autotree *)

 PROCEDURE mainio;  VAR nexe : INTEGER;
 BEGIN (* mainio *)
    freenode := NIL;
    FOR nexe := 1 TO maxexe DO
    BEGIN  numexe := nexe;
       IF maxexe > 1 THEN
          writeln(output,' * JOB :',numexe:4,' *');
       getinput;
       IF normal AND (numexe = 1) THEN maketransprob;
       IF normal THEN
       IF usertree THEN
          manutree
       ELSE
          autotree;
    (* IF maxexe > 1 THEN PAGE(output); *)
    END;
 END;  (* mainio *)

BEGIN (* PROTML *)

(* ASSIGN(seqfile,seqfname);      If TURBO Pascal, use *)
(* ASSIGN(tpmfile,tpmfname); *)
(* IF putlkoptn THEN
   ASSIGN (lklfile,lklfname); *)

   RESET (seqfile);            (* If SUN Pascal, don't use *)
(* RESET (tpmfile); *)
(* IF putlkoptn THEN
   REWRITE(lklfile); *)

(* RESET (seqfile,seqfname);      If SUN Pascal, use *)
(* RESET (tpmfile,tpmfname); *)
(* IF putlkoptn THEN
   REWRITE(lklfile,lklfname); *)

    mainio;

(* CLOSE (seqfile);               If TURBO Pascal, use *)
(* CLOSE (tpmfile); *)
(* IF putlkoptn THEN
   CLOSE (lklfile); *)

END. (* PROTML *)

