$ define lnk$library sys$library:vaxcrtl
$ write sys$output "Please wait, building VAX/VMS phylip V3.5c"
$
$ cc vaxfix.c
$
$ write sys$output "clique"
$ cc clique.c
$ link clique,vaxfix
$ write sys$output "consense"
$ cc consense.c
$ link consense,vaxfix
$ write sys$output "contml"
$ cc contml.c
$ link contml,vaxfix
$ write sys$output "contrast"
$ cc contrast.c
$ link contrast,vaxfix
$ write sys$output "dnacomp"
$ cc dnacomp.c
$ link dnacomp,vaxfix
$ write sys$output "dnadist"
$ cc dnadist.c
$ link dnadist,vaxfix
$ write sys$output "dnainvar"
$ cc dnainvar.c
$ link dnainvar,vaxfix
$ write sys$output "dnaml"
$ cc dnaml.c
$ cc dnaml2.c
$ link dnaml,dnaml2,vaxfix
$ write sys$output "dnamlk"
$ cc dnamlk.c
$ cc dnamlk2.c
$ link dnamlk,dnamlk2,vaxfix
$ write sys$output "dnamove"
$ cc dnamove.c
$ link dnamove,vaxfix
$ write sys$output "dnapars"
$ cc dnapars.c
$ link dnapars,vaxfix
$ write sys$output "dnapenny"
$ cc dnapenny.c
$ link dnapenny,vaxfix
$ write sys$output "dollop"
$ cc dollop.c
$ link dollop,vaxfix
$ write sys$output "dolmove"
$ cc dolmove.c
$ link dolmove,vaxfix
$ write sys$output "dolpenny"
$ cc dolpenny.c
$ link dolpenny,vaxfix
$ cc drawgraphics.c
$ write sys$output "drawgram"
$ cc drawgram.c
$ write sys$output "drawtree"
$ cc drawtree.c
$ link drawgram,drawgraphics,vaxfix
$ link drawtree,drawgraphics,vaxfix
$ write sys$output "factor"
$ cc factor.c
$ link factor,vaxfix
$ write sys$output "fitch"
$ cc fitch.c
$ link fitch,vaxfix
$ write sys$output "gendist"
$ cc gendist.c
$ link gendist,vaxfix
$ write sys$output "kitsch"
$ cc kitsch.c
$ link kitsch,vaxfix
$ write sys$output "mix"
$ cc mix.c
$ cc mix2.c
$ link mix,mix2,vaxfix
$ write sys$output "move"
$ cc move.c
$ link move,vaxfix
$ write sys$output "neighbor"
$ cc neighbor.c
$ link neighbor,vaxfix
$ write sys$output "penny"
$ cc penny.c
$ link penny,vaxfix
$ write sys$output "protdist"
$ cc protdist.c
$ link protdist,vaxfix
$ write sys$output "protpars"
$ cc protpars.c
$ link protpars,vaxfix
$ write sys$output "restml"
$ cc restml.c
$ cc restml2.c
$ link restml,restml2,vaxfix
$ write sys$output "retree"
$ cc retree.c
$ link retree,vaxfix
$ write sys$output "seqboot"
$ cc seqboot.c
$ link seqboot,vaxfix
