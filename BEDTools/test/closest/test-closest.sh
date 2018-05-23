BT=${BT-../../bin/bedtools}

check()
{
	if diff $1 $2; then
    	echo ok
		return 1
	else
    	echo fail
		return 0
	fi
}

# cat a.bed
# chr1	20	21

# cat b.bed
# chr1	20	21

# cat b-one-bp-closer.bed
# chr1	19	21

###########################################################
# test 1bp apart; checking for off-by-one errors
###########################################################
echo "    closest.t1...\c"
echo \
"chr1	10	20	chr1	20	21	1" > exp
$BT closest -a a.bed -b b.bed -d > obs
check obs exp
rm obs exp

###########################################################
# test reciprocal of t1
###########################################################
echo "    closest.t2...\c"
echo \
"chr1	20	21	chr1	10	20	1" > exp
$BT closest -a b.bed -b a.bed -d > obs
check obs exp
rm obs exp

###########################################################
# test 0bp apart; checking for off-by-one errors
###########################################################
echo "    closest.t3...\c"
echo \
"chr1	10	20	chr1	19	21	0" > exp
$BT closest -a a.bed -b b-one-bp-closer.bed -d > obs
check obs exp
rm obs exp

###########################################################
# test reciprocal of t3
###########################################################
echo "    closest.t4...\c"
echo \
"chr1	19	21	chr1	10	20	0" > exp
$BT closest -a b-one-bp-closer.bed -b a.bed -d > obs
check obs exp
rm obs exp

###########################################################
# test closest without forcing different names ( -N )
###########################################################
echo "    closest.t5...\c"
echo \
"chr1	10	20	break1	chr1	40	50	break1	21
chr1	55	58	break2	chr1	60	70	break2	3" > exp
$BT closest -a a.names.bed -b b.names.bed -d > obs
check obs exp
rm obs exp

###########################################################
# test closest with forcing different names ( -N )
###########################################################
echo "    closest.t6...\c"
echo \
"chr1	10	20	break1	chr1	60	70	break2	41
chr1	55	58	break2	chr1	40	50	break1	6" > exp
$BT closest -a a.names.bed -b b.names.bed -d -N > obs
check obs exp
rm obs exp

###########################################################
# test closest forcing -s yet no matching strands on chrom
###########################################################
echo "    closest.t7...\c"
echo \
"chr1	100	200	a	10	+	.	-1	-1	.	-1	." > exp
$BT closest -a strand-test-a.bed -b strand-test-b.bed -s > obs
check obs exp
rm obs exp

###########################################################
# test closest forcing -S with only an opp strands on chrom
###########################################################
echo "    closest.t8...\c"
echo \
"chr1	100	200	a	10	+	chr1	90	120	b	1	-" > exp
$BT closest -a strand-test-a.bed -b strand-test-b.bed -S > obs
check obs exp
rm obs exp


###########################################################
# test reproting of all overlapping features
###########################################################
echo "    closest.t9...\c"
echo \
"chr1	100	101	chr1	100	101
chr1	200	201	chr1	150	201
chr1	200	201	chr1	175	375
chr1	300	301	chr1	175	375
chr1	100000	100010	chr1	175	375
chr1	100020	100040	chr1	175	375
chr2	1	10	.	-1	-1
chr2	20	30	.	-1	-1" > exp
$BT closest -a close-a.bed -b close-b.bed > obs
check obs exp
rm obs exp

###########################################################
# test reproting of first overlapping feature
###########################################################
echo "    closest.t10...\c"
echo \
"chr1	100	101	chr1	100	101
chr1	200	201	chr1	150	201
chr1	300	301	chr1	175	375
chr1	100000	100010	chr1	175	375
chr1	100020	100040	chr1	175	375
chr2	1	10	.	-1	-1
chr2	20	30	.	-1	-1" > exp
$BT closest -a close-a.bed -b close-b.bed -t first > obs
check obs exp
rm obs exp

###########################################################
# test reproting of last overlapping feature
###########################################################
echo "    closest.t11...\c"
echo \
"chr1	100	101	chr1	100	101
chr1	200	201	chr1	175	375
chr1	300	301	chr1	175	375
chr1	100000	100010	chr1	175	375
chr1	100020	100040	chr1	175	375
chr2	1	10	.	-1	-1
chr2	20	30	.	-1	-1" > exp
$BT closest -a close-a.bed -b close-b.bed -t last > obs
check obs exp
rm obs exp



###########################################################
#
# TEST MULTIPLE DATABASES
#
###########################################################


###########################################################
# test 3 dbs, -each mode, which is the default
###########################################################
echo "    closest.t13...\c"
echo \
"chr1	80	100	q1	1	+	1	chr1	20	60	d1.2	2	-
chr1	80	100	q1	1	+	2	chr1	120	170	db2.2	2	-
chr1	80	100	q1	1	+	3	chr1	70	90	d3.1	3	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed mdb3.bed > obs
check obs exp
rm obs exp


###########################################################
# test 3 dbs with -names option
###########################################################
echo "    closest.t14...\c"
echo \
"chr1	80	100	q1	1	+	a	chr1	20	60	d1.2	2	-
chr1	80	100	q1	1	+	b	chr1	120	170	db2.2	2	-
chr1	80	100	q1	1	+	c	chr1	70	90	d3.1	3	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed mdb3.bed -names a b c > obs
check obs exp
rm obs exp


###########################################################
# test 3 dbs with -filenames option
###########################################################
echo "    closest.t15...\c"
echo \
"chr1	80	100	q1	1	+	mdb1.bed	chr1	20	60	d1.2	2	-
chr1	80	100	q1	1	+	mdb2.bed	chr1	120	170	db2.2	2	-
chr1	80	100	q1	1	+	mdb3.bed	chr1	70	90	d3.1	3	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed mdb3.bed -filenames > obs
check obs exp
rm obs exp

###########################################################
# test 3 dbs, -all mode
###########################################################
echo "    closest.t16...\c"
echo \
"chr1	80	100	q1	1	+	3	chr1	70	90	d3.1	3	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed mdb3.bed -mdb all > obs
check obs exp
rm obs exp


###########################################################
# test 2 dbs, tie mode = all
###########################################################
echo "    closest.t17...\c"
echo \
"chr1	80	100	q1	1	+	1	chr1	20	60	d1.2	2	-
chr1	80	100	q1	1	+	2	chr1	120	170	db2.2	2	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed -t all > obs
check obs exp
rm obs exp

###########################################################
# test 2 dbs, tie mode = first
###########################################################
echo "    closest.t18...\c"
echo \
"chr1	80	100	q1	1	+	1	chr1	20	60	d1.2	2	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed -mdb all -t first > obs
check obs exp
rm obs exp

###########################################################
# test 2 dbs, tie mode = last
###########################################################
echo "    closest.t19...\c"
echo \
"chr1	80	100	q1	1	+	2	chr1	120	170	db2.2	2	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed -mdb all -t last > obs
check obs exp
rm obs exp

###########################################################
# test same strand 
###########################################################
echo "    closest.t20...\c"
echo \
"chr1	80	100	q1	1	+	chr1	5	15	d1.1	1	+" > exp
$BT closest -a mq1.bed -b mdb1.bed -s> obs
check obs exp
rm obs exp

###########################################################
# test diff strand 
###########################################################
echo "    closest.t21...\c"
echo \
"chr1	80	100	q1	1	+	chr1	20	60	d1.2	2	-" > exp
$BT closest -a mq1.bed -b mdb1.bed -S> obs
check obs exp
rm obs exp


###########################################################
# test 2 dbs, tie mode = all, same strand
###########################################################
echo "    closest.t22...\c"
echo \
"chr1	80	100	q1	1	+	1	chr1	5	15	d1.1	1	+" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed -t all -mdb all -s > obs
check obs exp
rm obs exp


###########################################################
# test 2 dbs, tie mode = all, diff strand
###########################################################
echo "    closest.t23...\c"
echo \
"chr1	80	100	q1	1	+	1	chr1	20	60	d1.2	2	-
chr1	80	100	q1	1	+	2	chr1	120	170	db2.2	2	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed -t all -mdb all -S > obs
check obs exp
rm obs exp






###########################################################
#
# TEST -D OPTION
#
###########################################################


###########################################################
# hit on left, forward query, forward hit, ref mode
###########################################################
echo "    closest.t24...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	40	60	d1F.1	10	+	-21" > exp
$BT closest -a d_q1.bed -b d_d1F.bed -D ref > obs
check obs exp
rm obs exp


###########################################################
# hit on left, forward query, forward hit, a mode
###########################################################
echo "    closest.t25...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	40	60	d1F.1	10	+	-21" > exp
$BT closest -a d_q1.bed -b d_d1F.bed -D a > obs
check obs exp
rm obs exp

###########################################################
# hit on left, forward query, forward hit, b mode
###########################################################
echo "    closest.t26...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	40	60	d1F.1	10	+	21" > exp
$BT closest -a d_q1.bed -b d_d1F.bed -D b > obs
check obs exp
rm obs exp

###########################################################
# hit on left, forward query, reverse hit, ref mode
###########################################################
echo "    closest.t27...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	40	60	d1R.1	10	-	-21" > exp
$BT closest -a d_q1.bed -b d_d1R.bed -D ref > obs
check obs exp
rm obs exp


###########################################################
# hit on left, forward query, reverse hit, a mode
###########################################################
echo "    closest.t28...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	40	60	d1R.1	10	-	-21" > exp
$BT closest -a d_q1.bed -b d_d1R.bed -D a > obs
check obs exp
rm obs exp

###########################################################
# hit on left, forward query, reverse hit, b mode
###########################################################
echo "    closest.t29...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	40	60	d1R.1	10	-	-21" > exp
$BT closest -a d_q1.bed -b d_d1R.bed -D b > obs
check obs exp
rm obs exp

###########################################################
# hit on left, reverse query, forward hit, ref mode
###########################################################
echo "    closest.t30...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	40	60	d1F.1	10	+	-21" > exp
$BT closest -a d_q2.bed -b d_d1F.bed -D ref > obs
check obs exp
rm obs exp

###########################################################
# hit on left, reverse query, forward hit, a mode
###########################################################
echo "    closest.t31...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	40	60	d1F.1	10	+	21" > exp
$BT closest -a d_q2.bed -b d_d1F.bed -D a > obs
check obs exp
rm obs exp


###########################################################
# hit on left, reverse query, forward hit, b mode
###########################################################
echo "    closest.t32...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	40	60	d1F.1	10	+	21" > exp
$BT closest -a d_q2.bed -b d_d1F.bed -D b > obs
check obs exp
rm obs exp




###########################################################
# hit on left, reverse query, reverse hit, ref mode
###########################################################
echo "    closest.t33...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	40	60	d1R.1	10	-	-21" > exp
$BT closest -a d_q2.bed -b d_d1R.bed -D ref > obs
check obs exp
rm obs exp

###########################################################
# hit on left, reverse query, reverse hit, a mode
###########################################################
echo "    closest.t34...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	40	60	d1R.1	10	-	21" > exp
$BT closest -a d_q2.bed -b d_d1R.bed -D a > obs
check obs exp
rm obs exp


###########################################################
# hit on left, reverse query, reverse hit, b mode
###########################################################
echo "    closest.t35...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	40	60	d1R.1	10	-	-21" > exp
$BT closest -a d_q2.bed -b d_d1R.bed -D b > obs
check obs exp
rm obs exp



###########################################################
# hit on right, forward query, forward hit, ref mode
###########################################################
echo "    closest.t36...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	140	160	d2F.1	10	+	41" > exp
$BT closest -a d_q1.bed -b d_d2F.bed -D ref > obs
check obs exp
rm obs exp

###########################################################
# hit on right, forward query, forward hit, a mode
###########################################################
echo "    closest.t37...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	140	160	d2F.1	10	+	41" > exp
$BT closest -a d_q1.bed -b d_d2F.bed -D a > obs
check obs exp
rm obs exp

###########################################################
# hit on right, forward query, forward hit, b mode
###########################################################
echo "    closest.t38...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	140	160	d2F.1	10	+	-41" > exp
$BT closest -a d_q1.bed -b d_d2F.bed -D b > obs
check obs exp
rm obs exp


###########################################################
# hit on right, forward query, reverse hit, ref mode
###########################################################
echo "    closest.t39...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	140	160	d2R.1	10	-	41" > exp
$BT closest -a d_q1.bed -b d_d2R.bed -D ref > obs
check obs exp
rm obs exp

###########################################################
# hit on right, forward query, reverse hit, a mode
###########################################################
echo "    closest.t40...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	140	160	d2R.1	10	-	41" > exp
$BT closest -a d_q1.bed -b d_d2R.bed -D a > obs
check obs exp
rm obs exp

###########################################################
# hit on right, forward query, reverse hit, b mode
###########################################################
echo "    closest.t41...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	140	160	d2R.1	10	-	41" > exp
$BT closest -a d_q1.bed -b d_d2R.bed -D b > obs
check obs exp
rm obs exp

###########################################################
# hit on right, reverse query, forward hit, ref mode
###########################################################
echo "    closest.t42...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	140	160	d2F.1	10	+	41" > exp
$BT closest -a d_q2.bed -b d_d2F.bed -D ref > obs
check obs exp
rm obs exp

###########################################################
# hit on right, reverse query, forward hit, a mode
###########################################################
echo "    closest.t43...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	140	160	d2F.1	10	+	-41" > exp
$BT closest -a d_q2.bed -b d_d2F.bed -D a > obs
check obs exp
rm obs exp

###########################################################
# hit on right, reverse query, forward hit, b mode
###########################################################
echo "    closest.t44...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	140	160	d2F.1	10	+	-41" > exp
$BT closest -a d_q2.bed -b d_d2F.bed -D b > obs
check obs exp
rm obs exp

###########################################################
# hit on right, reverse query, reverse hit, ref mode
###########################################################
echo "    closest.t45...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	140	160	d2R.1	10	-	41" > exp
$BT closest -a d_q2.bed -b d_d2R.bed -D ref > obs
check obs exp
rm obs exp

###########################################################
# hit on right, reverse query, reverse hit, a mode
###########################################################
echo "    closest.t46...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	140	160	d2R.1	10	-	-41" > exp
$BT closest -a d_q2.bed -b d_d2R.bed -D a > obs
check obs exp
rm obs exp

###########################################################
# hit on right, reverse query, reverse hit, b mode
###########################################################
echo "    closest.t47...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	140	160	d2R.1	10	-	41" > exp
$BT closest -a d_q2.bed -b d_d2R.bed -D b > obs
check obs exp
rm obs exp

############################################################
# Make sure non-overlapping ties are reported 
############################################################
echo "    closest.t48...\c"
echo \
"chr1	10	20	a1	1	-	chr1	8	9	b1	1	+
chr1	10	20	a1	1	-	chr1	21	22	b2	1	-" > exp
$BT closest -a a2.bed -b b2.bed > obs
check obs exp
rm obs exp

############################################################
# Make sure non-overlapping ties are reported, but with -s
############################################################
echo "    closest.t49...\c"
echo \
"chr1	10	20	a1	1	-	chr1	21	22	b2	1	-" > exp
$BT closest -a a2.bed -b b2.bed -s > obs
check obs exp
rm obs exp

############################################################
# Make sure non-overlapping ties are reported, but with -S
############################################################
echo "    closest.t50...\c"
echo \
"chr1	10	20	a1	1	-	chr1	8	9	b1	1	+" > exp
$BT closest -a a2.bed -b b2.bed -S > obs
check obs exp
rm obs exp


echo "    closest.t51...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	140	160	d2R.1	10	-	41" > exp
$BT closest -a d_q2.bed -b d_d2R.bed -D b > obs
check obs exp
rm obs exp


###########################################################
# see that stranded sweep doesn't prematurely purge 
# records from the cache
###########################################################
echo "    closest.t52...\c"
echo \
"chr1	80	100	c1	20	+	chr1	50	60	d2	20	+
chr1	110	130	c2	20	-	chr1	20	40	d1	20	-" > exp
$BT closest -a strand-test-c.bed -b strand-test-d.bed -s > obs
check obs exp
rm obs exp

###########################################################
# check oppostite stranded sweep 
###########################################################
echo "    closest.t53...\c"
echo \
"chr1	80	100	c1	20	+	chr1	20	40	d1	20	-
chr1	110	130	c2	20	-	chr1	50	60	d2	20	+" > exp
$BT closest -a strand-test-c.bed -b strand-test-d.bed -S > obs
check exp obs
rm obs exp


###########################################################
# check -iu
###########################################################
echo "    closest.t54...\c"
echo \
"chr1	100	120	chr1	200	210	81" > exp
$BT closest -a d.bed -b d_iu.bed -D ref -iu > obs
check exp obs
rm exp obs


##########################################################
# check -id
###########################################################
echo "    closest.t55...\c"
echo \
"chr1	100	120	chr1	10	20	-81" > exp
$BT closest -a d.bed -b d_id.bed -D ref -id > obs
check exp obs
rm exp obs

##########################################################
# check ties, single db
###########################################################
echo "    closest.t56...\c"
echo \
"chr1	10	20	a1	1	-	chr1	8	9	b1	1	+
chr1	10	20	a1	1	-	chr1	21	22	b2	1	-" > exp
$BT closest -a bug157_a.bed -b bug157_b.bed > obs
check exp obs
rm exp obs


##########################################################
# check ties, single db, -iu
###########################################################
echo "    closest.t57...\c"
echo \
"chr1	10	20	a1	1	-	chr1	21	22	b2	1	-	2" > exp
$BT closest -a bug157_a.bed -b bug157_b.bed -D ref -iu > obs
check exp obs
rm exp obs

##########################################################
# check ties, single db, -id
###########################################################
echo "    closest.t58...\c"
echo \
"chr1	10	20	a1	1	-	chr1	8	9	b1	1	+	-2" > exp
$BT closest -a bug157_a.bed -b bug157_b.bed -D ref -id > obs 
check exp obs
rm exp obs

##########################################################
# check -header
###########################################################
echo "    closest.t59...\c"
echo \
"#Header for file a.bed
chr1	10	20	chr1	20	21" > exp
$BT closest -a a.bed -b b.bed -header > obs
check exp obs
rm exp obs

cd sortAndNaming
bash test-sort-and-naming.sh
cd ..

cd kclosest
bash test-kclosest.sh
cd ..

