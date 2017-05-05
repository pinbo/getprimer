#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  getprimer.py
#
#  Copyright 2016 Junli Zhang <zhjl86@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

### version update
# v1.2
# For primer 3 default version 1

# v.1.3
# Use output with version 2 and default settings from primer3web_v4
# added hairpin in primer class

# v.1.4
# added primer pair check with primer3
# added total number of different sequences of each primer (primer.difnum)

# 11/19/2016 v.1.5
# from getprimer.v1.4withpaircheck.py
# now input file only needs to provide the raw sequences

# 11/20/2016 v1.6
# added primer score and primer pair score
# and a few other changes

# 11/21/2016 v1.7
# scored now is (diffnumber*100/number_of_other_sequence + 100) / (dif_position + 1), so the dividing is from 2.
# added command line options

# 11/28/2016
# modified the score calculation
# add filter of primers to remove duplicated primers
# look for overlap primers in all possible primers (leftprimers and rightprimers)

# 12/26/2016
# modify to make primers more specific: 2 mismatches in the last 4 3'-termini bases or 5 mismatches in the last 15 3'-termini bases
# primer output now also includes base pair differences from all other sequences
# use primer3 default (instead of web_v4 settings) to loose the criteria

# 1/4/2017
# add support for windows 7.
# more consideration for primers across user-defined border

# 4/23/2017
# add support for MacOSX 64bit
# removed option "-p" for script path
# loosed the criteria to 1) 1 mismatch in the first 4 3'-termini bases or 4 mismatches in the last 15 3'-termini base

# 4/24/2017
# add an galaxy tool xml file getprimer.xml
# add more primer3 options

# 4/25/2017
# modified the fasta name processing so it can extract only the sequence identifier but not the comments after the first space.

# 4/26/2017
# add blast check options

# 5/3/2017: 1. set word_size of blastn to 7; 2. set potential target of primers in the blast to: less than 2 mismatches in the first 4 3' bases.

### Imported
from subprocess import call
import getopt, sys, os

usage="""
getprimer.py
	-i <sequence.fa>
	-s <product min size>
	-l <product max size>
	-g <seqID1,seqID2>
	-r <range limit: where the differences should be in the primer from 3' end, default 10>
	-o <output file name>
	-v <primer overlap region (such as intron): n1-n2,n3-n4>
	-f <1 or 0, default 1: filter the output>
	-a <primer_pair_compl_any_threshold: default 10>
	-e <primer_pair_compl_end_threshold: default 10>
	-c <primer_pair_score_threshold: default 50>
	-b <blast all produced primers against the genomes: 1 for YES and 0 for NO, default is NO>
	--mintm <primer min Tm, default 58>
	--maxtm <primer max Tm, default 62>
	--minsize <primer min size, default 18>
	--maxsize <primer max size, default 23>
	--maxtmdiff <max Tm difference between left and right primers, default 2>
	--maxhairpin <maximum hairpin score of each primer, default 24>
"""


#########################
# parameter or file names that need to be changed

seqfile = "sequence.fa"
targets = []
primer_pair_score_threshold = 50.0
primer_pair_compl_any_threshold = 10.0 # web_v4 default 45
primer_pair_compl_end_threshold = 10.0 # web_v4 default 35
product_min = 100
product_max = 1000
minTm = 58
maxTm = 62
maxTmdiff = 2
minSize = 18
maxSize = 23
maxhairpin = 35 # web_v4 default is 24. primer3 default is 47
blast = 0 # whether blast to check the primer specificity

getprimer_path = os.path.dirname(os.path.realpath(__file__))
rangelimit = 15 # only find difference in the first a few nt from 3' end of each primer
out = ""
msa = 1 # whether need to do multiple sequence alignment
overlap_region = [] # intron region
filter_flag = 0 # whether to filter the primers to remove primers in the same positions
# read command line options
print "Parsing command line options"

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:p:s:l:g:r:o:m:v:f:a:e:c:b:h", ["help", "mintm=", "maxtm=", "minsize=", "maxsize=", "maxtmdiff=", "maxhairpin="])
except getopt.GetoptError as err:
	# print help information and exit:
	print str(err)  # will print something like "option -a not recognized"
	print usage
	sys.exit(2)
for o, a in opts:
	if o == "-i":
		seqfile = a
	elif o in ("-h", "--help"):
		print usage
		sys.exit()
	elif o in ("-p"):
		getprimer_path = a
	elif o in ("-s"):
		product_min = int(a)
	elif o in ("-l"):
		product_max = int(a)
	elif o in ("-g"):
		targets = a.split(",")
		print targets
	elif o in ("-o"):
		out = a
	elif o in ("-m"):
		msa = int(a)
	elif o in ("-r"):
		rangelimit = int(a)
	elif o in ("-f"):
		filter_flag = int(a)
	elif o in ("-a"):
		primer_pair_compl_any_threshold = int(a)
	elif o in ("-e"):
		primer_pair_compl_end_threshold = int(a)
	elif o in ("-c"):
		primer_pair_score_threshold = int(a)
	elif o in ("-b"):
		blast = int(a)
	elif o in ("-v"):
		regions = a.split(",")
		for rr in regions:
			r1, r2 = rr.split("-")
			overlap_region += range(int(r1), int(r2)+1)
	elif o in ("--mintm"):
		minTm = int(a)
	elif o in ("--mintm"):
		minTm = int(a)
	elif o in ("--maxtm"):
		maxTm = int(a)
	elif o in ("--minsize"):
		minSize = int(a)
	elif o in ("--maxsize"):
		maxSize = int(a)
	elif o in ("--maxtmdiff"):
		maxTmdiff = int(a)
	elif o in ("--maxhairpin"):
		maxhairpin = int(a)
	else:
		assert False, "unhandled option"
print "OPtions done"
if not targets:
	print "Please give the sequence IDs for primer design!"
	sys.exit(1)
groupname = "-".join(targets)
getprimer_path = os.path.expanduser(getprimer_path)
#primer3_parameter_path = getprimer_path + "/primer3web_v4_JZ.txt"

#from sys import platform
if sys.platform.startswith('linux'): # linux
	primer3_path = getprimer_path + "/bin/primer3_core"
	muscle_path = getprimer_path + "/bin/muscle"
elif sys.platform == "win32" or sys.platform == "cygwin": # Windows...
	primer3_path = getprimer_path + "/bin/primer3_core.exe"
	muscle_path = getprimer_path + "/bin/muscle.exe"
elif sys.platform == "win32" or sys.platform == "darwin": # Windows...
        primer3_path = getprimer_path + "/bin/primer3_core_darwin64"
        muscle_path = getprimer_path + "/bin/muscle3.8.31_i86darwin64"
# other variables
if not out:
	out = 'selected_primers_for_' + groupname + ".txt"

#########################
# STEP 0: create alignment file and primer3output file

RawAlignFile = "alignment_raw.fa"
alignmentcmd = muscle_path + " -in " + seqfile + " -out " + RawAlignFile + " -quiet"
print "Alignment command: ", alignmentcmd
if msa:
	call(alignmentcmd, shell=True)

###################
# STEP 1: read alignment fasta file into a dictionary AND format it by removing leading and ending "-"
fasta = {} # dictionary for alignment

with open(RawAlignFile) as file_one:
	for line in file_one:
		line = line.strip()
		if line.startswith(">"):
			sequence_name = line.split()[0].lstrip(">")
		else:
			fasta.setdefault(sequence_name, "")
			fasta[sequence_name] += line.rstrip()

## calculate gap number
gap_left = 0
gap_right = 0

# gap left number
for k, v in fasta.items():
	n = 0
	for i in v:
		if i == "-":
			n += 1
		else:
			break
	if n > gap_left:
		gap_left = n
print "Gap left: ", gap_left

# right gap
for k, v in fasta.items():
	n = 0
	for i in reversed(v):
		if i == "-":
			n += 1
		else:
			break
	if n > gap_right:
		gap_right = n

print "Gap right: ", gap_right


### striping sequences to equal length
for k, v in fasta.items():
	if gap_right == 0:
		fasta[k] = v[gap_left:]
	else:
		fasta[k] = v[gap_left:-gap_right]

###########################
# STEP 2: get variation sites
ids = [] # all other sequence names
for kk in fasta.keys():
	if kk not in targets:
		ids.append(kk)

print "The other groups: ", ids

alignlen = len(fasta[targets[0]])
print "Alignment length: ", alignlen

# get the target ID template base coordinate in the alignment
t2a = {}
ngap = 0 # gaps
mainID = targets[0]
for i in range(alignlen):
	if fasta[mainID][i] == "-":
		ngap += 1
		continue
	t2a[i - ngap] = i

# function to get primer variations
def getvar(pp): # pp is a Primer object
	global fasta, t2a, targets, ids
	diffarray = [] # to record difference with each other seq in each site, for primer pair selection later
	ngap = 0 # gaps
	mainID = targets[0] # the one whose primer3 output will be used
	exclude = 0 # a tag to see whether the primer has difference within target group
	align_left = t2a[pp.start - 1] # pp.start and end is 1-based
	align_right = t2a[pp.end - 1]
	# check whether the primer is in all the targets
	#print "targets ", targets
	#print "pp.seq ", pp.seq
	#print "pp.start ", pp.start
	#print "pp.end ", pp.end
	#print "align_left is ", align_left
	#print "align_right ", align_right
	for j in targets:#check difference within targets
		targetseq = fasta[j][align_left:(align_right + 1)].replace("-","")
		#print "targetseq ", targetseq
		if targetseq != pp.seq:
			exclude = 1
			break
	if exclude == 1:
		return 0
	# if it is good for all targets
	for i in range(align_left, align_right + 1):
		if fasta[mainID][i] == "-":
			ngap += 1
			continue # if the target is a gap, just pass to next circle.
		nl1 = len(fasta[mainID][:i]) - len(fasta[mainID][:i].rstrip("-")) # number of gaps on the left of current site for the target
		nr1 = len(fasta[mainID][(i+1):]) - len(fasta[mainID][(i+1):].lstrip("-")) + 1 # number of gaps on the right of current site for the target
		da = [0] * len(ids) # differ array for each base
		m = 0 # counter of next loop
		for k in ids:
			# below nl1 to nr2 is for handling gaps
			nl2 = len(fasta[k][:i]) - len(fasta[k][:i].rstrip("-")) # number of gaps on the left of current site for the homeolog
			nr2 = len(fasta[k][(i+1):]) - len(fasta[k][(i+1):].lstrip("-")) + 1 # number of gaps on the right of current site for the homeolog
			if i - ngap > pp.length / 2.0: # determine whether to shift left or right
				b1 = fasta[mainID][i:].replace("-","")[nl1] # target non-gap base
				b2 = fasta[k][i:].replace("-","")[nl2] # homeolog non-gap bas
			else:
				b1 = fasta[mainID][:(i+1)].replace("-","")[-nr1] # target non-gap base
				b2 = fasta[k][:(i+1)].replace("-","")[-nr2] # homeolog non-gap bas
			if b1 != b2:
				da[m] = 1 # m sequence has variation from target
			m += 1
		diffarray.append(da)		
	# get differsite4 and differsite15
	pp.nvar =  sum([max(x) for x in zip(*diffarray)])
	diffnum = [sum(x) for x in diffarray] # diff number in each base
	#print "diffnum is ", diffnum
	#print "pp.length is ", pp.length
	if pp.direction == "LEFT_PRIMER":
		pp.difsite = [sum(x) for x in zip(*diffarray[(pp.length - 15):pp.length])]
		pp.difsite4 = [sum(x) for x in zip(*diffarray[(pp.length - 4):pp.length])]
		pp.score = sum([float(diffnum[i]) / len(ids) * 100 / (i + 1) for i in range(pp.length)[::-1]])
		if sum(diffarray[-1]) > 0:
			pp.difthree = "YES"
	else:
		pp.difsite = [sum(x) for x in zip(*diffarray[0:15])]
		pp.difsite4 = [sum(x) for x in zip(*diffarray[0:4])]
		pp.score = sum([float(diffnum[i]) / len(ids) * 100 / (i + 1) for i in range(pp.length)])
		if sum(diffarray[0]) > 0:
			pp.difthree = "YES"
	#print "pp.difsite is ", pp.difsite
	return pp


####################################
# STEP 3: create primer3.input and run it for output

# primer3 inputfile
primer3input = "primer3.input"
p3input = open(primer3input, 'w')

# I only need to use the output from the first target
line0 = "PRIMER_TASK=pick_primer_list"
line1 = "SEQUENCE_ID=" + mainID
line2 = "SEQUENCE_TEMPLATE=" + fasta[mainID].replace("-","") # remove "-" in the alignment seq
line3 = "PRIMER_FIRST_BASE_INDEX=1"
line4 = "PRIMER_PRODUCT_SIZE_RANGE=" + str(product_min) + "-" + str(product_max)
line5 = "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + getprimer_path + "/bin/primer3_config/"
line6 = "PRIMER_NUM_RETURN=10000"
line7 = "PRIMER_MIN_SIZE=" + str(minSize)
line8 = "PRIMER_MAX_SIZE=" + str(maxSize)
line9 = "PRIMER_MIN_TM=" + str(minTm)
line10 = "PRIMER_MAX_TM=" + str(maxTm)
line11 = "PRIMER_PAIR_MAX_DIFF_TM=" + str(maxTmdiff)
line12 = "PRIMER_MAX_HAIRPIN_TH=" + str(maxhairpin)
line13 = "="

p3input.write("\n".join([line0, line1, line2, line4, line5, line6, line7, line8, line9, line10, line11, line12, line13]) + "\n")

p3input.close()

# primer3 output file
primer3output = "primer3.output"
#p3cmd = primer3_path + " -default_version=2 -format_output -p3_settings_file=" + primer3_parameter_path + " -output=" + primer3output + " " + primer3input
p3cmd = primer3_path + " -default_version=2 -format_output -output=" + primer3output + " " + primer3input
print "Primer3 command 1st time: ", p3cmd
call(p3cmd, shell=True)

#####################################
# STEP 4: parse primer3 output

# function to get reverse complement
def ReverseComplement(seq):
	# too lazy to construct the dictionary manually, use a dict comprehension
	seq1 = 'ATCGTAGCatcgtagc'
	seq_dict = { seq1[i]:seq1[i+4] for i in range(16) if i < 4 or 8<=i<12 }
	return "".join([seq_dict[base] for base in reversed(seq)])

# Primer class
class Primers(object):
	"""A primer set designed by Primer3"""
	def __init__(self):
		self.name = ""
		self.direction = ""
		self.start = 0
		self.end = 0
		self.length = 0
		self.tm = 0.0
		self.gc = 0.0
		self.anys = 0.0
		self.three = 0.0
		self.hairpin = 0.0
		self.seq = ""
		self.difsite = [] # number of differences with other sequences in the first 15 bases from 3'
		self.difsite4 = [] # number of differences with other sequences in the first 4 bases from 3'
		self.nvar = 0 # number of variation from other sequences
		self.difthree = "NO" # whether 3' site is different
		self.difall = "NO" # whether one primer is enough to differ all the rest
		self.difnum = 0 # how many sequences can be differentiated within the rangelimit
		self.score = 0.0 # score based on number of different sites and differ position
		self.difsitedict = {} # difsite in each different position (just like 2-dim array): self.difsite would be sum of cols of self.difsitedict
		self.overlap = "NO" # whether across the user defined overlap region

	def formatprimer(self):
		if self.direction == "LEFT_PRIMER":
			formatout = "\t".join(str(x) for x in [self.name, self.direction, self.start, self.end, self.length, self.tm, self.gc, self.anys, self.three, self.hairpin, self.nvar, self.difthree, self.difall, self.difnum, self.score, self.seq, ""])
		else:
			formatout = "\t".join(str(x) for x in [self.name, self.direction, self.start, self.end, self.length, self.tm, self.gc, self.anys, self.three, self.hairpin, self.nvar, self.difthree, self.difall, self.difnum, self.score, ReverseComplement(self.seq), self.seq])
		return(formatout)

class PrimerPair(object):
	"""A pair of primers designed by Primer3"""
	def __init__(self):
		self.name = ""
		self.left = Primers()
		self.right = Primers()
		self.compl_any = "NA"
		self.compl_end = "NA"
		self.penalty = "NA"
		self.product_size = 0
		self.score = 0 # sum of scores of left and right primers
		self.difsite = []
		self.difsite4 = []

# parse primer3 output file
# use 0-based primer3 output
def parse(handle, direction):
	"""Iterate over primer3 output
	"""
	primerlist = []
	with open(handle) as infile:
		startread = 0
		for line in infile:
			if targets[0] in line:
				startread = 1
			if startread and direction in line:
				pp = line.strip().split() # split on white spaces
				primer = Primers()
				primer.direction = direction
				if direction == "LEFT_PRIMER":
					primer.name = "L" + pp[0]
					primer.start = int(pp[2]) + 1
					primer.end = int(pp[3]) + int(pp[2])
					primer.seq = pp[9]
				else:
					primer.name = "R" + pp[0]
					primer.end = int(pp[2]) + 1
					primer.start = int(pp[2]) + 1 - int(pp[3]) + 1
					primer.seq = ReverseComplement(pp[9])
				primer.length = int(pp[3])
				primer.tm = float(pp[4])
				primer.gc = float(pp[5])
				primer.anys = pp[6]
				primer.three = pp[7]
				primer.hairpin = pp[8]
				if set(overlap_region) & set(range(primer.start, primer.end-1)):
					primer.overlap = "YES"
				primerlist.append(primer)
			if startread and "PRIMER PICKING RESULTS FOR" in line and targets[0] not in line:
				break
	return(primerlist)

# get left and right primers
leftprimers = parse(primer3output, "LEFT_PRIMER")
rightprimers = parse(primer3output, "RIGHT_PRIMER")

print "Length of LEFT primers:", len(leftprimers)
print "Length of RIGHT primers:", len(rightprimers)



###################################
# STEP 5: filter primers with variations

# selected left primers
newleftprimers = [] # left primers that are different between target group and the others but cannot differ all.
alldifferenceleft = [] # primers that already can differentiate all other sequences
nodiffleft = [] # primers that have no differences among all the seq, will be used to pair with primers that can already differ all other sequences.

for pp in leftprimers:
	pp = getvar(pp)
	if pp != 0:
		if pp.nvar > 0: # if there is difference
			alldifsite15 = [int(x > 3) for x in pp.difsite] # at least 4 differences
			alldifsite4 = [int(x > 0) for x in pp.difsite4] # at least 1 differences
			alldifsite = [sum(x) for x in zip(alldifsite15, alldifsite4)]
			if min(alldifsite) > 0:
				alldifferenceleft.append(pp)
			else:
				newleftprimers.append(pp)
		else:
			nodiffleft.append(pp)

print "number of left primers that can diff all:",
print len(alldifferenceleft)

print "Number of selected LEFT primers", len(newleftprimers)

# selected right primers
newrightprimers = []
alldifferenceright = [] # primers that already can differenciate all seqences
nodiffright = []

for pp in rightprimers:
	pp = getvar(pp)
	if pp != 0:
		if pp.nvar > 0: # if there is difference
			alldifsite15 = [int(x > 3) for x in pp.difsite] # at least 4 differences
			alldifsite4 = [int(x > 0) for x in pp.difsite4] # at least 1 differences
			alldifsite = [sum(x) for x in zip(alldifsite15, alldifsite4)]
			if min(alldifsite) > 0:
				alldifferenceright.append(pp)
			else:
				newrightprimers.append(pp)
		else:
			nodiffright.append(pp)

print "number of right primers that can diff all:",
print len(alldifferenceright)

print "Number of selected RIGHT primers", len(newrightprimers)

#############################
# STEP 5: Select and Test Primer Pairs

# selected primers pairs
tempin = 'temp_primer_pair_check_input_' + groupname + ".txt"
p3temp = open(tempin, 'w') # output file
seqtemplate = fasta[mainID].replace("-","") # remove "-" in the alignment seq
primernumber = 0
primerpairs = {} # all the pairs with the right size and Tm differences

## define function to merge two dictionaries of difsite
def merge_dict(dx, dy):
	newdict = {}
	for key in dx.keys() + dy.keys():
		if key in dx:
			if key in dy:
				newdict[key] = [sum(x) for x in zip(dx[key], dy[key])]
			else:
				newdict[key] = dx[key]
		else:
			newdict[key] = dy[key]
	return(newdict)


def testpair(leftlist, rightlist):
	global product_max, product_min, primernumber, primerpairs, p3temp, maxTmdiff
	for pl in leftlist:
		if len(primerpairs) > 5000: # in case too many pairs of primers need to check.
			break
		for pr in rightlist:
			alldifsite15 = [int(max(x) > 3) for x in zip(pl.difsite, pr.difsite)] # at least 4 differences
			alldifsite4 = [int(max(x) > 0) for x in zip(pl.difsite4, pr.difsite4)] # at least 1 differences
			alldifsite = [sum(x) for x in zip(alldifsite15, alldifsite4)]
			if min(alldifsite) > 0 and pl.end < pr.start and abs(pl.tm - pr.tm) <= maxTmdiff:
				productsize = pr.end - pl.start + 1
				if productsize >= product_min and productsize <= product_max:
					primernumber += 1
					ppname = str(primernumber)
					pp = PrimerPair()
					pp.name = ppname
					pp.left = pl
					pp.right = pr
					pp.product_size = productsize
					pp.difsite = alldifsite15
					pp.difsite4 = alldifsite4
					# get score below
					diffmerge = merge_dict(pl.difsitedict, pr.difsitedict)
					for k, v in diffmerge.items():
						difnum = sum(i > 0 for i in v) # how many sequences this difpos can differ
						tt = sum(i > 1 for i in v) # to account for the same variations between left and right primers
						#pp.score += (difnum / len(ids) + 1) * 100 / (k + 1)
						pp.score += (float(difnum) / len(ids) * 100 + float(tt) / len(ids) * 50)/ k
					if pp.score >= primer_pair_score_threshold:
						primerpairs[ppname] = pp
						line1 = "SEQUENCE_ID=" + ppname
						line2 = "PRIMER_TASK=check_primers"
						line3 = "PRIMER_EXPLAIN_FLAG=1"
						line4 = "PRIMER_PRODUCT_SIZE_RANGE=" + str(product_min) + "-" + str(product_max)
						line5 = "SEQUENCE_TEMPLATE=" + seqtemplate
						line6 =  "SEQUENCE_PRIMER=" + pl.seq
						line7 = "SEQUENCE_PRIMER_REVCOMP=" + ReverseComplement(pr.seq)
						line8 = "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + getprimer_path + "/bin/primer3_config/"
						line9 = "="
						p3temp.write("\n".join([line1, line2, line4, line5, line6, line7, line8, line9]) + "\n")
	print "primernumber is ", primernumber
	print "Candidate primer pairs: ", len(primerpairs)

testpair(alldifferenceleft, alldifferenceright)
if len(primerpairs) < 1000:
	testpair(alldifferenceleft, newrightprimers)
if len(primerpairs) < 1000:
	testpair(newleftprimers, alldifferenceright)
if len(primerpairs) < 1000:
	testpair(alldifferenceleft, nodiffright)
if len(primerpairs) < 1000:
	testpair(nodiffleft, alldifferenceright)

p3temp.close()

# check to see whether no good primer pairs found
#print "Candidate primer pairs: ", len(primerpairs)
if not primerpairs:
	print "\nNo GOOD primers found!"
	sys.exit(1)

## use primer3 to check the primer pair quality
tempout = "temp_primer_pair_test_out_" + groupname + ".txt"
#p3cmd = primer3_path + " -default_version=2 -p3_settings_file=" + primer3_parameter_path + " -output=" + " ".join([tempout, tempin])
p3cmd = primer3_path + " -default_version=2 -output=" + " ".join([tempout, tempin])
print "Primer3 command 2nd time: ", p3cmd
call(p3cmd, shell=True)

### parse primer 3 primer check output

with open(tempout) as infile:
	for line in infile:
		line = line.strip()
		if "SEQUENCE_ID" in line:
			seqid = line.split("=")[1]
		if "PRIMER_PAIR_0_PENALTY" in line:
			primerpairs[seqid].penalty = line.split("=")[1]
		if "PRIMER_PAIR_0_COMPL_ANY" in line:
			primerpairs[seqid].compl_any = line.split("=")[1]
		if "PRIMER_PAIR_0_COMPL_END" in line:
			primerpairs[seqid].compl_end = line.split("=")[1]

outfile = open(out, 'w') # output file
outfile.write("index\tproduct_size\tprimerID\ttype\tstart\tend\tlength\tTm\tGCcontent\tany\t3'\thairpin\tprimer_nvar\t3'Diff\tDiffAll\tDifNumber\tprimer_score\tprimer_seq\tReverseComplement\tpenalty\tcompl_any\tcompl_end\tprimerpair_score\tprimer_diff15\tprimer_diff4\tacross_overlap\tmatched_chromosomes\n")

print "primer_pair_compl_any_threshold ", primer_pair_compl_any_threshold
print "primer_pair_compl_end_threshold ", primer_pair_compl_end_threshold

#for pp in primerpairs:
pp_vector = primerpairs.values()
pp_vector.sort(key=lambda x: x.score, reverse=True)
exist_left = []
exist_right = []
final_primers = [] # for final output
primer_for_blast = {} # primer sequences for blast
for pp in pp_vector:
	if len(final_primers) > 50: # only get 51 primer pairs maximum. Otherwise too many.
		break
	if pp.compl_any != "NA" and float(pp.compl_any) <= primer_pair_compl_any_threshold and float(pp.compl_end) <= primer_pair_compl_end_threshold:
		pl = pp.left
		pr = pp.right
		if filter_flag:
			if pl.end not in exist_left or pr.start not in exist_right:
				#forblast.write(">" + pl.name + "\n" + pl.seq + "\n>" + pr.name + "\n" + pr.seq + "\n")
				primer_for_blast[pl.name] = pl.seq
				primer_for_blast[pr.name] = ReverseComplement(pr.seq) # right prmer sequences is actually its RC
				final_primers.append(pp)
				#outfile.write("\t".join([pp.name, str(pp.product_size), pl.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end, str(pp.score), str(pl.difsite), str(pl.difsite4), pl.overlap]) + "\n")
				#outfile.write("\t".join([pp.name, str(pp.product_size), pr.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end, str(pp.score), str(pr.difsite), str(pr.difsite4), pr.overlap]) + "\n")
				exist_left += range(pl.end - 5, pl.end + 6)
				exist_right += range(pr.start -5, pr.start + 6)
		else:
			#forblast.write(">" + pl.name + "\n" + pl.seq + "\n>" + pr.name + "\n" + pr.seq + "\n")
			primer_for_blast[pl.name] = pl.seq
			primer_for_blast[pr.name] = pr.seq
			final_primers.append(pp)
			#outfile.write("\t".join([pp.name, str(pp.product_size), pl.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end, str(pp.score), str(pl.difsite), str(pl.difsite4), pl.overlap]) + "\n")
			#outfile.write("\t".join([pp.name, str(pp.product_size), pr.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end, str(pp.score), str(pr.difsite), str(pr.difsite4), pr.overlap]) + "\n")

#########
#outfile.close()
forblast = open("for_blast.fa", 'w') # for blast against the gnome
for k, v in primer_for_blast.items():
	forblast.write(">" + k + "\n" + v + "\n")
forblast.close()

blast_hit = {} # matched chromosomes for primers: at least perfect match for the first 5 bp from 3'
### for blast

# function to count mismtaches
def mismatchn (s1, s2):
	return sum(c1!=c2 for c1,c2 in zip(s1,s2))

reference = "/Library/WebServer/Documents/blast/db/nucleotide/161010_Chinese_Spring_v1.0_pseudomolecules.fasta"
if blast:
	cmd2 = 'blastn -task blastn -db ' + reference + ' -query for_blast.fa -outfmt "6 std qseq sseq qlen slen" -num_threads 3 -word_size 7 -out blast_out.txt'
	print "Step 2: Blast command:\n", cmd2
	call(cmd2, shell=True)
	# process blast file
	# blast fields
	# IWB50236_7A_R	IWGSC_CSS_7DS_scaff_3919748	98.718	78	1	0	24	101	4891	4968	1.55e-30	138	CTCATCAAATGATTCAAAAATATCGATRCTTGGCTGGTGTATCGTGCAGACGACAGTTCGTCCGGTATCAACAGCATT	CTCATCAAATGATTCAAAAATATCGATGCTTGGCTGGTGTATCGTGCAGACGACAGTTCGTCCGGTATCAACAGCATT	101 5924
	# Fields: 
	# 1: query id, subject id, % identity, alignment length, mismatches, gap opens, 
	# 7: q. start, q. end, s. start, s. end, evalue, bit score
	# 13: q. sequence, s. sequence, q. length s. length
	for line in open("blast_out.txt"):
		if line.startswith('#'):
			continue
		fields = line.split("\t")
		query, subject, pct_identity, align_length= fields[:4]
		qstart, qstop, sstart, sstop = [int(x) for x in fields[6:10]]
		qseq, sseq = fields[12:14]
		qlen = int(fields[14])
		n1 = qlen - qstop
		if n1 < 2 and mismatchn(qseq[(n1 - 4):], sseq[(n1 - 4):]) + n1 < 2: # if less than 2 mismtaches in the first 4 bases from the 3' end of the primer
			blast_hit[query] = blast_hit.setdefault(query, "") + ";" + subject + ":" + str(sstart)

## write output
for pp in final_primers:
	pl = pp.left
	pr = pp.right
	outfile.write("\t".join([pp.name, str(pp.product_size), pl.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end, str(pp.score), str(pl.difsite), str(pl.difsite4), pl.overlap, blast_hit.setdefault(pl.name, "")]) + "\n")
	outfile.write("\t".join([pp.name, str(pp.product_size), pr.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end, str(pp.score), str(pr.difsite), str(pr.difsite4), pr.overlap, blast_hit.setdefault(pr.name, "")]) + "\n")

outfile.close()

## remove all tempotary files
#call("rm temp_*", shell=True)
