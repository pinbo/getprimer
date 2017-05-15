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

# 5/9/2017: change a lot of steps to functions and add main function


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

###########################
# functions to use

#from sys import platform
def get_software_path(base_path):
	if sys.platform.startswith('linux'): # linux
		primer3_path = base_path + "/bin/primer3_core"
		muscle_path = base_path + "/bin/muscle"
	elif sys.platform == "win32" or sys.platform == "cygwin": # Windows...
		primer3_path = base_path + "/bin/primer3_core.exe"
		muscle_path = base_path + "/bin/muscle.exe"
	elif sys.platform == "darwin": # MacOSX
		primer3_path = base_path + "/bin/primer3_core_darwin64"
		muscle_path = base_path + "/bin/muscle3.8.31_i86darwin64"
	return primer3_path, muscle_path

# simple Tm calculator
def Tm(seq):
	t=0
	for a in seq:
		if a=='A' or a=='T':
			t=t+2
		if a=='C' or a=='G':
			t=t+4
	return t

# calculate GC content of a sequence
def Calc_GC(seq):
	t = 0.0 # float
	for a in seq:
		if a=='C' or a=='G':
			t += 1
	return t / len(seq) * 100

# function to find the segments with largest Tm
def FindLongestSubstring(s1, s2):
	longestStart = 0
	longestEnd = 0
	largestTm = 0
	start = 0
	gaps = [i for i, c in enumerate(s1) if c=='-' or s2[i]=='-']
	gaps.append(len(s1))
	for gap in gaps:
		end = gap
		tm = Tm(s1[start:end])
		if  tm > largestTm:
			longestStart = start
			longestEnd = end
			largestTm = tm
		start = gap + 1
	nL = len(s1[:longestStart].replace("-","")) # number of bases on the left of the longest segments
	nR = len(s1[longestEnd:].replace("-",""))  # number of bases on the right of the longest segments
	return [longestStart, longestEnd, nL, nR]

# function to creat alignment file
def get_alignment(muscle_path, infile, outfile):
	alignmentcmd = muscle_path + " -in " + infile + " -out " + outfile + " -quiet"
	print "Alignment command: ", alignmentcmd
	call(alignmentcmd, shell=True)

# function to extract sequences from a fasta file 
def get_fasta(infile):
	fasta = {} # dictionary for alignment
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if line.startswith(">"):
				sequence_name = line.split()[0].lstrip(">")
			else:
				fasta.setdefault(sequence_name, "")
				fasta[sequence_name] += line.rstrip()
	return fasta

# format alignment sequences by removing leading and ending "-"
def format_alignment(fasta):
	max_gap_left = 0
	max_gap_right = 0
	for k, v in fasta.items():
		nL = len(v) - len(v.lstrip("-"))
		nR = len(v) - len(v.rstrip("-"))
		if nL > max_gap_left:
			max_gap_left = nL
		if nR > max_gap_right:
			max_gap_right = nR
	print "Gap left: ", max_gap_left
	print "Gap right: ", max_gap_right
	### striping sequences to equal length
	for k, v in fasta.items():
		fasta[k] = v[max_gap_left:(len(v) - max_gap_right)]
	return fasta

# Creat primer3 input and run primer3 to create primer lists
def get_primers(primer3_param, primer3_path, primer3output):
	primer3input = "primer3.input"
	p3input = open(primer3input, 'w')
	# I only need to use the output from the first target
	lines = []
	lines.append("PRIMER_TASK=pick_primer_list")
	lines.append("SEQUENCE_ID=" + primer3_param["seqID"])
	lines.append("SEQUENCE_TEMPLATE=" + primer3_param["seq"])
	lines.append("PRIMER_FIRST_BASE_INDEX=1")
	lines.append("PRIMER_PRODUCT_SIZE_RANGE=" + primer3_param["product_range"])
	lines.append("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + primer3_param["TH_param_path"])
	lines.append("PRIMER_NUM_RETURN=10000")
	lines.append("PRIMER_MIN_SIZE=" + primer3_param["primer_size_min"])
	lines.append("PRIMER_MAX_SIZE=" + primer3_param["primer_size_max"])
	lines.append("PRIMER_MIN_TM=" + primer3_param["primer_tm_min"])
	lines.append("PRIMER_MAX_TM=" + primer3_param["primer_tm_max"])
	lines.append("PRIMER_PAIR_MAX_DIFF_TM=" + primer3_param["primer_tm_diff_max"])
	lines.append("PRIMER_MAX_HAIRPIN_TH=" + primer3_param["primer_hairpin_max"])
	lines.append("=\n")

	p3input.write("\n".join(lines))
	p3input.close()
	# primer3 output file
	#p3cmd = primer3_path + " -default_version=2 -format_output -p3_settings_file=" + primer3_parameter_path + " -output=" + primer3output + " " + primer3input
	p3cmd = primer3_path + " -default_version=2 -format_output -output=" + primer3output + " " + primer3input
	print "Primer3 command 1st time: ", p3cmd
	call(p3cmd, shell=True)

# function to get reverse complement
def ReverseComplement(seq):
	# too lazy to construct the dictionary manually, use a dict comprehension
	seq1 = 'ATCGTAGCatcgtagc'
	seq_dict = { seq1[i]:seq1[i+4] for i in range(16) if i < 4 or 8<=i<12 }
	return "".join([seq_dict[base] for base in reversed(seq)])

# parse primer3 output file
# use 0-based primer3 output
def parse_primer3_output(handle, direction, overlap_region):
	"""Iterate over primer3 output
	"""
	primerlist = []
	with open(handle) as infile:
		for line in infile:
			if direction in line:
				pp = line.strip().split() # split on white spaces
				primer = Primers()
				primer.direction = direction
				if direction == "LEFT_PRIMER":
					primer.name = "L" + pp[0]
					primer.start = int(pp[2])
					primer.end = int(pp[3]) + int(pp[2]) - 1
					primer.seq = pp[9]
				else:
					primer.name = "R" + pp[0]
					primer.end = int(pp[2])
					primer.start = int(pp[2]) - int(pp[3]) + 1
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
	return(primerlist)

# see whether the primer has difference within target group
def exclude_primer(pp, targets, fasta, align_left, align_right):
	exclude = 0
	for j in targets:#check difference within targets
		targetseq = fasta[j][align_left:(align_right + 1)].replace("-","")
		if targetseq != pp.seq:
			exclude = 1
			break
	return exclude

# get the list of sequences in the homeolog groups for comparison with current primer
def get_homeo_seq(fasta, mainID, ids, align_left, align_right):
	s1 = fasta[mainID] # primer sequence in the template with gaps
	seq2comp = [] # final sequence to compare for each other homeolog
	for k in ids:
		s2 = fasta[k]
		targetSeq = s1[align_left:(align_right + 1)]
		homeoSeq = s2[align_left:(align_right + 1)]
		# Get the sequences for comparison
		indexL, indexR, nL, nR = FindLongestSubstring(targetSeq, homeoSeq)
		indexL += align_left
		indexR += align_left
		seqL = s2[:indexL].replace("-","")
		seqR = s2[indexR:].replace("-","")
		print "indexL, indexR, nL, nR ", indexL, indexR, nL, nR
		print "s2[indexL:indexR] ", s2[indexL:indexR]
		if len(seqL) < nL: # in case it does not have enough bases
			seqL = "-" * (nL - len(seqL)) + seqL
		if len(seqR) < nR:
			seqL = seqR + "-" * (nR - len(seqR))
		seqk = seqL[::-1][:nL][::-1] + s2[indexL:indexR] + seqR[:nR]
		seq2comp.append(seqk)
		print k, "\t", seqk
	return seq2comp

# group primers into different category: differ all, differ some, or no differ from the homeolog
def group_primers(primerlist, fasta, t2a, targets, ids):
	new_primers = [] # primers that are different between target group and the others but cannot differ all.
	primer_diff_all = [] # primers that already can differentiate all other sequences
	primer_no_diff = [] # primers cannot differ any homeolog sequences
	for pp in primerlist:
		pp = getvar(pp, fasta, t2a, targets, ids)
		if pp != 0:
			if pp.nvar > 0: # if there is difference
				alldifsite15 = [int(x > 3) for x in pp.difsite] # at least 4 differences
				alldifsite4 = [int(x > 0) for x in pp.difsite4] # at least 1 differences
				alldifsite = [sum(x) for x in zip(alldifsite15, alldifsite4)]
				if min(alldifsite) > 0:
					primer_diff_all.append(pp)
				else:
					new_primers.append(pp)
			else:
				primer_no_diff.append(pp)
	return new_primers, primer_diff_all, primer_no_diff

# function to get primer variations
def getvar(pp, fasta, t2a, targets, ids): # pp is a Primer object
	print "pp.seq ", pp.seq
	pp.difsite = [0] * len(ids) # initialize the number of differnt bases from other sequences
	pp.difsite4 = [0] * len(ids) # initialize the number of differnt bases from other sequences
	mainID = targets[0] # the one whose primer3 output will be used
	align_left = t2a[pp.start - 1] # pp.start and end is 1-based
	align_right = t2a[pp.end - 1]
	if exclude_primer(pp, targets, fasta, align_left, align_right): # if not match all targets
		return 0
	# if it is good for all targets
	print "pp.seq\t", pp.seq
	### Idea: split the sequences by "-", find the piece that has the largest Tm, then extended that.
	seq2comp = get_homeo_seq(fasta, mainID, ids, align_left, align_right) # list of sequences for comparison
	diffarray = [] # to record difference with each other seq in each site, for primer pair selection later
	for i in range(pp.length):
		da = [0] * len(ids) # differ array for each base
		m = 0 # counter of next loop
		for k in seq2comp:
			# below nl1 to nr2 is for handling gaps
			b1 = pp.seq[i] # target non-gap base
			b2 = k[i] # homeolog non-gap bas
			if b1 != b2:
				da[m] = 1 # m sequence has variation from target
			m += 1
		# add the da to diffarray, diffarray now a list of difference for each homeolog	
		diffarray.append(da)
		if min(da) > 0: # if there is variation between targets and other sequences
			if pp.direction == "LEFT_PRIMER":
				difpos = pp.length - i + 1
				pp.difsitedict[difpos] = da
			else:
				difpos = i + 1
				pp.difsitedict[difpos] = da
			pp.nvar += 1
			difnum = sum(i > 0 for i in da)
			pp.difsite = [sum(x) for x in zip(pp.difsite, da)]
			difnum_new = sum(i > 0 for i in pp.difsite)
			dif_add = difnum_new - difnum # whether it can increase the overal variation
			pp.score += (float(difnum) / len(ids) * 100 + float(dif_add) / len(ids) * 50) / difpos
	# get differsite4 and differsite15
	pp.nvar =  sum([max(x) for x in zip(*diffarray)])
	pp.difnum = sum(i > 0 for i in pp.difsite)
	if min(pp.difsite) > 0:
		pp.difall = "YES"
	if pp.direction == "LEFT_PRIMER":
		pp.difsite = [sum(x) for x in zip(*diffarray[(pp.length - 15):pp.length])]
		pp.difsite4 = [sum(x) for x in zip(*diffarray[(pp.length - 4):pp.length])]
		if sum(diffarray[-1]) > 0:
			pp.difthree = "YES"
	else:
		pp.difsite = [sum(x) for x in zip(*diffarray[0:15])]
		pp.difsite4 = [sum(x) for x in zip(*diffarray[0:4])]
		if sum(diffarray[0]) > 0:
			pp.difthree = "YES"
	return pp

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

# test primer pairs
def testpair(leftlist, rightlist, primer3_param, primerpairs):
	ids = primer3_param["homeologs"] # list of other homeolog names
	primer_pair_score_threshold = primer3_param["primer_pair_score_threshold"]
	template = primer3_param["seq"]
	if primerpairs:
		primernumber = max(int(k) for k in primerpairs) + 1
	else:
		primernumber = 0
	product_min, product_max = [int(x) for x in primer3_param["product_range"].split("-")]
	maxTmdiff = primer3_param["primer_tm_diff_max"]
	for pl in leftlist:
		if len(primerpairs) > 1000: # in case too many pairs of primers need to check.
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
					amplicon = template[pl.start - 1:pr.end]
					pp.amplicon = amplicon
					pp.ampliconGC = Calc_GC(amplicon)
					# get score below
					diffmerge = merge_dict(pl.difsitedict, pr.difsitedict)
					for k, v in diffmerge.items():
						difnum = sum(i > 0 for i in v) # how many sequences this difpos can differ
						tt = sum(i > 1 for i in v) # to account for the same variations between left and right primers
						pp.score += (float(difnum) / len(ids) * 100 + float(tt) / len(ids) * 50)/ k
					if pp.score >= primer_pair_score_threshold:
						primerpairs[ppname] = pp
	print "primernumber is ", primernumber
	print "Candidate primer pairs: ", len(primerpairs)
	return primerpairs

def prepare_primerpair_check_input(primerpairs, primer3_param, tempin):
	p3temp = open(tempin, 'w') # output file
	for nn, pp in primerpairs.items():
		pl = pp.left
		pr = pp.right
		line1 = "SEQUENCE_ID=" + pp.name
		line2 = "PRIMER_TASK=check_primers"
		line3 = "PRIMER_EXPLAIN_FLAG=1"
		line4 = "PRIMER_PRODUCT_SIZE_RANGE=" + primer3_param["product_range"]
		line5 = "SEQUENCE_TEMPLATE=" + primer3_param["seq"]
		line6 =  "SEQUENCE_PRIMER=" + pl.seq
		line7 = "SEQUENCE_PRIMER_REVCOMP=" + ReverseComplement(pr.seq)
		line8 = "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + primer3_param["TH_param_path"]
		line9 = "="
		p3temp.write("\n".join([line1, line2, line4, line5, line6, line7, line8, line9]) + "\n")
	p3temp.close()

def filter_primerpairs_for_blast(primerpairs, primer_pair_compl_any_threshold, primer_pair_compl_end_threshold, filter_flag):
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
					primer_for_blast[pl.name] = pl.seq
					primer_for_blast[pr.name] = ReverseComplement(pr.seq) # right prmer sequences is actually its RC
					final_primers.append(pp)
					exist_left += range(pl.end - 2, pl.end + 4) # only filter out closest 3 bases
					exist_right += range(pr.start - 3, pr.start + 3)
			else:
				primer_for_blast[pl.name] = pl.seq
				primer_for_blast[pr.name] = pr.seq
				final_primers.append(pp)
	return final_primers, primer_for_blast

# function to count mismtaches
def mismatchn (s1, s2):
	return sum(c1!=c2 for c1,c2 in zip(s1,s2))

# function to blast and parse the output of each primer in the wheat genome
def primer_blast(primer_for_blast):
	forblast = open("for_blast.fa", 'w') # for blast against the gnome
	for k, v in primer_for_blast.items():
		forblast.write(">" + k + "\n" + v + "\n")
	forblast.close()
	blast_hit = {} # matched chromosomes for primers: less than 2 mismatches in the first 4 bps from 3'
	### for blast
	reference = "/Library/WebServer/Documents/blast/db/nucleotide/161010_Chinese_Spring_v1.0_pseudomolecules.fasta"
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
	return blast_hit

##############################
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
		self.amplicon = ""
		self.ampliconGC = 0

#########################

def main():
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
	out = "" # output file
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
	print "Options done"
	if not targets:
		print "Please give the sequence IDs for primer design!"
		sys.exit(1)

	groupname = "-".join(targets)
	if not out:
		out = 'selected_primers_for_' + groupname + ".txt"

	primer3_path, muscle_path = get_software_path(getprimer_path)
	#########################
	# STEP 0: create alignment file and primer3output file
	RawAlignFile = "alignment_raw.fa"
	if msa:
		get_alignment(muscle_path, seqfile, RawAlignFile)

	###################
	# STEP 1: read alignment fasta file into a dictionary AND format it by removing leading and ending "-"
	fasta = get_fasta(RawAlignFile) # read alignment file
	fasta = format_alignment(fasta) # remove leading and ending "-"

	####################################
	# STEP 2: design primers
	mainID = targets[0]
	primer3_param = {
		"seqID": mainID,
		"seq": fasta[mainID].replace("-",""), # remove "-" in the alignment seq
		"product_range": str(product_min) + "-" + str(product_max),
		"TH_param_path": getprimer_path + "/bin/primer3_config/",
		"primer_size_min": str(minSize),
		"primer_size_max": str(maxSize),
		"primer_tm_min":  str(minTm),
		"primer_tm_max": str(maxTm),
		"primer_tm_diff_max": str(maxTmdiff),
		"primer_hairpin_max": str(maxhairpin),
		"primer_pair_score_threshold": primer_pair_score_threshold
		}
	# Creat primer3 input and run primer3 to create primer lists
	primer3output = "primer3.output"
	get_primers(primer3_param, primer3_path, primer3output)

	#####################################
	# STEP 3: parse primer3 output
	leftprimers = parse_primer3_output(primer3output, "LEFT_PRIMER", overlap_region)
	rightprimers = parse_primer3_output(primer3output, "RIGHT_PRIMER", overlap_region)

	print "Length of LEFT primers:", len(leftprimers)
	print "Length of RIGHT primers:", len(rightprimers)

	###########################
	# STEP 4: get the primer start and end positions on the alignment
	ids = [] # all other sequence names
	for key in fasta:
		if key not in targets:
			ids.append(key)

	print "The other groups: ", ids
	primer3_param["homeologs"] = ids

	alignlen = len(fasta[targets[0]])
	print "Alignment length: ", alignlen

	# get the target ID template base coordinate in the alignment
	t2a = {} # template to alignment
	ngap = 0 # gaps
	for i in range(alignlen):
		if fasta[mainID][i] == "-":
			ngap += 1
			continue
		t2a[i - ngap] = i

	print "last key of t2a", i - ngap

	###################################
	# STEP 5: filter primers with variations

	# selected left primers
	newleftprimers, alldifferenceleft, nodiffleft = group_primers(leftprimers, fasta, t2a, targets, ids)
	print "number of left primers that can diff all:",
	print len(alldifferenceleft)
	print "Number of selected LEFT primers", len(newleftprimers)
	# selected right primers
	newrightprimers, alldifferenceright, nodiffright = group_primers(rightprimers, fasta, t2a, targets, ids)
	print "number of right primers that can diff all:",
	print len(alldifferenceright)
	print "Number of selected RIGHT primers", len(newrightprimers)

	#############################
	# STEP 6: Select Primer Pairs

	# selected primers pairs
	primerpairs = {} # all the pairs with the right size and Tm differences
	# test primer pairs
	primerpairs = testpair(alldifferenceleft, alldifferenceright, primer3_param, primerpairs)
	if len(primerpairs) < 1000:
		primerpairs = testpair(alldifferenceleft, newrightprimers, primer3_param, primerpairs)
	if len(primerpairs) < 1000:
		primerpairs = testpair(newleftprimers, alldifferenceright, primer3_param, primerpairs)
	if len(primerpairs) < 1000:
		primerpairs = testpair(alldifferenceleft, nodiffright, primer3_param, primerpairs)
	if len(primerpairs) < 1000:
		primerpairs = testpair(nodiffleft, alldifferenceright, primer3_param, primerpairs)
	# check to see whether no good primer pairs found
	if not primerpairs:
		print "\nNo GOOD primers found!"
		sys.exit(1)

	#############################
	# STEP 7: Check Primer Pairs quality
	tempin = "temp_primer_pair_check_input_" + groupname + ".txt"
	prepare_primerpair_check_input(primerpairs, primer3_param, tempin)
	tempout = "temp_primer_pair_test_out_" + groupname + ".txt"
	p3cmd = primer3_path + " -default_version=2 -output=" + " ".join([tempout, tempin])
	print "Primer3 command 2nd time: ", p3cmd
	call(p3cmd, shell=True)

	# parse primer 3 primer check output
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

	#############################
	# STEP 8: filter and write output

	# filter primer pairs if filter_flag == 1 and get the final list of primers
	final_primers, primer_for_blast = filter_primerpairs_for_blast(primerpairs, primer_pair_compl_any_threshold, primer_pair_compl_end_threshold, filter_flag)

	#############################
	# STEP 9: blast the primers in the wheat genome if blast parameter
	blast_hit = {}
	if blast:
		blast_hit = primer_blast(primer_for_blast) # chromosome hit for each primer

	#############################
	# STEP 10: write output
	outfile = open(out, 'w') # output file
	outfile.write("index\tproduct_size\tprimerID\ttype\tstart\tend\tlength\tTm\tGCcontent\tany\t3'\thairpin\tprimer_nvar\t3'Diff\tDiffAll\tDifNumber\tprimer_score\tprimer_seq\tReverseComplement\tpenalty\tcompl_any\tcompl_end\tprimerpair_score\tprimer_diff15\tprimer_diff4\tacross_overlap\tampliconGC\tmatched_chromosomes\n")

	for pp in final_primers:
		pl = pp.left
		pr = pp.right
		outfile.write("\t".join([pp.name, str(pp.product_size), pl.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end, str(pp.score), str(pl.difsite), str(pl.difsite4), pl.overlap, str(pp.ampliconGC), blast_hit.setdefault(pl.name, "")]) + "\n")
		outfile.write("\t".join([pp.name, str(pp.product_size), pr.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end, str(pp.score), str(pr.difsite), str(pr.difsite4), pr.overlap, str(pp.ampliconGC), blast_hit.setdefault(pr.name, "")]) + "\n")

	outfile.close()

	print "\n\nPrimer design is finished!\n\n"
	## remove all tempotary files
	#call("rm temp_*", shell=True)
	return 0

if __name__ == '__main__':
    sys.exit(main())
