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
### Imported
from subprocess import call




#########################
# parameter or file names that need to be changed
# Gene groups (5 here)
# ["Chr-B2.2"]
# ["Chr-B2.3"]
# ["Chr-D2.2", "Chr-A2.2"]
# ["Chr-D2.3", "Chr-A2.3"]
# ["Chr-B2.1", "Chr-D2.1", "Chr-A2.1"]

seqfile = "sequence.fa" # CHANGE HERE
targets = ["Chr-D2.2", "Chr-A2.2"] # CHAGE HERE
groupname = "Chr2.2AD" # CHANGE HERE
primer_pair_score_threshold = 200.0
primer_pair_compl_any_threshold = 1.0
primer_pair_compl_end_threshold = 1.0
product_min = 50
product_opt = 100
product_max = 150
muscle_path = "./bin/muscle"
primer3_path = "./bin/primer3_core"
rangelimit = 10 # only find difference in the first a few nt from 3' end of each primer
outfile = open('selected_primers_for_' + groupname + ".txt", 'w') # output file
#########################
# STEP 0: create alignment file and primer3output file

RawAlignFile = "alignment_raw.fa"
alignmentcmd = muscle_path + " -in " + seqfile + " -out " + RawAlignFile + " -quiet"
print "Alignment command: ", alignmentcmd
call(alignmentcmd, shell=True)

###################
# STEP 1: read alignment fasta file into a dictionary AND format it by removing leading and ending "-"
fasta = {} # dictionary for alignment

with open(RawAlignFile) as file_one:
	for line in file_one:
		if line.startswith(">"):
			sequence_name = line.rstrip().lstrip(">")
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


variation = [] # variation sites
varexclude = [] # varaition among targets, any primers that have these should be excluded. 
diffarray = {} # to record difference with each other seq in each site, for primer pair selection later
ngap = 0 # gaps
mainID = targets[0] # the one whose primer3 output will be used
for i in range(alignlen):
	nd = 0 # number of difference
	dd = [] # to record site difference for each seq
	if fasta[mainID][i] == "-":
		ngap += 1
	for j in targets:#check difference within targets
		if fasta[j][i] == fasta[mainID][i]: nd+=1
	if nd == len(targets): # if the site i has no differences within targets
		for k in ids:
			if fasta[k][i] == fasta[mainID][i]:
				nd+=1
				dd.append(0)
			else:
				dd.append(1)
		if nd < len(targets) + len(ids): # if there is variation
			if fasta[mainID][i] == "-":
				coordinates = (i - ngap,i - ngap + 1) # coordinates
				variation.append(coordinates)
				if coordinates not in diffarray:
					diffarray[coordinates] = dd
				else:
					diffarray[coordinates] = [sum(x) for x in zip(dd,diffarray[coordinates])]
			else:
				variation.append((i - ngap,)) # shift due to gap
				diffarray[(i - ngap,)] = dd
	else: # if there is difference among the targets
		if fasta[mainID][i] == "-":
			coordinates = (i - ngap,i - ngap + 1) # coordinates
			varexclude.append(coordinates)
		else:
			varexclude.append((i - ngap,)) # shift due to gap	

####################################
# STEP 3: create primer3.input and run it for output

# primer3 inputfile
primer3input = "primer3.input"
p3input = open(primer3input, 'w')
"""
for k, v in fasta.items():
	line1 = "SEQUENCE_ID=" + k
	line2 = "SEQUENCE_TEMPLATE=" + v.replace("-","") # remove "-" in the alignment seq
	line3 = "PRIMER_PRODUCT_OPT_SIZE=" + str(product_opt)
	line4 = "PRIMER_PRODUCT_SIZE_RANGE=" + str(product_min) + "-" + str(product_max)
	line5 = "="
	p3input.write("\n".join([line1, line2, line3, line4, line5]) + "\n")
"""
# I only need to use the output from the first target
line1 = "SEQUENCE_ID=" + mainID
line2 = "SEQUENCE_TEMPLATE=" + fasta[mainID].replace("-","") # remove "-" in the alignment seq
line3 = "PRIMER_PRODUCT_OPT_SIZE=" + str(product_opt)
line4 = "PRIMER_PRODUCT_SIZE_RANGE=" + str(product_min) + "-" + str(product_max)
line5 = "="
p3input.write("\n".join([line1, line2, line3, line4, line5]) + "\n")

p3input.close()

# primer3 output file
primer3output = "primer3.output"
p3cmd = primer3_path + " -default_version=2 -format_output -p3_settings_file=./primer3web_v4_JZ.txt -output=" + primer3output + " " + primer3input
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
		self.difsite = []
		self.nvar = 0 # number of variation from other sequences
		self.difthree = "NO" # whether 3' site is different
		self.difall = "NO" # whether one primer is enough to differ all the rest
		self.difnum = 0 # how many sequences can be differentiated within the rangelimit
		self.score = 0.0 # score based on number of different sites and differ position
		self.difsitedict = {} # difsite in each different position (just like 2-dim array): self.difsite would be sum of cols of self.difsitedict

	def __len__(self):
		"""Length of the primer"""
		return self.length
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
newleftprimers = [] # left primers that are different between target group and the others
alldifferenceleft = [] # primers that already can differentiate all other sequences
nodiffleft = [] # primers that have no differences among all the seq, will be used to pair with primers that can already differ all other sequences.

for pp in leftprimers:
	exclude = 0 # a tag to see whether the primer has difference within target group
	for var in varexclude: # check whether this primer has sites different within targets
		if set(var) < set(range(pp.start-1, pp.end)):
			exclude += 1
	if exclude == 0: # only if the primer does not have any exlude sites
		for var in variation:
			#if set([var]) < set(range(pp.start-1, pp.end)):
			if set(var) < set(range(pp.end - rangelimit, pp.end)): # pp.start and pp.end is 1-based, not 0-based
				difpos = pp.end - var[-1] # position of different site from the end: for calculating score
				difnum = sum(i > 0 for i in diffarray[var]) # how many sequences this difpos can differ
				pp.score += (difnum / len(ids) + 1) * 100 / difpos # consider both the difpos and how many it can differ
				pp.difsitedict[difpos] = diffarray[var]
				pp.nvar += 1
				if pp.difsite:
					pp.difsite = [sum(x) for x in zip(pp.difsite, diffarray[var])]
				else:
					pp.difsite = diffarray[var]
			if var == (pp.end-1,):
				pp.difthree = "YES"
		if pp.nvar > 0: # if there is difference
			newleftprimers.append(pp)
			pp.difnum = sum(i > 0 for i in pp.difsite) # count how many sequences can be differentiated within the rangelimit
			if min(pp.difsite) > 0:
				alldifferenceleft.append(pp)
				pp.difall = "YES"
		else:
			nodiffleft.append(pp)

print "number of left primers that can diff all:",
print len(alldifferenceleft)

#print newleftprimers[0].seq
#print newleftprimers[0].nvar
#print newleftprimers[0].difsite
#print newleftprimers[0].end
print "Number of selected LEFT primers", len(newleftprimers)

# selected right primers
newrightprimers = []
alldifferenceright = [] # primers that already can differenciate all seqences
nodiffright = []

for pp in rightprimers:
	exclude = 0
	for var in varexclude: # check whether this primer has sites different within targets
		if set(var) < set(range(pp.start-1, pp.end)):
			exclude += 1
	if exclude == 0: # only if the primer does not have any exlude sites
		for var in variation:
			#if set([var]) < set(range(pp.start-1, pp.end)):
			if set(var) < set(range(pp.start-1, pp.start-1 + rangelimit)):
				difpos = var[0] - pp.start + 2  # position of different site from the end: for calculating score
				difnum = sum(i > 0 for i in diffarray[var]) # how many sequences this difpos can differ
				pp.score += (difnum / len(ids) + 1) * 100 / difpos # consider both the difpos and how many it can differ
				pp.difsitedict[difpos] = diffarray[var]
				pp.nvar += 1
				if pp.difsite:
					pp.difsite = [sum(x) for x in zip(pp.difsite, diffarray[var])]
				else:
					pp.difsite = diffarray[var]
			if var == (pp.start-1,): # see whether the variation site is the 3' first site
				pp.difthree = "YES"
		if pp.nvar > 0:
			newrightprimers.append(pp)
			pp.difnum = sum(i > 0 for i in pp.difsite)
			if min(pp.difsite) > 0:
				alldifferenceright.append(pp)
				pp.difall = "YES"
		else:
			nodiffright.append(pp)

print "number of right primers that can diff all:",
print len(alldifferenceright)

#print newrightprimers[0].seq
#print newrightprimers[0].nvar
#print newrightprimers[0].difsite
#print newrightprimers[0].start



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

for pl in newleftprimers:
	for pr in newrightprimers:
		alldifsite = [sum(x) for x in zip(pl.difsite, pr.difsite)]
		if min(alldifsite) > 0 and pl.end < pr.start and abs(pl.tm - pr.tm) < 1.0:
			#print primernumber
			productsize = pr.end - pl.start + 1
			if productsize >= product_min and productsize <= product_max:
				primernumber += 1
				ppname = str(primernumber)
				pp = PrimerPair()
				pp.name = ppname
				pp.left = pl
				pp.right = pr
				pp.product_size = productsize
				# get score below
				diffmerge = merge_dict(pl.difsitedict, pr.difsitedict)
				for k, v in diffmerge.items():
					difnum = sum(i > 0 for i in v) # how many sequences this difpos can differ
					pp.score += (difnum / len(ids) + 1) * 100 / k
				if pp.score >= primer_pair_score_threshold:
					primerpairs[ppname] = pp
					line1 = "SEQUENCE_ID=" + ppname
					line2 = "PRIMER_TASK=check_primers"
					line3 = "PRIMER_PRODUCT_OPT_SIZE=100"
					line4 = "PRIMER_PRODUCT_SIZE_RANGE=50-150"
					line5 = "SEQUENCE_TEMPLATE=" + seqtemplate
					line6 =  "SEQUENCE_PRIMER=" + pl.seq
					line7 = "SEQUENCE_PRIMER_REVCOMP=" + ReverseComplement(pr.seq)
					line8 = "="
					p3temp.write("\n".join([line1, line2, line3, line4, line5, line6, line7, line8]) + "\n")

p3temp.close()
## use primer3 to check the primer pair quality
tempout = "temp_primer_pair_test_out_" + groupname + ".txt"
p3cmd = primer3_path + " -default_version=2 -p3_settings_file=./primer3web_v4_JZ.txt -output=" + " ".join([tempout, tempin])
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

outfile.write("index\tproduct_size\tprimerID\ttype\tstart\tend\tlength\tTm\tGCcontent\tany\t3'\thairpin\tprimer_nvar\t3'Diff\tDiffAll\tDifNumber\tprimer_score\tprimer_seq\tReverseComplement\tpenalty\tcompl_any\tcompl_end\tprimerpair_score\n")

#for pp in primerpairs:
for np, pp in primerpairs.iteritems():
	if pp.compl_any != "NA" and float(pp.compl_any) <= primer_pair_compl_any_threshold and float(pp.compl_end) <= primer_pair_compl_end_threshold:
		pl = pp.left
		pr = pp.right
		outfile.write("\t".join([pp.name, str(pp.product_size), pl.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end, str(pp.score)]) + "\n")
		outfile.write("\t".join([pp.name, str(pp.product_size), pr.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end, str(pp.score)]) + "\n")



"""
#outfile = open('selected_primers_for_Chr-B2.2.txt', 'w')
outfile.write("index\tprimerID\ttype\tstart\tend\tlength\tproduct_size\tTm\tGCcontent\tany\t3'\thairpin\tprimer_nvar\t3'Diff\tDiffAll\tprimer_seq\tReverseComplement\n")
primernumber = 0
for pl in newleftprimers:
	for pr in newrightprimers:
		alldifsite = [sum(x) for x in zip(pl.difsite, pr.difsite)]
		if min(alldifsite) > 0 and pl.end < pr.start and abs(pl.tm - pr.tm) < 1.0:
			#print primernumber
			productsize = pr.end - pl.start + 1
			#if productsize >= 50 and productsize <= 150:
			if productsize >= 50 and productsize <= 150:
				primernumber += 1
				outfile.write("\t".join([str(primernumber), pl.formatprimer()]) + "\n")
				outfile.write("\t".join([str(primernumber), pr.formatprimer()]) + "\n")
		
outfile.write("\nLeft that can differ all\n\n")

# for LEFT
for pl in alldifferenceleft:
	for pr in nodiffright:
		if pl.end < pr.start and abs(pl.tm - pr.tm) < 1.0:
			productsize = pr.end - pl.start + 1
			if productsize >= 50 and productsize <= 150:
				primernumber += 1
				outfile.write("\t".join([str(primernumber), pl.formatprimer()]) + "\n")
				outfile.write("\t".join([str(primernumber), pr.formatprimer()]) + "\n")

# For RIGHT	
outfile.write("\nRight that can differ all\n\n")

for pr in alldifferenceright:
	for pl in nodiffleft:
		if pl.end < pr.start and abs(pl.tm - pr.tm) < 1.0:
			productsize = pr.end - pl.start + 1
			if productsize >= 50 and productsize <= 150:
				primernumber += 1
				outfile.write("\t".join([str(primernumber), pl.formatprimer()]) + "\n")
				outfile.write("\t".join([str(primernumber), pr.formatprimer()]) + "\n")
"""


######### Get all primers across intron-exon border
outfile.write("\nLeft that across border\n\n")

# for LEFT
primernumber = 0
for pp in newleftprimers:
	if set((59,)) < set(range(pp.start-1, pp.end)) or set((505,)) < set(range(pp.start-1, pp.end)):
		primernumber += 1
		outfile.write("\t".join([str(primernumber), pp.formatprimer()]) + "\n")
print "Left primer that across border:", primernumber

# For RIGHT	
outfile.write("\nRight that across border\n\n")
primernumber = 0
for pr in newrightprimers:
	if set((59,)) < set(range(pp.start-1, pp.end)) or set((505,)) < set(range(pp.start-1, pp.end)):
		primernumber += 1
		outfile.write("\t".join([str(primernumber), pp.formatprimer()]) + "\n")
print "Right primer that across border:", primernumber

#########
outfile.close()

## remove all tempotary files
call("rm temp_*", shell=True)
