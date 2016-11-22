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

### Imported
from subprocess import call
import getopt, sys, os

usage="""
getprimer.py 
	-i <sequence.fa>
	-p <GetPrimer path> 
	-s <product min size> 
	-l <product max size> 
	-o <output file name>"
"""


#########################
# parameter or file names that need to be changed

seqfile = "sequence.fa" # CHANGE HERE
#targets = ["Chr-D2.2", "Chr-A2.2"] # CHAGE HERE
primer_pair_score_threshold = 100.0
primer_pair_compl_any_threshold = 1.0
primer_pair_compl_end_threshold = 1.0
product_min = 50
product_opt = 100
product_max = 150
getprimer_path = "~/Research/Software/git/github/getprimer"
out = ""
msa = 1 # whether need to do multiple sequence alignment
overlap_region = range(60,159+1) + range(605,907+1) # intron region

# read command line options
print "Parsing command line optinos"

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:p:s:l:r:o:m:h", ["help"])
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
		primer3_path = a
	elif o in ("-s"):
		product_min = int(a)
	elif o in ("-l"):
		product_max = int(a)
	elif o in ("-o"):
		out = a
	elif o in ("-m"):
		msa = int(a)
	else:
		assert False, "unhandled option"

print "OPtions done"

getprimer_path = os.path.expanduser(getprimer_path)
muscle_path = getprimer_path + "/bin/muscle"
primer3_path = getprimer_path + "/bin/primer3_core"
primer3_parameter_path = getprimer_path + "/primer3web_v4_JZ.txt"
# other variables

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
	ids.append(kk)

print "Sequence IDs: ", ids

#groupname = "-".join(ids)

if not out:
	out = "selected_common_primers.txt"
outfile = open(out, 'w') # output file
	
alignlen = len(fasta[ids[0]])
print "Alignment length: ", alignlen

varexclude = [] # varaition among targets, any primers that have these should be excluded. 
ngap = 0 # gaps
mainID = ids[0] # the one whose primer3 output will be used
print "Sequence used for creating primers:", mainID

for i in range(alignlen):
	nd = 0 # number of difference
	dd = [] # to record site difference for each seq
	if fasta[mainID][i] == "-":
		ngap += 1
	for j in ids:#check difference within targets
		if fasta[j][i] == fasta[mainID][i]: nd+=1
	if nd < len(ids): # if the site i has no differences within targets
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

# I only need to use the output from the first target
line1 = "SEQUENCE_ID=" + mainID
line2 = "SEQUENCE_TEMPLATE=" + fasta[mainID].replace("-","") # remove "-" in the alignment seq
#line3 = "PRIMER_PRODUCT_OPT_SIZE=" + str(product_opt)
line3 = "PRIMER_PRODUCT_SIZE_RANGE=" + str(product_min) + "-" + str(product_max)
line4 = "="
p3input.write("\n".join([line1, line2, line3, line4]) + "\n")

p3input.close()

# primer3 output file
primer3output = "primer3.output"
p3cmd = primer3_path + " -default_version=2 -format_output -p3_settings_file=" + primer3_parameter_path + " -output=" + primer3output + " " + primer3input
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

	def formatprimer(self):
		if self.direction == "LEFT_PRIMER":
			formatout = "\t".join(str(x) for x in [self.name, self.direction, self.start, self.end, self.length, self.tm, self.gc, self.anys, self.three, self.hairpin, self.seq, ""])
		else:
			formatout = "\t".join(str(x) for x in [self.name, self.direction, self.start, self.end, self.length, self.tm, self.gc, self.anys, self.three, self.hairpin, ReverseComplement(self.seq), self.seq])
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


# parse primer3 output file
# use 0-based primer3 output
def parse(handle, direction):
	"""Iterate over primer3 output
	"""
	primerlist = []
	with open(handle) as infile:
		startread = 0
		for line in infile:
			if mainID in line:
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
			if startread and "PRIMER PICKING RESULTS FOR" in line and mainID not in line:
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

for pp in leftprimers:
	exclude = 0 # a tag to see whether the primer has difference within target group
	for var in varexclude: # check whether this primer has sites different within targets
		if set(var) < set(range(pp.start-1, pp.end)):
			exclude += 1
	if exclude == 0: # only if the primer does not have any exlude sites
		newleftprimers.append(pp)

print "Number of selected LEFT primers", len(newleftprimers)

if len(newleftprimers) == 0:
	print "No COMMON left primers found!"
	sys.exit(1)

# selected right primers
newrightprimers = []


for pp in rightprimers:
	exclude = 0
	for var in varexclude: # check whether this primer has sites different within targets
		if set(var) < set(range(pp.start-1, pp.end)):
			exclude += 1
	if exclude == 0: # only if the primer does not have any exlude sites
		newrightprimers.append(pp)

print "Number of selected RIGHT primers", len(newrightprimers)

if len(newrightprimers) == 0:
	print "No COMMON RIGHT primers found!"
	sys.exit(1)
	
#############################
# STEP 5: Select and Test Primer Pairs

# selected primers pairs
tempin = "temp_primer_pair_check_input.txt"
p3temp = open(tempin, 'w') # output file
seqtemplate = fasta[mainID].replace("-","") # remove "-" in the alignment seq
primernumber = 0
primerpairs = {} # all the pairs with the right size and Tm differences


for pl in newleftprimers:
	for pr in newrightprimers:
		productsize = pr.end - pl.start + 1
		if productsize >= product_min and productsize <= product_max:
			primernumber += 1
			ppname = str(primernumber)
			pp = PrimerPair()
			pp.name = ppname
			pp.left = pl
			pp.right = pr
			pp.product_size = productsize
			primerpairs[ppname] = pp
			line1 = "SEQUENCE_ID=" + ppname
			line2 = "PRIMER_TASK=check_primers"
			#line3 = "PRIMER_PRODUCT_OPT_SIZE=100" + str(product_opt)
			line4 = "PRIMER_PRODUCT_SIZE_RANGE=" + str(product_min) + "-" + str(product_max)
			line5 = "SEQUENCE_TEMPLATE=" + seqtemplate
			line6 =  "SEQUENCE_PRIMER=" + pl.seq
			line7 = "SEQUENCE_PRIMER_REVCOMP=" + ReverseComplement(pr.seq)
			line8 = "="
			p3temp.write("\n".join([line1, line2, line4, line5, line6, line7, line8]) + "\n")

p3temp.close()

if primernumber == 0:
	print "No COMMON primer pairs found!"
	sys.exit(1)

## use primer3 to check the primer pair quality
tempout = "temp_primer_pair_test_out.txt"
p3cmd = primer3_path + " -default_version=2 -p3_settings_file=" + primer3_parameter_path + " -output=" + " ".join([tempout, tempin])
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
		outfile.write("\t".join([pp.name, str(pp.product_size), pl.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end]) + "\n")
		outfile.write("\t".join([pp.name, str(pp.product_size), pr.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end]) + "\n")



######### Get all primers across intron-exon border
outfile.write("\nLeft that across border\n\n")

# for LEFT
primernumber = 0
for pp in newleftprimers:
	if set(overlap_region)|set(range(pp.start, pp.end-1)):
		primernumber += 1
		outfile.write("\t".join([str(primernumber), pp.formatprimer()]) + "\n")
print "Left primer that across border:", primernumber

# For RIGHT	
outfile.write("\nRight that across border\n\n")
primernumber = 0
for pr in newrightprimers:
	if set(overlap_region)|set(range(pp.start, pp.end-1)):
		primernumber += 1
		outfile.write("\t".join([str(primernumber), pp.formatprimer()]) + "\n")
print "Right primer that across border:", primernumber

#########
outfile.close()

## remove all tempotary files
#call("rm temp_*", shell=True)
