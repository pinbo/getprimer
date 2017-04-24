# GetPrimer
Python script to design genome or group specific primers using Muscle and Primer3.

For example:

I have a gene that has both homeologs and paralogs in wheat. Let's say the gene is CG. It has 3 paralogs in each wheat genome, and 3 homeologs across genomes:  
A genome: CG-A1, CG-A2, CG-A3  
B genome: CG-B1, CG-B2, CG-B3  
D genome: CG-D1, CG-D2, CG-D3

I want primers to amplify homeologs but not paralogs, e.g. one pair of primers to only amplify CG-A1, CG-B1 and CG-D1, one pair of primers to only amplify CG-A2, CG-B2 and CG-D2, and one pair of primers to only amplify CG-A3, CG-B3 and CG-D3.

This python script with Primer3 can design primers that can only amplify one group of sequences but not the others.

# Changes
- 01/01/2017 Added **getprimer2.py**, which looses the primer criteria but strengths the primer specific (2-base differences in the first 4 3'-termini bases or 5-base differences in the first 15 3'termini bases)
- 01/04/2017 Added support for windows 7 users (but still need to use terminal either built-in or Cygwin)
- 04/23/2017 Added support for MacOSX system; removed the option "-p" for script path.
- 04/23/2017 Loosed the primer criteria: 1 base differences in the first 4 3'-termini bases or 4 base differences in the first 15 3'termini bases)

# Dependencies

GetPrimer use "muscle" to create multiple sequence alignment and "Primer3" to design primers.

1. Muscle: Multiple sequence alignment program (http://www.drive5.com/muscle/).
2. Primer3: program for designing PCR primers (http://primer3.sourceforge.net/).

# How it works
1. Use muscle to do multiple sequence alignment of user provided fasta sequences;
2. Find all the different sites;
3. Use Primer3 to design all possible left and right primers in one of the target sequences;
4. Only keep primers with no variations among target sequences but with variations from all other sequences. Score them based on variation positions and differentiation number.
5. Keep pairs of primers with required size and other criteria.
6. Use Primer3 again to score the selected pairs for primer dimers.

# Usage
```
getprimer.py
	-i <sequence.fa> : a fasta file with all sequences.
	-s <product min size> : default now is 50. Modify the script to change the default value.
	-l <product max size> : default now is 150.
	-g <seqID1,seqID2> : target sequences, the sequence IDs in the fasta file which you want the primers to amplify.
	-r <range limit> : where the differences should be in the primer from 3' end, 
	   default 10, because variations in the 5' do not count much for primer specification.
	-o <output file name>
	-v <n1-n2,n3-n4> : primer overlap region (such as intron or exon-exon junction)
	   to design genomic DNA specific or RNA specific primers
	-f <1 or 0, default 0: no filtering of the output>
	-a <primer_pair_compl_any_threshold: default 10>
	-e <primer_pair_compl_end_threshold: default 10>
	-c <primer_pair_score_threshold: default 50>
	-h help
```
1. Download the full depository to your computer (it has 64bit muscle and primer3_core working in Ubuntu 16.04).
2. Download and/or Recompile muscle or primer3_core from source the if they do not work in your computer. Rename them as "muscle" and "primer3_core" (without .exe) and put them in the bin sub-folder (replace the existed ones).
3. You should be able to use it.
4. I found a lot of primers are around one same sequence variaiton site in the output file. You can sort them by primer pair score to find a few good primers. I added a filter option: "-f 1" will remove similar primers in the same region.

**Example**

./getprimer2.py -i sequence.fa -s 100 -l 500 -g Chr-B2.2,Chr-B2.3 -v 59-60,300-400 -f 1 -o primers_for_B2.2_2.3.txt

# To-do
1. Clean and optimize code


# Credits
I borrowed ideas from GSP (https://github.com/bioinfogenome/GSP), a great tool for Genome Specific Primers design in polyploid species. I only rewrote it with python (not complete copy of GSP) with my own needs. Thanks to the author of GSP.

I also borrowed some codes from biopython (https://github.com/biopython/biopython/blob/master/Bio/Emboss/Primer3.py). Thanks to them too.

Thanks to the open source software **Muscle** (http://www.drive5.com/muscle/) and **Primer3** (http://primer3.sourceforge.net/).
