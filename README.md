# getprimer
Python script to design genome specific or homeolog specific primers using Primer3

For example:

I have a gene that has both homeologs and paralogs in wheat. Let's say the gene is CG. It has 3 paralogs in each wheat genome, and 3 homeologs across genomes:  
A genome: CG-A1, CG-A2, CG-A3  
B genome: CG-B1, CG-B2, CG-B3  
D genome: CG-D1, CG-D2, CG-D3

I want primers to amplify homeologs but not paralogs, e.g. one pair of primers to only amplify CG-A1, CG-B1 and CG-D1, one pair of primers to only amplify CG-A2, CG-B2 and CG-D2, and one pair of primers to only amplify CG-A3, CG-B3 and CG-D3.

This python script with Primer3 can design primers that can only amplify one group of sequences but not the others.

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
	-i <sequence.fa>
	-p <GetPrimer path>
	-s <product min size>
	-l <product max size>
	-g <target sequences: seqID1,seqID2>
	-r <range limit: where the differences should be in the primer from 3' end, default 10>
	-o <output file name>"
	-v <primer overlap region (such as intron): n1-n2,n3-n4>
  -h help
```
1. Download the full depository to your computer (it has 64bit muscle and primer3_core working in Ubuntu 16.04).
2. Recompile muscle or primer3_core from source if they do not work in your computer. Put them in the bin subfolder.
3. You should be able to use it.


# To-do
1. Clean and optimize code


# Credits
I borrowed ideas from GSP (https://github.com/bioinfogenome/GSP), a great tool for Genome Specific Primers design in polyploid species. I only rewrote it with python (not complete copy of GSP) with my own needs. Thanks to the author of GSP.

I also borrowed some codes from biopython (https://github.com/biopython/biopython/blob/master/Bio/Emboss/Primer3.py). Thanks to them too.
