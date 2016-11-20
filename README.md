# getprimer
Python script to design genome specific or homeolog specific primers using Primer3

For example:

I have a gene that has both homeologs and paralogs in wheat. Let's say the gene is CG. It has 3 paralogs in each wheat genome, and 3 homeologs across genomes:
A genome: CG-A1, CG-A2, CG-A3
B genome: CG-B1, CG-B2, CG-B3
D genome: CG-D1, CG-D2, CG-D3.

I want primers to amplify homeologs but not paralogs, e.g. one pair of primers to only amplify CG-A1, CG-B1 and CG-D1, one pair of primers to only amplify CG-A2, CG-B2 and CG-D2, and one pair of primers to only amplify CG-A3, CG-B3 and CG-D3.

This python script with Primer3 should be able to design that.

# Dependencies

GetPrimer use "muscle" to create multiple sequence alignment and "Primer3" to design primers.

1. Muscle: Multiple sequence alignment program (http://www.drive5.com/muscle/).
2. Primer3: program for designing PCR primers (http://primer3.sourceforge.net/).

# How to use it
1. Download the full depository to your computer (it has 64bit muscle and primer3_core working in Ubuntu 16.04).
2. Right now you need to make a fasta file that include all the similar sequences (homeologs and paralogs or any other sequences) for an input file.
3. Edit the first part of the file to change some parameters.
4. Run getprimer.py

# To-do
1. Give scores for each pair of primers for easy selection
2. Clean code

# Credits
I borrowed ideas from GSP (https://github.com/bioinfogenome/GSP), a great tool for Genome Specific Primers design in polyploid species. I only rewrote it with python (not complete copy of GSP) with my own needs. Thanks to the author of GSP.

I also borrowed some codes from biopython (https://github.com/biopython/biopython/blob/master/Bio/Emboss/Primer3.py). Thanks to them too.