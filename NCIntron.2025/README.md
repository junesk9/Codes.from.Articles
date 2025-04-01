***
### "The noncanonical splicing code in Euglena"
T Nomura et al. (in prep.)

1. ***[1_check_splice_boundary1.pl](https://github.com/junesk9/Codes.from.Articles/blob/main/2023nomura2/1_check_splice_boundary1.pl)*** - Perl script to extract the intron information from inputs of the genome fasta and the GFF file.
2. ***[2_check_splice_boundary2.2.pl](https://github.com/junesk9/Codes.from.Articles/blob/main/2023nomura2/2_check_splice_boundary2.2.pl)*** - Perl script to parse the #1 output and check the sequence integrity of the exon-intron boarder. STDOUT for next process.
3. ***[3_check_splice_boundary3.pl](https://github.com/junesk9/Codes.from.Articles/blob/main/2023nomura2/3_check_splice_boundary3.pl)*** - Perl script to parse the #2 output. STDOUT for next process.
4. ***[4_check_boundary_shift.pl](https://github.com/junesk9/Codes.from.Articles/blob/main/2023nomura2/4_check_boundary_shift.pl)*** - Perl script to analyze the #3 output to check the slidalbility of each exon-intron boarder. STDOUT for next process.
5. ***[5_collect.nonSlidables.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023nomura2/5_collect.nonSlidables.py)*** - Python3 script to collect the non-slidable introns from the #4 output. 

