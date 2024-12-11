use Data::Dumper;

#my $gff="Eag.stringtie-merge.edited.gff3";
#my $gff="Eag.draft-genome.20160925-2.gene.id_edited.cds_edited.removeEugag004566s00025182.gff3";
#my $gff="Eag.draft-genome.20180531.gene.id_edited.cds_edited.gff3";
#my $fasta="Eag.draft-genome.20160925.dna.fa";

###2021.0901 Junesk9
#my $gff="Eag_consensusScaffold.1kb.rename.gff3";
#my $fasta="Eag_consensusScaffold.1kb.fa";
#2021.0904 Junesk9
#my $gff=" Eag.ref.gff3";
#my $fasta="Eag.ref.fa";
my $fasta="Eg-asm.chen2024.fa";
#my $gff="gmap.Eag_old-r2_comp.gff";
my $gff="Euglena.gene.change.chr.gff";

#my $bedtools_path="/data/inoue/tools/bedtools2/bin/bedtools";
my $bedtools_path="~/anaconda3/bin/bedtools";


open(GFF,$gff);

my %hash=();
while(<GFF>){
	s/\r//;
        chomp;
        my @tmp=split(/\t/,$_);
#	if(($_!~/^#/)and($tmp[2]eq "exon")){
	if(($_!~/^#/)and($tmp[2]eq "exon")){
		my @tmp2=split(";",$tmp[8]);
                my $parent="";
		my $id="";
		foreach my $each3(@tmp2){
			if($each3=~/^Parent=(\S+)/){
				$parent=$1
			}elsif($each3=~/^ID=(\S+)/){
				$id=$1;
			}else{
				
			}
		}
		if(($parent eq "")and($id eq "")){
			die();
		}
		push(@{$hash{$parent}},"$tmp[3] $tmp[4]	$tmp[6] $tmp[0] $id");
                	
	}
}
close(GFF);
#print Dumper %hash;
#die();	

my %intron_cds=();
my $intron_bed=();
foreach my $each2(keys(%hash)){
	@{$hash{$each2}}=sort{(split(/\s/,$a))[0] <=> (split(/\s/,$b))[0] }@{$hash{$each2}};
	for(my $n=0;$n<@{$hash{$each2}};$n++){
		if($n>=1){
			my $m=$n-1;
			my @mae=split(/\s/,${$hash{$each2}}[$m]);
			my @ato=split(/\s/,${$hash{$each2}}[$n]);
			
			#intron
			my $intron_start=$mae[1]+1;
			my $intron_end=$ato[0]-1;
			my $intron_length=$intron_end-$intron_start+1;
			
			#exon
			my $left_exon_end=$mae[1];
			my $right_exon_start=$ato[0];

			#seqeunce_area
			my $a_from=$mae[1]-10;
			my $a_to=$mae[1]+10;
			my $a_direc=$mae[2];
			my $b_from=$ato[0]-10;
			my $b_to=$ato[0]+10;
			my $b_direc=$ato[2];

			#print "$mae[3]	$a_from $mae[1] ($intron_start) $a_to $a_direc .. $b_from $ato[0] ($intron_end) $b_to $b_direc $intron_length\n";
			
			#bed
			my $bed_a_from=$a_from-1;
                        my $bed_b_from=$b_from-1;
                        my $bed_a_to=$a_to;
                        my $bed_b_to=$b_to;
			my $bed_intron_from=$intron_start-1;
			my $bed_intron_to=$intron_end;                        
			
			my $intron_id="$mae[3].$intron_start-$intron_end.$mae[2]";
	
			
			if($a_direc eq "+"){
                                #print "-10bp|+10bp donner site BED (eeeeeeeeeeEDDiiiiiiii) [+]\n";
				#print "[DONER_BED]	$mae[3]	$bed_a_from	$bed_a_to	$mae[4].$a_direc.$mae[1].don	.	$a_direc\n";
                                #print "-10bp|+10bp accepter site BED (iiiiiiiiAAEeeeeeeeeee) [+]\n";
				#print "[ACCEPTER_BED]	$ato[3]	$bed_b_from	$bed_b_to	$ato[4].$b_direc.$ato[0].acc	.	$b_direc\n";
                        	#print "intron fromto [+]\n";
				#print "[INTRON_BED]	$mae[3]	$bed_intron_from	$bed_intron_to	$mae[4].$intron_start-$intron_end.$ato[4].$a_direc	intron	.	$a_direc\n";
				push(@{$intron_cds{$intron_id}},"$mae[4].$intron_start-$intron_end.$ato[4].$a_direc");			
				${$intron_bed{$intron_id}}{"doner"}="$mae[3]	$bed_a_from	$bed_a_to	$intron_id.don	.	$a_direc";
				${$intron_bed{$intron_id}}{"accepter"}="$ato[3]	$bed_b_from	$bed_b_to	$intron_id.acc	.	$b_direc";
			}elsif($a_direc eq "-"){
                                #print "-10bp|+10bp donner site BED (eeeeeeeeeeEDDiiiiiiii) [-]\n";
				#print "[ACCEPTER_BED]	$mae[3]	$bed_a_from	$bed_a_to	$mae[4].$a_direc.$mae[1].acc	.	$a_direc\n";
                                #print "-10bp|+10bp accepter site BED (iiiiiiiAAEeeeeeeeeee) [+]\n";
				#print "[DONER_BED]	$ato[3]	$bed_b_from	$bed_b_to	$ato[4].$b_direc.$ato[0].don	.	$b_direc\n";
				#print "intron_fromto [-]\n";
                        	#print "[INTRON_BED]	$mae[3]	$bed_intron_from	$bed_intron_to	$mae[4].$intron_start-$intron_end.$ato[4].$a_direc	intron	.	$a_direc\n";
				push(@{$intron_cds{$intron_id}},"$mae[4].$intron_start-$intron_end.$ato[4].$a_direc");
				${$intron_bed{$intron_id}}{"accepter"}="$mae[3]	$bed_a_from	$bed_a_to	$intron_id.acc	.	$a_direc";		
				${$intron_bed{$intron_id}}{"doner"}="$ato[3]	$bed_b_from	$bed_b_to	$intron_id.don	.	$b_direc";
			}
			#print "\n";
			
		}
	}
}
close(GFF);

#print Dumper %intron_cds;
#print Dumper %intron_bed;

print "Annotation gff:	$gff\n";
open(BED,"+>$gff.intron.bed");
foreach my $intron_id(sort(keys(%intron_cds))){	
	my $don_bed=${$intron_bed{$intron_id}}{"doner"};
	my $acc_bed=${$intron_bed{$intron_id}}{"accepter"};
	print BED "#ID	$intron_id\n";
	print BED "#Exon	@{$intron_cds{$intron_id}}\n";
	print BED "$don_bed\n";
	print BED "$acc_bed\n";	
}
close(BED);
print "Intron BED:	$gff.intron.bed (-10bps|+10bps)\n";

###21.0901 "-nameOnly" for further parsing after bedtools v2.30.0 
#system("$bedtools_path getfasta -tab -name -s -fi $fasta -bed $gff.intron.bed -fo $gff.intron.seq");
system("$bedtools_path getfasta -tab -nameOnly -s -fi $fasta -bed $gff.intron.bed -fo $gff.intron.seq");
print "Intron SEQ:	$gff.intron.seq\n";



