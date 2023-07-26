use Data::Dumper;
use strict;

my $gff="Eag.draft-genome.20180531.gene.id_edited.cds_edited.gt.gff3";
my $fasta="Eag.draft-genome.20160925.dna.fa";
my $bedtools_path="/usr/local/bio/bedtools2-2.26.0/bin/bedtools";

my $exon_overhang=11;
my $intron_overhang=20;
my $slidable_range=5; 
my $N_range=12;

if(($N_range+$slidable_range)>$intron_overhang){
	print "ERROR	intron length error\n";
	die();
}
		



print "#GFF:	$gff\n";
print "#FASTA:	$fasta\n";
print "#length of exonic regions:	$exon_overhang\n";
print "#length of intronic regions	$intron_overhang\n";
print "#length of slidable range checked	$slidable_range\n";
print "#length of introns containing 'N' checked	$N_range\n";



=pod
my ($intron_bed,$intron_parents)=get_intron_bed($gff,$exon_overhang,$intron_overhang);
#print Dumper %intron_bed;
open(BED,"+>$gff.intron_boundary.bed");
print BED "#BED file for $gff.intron_boundary, [$exon_overhang bps]/[$intron_overhang bps]......[$intron_overhang bps]/[$exon_overhang bps]\n";
foreach my $intron_id(sort(keys(%$intron_bed))){
	my $don_bed=${$$intron_bed{$intron_id}}{"doner"};
        my $acc_bed=${$$intron_bed{$intron_id}}{"accepter"};
       	my $don_exon_bed=${$$intron_bed{$intron_id}}{"doner_exon"};
	my $acc_exon_bed=${$$intron_bed{$intron_id}}{"accepter_exon"};	
	print BED "#ID	$intron_id\n";
        print BED "#mRNAs	@{$$intron_parents{$intron_id}}\n";
       	print BED "$don_bed\n";
	print BED "$don_exon_bed\n";
        print BED "$acc_bed\n";
	print BED "$acc_exon_bed\n";	
}
close(BED);
print "$gff.intron_boundary.bed was produced.\n";
system("$bedtools_path getfasta -name -tab -s -fi $fasta -bed $gff.intron_boundary.bed -fo $gff.intron_boundary.seq");
print "$gff.intron_boundary.seq was produced.\n";
=cut

my %intron_seq=get_intron_seq("$gff.intron_boundary.seq");
#print Dumper %intron_seq;
my %intron_slidability=check_slidablity_intron_boundary(\%intron_seq);
#print Dumper %intron_slidability;
open(SLD,"+>$gff.intron_boundary.slided");
foreach my $each(sort(keys(%intron_slidability))){
	print SLD "$each\n";
	foreach my $each2(@{$intron_slidability{$each}}){
		print SLD "$each2\n";
	}
	print SLD "//\n";
}
close(SLD);
print "$gff.intron_boundary.slided was generated\n";

my %intron_ambiguity=check_ambiguity_intron_boundary(\%intron_slidability, $N_range);
#print Dumper %intron_ambiguity;
open(AMB,"+>$gff.intron_boundary.ambiguous");
foreach my $each(sort(keys(%intron_ambiguity))){
        print AMB "$each\n";
        foreach my $each2(@{$intron_ambiguity{$each}}){
                print AMB "$each2\n";
        }
        print AMB "//\n";
}
close(AMB);
print "$gff.intron_boundary.ambiguity was generated\n";

my %intron_class=check_intron_class(\%intron_ambiguity);


sub check_intron_class{
	my $ref=$_[0];
	my %intron_ambiguity=%$ref;
	my %slidablity_introns=();	
	
	foreach my $each(sort(keys(%intron_ambiguity))){
		my $count_slide=0;
		foreach my $each2(@{$intron_ambiguity{$each}}){
			my @tmp=split(/\t/,$each2);
			if(($tmp[3]=~/\d+/)and($tmp[3]==1)){
				$count_slide++;
				#print "$each2	$tmp[3]\n";
			}
		}
		if($count_slide ==0){
			push(@{$slidablity_introns{"non_slidable"}},$each);
			#print Dumper @{$intron_ambiguity{$each}};
		}elsif($count_slide >=1){
			push(@{$slidablity_introns{"slidable"}},$each);
			#print Dumper @{$intron_ambiguity{$each}};
		}else{
			die();
		}
	}
	my $num_slidable=scalar(@{$slidablity_introns{"slidable"}});
	my $num_nonslidable=scalar(@{$slidablity_introns{"non_slidable"}});
	print "$num_slidable slidable inrons were identified\n";
	print "$num_nonslidable non_slidable inrons were identified\n";
				
				
}





sub check_ambiguity_intron_boundary{
	my $ref=$_[0];
	my %intron_slidability=%$ref;
	my $length_N=$_[1];
	my %intron_ambiguity=();
	foreach my $each(sort(keys(%intron_slidability))){
		#print "$each\n";
		foreach my $each2(@{$intron_slidability{$each}}){
			my @tmp=split(/\t/,$each2);
			if($tmp[1]=~/\S+\/(\S{$length_N}).*\s\.\.\s.*(\S{$length_N})\/\S+/){
				my $n_flag=0;
				my $intron_don=$1;
				my $intron_acc=$2;
				if(($intron_don=~/N/)or($intron_acc=~/N/)){
					$n_flag++;
				}
				#print "$each2	$n_flag\n";
				push(@{$intron_ambiguity{$each}},"$each2	$n_flag");
			}else{
				print "ERROR	intron boundary format error\n";
				die();	
			}	
		
		}
	}
	return(%intron_ambiguity);
}

sub check_slidablity_intron_boundary{
	my $ref=$_[0];
	my %intron_seq=%$ref;
	my $num_introns=scalar(keys(%intron_seq));
	print "$num_introns introns and thier boundary sequences were identified.\n";
	my %intron_slidability=();
	foreach my $intron_id(sort(keys(%intron_seq))){
		my $don_intron_seq=${$intron_seq{$intron_id}}{"don_intron"};
                my $acc_intron_seq=${$intron_seq{$intron_id}}{"acc_intron"};
		my $don_exon_seq=${$intron_seq{$intron_id}}{"don_exon"};
		my $acc_exon_seq=${$intron_seq{$intron_id}}{"acc_exon"};

		#print "#ID	$intron_id\n";
		#print "$don_exon_seq/$don_intron_seq .. $acc_intron_seq/$acc_exon_seq\n";
		my @slided=check_slidability($don_intron_seq,$acc_intron_seq,$don_exon_seq,$acc_exon_seq);
		@{$intron_slidability{$intron_id}}=@slided;
	}
	return(%intron_slidability);
}

sub check_slidability{
	my $don_intron_seq=$_[0];
	my $acc_intron_seq=$_[1];
	my $don_exon_seq=$_[2];
	my $acc_exon_seq=$_[3];

	my %hash=();
	#print "0	$don_exon_seq/$don_intron_seq .. $acc_intron_seq/$acc_exon_seq\n";
	my $don=substr($don_intron_seq,0,2);
        my $acc=substr($acc_intron_seq,-2,2);
	my @slided=();
	push(@slided,"0	$don_exon_seq/$don_intron_seq .. $acc_intron_seq/$acc_exon_seq	$don $acc	_");
	for(my $n=1;$n<=$slidable_range;$n++){
		my $don_shifted=substr($don_intron_seq,0,$n);
		my $don_exon_seq_shift=$don_exon_seq.$don_shifted;
		my $don_intron_seq_shift=substr($don_intron_seq,$n,length($don_intron_seq)-$n);
			
		my $acc_shifted=substr($acc_exon_seq,0,$n);
		my $acc_intron_seq_shift=$acc_intron_seq.$acc_shifted;
		my $acc_exon_seq_shift=substr($acc_exon_seq,$n,length($acc_exon_seq)-$n);
		my $slide_flag=0;
		if($don_shifted eq $acc_shifted){
			$slide_flag=1;
		}
		my $shifted_don=substr($don_intron_seq_shift,0,2);
		my $shifted_acc=substr($acc_intron_seq_shift,-2,2);
		#print "+$n	$don_exon_seq_shift/$don_intron_seq_shift .. $acc_intron_seq_shift/$acc_exon_seq_shift	$shifted_don $shifted_acc	$slide_flag\n";
		push(@slided,"+$n	$don_exon_seq_shift/$don_intron_seq_shift .. $acc_intron_seq_shift/$acc_exon_seq_shift	$shifted_don $shifted_acc	$slide_flag");
	}
	
	for(my $n=1;$n<=$slidable_range;$n++){
                my $don_shifted=substr($don_exon_seq,-1*$n);
                my $don_intron_seq_shift=$don_shifted.$don_intron_seq;
                my $don_exon_seq_shift=substr($don_exon_seq,0,length($don_exon_seq)-$n);

                my $acc_shifted=substr($acc_intron_seq,-1*$n);
                my $acc_exon_seq_shift=$acc_shifted.$acc_exon_seq;
                my $acc_intron_seq_shift=substr($acc_intron_seq,0,length($acc_intron_seq)-$n);
                my $slide_flag=0;
                if($don_shifted eq $acc_shifted){
                        $slide_flag=1;
                }
                my $shifted_don=substr($don_intron_seq_shift,0,2);
                my $shifted_acc=substr($acc_intron_seq_shift,-2,2);

		#print "-$n	$don_exon_seq_shift/$don_intron_seq_shift .. $acc_intron_seq_shift/$acc_exon_seq_shift	$shifted_don $shifted_acc	$slide_flag\n";
        	push(@slided,"-$n	$don_exon_seq_shift/$don_intron_seq_shift .. $acc_intron_seq_shift/$acc_exon_seq_shift	$shifted_don $shifted_acc	$slide_flag")
	}
	return(@slided);	
}

sub get_intron_seq{
	my $seq=$_[0];
	my %hash=();
	open(SEQ,$seq);
	while(<SEQ>){
        	s/\r//;
        	chomp;
       		my @tmp=split(/\t/,$_);
        	my $id=$tmp[0];
		$id=~s/::\S+$//;
        	if($id=~/^(\S+)\.(don_intron|acc_intron|don_exon|acc_exon)$/){
                	my $id2=$1;
	                my $type=$2;
        	        #print "$id2 $type $tmp[1]\n";
			${$hash{$id2}}{$type}=$tmp[1];
        	}
	}
	close(SEQ);	
	return(%hash);
}

sub get_intron_bed{
	my $gff=$_[0];
	my %intron_parents=();
	my $exon_overhang=$_[1];
	my $intron_overhang=$_[2];

	open(GFF,$gff);
	while(<GFF>){
		if($_!~/^#/){
			s/\r//;
			chomp;
			my @tmp=split(/\t/,$_);
			if($tmp[2]eq"intron"){
				my $intron_id="$tmp[0].$tmp[3]-$tmp[4].$tmp[6]";
				my $parent=$tmp[8];
				$parent=~s/^Parent=//;
				push(@{$intron_parents{$intron_id}},$parent);
			}
		}
	}
	close(GFF);
	#print Dumper %intron_parents;
	my $num_introns=scalar(keys(%intron_parents));
	print "$num_introns introns were identified in $gff\n";

	my %intron_bed=();
	foreach my $intron_id(sort(keys(%intron_parents))){
		my $chr="";
		my $intron_start="";
		my $intron_end="";
		my $direc="";
		if($intron_id=~/^(\S+)\.(\S+)-(\S+)\.([+-])/){
			$chr=$1;
			$intron_start=$2;
			$intron_end=$3;
			$direc=$4;
		}else{
			die();
		}
		#print "$intron_id -> $chr $intron_start $intron_end $direc\n";
		my $a_intron_from=$intron_start;
		my $a_intron_to=$intron_start+$intron_overhang-1;
		my $a_exon_to=$a_intron_from-1;
		my $a_exon_from=$a_exon_to-$exon_overhang+1;
		my $a_direc=$direc;

		my $b_intron_from=$intron_end-$intron_overhang+1;
		my $b_intron_to=$intron_end;
                my $b_exon_from=$b_intron_to+1;
		my $b_exon_to=$b_exon_from+$exon_overhang-1;
		my $b_direc=$direc;

		my $bed_a_intron_from=$a_intron_from-1;
                my $bed_b_intron_from=$b_intron_from-1;
                my $bed_a_exon_from=$a_exon_from-1;
                my $bed_b_exon_from=$b_exon_from-1;

		if($a_direc eq "+"){
			${$intron_bed{$intron_id}}{"doner"}="$chr	$bed_a_intron_from	$a_intron_to	$intron_id.don_intron	.	$a_direc";
                        ${$intron_bed{$intron_id}}{"accepter"}="$chr	$bed_b_intron_from	$b_intron_to	$intron_id.acc_intron	.	$b_direc";
			${$intron_bed{$intron_id}}{"doner_exon"}="$chr	$bed_a_exon_from	$a_exon_to	$intron_id.don_exon	.	$a_direc";
			${$intron_bed{$intron_id}}{"accepter_exon"}="$chr	$bed_b_exon_from	$b_exon_to	$intron_id.acc_exon	.	$b_direc";
                }elsif($a_direc eq "-"){
                        ${$intron_bed{$intron_id}}{"accepter"}="$chr	$bed_a_intron_from	$a_intron_to	$intron_id.acc_intron	.	$a_direc";           
                        ${$intron_bed{$intron_id}}{"doner"}="$chr	$bed_b_intron_from	$b_intron_to	$intron_id.don_intron	.	$b_direc";
			${$intron_bed{$intron_id}}{"accepter_exon"}="$chr	$bed_a_exon_from	$a_exon_to	$intron_id.acc_exon	.	$a_direc";   
                        ${$intron_bed{$intron_id}}{"doner_exon"}="$chr	$bed_b_exon_from	$b_exon_to	$intron_id.don_exon	.	$b_direc";

		}
	}
	return(\%intron_bed,\%intron_parents);
}
	



