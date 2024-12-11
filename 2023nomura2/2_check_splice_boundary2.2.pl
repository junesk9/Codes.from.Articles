use strict;
use Data::Dumper;

#my $fasta="Eag.draft-genome.20160925.dna.fa";
#my $bedtools_path="/data/inoue/tools/bedtools2/bin/bedtools";

#my $intron_seq="Eag.stringtie-merge.edited.gff3.intron.seq";
#my $bed="Eag.stringtie-merge.edited.gff3.intron.bed";
#my $intron_seq="Eag.draft-genome.20180531.gene.id_edited.cds_edited.gff3.intron.seq";
#my $bed="Eag.draft-genome.20180531.gene.id_edited.cds_edited.gff3.intron.bed";

#### input adjusted by Junesk9
#my $bedtools_path="~/anaconda3/bin/bedtools";
#my $fasta="Eag_consensusScaffold.1kb.rename.fa";
#my $intron_seq="Eag_consensusScaffold.1kb.rename.gff3.intron.seq";
#my $bed="Eag_consensusScaffold.1kb.rename.gff3.intron.bed";

### 21.0904 Junesk9
my $bedtools_path="~/anaconda3/bin/bedtools";
my $fasta="Eg-asm.chen2024.fa";
my $intron_seq="Euglena.gene.change.chr.gff.intron.seq";
my $bed="Euglena.gene.change.chr.gff.intron.bed";



open(SEQ,$intron_seq);
my $id="";
my %hash=();
while(<SEQ>){
	s/\r//;
	chomp;
	my @tmp=split(/\t/,$_);
	my $id=$tmp[0];
	if($id=~/^(\S+)\.(\S\S\S)\(\S\)$/){
		my $id2=$1;
		my $type=$2;
		${$hash{$id2}}{$type}=$tmp[1];
	}
}
close(SEQ);
#print Dumper %hash;

	
open(BED,$bed);
my %hash2=();
my $comment2="";

my %intron_class=();
my %intron_seq=();
my %intron_ddaa=();

while(<BED>){
	if(/^#ID	(\S+)/){
		s/\r//;
		chomp;
		print "$_\n";
		my $id=$1;
		my $don_seq=${$hash{$id}}{"don"};
		my $acc_seq=${$hash{$id}}{"acc"};
		print "eeeeeeeeeeEDDIiiiiiii .. iiiiiiiIAAEeeeeeeeeee\n";
		print "$don_seq .. $acc_seq\n"; 
		my ($comment1,$ddaa)=check_intron_boundary($don_seq,$acc_seq,$id);
		print "#Comment1	$comment1\n";
		$comment2="";
		if($comment1=~/^Non-conventional/){
			if($ddaa=~/AT-A./){
				print "#Comment2	Non-conventional(U12-type)\n";
				$comment2="Non-conventional(U12-type)";
			}elsif(($ddaa=~/GT-/)or($ddaa=~/-AG/)){
				print "#Comment2	Non-conventional(Intermediate)\n";
				$comment2="Non-conventional(Intermediate)";
			}else{
				print "#Comment2	Non-conventional(Unclassified)\n";
				$comment2="Non-conventional(Unclassified)";
			}
		}else{
			print "#Comment2	$comment1\n";
			$comment2=$comment1;	
		}	
		$hash2{$comment2}++;
		push(@{$intron_class{$id}},$comment2);		
		push(@{$intron_ddaa{$id}},$ddaa);
		push(@{$intron_seq{$id}},"$don_seq .. $acc_seq");
		print "//\n\n";
	}
}
close(BED);

print Dumper %hash2;

my $out=$intron_seq;
$out=~s/\.seq$/\.filtered/;
open(OUT,"+>$out");

foreach my $each(sort(keys(%intron_class))){
	print OUT "$each	@{$intron_class{$each}}	${$intron_ddaa{$each}}[0]	${$intron_seq{$each}}[0]\n";
}
close(OUT);




sub check_intron_boundary{
	my $don_seq=$_[0];
	my $acc_seq=$_[1];
	my $id=$_[2];

	my @don=split("",$don_seq);
	my @acc=split("",$acc_seq);
	
	my $comment="";
	my %hash=();
	

	$hash{"DD-AA"}="$don[11]$don[12]-$acc[8]$acc[9]";
	print "$don[11]$don[12]-$acc[8]$acc[9]\n";
	$hash{"don_dinuc"}="$don[11]$don[12]";
	$hash{"acc_dinuc"}="$acc[8]$acc[9]";


	if($hash{"DD-AA"}=~/N/){
		print "$hash{'DD-AA'} -> The inton include ambiguous nucleotide(s)in its terminal dinucleotides.\n";
		$comment="Ambiguous Nucleotides";
	}else{
		if($hash{'DD-AA'}eq'GT-AG'){
			print "$hash{'DD-AA'} -> The intron is possessing the canonical GT-AG splice sites rocessed U2 snRNA-containing spliceosome.\n";
			$comment="Conventional (GT-AG)";
		}elsif($hash{'DD-AA'}eq'GC-AG'){
			print "$hash{'DD-AA'} -> The intron is possessing the canonical GC-AG splice sites rocessed U2 snRNA-containing spliceosome.\n";
			$comment="Conventional (GC-AG)";
		}else{
			my $seq="$don_seq"."xxxx"."$acc_seq";
			
			print "$hash{'DD-AA'} -> The intron is possessing the non-canonical splice sites (No GT/C-AG found in the 10 bps franking)\n";
			$comment="Non-conventional(No GT/C-AG found)";
			
		}
	}
	return($comment,$hash{'DD-AA'});	
	
}


sub get_intron_boundary{
	my $id=$_[0];
	my $scaffold="";
	my $direc="";
	my $intron_start="";
	my $intron_end="";
	if($id=~/^(Scaffold\S+)\.(\d+)\-(\d+)\.([+-])/){
	#if($id=~/^(scaffold\S+)\.(\d+)\-(\d+)\.([+-])/){
		$scaffold=$1;
		$intron_start=$2;
		$intron_end=$3;
		$direc=$4;
	}	

	my $boundary_seq="";
	if($direc eq "+"){
		my $don_start=$intron_start-10-1;
		my $don_end=$intron_start+10-1;	
		my $acc_start=$intron_end-10+1;
		my $acc_end=$intron_end+10+1;

		my $don_start_bed=$don_start-1;		
		my $acc_start_bed=$acc_start-1;
		
		my $don_bed="$scaffold	$don_start_bed	$don_end	$id.don	.	$direc";
		my $acc_bed="$scaffold	$acc_start_bed	$acc_end	$id.acc	.	$direc";
		print "$don_bed\n";
		print "$acc_bed\n";

		my ($don_seq,$acc_seq)=call_bed_get_fasta($don_bed,$acc_bed);
		return($don_seq,$acc_seq);
	}elsif($direc eq "-"){
		my $don_start=$intron_end-10+1;
		my $don_end=$intron_end+10+1;
		my $acc_start=$intron_start-10-1;
		my $acc_end=$intron_start+10-1;
		
		my $don_start_bed=$don_start-1;
                my $acc_start_bed=$acc_start-1;

		my $don_bed="$scaffold	$don_start_bed	$don_end	$id.don	.	$direc";
		my $acc_bed="$scaffold	$acc_start_bed	$acc_end	$id.acc	.	$direc";	
		print "$don_bed\n";
		print "$acc_bed\n";		

		my ($don_seq,$acc_seq)=call_bed_get_fasta($don_bed,$acc_bed);		
		return($don_seq,$acc_seq);
	}

}


sub call_bed_get_fasta{
	my $don_bed=$_[0];
	my $acc_bed=$_[1];
	open(OUT,"+>tmp.bed");	
	print OUT "$don_bed\n";
	print OUT "$acc_bed\n";
	close(OUT);

	system("$bedtools_path getfasta -tab -name -s -fi $fasta -bed tmp.bed -fo tmp.intron.seq");	
	
	my $don_seq="";
	my $acc_seq="";

	open(SEQ,"tmp.intron.seq");
	while(<SEQ>){
		s/\r//;
		chomp;
		my @tmp=split(/\t/,$_);
		if($tmp[0]=~/\.don/){
			$don_seq=$tmp[1];
		}elsif($tmp[0]=~/\.acc/){
			$acc_seq=$tmp[1];
		}else{
			die();
		}
	}
	close(SEQ);
	return($don_seq,$acc_seq);	

}

