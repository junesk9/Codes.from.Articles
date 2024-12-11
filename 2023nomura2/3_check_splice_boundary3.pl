use Data::Dumper;

#my $filtered="Eag.stringtie-merge.edited.gff3.intron.filtered";
#my $filtered="Eag.draft-genome.20180531.gene.id_edited.cds_edited.gff3.intron.filtered";

##by Junesk9 21.0901
#my $filtered="Eag_consensusScaffold.1kb.rename.gff3.intron.filtered";
##by Junesk9 21.0904
#my $filtered="Eag.ref.gff3.intron.filtered";
my $filtered="Euglena.gene.change.chr.gff.intron.filtered";


my %don_intron=();
my %acc_intron=();

my $mae=5;
my $ato=5;

my %intron=();
my %intron2=();
open(FILTERED,$filtered);
while(<FILTERED>){
	s/\r//;
	chomp;
	my @tmp=split(/\t/,$_);
	my $intron_id=$tmp[0];
	#if($intron_id=~/(scaffold\S+)\.(\d+)\-(\d+)\.([\+\-])/){
	if($intron_id=~/(Chromosome\d+)\.(\d+)\-(\d+)\.([\+\-])/){
	#if($intron_id=~/(EugagSCF\d+)\.(\d+)\-(\d+)\.([\+\-])/){
		my $scaffold=$1;
		my $start=$2;
		my $end=$3;
		my $direc=$4;	
		#print "$intron_id	$scaffold	$start $end	$direc\n";	
		if($direc eq "+"){
			push(@{${$don_intron{$scaffold}}{"$scaffold $start $direc"}},"$intron_id");
			push(@{${$acc_intron{$scaffold}}{"$scaffold $end $direc"}},"$intron_id");
		}elsif($direc eq "-"){
			push(@{${$acc_intron{$scaffold}}{"$scaffold $start $direc"}},"$intron_id");
                        push(@{${$don_intron{$scaffold}}{"$scaffold $end $direc"}},"$intron_id");
		}
		$intron{"$intron_id"}="$scaffold $start $end $direc\n";
		$intron2{$intron_id}=$_;
	}
}
close(FILTERED);

#print Dumper %acc_intron

my %don_amb=();
my %acc_amb=();

my $num_introns=scalar(keys(%intron));
foreach my $each(sort(keys(%intron))){
	my @tmp=split(/\s/,$intron{$each});
	my $scaffold=$tmp[0];
	my $start=$tmp[1];
	my $end=$tmp[2];
	my $direc=$tmp[3];

	my $don="";
	my $acc="";
	if($direc eq "+"){
		$don=$start;
		$acc=$end;
	}elsif($direc eq "-"){
		$acc=$start;	
		$don=$end;
	}
	#print "$each	$scaffold $start $end $direc\n";
	
	for($n=$don-$mae;$n<=$don+$ato;$n++){
		if($n != $don){
			if((scalar(@{${$don_intron{$scaffold}}{"$scaffold $n $direc"}})>=1)){	
				my @array=@{${$don_intron{$scaffold}}{"$scaffold $n $direc"}};
				#print "[don]$n	@array\n";
				$don_amb{"$scaffold $n $direc"}="ambiguous";
				#die();
			}else{
				#print "[don]$n\n";
			}	
		}else{
			#print "([don]$n	$each)\n";
		}
	}
	#print "+\n";
	for($m=$acc-$mae;$m<=$acc+$ato;$m++){
		if($m != $acc){
                        if((scalar(@{${$acc_intron{$scaffold}}{"$scaffold $m $direc"}})>=1)){
        			my @array=@{${$acc_intron{$scaffold}}{"$scaffold $m $direc"}};	 
	                	#print "[acc]$m	@array\n";
                        	$acc_amb{"$scaffold $m $direc"}="ambiguous";
				#die();
			}else{
				#print "[acc]$m\n";
			}
                }else{
			#print "([acc]$m	$each)\n";
		}
	}
	#print "\n";
}	

print "Number of introns $num_introns in $filtered\n";
my $num_don_amb=scalar(keys(%don_amb));
my $num_acc_amb=scalar(keys(%acc_amb));
print "Number of ambiguous donor sites $num_don_amb\n";
print "Number of ambigupus accepter sites $num_acc_amb\n";


open(OUT,"+>$filtered.ambiguity_removed.2");
my $num_amb_intron=0;
my $num_filt_intron=0;
foreach my $each(sort(keys(%intron))){
        my @tmp=split(/\s/,$intron{$each});
        my $scaffold=$tmp[0];
        my $start=$tmp[1];
        my $end=$tmp[2];
        my $direc=$tmp[3];

        my $don="";
        my $acc="";
        if($direc eq "+"){
                $don=$start;
                $acc=$end;
        }elsif($direc eq "-"){
                $acc=$start;
                $don=$end;
        }

	if(($don_amb{"$scaffold $don $direc"}eq"ambiguous")or($acc_amb{"$scaffold $acc $direc"}eq"ambiguous")){
		print OUT "$intron2{$each}	... ambiguous\n";
		$num_amb_intron++;
	}else{
		print OUT "$intron2{$each}\n";
		$num_filt_intron++;
	}
}
print "Number of ambiguous introns $num_amb_intron\n";
print "Number of filtered introns $num_filt_intron\n";		
close(OUT);	


sub array_nr{
        my $array=$_[0];
        my %count=();
        my @nr_array=grep(!$count{$_}++,@$array);
        #print "@nr_array\n";
        return(@nr_array);
}
        
