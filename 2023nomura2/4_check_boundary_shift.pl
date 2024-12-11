use Data::Dumper;

#my $file="Eag.draft-genome.20160925-3.gene.id_edited.cds_edited.gff3.intron.filtered.ambiguity_removed.2";
#my $file="Eag.draft-genome.20180531.gene.id_edited.cds_edited.gff3.intron.filtered.ambiguity_removed.2";
#my $file="Eag_consensusScaffold.1kb.rename.gff3.intron.filtered.ambiguity_removed.2";
#my $file="Eag.ref.gff3.intron.filtered.ambiguity_removed.2";
#my $file="gmap.Eag_new-r2_comp.gff.intron.filtered.ambiguity_removed.2";
my $file="Euglena.gene.change.chr.gff.intron.filtered.ambiguity_removed.2";
open(OUT,"+>$file.shift.check");

print OUT "#Intron	Type1	DD-AA	EEEEEEEEEEEDDIIIIIIII .. IIIIIIIIAAEEEEEEEEEEE	Type2	No.Shift	Shift DD-AA	Shift bps	Conventional	U12	Intermediate	Nonconventional\n";

open(FILE,$file);
my $all=0;
my $shift=0;
my $conv=0;
my $u12=0;
my $nonconv=0;
my $intermed=0;


while(<FILE>){
	s/\r//;
	chomp;
	my @tmp=split(/\t/,$_);
	my $seq=$tmp[3];
	print "#@tmp\n";
	#print "#$seq\n";

	my %hash=();
	
	my @donacc=();
	my @shift=();
	if($seq=~/^(...........)(.)(.........) \.\. (..........)(.)(..........)$/){
                print "+1 $1 $2 $3 .. $4 $5 $6";
		my @tmp=split("","$1$2$3$4$5$6");
		$don="$tmp[12]$tmp[13]";
		$acc="$tmp[30]$tmp[31]";
		print "	$don $acc";
		${$hash{"+1"}}{"seq"}="$1 $2 $3 .. $4 $5 $6";
		${$hash{"+1"}}{"don-acc"}="$don-$acc";
		if($2 eq $5){
			push(@donacc,"$don-$acc");
			push(@shift,"+1");
			${$hash{"+1"}}{"shift"}="yes";
			print "	*\n";
		}else{
			${$hash{"+1"}}{"shift"}="no";
			print "	\n";
		}	
	}
	if($seq=~/^(...........)(..)(........) \.\. (..........)(..)(.........)$/){
                print "+2 $1 $2 $3 .. $4 $5 $6";
        	my @tmp=split("","$1$2$3$4$5$6");
                $don="$tmp[13]$tmp[14]";
                $acc="$tmp[31]$tmp[32]";
                print "	$don $acc";
		${$hash{"+2"}}{"seq"}="$1 $2 $3 .. $4 $5 $6";
                ${$hash{"+2"}}{"don-acc"}="$don-$acc";
		if($2 eq $5){
                        push(@donacc,"$don-$acc");
                	push(@shift,"+2");
			${$hash{"+2"}}{"shift"}="yes";
			print "	*\n";
		}else{
			${$hash{"+2"}}{"shift"}="no";
			print "	\n";
		}

	}
	if($seq=~/^(...........)(...)(.......) \.\. (..........)(...)(........)$/){
                print "+3 $1 $2 $3 .. $4 $5 $6";
        	my @tmp=split("","$1$2$3$4$5$6");
                $don="$tmp[14]$tmp[15]";
                $acc="$tmp[32]$tmp[33]";
                print "	$don $acc";
		${$hash{"+3"}}{"seq"}="$1 $2 $3 .. $4 $5 $6";
                ${$hash{"+3"}}{"don-acc"}="$don-$acc";
		if($2 eq $5){
                        push(@donacc,"$don-$acc");
			push(@shift,"+3");
			${$hash{"+3"}}{"shift"}="yes";
                	print "	*\n";
		}else{
			${$hash{"+3"}}{"shift"}="no";
			print "	\n";
		}
	}
	if($seq=~/^(...........)(....)(......) \.\. (..........)(....)(.......)$/){
                print "+4 $1 $2 $3 .. $4 $5 $6";
        	my @tmp=split("","$1$2$3$4$5$6");
                $don="$tmp[15]$tmp[16]";
                $acc="$tmp[33]$tmp[34]";
                print "	$don $acc";
		${$hash{"+4"}}{"seq"}="$1 $2 $3 .. $4 $5 $6";
                ${$hash{"+4"}}{"don-acc"}="$don-$acc";

		if($2 eq $5){
                        push(@donacc,"$don-$acc");
                	${$hash{"+4"}}{"shift"}="yes";
			push(@shift,"+4");
			print "	*\n";
		}else{
			${$hash{"+4"}}{"shift"}="no";
			print "	\n";
		}
	}
	if($seq=~/^(...........)(.....)(.....) \.\. (..........)(.....)(......)$/){
                print "+5 $1 $2 $3 .. $4 $5 $6";
        	my @tmp=split("","$1$2$3$4$5$6");
                $don="$tmp[16]$tmp[17]";
                $acc="$tmp[34]$tmp[35]";
                print "	$don $acc";
		${$hash{"+5"}}{"seq"}="$1 $2 $3 .. $4 $5 $6";
                ${$hash{"+5"}}{"don-acc"}="$don-$acc";

		if($2 eq $5){
                        push(@donacc,"$don-$acc");
                	${$hash{"+5"}}{"shift"}="yes";
			push(@shift,"+5");
			print "	*\n";
		}else{
			${$hash{"+5"}}{"shift"}="no";
			print "	\n";
		}
	}


	if($seq=~/^(..........)(.)(..........) \.\. (.........)(.)(...........)$/){
                print "-1 $1 $2 $3 .. $4 $5 $6";
        	my @tmp=split("","$1$2$3$4$5$6");
                $don="$tmp[10]$tmp[11]";
                $acc="$tmp[28]$tmp[29]";
                print "	$don $acc";
		${$hash{"-1"}}{"seq"}="$1 $2 $3 .. $4 $5 $6";
                ${$hash{"-1"}}{"don-acc"}="$don-$acc";

		if($2 eq $5){
                        push(@donacc,"$don-$acc");
                	${$hash{"-1"}}{"shift"}="yes";
			push(@shift,"-1");
			print "	*\n";
		}else{
			${$hash{"-1"}}{"shift"}="no";
			print "	\n";
		}
	}
        if($seq=~/^(.........)(..)(..........) \.\. (........)(..)(...........)$/){
                print "-2 $1 $2 $3 .. $4 $5 $6";
        	my @tmp=split("","$1$2$3$4$5$6");
                $don="$tmp[9]$tmp[10]";
                $acc="$tmp[27]$tmp[28]";
                print "	$don $acc";
		${$hash{"-2"}}{"seq"}="$1 $2 $3 .. $4 $5 $6";
                ${$hash{"-2"}}{"don-acc"}="$don-$acc";

		if($2 eq $5){
                        push(@donacc,"$don-$acc");
                	${$hash{"-2"}}{"shift"}="yes";
			push(@shift,"-2");
			print "	*\n";
		}else{
			${$hash{"-2"}}{"shift"}="no";
			print "	\n";
		}


	}
        if($seq=~/^(........)(...)(..........) \.\. (.......)(...)(...........)$/){
                print "-3 $1 $2 $3 .. $4 $5 $6";
        	my @tmp=split("","$1$2$3$4$5$6");
                $don="$tmp[8]$tmp[9]";
                $acc="$tmp[26]$tmp[27]";
                print "	$don $acc";
		${$hash{"-3"}}{"seq"}="$1 $2 $3 .. $4 $5 $6";
                ${$hash{"-3"}}{"don-acc"}="$don-$acc";

		if($2 eq $5){
                        push(@donacc,"$don-$acc");
                	${$hash{"-3"}}{"shift"}="yes";
			push(@shift,"-3");
			print "	*\n";
		}else{
			${$hash{"-3"}}{"shift"}="no";
			print "	\n";
		}

	}
        if($seq=~/^(.......)(....)(..........) \.\. (......)(....)(...........)$/){
                print "-4 $1 $2 $3 .. $4 $5 $6";
        	my @tmp=split("","$1$2$3$4$5$6");
                $don="$tmp[7]$tmp[8]";
                $acc="$tmp[25]$tmp[26]";
                print "	$don $acc";
		${$hash{"-4"}}{"seq"}="$1 $2 $3 .. $4 $5 $6";
                ${$hash{"-4"}}{"don-acc"}="$don-$acc";

		if($2 eq $5){
                        push(@donacc,"$don-$acc");
                	${$hash{"-4"}}{"shift"}="yes";
			push(@shift,"-4");
			print "	*\n";
		}else{
			${$hash{"-4"}}{"shift"}="no";
			print "	\n";
		}
	}
        if($seq=~/^(......)(.....)(..........) \.\. (.....)(.....)(...........)$/){
                print "-5 $1 $2 $3 .. $4 $5 $6";
        	my @tmp=split("","$1$2$3$4$5$6");
                $don="$tmp[6]$tmp[7]";
                $acc="$tmp[24]$tmp[25]";
                print "	$don $acc";
		${$hash{"-5"}}{"seq"}="$1 $2 $3 .. $4 $5 $6";
                ${$hash{"-5"}}{"don-acc"}="$don-$acc";

		if($2 eq $5){
                        push(@donacc,"$don-$acc");
                	${$hash{"-5"}}{"shift"}="yes";
			push(@shift,"-5");
			print "	*\n";
		}else{
			${$hash{"-5"}}{"shift"}="no";
			print "	\n";
		}

	}

	#print Dumper %hash;
	my $type="";
	if(($_=~/ambiguous/)or($_=~/Amb/)){
		$type="ambiguous";
	}else{
		if($tmp[1]=~/Non-conventional/){
                        $type="$tmp[1] $tmp[2]";
                }elsif($tmp[1]=~/Conventional/){
                        $type="$tmp[1]";
                }else{
                        #print "$tmp[1]\n";
                        die();
                }
	}
	my $motodata=join("\t",@tmp[0..3]);
	my $num_shift=0;
	$num_shift=scalar(@donacc);
	my @conv=();
	my @U12=();
	my @nonconv=();	
	my @intermed=();
	
	foreach my $each(@donacc){
		if(($each =~/GT-AG/)or($each=~/GC-AG/)){
			push(@conv,$each);
		}elsif($each =~/AT-A./){
			push(@U12,$each);
		}elsif(($each=~/GT-/)or($each=~/-AG/)){
			push(@intermed,$each);
		}else{
			push(@nonconv,$each);
		}
	}
	print OUT "$motodata	$type	$num_shift	@donacc	@shift	@conv	@U12	@intermed	@nonconv\n";
	print "$motodata	$type	$num_shift	@donacc	@shift	@conv	@U12	@intermed	@nonconv\n";

	print "//\n";
	if(scalar(@donacc)>=1){
		$shift++;
		if(scalar(@conv)>=1){
			$conv++;
		}
		if(scalar(@U12)>=1){
			$u12++;
		}
		if(scalar(@intermed)>=1){
			$intermed++;
		}
		if(scalar(@nonconv)>=1){
			$nonconv++;
		}
		if(($type=~/Non-conventional/)and(scalar(@conv)>=1)){
			$revertant++;
		}

	}	
	$all++;

	
}
close(FILE);

print "all $all	shift $shift	conv $conv u12 $u12 int $intermed nonconv $nonconv	revertant $revertant\n";



