#!/usr/bin/perl -w
use strict;
use warnings;
use Time::HiRes;
if (scalar(@ARGV) < 2){
	print "usage:\n";
	print "elimination_sam2mismatches.pl MSA.fasta alignment.sam output.txt\n";
	print "outputs: output.txt\n";
	exit;
}
my $msa_file = $ARGV[0];
my $sam_file = $ARGV[1];
my @splitname = split(/\./,$sam_file);
my $outfile = $ARGV[2];
my %chrom_sizes = ();
$chrom_sizes{"NC_045512v2"}=29903;
my $ref_seq=0;
my %reference_map;
my %reference_genomes;
my %allele_frequencies = ();
my $name;
my $msa_length=0;
my $total_start = Time::HiRes::gettimeofday();
my $start_time = Time::HiRes::gettimeofday();
print "Reading in MSA... \n";
open(MSA,$msa_file) || die("cannot open $msa_file!");
while(<MSA>){
	my $line = $_;
	chomp($line);
	if ($line =~ /^\>/){
		my @spl = split(/\>/,$line);
		$name = $spl[1];
		if ( $name eq "NC_045512v2" ){
			$ref_seq = 1;
		}
	}elsif($ref_seq==1){
		my @spl = split(//,$line);
		$msa_length=scalar(@spl);
		for(my $i=0; $i<scalar(@spl); $i++){
			push(@{$reference_genomes{$name}},$spl[$i]);
		}
		my $ref_pos=0;
		for(my $i=0; $i<scalar(@spl); $i++){
			if ($spl[$i] ne "-"){
				$reference_map{$ref_pos}=$i;
				$ref_pos++;		
			}
		}
		$ref_seq=0;
	}else{
		my @spl = split(//,$line);
		for(my $i=0; $i<scalar(@spl); $i++){
			push(@{$reference_genomes{$name}},$spl[$i]);
		}
	}
}
close(MSA);
my $stop_time = Time::HiRes::gettimeofday();
printf("Took %.2f\n",$stop_time - $start_time);
$start_time = Time::HiRes::gettimeofday();
print "Calculating allele frequencies... $start_time\n";
open (SAM, $sam_file) || die ("Could not open sam file!");
while (<SAM>){
        my $line = $_;
        chomp($line);
	if ($line !~ m/^\@/){
		my @spl = split('\t',$line);
		my $chrom = $spl[2];
		my $tSize = $chrom_sizes{$chrom};
		#if there is a match, process line
		if ($chrom ne "*"){
			my $qName = $spl[0];
			my $flag = $spl[1];
			my $tStart = $spl[3] - 1;
			my $cigar = $spl[5];
			my $sequence = $spl[9];
			my @sequence_chars = split(//,$sequence);
			my $misMatches = $spl[16];
			my @spl5 = split("NM:i:",$misMatches);
			$misMatches = $spl5[1];
			my $qSize = seqLength($sequence);
			my $totalmatches = 0;
			my $blockCount = 0;
			#process cigar string
			if ( $cigar =~ m/M/ ){
				my @counts = getCount($cigar,"M");
				$totalmatches = $counts[0];
				$blockCount = $counts[1];
			}
			my $nCount = 0;
			if ($cigar =~ m/N/ ){
				my @counts = getCount($cigar,"N");
				$nCount = $counts[0];
			}
			my $tBaseInsert = 0;
			my $tNumInsert = 0;
			if ($cigar =~ m/I/ ){
				my @counts = getCount($cigar,"I");
				$tBaseInsert = $counts[0];
				$tNumInsert = $counts[1];
			}
			my $qBaseInsert = 0;
			my $qNumInsert = 0;
			if ($cigar =~ m/D/ ){
				my @counts = getCount($cigar,"D");
				$qBaseInsert = $counts[0];
				$qNumInsert = $counts[1];
			}
			my $qStart = 0;
			if ($cigar =~ m/^[0-9]+S/ ){
				my @counts = getCount($cigar,"S");
				$qStart = $counts[2];
			}elsif( $cigar =~ m/^[0-9]+H/ ){
				my @counts = getCount($cigar,"H");
				$qStart = $counts[2];
			}
			my $qEnd = 0;
			my $matches = $totalmatches - $misMatches;
			my @qStart_blocks;
			$qStart_blocks[$blockCount] = undef;
			$qStart_blocks[0] = $qStart;
			my @blockSizes;
			my @qStarts;
			my $qStart_sum=$qStart;
			my @tStarts;
			my $tStart_sum = $tStart;
			if ($cigar =~ m/M/){
				my $tmp = $cigar;
				$tmp =~ s/[DNM=XISHP]/|/g;
				my @spl = split(/\Q|/,$tmp);	
				my $tmp2 = $cigar;
				$tmp2 =~ s/[0-9]+/|/g;
				my @spl2 = split(/\Q|/,$tmp2);
				for(my $i=0; $i<scalar(@spl2); $i++){
					if ($spl2[$i] eq 'M' or $spl2[$i] eq 'I'){
						if ($spl2[$i] eq 'M'){
							push(@blockSizes,$spl[$i-1]);
							my $scale = scalar(@blockSizes);
							if ( scalar(@blockSizes) > 1 ){
								push(@qStarts,$qStart_sum);
							}
							else{ push(@qStarts, $qStart); }
							$qStart_sum = $qStart_sum + $spl[$i-1];
						}else{
							$qStart_sum = $qStart_sum + $spl[$i-1];
						}	
					}
					if ($spl2[$i] eq 'M' or $spl2[$i] eq 'D' or $spl2[$i] eq 'N'){
						if ($spl2[$i] eq 'M'){
							if ( scalar(@tStarts) > 0 ){
								push(@tStarts,$reference_map{$tStart_sum});
							}else{ 
								push(@tStarts,$reference_map{$tStart}); 
							}
							$tStart_sum = $tStart_sum + $spl[$i-1];
						}else{
							$tStart_sum = $tStart_sum + $spl[$i-1];
						}

	

					}	
				}
			
			}
			$qEnd = $qStarts[scalar(@qStarts)-1]+$blockSizes[scalar(@qStarts)-1];
			if ($cigar =~ m/[0-9]+H$/ ){
                                my @counts = getCount($cigar,"H");
                                $qEnd = $qSize - $counts[3];
                        }elsif ($cigar =~ m/[0-9]+S$/ ){
                                my @counts = getCount($cigar, "S");
                                $qEnd = $qSize - $counts[3];
                        }
			$cigar =~ s/([0-9]+[ISHP])+//g;
			my @spl2 = split('[MDN=X]',$cigar);
			my $tAlignment=0;
			for(my $i=0; $i<scalar(@spl2); $i++){
				$tAlignment = $tAlignment + $spl2[$i];
			}
			#my $tEnd = $tStart + $tAlignment - 2;
			my $tEnd = $tStart + $tAlignment;
			my @positions;
			for(my $i=$tStart; $i<$tEnd; $i++){
				push(@positions,$reference_map{$i});
			}
			my @spl3 = split('',$flag);
			my $binary = dec2bin($flag);
			my @spl4 = split('',$binary);
			my $strand_code = $spl4[7];
			my $strand;
			if ($strand_code == 1){
				$strand = "-";
			}elsif($strand_code == 0){
				$strand = "+";
			}
			my $repMatches = 0;
			#print OUTFILE "$qName\t$blockSizes[0]";
			my $blockSizes_print = join(",",@blockSizes);
			my $qStarts_print = join(",",@qStarts);
			my $tStarts_print = join(",",@tStarts);
			my $j;
			if (scalar(@blockSizes) > 1){
				print "$qName has more than 1 block\n";
				exit(1);
			}
			for(my $i=0; $i<scalar(@positions); $i++){
				$j=$i+$qStart;
				if (uc($sequence_chars[$j]) eq "A" ){
					$allele_frequencies{$positions[$i]}{"A"}++;
				}elsif( uc($sequence_chars[$j]) eq "G" ){
					$allele_frequencies{$positions[$i]}{"G"}++;
				}elsif( uc($sequence_chars[$j]) eq "C" ){
					$allele_frequencies{$positions[$i]}{"C"}++;
				}elsif( uc($sequence_chars[$j]) eq "T" ){
					$allele_frequencies{$positions[$i]}{"T"}++;
				}elsif ( uc($sequence_chars[$j]) ne "N"){
					print "strange: $sequence_chars[$j]\n";
				}
			}
		}
	}
}
close(SAM);
foreach my $position ( sort {$a cmp $b} (keys %allele_frequencies) ){
	my $total = 0;
	my $A_freq=0;
	my $G_freq=0;
	my $C_freq=0;
	my $T_freq=0;
	foreach my $base (keys %{$allele_frequencies{$position}}){
		$total += $allele_frequencies{$position}{$base};
	}
	foreach my $base (keys %{$allele_frequencies{$position}}){
		if ( $base eq "A" ){
			$A_freq = $allele_frequencies{$position}{"A"}/$total;
		}elsif( $base eq "G" ){
			$G_freq = $allele_frequencies{$position}{"G"}/$total;
		}elsif( $base eq "C" ){
			$C_freq = $allele_frequencies{$position}{"C"}/$total;
		}elsif( $base eq "T" ){
			$T_freq = $allele_frequencies{$position}{"T"}/$total;
		}
	}
	$allele_frequencies{$position}{"A"} = $A_freq;
	$allele_frequencies{$position}{"G"} = $G_freq;
	$allele_frequencies{$position}{"C"} = $C_freq;
	$allele_frequencies{$position}{"T"} = $T_freq;
}
$stop_time = Time::HiRes::gettimeofday();
printf("Took %.2f\n",$stop_time - $start_time);
$start_time = Time::HiRes::gettimeofday();
print "Eliminating sequences...\n";
my $num_eliminated = 0;
foreach my $i ( sort {$a <=> $b} (keys %allele_frequencies) ){
	foreach my $key (sort keys %reference_genomes){
		my $nuc = uc(@{$reference_genomes{$key}}[$i]);
		if ( $nuc ne "A" and $nuc ne "G" and $nuc ne "C" and $nuc ne "T"){ next; }
		if ( $allele_frequencies{$i}{$nuc} < 0.03 ){
			delete($reference_genomes{$key});
			$num_eliminated++;
		}
	}
}
my $ref_left = (keys %reference_genomes);
print "Number eliminated: $num_eliminated\n";
print "Number of reference genomes left: $ref_left\n";
$stop_time = Time::HiRes::gettimeofday();
printf("Took %.2f\n",$stop_time - $start_time);
if ( $ref_left == 1){
	my $key_1 = (sort { $a cmp $b } keys %reference_genomes )[0];
	print "only 1 reference left: $key_1\n";
	exit(1);
}
$start_time = Time::HiRes::gettimeofday();
print "Finding indifferent positions... \n";
my %indifferent_positions = ();
my $first_key = (sort { $a cmp $b } keys %reference_genomes )[0];
print "First key is $first_key\n";
for(my $i=0; $i<$msa_length; $i++){
	my $different = 0;
	foreach my $key (sort { $a cmp $b } keys %reference_genomes){
		if ( @{$reference_genomes{$first_key}}[$i] ne "-" and @{$reference_genomes{$key}}[$i] ne "-" and @{$reference_genomes{$first_key}}[$i] ne @{$reference_genomes{$key}}[$i] ){
			$different=1;
			last;
		}
	}
	if ($different == 0){
		$indifferent_positions{$i}=1;
	}
}
my $num_indifferent = (keys %indifferent_positions);
print "Number of indifferent sites: $num_indifferent\n";
$stop_time = Time::HiRes::gettimeofday();
printf("Took %.2f\n",$stop_time - $start_time);
$start_time = Time::HiRes::gettimeofday();
print "Creating mismatch matrix...\n";
open (OUTFILE, ">$outfile") || die ("Could not open $outfile!");
print OUTFILE "qName\tblockSizes";
foreach my $key (sort { $a cmp $b } keys %reference_genomes){
	print OUTFILE "\t$key";
}
print OUTFILE "\n";
open (SAM, $sam_file) || die ("Could not open sam file!");
while (<SAM>){
        my $line = $_;
        chomp($line);
	if ($line !~ m/^\@/){
		my @spl = split('\t',$line);
		my $chrom = $spl[2];
		my $tSize = $chrom_sizes{$chrom};
		#if there is a match, process line
		if ($chrom ne "*"){
			my $qName = $spl[0];
			my $flag = $spl[1];
			my $tStart = $spl[3] - 1;
			my $cigar = $spl[5];
			my $sequence = $spl[9];
			my @sequence_chars = split(//,$sequence);
			my $misMatches = $spl[16];
			my @spl5 = split("NM:i:",$misMatches);
			$misMatches = $spl5[1];
			my $qSize = seqLength($sequence);
			my $totalmatches = 0;
			my $blockCount = 0;
			#process cigar string
			if ( $cigar =~ m/M/ ){
				my @counts = getCount($cigar,"M");
				$totalmatches = $counts[0];
				$blockCount = $counts[1];
			}
			my $nCount = 0;
			if ($cigar =~ m/N/ ){
				my @counts = getCount($cigar,"N");
				$nCount = $counts[0];
			}
			my $tBaseInsert = 0;
			my $tNumInsert = 0;
			if ($cigar =~ m/I/ ){
				my @counts = getCount($cigar,"I");
				$tBaseInsert = $counts[0];
				$tNumInsert = $counts[1];
			}
			my $qBaseInsert = 0;
			my $qNumInsert = 0;
			if ($cigar =~ m/D/ ){
				my @counts = getCount($cigar,"D");
				$qBaseInsert = $counts[0];
				$qNumInsert = $counts[1];
			}
			my $qStart = 0;
			if ($cigar =~ m/^[0-9]+S/ ){
				my @counts = getCount($cigar,"S");
				$qStart = $counts[2];
			}elsif( $cigar =~ m/^[0-9]+H/ ){
				my @counts = getCount($cigar,"H");
				$qStart = $counts[2];
			}
			my $qEnd = 0;
			my $matches = $totalmatches - $misMatches;
			my @qStart_blocks;
			$qStart_blocks[$blockCount] = undef;
			$qStart_blocks[0] = $qStart;
			my @blockSizes;
			my @qStarts;
			my $qStart_sum=$qStart;
			my @tStarts;
			my $tStart_sum = $tStart;
			if ($cigar =~ m/M/){
				my $tmp = $cigar;
				$tmp =~ s/[DNM=XISHP]/|/g;
				my @spl = split(/\Q|/,$tmp);	
				my $tmp2 = $cigar;
				$tmp2 =~ s/[0-9]+/|/g;
				my @spl2 = split(/\Q|/,$tmp2);
				for(my $i=0; $i<scalar(@spl2); $i++){
					if ($spl2[$i] eq 'M' or $spl2[$i] eq 'I'){
						if ($spl2[$i] eq 'M'){
							push(@blockSizes,$spl[$i-1]);
							my $scale = scalar(@blockSizes);
							if ( scalar(@blockSizes) > 1 ){
								push(@qStarts,$qStart_sum);
							}
							else{ push(@qStarts, $qStart); }
							$qStart_sum = $qStart_sum + $spl[$i-1];
						}else{
							$qStart_sum = $qStart_sum + $spl[$i-1];
						}	
					}
					if ($spl2[$i] eq 'M' or $spl2[$i] eq 'D' or $spl2[$i] eq 'N'){
						if ($spl2[$i] eq 'M'){
							if ( scalar(@tStarts) > 0 ){
								push(@tStarts,$reference_map{$tStart_sum});
							}else{ 
								push(@tStarts,$reference_map{$tStart}); 
							}
							$tStart_sum = $tStart_sum + $spl[$i-1];
						}else{
							$tStart_sum = $tStart_sum + $spl[$i-1];
						}

	

					}	
				}
			
			}
			$qEnd = $qStarts[scalar(@qStarts)-1]+$blockSizes[scalar(@qStarts)-1];
			if ($cigar =~ m/[0-9]+H$/ ){
                                my @counts = getCount($cigar,"H");
                                $qEnd = $qSize - $counts[3];
                        }elsif ($cigar =~ m/[0-9]+S$/ ){
                                my @counts = getCount($cigar, "S");
                                $qEnd = $qSize - $counts[3];
                        }
			$cigar =~ s/([0-9]+[ISHP])+//g;
			my @spl2 = split('[MDN=X]',$cigar);
			my $tAlignment=0;
			for(my $i=0; $i<scalar(@spl2); $i++){
				$tAlignment = $tAlignment + $spl2[$i];
			}
			#my $tEnd = $tStart + $tAlignment - 2;
			my $tEnd = $tStart + $tAlignment;
			my @positions;
			for(my $i=$tStart; $i<$tEnd; $i++){
				push(@positions,$reference_map{$i});
			}
			my @spl3 = split('',$flag);
			my $binary = dec2bin($flag);
			my @spl4 = split('',$binary);
			my $strand_code = $spl4[7];
			my $strand;
			if ($strand_code == 1){
				$strand = "-";
			}elsif($strand_code == 0){
				$strand = "+";
			}
			my $repMatches = 0;
			print OUTFILE "$qName\t$blockSizes[0]";
			my $blockSizes_print = join(",",@blockSizes);
			my $qStarts_print = join(",",@qStarts);
			my $tStarts_print = join(",",@tStarts);
			my $j;
			$first_key = (sort { $a cmp $b } keys %reference_genomes )[0];
			foreach my $key ( sort { $a cmp $b } keys %reference_genomes ){
				my $number_of_mismatches=0;
				for(my $i=0; $i<scalar(@positions); $i++){
					$j=$i+$qStart;
					#if ( uc(@{$reference_genomes{$key}}[$positions[$i]]) ne "A" and uc(@{$reference_genomes{$key}}[$positions[$i]]) ne "G" and uc(@{$reference_genomes{$key}}[$positions[$i]]) ne "C" and uc(@{$reference_genomes{$key}}[$positions[$i]]) ne "T" ){ next;}
					#if ( uc($sequence_chars[$j]) ne "A" and uc($sequence_chars[$j]) ne "G" and uc($sequence_chars[$j]) ne "C" and uc($sequence_chars[$j]) ne "T" ){ next; }
					#if ( not exists $indifferent_positions{$positions[$i]} and uc(@{$reference_genomes{$key}}[$positions[$i]]) ne uc($sequence_chars[$j]) and @{$reference_genomes{$key}}[$positions[$i]] ne "-" and uc(@{$reference_genomes{$key}}[$positions[$i]]) ne "N"){
					if ( @{$reference_genomes{$key}}[$positions[$i]] =~ /[AaGgCcTt]/ and $sequence_chars[$j] =~ /[AaGgCcTt]/ and uc(@{$reference_genomes{$key}}[$positions[$i]]) ne uc($sequence_chars[$j]) ){
						$number_of_mismatches++;
					}elsif ( @{$reference_genomes{$first_key}}[$positions[$i]] ne "-" and exists $indifferent_positions{$positions[$i]} and uc(@{$reference_genomes{$first_key}}[$positions[$i]]) ne uc($sequence_chars[$j]) ){
						$number_of_mismatches++;
						#print "$sequence_chars[$j]\t@{$reference_genomes{$first_key}}[$positions[$i]]\n";
					}
				}
				print OUTFILE "\t$number_of_mismatches";
			}
			print OUTFILE "\n";
		}
	}
}
close(SAM);
$stop_time = Time::HiRes::gettimeofday();
printf("Took %.2f\n",$stop_time - $start_time);
my $total_stop = Time::HiRes::gettimeofday();
printf("Entire script took %.2f\n",$total_stop-$total_start);
#converts to binary
sub dec2bin {
	my $str = unpack("B32", pack("N", shift));
	my @spl = split('',$str);
	my $binary = $spl[20];
	for (my $i=21; $i<=31; $i++){
		$binary .= $spl[$i];
	}	
	return $binary;
}
#get length of sequence
sub seqLength {
	my @spl = split('',$_[0]);
	return scalar(@spl);
}
sub getCount {
	my $cigar = $_[0];
	my $symbol = $_[1];
	my @spl;
	my $beg=0;
	my $end=0;
	if ($symbol eq 'M'){
		$cigar =~ s/([0-9]+[ISHPDN=X])+//g;
		@spl = split('M', $cigar);
	}
	if ($symbol eq 'N'){
		$cigar =~ s/([0-9]+[ISHPDM=X])+//g;
		@spl = split('N',$cigar);
	}
	if ($symbol eq 'I'){
		$cigar =~ s/([0-9]+[NSHPDM=X])+//g;
		@spl = split('I',$cigar);
	}
	if ($symbol eq 'D'){
		$cigar =~ s/([0-9]+[NSHPIM=X])+//g;
		@spl = split('D',$cigar);
	}
	if ($symbol eq 'S'){
		if ( $cigar =~ m/^[0-9]+S/){
		     	$cigar =~ s/([0-9]+[NDHPIM=X])+//g;
			@spl = split('S',$cigar);
			if ( scalar(@spl) > 1){
                        	$beg = $spl[0];
				$end = $spl[1];
			}else{
				$beg = $spl[0];
			}
                }else{
			$cigar =~ s/([0-9]+[NDHPIM=X])+//g;
			@spl = split('S',$cigar);
			$end = $spl[0];
		}
	}
	if ($symbol eq 'H'){
		if ( $cigar =~ m/^[0-9]+H/){
			$cigar =~ s/([0-9]+[NDSPIM=X])+//g;
			@spl = split('H',$cigar);
			if (scalar(@spl) > 1){
				$beg = $spl[0];
				$end = $spl[1];
			}else{
				$beg=$spl[0];
			}
		}else{
			$cigar =~ s/([0-9]+[NDSPIM=X])+//g;
			@spl = split('H',$cigar);
			$end = $spl[0];
		}
	}
	my $count=0;
	for(my $i=0; $i<scalar(@spl); $i++){
		$count = $count + $spl[$i];
	}
	return $count,scalar(@spl),$beg,$end;
}
