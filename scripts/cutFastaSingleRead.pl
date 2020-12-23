#!/usr/bin/perl
use strict;
use warnings;
	my $name;
	my $error_rate = 0;
	my $number_of_seq = 1;
	my $read_length = 300;
	my $insert_size = 0;
	my %hash = ();
	my $acc;
	my $hash_index=0;
open(FASTA,$ARGV[0]) || die("cannot open file!");
while(<FASTA>){
	my $line = $_;
	chomp($line);
	if ( $line =~ /^\>/ ){
		my @spl = split(/\>/,$line);
		$acc = $spl[1];
	}else{
		#@{$hash{$hash_index}}[0] = $acc;
		#@{$hash{$hash_index}}[1] = $line;
		#$hash_index++;
		my @spl = split(//,$line);
		my $sequence;
		for(my $i=0; $i<scalar(@spl); $i++){
			if ($spl[$i] ne "-"){
				$sequence .= $spl[$i];
			}
		}
		$hash{$acc} = $sequence;
	}
}
close(FASTA);

my $random_index=0;
my %random_sequences = ();

#for(my $i=0; $i<$number_of_seq; $i++){
#foreach my $key (sort keys %hash){
#	my $submitted = 0;
#	while ( $submitted==0 ){
#		$random_index = int(rand(keys %hash));
#		if ( not exists $random_sequences{@{$hash{$random_index}}[0]} ){
#		$random_sequences{@{$hash{$random_index}}[0]} = @{$hash{$random_index}}[1];
#		$random_index++;
#			$submitted = 1;
#		}
#	}
#}
my $numreads = 0;
open(READ1,'>', $ARGV[1] . '.simulated_300_0.fasta') || die("Cannot open file!");
foreach my $key(keys %hash){
	my @spl = split(//,$hash{$key});
	my $index=0;
	while($index<=scalar(@spl)-$read_length){
		#if ( $numreads < 28348 ){
		print READ1 ">$key\_$index\n";
		for (my $i=$index; $i<$index+$read_length; $i++){
			my $base=$spl[$i];
			my $random_number = rand();
			if ($random_number < $error_rate){
				my $random_nuc = rand();
				if ( $base eq "A" and $random_nuc <= 0.4918 ){
				print READ1 "C";
				}elsif ( $base eq "A" and $random_nuc > 0.4918 and $random_nuc < 0.8295 ){
					print READ1 "G";
				}elsif ( $base eq "A" and $random_nuc >= 0.8295 ){
					print READ1 "T";
				}elsif ( $base eq "C" and $random_nuc <=0.5238 ){
					print READ1 "A";
				}elsif ( $base eq "C" and $random_nuc > 0.5238 and $random_nuc < 0.7899 ){
					print READ1 "G";
				}elsif ( $base eq "C" and $random_nuc >= 0.7899 ){
					print READ1 "T";
				}elsif ( $base eq "G" and $random_nuc <= 0.3754 ){
					print READ1 "A";
				}elsif ( $base eq "G" and $random_nuc > 0.3754 and $random_nuc < 0.6109 ){
					print READ1 "C";
				}elsif ( $base eq "G" and $random_nuc >= 0.6109 ){
					print READ1 "T";
				}elsif ( $base eq "T" and $random_nuc <= 0.2505 ){
					print READ1 "A";
				}elsif ( $base eq "T" and $random_nuc > 0.2505 and $random_nuc < 0.5057 ){
					print READ1 "C";
				}elsif ( $base eq "T" and $random_nuc >= 0.5057 ){
					print READ1 "G";
				}
			}else{
				print READ1 "$base";
			}
		}
		print READ1 "\n";
		$index++;
		$numreads++;
		#}else{
		#	last;
		#}
	}
}
close(READ1);

