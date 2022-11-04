#!/usr/bin/perl
use strict;
use warnings;

my %seqs = ();
open(META,$ARGV[0]) || die("cannot open file!");
while(<META>){
	my $line = $_;
	chomp($line);
	if ($line !~ /^accession/){
	my @spl = split(/\,/,$line);
	#$seqs{$spl[1]}=$spl[0];
	$seqs{$spl[0]}=$spl[0];
	}
}
close(META);
my $name;
my $sequence;
my $first=0;
my $tree_size = keys %seqs;
my %msa_ids = ();
open(OUT_MSA,'>','seqs.fasta') || die("cannot open file!");
open(MSA,$ARGV[1]) || die("cannot open file!");
while(<MSA>){
	my $line = $_;
	chomp($line);
	if ($line =~/^\>/){
		if ( $first==1){
			if ( exists $seqs{$name} ){
				print OUT_MSA ">$seqs{$name}\n$sequence\n";
				$msa_ids{$name}=0;
			}
		}
		if ( $first==0){
			$first=1;
		}
		my @spl = split(/\>/,$line);
		#@spl = split(/\|/,$spl[1]);
		#$name=$spl[0];
		$name = $spl[1];
		$sequence = undef;
	}else{
		if ( not defined $sequence ){
			$sequence =  $line;
		}else{
			$sequence = $sequence . $line;
		}
	}
}
close(MSA);
close(OUT_MSA);
my $msa_size = keys %msa_ids;
if ($msa_size != $tree_size){
	open(MISSING,'>','missing.txt') || die("cannot open file!");
	foreach my $key (keys %seqs){
		if ( not exists $msa_ids{$key} ){
			print MISSING "$key\n";
		}
	}
	close(MISSING);
}else{
	print "success!\n";
}
