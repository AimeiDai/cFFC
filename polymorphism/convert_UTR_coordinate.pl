#!/usr/bin/perl -w
use strict;

die "Usage:perl $0 motif_with_redundancy.txt TSFly_r6.19_3UTRs.gff out" if scalar @ARGV !=3;

open UTR, $ARGV[1] || die "open UTR failed!\n";
my $utr;
while (my $line = <UTR>) {
	# body...
	chomp $line;
	my @arr = split "\t",$line;
	my ($fbtr) = ($arr[8] =~ /(FBtr\d+)?/sg);
	#$arr[0] =~ s/chr//sg;
	$utr -> {$fbtr} -> {"chr"} = $arr[0];
	$utr -> {$fbtr} -> {"strand"} = $arr[6];
	my @poses = $arr[3]..$arr[4];
	if (exists $utr -> {$fbtr} -> {"pos"}){
		@poses = (@poses, @{$utr -> {$fbtr} -> {"pos"}});
	}
	$utr -> {$fbtr} -> {"pos"} = \@poses;
}

open OUT, ">$ARGV[2]" || die "open OUT failed!\n";
open RED, $ARGV[0] || die "open RED failed!\n";
#my $no;
while (my $line = <RED>) {
	# body...
	chomp $line;
	my @arr = split "\t",$line;
	next if ! exists $utr -> {$arr[3]};
	my @pos;
	if ($utr -> {$arr[3]} -> {"strand"} eq "+"){
		@pos = sort{$a<=>$b} @{$utr -> {$arr[3]} -> {"pos"}};
		}else{
			@pos = sort{$b <=> $a} @{$utr -> {$arr[3]} -> {"pos"}};
		}
	foreach my $s (@pos[($arr[5]-1)..($arr[6]-1)]){
		if(!defined $s){
			#$no++;
			print $arr[3],"\n";
			next;
		}
		print OUT $utr -> {$arr[3]} -> {"chr"},"\t",$s,"\t",$arr[12],"\t",$arr[9],"\n";
	}
}

#print "no coordinate: ",$no,"\n";
