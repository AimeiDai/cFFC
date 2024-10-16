#!/usr/bin/perl -w
use strict;

die "Usage: perl $0 TF_TG miR_TG expressed_TFs expressed_miRs miR_Family_Info.txt out" if scalar @ARGV != 6;

my $expressed_tfs;
open TFS,$ARGV[2] || die "open expressed_TFs failed!\n";
while (my $line = <TFS>){
	chomp $line;
	$expressed_tfs -> {$line} = 1;
}

my $expressed_mirs;
open MIRS,$ARGV[3] || die "open expressed_miRs failed!\n";
while (my $line = <MIRS>){
	chomp $line;
	$expressed_mirs -> {$line} = 1;
}

my $seeds;
open SEEDS, $ARGV[4] || die "open miR_Family_Info.txt failed!\n";
while (my $line = <SEEDS>){
	chomp $line;
	my @arr = split "\t",$line;
	next if ! defined $arr[4] || $arr[4] eq "-1";
	$seeds -> {$arr[1]} = $arr[0];
}

my $mir_tg;
open MIR,$ARGV[1] || die "open MIR failed!\n";
while (my $line = <MIR>) {
	# body...
	chomp $line;
	my @arr = split "\t",$line;
	next if $arr[4] ne "7227";
	if (! defined $seeds -> {$arr[0]}){
		next;
	}
	next if ! exists $expressed_mirs -> {$seeds->{$arr[0]}};
	$mir_tg -> {$arr[1]} -> {$seeds->{$arr[0]}} = 1;
		# target,miR,transcript
}


open TFTG,$ARGV[0] || die "open TFTG failed!\n";
open OUT,">$ARGV[5]" || die "open OUT failed!\n";
while (my $line = <TFTG>) {
	# body...
	chomp $line;
	my @arr = split "\t",$line;
	next if $arr[0] eq $arr[1];
	my @overlap_mirs;
	foreach my $mir (keys %{$mir_tg -> {$arr[1]}}){
		if( exists $mir_tg -> {$arr[2]} -> {$mir}){
			push @overlap_mirs,$mir;
		}
	}
	if (scalar @overlap_mirs != 0){
		my $num = scalar @overlap_mirs;
		foreach my $mir (@overlap_mirs){
			# my $x1 = scalar keys %{$mir_tg->{$arr[0]}->{$mir}};
			# my $x2 = scalar keys %{$mir_tg->{$arr[1]}->{$mir}};
			# $redundancy += $x1*$x2;
			print OUT $mir,"\t",$arr[1],"\t",$arr[2],"\t",$num,"\n";
		}
	}

}




