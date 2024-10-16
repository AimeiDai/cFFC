#!/usr/bin/perl -w
use strict;

die "Usage: perl $0 divergence.file site_type\n" if scalar @ARGV != 2;
my ($div_file, $site_type) = @ARGV;

my $div;
open DIV ,$div_file || die "open $div_file failed!\n";
while (my $line = <DIV>) {
	# body...
	chomp $line;
	my @arr = split /\t/,$line;
	next if $arr[2] eq "chr4";
	$div -> {$arr[2]} -> {$arr[3]} = 1;
	#$x_linked ++ if exists $sites -> {$arr[0]} -> {$arr[1]} && $arr[0] eq "X";
}
close DIV;


open TP, $site_type || die "open $site_type failed!\n";
my $sites;
while (my $line = <TP>) {
	chomp $line;
	next if $line eq "";
	my @arr = split /\t/,$line;
	#$arr[0] =~ s/^chr//sg;
	next if $arr[0] eq "chr4";
	if (exists $div -> {$arr[0]} -> {$arr[1]}){
		$sites -> {$arr[2]} ++;
	}
}
close TP;

print $site_type,"\n";
foreach my $type (keys %{$sites}){
	print $type,":\t",$sites -> {$type},"\n";
}

