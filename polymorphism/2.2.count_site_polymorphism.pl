#!/usr/bin/perl -w
use strict;

die "Usage: perl $0 fourfold.site snp.file pop_size:30 out\n" if scalar @ARGV != 4;
my ($fourfold, $snp_file, $pop_size, $out) = @ARGV;

open SNP ,$snp_file || die "open $snp_file failed!\n";
my $snp;
<SNP>;
while (my $line = <SNP>) {
	# body...
	chomp $line;
	my @arr = split /\t/,$line;
	next if $arr[1] eq "chr4" || $arr[5] eq "u" || $arr[8] == 0 || $arr[8] == 1 || $arr[6]/$arr[8] < $pop_size;
	$snp -> {$arr[1]} -> {$arr[2]} -> {"AC"} = $arr[6];
	$snp -> {$arr[1]} -> {$arr[2]} -> {"MLAC"} = $arr[7];
	# next if $arr[3] + $arr[5] < $pop_size*2 || $arr[3] == 0 || $arr[5] == 0 || $arr[0] eq "4";

	# my $daf = $arr[3]/($arr[3] + $arr[5]);
	# my $maf = $daf > 1 - $daf ? 1-$daf : $daf;

	# if (exists $sites -> {$arr[0]} -> {$arr[1]}){
	# 	$counts ++;
	# 	$maf5 ++ if $maf >= 0.05;
		
	# }
	#$neutral_MAF5 ++ if $daf > 0.05;
}
close SNP;

open OUT, ">$out" || die "open OUT failed!\n";
open FOLD, $fourfold || die "open $fourfold failed!\n";
my ($sites, $sites_MAF5);
while (my $line = <FOLD>) {
	chomp $line;
	next if $line eq "";
	my @arr = split /\t/,$line;
	#$arr[0] =~ s/^chr//sg;
	next if $arr[0] eq "chr4";
	next if ! $snp -> {$arr[0]} -> {$arr[1]};
	#my @freq = split "\t", $snp -> {$arr[0]} -> {$arr[1]};
	#next if $freq[1] + $freq[3] < $pop_size || $freq[1] == 0 || $freq[3] == 0;
	print OUT $line,"\t",$snp -> {$arr[0]} -> {$arr[1]} -> {"AC"},"\n";
	$sites -> {$arr[2]} ++;

	# my $daf = $freq[1]/($freq[1] + $freq[3]);
	# my $maf = $daf > 1 - $daf ? 1-$daf : $daf;
	# $sites_MAF5 -> {$arr[2]} ++ if $maf >= 0.05;
	$sites_MAF5 -> {$arr[2]} ++ if $snp -> {$arr[0]} -> {$arr[1]} -> {"MLAC"} > 2;
}

close FOLD;
close OUT;

print $fourfold,"\n";
foreach my $type (keys %{$sites}){
	print $type,"\n";
	print "total polymorphism: ",$sites -> {$type},"\n";
	print "polymorphism with MAF>=5%: ", $sites_MAF5 -> {$type},"\n";
}
print "\n";
