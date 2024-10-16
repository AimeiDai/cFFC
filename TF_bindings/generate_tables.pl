#!/usr/bin/perl -w
use strict;

die "Usage: perl $0 metatab gtf_file outtab!\n" if scalar @ARGV != 3;

open GTF, $ARGV[1] || die "open $ARGV[1] failed!\n";
my $symbol2fbgn;
while (my $line = <GTF>) {
	# body...
	next if $line =~ /^#/sg;
	my @arr = split "\t",$line;
	next if $arr[2] ne "gene";
	my ($gene_id) = ($arr[8] =~ /gene_id \"(.*?)\"/sg);
	my ($gene_name) = ($arr[8] =~ /gene_name \"(.*?)\"/sg);
	next if ! defined $gene_name;
	$symbol2fbgn -> {$gene_name} = $gene_id;
}

open OUT, ">".$ARGV[2] || die "open $ARGV[2] failed!\n";
print OUT "geneSymbol,FBgnID,Accession,GEOacc,stage,BiologicalReplicate,Lab\n";
open META, $ARGV[0] || die "open $ARGV[0] failed!\n";
<META>;
<META>;
while (my $line = <META>){
	my @arr = split "\t",$line;
	my ($geo_acc) = ($arr[10] =~ /(GEO:\w+)/sg); 
	$arr[18] =~ s/,/|/g;
	print OUT $arr[7],",",$symbol2fbgn -> {$arr[7]},",",$arr[1],",",$geo_acc,",",$arr[22]," ",$arr[23],",",$arr[18],",",$arr[12],"\n";
}

