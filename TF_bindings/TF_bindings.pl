#!/usr/bin/perl -w

use strict;

die "Usage: perl $0 ChIP.bed STable gtf_file output!\n" if scalar @ARGV != 4;

my $gene;  #coordinates of genes
open GTF,$ARGV[2]||die "open $ARGV[1] failed!!\n";
while(my $line = <GTF>){
	chomp $line;
	next if $line =~/^#/sg || $line eq "";
	my @arr = split "\t",$line;
	next if $arr[2] ne "gene";
	my ($gene_id) = $arr[8] =~ /gene_id \"(FBgn\d+)\"/sg;
	next if ! defined $gene_id;
	my ($start,$end) = ($arr[3],$arr[4]);
	if ($arr[6] eq "+"){
		$start -= 1500;
		$end += 500;
	}
	if ($arr[6] eq "-"){
		$start -= 500;
		$end += 1500;
	}
	$gene -> {$arr[0]} -> {$start} -> {$end} = $gene_id;
}

my $tfs;
my $cgids;
open STAB, $ARGV[1] || die "open $ARGV[1] failed!\n";
<STAB>;
while (my $line = <STAB>){
	chomp $line;
	my @arr = split ",",$line;
	my ($cgid,$geneSymbol) = (split /\|/,$arr[0]);
	if (defined $geneSymbol){
		$tfs -> {$geneSymbol} = $arr[1];
		$cgids -> {$cgid} = $geneSymbol;
	}else{
		$tfs -> {$cgid} = $arr[1];
	}
}

#my $tf_bindings;
#opendir DIR, $ARGV[0] || die "open $ARGV[0] failed!\n";
# foreach my $bed_file (readdir DIR){
# 	next if $bed_file eq "." || $bed_file eq "..";
	my ($bed_file)= ($ARGV[0] =~ /\/(.*$)/);
	my ($tf) = ($bed_file =~ /(.*?)_/);
	open BED, $ARGV[0];
	open OUT, ">".$ARGV[3]."/".$bed_file;
	while (my $line = <BED>){
		chomp $line;
		my @arr = split "\t",$line;
		$arr[0] =~ s/^chr//;
		foreach my $start (sort{$a<=>$b} keys %{$gene -> {$arr[0]}}){
			foreach my $end (sort{$a<=>$b} keys %{$gene -> {$arr[0]} -> {$start}}){
				my $l = $arr[1]>$start?$arr[1]:$start;
				my $s = $arr[2]<$end?$arr[2]:$end;
				if ($s - $l >= 0){
					if (! exists $tfs -> {$tf}){
						$tf = $cgids -> {$tf};
					}
					print OUT $tf,"\t",$tfs->{$tf},"\t",$gene -> {$arr[0]} -> {$start} -> {$end},"\n";
					#$tf_bindings -> {$tfs -> {$tf}} -> {$gene -> {$arr[0]} -> {$start} -> {$end}} = $tf;
				}
			}
		}
	}
	close OUT;
	close BED;
#}

# open OUT,">$ARGV[3]" || die "open $ARGV[3] failed!\n";
# foreach my $tf (keys %{$tf_bindings}){
# 	foreach my $tg (keys %{$tf_bindings -> {$tf}}){
# 		print OUT $tf,"\t",$tg,"\t",$tf_bindings->{$tf}->{$tg},"\n";
# 	}
# }
