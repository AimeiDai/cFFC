#!/usr/bin/perl -w

use strict;

die "Usage: perl $0 IDRurl metadata outdir download.list \n" if scalar @ARGV != 4;

open FILE, $ARGV[1] || die "open $ARGV[1] failed!\n";
<FILE>;
<FILE>;
my $gene2file;
while(my $line = <FILE>){
	chomp $line;
	my @arr = split "\t",$line;
	my (@files) = ($arr[15] =~ /\/files\/(.*?)\//g);
	foreach my $file (@files){
		$gene2file -> {$file} = $arr[7];
	}
}


open FILE, $ARGV[0] || die "open $ARGV[0] failed!\n";
my $files;
<FILE>;
while (my $line = <FILE>) {
	# body...
	chomp $line;
	my ($file_name) = ($line =~ /.*\/(.*?\.gz)/);
	$files -> {$file_name} = $line;
}


opendir DIR, $ARGV[2] || die "open $ARGV[2] failed!\n";
my $downloaded;
while (my $file_name = readdir DIR){
	next if $file_name eq "." || $file_name eq "..";
	my ($original) = ($file_name =~ /.*?_(.*?\.bed\.gz)/);
	if (! exists $files -> {$original} ){
		system "rm ".$ARGV[2]."/"."\"".$file_name."\"";
		#print "$ARGV[1]/$file_name\n";
	}else{
		$downloaded -> {$original} = 1;
	}
}

open OUT,">$ARGV[3]" || die "open $ARGV[3] failed!\n";
my $ix;
foreach my $file_name (keys %{$files}){
	if (! exists $downloaded -> {$file_name}){
		my ($sur) = ($file_name =~ /(.*?)\.bed\.gz/);
		my $shell = "wget -c ".$files -> {$file_name}." -O "."\"".$ARGV[2]."/".$gene2file -> {$sur}."_".$file_name."\"";
		print OUT $shell," &\n";
		$ix ++;
		if ($ix % 20 == 0){
			print OUT "wait\n";
			$ix=0;
		}
		#system $shell;
	}
}


