#!/usr/bin/perl -w
use strict;

# covert pbsv output files (*.vcf) to vcf

my $file = shift @ARGV;

open (FILE, $file) or die "$file is not found: $!\n";
while (my $line = <FILE>){
	chomp $line;
	next if ($line =~ /^#|^$/);
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	my $pos = $line[1];
	my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
	my $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
	my $read = 3;
	my $gt = './.';
	$gt = $1 if ($line[9] =~ /^([01]\/[01]):/);
	print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$read;GT=$gt\n";
}
close (FILE);
