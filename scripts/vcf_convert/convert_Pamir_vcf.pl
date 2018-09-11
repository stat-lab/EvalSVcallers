#!/usr/bin/perl -w
use strict;

# covert pamir output files (insertions_setcover.vcf) to vcf

my $file = shift @ARGV;

my $min_len = 30;

my %vcf;

open (FILE, $file) or die "$file is not found: $!\n";
while (my $line = <FILE>){
	chomp $line;
	next if ($line =~ /^#|^$/);
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	my $pos = $line[1];
	my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
	my $len = $1 if ($line[7] =~ /SVLEN=(\d+)/);
	next if ($len < $min_len);
	my $read = $1 if ($line[7] =~ /Support=(\d+)/);
	print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$read\n";
}
close (FILE);
