#!/usr/bin/perl -w
use strict;

my %vcf;

foreach my $file (@ARGV){
    open (FILE, $file) or die "$file is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#/){
	    next;
	}
	my @line = split (/\t/, $line);
	next if ($line[3] =~ /NB=Y/);
	my $chr = $line[0];
	my $pos = $line[1];
	my $end = $line[2];
	my $type = $1 if ($line[3] =~ /NM=(.+?);/);
	my $len = $end - $pos + 1;
	my $reads = 0;
	$reads = $1 if ($line[3] =~ /SR=(\d+),/);
	${$vcf{$chr}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=MEI;SVLEN=$len;READS=$reads\n";
    }
    close (FILE);
}

foreach my $chr (keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	print ${$vcf{$chr}}{$pos};
    }
}
