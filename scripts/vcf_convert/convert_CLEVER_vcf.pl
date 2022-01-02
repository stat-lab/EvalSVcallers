#!/usr/bin/perl -w
use strict;

my $var_file = shift @ARGV;

my %vcf;

open (FILE, $var_file) or die "$var_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
	print $line, "\n";
	next;
    }	
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    my $len = 0;
    $len = $1 if ($line[7] =~ /SVLEN=-*(\d+);/);
    my $end = $pos + $len - 1;
    my $reads = 0;
    $reads = $1 if ($line[-1] =~ /:(\d+)$/);
    my $chr_02d = $chr;
    $chr_02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
    ${${$vcf{$chr_02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads";
}
close (FILE);

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
	    print ${${$vcf{$chr}}{$pos}}{$type}, "\n";
	}
    }
}
