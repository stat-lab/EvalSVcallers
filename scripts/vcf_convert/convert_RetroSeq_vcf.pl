#!/usr/bin/perl -w
use strict;

my $file = shift @ARGV;

my %vcf;

open (FILE, $file) or die "$file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
	next;
    }
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /MEINFO=(.+?),/);
    my $len = 0;
    my $reads = 3;
    my @info = split (/:/, $line[9]);
    $reads = $info[2] if ($info[2] > 3);
    my $gt = $1 if ($line[9] =~ /^(.+?):/);
    ${$vcf{$chr}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=MEI;SVLEN=$len;READS=$reads;GT=$gt";
}
close (FILE);

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	print ${$vcf{$chr}}{$pos}, "\n";
    }
}
