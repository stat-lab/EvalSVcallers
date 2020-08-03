#!/usr/bin/perl -w
use strict;

# covert Delly output files to vcf

my $var_file = shift @ARGV;

my $min_sv_len = 30;

my $min_reads = 2;

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
    next if ($len < $min_sv_len) and ($len > 0);
    my $end = $pos + $len - 1;
    my $reads = 0;
    $reads = $1 if ($line[-1] =~ /:(\d+)$/);
    my $chr_02d = $chr;
    $chr_02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
    next if ($chr !~ /^chr/) and ($chr !~ /^\d+$|[XY]/);
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
