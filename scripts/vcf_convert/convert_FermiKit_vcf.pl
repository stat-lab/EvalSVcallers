#!/usr/bin/perl -w
use strict;

my $var_file = shift @ARGV;

my %vcf;
my $count = 0;

open (FILE, $var_file) or die "$var_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
	print $line, "\n" if ($count == 1);
	next;
    }
    my @line = split (/\s+/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    next if ($type eq 'COMPLEX');
    my $end = $1 if ($line[7] =~ /END=(\d+);/);
    my $mintipq = $1 if ($line[7] =~ /MINTIPQ=(\d+)/);
    my $reads = 0;
    if ($mintipq >= 20){
	$reads = 15;
    }
    elsif ($mintipq >= 15){
	$reads = 10;
    }
    elsif ($mintipq >= 10){
	$reads = 7;
    }
    elsif ($mintipq >= 5){
	$reads = 5;
    }
    elsif ($mintipq >= 1){
	$reads = 3;
    }
    my $len = $end - $pos + 1;
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
