#!/usr/bin/perl -w
use strict;

my $var_file = shift @ARGV;

my %vcf;

open (FILE, $var_file) or die "$var_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
	next;
    }
    my @line = split (/\s+/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /EVENT=([A-Z]+)/);
    $type = 'DEL' if ($type eq 'DELETION');
    $type = 'DUP' if ($type eq 'DUPLICATION');
    $type = 'TRA' if ($type eq 'TRANSLOCATION');
    $type = 'INV' if ($type eq 'INVERSION');
$type = 'INS' if ($type eq 'INSERTION');
    my $end = 0;
    my $len = 0;
    my $pos2 = 0;
    my $chr2 = '';
    my $reads = 0;
    my $start = 0;
    my @item = split (/:/, $line[9]);
    $reads = $item[1];
    if ($line[4] =~ /[\[\]](\d+):(\d+)[\[\]]/){
	$chr2 = $1;
	$pos2 = $2;
    }
    if ($type ne 'TRA'){
	$start = $pos if ($pos <= $pos2);
	$start = $pos2 if ($pos2 < $pos);
	$len = abs ($pos - $pos2 + 1);
    }
    my $chr_02d = $chr;
    $chr_02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
    ${${$vcf{$chr_02d}}{$start}}{$type} = "$chr\t$start\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads" if ($type ne 'TRA');
    ${${$vcf{$chr_02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads\tCHR2=$chr2;POS2=$pos2" if ($type eq 'TRA');
}


foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
	    print ${${$vcf{$chr}}{$pos}}{$type}, "\n";
	}
    }
}
