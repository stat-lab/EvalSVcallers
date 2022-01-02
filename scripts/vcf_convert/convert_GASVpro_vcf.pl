#!/usr/bin/perl -w
use strict;

my $var_file = shift @ARGV;

open (FILE, $var_file) or die "$var_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
	print "$line\n";
	next;
    }
    my @line = split (/\t/, $line);
    my $id = $line[0];
    my $chr1 = $line[1];
    my ($start1, $end1) = split (/,/, $line[2]);
    my $start = int (($end1 + $start1) / 2);
    my $chr2 = $line[3];
    my ($start2, $end2) = split (/,/, $line[4]);
    my $end = int (($end2 + $start2) / 2);
    my $svlen = $end - $start + 1;
    my $type = $line[7];
    if ($type eq 'D'){
	$type = 'DEL';
    }
    elsif (($type eq 'IR') or ($type eq 'I+') or ($type eq 'I-')){
	$type = 'INV';
    }
    elsif (($type eq 'TR') or ($type eq 'TN')){
	$type = 'TRA';
	$svlen = 0;
    }
    my $reads = $line[5];
    if ($type ne 'TRA'){
	print "$chr1\t$start1\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$svlen;READS=$reads\n";
    }
    else{
	print "$chr1\t$start\t$type\t.\t.\t.\tPASS\t;SVTYPE=$type;READS=$reads\tCHR2=$chr2;POS2=$end\n";
    }
}
close (FILE);
