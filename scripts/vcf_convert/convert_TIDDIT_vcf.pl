#!/usr/bin/perl -w
use strict;

my $file = shift @ARGV;

open (FILE, $file) or die "$file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if (($line =~ /^#/) or ($line =~ /^sample\tchr/)){
	next;
    }
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    $type = 'DUP' if ($type eq 'TDUP');
	next if ($type !~ /DEL|^DUP$|INV/);
	my $end = $1 if ($line[7] =~ /END=(\d+)/);
	my $len = $end - $pos + 1;
	my $reads = 0;
    $reads = $1 if ($line[7] =~ /LTE=(\d+)/);
	$reads = 10 if ($reads == 0);
	my $gt = './.';
	$gt = $1 if ($line[9] =~ /^(.+?):/);
	$gt = '0/1' if ($gt eq './1');
next if ($line[6] ne 'PASS');   
    print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads;GT=$gt\n";
}
close (FILE);
