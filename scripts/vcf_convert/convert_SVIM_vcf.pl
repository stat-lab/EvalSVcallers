#!/usr/bin/perl -w
use strict;

my $file = shift @ARGV;

open (FILE, $file) or die "$file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#/);
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?)[:;]/);
    next if ($type !~ /^DEL|^INS|^INV|^DUP/);
    my $len = 0;
    $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
    if ($len == 0){
        my $end = 0;
        $end = $1 if ($line[7] =~ /END=(\d+)/);
        $len = $end - $pos + 1;
    }
    my $reads = $1 if ($line[7] =~ /SUPPORT=(\d+)/);
    print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads\n";
}
close (FILE);
