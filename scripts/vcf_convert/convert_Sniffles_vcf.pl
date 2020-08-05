#!/usr/bin/perl -w
use strict;

my $var_file = shift @ARGV;

open (FILE, $var_file) or die "$var_file is not found: \n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#/);
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    next if ($type eq 'TRA');
    my $len = 0;
    $len = $1 if ($line[7] =~ /SVLEN=(\d+)/);
    my $gt = $1 if ($line[9] =~ /^([01\.]\/[01\.])/);
    my $reads = $1 if ($line[9] =~ /:(\d+)$/);
    print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads;GT=$gt\n";
}
close (FILE);
