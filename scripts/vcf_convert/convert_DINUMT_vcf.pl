#!/usr/bin/perl -w
use strict;

my $var_file = shift @ARGV;

open (FILE, $var_file) or die "$var_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
	next;
    }
    my @line = split (/\s+/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $reads = 10;
    my $len = $1 if ($line[7] =~ /MLEN=(\d+)/);
    print "$chr\t$pos\tNUMT\t.\t.\t.\tPASS\tSVTYPE=INS;SVLEN=$len;READS=$reads\n";
}
close (FILE);
