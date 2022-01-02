#!/usr/bin/perl -w
use strict;

my $vcf_file = shift @ARGV;

open (FILE, $vcf_file) or die "$vcf_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#/);
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?)[:;]/);
    my $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
    next if ($type eq 'INS') or ($type eq 'INV');
    my $gt = './.';
    $gt = $1 if ($line[9] =~ /^([01]\/[01])/);
    $gt = './.' if ($gt eq '0/0');
    print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=5;GT=$gt\n";
}
close (FILE);

