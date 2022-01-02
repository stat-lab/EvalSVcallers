#!/usr/bin/perl -w
use strict;

foreach my $file (@ARGV){
    open (FILE, $file) or die "$file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        if ($line =~ /^#/){
            next;
        }
        my @line = split (/\s+/, $line);
        my $chr = $line[0];
        my $pos = $line[1];
        my $len = $1 if ($line[7] =~ /SVLEN=(\d+)/);
        my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
        $type = 'DUP' if ($type eq 'DUP:TANDEM');
        my $reads = 10;
        my $gt = './.';
        $gt = '0/1' if ($line[9] =~ /^1\/0/) or ($line[9] =~ /^1\/\./);
        $gt = '1/1' if ($line[9] =~ /^1\/1/);
        print "$chr\t$pos\t.\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads;GT=$gt\n";
    }
    close (FILE);
}
