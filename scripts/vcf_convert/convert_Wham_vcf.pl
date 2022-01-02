#!/usr/bin/perl -w
use strict;

my $file = shift @ARGV;

open (FILE, $file) or die "$file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
        print "$line\n";
        next;
    }
    my @line = split (/\s+/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $type = '';
    $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    next if ($type eq 'BND');
    my $len = 0;
    $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
    my $read = 0;
    $read = $1 if ($line[7] =~ /A=(\d+);/);
    my $cw = $1 if ($line[7] =~ /CW=(.+?);/);
    my @cw = split (/,/, $cw);
    next if ($cw[4] > 0.2);
    if ($type eq 'DEL'){
        next if ($cw[0] < 0.2);
    }
    elsif ($type eq 'DUP'){
        next if ($cw[1] < 0.2);
    }
    elsif ($type eq 'INV'){
        next if ($cw[2] < 0.2);
    }
    elsif ($type eq 'INS'){
        next if ($cw[3] < 0.2);
    }
    print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$read\n";
}
close (FILE);
