#!/usr/bin/perl -w
use strict;

# covert Wham output files to vcf

my $file = shift @ARGV;

my $min_del_size = 30;
my $min_sv_size = 50;
my $max_sv_size = 2000000;

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
    next if ($len > $max_sv_size);
    next if ($type eq 'DEL') and ($len < $min_del_size);
    next if ($type ne 'DEL') and ($len < $min_sv_size);
    my $read = 0;
    $read = $1 if ($line[7] =~ /A=(\d+);/);
    next if ($read < 3);
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
