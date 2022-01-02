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
    my $end = $line[2];
    my $type = $line[3];
    my $len = $line[4];
    $len = 0 - $len if ($len < 0);
    my $RP_info = $line[8];
    my $SR_info = $line[10];
    my @RP_info = split (/;/, $RP_info);
    my $RP_num = $RP_info[2];
    $RP_num = $1 if ($RP_num =~ /^-*(\d+):/);
    $RP_num = 0 if ($RP_num < 0);
    my @SR_info = split (/;/, $SR_info);
    my $SR_num = $SR_info[1];
    $SR_num = 0 if ($SR_num < 0);
    my $signal = $RP_num + $SR_num;
    my $reads = 2;
    if ($signal <= 2){
        $reads = 2;
    }
    elsif ($signal <= 10){
        $reads = 3;
    }
    elsif ($signal <= 20){
        $reads = 4;
    }
    elsif ($signal <= 30){
        $reads = 5;
    }
    elsif ($signal <= 50){
        $reads = 6;
    }
    elsif ($signal <= 80){
        $reads = 7;
    }
    elsif ($signal <= 130){
        $reads = 8;
    }
    elsif ($signal <= 200){
        $reads = 9;
    }
    else{
        $reads = 10;
    }
    print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads\n";
}
close (FILE);
