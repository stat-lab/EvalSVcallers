#!/usr/bin/perl -w
use strict;

my $file = shift @ARGV;

my %vcf;

open (FILE, $file) or die "$file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^GeneName/);
    my @line = split (/\t/, $line);
    my $chr = $line[2];
    my $chr02d = $chr;
    $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
    my $pos = $line[3];
    my $len = $line[6] - $line[5] + 1;
    my $score = $line[-1];
    my $read = 2;
    if ($score < -1){
        $read = 2;
    }
    elsif ($score < 0){
        $read = 3;
    }
    elsif ($score < 2){
        $read = 4;
    }
    elsif ($score < 5){
        $read = 5;
    }
    elsif ($score < 10){
        $read = 6;
    }
    elsif ($score < 20){
        $read = 7;
    }
    elsif ($score < 50){
        $read = 8;
    }
    else{
        $read = 9;
    }
    ${$vcf{$chr02d}}{$pos} = "$chr\t$pos\tINS\t.\t.\t.\tPASS\tSVTYPE=INS;SVLEN=$len;READS=$read";
}
close (FILE);

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
        print ${$vcf{$chr}}{$pos}, "\n";
    }
}
