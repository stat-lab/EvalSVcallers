#!/usr/bin/perl -w
use strict;

my $var_file = shift @ARGV;

my %vcf;
my $count = 0;

open (FILE, $var_file) or die "$var_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#|^$|^\s+/){
        next;
    }
    my @line = split (/\s+/, $line);
    my $chr = $line[0];
    my $chr02d = $chr;
    $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    next if ($type eq 'COMPLEX');
    my $end = $1 if ($line[7] =~ /END=(\d+);/);
    my $mintipq = $1 if ($line[7] =~ /MINTIPQ=(\d+)/);
    my $reads = 3;
    $reads = 5 if ($line[7] =~ /^PRECISE/);
    my $gt = './.';
    if ($type eq 'DEL'){
        $gt = '0/1' if ($line[9] eq '2:1');
        $gt = '1/1' if ($line[9] eq '2:0');
    }
    elsif ($type eq 'DUP'){
        my $cr = 0;
        if ($line[9] =~ /(\d+):(\d+)$/){
            $cr = $2 / $1;
        }
        if ($cr < 2){
            $gt = '0/1';
        }
        else{
            $gt = '1/1';
        }
    }
    $line[7] .= ";READS=$reads;GT=$gt";
    my $new_line = join ("\t", @line);
    ${$vcf{$chr02d}}{$pos} = $new_line;
}
close (FILE);

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	    print ${$vcf{$chr}}{$pos}, "\n";
    }
}
