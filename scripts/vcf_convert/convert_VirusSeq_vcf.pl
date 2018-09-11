#!/usr/bin/perl -w
use strict;

# covert VirusSeq output file (*.integration-sites.txt) to vcf

my $file = shift @ARGV;

my %vcf;

open (FILE, $file) or die "$file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^Viral_Transcript/);
    my @line = split (/\s+/, $line);
    my $chr = $line[-2];
    my $chr02d = $chr;
    $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
    my $pos = $line[-1];
    my $disc_read = $line[2];
    my $read = 2;
    if ($disc_read <= 2){
        $read = 2;
    }
    elsif ($disc_read <= 4){
        $read = 3;
    }
    elsif ($disc_read <= 7){
        $read = 4;
    }
    elsif ($disc_read <= 10){
        $read = 5;
    }
    elsif ($disc_read <= 15){
        $read = 6;
    }
    elsif ($disc_read <= 20){
        $read = 7;
    }
    elsif ($disc_read <= 30){
        $read = 8;
    }
    elsif ($disc_read <= 40){
        $read = 9;
    }
    else{
        $read = 10;
    }
    ${$vcf{$chr02d}}{$pos} = "$chr\t$pos\tINS\t.\t.\t.\tPASS\tSVTYPE=INS;SVLEN=0;READS=$read\n";
}
close (FILE);

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
        print ${$vcf{$chr}}{$pos};
    }
}
