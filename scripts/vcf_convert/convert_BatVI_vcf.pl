#!/usr/bin/perl -w
use strict;

my %vcf;

foreach my $file (@ARGV){
    open (FILE, $file) or die "$file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^LIB/);
        my @line = split (/\t/, $line);
        my $chr = $line[1];
        my $chr02d = $chr;
        $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
        my $pos = $line[2];
        my $split_read = 2;
        $split_read = $line[7] if (@line > 7);
        next if ($line[5] == 0);
        my $read = 1;
        if ($split_read <= 1){
            $read = 1;
        }
        elsif ($split_read <= 2){
            $read = 2;
        }
        elsif ($split_read <= 5){
            $read = 3;
        }
        elsif ($split_read <= 10){
            $read = 4;
        }
        elsif ($split_read <= 20){
            $read = 5;
        }
        elsif ($split_read <= 30){
            $read = 6;
        }
        elsif ($split_read <= 50){
            $read = 7;
        }
        elsif ($split_read <= 80){
            $read = 8;
        }
        else{
            $read = 9;
        }
        ${$vcf{$chr02d}}{$pos} = "$chr\t$pos\tINS\t.\t.\t.\tPASS\tSVTYPE=INS;SVLEN=0;READS=$read\n";
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
        print ${$vcf{$chr}}{$pos};
    }
}
