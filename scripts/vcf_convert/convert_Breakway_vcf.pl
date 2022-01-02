#!/usr/bin/perl -w
use strict;

my %vcf;

foreach my $var_file (@ARGV){
    open (FILE, $var_file) or die "$var_file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        if ($line =~ /^chrom1:/){
            next;
        }	
        my @line = split (/\t/, $line);
        my $chr = '';
        my $pos = 0;
        my $end = 0;
        my $bp1_s = 0;
        my $bp1_e = 0;
        my $bp2_s = 0;
        my $bp2_e = 0;
        if ($line[0] =~ /(.+?):(\d+)-(\d+)/){
            $chr = $1;
            $bp1_s = $2;
            $bp1_e = $3;
            $pos = int (($bp1_s + $bp1_e) / 2);
        }
        if ($line[1] =~ /(.+?):(\d+)-(\d+)/){
            $bp2_s = $2;
            $bp2_e = $3;
            $end = int (($bp2_s + $bp2_e) / 2);
        }
        if ($pos > $end){
            my $end2 = $end;
            $end = $pos;
            $pos = $end2;
        }
        my $type = $line[2];
        $pos = int (($end + $pos) / 2) if ($type eq 'INS');
        my $len = 0;
        $len = $end - $pos + 1 if ($type ne 'INS');
        my $reads = $line[3];
        my $chr02d = $chr;
        $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
        ${$vcf{$chr02d}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads";
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
        print "${$vcf{$chr}}{$pos}\n";
    }
}
