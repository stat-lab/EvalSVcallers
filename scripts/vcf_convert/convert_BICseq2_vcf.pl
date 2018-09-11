#!/usr/bin/perl -w
use strict;

# covert to vcf format

my @out_files = (@ARGV);

my $min_dup_ratio = 0.3;
my $max_del_ratio = -0.5;

foreach my $out_file (@out_files){
    open (FILE, $out_file) or die "$out_file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^chrom|^$/);
        my @line = split (/\t/, $line);
        my $chr = $line[0];
        my $start = $line[1];
        my $end = $line[2];
        my $log2ratio = $line[6];
        my $pvalue = $line[7];
        my $len = $end - $start + 1;
        next if ($log2ratio > $max_del_ratio) and ($log2ratio < $min_dup_ratio);
        my $type = 'DEL';
        $type = 'DUP' if ($log2ratio > 0);
        my $gt = '0/1';
        $gt = '1/1' if ($type eq 'DEL') and ($log2ratio < -2);
        $gt = '1/1' if ($type eq 'DUP') and ($log2ratio > 0.8);
        my $read = 3;
        $pvalue = 1e-100 if ($pvalue == 0);
        my $log_p = 0 - (log ($pvalue) / log (10));
        if ($log_p <= 1){
            $read = 3;
        }
        elsif ($log_p <= 2){
            $read = 4;
        }
        elsif ($log_p <= 3){
            $read = 5;
        }
        elsif ($log_p <= 5){
            $read = 6;
        }
        elsif ($log_p <= 10){
            $read = 7;
        }
        elsif ($log_p <= 30){
            $read = 8;
        }
        elsif ($log_p <= 60){
            $read = 9;
        }
        else{
            $read = 10;
        }
        print "$chr\t$start\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$read;GT=$gt\n";
    }
    close (FILE);
}
