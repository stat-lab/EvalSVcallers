#!/usr/bin/perl -w
use strict;

my %vcf;

foreach my $var_file (@ARGV){
    open (FILE, $var_file) or die "$var_file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        if ($line =~ /^#/){
        	next;
        }	
        my @line = split (/\t/, $line);
        my $chr = $line[0];
        my $chr02d = $chr;
	$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
        my $bp1_s = $line[1];
        my $bp1_e = $line[2];
        my $bp2_s = $line[4];
        my $bp2_e = $line[5];
        my $pos = int (($bp1_e + $bp1_s) / 2);
        my $end = int (($bp2_e + $bp2_s) / 2);
        my $type = 'DEL';
        my $len = $end - $pos;
        $len = abs ($len);
        $type = 'INS' if ($len < 10);
        my $reads = $line[7];
        ${${$vcf{$chr02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads";
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
	    print ${${$vcf{$chr}}{$pos}}{$type}, "\n";
	}
    }
}
