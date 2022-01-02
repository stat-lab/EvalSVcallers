#!/usr/bin/perl -w
use strict;

my %vcf;

my $count = 0;

foreach my $var_file (@ARGV){
    if (!-f $var_file){
	print STDERR "$var_file is not found: \n";
	next;
    }
    $count ++;
    open (FILE, $var_file) or die "$var_file is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#/){
	    print "$line\n" if ($count == 1);
	    next;
	}
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	my $start = $line[1];
	my $end = $line[4];
	my $chr2 = $line[3];
	my $type = $line[6];
	my $svsize = $line[7];
	$svsize = 0 - $svsize if ($svsize < 0);
	my $bpsize = $end - $start + 1;
	my $reads = $line[9];
	if ($type eq 'ITX'){
	    $type = 'DEL' if ($line[7] > 0);
	    $type = 'INS' if ($line[7] < 0);
	}
	my $chr_02d = $chr;
	$chr_02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	if ($type ne 'CTX'){
	    ${$vcf{$chr_02d}}{$start} = "$chr\t$start\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$svsize;BPSIZE=$bpsize;READS=$reads\n";
    #	print "$chr\t$start\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$svsize;BPSIZE=$bpsize;READS=$reads\n";
	}
	else{
	    $type = 'TRA';
	    ${$vcf{$chr_02d}}{$start} = "$chr\t$start\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$svsize;READS=$reads\tCHR2=$chr2;POS2=$end\n";
    #	print "$chr\t$start\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$svsize;READS=$reads\tTR-chr=$chr2;TR-pos=$end\n";
	}
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	print ${$vcf{$chr}}{$pos};
    }
}
