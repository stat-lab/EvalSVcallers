#!/usr/bin/perl -w
use strict;

my %vcf;

foreach my $var_file (@ARGV){
    if (!-e $var_file){
	print STDERR "$var_file is not found\n";
	next;
    }
    open (FILE, $var_file);
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#/){
    #	print $line, "\n";
	    next;
	}	
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	my $pos = $line[1];
	my $type = 'INS';
	my $len = 0;
	$len = $1 if ($line[4] =~ /length_(\d+)/);
	my $reads = 0;
	$reads = $1 if ($line[7] =~ /RP=(\d+);/);
	my $gt = '';
	$gt = $1 if (@line > 9) and ($line[9] =~ /^(.+?):/);
	my $chr_02d = $chr;
	$chr_02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	${$vcf{$chr_02d}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads;GT=$gt" if ($gt ne '');
	${$vcf{$chr_02d}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads" if ($gt eq '');
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	print ${$vcf{$chr}}{$pos}, "\n";
    }
}
