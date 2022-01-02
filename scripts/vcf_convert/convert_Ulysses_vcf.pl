#!/usr/bin/perl -w
use strict;

my %vcf;
my $count = 0;

foreach my $file (@ARGV){
    if (!-f $file){
	print STDERR "$file is not found: \n";
	next;
    }
    open (FILE, $file);
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#/){
    #	    print $line, "\n" if ($count == 1);
	    next;
	}	
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	my $pos = $line[1];
	my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
	my $len = 0;
	$type = 'INS' if ($type eq 'sINS');
    #    $type = 'TRA' if ($type eq 'RT') or ($type eq 'NRT');
	next if ($type eq 'RT') or ($type eq 'NRT');
	my $len1 = $1 if ($line[7] =~ /SVMI=-*(\d+);/);
	my $len2 = $1 if ($line[7] =~ /SVM=-*(\d+);/);
	$len = int (($len1 + $len2) / 2);
	my $end = $pos + $len - 1;
	my $pvalue = 0;
	$pvalue = $1 if ($line[7] =~ /PVAL=[\d\.]+e-(\d+)/);
	$pvalue =~ s/^0*// if ($pvalue =~ /^0\d+/);
	my $reads = $pvalue + 3;
	my $chr_02d = $chr;
	$chr_02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	${${$vcf{$chr_02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads";
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
