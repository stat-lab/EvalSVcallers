#!/usr/bin/perl -w
use strict;

my %vcf;

foreach my $file (@ARGV){
    if (!-e $file){
	print STDERR "$file is not found\n";
	next;
    }
    open (FILE, $file);
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#/){
	    next;
	}
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	my $pos = $line[1];
	my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
	$type = 'NUMT' if ($type eq 'MT');
	my $len = 0;
	$len = $1 if ($line[7] =~ /SVLEN=(\d+)/);
	my $reads = 0;
	my $read1 = $1 if ($line[7] =~ /LP=(\d+)/);
	my $read2 = $1 if ($line[7] =~ /RP=(\d+)/);
	my $gt = $1 if ($line[9] =~ /^(.+?):/);
	$reads = int (($read1 + $read2) * 0.5);
	next if ($type eq 'SVA') and ($reads < 8);
	my $class = $type;
	if (($type ne 'NUMT') and ($type ne 'VEI')){
	    $class = 'MEI';
	}
	my $chr02d = $chr;
	$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	${$vcf{$chr02d}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$class;SVLEN=$len;READS=$reads;GT=$gt\n";
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	print ${$vcf{$chr}}{$pos};
    }
}
