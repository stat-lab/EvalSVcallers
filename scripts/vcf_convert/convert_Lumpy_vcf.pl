#!/usr/bin/perl -w
use strict;

my %vcf;
my $count = 0;

foreach my $file (@ARGV){
    if (!-e $file){
	print STDERR "$file is not found\n";
	next;
    }
    $count ++;
    open (FILE, $file) if ($file !~ /\.gz$/);
    open (FILE, "gzip -dc $file |") if ($file =~ /\.gz$/);
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#/){
#	    print $line, "\n" if ($count == 1);
	    next;
	}	
	my @line = split (/\t/, $line);
	    if (@line >= 11){
		if ($line[9] =~ /:0:0:0$/){
		    next;
		}
	    }
	my $chr = $line[0];
	my $pos = $line[1];
	my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
	next if ($type eq 'BND');
	my $len = 0;
	$len = $1 if ($line[7] =~ /SVLEN=-*(\d+);/);
	my $end = $pos + $len - 1;
	my $reads = 0;
	$reads = $1 if ($line[7] =~ /SU=(\d+);/);
	my $gt = '';
	$gt = $1 if ($line[-1] =~ /^(.+?):/);
	my $chr_02d = $chr;
	$chr_02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	${${$vcf{$chr_02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads;GT=$gt" if ($gt ne '');
	${${$vcf{$chr_02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads" if ($gt eq '');
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
