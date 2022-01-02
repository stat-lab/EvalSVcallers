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
    open (FILE, $file);
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#/){
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
	my $type = $1 if ($line[7] =~ /SVTYPE=(.+)/);
    	$type = substr $type, 0, 3;
	next if ($type eq 'BND');
	my $len = 0;
	$len = $1 if ($line[7] =~ /SVLEN=-*(\d+);/);
	my $reads = 0;
	$reads = $1 if ($line[7] =~ /SE=(\d+);/);
	$reads += $1 if ($line[7] =~ /PE=(\d+);/);
	$reads += $1 if ($line[7] =~ /CE=(\d+);/);
	my $gt = '';
	print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads\n";
    }
    close (FILE);
}
