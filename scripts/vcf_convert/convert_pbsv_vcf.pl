#!/usr/bin/perl -w
use strict;

my $file = shift @ARGV;

open (FILE, $file) or die "$file is not found: $!\n";
while (my $line = <FILE>){
	chomp $line;
	next if ($line =~ /^#|^$/);
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	my $pos = $line[1];
	my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
	my $len = 0;
	$len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
	if ($line[7] !~ /SVLEN=/){
		my $end = $1 if ($line[7] =~ /END=(\d+)/);
		$len = $end - $pos + 1;
	}
	my $read = 3;
	my $dp = 3;
	my $gt = './.';
	$gt = $1 if ($line[9] =~ /^([01\.]\/[01\.]):/);
	if ($line[9] =~ /^[01\.]\/[01\.]:\d+,(\d+):/){
		$read = $1;
	}
	elsif ($line[9] =~ /^[01\.]\/[01\.]:(\d+):/){
		$read = $1;
	}
	if ($line[9] =~ /^[01\.]\/[01\.]:\d+,\d+:(\d+)/){
		$dp = $1;
	}
	elsif ($line[9] =~ /^[01\.]\/[01\.]:\d+:(\d+)/){
		$dp = $1;
	}
	my $vrr = 0;
	$vrr = int ($read / $dp * 100 + 0.5) / 100 if ($dp > 0);
	print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$read;VRR=$vrr;GT=$gt\n";
}
close (FILE);
