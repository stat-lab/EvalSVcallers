#!/usr/bin/perl -w
use strict;

my $file = shift @ARGV;

open (FILE, $file) or die "$file is not found: $!\n";
while (my $line = <FILE>){
	chomp $line;
	next if ($line =~ /^#|^$/);
	my @line = split (/\s+/, $line);
	my $chr = $line[0];
	my $pos = $line[1];
	my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
	$type = 'DEL' if ($type eq 'del');
	$type = 'DUP' if ($type eq 'tandup');
	$type = 'INV' if ($type eq 'inv');
	my $end = $1 if ($line[7] =~ /END=(\d+)/);
	my $len = $end - $pos + 1;
	my $cn = '';
	$cn = $1 if ($line[7] =~ /CN=(\d+)/);
	my $gt = $line[9];
	my $read = 3;
	if ($line[5] == 100){
		$read = 9;
	}
	elsif ($line[5] >= 80){
		$read = 8;
	}
	elsif ($line[5] >= 60){
		$read = 7;
	}
	elsif ($line[5] >= 40){
		$read = 6;
	}
	elsif ($line[5] >= 20){
		$read = 5;
	}
	elsif ($line[5] >= 1){
		$read = 4;
	}
	print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$read;GT=$gt;CN=$cn\n" if ($cn ne '');
	print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$read;GT=$gt\n" if ($cn eq '');
	
}
close (FILE);
