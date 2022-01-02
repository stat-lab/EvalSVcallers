#!/usr/bin/perl -w
use strict;

my $file = shift @ARGV;

open (FILE, $file) or die "$file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^Chr\s+/);
	my @line = split (/\s+/, $line);
	my $chr = $line[0];
	my $pos = $line[1];
	my $end = $line[2];
	my $len = $end - $pos + 1;
	my $type = $line[3];
	my $read = $line[6];
	$chr =~ s/^chr//;
	print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=MEI;SVLEN=$len;READS=$read\n";
}
close (FILE);
