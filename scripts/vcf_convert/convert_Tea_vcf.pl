#!/usr/bin/perl -w
use strict;

# covert Mobster output files to vcf

my $file = shift @ARGV;

open (FILE, $file) or die "$file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if (($line =~ /^#/) or ($line =~ /^sample\tchr/)){
	next;
    }
    my @line = split (/\t/, $line);
    my $chr = $line[1];
    $chr =~ s/^chr// if ($chr =~ /^chr/);
    my $pos = $line[6];
    my $type = $line[9];
    $type = 'ALU' if ($type eq 'Alu');
    my $len = $line[4];
    my $reads = $line[11];
    next if ($chr !~ /^\d+$|[XY]/);
    print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=MEI;SVLEN=$len;READS=$reads\n";
}
close (FILE);
