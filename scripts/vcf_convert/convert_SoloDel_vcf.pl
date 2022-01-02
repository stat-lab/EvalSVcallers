#!/usr/bin/perl -w
use strict;

my %vcf;

foreach my $var_file (@ARGV){
    open (FILE, $var_file) or die "$var_file is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^Chr\t/){
	    next;
	}
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	my $pos = $line[1];
	my $type = 'DEL';
	my $end = $line[2];
	my $reads = $line[4];
	my $len = $line[3];
	my $chr_02d = $chr;
	$chr_02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	${$vcf{$chr_02d}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads\n";
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	print ${$vcf{$chr}}{$pos};
    }
}
