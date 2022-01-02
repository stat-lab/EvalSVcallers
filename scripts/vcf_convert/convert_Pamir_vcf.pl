#!/usr/bin/perl -w
use strict;

my %vcf;

foreach my $file (@ARGV){
	open (FILE, $file) or die "$file is not found: $!\n";
	while (my $line = <FILE>){
		chomp $line;
		next if ($line =~ /^#|^$/);
		my @line = split (/\t/, $line);
		my $chr = $line[0];
		my $chr02d = $chr;
		$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
		my $pos = $line[1];
		my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
		my $len = $1 if ($line[7] =~ /SVLEN=(\d+)/);;
		my $read = $1 if ($line[7] =~ /Support=(\d+)/);
		${${$vcf{$chr02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$read\n";
	}
	close (FILE);
}

foreach my $chr (sort keys %vcf){
	foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
		foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
			print ${${$vcf{$chr}}{$pos}}{$type};
		}
	}
}
