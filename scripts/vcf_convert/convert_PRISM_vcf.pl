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
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	my $type = $line[4];
	my $bp1 = $line[1];
	my $bp2 = $line[2];
	my $svsize = $line[3];
	my $reads = $line[5];
	if (($svsize >= $min_sv_len) and ($reads >= $min_reads)){
	    my $chr02d = $chr;
	    $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	    ${${$vcf{$chr02d}}{$bp1}}{$type} = "SVTYPE=$type;SVLEN=$svsize;BP2=$bp2;READS=$reads";
	}
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
	    my $chr2 = $chr;
	    $chr2 =~ s/^0*//;
	    print "$chr2\t$pos\t$type\t.\t.\t.\tPASS\t${${$vcf{$chr}}{$pos}}{$type}\n";
	}
    }
}
