#!/usr/bin/perl -w
use strict;

my %vcf;

foreach my $file (@ARGV){
    open (FILE, $file) or die "$file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        if (($line =~ /^Chromosome\s/) or ($line =~ /^chrVirus/)){
        	next;
        }
        my @line = split (/\s+/, $line);
        my $chr = $line[0];
        my $chr02d = $chr;
        $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
        my $pos = $line[1];
        my $len = 0;
        my $type = 'VEI';
        my $reads = $line[6];
        if ($reads =~ /\+/){
		my ($read_s, $read_p) = split (/\+/, $reads);
		if ($read_s >= 3){
		    $reads = $read_s;
		}
		else{
		    $reads = $read_p;
		}
        }
        ${$vcf{$chr02d}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads\n";
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	print ${$vcf{$chr}}{$pos};
    }
}
