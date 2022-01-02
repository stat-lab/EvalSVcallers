#!/usr/bin/perl -w
use strict;

my %vcf;

my @calls;

if (@ARGV > 0){
    @calls = (@ARGV);
}
else{
    @calls = <./*.calls>;
}

foreach my $file (@calls){
    if (!-e $file){
	print STDERR "$file is not found\n";
	next;
    }
    open (FILE, $file);
    while (my $line = <FILE>){
	chomp $line;
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	my $chr2 = $line[3];
	next if ($chr ne $chr2);
	my $pos1 = $line[1];
	my $pos2 = $line[2];
	my $pos3 = $line[4];
	my $pos4 = $line[5];
	my $pos = int (($pos1 + $pos2) / 2 + 0.5);
	my $end = int (($pos3 + $pos4) / 2 + 0.5);
	my $len = $end - $pos + 1;
	my $type = $1 if ($line[6] =~ /(\w+?)\./);
	my $reads = 3;
	my $chr02d = $chr;
	$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	${$vcf{$chr02d}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads";
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	print "${$vcf{$chr}}{$pos}\n";
    }
}
