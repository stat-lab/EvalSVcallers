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
	next if ($line =~ /^#|^$/);
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	my $pos = $line[1];
	my $type = $line[3];
	my $len = $line[4];
	my $reads = $1 if ($line[5] =~ /szCount=(\d+)/);
	my $gt = 'GT=./.';
	$gt = $1 if ($line[5] =~ /(GT=[01\.]\/[01\.])/);
	my $chr02d =$chr;
	$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	${$vcf{$chr02d}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads;$gt";
    }
    close (FILE);
}

foreach my $chr (sort %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	print ${$vcf{$chr}}{$pos}, "\n";
    }
}
