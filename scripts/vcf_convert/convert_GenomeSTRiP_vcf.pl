#!/usr/bin/perl -w
use strict;

my @var_file = (@ARGV);

my %vcf;

my $target_sample = 'NA12878';

foreach my $var_file (@var_file){
	open (FILE, $var_file) or die "$var_file is not found: $!\n" if ($var_file !~ /\.gz$/);
	open (FILE, "gzip -dc $var_file |") or die "$var_file is not found: $!\n" if ($var_file =~ /\.gz$/);
	while (my $line = <FILE>){
		chomp $line;
		if ($line =~ /^#/){
			next;
		}	
		my @line = split (/\t/, $line);
		my $chr = $line[0];
		my $chr02d = $chr;
		$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
		my $pos = $line[1];
		my $end = 0;
		my $type = $1 if ($line[7] =~ /SVTYPE=([^;]+)/);
		$end = $1 if ($line[7] =~ /END=(\d+);/);
		my $len = 0;
		$len = $1 if ($line[7] =~ /SVLEN=-*(\d+);/);
		my $reads = $1 if ($line[7] =~ /GSNPAIRS=(\d+);/);
		my $gt = './.';
		$gt = $1 if (@line > 9) and ($line[9] =~ /^(.+?):/);
		if ($line =~ /GSSAMPLES=([^;]+)/){
			my $sample_str = $1;
			if ($sample_str !~ /$target_sample/){
				next;
			}
		}
		${$vcf{$chr02d}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads;GT=$gt\n";
	}
	close (FILE);
}

foreach my $chr (sort keys %vcf){
	foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
		print ${$vcf{$chr}}{$pos};
	}
}
