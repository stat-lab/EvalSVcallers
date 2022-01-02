#!/usr/bin/perl -w
use strict;

# covert Delly output files to vcf

my $min_gain = 1.1;
my $max_loss = 0.9;

my %vcf;
my $count = 0;

my @files = (@ARGV) if (@ARGV > 0);
@files = <./*_segment.dat> if (@ARGV == 0);

print STDERR "@files\n";

foreach my $var_file (@files){
    open (FILE, $var_file) or die "$var_file is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#/){
	    next;
	}
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	next if ($chr eq 'chr');
	$chr = $1 if ($chr =~ /^chr(\.+)/);
	my $pos = $line[1];
	my $end = $line[2];
	my $len = $end - $pos + 1;
	my $type = 'DUP' if ($line[4] >= 2);
	$type = 'DEL' if ($line[4] < 2);
	my $score = $line[6];
	my $reads = 0;
	if ($score <= 1){
	    $reads = 3;
	}
	elsif ($score <= 5){
	    $reads = 4;
	}
	elsif ($score <= 10){
	    $reads = 5;
	}
	elsif ($score <= 20){
	    $reads = 6;
	}
	elsif ($score <= 30){
	    $reads = 7;
	}
	else{
	    $reads = 8;
	}
	my $cn = $line[4];
	my $gt = './.';
	$gt = '1/0' if ($cn == 1) or ($cn == 3);
	$gt = '1/1' if ($cn == 0) or ($cn >= 4);
	my $chr_02d = $chr;
	$chr_02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	${${$vcf{$chr_02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads;GT=$gt";
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
	    print ${${$vcf{$chr}}{$pos}}{$type}, "\n";
	}
    }
}
