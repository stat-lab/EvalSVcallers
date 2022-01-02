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
	next if ($line =~ /^#/);
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	my $pos = $line[1];
	my $end = $line[2];
	my $len = $end - $pos + 1;
	my $type = 'DUP';
	my $cn = $line[4];
	my $read = 1;
	if ($cn >= 1){
	    next if ($cn < 1.1);
	    if ($cn < 1.5){
		$read = 3;
	    }
	    elsif ($cn < 2.0){
		$read = 4;
	    }
	    elsif ($cn < 2.5){
		$read = 5;
	    }
	    elsif ($cn < 3.0){
		$read = 6;
	    }
	    elsif ($cn < 4.0){
		$read = 7;
	    }
	    elsif ($cn >= 4.0){
		$read = 8;
	    }
	}
	else{
	    $type = 'DEL';
	    next if ($cn > 0.9);
	    if ($cn > 0.5){
		$read = 3;
	    }
	    elsif ($cn > 0.2){
		$read = 4;
	    }
	    elsif ($cn > 0.1){
		$read = 5;
	    }
	    elsif ($cn > 0.05){
		$read = 6;
	    }
	    elsif ($cn > 0.01){
		$read = 7;
	    }
	    elsif ($cn <= 0.01){
		$read = 8;
	    }
	}
	my $gt = './.';
	$gt = '1/1' if ($cn < 0.5) or ($cn >= 3.5);
	$gt = '1/0' if (($cn >= 0.5) and ($cn <= 1.75)) or (($cn >= 2.5) and ($cn < 3.5));
	my $chr02d = $chr;
	$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	${${$vcf{$chr02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$read;GT=$gt";
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
	    print "${${$vcf{$chr}}{$pos}}{$type}\n";
	}
    }
}
