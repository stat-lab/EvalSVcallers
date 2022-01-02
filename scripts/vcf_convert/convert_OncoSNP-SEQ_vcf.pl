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
	next if ($line =~ /^#|^Chromosome/);
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	my $pos = $line[1];
	my $end = $line[2];
	my $len = $end - $pos + 1;
	my $copy_num = $line[3];
	my $type = 'DUP' if ($copy_num >= 2);
	$type = 'DEL' if ($copy_num < 2);
	my $log_like = $line[6];
	my $reads = 1;
	if ($log_like < -10000){
	    $reads = 3;
	}
	elsif ($log_like < -1000){
	    $reads = 5;
	}
	elsif ($log_like < -400){
	    $reads = 7;
	}
	elsif ($log_like < -100){
	    $reads = 9;
	}
	else{
	    $reads = 11;
	}
	my $chr02d = $chr;
	$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	${$vcf{$chr02d}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;COPYNUM=$copy_num;READS=$reads\n";
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	print ${$vcf{$chr}}{$pos};
    }
}
