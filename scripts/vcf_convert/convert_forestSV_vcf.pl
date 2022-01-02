#!/usr/bin/perl -w
use strict;

# score <= 0.2, reads 2
# score <= 0.3, reads 3
# score <= 0.8, reads 8
# score <= 0.9, reads 9

my %vcf;

foreach my $file (@ARGV){
    if (!-e $file){
	print STDERR "$file is not found\n";
	next;
    }
    open (FILE, $file);
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^chr\tstart/){
	    next;
	}	
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	my $pos = $line[1];
	my $end = $line[2];
	my $len = $line[3];
	my $type = $line[5];
	my $score = $line[4];
	my $reads = 0;
	if ($score > 0.9){
	    $reads = 10;
	}
	elsif ($reads > 0.8){
	    $reads = 9;
	}
	elsif ($score > 0.7){
	    $reads = 8;
	}
	elsif ($score > 0.6){
	    $reads = 7;
	}
	elsif ($score > 0.5){
	    $reads = 6;
	}
	elsif ($score > 0.4){
	    $reads = 5;
	}
	elsif ($score > 0.3){
	    $reads = 4;
	}
	elsif ($score > 0.2){
	    $reads = 3
	}
	else{
	    $reads = 2;
	}
	my $chr_02d = $chr;
	$chr_02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	${${$vcf{$chr_02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads";
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
