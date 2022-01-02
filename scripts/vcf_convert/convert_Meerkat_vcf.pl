#!/usr/bin/perl -w
use strict;

my %vcf;

foreach my $var_file (@ARGV){
    if (!-e $var_file){
	print STDERR "$var_file is not found\n";
	next;
    }
    open (FILE, $var_file);
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#/){
    #	print $line, "\n";
	    next;
	}	
	my @line = split (/\t/, $line);
	my $type = $line[0];
	$type = 'DEL' if ($type =~ /^del/);
	$type = 'INS' if ($type =~ /^ins/);
	$type = 'DUP' if ($type =~ /^tandem/);
	$type = 'INV' if ($type =~ /^inv/);
	$type = 'TRA' if ($type =~ /^transl/);
	my $chr = $line[5];
	my $reads = $line[3];
	$reads = int (($1 + $2) / 2 + 0.5) if ($reads =~ /(\d+)\/(\d+)/);
	my $pos = 0;
	my $len = 0;
	my $chr2 = '';
	my $pos2 = 0;
	if ($type eq 'INS'){
	    $chr = $line[9];
	    $pos = $line[10];
	    $len = $line[12];
	    if ($len >= 10){
		$pos = int (($line[10] + $line[11]) / 2);
	    }
	    $len = 0 if ($len < 0);
	}
	elsif ($type eq 'TRA'){
	    $pos = $line[6];
	    $chr2 = $line[8];
	    $pos2 = $line[9];
	    $len = 0;
	}
	else{
	    $pos = $line[6];
	    $len = $line[8];
	}
	my $chr02d = $chr;
	$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	if ($pos > 0){
	    ${${$vcf{$chr02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads" if ($type ne 'TRA');
	    ${${$vcf{$chr02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads\tCHR2=$chr2;POS2=$pos2" if ($type eq 'TRA');
	}
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
