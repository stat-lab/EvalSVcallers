#!/usr/bin/perl -w
use strict;

my %vcf;

foreach my $var_file (@ARGV){
    open (FILE, $var_file) or die "$var_file is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	if (($line =~ /^#/) or ($line =~ /^Chromosome/)){
	    next;
	}	
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	$chr = $1 if ($chr =~ /chr(.+)$/);
	my $pos = $line[1];
	my $len = $line[3];
	my $type = '';
	if ($var_file =~ /deletion/){
	    $type = 'DEL';
	}
	elsif ($var_file =~ /insertion/){
	    $type = 'INS';
	}
	elsif ($var_file =~ /inversion/){
	    $type = 'INV';
	}
	elsif ($var_file =~ /tandem/){
	    $type = 'DUP';
	}
	elsif ($var_file =~ /translocation/){
	    $type = 'TRA';
	}
	my $reads = $line[5];
	my $end = $pos + $len - 1;
	my $chr2 = '';
	my $pos2 = 0;
	my $subtype = '';
	if ($type eq 'TRA'){
	    $len = 0;
	    $chr2 = $line[2];
	    $chr2 = $1 if ($chr2 =~ /chr(.+)$/);
	    $pos2 = $line[3];
	}
	my $chr02d = $chr;
	$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	${${$vcf{$chr02d}}{$pos}}{$type} = "SVTYPE=$type;SVLEN=$len;READS=$reads" if ($type ne 'TRA');
	${${$vcf{$chr02d}}{$pos}}{$type} = "SVTYPE=$type;SVLEN=$len;READS=$reads\tCHR2=$chr2;POS2=$pos2" if ($type eq 'TRA');
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
