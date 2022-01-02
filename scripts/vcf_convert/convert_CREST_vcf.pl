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
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	my $type = $line[8];
	next if ($type eq 'ITX');
	my $pos = $line[1];
	my $end = $line[5];
	my $len = $end - $pos + 1;
	my $reads = int (($line[3] + $line[7]) / 2 + 0.5);
	my $chr2 = '';
	my $pos2 = 0;
	my $chr_02d = $chr;
	$chr_02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	if ($type eq 'CTX'){
	    $type = 'TRA';
	    $chr2 = $line[4];
	    $pos2 = $line[5];
	    $end = 0;
	    $len = 0;
	    my $chr2_02d = $chr2;
	    $chr2_02d = sprintf ("%02d", $chr2) if ($chr2 =~ /^\d+$/);
	    next if (exists ${${$vcf{$chr2_02d}}{$pos2}}{'TRA'});
	    ${${$vcf{$chr_02d}}{$pos}}{$type} = "SVTYPE=$type;SVLEN=$len;READS=$reads\tCHR2=$chr2;POS2=$pos2";
	}
	else{
	    ${${$vcf{$chr_02d}}{$pos}}{$type} = "SVTYPE=$type;SVLEN=$len;READS=$reads";
	}
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
	    print "$chr\t$pos\t$type\t.\t.\t.\tPASS\t${${$vcf{$chr}}{$pos}}{$type}\n";
	}
    }
}
