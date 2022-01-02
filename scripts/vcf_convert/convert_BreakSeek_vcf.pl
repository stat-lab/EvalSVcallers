#!/usr/bin/perl -w
use strict;

my %vcf;

foreach my $file (@ARGV){
    open (FILE, $file) or die "$file is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^ST/){
	    print $line, "\n";
	    next;
	}	
	my @line = split (/\t/, $line);
	my $chr = $line[3];
	my $pos = $line[0];
	my $end = $line[1];
	my $type = $line[2];
	$type = 'DEL' if ($type eq 'D');
	$type = 'INS' if ($type eq 'I');
	my $len = $end - $pos + 1;
	my $reads = 0;
	my @readL = ();
	my @readR = ();
	@readL = split (/\//, $line[5]) if ($line[5] ne '-');
	@readR = split (/\//, $line[6]) if ($line[6] ne '-');
	sort {$b <=> $a} @readL;
	sort {$b <=> $a} @readR;
	if ((@readL > 0) and (@readR > 0)){
	    $reads = int (($readL[0] + $readR[0]) / 2 + 0.5);
	}
	elsif (@readL > 0){
	    $reads = $readL[0];
	}
	elsif (@readR > 0){
	    $reads = $readR[0];
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
