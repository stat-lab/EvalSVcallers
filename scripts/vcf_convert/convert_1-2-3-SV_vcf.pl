#!/usr/bin/perl -w
use strict;

# covert 1-2-3-SV output files to vcf

my $target_chr = '';

my %vcf;
my $type = '';

foreach my $file (@ARGV){
    if (!-e $file){
	print STDERR "$file is not found\n";
	next;
    }
    open (FILE, $file);
    if ($file =~ /deletion/){
	$type = 'DEL';
    }
    elsif ($file =~ /insertion/){
	$type = 'INS';
    }
    elsif ($file =~ /inversion/){
	$type = 'INV';
    }
    elsif ($file =~ /translocation/){
	$type = 'TRA';
    }
    while (my $line = <FILE>){
	chomp $line;
	next if ($line =~ /^#/);
	my @line = split (/\s+/, $line);
	my $chr = $line[0];
	my $bp1_s = $line[1];
	my $bp1_e = $line[2];
	my $bp2_s = $line[4];
	my $bp2_e = $line[5];
	my $pos = int (($bp1_s + $bp1_e) / 2);
	my $end = int (($bp2_s + $bp2_e) / 2);
	my $len = $end - $pos + 1;
	my $chr2 = $line[3];
	if ($len < 0){
	    my $bp2_s2 = $bp2_s;
	    my $bp2_e2 = $bp2_e;
	    my $end2 = $end;
	    $bp2_s = $bp1_s;
	    $bp2_e = $bp1_e;
	    $bp1_s = $bp2_s2;
	    $bp1_e = $bp2_e2;
	    $end = $pos;
	    $pos = $end2;
	    $len = 0 if ($type eq 'INS');
	    $len = 0 - $len;
	}
	my $pos2 = 0;
	if ($type eq 'TRA'){
	    $pos2 = int (($bp2_s + $bp2_e) / 2);
	    $len = 0;
	    $end = 0;
	}
	my $reads = $line[8];
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
