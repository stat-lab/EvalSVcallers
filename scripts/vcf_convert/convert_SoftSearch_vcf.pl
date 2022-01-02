#!/usr/bin/perl -w
use strict;

my %vcf;
my %used;

foreach my $var_file (@ARGV){
    open (FILE, $var_file) or die "$var_file is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#/){
	    next;
	}	
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	$chr = $1 if ($chr =~ /chr(.+)$/);
	my $pos = $line[1];
	my $type = $1 if ($line[7] =~ /EVENT=(.+?);/);
	my @rn = split (/:/, $line[9]);
	my $reads = 0;
	$reads = $rn[1] if ($type eq 'CTX');
	$reads = $rn[2] if ($type eq 'DEL');
	$reads = $rn[3] if ($type eq 'INS');
	$reads = $rn[4] if ($type eq 'INV');
	$reads = $rn[5] if ($type eq 'NOV_INS');
	$reads = $rn[6] if ($type eq 'TDUP');
	$type = 'INS' if ($type eq 'NOV_INS');
	$type = 'DUP' if ($type eq 'TDUP');
	my $len = $1 if ($line[7] =~ /ISIZE=(\d+);/);
	$len = 0 if ($len < 20) and ($type eq 'INS');
	my $end = $pos + $len - 1;
	my $chr2 = '';
	my $pos2 = 0;
	my $subtype = '';
	if ($type eq 'CTX'){
	    $line[4] =~ /[\[\]](.+?):(\d+)[\[\]]/;
	    $chr2 = $1;
	    $pos2 = $2;
	    if ($chr eq $chr2){
		$type = 'INS';
		$subtype = 'CTX';
	    }
	    else{
		$type = 'TRA';
		$len = 0;
		$end = 0;
	    }
	}
	next if ($type eq 'TRA') and (exists ${${$vcf{$chr2}}{$pos2}}{'TRA'});
	my $chr_second = '';
	my $pos_second = 0;
	if ($line[4] =~ /([^\[\]]+):(\d+)/){
		$chr_second = $1;
		$pos_second = $2;
	}
#print STDERR "$chr\t$pos\t$chr_second:$pos_second\n";
	next if (exists ${${$used{$type}}{$chr}}{$pos});
	my $chr02d = $chr;
	$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	${${$vcf{$chr02d}}{$pos}}{$type} = "SVTYPE=$type;SVLEN=$len;READS=$reads" if ($type ne 'TRA');
	${${$vcf{$chr02d}}{$pos}}{$type} = "SVTYPE=$type;SVLEN=$len;READS=$reads\tCHR2=$chr2;POS2=$pos2" if ($type eq 'TRA');
	${${$used{$type}}{$chr_second}}{$pos_second} = 1;
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    my $chr2 = $chr;
    $chr2 =~ s/^0*//;
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
	    print "$chr2\t$pos\t$type\t.\t.\t.\tPASS\t${${$vcf{$chr}}{$pos}}{$type}\n";
	}
    }
}
