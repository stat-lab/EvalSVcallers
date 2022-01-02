#!/usr/bin/perl -w
use strict;

my %vcf;

foreach my $file (@ARGV){
    open (FILE, $file) or die "$file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        if (($line =~ /^#/) or ($line =~ /^sample\tchr/)){
        	next;
        }
        my @line = split (/\t/, $line);
        my $chr = $line[0];
        my $chr02d = $chr;
        $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
        my $pos = $line[1];
        my $ref = $line[3];
        my $type = $1 if ($line[7] =~ /TYPE=(.+?);/);
        if ($type eq 'ALU'){
        	$type = 'ALU';
        }
        elsif ($type eq 'SV'){
        	$type = 'SVA';
        }
        elsif ($type eq 'HE'){
        	$type = 'HERVK';
        }
        elsif ($type eq 'MT'){
        	$type = 'NUMT';
        }
        elsif (($type ne 'L1') and ($type ne 'LINE1')){
        	$type = 'VEI';
        }
        my $len = 0;
        $len = $1 if ($line[7] =~ /MEILEN=(\d+)/);
        my $reads = 2;
        my $read5 = $1 if ($line[7] =~ /SR5=(\d+);/);
        my $read3 = $1 if ($line[7] =~ /SR3=(\d+);/);
        $reads = $read5 if ($read5 >= $read3) and ($read5 > 2);
        $reads = $read3 if ($read5 < $read3) and ($read3 > 2);
        my $strand = $1 if ($line[7] =~ /STRAND=([+\-])/);
        my $class = $type;
        if (($type ne 'NUMT') and ($type ne 'VEI')){
        	$class = 'MEI';
        }
        my $gt = $1 if ($line[9] =~ /^(.+?):/);
        $gt = './.' if ($gt ne '1/1') and ($gt ne '0/1') and ($gt ne '1/0');
        ${$vcf{$chr02d}}{$pos} = "$chr\t$pos\t$type\t$ref\t.\t.\tPASS\tSVTYPE=$class;SVLEN=$len;READS=$reads;STRAND=$strand;GT=$gt\n";
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    my $pre_pos = 0;
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	if (($pre_pos > 0) and ($pos - $pre_pos + 1 < $min_sv)){
	    next;
	}
	print ${$vcf{$chr}}{$pos};
	$pre_pos = $pos;
    }
}
