#!/usr/bin/perl -w
use strict;

my $var_file = shift @ARGV;

my %vcf;
my $count = 0;

open (FILE, $var_file) or die "$var_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^chr_type/){
	next;
    }	
    my @line = split (/\t/, $line);
    my $chr = $line[3];
    $chr = $1 if ($chr =~ /chr(.+)/);
    my ($bp1_s, $bp1_e) = split (/-/, $line[4]);
    my ($bp2_s, $bp2_e) = split (/-/, $line[7]);
    my $pos = int (($bp1_s + $bp1_e) / 2);
    my $end = int (($bp2_s + $bp2_e) / 2);
    my $len = $line[5];
    my $type = $line[1];
    my $chr2 = '';
    my $pos2 = 0;
    next if ($type eq 'UNDEFINED');
    if ($line[0] eq 'INTRA'){
	if ($type eq 'DELETION'){
	    $type = 'DEL';
	}
	elsif (($type eq 'INSERTION') or ($type eq 'INS_FRAGMT')){
	    $type = 'INS';
	}
	elsif (($type eq 'DUPLICATION') or ($type eq 'SMALL_DUPLI') or ($type eq 'LARGE_DUPLI')){
	    $type = 'DUP';
	}
	elsif ($type =~ /^INV/){
	    $type = 'INV';
	}
	else{
	    next;
	}
    }
    elsif (($line[0] eq 'INTER') or ($type eq 'TRANSLOC') or ($type eq 'INV_TRANSLOC')){
	$type = 'TRA';
	$chr2 = $line[6];
	$chr2 = $1 if ($chr2 =~ /chr(.+)/);
	next if ($chr eq $chr2);
	$pos2 = $end;
	$end = 0;
	$len = 0;
    }
    else{
	next;
    }
#    next if ($len < $min_sv_len) and ($len > 0);
    my $score = $line[12];
    my $reads = 0;
    if ($score == 1){
	$reads = 15;
    }
    elsif ($score >= 0.98){
	$reads = 10;
    }
    elsif ($score >= 0.95){
	$reads = 7;
    }
    elsif ($score >= 0.9){
	$reads = 5;
    }
    elsif ($score >= 0.85){
	$reads = 3;
    }
    my $chr_02d = $chr;
    $chr_02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
    my $chr2_02d = $chr2;
    $chr2_02d = sprintf ("%02d", $chr2) if ($chr2 =~ /^\d+$/);
    next if ($type eq 'TRA') and (exists ${${$vcf{$chr2_02d}}{$pos2}}{$type});
    ${${$vcf{$chr_02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;END=$end;READS=$reads" if ($type ne 'TRA');
    ${${$vcf{$chr_02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;END=$end;READS=$reads\tCHR2=$chr2;POS2=$pos2" if ($type eq 'TRA');
}
close (FILE);

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
	    print ${${$vcf{$chr}}{$pos}}{$type}, "\n";
	}
    }
}
