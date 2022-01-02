#!/usr/bin/perl -w
use strict;

my $var_file = shift @ARGV;

my %vcf;

open (FILE, $var_file) or die "$var_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
	next;
    }	
    my @line = split (/\t/, $line);
    my ($chr, $pos) = split (/:/, $line[0]);
    next if (!defined $pos) or ($pos eq '');
    my $end = $1 if ($line[12] =~ /$chr:(\d+)/);
    next if (!defined $end) or ($end eq '');
    if ($pos > $end){
	my $pos2 = $pos;
	$pos = $end;
	$end = $pos2;
    }
    my $len = $end - $pos + 1;
    my $type = '';
    if ($len >= $min_sv_len){
	$type = 'DEL';
    }
    else{
	$type = 'INS';
    }
    $len = 0 if ($type eq 'INS');
    my $reads = $line[6];
    my $chr_02d = $chr;
    $chr_02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
    ${$vcf{$chr}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads";
}
close (FILE);

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	print ${$vcf{$chr}}{$pos}, "\n";
    }
}
