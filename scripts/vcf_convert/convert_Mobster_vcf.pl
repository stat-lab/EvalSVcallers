#!/usr/bin/perl -w
use strict;

my $file = shift @ARGV;
die "The second arg (MEI|VEI|NUMT) is missing:\n" if (@ARGV == 0);
my $tag = shift @ARGV;

open (FILE, $file) or die "$file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if (($line =~ /^#/) or ($line =~ /^Chr\tMobile/)){
	next;
    }
    my @line = split (/\t/, $line);
    my $NA_count = 0;
    if ($line[8] eq 'NA'){
	$NA_count ++;
    }
    if ($line[9] eq 'NA'){
	$NA_count ++;
    }
    if ($line[10] eq 'NA'){
	$NA_count ++;
    }
    if ($line[11] eq 'NA'){
	$NA_count ++;
    }
    if (($tag ne 'NUMT') and ($line[-1] eq 'unknown') and ($line[-2] == -1) and ($line[-3] == -1) and ($line[-4] == -1) and ($NA_count >= 2)){
	next;
    }
    my $chr = $line[0];
    $chr =~ s/^chr// if ($chr =~ /^chr/);
    my $pos = $line[2];
    my $type = $line[1];
    if (($type eq 'HERV') and (($tag eq 'NUMT') or ($tag eq 'MT'))){
	$type = 'NUMT';
    }
    elsif (($type eq 'HERV') and ($tag eq 'VEI')){
	$type = 'VEI';
    }
    my $len = 0;
    $len = $line[8] if ($line[8] ne 'NA') and ($line[9] eq 'NA');
    $len = $line[9] if ($line[9] ne 'NA') and ($line[8] eq 'NA');
    $len = int (($line[8] + $line[9]) / 2) if ($line[9] ne 'NA') and ($line[8] ne 'NA');
    next if ($len < 150);
    my $reads = 0;
    $reads = $1 if ($line[7] =~ /.+=(\d+)$/);
    my $split_reads = $line[12] + $line[13];
    next if ($split_reads < 0);
    my $class = $type;
    if (($type ne 'NUMT') and ($type ne 'VEI')){
	$class = 'MEI';
    }
    print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$class;SVLEN=$len;READS=$reads;SP_READS=$split_reads\n";
}
close (FILE);
