#!/usr/bin/perl -w
use strict;

my $var_file = shift @ARGV;

my %vcf;

open (FILE, $var_file) or die "$var_file is not found: $!\n" if ($var_file !~ /\.gz$/);
open (FILE, "gzip -dc $var_file |") or die "$var_file is not found: $!\n" if ($var_file =~ /\.gz$/);
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
	print $line, "\n";
	next;
    }	
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    my $len = 0;
    $len = $1 if ($line[7] =~ /SVLEN=-*(\d+);/);
    my $end = $pos + $len - 1;
    my $reads = 30;
    if ($line[7] =~ /BD_SUPPORTING_READ_PAIRS=(\d+)/){
	$reads = $1;
    }
    elsif($line[7] =~ /PD_UNIQ_READ_SUPP=(\d+)/){
	$reads = $1;
    }
    elsif($line[7] =~ /PD_READ_SUPP=(\d+)/){
	$reads = $1;
    }
    my $chr2 = '';
    my $pos2 = 0;
    if (($type eq 'CTX') or ($type eq 'CTX')){
	$type = 'TRA';
	$chr2 = $1 if ($line[7] =~ /CHR2=(.+?);/);
	$pos2 = $1 if ($line[7] =~ /POS2=(.+?);/);
	$len = 0;
	$end = 0;
    }
    elsif ($type eq 'CNV'){
	if ($line[7] =~ /SVLEN=-\d+;/){
	    $type = 'DEL';
	}
	else{
	    $type = 'DUP';
	}
    }
    my $qual = $line[6];
    if (($type eq 'DEL') or ($type eq 'INV')){
	next if ($qual ne 'PASS');				# for DEL and INV, low quality sites are removed
    }
    print "$chr\t$pos\t$type\t.\t.\t.\t$qual\tSVTYPE=$type;SVLEN=$len;READS=$reads\n" if ($type ne 'TRA');
    print "$chr\t$pos\t$type\t.\t.\t.\t$qual\tSVTYPE=$type;SVLEN=$len;READS=$reads\tCHR2=$chr2;POS2=$pos2\n" if ($type eq 'TRA');
}
close (FILE);
