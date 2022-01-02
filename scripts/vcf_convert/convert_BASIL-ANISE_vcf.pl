#!/usr/bin/perl -w
use strict;

# covert Basil-Anise output files to vcf

my $var_file = shift @ARGV;

my $ins_file = '';
$ins_file = shift @ARGV if (@ARGV > 0) and (-f $ARGV[0]);

my %ins_seq;

my $seq = '';
my $header = '';
if ($ins_file ne ''){
    open (FILE, $ins_file) or die "$ins_file is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^>/){
	    if ($seq ne ''){
		my $chr = $1 if ($header =~ /REF=(\S+)/);
		my $pos = $1 if ($header =~ /POS=(\d+)/);
		${$ins_seq{$chr}}{$pos} = $seq;
	    }
	    $seq = '';
	    $header = $1 if ($line =~ /^>(.+)/);
	}
	else{
	    $seq .= uc $line;
	}
    }
    if ($seq ne ''){
	my $chr = $1 if ($header =~ /REF=(\S+)/);
	my $pos = $1 if ($header =~ /POS=(\d+)/);
	${$ins_seq{$chr}}{$pos} = $seq;
    }
    close (FILE);
}

open (FILE, $var_file) or die "$var_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
	print $line, "\n";
	next;
    }	
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $type = $1 if ($line[4] =~ /<(.+)>/);
    my $len = 0;
    if (exists ${$ins_seq{$chr}}{$pos}){
	$len = length ${$ins_seq{$chr}}{$pos};
    }
    my $end = $pos + $len - 1;
    my $oea_read_1 = 0;
    my $oea_read_2 = 0;
    $oea_read_1 = $1 if ($line[9] =~ /:(\d+):\d+$/);
    $oea_read_2 = $1 if ($line[9] =~ /:\d+:(\d+)$/);
    my $reads = 0;
    $reads = int (($oea_read_1 + $oea_read_2) / 2 + 0.5);
    print"$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads\n";
}
close (FILE);
