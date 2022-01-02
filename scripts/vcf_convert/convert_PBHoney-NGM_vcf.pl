#!/usr/bin/perl -w
use strict;

my $min_supp_read_rate = 0.2;

my %vcf;

my $spot_del_file = shift @ARGV;

my $spot_ins_file = shift @ARGV;

open (FILE, $spot_del_file) or die "$spot_del_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#|^$/);
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $type = $line[3];
    my $len = $line[4];
    next if ($type ne 'DEL');
    my $reads = $1 if ($line[5] =~ /szCount=(\d+)/);
    my $gt = 'GT=./.';
    $gt = $1 if ($line[5] =~ /(GT=[01\.]\/[01\.])/);
    my $coverage = $1 if ($line[5] =~ /coverage=([\d\.]+)/);
    $coverage = int ($coverage + 0.5);
    my $supp_rate = $reads / $coverage;
    next if ($supp_rate < $min_supp_read_rate);
    my $chr02d =$chr;
    $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
    ${$vcf{$chr02d}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads;$gt";
}
close (FILE);

open (FILE, $spot_ins_file) or die "$spot_ins_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#|^$/);
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $type = $line[3];
    my $len = $line[4];
    next if ($type ne 'INS');
    my $reads = $1 if ($line[5] =~ /szCount=(\d+)/);
    my $gt = 'GT=./.';
    $gt = $1 if ($line[5] =~ /(GT=[01\.]\/[01\.])/);
    my $coverage = $1 if ($line[5] =~ /coverage=([\d\.]+)/);
    $coverage = int ($coverage + 0.5);
    my $supp_rate = $reads / $coverage;
    next if ($supp_rate < $min_supp_read_rate);
    my $chr02d =$chr;
    $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
    ${$vcf{$chr02d}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads;$gt";
}
close (FILE);

foreach my $chr (sort %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	print ${$vcf{$chr}}{$pos}, "\n";
    }
}
