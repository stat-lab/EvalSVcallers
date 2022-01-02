#!/usr/bin/perl -w
use strict;
use File::Basename;

# covert PennCNV-Seq output file (*.rawcnv) to vcf

my $file = shift @ARGV;

my $min_probe = 3;

my %vcf;

open (FILE, $file);
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^WARNING/);
    my @line = split (/\s+/, $line);
    my $chr_pos = $line[0];
    my $snp = $1 if ($line[1] =~ /numsnp=([\d,]+)/);
    next if ($snp < $min_probe);
    my ($chr, $range) = split (/:/, $chr_pos);
    $chr =~ s/^chr// if ($chr =~ /^chr/);
    my $chr02d = $chr;
    $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
    my ($pos, $end) = split (/-/, $range);
    my $len = $end - $pos + 1;
    my $cn = $1 if ($line[3] =~ /cn=(\d+)/);
    next if ($cn == 2);
    my $type = 'DEL';
    $type = 'DUP' if ($cn > 2);
    my $gt = '0/1';
    $gt = '1/1' if ($type eq 'DEL') and ($cn == 0);
    $gt = '1/1' if ($type eq 'DUP') and ($cn >= 4);
    ${${$vcf{$chr02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t\.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=3;CN=$cn;GT=$gt";
}
close (FILE);

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
        foreach my $type (keys %{${vcf{$chr}}{$pos}}){
            print ${${vcf{$chr}}{$pos}}{$type}, "\n";
            
        }
    }
}
