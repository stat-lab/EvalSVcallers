#!/usr/bin/perl -w
use strict;

my $vcf_file = shift @ARGV;

my %vcf;
my %vcf2;

open (FILE, $vcf_file) or die "$vcf_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#/);
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $qual = $line[5];
    my $len = $1 if ($line[7] =~ /SPAN=(-*\d+)/);
    my $alt = $line[4];
    my $gt = './.';
    $gt = $1 if ($line[9] =~ /^([01]\/[01])/);
    $gt = './.' if ($gt eq '0/0');
    my $mate_chr = $1 if ($alt =~ /([XY\d]+):\d+/);
    my $mate_pos = $1 if ($alt =~ /[XY\d]+:(\d+)/);
    if (exists ${$vcf{$mate_chr}}{$mate_pos}){
        my ($qual1, $len1, $gt1, $alt1) = split (/=/, ${$vcf{$mate_chr}}{$mate_pos});
        my $dir1 = '';
        my $dir2 = '';
        $dir1 = 'D' if ($alt1 =~ /^[ACGT]\[/);
        $dir1 = 'U' if ($alt1 =~ /.+\][ACGT]/);
        $dir1 = 'V' if ($alt1 =~ /.+\[[ACGT]/);
        $dir2 = 'D' if ($alt =~ /.+\][ACGT]/);
        $dir2 = 'U' if ($alt =~ /^[ACGT]\[/);
        $dir2 = 'V' if ($alt =~ /.+\[[ACGT]/);
        my $type = '';
        if (($dir1 eq 'D') and ($dir2 eq 'D')){
            $type = 'DEL';
        }
        elsif (($dir1 eq 'U') and ($dir2 eq 'U')){
            $type = 'DUP';
        }
        elsif (($dir1 eq 'V') and ($dir2 eq 'V')){
            $type = 'INV';
        }
        if ($type ne ''){
            my $chr02d = $mate_chr;
            $chr02d = sprintf ("%02d", $mate_chr) if ($mate_chr =~ /^\d+$/);
            ${$vcf2{$chr02d}}{$mate_pos} = "${$vcf{$mate_chr}}{$mate_pos}=$type";
        }
    }
    else{
        if ($len == -1){
            my $chr02d = $chr;
            $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
            ${$vcf2{$chr02d}}{$pos} = "$qual=0=$gt=$alt=INS";
        }
        else{
            ${$vcf{$chr}}{$pos} = "$qual=$len=$gt=$alt";
        }
    }
}
close (FILE);

foreach my $chr (sort keys %vcf2){
    my $chr2 = $chr;
    $chr2 =~ s/^0*//;
    foreach my $pos (sort {$a <=> $b} keys %{$vcf2{$chr}}){
        my ($qual, $len, $gt, $alt, $type) = split (/=/, ${$vcf2{$chr}}{$pos});
        my $read = 2;
        $read = int ($qual / 10 + 0.5) if ($qual > 20);
        print "$chr2\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$read;GT=$gt\n";
    }
}
