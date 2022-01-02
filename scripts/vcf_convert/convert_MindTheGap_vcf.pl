#!/usr/bin/perl -w
use strict;

my $var_file = shift @ARGV;

my $var_sd = 125;

my %ins;


my $pre_chr = 0;
my $pre_pos = 0;
my $pre_len = 0;
open (FILE, $var_file) or die "$var_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line !~ /^>/);
    my $type = 'INS';
    my $chr = '';
    my $pos = 0;
    my $len = 0;
    if ($line =~ /_contig_\d+_(.+?)_pos_(\d+)/){
        $chr = $1;
        $pos = $2;
        next if ($chr !~ /^[\dXY]+$/);
    }
    elsif ($line =~ /bkpt\d+_(.+?)_pos_(\d+)/){
        $chr = $1;
        $pos = $2;
        next if ($chr !~ /^[\dXY]+$/);
    }
    else{
        next;
    }
    $len = $1 if ($line =~ /\(\s+len=\s(\d+)\s+\)/);
    $len = $1 if ($line =~ /len_(\d+)/);
    if ($pre_chr eq $chr){
        if ($pos > $pre_pos + $var_sd){
            ${$ins{$chr}}{$pos} = $len;
        }
        else{
            if ($len >= $pre_len){
                if (exists ${$ins{$pre_chr}}{$pre_pos}){
                    delete ${$ins{$pre_chr}}{$pre_pos};
                }
                ${$ins{$chr}}{$pos} = $len;
            }
            else{
                next;
            }
        }
    }
    $pre_pos = $pos;
    $pre_len = $len;
    $pre_chr = $chr;
}
close (FILE);

foreach my $chr (keys %ins){
    foreach my $pos (sort {$a <=> $b} keys %{$ins{$chr}}){
        my $type = 'INS';
        my $len = ${$ins{$chr}}{$pos};
        print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=15\n";
    }
}
