#!/usr/bin/perl -w
use strict;

my %var;

foreach my $var_file (@ARGV){
    open (FILE, $var_file) or die "$var_file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        if ($line =~ /^type/){
            next;
        }
        my @line = split (/\t/, $line);
        my $chr = $line[2];
        $chr = $1 if ($chr =~ /^chr(.+)/);
        my $type = $line[0];
        my $bp1_s = $line[3];
        my $bp1_e = $line[4];
        my $bp2_s = $line[7];
        my $bp2_e = $line[8];
        my $bp1 = int (($bp1_s + $bp1_e) / 2);
        my $bp2 = int (($bp2_s + $bp2_e) / 2);
        my $pos = $bp1;
        my $end = $bp2;
        my $chr2 = '';
        my $len = 0;
        if ($type eq 'Inversion'){
            $type = 'INV';
        }
        elsif ($type eq 'Deletion'){
            $type = 'DEL';
        }
        elsif ($type eq 'Insertion'){
            $type = 'INS';
        }
        elsif ($type eq 'Intrachromosomal_translocation'){
            $chr2 = $line[7];
            if ($chr eq $chr2){
                $type = 'INS-T';
            }
            else{
                $type = 'TRA';
            }
        }
            elsif ($type eq 'Interchromosomal_translocation'){
                $chr2 = $line[7];
                $type = 'TRA';
            }
            if (($bp1 > $bp2) and ($type ne 'INS-T')){
                $pos = $bp2;
                $end = $bp1;
                $bp1_s = $line[7];
                $bp1_e = $line[8];
                $bp2_s = $line[3];
                $bp2_e = $line[4];
            }

            $len = $end - $pos + 1 if ($type ne 'INS-T');
            $type = 'INS' if ($type eq 'INS-T');
            my $pos2 = 0;
            if ($type eq 'TRA'){
                $end = 0;
                $len = 0;
            }
        my $reads = $line[1];
        ${${$var{$chr}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads" if ($type ne 'TRA');
        ${${$var{$chr}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads\tCHR2=$chr2;POS2=$pos2" if ($type eq 'TRA');
    }
    close (FILE);
}

foreach my $chr (sort keys %var){
    foreach my $pos1 (sort {$a <=> $b} keys %{$var{$chr}}){
        foreach my $type (keys %{${$var{$chr}}{$pos1}}){
            print ${${$var{$chr}}{$pos1}}{$type}, "\n";
        }
    }
}
