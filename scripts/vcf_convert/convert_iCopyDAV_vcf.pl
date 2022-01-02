#!/usr/bin/perl -w
use strict;

my %vcf;

foreach my $var_file (@ARGV){
    open (FILE, $var_file) or die "$var_file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        my ($chr, $start, $end, $tag, $cn) = split (/\t/, $line);
        my $chr02d = $chr;
        $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
        my $type = 'DEL';
        $type = 'DUP' if ($tag == 1);
        my $len = $end - $start + 1;
        $cn = int ($cn * 100 + 0.5) / 100;
        ${$vcf{$chr02d}}{$start} = "$chr\t$start\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;CN=$cn;READS=3\n";
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
        print ${$vcf{$chr}}{$pos};
    }
}
