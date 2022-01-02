#!/usr/bin/perl -w
use strict;

my %vcf;

foreach my $file (@ARGV){
    open (FILE, $file) or die "$file is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#/){
	    next;
	}
	if ($line =~ /^Chr:/){
	    my @line = split (/\s+/, $line);
	    my $chr = $line[0];
	    $chr =~ s/^Chr:// if ($chr =~ /^Chr:/);
	    my $type = $line[5];
	    if ($type eq 'SVtype:D'){
		$type = 'DEL';
	    }
	    elsif ($type eq 'SVtype:V'){
		$type = 'INV';
	    }
	    elsif ($type eq 'SVtype:I'){
		$type = 'INS';
	    }
	    my $pos = 0;
	    my $len = 0;
	    my $start_1 = $1 if ($line[1] =~ /Start_Outer:(\d+)/);
	    my $start_2 = $1 if ($line[2] =~ /Start_Inner:(\d+)/);
	    my $end_1 = $1 if ($line[3] =~ /End_Inner:(\d+)/);
	    my $end_2 = $1 if ($line[4] =~ /End_Outer:(\d+)/);
	    if ($type eq 'INS'){
		$pos = int (($start_1 + $start_2 + $end_1 + $end_2) / 4);
	    }
	    else{
		$pos = int (($start_1 + $start_2) / 2);
		my $end = int (($end_1 + $end_2) / 2);
		$len = $end - $pos + 1;
	    }
	    my $reads = 0;
	    $reads = $1 if ($line[6] =~ /sup:(\d+)/);
	    ${$vcf{$chr}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads\n";
	}
    }
    close (FILE);
}

foreach my $chr (keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	print ${$vcf{$chr}}{$pos};
    }
}
