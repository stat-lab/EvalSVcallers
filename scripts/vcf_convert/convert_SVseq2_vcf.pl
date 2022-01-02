#!/usr/bin/perl -w
use strict;

my %vcf;
my $flag = 0;

foreach my $file (@ARGV){
    if (!-e $file){
	print STDERR "$file is not found\n";
	next;
    }
    open (FILE, $file);
    my $type = '';
    my $pre_line = '';
    my $reads = 0;
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#+$/){
	    if ($pre_line ne ''){
		my @line = split (/\s+/, $pre_line);
		my $len = 0;
		my $pos = 0;
		my $chr = '';
		if ($type eq 'DEL'){
		    $chr = $line[1];
		    $pos = $line[2];
		    $len = $line[4] - $pos + 1;
		}
		elsif ($type eq 'INS'){
		    $chr = $line[0];
		    $pos = $line[1];
		}
		my $chr02d = $chr;
		$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
		${${$vcf{$chr02d}}{$pos}}{$type} = "$len=$reads";
	    }
	    $pre_line = '';
	    $type = '';
	    $flag = 0;
	    $reads = 0;
	}
	elsif ($flag == 0){
	    $pre_line = $line;
	    $flag = 1;
	    if ($line =~ /^range/){
		$type = 'DEL';
	    }
	    else{
		$type = 'INS';
	    }
	}
	elsif ($flag == 1){
	    if ($type eq 'DEL'){
		$reads ++;
	    }
	    elsif ($type eq 'INS'){
		$reads ++ if ($line =~ /^\d+/);
	    }
	}
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
	    my ($len, $reads) = split (/=/, ${${$vcf{$chr}}{$pos}}{$type});
	    my $chr2 = $chr;
	    $chr2 =~ s/^0*//;
	    print "$chr2\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads\n";
	}
    }
}
