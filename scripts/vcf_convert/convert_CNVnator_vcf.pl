#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);

my $min_del_len = 3000;
my $max_del_len = 15000000;
my $max_dup_len = 30000000;

my $min_gain = 1.1;
my $max_loss = 0.9;

my $var_file = shift @ARGV;

my $long_read_flag = 0;

$long_read_flag = shift @ARGV if (@ARGV > 0);

my %vcf;
my $count = 0;

my $gap_bed = "$Bin/../../Ref_SV/gap.bed";

my %gap;

open (FILE, $gap_bed) or die "$gap_bed is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    my @line = split (/\s+/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $end = $line[2];
    ${$gap{$chr}}{$pos} = $end;
}
close (FILE);

open (FILE, $var_file) or die "$var_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
	next;
    }
    my @line = split (/\t/, $line);
    my $chr = '';
    my $pos = 0;
    my $end = 0;
    if ($line[1] =~ /^(.+?):(\d+)-(\d+)$/){
	$chr = $1;
	$pos = $2;
	$end = $3;
    }
    my $len = $line[2];
    if ($len =~ /^([\d\.]+)e\+(\d+)$/){
	my $const = $1;
	my $index = $2;
	$index =~ s/^0+//;
	$len = $const * (10 ** $index);
    }
    my $type = 'DUP' if ($line[0] eq 'duplication');
    $type = 'DEL' if ($line[0] eq 'deletion');
    next if ($type eq 'DEL') and ($len < $min_del_len);
    next if ($len > $max_del_len) and ($type eq 'DEL');
    next if ($len > $max_dup_len) and ($type eq 'DUP');
    if ($type eq 'DUP'){
	next if ($line[3] < $min_gain);
    }
    elsif ($type eq 'DEL'){
	next if ($line[3] > $max_loss);
    }
    my $chr2 = $chr;
    $chr = $1 if ($chr =~ /^chr(.+)/);
    my $gap_overlap = 0;
    foreach my $gstart (sort {$a <=> $b} keys %{$gap{$chr}}){
	my $gend = ${$gap{$chr}}{$gstart};
	if (($pos <= $gstart) and ($end >= $gend)){
	    if ($gend - $gstart + 1 >= $len * 0.5){
		$gap_overlap = 1;
		last;
	    }
	}
	elsif (($pos >= $gstart) and ($pos <= $gend)){
	    if ($gend - $pos + 1 >= $len * 0.5){
		$gap_overlap = 1;
		last;
	    }
	}
	elsif (($gstart >= $pos) and ($gstart <= $end)){
	    if ($end - $gstart + 1 >= $len * 0.5){
		$gap_overlap = 1;
		last;
	    }
	}
	last if ($gstart > $end);
    }
    next if ($gap_overlap == 1);
    my @evalue;
    push @evalue, $line[4], $line[5], $line[6], $line[7];
    @evalue = sort {$a <=> $b} @evalue;
    next if ($evalue[0] != 0) and ($evalue[0] > 1e-7);
    my $Q0_rate = $line[8];
    my $reads = 0;
    if ($type eq 'DEL'){
	if ($Q0_rate <= 0.025){
	    $reads = 12;
	    if ($len >= 500000){
		$reads = 14;
	    }
	}
	elsif ($Q0_rate <= 0.5){
	    $reads = 10;
	}
	elsif ($Q0_rate <= 0.8){
	    $reads = 7;
	}
	elsif ($Q0_rate <= 0.85){
	    $reads = 5;
	}
	elsif ($Q0_rate <= 0.9){
	    $reads = 4;
	}
	else{
	    $reads = 3;
	}
    }
    elsif ($type eq 'DUP'){
	next if ($Q0_rate > 0.7) and ($long_read_flag == 0);
        next if ($Q0_rate > 0.9) and ($long_read_flag == 1);
	if ($line[3] >= 2){
	    $reads = 10
	}
	elsif ($line[3] >= 1.8){
	    $reads = 7;
	}
	elsif ($line[3] >= 1.6){
	    $reads = 5;
	}
	elsif ($line[3] >= 1.4){
	    $reads = 4;
	}
	else{
	    $reads = 3;
	}
	$reads = 3 if ($len > 1000000) and ($len <= 2000000) and ($Q0_rate <= 0.01);
	$reads = 2 if ($len > 2000000) and ($Q0_rate <= 0.01)
    }
    my $chr_02d = $chr2;
    $chr_02d = sprintf ("%02d", $chr2) if ($chr2 =~ /^\d+$/);
    my $gt = './.';
#    $gt = '1/1' if ($line[3] < 0.25) or (($line[3] > 1.75) and ($line[3] <= 2.5));
#    $gt = '1/0' if (($line[3] >= 0.25) and ($line[3] < 1)) or (($line[3] > 1) and ($line[3] <= 1.75));
    $gt = '1/1' if ($line[3] < 0.25) or (($line[3] > 1.85) and ($line[3] <= 2.2));
    $gt = '1/0' if (($line[3] >= 0.25) and ($line[3] < 1)) or (($line[3] > 1.35) and ($line[3] <= 1.65));
    ${${$vcf{$chr_02d}}{$pos}}{$type} = "$chr2\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;RDRATIO=$line[3];READS=$reads;GT=$gt";
}
close (FILE);

foreach my $chr (sort keys %vcf){
    my $pre_pos = 0;
    my $pre_end = 0;
    my $pre_len = 0;
    my $pre_type = '';
    my $pre_info = '';
    my @read = ();
    my @RDrate = ();
    my %GT;
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
	    my @line = split (/\t/, ${${$vcf{$chr}}{$pos}}{$type});
	    my $len = $1 if ($line[7] =~ /SVLEN=(\d+)/);
	    my $end = $pos + $len -1;
	    my $read = $1 if ($line[7] =~ /READS=(\d+)/);
	    my $RD = $1 if ($line[7] =~ /RDRATIO=(.+?);/);
	    my $GT = $1 if ($line[7] =~ /GT=(.+)/);
	    if ((($pre_type ne '') and ($pre_type ne $type)) or (($pre_pos > 0) and ($pos - $pre_end > 2010))){
		if (@RDrate > 0){
		    if (@RDrate == 1){
			print "$pre_info\n";
		    }
		    else{
			@read = sort {$b <=> $a} @read;
			my $max_read = $read[0];
			my $sum_RD = 0;
			map{$sum_RD += $_} @RDrate;
			my $ave_RD = int ($sum_RD / @RDrate * 100000) / 100000;
			my $GTcons = '';
			foreach my $gt (sort {$GT{$b} <=> $GT{$a}} keys %GT){
			    $GTcons = $gt;
			    last;
			}
			my $chr2 = $chr;
			$chr2 =~ s/^0*// if ($chr2 =~ /^0/);
			print "$chr2\t$pre_pos\t$pre_type\t.\t.\t.\tPASS\tSVTYPE=$pre_type;SVLEN=$pre_len;RDRATIO=$ave_RD;READS=$max_read;GT=$GTcons\n";
		    }
		    @RDrate = ();
		    @read = ();
		    %GT = ();
		}
	    }
	    elsif (($pre_type ne '') and ($pre_type eq $type) and ($pre_pos > 0) and ($pos - $pre_end <= 2010)){
		if (($len >= 10000) or ($pre_len >= 10000)){
		    $pre_len = $end - $pre_pos + 1;
		    $pre_end = $end;
		    push @RDrate, $RD;
		    push @read, $read;
		    $GT{$GT} ++;
		    next;
		}
	    }
	    push @RDrate, $RD;
	    push @read, $read;
	    $GT{$GT} ++;
	    $pre_pos = $pos;
	    $pre_info = ${${$vcf{$chr}}{$pos}}{$type};
	    $pre_end = $end;
	    $pre_len = $len;
	    $pre_type = $type;
	}
    }
    if (@RDrate > 0){
	if (@RDrate == 1){
	    print "$pre_info\n";
	}
	else{
	    @read = sort {$b <=> $a} @read;
	    my $max_read = $read[0];
	    my $sum_RD = 0;
	    map{$sum_RD += $_} @RDrate;
	    my $ave_RD = int ($sum_RD / @RDrate * 100000) / 100000;
	    my $GTcons = '';
	    foreach my $gt (sort {$GT{$b} <=> $GT{$a}} keys %GT){
		$GTcons = $gt;
		last;
	    }
	    my $chr2 = $chr;
	    $chr2 =~ s/^0*// if ($chr2 =~ /^0/);
	    print "$chr2\t$pre_pos\t$pre_type\t.\t.\t.\tPASS\tSVTYPE=$pre_type;SVLEN=$pre_len;RDRATIO=$ave_RD;READS=$max_read;GT=$GTcons\n";
	}
	@RDrate = ();
	@read = ();
	%GT = ();
    }
}
