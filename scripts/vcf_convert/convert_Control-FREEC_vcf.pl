#!/usr/bin/perl -w
use strict;

my $var_file = shift @ARGV;

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
    my $end = $line[2];
    my $len = $end - $pos + 1;
    my $copy = $line[3];
    my $type = 'DEL' if ($line[4] eq 'loss');
    $type = 'DUP' if ($line[4] eq 'gain');
    next if ($line[4] eq 'normal');
    my $reads = 0;
    my $score = $line[6];
    if ($score < 0){
	$reads = 3;
    }
    elsif ($score < 1){
	$reads = 4;
    }
    elsif ($score < 2){
	$reads = 5;
    }
    elsif ($score < 3){
	$reads = 6;
    }
    elsif ($score < 5){
	$reads = 7;
    }
    elsif ($score < 10){
	$reads = 8;
    }
    else{
	$reads = 10;
    }
	if (@line > 5){
		my $gt = $line[5];
		$gt = '1/1' if ($gt eq '-') or ($gt eq 'AABB') or ($gt eq 'AAABB') or ($gt eq 'AABBB') or ($gt eq 'AAABBB');
		$gt = '1/0' if ($gt eq 'A') or ($gt eq 'B') or ($gt eq 'AAB') or ($gt eq 'ABB') or ($gt =~ /^A{2,}$/) or ($gt =~ /^B{2,}$/) or ($gt eq 'AAAB') or ($gt eq 'ABBB');
		$gt = './.' if ($gt !~ /\//);
		print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;COPYNUM=$copy;READS=$reads;GT=$gt\n";
	}
	else{
		print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;COPYNUM=$copy;READS=$reads\n";
	}
}
close (FILE);
