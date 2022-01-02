#!/usr/bin/perl -w
use strict;

my $min_qual = 0;

my $min_ins_dist = 170;

my %var;
my @file = (@ARGV);
if ($ARGV[0] eq '-'){
    @file = ();
    while (my $file = <STDIN>){
	chomp $file;
	my @str = split (/\s+/, $file);
	push @file, @str;
    }
}

foreach my $var_file (@file){
    my $pre_line = '';
    my $pre_chr = '';
    my $pre_ins = 0;
    if (!-e $var_file){
	print STDERR "$var_file is not found\n";
	next;
    }
    open (FILE, $var_file);
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#/){
	    next;
	}
	if ($line =~ /^\w+/){
	    if ($pre_line ne ''){
		my @pre_line = split (/\t/, $pre_line);
		my $type = $pre_line[0];
		my $chr = $pre_line[3];
		my $pos = $pre_line[4];
		my $size = $pre_line[5] - $pos + 1;
		my $reads = 100;
		if ($type eq 'deletion'){
		    $type = 'DEL';
		}
		elsif ($type eq 'insertion'){
		    $type = 'INS';
		}
		elsif ($type eq 'duplication'){
		    $type = 'DUP';
		}
		elsif ($type eq 'inversion'){
		    $type = 'INV';
		}
		elsif ($type eq 'transposition'){
		    if ($size >= 170){
			$type = 'DEL';
		    }
		    else{
			$type = 'INS';
		    }
		}
		my $qual = $pre_line[2];
		if ($type eq 'unknown'){
		    $pre_chr = $chr;
		    $pre_line = $line;
		    next;
		}
		my $chr02d = $chr;
		$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
		if (($type ne 'INS') or (($type eq 'INS') and (abs ($pre_ins - $pos) >= $min_ins_dist) and ($qual >= $min_qual))){
		    ${$var{$chr02d}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$size;QUAL=$qual;READS=$reads\n";
		}
		$pre_chr = $chr;
		$pre_ins = $pos if ($type eq 'INS');
	    }
	    $pre_line = $line;
	}
	elsif ($pre_line ne ''){
	    my @pre_line = split (/\t/, $pre_line);
	    my $type = $pre_line[0];
	    my $chr = $pre_line[3];
	    my $pos = $pre_line[4];
	    my $size = $pre_line[5] - $pos + 1;
	    my $reads = $1 if ($line =~ /num=(\d+)/);
	    if ($type eq 'deletion'){
		$type = 'DEL';
	    }
	    elsif ($type eq 'insertion'){
		$type = 'INS';
	    }
	    elsif ($type eq 'duplication'){
		$type = 'DUP';
	    }
	    elsif ($type eq 'inversion'){
		$type = 'INV';
	    }
	    elsif ($type eq 'transposition'){
		if ($size >= 170){
		    $type = 'DEL';
		}
		else{
		    $type = 'INS';
		}
	    }
	    my $qual = $pre_line[2];
	    if ($type eq 'unknown'){
		$pre_chr = $chr;
		$pre_line = '';
		next;
	    }
	    my $chr02d = $chr;
	    $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	    if (($type ne 'INS') or (($type eq 'INS') and (abs ($pre_ins - $pos) >= $min_ins_dist) and ($qual >= $min_qual))){
		${$var{$chr02d}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$size;QUAL=$qual;READS=$reads\n";
	    }
	    $pre_line = '';
	    $pre_chr = $chr;
	    $pre_ins = $pos if ($type eq 'INS');
	}
    }
    if ($pre_line ne ''){
	my @pre_line = split (/\t/, $pre_line);
	my $type = $pre_line[0];
	my $chr = $pre_line[3];
	my $pos = $pre_line[4];
	my $size = $pre_line[5] - $pos + 1;
	my $reads = 100;
	if ($type eq 'deletion'){
	    $type = 'DEL';
	}
	elsif ($type eq 'insertion'){
	    $type = 'INS';
	}
	elsif ($type eq 'duplication'){
	    $type = 'DUP';
	}
	elsif ($type eq 'inversion'){
	    $type = 'INV';
	}
	elsif ($type eq 'transposition'){
	    if ($size >= 170){
		$type = 'DEL';
	    }
	    else{
		$type = 'INS';
	    }
	}
	my $qual = $pre_line[2];
	next if ($type eq 'unknown');
	my $chr02d = $chr;
	$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	if (($type ne 'INS') or (($type eq 'INS') and (abs ($pre_ins - $pos) >= $min_ins_dist) and ($qual >= $min_qual))){
	    ${$var{$chr02d}}{$pos} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$size;QUAL=$qual;READS=$reads\n";
	}
    }
    close (FILE);
}

foreach my $chr (sort keys %var){
    foreach my $pos (sort {$a <=> $b} keys %{$var{$chr}}){
	print ${$var{$chr}}{$pos};
    }
}
