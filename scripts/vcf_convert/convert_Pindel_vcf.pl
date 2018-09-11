#!/usr/bin/perl -w
use strict;

# covert Pindel output files to vcf

my @del_file;
my @dup_file;
my @ins_file;
my @large_ins_file;
my @inv_file;
my @tra_file;
my @out_prefix;
my %gt;
my $gt_flag = 0;

if (@ARGV == 0){
    @ARGV = <./*_*>;
}

foreach my $file (@ARGV){
    if ((-f $file) and ($file =~ /(.+)_D$/)){
	push @del_file, $file;
    push @out_prefix, $1;
    }
    elsif ((-f $file) and ($file =~ /_TD$/)){
	push @dup_file, $file;
    }
    elsif ((-f $file) and ($file =~ /_SI$/)){
	push @ins_file, $file;
    }
    elsif ((-f $file) and ($file =~ /_LI$/)){
	push @large_ins_file, $file;
    }
    elsif ((-f $file) and ($file =~ /_INV$/)){
	push @inv_file, $file;
    }
    elsif ((-f $file) and ($file =~ /_INT_final$/)){
	push @tra_file, $file;
    }
}
=pod
foreach my $prefix (@out_prefix){
    my $chr = '';
    $chr = $1 if ($prefix =~ /\.(chr[\dXY]+)/);
    my $convert_command = "/share/project/sg/kosugi/tools/pindel-master/pindel2vcf -P $prefix -r /home/lustre/sg/kosugi/SV_calls/Ref/hs37d5.$chr.fasta -R GRCh37.p5 -d 20110817 2>pindel2vcf.log" if ($chr ne '');
    $convert_command = "/share/project/sg/kosugi/tools/pindel-master/pindel2vcf -P $prefix -r /home/lustre/sg/kosugi/SV_calls/Ref/hs37d5.fa -R GRCh37.p5 -d 20110817 2>pindel2vcf.log" if ($chr eq '');
    system ("$convert_command");
    my $vcf_file = "$prefix.vcf";
    open (FILE, $vcf_file) or die "$vcf_file is not found:$!\n";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^#|^$/);
        my @line = split (/\t/, $line);
        my $chr = $line[0];
        my $pos = $line[1];
        my $gt = $1 if ($line[9] =~ /(.+?):/);
        my $type = $1 if ($line[7] =~ /SVTYPE=([^:;]+)/);
        ${${$gt{$chr}}{$pos}}{$type} = $gt;
        if (($gt ne '0/0') and ($gt ne './.')){
            $gt_flag = 1;
        }
    }
    close (FILE);
}
=cut
my $min_sv_len = 80;

my $min_reads = 2;

my %vcf;

my $flag = 0;

foreach my $del_file (@del_file){
    open (FILE, $del_file) or die "$del_file is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#+$/){
	    $flag = 1;
	    next;
	}
	if ($flag == 1){
	    my @line = split (/\t/, $line);
	    my $id = $line[0];
	    my $type = $line[1];
	    my $svsize = $1 if ($type =~ /D\s(\d+)/);
	    my $chr = $1 if ($line =~ /ChrID\s(\S+)/);
	    my $bp1 = $1 if ($line =~ /BP\s(\d+)/);
	    my $bp2 = $1 if ($line =~ /BP\s\d+\t(\d+)/);
	    my $reads = $1 if ($line =~ /Supports\s(\d+)/);
	    if (($svsize >= $min_sv_len) and ($reads >= $min_reads)){
		${${$vcf{$chr}}{$bp1}}{'DEL'} = "SVTYPE=DEL;SVLEN=$svsize;BP2=$bp2;READS=$reads;ID=D$id";
	    }
	    $flag = 0;
	}
    }
    close (FILE);
}

foreach my $dup_file (@dup_file){
    open (FILE, $dup_file) or die "$dup_file is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#+$/){
	    $flag = 1;
	    next;
	}
	if ($flag == 1){
	    my @line = split (/\t/, $line);
	    my $id = $line[0];
	    my $type = $line[1];
	    my $svsize = $1 if ($type =~ /TD\s(\d+)/);
	    my $chr = $1 if ($line =~ /ChrID\s(\S+)/);
	    my $bp1 = $1 if ($line =~ /BP\s(\d+)/);
	    my $bp2 = $1 if ($line =~ /BP\s\d+\t(\d+)/);
	    my $reads = $1 if ($line =~ /Supports\s(\d+)/);
	    if (($svsize >= $min_sv_len) and ($reads >= $min_reads)){
		${${$vcf{$chr}}{$bp1}}{'DUP'} = "SVTYPE=DUP;SVLEN=$svsize;BP2=$bp2;READS=$reads;ID=TD$id";
	    }
	    $flag = 0;
	}
    }
    close (FILE);
}

foreach my $inv_file (@inv_file){
    open (FILE, $inv_file) or die "$inv_file is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#+$/){
	    $flag = 1;
	    next;
	}
	if ($flag == 1){
	    my @line = split (/\t/, $line);
	    my $id = $line[0];
	    my $type = $line[1];
	    my $svsize = $1 if ($type =~ /INV\s(\d+)/);
	    my $chr = $1 if ($line =~ /ChrID\s(\S+)/);
	    my $bp1 = $1 if ($line =~ /BP\s(\d+)/);
	    my $bp2 = $1 if ($line =~ /BP\s\d+\t(\d+)/);
	    my $reads = $1 if ($line =~ /Supports\s(\d+)/);
	    if (($svsize >= $min_sv_len) and ($reads >= $min_reads)){
		${${$vcf{$chr}}{$bp1}}{'INV'} = "SVTYPE=INV;SVLEN=$svsize;BP2=$bp2;READS=$reads;ID=IV$id";
	    }
	    $flag = 0;
	}
    }
    close (FILE);
}

foreach my $ins_file (@ins_file){
    open (FILE, $ins_file) or die "$ins_file is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#+$/){
	    $flag = 1;
	    next;
	}
	if ($flag == 1){
	    my @line = split (/\t/, $line);
	    my $id = $line[0];
	    my $type = $line[1];
	    my $svsize = $1 if ($type =~ /I\s(\d+)/);
	    my $chr = $1 if ($line =~ /ChrID\s(\S+)/);
	    my $bp1 = $1 if ($line =~ /BP\s(\d+)/);
	    my $bp2 = $1 if ($line =~ /BP\s\d+\t(\d+)/);
	    my $reads = $1 if ($line =~ /Supports\s(\d+)/);
	    if (($svsize >= $min_sv_len) and ($reads >= $min_reads)){
		${${$vcf{$chr}}{$bp1}}{'INS'} = "SVTYPE=INS;SVLEN=$svsize;BP2=$bp2;READS=$reads;ID=I$id";
	    }
	    $flag = 0;
	}
    }
    close (FILE);
}

foreach my $large_ins_file (@large_ins_file){
    open (FILE, $large_ins_file) or die "$large_ins_file is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#+$/){
	    $flag = 1;
	    next;
	}
	if ($flag == 1){
	    my @line = split (/\t/, $line);
		my $id = $line[0];
		my $type = $line[1];
		my $chr = $1 if ($line =~ /ChrID\s(\S+)/);
		my $bp1 = $line[3];
		my $item = $line[7];
		my $read1 = $1 if ($item =~ /\-\s(\d+)$/);
		my $read2 = $1 if ($item =~ /\+\s(\d+)\s/);
		my $reads = $read1 + $read2;
		${${$vcf{$chr}}{$bp1}}{'INS'} = "SVTYPE=INS;READS=$reads;ID=LI$id" if ($reads >= $min_reads);
	    $flag = 0;
	    next;
	}
    }
    close (FILE);
}

if (@tra_file > 0){
    foreach my $tra_file (@tra_file){
	if ($tra_file ne ''){
	    open (FILE, $tra_file) or die "$tra_file is not found: $!\n";
	    while (my $line = <FILE>){
		chomp $line;
		my @line = split (/\s+/, $line);
		my $type = 'TRA';
		my $chr = $line[1];
		my $bp1 = $line[3];
		my $chr2 = $line[5];
		my $pos2 = $line[7];
		my $reads = $line[11];
		next if ($reads < $min_reads);
		next if (exists ${${$vcf{$chr2}}{$pos2}}{'TRA'});
		${${$vcf{$chr}}{$bp1}}{'TRA'} = "SVTYPE=TRA;READS=$reads\tCHR2=$chr2;POS2=$pos2" if ($reads >= $min_reads);
		next;
	    }
	    close (FILE);
	}
    }
}

if ($gt_flag == 1){
    foreach my $chr (sort keys %vcf){
        next if ($chr !~ /^\d+$|[XY]/);
        foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
        foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
            my $gt = './.';
            if (exists ${${$gt{$chr}}{$pos}}{$type}){
                my $gt1 = ${${$gt{$chr}}{$pos}}{$type};
                if ($gt1 eq '0/1'){
                    $gt = $gt1;
                }
                elsif ($gt1 eq '1/1'){
                    $gt = $gt1;
                }
            }
            print "$chr\t$pos\t$type\t.\t.\t.\tPASS\t${${$vcf{$chr}}{$pos}}{$type};GT=$gt\n";
        }
        }
    }
}
else{
    foreach my $chr (sort keys %vcf){
        next if ($chr !~ /^\d+$|[XY]/);
        foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
        foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
            print "$chr\t$pos\t$type\t.\t.\t.\tPASS\t${${$vcf{$chr}}{$pos}}{$type}\n";
        }
        }
    }
}
