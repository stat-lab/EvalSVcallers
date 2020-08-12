#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use FindBin qw/$Bin $Script/;

# First, select SVs (shared set) shared between 2 SV sets (A AND B => AB, C AND D => CD), second merge different two shared sets (AB OR CD)

my $eval_file = '';

my $sv_type = 'ALL';

my $sv_size = 'ALL';

my $target_chr = 'all';

my $min_sv_len = 50;

my $max_sv_len = 2000000;

my $max_inv_size = 200000;

my $read_length = 125;
my $var_sd = 125;
my $ins_sd = 200;

my $min_qual = 90;
my $min_del_qual = 0;

my $min_reads_parent = 2;

my $min_overlap_ratio = 0.6;        # for overlap between SV call sets

my $min_overlap_ratio2 = 0.5;       # for overlap with ref SVs

my $min_overlap = 0;
my $min_overlap2 = 0;

my $help;

my $ref_sv = 'N';

my $script_path = $Bin;

my $ref_sv_simA = "$script_path/../Ref_SV/Sim-A.SV.vcf";

my $ref_sv_NA12878 = "$script_path/../Ref_SV/NA12878_DGV-2016_LR-assembly.vcf";

my $SD_len_file = "$script_path/SV_SDlen_list_2.txt";

my $gap_bed = "$script_path/../Ref_SV/gap.bed";

my $out_prefix = '';

my $input_dir = './';

my $parent_dir = '';

my $include_y = 0;

my $output_vcf = 0;

my @tools;

my @parent_vcf = ();

my $merge_2set = 0;
my $merge_3set = 0;


GetOptions(
    'ref|r=s' => \$ref_sv,
    'tools|t=s{,}' => \@tools,
    'parent_vcf|pv=s{,}' => \@parent_vcf,
    'parent_dir|pd=s' => \$parent_dir,
    'sv_type|st=s' => \$sv_type,
    'sv_size|ss=s' => \$sv_size,
    'min_ins|mins=i' => \$ins_sd,
    'min_ovl|mo=f' => \$min_overlap,
    'min_ovl2|mo2=f' => \$min_overlap2,
    'read_len|rl=i' => \$read_length,
    'dir|d=s' => \$input_dir,
    'prefix|p=s' => \$out_prefix,
    'chr|c=s' => \$target_chr,
    'in_y|y' => \$include_y,
    'vcf_out|vcf' => \$output_vcf,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  evaluate_SV_overlap_calls.pl <option> [vcf file] 

  Options:
   --ref or -r <STR>        reference SV type  (A|N), A: Sim-A, N: NA12878 [default: N]
   --tools or -t <STR>      list of tool:read (e.g., Pindel:3 Lumpy:5)
   --sv_type or -st <STR>   SV type (ALL|DEL|DUP|INS|INV|TRA) [default: ALL]
   --sv_size or -ss <STR>   SV size, SS, S, M, L, Lo, or ALL, SS: 30-100 bp, S: 50-1000 bp, M: 1000-100,000 bp, Lo: 500-2,000,000 bp, L: 100,000-2,000,000 bp [default: ALL]
   --min_size or -s <INT>   minimum size (bp) of SV [default: 50]
   --chr or -c <STR>        target chromosome to be analyzed [all or chr name(s), e.g., 4,5,6,X; default: all]
   --min_ovl or -mo <FLOAT> minimum rate of reciprocal overlap between called non-INS-SVs and reference non-INS-SVs [default: 0.5 for NA12878, 0.8 for Sim-A > 1Kb SVs, 0.6 for Sim-A <= 1 Kb SVs]
   --min_ovl2 or -mo2 <FLOAT> minimum rate of reciprocal overlap between non-INS-SVs from 2 datasets [default: 0.6]
   --min_ins or -mins <INT> maximum allowable length between overlap INS breakpoints [default: 200]
   --read_len or -rl <INT>  read length [default: 125]
   --dir or -d <STR>        directory of input vcf files [default: ./]
   --parent_vcf or -pv <STR> list of parent vcf files (e.g., Pindel.mother.vcf Pindel.father.vcf Lumpy.mother.vcf Lumpy.father.vcf) [either specify -pd or -pv when using trio data]
   --parent_dir or -pd <STR> directory containing parent vcf files (should not include child vcf files in this directory)
   --in_y or -y             include chrY [default: false]
   --prefix or -p <STR>     prefix of output
   --vcf_out or -vcf        output vcf file for overlapping calls between tool pairs [default: false]
   --help or -h             output help message
   
=cut

die "eval file (--eval_file) or tool list (--tools) not specified: \n" if ($eval_file eq '') and (@tools == 0);
die "Double-specification of both -m2 and -m3 is not allowed: \n" if ($merge_2set == 1) and ($merge_3set == 1);

my $ref_type = $ref_sv;

if ($ref_sv eq 'A'){
    $ref_sv = $ref_sv_simA;
}
elsif ($ref_sv eq 'N'){
    $ref_sv = $ref_sv_NA12878;
#    $min_overlap_ratio = 0.01 if ($sv_type eq 'DUP');
}

my $min_overlap_ratio3 = 0.4;           # for overlap between parent-child SVs

if ($min_overlap2 == 0){
    $min_overlap_ratio2 = 0.5 if ($ref_type eq 'N');
    $min_overlap_ratio2 = 0.6 if ($ref_type eq 'A') and (($sv_size eq 'SS') or ($sv_size eq 'S'));
    $min_overlap_ratio2 = 0.8 if ($ref_type eq 'A') and (($sv_size eq 'M') or ($sv_size eq 'L') or ($sv_size eq 'ALL'));
}
else{
    $min_overlap_ratio2 = $min_overlap2;
}

if ($min_overlap == 0){
    $min_overlap_ratio = 0.6;
}
else{
    $min_overlap_ratio = $min_overlap;
}

my %dup_tool = ('AS-GENSENG' => 1, 'CNVnator' => 1, 'Control-FREEC' => 1, 'Delly' => 1, 'ERDS' => 1, 'forestSV' => 1, 'inGAP' => 1, 'laSV' => 1, 'Lumpy' => 1, 'Manta' => 1, 'Meerkat' => 1, 'metaSV' => 1, 'OncoSNP-Seq' => 1,
		'Pindel' => 1, 'RAPTR' => 1, 'readDepth' => 1, 'Sniffles' => 1, 'SoftSearch' => 1, 'SoftSV' => 1, 'SVDetect' => 1, 'Ulysses' => 1, 'TIDDIT' => 1, 'PennCNV-Seq' => 1, 'BICseq2' => 1, 'SVelter' => 1, 'Wham' => 1, 'GRIDSS' => 1);
my %ins_tool = ('123SV' => 1, 'Basil' => 1, 'BreakDancer' => 1, 'BreakSeek' => 1, 'BreakSeq2' => 1, 'CLEVER' => 1, 'CREST' => 1, 'FermiKit' => 1, 'inGAP' => 1, 'laSV' => 1, 'Manta' => 1, 'Meerkat' => 1, 'metaSV' => 1, 'MELT' => 1, 
		'MindTheGap' => 1, 'Mobster' => 1, 'PBHoney' => 1, 'PBSV' => 1, 'Popins' => 1, 'PRISM' => 1, 'RAPTR' => 1, 'Sniffles' => 1, 'Socrates' => 1, 'SoftSearch' => 1, 'SVDetect' => 1, 'SVfinder' => 1, 'SVseq2' => 1, 'laSV' => 1, 'Pamir' => 1, 'Wham' => 1, 'GRIDSS' => 1);

my %SDlen;
my %SDbp;

open (FILE, $SD_len_file) or die ("SD_len_file is not found: $!\n");
while (my $line = <FILE>){
    chomp $line;
    my @line = split (/\s+/, $line);
    my $id = lc $line[0];
    my $class = $line[1];
    my $type = $line[2];
    my $sd_len = $line[3];
    my $sd_bp = $line[4];
    if ($type eq 'INS'){
	${${$SDlen{$type}}{'S'}}{$id} = $sd_len;
	${${$SDbp{$type}}{'S'}}{$id} = $sd_bp;
	${${$SDlen{$type}}{'M'}}{$id} = $sd_len;
	${${$SDbp{$type}}{'M'}}{$id} = $sd_bp;
	${${$SDlen{$type}}{'L'}}{$id} = $sd_len;
	${${$SDbp{$type}}{'L'}}{$id} = $sd_bp;
    }
    else{
	${${$SDlen{$type}}{$class}}{$id} = $sd_len;
	${${$SDbp{$type}}{$class}}{$id} = $sd_bp;
    }
    if (($type eq 'DEL') or ($type eq 'INV')){
	if ($class eq 'M'){
	    ${${$SDlen{$type}}{'L'}}{$id} = $sd_len;
	    ${${$SDbp{$type}}{'L'}}{$id} = $sd_bp;
	}
    }
}
close (FILE);

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

my $call_ctx_num = 0;
my $call_ctx_s = 0;
my $call_ctx_m = 0;
my $call_ctx_l = 0;

#my %tool_set;
my %tool_single;
my @tool_set;
my %var_file;
my %parent1_file;
my %parent2_file;
my %var_read;
my @var_file;

my %call;
my %sv_read;
my %match;
my %call_count;
my %match_count;
my %dup_flag;
my @var_id2 = ();

my @work_files = <$input_dir/*.vcf>;
my @work_files2;
foreach (@work_files){
    push @work_files2, $_ if ($_ !~ /\.TP-FP\.vcf$|\.FP\.vcf$|\.FN\.vcf$/);
}

if (@tools > 0){
    $eval_file = "$out_prefix.tool.list" if ($out_prefix ne '');
    $eval_file = "$sv_type.$sv_size.tool.list" if ($out_prefix eq '');
    my $count_1 = 0;
    my @tool_comb;
    my @tool_comb2;
    my @tool_comb3;
    foreach my $tool_1 (@tools){
	push @tool_comb, "$tool_1=$tool_1";
    }
    foreach my $tool_1 (@tools){
	$count_1 ++;
	my $count_2 = 0;
	foreach my $tool_2 (@tools){
	    $count_2 ++;
	    next if ($count_2 <= $count_1);
	    my $set = "$tool_1=$tool_2";
	    push @tool_comb, $set;
	}
    }
    if (($merge_2set == 0) and ($merge_3set == 0)){
	open (OUT, "> $eval_file");
	foreach my $set (@tool_comb){
	    print OUT "$set\n";
	}
	close (OUT);
    }
    elsif (($merge_2set == 1) or ($merge_3set == 1)){
	$count_1 = 0;
	foreach my $set1 (@tool_comb){
	    $count_1 ++;
	    my $count_2 = 0;
	    foreach my $set2 (@tool_comb){
		$count_2 ++;
		next if ($count_2 <= $count_1);
		my $set3 = "$set1~$set2";
		push @tool_comb2, $set3;
	    }
	}
	if ($merge_2set == 1){
	    open (OUT, "> $eval_file");
	    foreach my $set (@tool_comb2){
		print OUT "$set\n";
	    }
	    close (OUT);
	}
    }
    if ($merge_3set == 1){
	$count_1 = 0;
	foreach my $set1 (@tool_comb2){
	    foreach my $set2 (@tool_comb){
		next if ($set1 =~ $set2);
		my $set3 = "$set1~$set2";
		push @tool_comb3, $set3;
	    }
	}
	open (OUT, "> $eval_file");
	foreach my $set (@tool_comb3){
	    print OUT "$set\n";
	}
	close (OUT);
    }
}

$eval_file =~ s/\'//g if ($eval_file =~ /\'/);

@parent_vcf = <$parent_dir/*.vcf> if ($parent_dir ne '') and (@parent_vcf == 0);

open (FILE, $eval_file) or die "$eval_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#/);
    next if ($line =~ /^$/);
    my $set = $line;
    $set =~ s/~/|/g if ($set =~ /~/);
    my @pairs = split (/\|/, $set);
    push @tool_set, $set;
    foreach my $pair (@pairs){
        my ($id1, $id2) = ('', '');
        ($id1, $id2) = split (/=/, $pair) if ($pair =~ /=/);
        $id1 = $pair if ($pair !~ /=/);
        my ($id1_base, $read1) = split (/:/, $id1);
        ${$var_read{$id1_base}}{$read1} = 1;
        if (!exists $var_file{$id1_base}){
            my $exist_flag = 0;
            foreach my $file (@work_files2){
                if (($file =~ /\/$id1_base\./i) and ($file !~ /=|~/)){
                    $var_file{$id1_base} = $file;
                    $exist_flag = 1;
                    last;
                }
            }
            if ($exist_flag == 0){
                die "vcf file for $id1_base is not found:\n";
            }
        }
        if (@parent_vcf > 0){
            my $parent1_id1 = '';
            my $parent2_id1 = '';
            my $parent1_vcf = '';
            my $parent2_vcf = '';
            foreach my $vcf (@parent_vcf){
                if ($vcf =~ /.+\/($id1_base.+)/){
                    my $vcf_base = $1;
                    if ($parent1_id1 eq ''){
                        $parent1_id1 = $vcf_base;
                        $parent1_vcf = $vcf;
                    }
                    elsif ($parent2_id1 eq ''){
                        $parent2_id1 = $vcf_base;
                        $parent2_vcf = $vcf;
                        last;
                    }
                }
            }

            if ($parent2_id1 ne ''){
                $parent1_file{$id1_base} = $parent1_vcf;
                $parent2_file{$id1_base} = $parent2_vcf;
            }
        }
        if ($id2 ne ''){
            my ($id2_base, $read2) = split (/:/, $id2);
            ${$var_read{$id2_base}}{$read2} = 1;
            if (!exists $var_file{$id2_base}){
                my $exist_flag = 0;
                foreach my $file (@work_files2){
                    if (($file =~ /\/$id2_base\./i) and ($file !~ /=|~/)){
                        $var_file{$id2_base} = $file;
                        $exist_flag = 1;
                        last;
                    }
                }
                if ($exist_flag == 0){
                    die "vcf file for $id2_base is not found:\n";
                }
            }
            if (@parent_vcf > 0){
                my $parent1_id2 = '';
                my $parent2_id2 = '';
                my $parent1_vcf = '';
                my $parent2_vcf = '';
                foreach my $vcf (@parent_vcf){
                    if ($vcf =~ /.+\/($id2_base.+)/){
                        my $vcf_base = $1;
                        if ($parent1_id2 eq ''){
                            $parent1_id2 = $vcf_base;
                            $parent1_vcf = $vcf;
                        }
                        elsif ($parent2_id2 eq ''){
                            $parent2_id2 = $vcf_base;
                            $parent2_vcf = $vcf;
                            last;
                        }
                    }
                }
                if ($parent2_id2 ne ''){
                    $parent1_file{$id2_base} = $parent1_vcf;
                    $parent2_file{$id2_base} = $parent2_vcf;
                }
            }
        }
        else{
            $tool_single{$pair} = 1;
        }
    }
}
close (FILE);

my %ref;
my %ref_ins;
my %ref_dup;
my %ref_num;

my $min_ref_len = 30;
my $max_ref_len = 2000000;

open (FILE, $ref_sv) or die "$ref_sv is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#/);
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    next if ($target_chr ne 'all') and (($chr ne $target_chr) and ($target_chr !~ /,$chr,|,$chr$|^$chr,/));
    next if (($chr eq 'Y')) and ($include_y == 0);
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    my $type2 = $line[2];
#    $type = 'DEL' if ($type eq 'INDEL');
    $type = 'DUP' if ($type2 eq 'tandem');
    my $svlen = 0;
    $svlen = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
    next if ($svlen < $min_ref_len) and ($type ne 'INS');
    next if ($svlen > $max_ref_len);
    next if ($type2 eq 'INS-DUP');
    my $end = $pos + $svlen;
    my $gap_overlap = 0;
    foreach my $gstart (sort {$a <=> $b} keys %{$gap{$chr}}){
	my $gend = ${$gap{$chr}}{$gstart};
	if (($pos >= $gstart) and ($pos <= $gend)){
	    if ($type eq 'INS'){
		$gap_overlap = 1;
		last;
	    }
	    else{
		if ($gend - $pos + 1 >= $svlen * 0.5){
		    $gap_overlap = 1;
		    last;
		}
	    }
	}
	elsif (($type ne 'INS') and ($gstart >= $pos) and ($gstart <= $end)){
	    if ($end - $gstart + 1 >= $svlen * 0.5){
		$gap_overlap = 1;
		last;
	    }
	}
	last if ($gstart > $end);
    }
    if ($gap_overlap == 1){
	next;
    }
    if ($type eq 'DEL'){
	if ($sv_size eq 'SS'){
	    $ref_num{$type} ++ if ($svlen <= 100);
	}
	elsif ($sv_size eq 'S'){
	    $ref_num{$type}++ if ($svlen > 100) and ($svlen <= 1000);
	}
	elsif ($sv_size eq 'M'){
	    $ref_num{$type} ++ if ($svlen > 1000) and ($svlen <= 100000);
	}
	elsif ($sv_size eq 'L'){
	    $ref_num{$type} ++ if ($svlen > 100000);
	}
	else{
	    $ref_num{$type} ++;
	}
    }
    elsif ($type ne 'INS'){
	if ($sv_size eq 'S'){
	    $ref_num{$type} ++ if ($svlen <= 1000);
	}
	elsif ($sv_size eq 'M'){
	    $ref_num{$type} ++ if ($svlen > 1000) and ($svlen <= 100000);
	}
	elsif ($sv_size eq 'L'){
	    $ref_num{$type} ++ if ($svlen > 100000);
	}
	else{
	    $ref_num{$type} ++;
	}
    }
    $ref_num{$type} ++ if ($type eq 'INS');
    ${${$ref{$type}}{$chr}}{$pos} = $svlen;
    
    if (($type eq 'INS') or ($type eq 'MEI')){
        ${$ref_ins{$chr}}{$pos} = $svlen;
        ${$ref_ins{$chr}}{$pos - $svlen} = $svlen;
    }
    elsif ($type eq 'DUP'){
        ${$ref_dup{$chr}}{$pos} = $svlen;
        ${$ref_dup{$chr}}{$pos - $svlen} = $svlen;
    }
}
close (FILE);

my %parent1;
my %parent2;

foreach my $id (keys %parent1_file){
    my $vcf = $parent1_file{$id};
    my $vcf_base = basename ($vcf);
    my %min_reads_parent;
#print STDERR "$id\tparent vcf: $vcf\n";
    if ($vcf_base =~ /PBHoney|Sniffles/){
        $min_reads_parent = 2;
    }
    elsif ($vcf_base =~ /RAPTR/){
        $min_reads_parent = 0;
        $min_reads_parent{'DEL'} = 8;
        $min_reads_parent{'DUP'} = 4;
        $min_reads_parent{'INS'} = 4;
    }
    elsif ($vcf_base =~ /SVfinder/){
        $min_reads_parent = 0;
        $min_reads_parent{'DEL'} = 5;
        $min_reads_parent{'INS'} = 4;
        $min_reads_parent{'INV'} = 5;
    }
    open (FILE, $vcf) or die "$vcf is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        if ($line =~ /^#|^$/){
            next;
        }
        my @line = split (/\t/, $line);
        my $chr = $line[0];
	$chr =~ s/^chr// if ($chr =~ /^chr/);
        next if ($target_chr ne 'all') and (($chr ne $target_chr) and ($target_chr !~ /,$chr,|,$chr$|^$chr,/));
        next if ($chr !~ /^\d+$|[XY]/);
        next if ($chr eq 'Y') and ($ref_type eq 'N');
        my $pos = $line[1];
        my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
        $type = 'INS' if ($type =~ /MEI|NUMT|VEI/);
        if ($sv_type ne 'ALL'){
            next if ($type ne $sv_type);
        }
        my $len = 0;
        $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
        next if ($type eq 'INV') and ($len > $max_inv_size);
        if (($type eq 'DEL') or ($type eq 'DUP') or ($type eq 'INV')){
            if ($sv_size eq 'SS'){
                next if ($len > 120);
            }
            elsif ($sv_size eq 'S'){
                next if ($len > 1200);
                next if ($len < 80);
            }
            elsif ($sv_size eq 'M'){
                next if ($len < 800);
                next if ($len > 120000);
            }
            elsif ($sv_size eq 'Lo'){
                next if ($len < 400);
            }
            elsif ($sv_size eq 'L'){
                next if ($len < 80000);
            }
        }
        if (($type eq 'INS') and ($line[7] =~ /QUAL=(\d+)/)){
            next if ($1 < $min_qual);
        }
        if (($type eq 'DEL') and ($line[7] =~ /QUAL=(\d+)/)){
            next if ($1 < $min_del_qual);
        }
        my $end = 0;
        if ($line[7] =~ /END=(\d+);/){
            $end = $1;
        }
        else{
            $end = $pos + $len - 1;
        }
        my $reads = $1 if ($line[7] =~ /READS=(\d+)/);
        
        next if ($reads < $min_reads_parent);
        next if (exists $min_reads_parent{$type}) and ($reads < $min_reads_parent{$type});
        ${${${$parent1{$id}}{$type}}{$chr}}{$pos} = $len;
    }
    close (FILE);
}
foreach my $id (keys %parent2_file){
    my $vcf = $parent2_file{$id};
    my $vcf_base = basename ($vcf);
    my %min_reads_parent;
    if ($vcf_base =~ /PBHoney|Sniffles/){
        $min_reads_parent = 2;
    }
    elsif ($vcf_base =~ /RAPTR/){
        $min_reads_parent = 0;
        $min_reads_parent{'DEL'} = 8;
        $min_reads_parent{'DUP'} = 4;
        $min_reads_parent{'INS'} = 4;
    }
    elsif ($vcf_base =~ /SVfinder/){
        $min_reads_parent = 0;
        $min_reads_parent{'DEL'} = 5;
        $min_reads_parent{'INS'} = 4;
        $min_reads_parent{'INV'} = 5;
    }
    open (FILE, $vcf) or die "$vcf is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        if ($line =~ /^#|^$/){
            next;
        }
        my @line = split (/\t/, $line);
        my $chr = $line[0];
	$chr =~ s/^chr// if ($chr =~ /^chr/);
        next if ($target_chr ne 'all') and (($chr ne $target_chr) and ($target_chr !~ /,$chr,|,$chr$|^$chr,/));
        next if ($chr !~ /^\d+$|[XY]/);
        next if ($chr eq 'Y') and ($ref_type eq 'N');
        my $pos = $line[1];
        my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
        $type = 'INS' if ($type =~ /MEI|NUMT|VEI/);
        if ($sv_type ne 'ALL'){
            next if ($type ne $sv_type);
        }
        my $len = 0;
        $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
        next if ($type eq 'INV') and ($len > $max_inv_size);
        if (($type eq 'DEL') or ($type eq 'DUP') or ($type eq 'INV')){
            if ($sv_size eq 'SS'){
                next if ($len > 120);
            }
            elsif ($sv_size eq 'S'){
                next if ($len > 1200);
                next if ($len < 80);
            }
            elsif ($sv_size eq 'M'){
                next if ($len < 800);
                next if ($len > 120000);
            }
            elsif ($sv_size eq 'Lo'){
                next if ($len < 400);
            }
            elsif ($sv_size eq 'L'){
                next if ($len < 80000);
            }
        }
        if (($type eq 'INS') and ($line[7] =~ /QUAL=(\d+)/)){
            next if ($1 < $min_qual);
        }
        if (($type eq 'DEL') and ($line[7] =~ /QUAL=(\d+)/)){
            next if ($1 < $min_del_qual);
        }
        my $end = 0;
        if ($line[7] =~ /END=(\d+);/){
            $end = $1;
        }
        else{
            $end = $pos + $len - 1;
        }
        my $reads = $1 if ($line[7] =~ /READS=(\d+)/);
        next if ($reads < $min_reads_parent);
        next if  (exists $min_reads_parent{$type}) and ($reads < $min_reads_parent{$type});
        ${${${$parent2{$id}}{$type}}{$chr}}{$pos} = $len;
    }
    close (FILE);
}

my $hg19_flag = 0;

foreach my $id (keys %var_file){
    my $var_file = $var_file{$id};
    my %pre_info;
    my %overlap;
    my %match_ref;
    $overlap{'DEL'} = 0;
    $overlap{'DUP'} = 0;
    $overlap{'INV'} = 0;
    $overlap{'INS'} = 0;
    $overlap{'TRA'} = 0;
    foreach my $min_read (sort {$a <=> $b} keys %{$var_read{$id}}){
	my $id2 = "$id:$min_read";
	push @var_id2, $id2;
	open (FILE, $var_file) or die "$var_file is not found: $!\n";
	while (my $line = <FILE>){
	    chomp $line;
	    next if ($line =~ /^#/);
	    my @line = split (/\t/, $line);
	    my $chr = $line[0];
	    $hg19_flag = 1 if ($chr =~ /^chr[\dXY]$/);
	    $chr =~ s/^chr// if ($chr =~ /^chr/);
	    
	    next if ($chr !~ /^\d+$|[XY]/);
	    next if ($target_chr ne 'all') and (($chr ne $target_chr) and ($target_chr !~ /,$chr,|,$chr$|^$chr,/));
            next if (($chr eq 'Y')) and ($include_y == 0);
	    my $pos = $line[1];
	    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
	    $type = 'INS' if ($type =~ /MEI|NUMT|VEI/);
	    if ($sv_type ne 'ALL'){
		next if ($type ne $sv_type);
	    }
            next if ($type !~ /^DEL$|^DUP$|^INS$|^INV$/);
	    my $len = 0;
	    $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
	    next if ($type eq 'INV') and ($len > $max_inv_size);
	    if (($type eq 'DEL') or ($type eq 'DUP') or ($type eq 'INV')){
		if ($sv_size eq 'SS'){
		    next if ($len > 100);
		}
		elsif ($sv_size eq 'S'){
		    next if ($len > 1000);
		    next if ($len < 100);
		}
		elsif ($sv_size eq 'M'){
		    next if ($len < 1000);
		    next if ($len > 100000);
		}
		elsif ($sv_size eq 'Lo'){
		    next if ($len < 500);
		}
		elsif ($sv_size eq 'L'){
		    next if ($len < 100000);
		}
	    }
            if (($type eq 'INS') and ($line[7] =~ /QUAL=(\d+)/)){
                next if ($1 < $min_qual);
            }
            if (($type eq 'DEL') and ($line[7] =~ /QUAL=(\d+)/)){
                next if ($1 < $min_del_qual);
            }
	    my $bplen = 0;
	    $bplen = $1 if ($line[7] =~ /BPSIZE=(\d+);/);
	    my $end = 0;
	    if ($bplen > 0){
		$end = $pos + $bplen - 1;
	    }
	    else{
		$end = $pos + $len - 1;
	    }
	    my $gap_overlap = 0;
	    foreach my $gstart (sort {$a <=> $b} keys %{$gap{$chr}}){
		my $gend = ${$gap{$chr}}{$gstart};
		if (($pos >= $gstart) and ($pos <= $gend)){
		    if ($type eq 'INS'){
			$gap_overlap = 1;
			last;
		    }
		    else{
			if ($gend - $pos + 1 >= $len * 0.5){
			    $gap_overlap = 1;
			    last;
			}
		    }
		}
		elsif (($type ne 'INS') and ($gstart >= $pos) and ($gstart <= $end)){
		    if ($end - $gstart + 1 >= $len * 0.5){
			$gap_overlap = 1;
			last;
		    }
		}
		last if ($gstart > $end);
	    }
	    if ($gap_overlap == 1){
		next;
	    }
	    
	    ${$call_count{$id2}}{$type} = 0 if (!exists ${$call_count{$id2}}{$type});
	    my $reads = 3;
	    $reads = $1 if ($line[7] =~ /READS=(\d+)/);
	    next if ($reads < $min_read);
	    my $pos2 = 0;
	    next if ($len < $min_sv_len) and ($len > 0) and ($type ne 'INS') and ($type ne 'TRA');
	    next if ($len > $max_sv_len);
	    my ($pre_chr, $pre_pos, $pre_size, $pre_hit) = ('', 0, 0, 0);
	    ($pre_chr, $pre_pos, $pre_size, $pre_hit) = split (/=/, $pre_info{$type}) if (exists $pre_info{$type});
	    if (($type eq 'INS') and ($pre_chr eq $chr) and ($pos - $pre_pos < $var_sd)){
		if ((($len > 0) and ($pre_size > $var_sd * 0.5) and ($len > $var_sd * 0.5)) or ($len == 0)){
		    $overlap{$type} ++;
		    if ($pre_hit == 1){
			$pre_info{$type} = "$chr=$pos=$len=1" if ($pos + $len > $pre_pos + $pre_size);
			next;
		    }
		}
	    }
	    elsif (($type ne 'INS') and ($pre_chr eq $chr) and ($pre_pos + $pre_size - $pos > $len * 0.8) and ($pre_pos + $pre_size - $pos > $pre_size * 0.8) and ($pre_size > $var_sd * 0.5) and ($len > $var_sd * 0.5)){
		$overlap{$type} ++;
		if ($pre_hit == 1){
		    $pre_info{$type} = "$chr=$pos=$len=1" if ($pos + $len > $pre_pos + $pre_size);
		    next;
		}
	    }
	    else{
		${${${$call{$id2}}{$type}}{$chr}}{$pos} = $len;
		${${${$sv_read{$id2}}{$type}}{$chr}}{$pos} = $reads;
		${$call_count{$id2}}{$type} ++;
	    }
	}
    }
}

my %call_share_2;

my %result_2p;

my %max_or;
my %max_and;
my %sv_and;
my %sv_and_or;
my %sv_and_or_or;

my %dup_and_flag;

my %classify_1;
my %classify_2;

foreach my $id1 (keys %call_count){
    my ($id2, $min_read) = split (/:/, $id1);
    foreach my $type (keys %{$call_count{$id1}}){
	if (exists $tool_single{$id1}){
	    foreach my $chr (keys %{${$call{$id1}}{$type}}){
		foreach my $pos (keys %{${${$call{$id1}}{$type}}{$chr}}){
		    my $len = ${${${$call{$id1}}{$type}}{$chr}}{$pos};
		    ${${${$sv_and{$type}}{$id1}}{$chr}}{$pos} = "$pos|$len";
		}
	    }
	}
    }
}

my $vcf_out = '';
$vcf_out = "$out_prefix.$ref_type.out" if ($out_prefix ne '');
open (OUT, "> $vcf_out") if ($vcf_out ne '');
print "#SV_type (size)\tRecall(%)\tPrecision(%)\tNum_overlapped_SVs\n";
print OUT "#SV_type (size)\tRecall(%)\tPrecision(%)\tNum_overlapped_SVs\n" if ($vcf_out ne '');

foreach my $set (@tool_set){
    my %vcf;
    my $set2 = $set;
    $set2 =~ s/\|/~/g;
    if ($set2 !~ /~/){
	my ($id1, $id2) = split (/=/, $set2);
	if ($id1 eq $id2){
	    $set2 = $id1;
	}
    }
    my $overlap_vcf = '';
    if (($output_vcf == 1) and ($set2 =~ /=/)){
        my $set3 = $set2;
        $set3 =~ s/:\d+//g;
        $overlap_vcf = "$out_prefix.$set3.vcf" if ($out_prefix ne '');
        $overlap_vcf = "$set3.vcf" if ($out_prefix eq '');
        open (OUT2, "> $overlap_vcf");
    }
    my %call_share;
    %sv_and_or = ();
    %sv_and_or_or = ();
    my %pair;
    $set =~ s/^\.\/// if ($set =~ /^\.\//);
    my @pairs = split (/\|/, $set);
    my $count = 0;
    map{$pair{$_} = 1} @pairs;
    foreach my $pair (@pairs){
	my ($id1, $id2) = ('', '');
	($id1, $id2) = split (/=/, $pair) if ($pair =~ /=/);
	$id1 = $pair if ($pair !~ /=/);
	my ($id1_base, $min_read1) = split (/:/, $id1);
	my ($id2_base, $min_read2) = split (/:/, $id2);
	foreach my $type (keys %{$call{$id1}}){
	    if (($id2 ne '') and (exists ${$call{$id2}}{$type})){
		foreach my $chr (keys %{${$call{$id1}}{$type}}){
		    if (exists ${${$call{$id2}}{$type}}{$chr}){
			foreach my $pos1 (sort {$a <=> $b} keys %{${${$call{$id1}}{$type}}{$chr}}){
			    my $len1 = ${${${$call{$id1}}{$type}}{$chr}}{$pos1};
			    my $end1 = $pos1 + $len1 - 1;
			    $end1 = $pos1 if ($type eq 'INS');
			    my @match = ();
			    foreach my $pos2 (sort {$a <=> $b} keys %{${${$call{$id2}}{$type}}{$chr}}){
				last if ($pos2 > $pos1 + $ins_sd);
				my $len2 = ${${${$call{$id2}}{$type}}{$chr}}{$pos2};
				my $end2 = $pos2 + $len2 - 1;
				$end2 = $pos2 if ($type eq 'INS');
				next if ($end2 + $ins_sd < $pos1);
				my $ovlrate = 0;
				if ($type eq 'INS'){
				    if (abs ($pos1 - $pos2) <= $ins_sd){
					$ovlrate = 1 / abs ($pos1 - $pos2) if (abs ($pos1 - $pos2) > 0);
					$ovlrate = 1 if (abs ($pos1 - $pos2) == 0);
				    }
				}
				else{
				    next if ($len1 == 0) or ($len2 == 0);
				    if ((abs ($pos1 - $pos2) <= $var_sd) and (abs ($end1 - $end2) <= $var_sd)){
					if ($len1 >= $len2){
					    $ovlrate = $len2 / $len1;
					}
					else{
					    $ovlrate = $len1 / $len2;
					}
				    }
				    if (($pos1 >= $pos2) and ($pos1 <= $end2)){
					if ($end1 <= $end2){
					    if ($len1 >= $len2 * $min_overlap_ratio){
						$ovlrate = $len1 / $len2;
					    }
					}
					else{
					    my $overlen = $end2 - $pos1 + 1;
					    if (($overlen >= $len2 * $min_overlap_ratio) or ($overlen >= $len1 * $min_overlap_ratio)){
						my $ovlrate1 = $overlen / $len1;
						my $ovlrate2 = $overlen / $len2;
						if (($overlen >= $len2 * $min_overlap_ratio) and ($overlen >= $len1 * $min_overlap_ratio)){
						    $ovlrate = ($ovlrate1 + $ovlrate2) / 2;
						}
						elsif ((($len1 >= $len2) and ($len1 < $len2 * 5)) or (($len1 < $len2) and ($len2 < $len1 * 5))){
						    $ovlrate = ($ovlrate1 + $ovlrate2) / 2;
						}
					    }
					}
				    }
				    elsif (($pos2 >= $pos1) and ($pos2 <= $end1)){
					if ($end2 <= $end1){
					    if ($len2 >= $len1 * $min_overlap_ratio){
						$ovlrate = $len2 / $len1;
					    }
					}
					else{
					    my $overlen = $end1 - $pos2 + 1;
					    if (($overlen >= $len2 * $min_overlap_ratio) or ($overlen >= $len1 * $min_overlap_ratio)){
						my $ovlrate1 = $overlen / $len1;
						my $ovlrate2 = $overlen / $len2;
						if (($overlen >= $len2 * $min_overlap_ratio) and ($overlen >= $len1 * $min_overlap_ratio)){
						    $ovlrate = ($ovlrate1 + $ovlrate2) / 2;
						}
						elsif ((($len1 >= $len2) and ($len1 < $len2 * 5)) or (($len1 < $len2) and ($len2 < $len1 * 5))){
						    $ovlrate = ($ovlrate1 + $ovlrate2) / 2;
						}
					    }
					}
				    }
				}
				if ($ovlrate > 0){
				    push @match, "$pos2=$len2=$ovlrate";
				}
			    }
			    if (@match > 0){
				if (@match == 1){
				    my ($mpos, $mlen) = split (/=/, $match[0]);
				    ${$call_share{$pair}}{$type} .= "$chr=$pos1=$mpos=$len1=$mlen||";
				}
				else{
				    my $max_rate = 0;
				    my $best_pos = 0;
				    my $best_len = 0;
				    foreach (@match){
					my @info = split (/=/, $_);
					my $rate = $info[2];
					if ($rate > $max_rate){
					    $max_rate = $rate;
					    $best_pos = $info[0];
					    $best_len = $info[1];
					}
				    }
				    ${$call_share{$pair}}{$type} .= "$chr=$pos1=$best_pos=$len1=$best_len||";
				}
			    }
			}
		    }
		}
		if (!exists ${$sv_and{$type}}{$pair}){
		    foreach my $type (sort keys %{$call_share{$pair}}){
			my $call_share = 0;
			my @call_share = ();
			${$call_share{$pair}}{$type} =~ s/\|\|$// if (${$call_share{$pair}}{$type} =~ /\|\|$/);
			@call_share = split (/\|\|/, ${$call_share{$pair}}{$type});
                        foreach my $info (@call_share){
                            my ($chr, $pos1, $pos2, $len1, $len2) = split (/=/, $info);
                            if ($type eq 'INS'){
                                $len1 = 0;
                                $len2 = 0;
                            }
                            my @pos = ($pos1, $pos2);
                            my @len = ($len1, $len2);
                            my $pos = int (($pos1 + $pos2) / 2 + 0.5);
                            ${${${$sv_and{$type}}{$pair}}{$chr}}{$pos} = "$pos1=$pos2|$len1=$len2";
                        }
		    }
		}
	    }
	}
    }

    my $first_pair = $pairs[0];
    my $second_pair = '';
    $second_pair = $pairs[1] if (@pairs >= 2);
    my ($id1, $id2) = ('', '');
    my ($id3, $id4) = ('', '');
    ($id1, $id2) = split (/=/, $first_pair) if ($first_pair =~ /=/);
    $id1 = $first_pair if ($first_pair !~ /=/);
    ($id3, $id4) = split (/=/, $second_pair) if ($second_pair =~ /=/);
    $id3 = $second_pair if ($second_pair !~ /=/);
    my ($id1_base) = split (/:/, $id1);
    my ($id2_base) = split (/:/, $id2) if ($id2 ne '');
    my ($id3_base) = split (/:/, $id3) if ($id3 ne '');
    my ($id4_base) = split (/:/, $id4) if ($id4 ne '');
    my @ids;
    push @ids, $id1_base, $id3_base if ($id3 ne '');
    push @ids, $id2_base if ($id2 ne '');
    push @ids, $id4_base if ($id4 ne '');
    if (@pairs == 1){
	foreach my $type (keys %sv_and){
	    foreach my $pair (keys %{$sv_and{$type}}){
		next if ($pair ne $pairs[0]);
		my ($id1, $id2) = split (/=/, $pair);
		my $id1_base = $1 if ($id1 =~ /(.+):\d+/);
		my $id2_base = $1 if ($id2 =~ /(.+):\d+/);
		my $id1_lc = lc $id1_base;
		my $id2_lc = lc $id2_base;
		foreach my $chr (keys %{${$sv_and{$type}}{$pair}}){
		    foreach my $pos (sort {$a <=> $b} keys %{${${$sv_and{$type}}{$pair}}{$chr}}){
			my ($pos_g, $len_g) = split (/\|/, ${${${$sv_and{$type}}{$pair}}{$chr}}{$pos});
			my %tool_pos;
			my %tool_read;
			my @pos_g = split (/=/, $pos_g);
			my @len_g = split (/=/, $len_g);
			my $av_len = int (($len_g[0] + $len_g[1]) / 2 + 0.5);
			my $class = 'S';
			if ($type eq 'DUP'){
			    $class = 'S' if ($av_len < 1000);
			    $class = 'M' if ($av_len >= 1000);
			    $class = 'L' if ($av_len >= 100000);
			}
			else{
			    $class = 'S' if ($av_len < 1000);
			    $class = 'M' if ($av_len >= 1000);
			}
			my $av_pos = int (($pos_g[0] + $pos_g[1]) / 2 + 0.5);			
			my $read_1 = ${${${$sv_read{$id1}}{$type}}{$chr}}{$pos_g[0]};
			my $read_2 = ${${${$sv_read{$id2}}{$type}}{$chr}}{$pos_g[1]};
			my $ave_read = int ($read_1 + $read_2 / 2 + 0.5);
			my $info = "$id1_base-$pos_g[0],$id2_base-$pos_g[1]";
			${$vcf{$chr}}{$av_pos} = "$chr\t$av_pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$av_len;READS=$ave_read;POS=$pos_g\tTOOLS=$info";
		    }
		}
	    }
	}
    }
    
    my %match_ref;
    my %pre_info;
    my %overlap;
    $overlap{'DEL'} = 0;
    $overlap{'DUP'} = 0;
    $overlap{'INV'} = 0;
    $overlap{'INS'} = 0;
    
    my %call_num;
    my %match_num;
    my %uncall_num;
    my %recal_num;
    
    my $under_del_ref = 0;
    my $over_del_ref = 0;
    my $sum_under_ref = 0;
    
    my $ID1 = '';
    my $ID2 = '';
    my %ID1_MIE_TP;
    my %ID2_MIE_TP;
    my %ID1_MIE_TP2;
    my %ID2_MIE_TP2;
    my %ID1_MIE_recal;
    my %ID2_MIE_recal;
    my $ID1_parent_flag = 1;
    my $ID2_parent_flag = 1;
    if ($set2 !~ /=/){
        $ID1 = $1 if ($set2 =~ /(.+):\d+$/);
    }
    else{
        ($ID1, $ID2) = split (/=/, $set2);
        $ID1 = $1 if ($ID1 =~ /(.+):\d+$/);
        $ID2 = $1 if ($ID2 =~ /(.+):\d+$/);
    }
    if ((!exists $parent1_file{$ID1}) and (!exists $parent2_file{$ID1})){
#        print STDERR "$ID1 parent vcf files not found: \n";
        $ID1_parent_flag = 0;
    }
    if (($ID2 ne '') and (!exists $parent1_file{$ID2}) and (!exists $parent2_file{$ID2})){
#        print STDERR "$ID2 parent vcf files not found: \n";
        $ID2_parent_flag = 0;
    }
    foreach my $chr (sort keys %vcf){
	foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
	    my $len = 0;
	    $len = $1 if (${$vcf{$chr}}{$pos} =~ /SVLEN=(\d+)/);
	    my $type = $1 if (${$vcf{$chr}}{$pos} =~ /SVTYPE=(.+?);/);
	    my $end = $pos + $len;
	    my ($pre_chr, $pre_pos, $pre_size, $pre_hit) = ('', 0, 0, 0);
	    ($pre_chr, $pre_pos, $pre_size, $pre_hit) = split (/=/, $pre_info{$type}) if (exists $pre_info{$type});
	    if (($type eq 'INS') and ($pre_chr eq $chr) and ($pos - $pre_pos < $var_sd)){
		if ((($len > 0) and ($pre_size > $var_sd * 0.5) and ($len > $var_sd * 0.5)) or ($len == 0)){
		    $overlap{$type} ++;
		    if ($pre_hit == 1){
			$pre_info{$type} = "$chr=$pos=$len=1" if ($pos + $len > $pre_pos + $pre_size);
			next;
		    }
		}
	    }
	    elsif (($type ne 'INS') and ($pre_chr eq $chr) and ($pre_pos + $pre_size - $pos > $len * 0.8) and ($pre_pos + $pre_size - $pos > $pre_size * 0.8) and ($pre_size > $var_sd * 0.5) and ($len > $var_sd * 0.5)){
		$overlap{$type} ++;
		if ($pre_hit == 1){
		    $pre_info{$type} = "$chr=$pos=$len=1" if ($pos + $len > $pre_pos + $pre_size);
		    next;
		}
	    }
	    else{
                $call_num{$type} ++;
	    }
            my $trio_flag1 = 0;
            my $trio_flag2 = 0;
            if (exists ${$parent1{$ID1}}{$type}){
                foreach my $fpos (sort {$a <=> $b} keys %{${${$parent1{$ID1}}{$type}}{$chr}}){
                    last if ($fpos > $end + $ins_sd);
                    my $flen = ${${${$parent1{$ID1}}{$type}}{$chr}}{$fpos};
                    my $fend = $fpos + $flen - 1;
                    next if ($fend < $pos - $ins_sd);
                    if ($type eq 'INS'){
                        if (abs ($pos - $fpos) <= $ins_sd){
                            $trio_flag1 = 1;
                            last;
                        }
                    }
                    else{
                        if (($pos <= $fpos) and ($end >= $fend)){
                            if ($flen >= $len * $min_overlap_ratio3){
                                $trio_flag1 = 1;
                            }
                        }
                        elsif (($pos >= $fpos) and ($end <= $fend)){
                            if ($len >= $flen * $min_overlap_ratio3){
                                $trio_flag1 = 1;
                            }
                        }
                        elsif (($pos >= $fpos) and ($pos <= $fend)){
                            my $overlap = $fend - $pos + 1;
                            if (($overlap >= $len * $min_overlap_ratio3) or ($overlap >= $flen * $min_overlap_ratio3)){
                                $trio_flag1 = 1;
                            }
                        }
                        elsif (($end >= $fpos) and ($end <= $fend)){
                            my $overlap = $end - $fpos + 1;
                            if (($overlap >= $len * $min_overlap_ratio3) or ($overlap >= $flen * $min_overlap_ratio3)){
                                $trio_flag1 = 1;
                            }
                        }
                    }
                    last if ($trio_flag1 == 1);
                }
                if ($trio_flag1 == 0){
                    foreach my $fpos (sort {$a <=> $b} keys %{${${$parent2{$ID1}}{$type}}{$chr}}){
                        last if ($fpos > $end + $ins_sd);
                        my $flen = ${${${$parent2{$ID1}}{$type}}{$chr}}{$fpos};
                        my $fend = $fpos + $flen - 1;
                        next if ($fend < $pos - $ins_sd);
                        if ($type eq 'INS'){
                            if (abs ($pos - $fpos) <= $ins_sd){
                                $trio_flag1 = 1;
                                last;
                            }
                        }
                        else{
                            if (($pos <= $fpos) and ($end >= $fend)){
                                if ($flen >= $len * $min_overlap_ratio3){
                                    $trio_flag1 = 1;
                                }
                            }
                            elsif (($pos >= $fpos) and ($end <= $fend)){
                                if ($len >= $flen * $min_overlap_ratio3){
                                    $trio_flag1 = 1;
                                }
                            }
                            elsif (($pos >= $fpos) and ($pos <= $fend)){
                                my $overlap = $fend - $pos + 1;
                                if (($overlap >= $len * $min_overlap_ratio3) or ($overlap >= $flen * $min_overlap_ratio3)){
                                    $trio_flag1 = 1;
                                }
                            }
                            elsif (($end >= $fpos) and ($end <= $fend)){
                                my $overlap = $end - $fpos + 1;
                                if (($overlap >= $len * $min_overlap_ratio3) or ($overlap >= $flen * $min_overlap_ratio3)){
                                    $trio_flag1 = 1;
                                }
                            }
                        }
                        last if ($trio_flag1 == 1);
                    }
                }
            }
            if ($ID2 ne ''){
                if (exists ${$parent1{$ID2}}{$type}){
                    foreach my $fpos (sort {$a <=> $b} keys %{${${$parent1{$ID2}}{$type}}{$chr}}){
                        last if ($fpos > $end + $ins_sd);
                        my $flen = ${${${$parent1{$ID2}}{$type}}{$chr}}{$fpos};
                        my $fend = $fpos + $flen - 1;
                        next if ($fend < $pos - $ins_sd);
                        if ($type eq 'INS'){
                            if (abs ($pos - $fpos) <= $ins_sd){
                                $trio_flag2 = 1;
                                last;
                            }
                        }
                        else{
                            if (($pos <= $fpos) and ($end >= $fend)){
                                if ($flen >= $len * $min_overlap_ratio3){
                                    $trio_flag2 = 1;
                                }
                            }
                            elsif (($pos >= $fpos) and ($end <= $fend)){
                                if ($len >= $flen * $min_overlap_ratio3){
                                    $trio_flag2 = 1;
                                }
                            }
                            elsif (($pos >= $fpos) and ($pos <= $fend)){
                                my $overlap = $fend - $pos + 1;
                                if (($overlap >= $len * $min_overlap_ratio3) or ($overlap >= $flen * $min_overlap_ratio3)){
                                    $trio_flag2 = 1;
                                }
                            }
                            elsif (($end >= $fpos) and ($end <= $fend)){
                                my $overlap = $end - $fpos + 1;
                                if (($overlap >= $len * $min_overlap_ratio3) or ($overlap >= $flen * $min_overlap_ratio3)){
                                    $trio_flag2 = 1;
                                }
                            }
                        }
                        last if ($trio_flag2 == 1);
                    }
                    if ($trio_flag2 == 0){
                        foreach my $fpos (sort {$a <=> $b} keys %{${${$parent2{$ID2}}{$type}}{$chr}}){
                            last if ($fpos > $end + $ins_sd);
                            my $flen = ${${${$parent2{$ID2}}{$type}}{$chr}}{$fpos};
                            my $fend = $fpos + $flen - 1;
                            next if ($fend < $pos - $ins_sd);
                            if ($type eq 'INS'){
                                if (abs ($pos - $fpos) <= $ins_sd){
                                    $trio_flag2 = 1;
                                    last;
                                }
                            }
                            else{
                                if (($pos <= $fpos) and ($end >= $fend)){
                                    if ($flen >= $len * $min_overlap_ratio3){
                                        $trio_flag2 = 1;
                                    }
                                }
                                elsif (($pos >= $fpos) and ($end <= $fend)){
                                    if ($len >= $flen * $min_overlap_ratio3){
                                        $trio_flag2 = 1;
                                    }
                                }
                                elsif (($pos >= $fpos) and ($pos <= $fend)){
                                    my $overlap = $fend - $pos + 1;
                                    if (($overlap >= $len * $min_overlap_ratio3) or ($overlap >= $flen * $min_overlap_ratio3)){
                                        $trio_flag2 = 1;
                                    }
                                }
                                elsif (($end >= $fpos) and ($end <= $fend)){
                                    my $overlap = $end - $fpos + 1;
                                    if (($overlap >= $len * $min_overlap_ratio3) or ($overlap >= $flen * $min_overlap_ratio3)){
                                        $trio_flag2 = 1;
                                    }
                                }
                            }
                            last if ($trio_flag2 == 1);
                        }
                    }
                }
            }

	    my $flag = 0;
	    my $hit_len = 0;
	    if ($type eq 'DEL'){
		my $dellen = 0;
		my $hit_bp = 0;
		foreach my $bp (sort {$a <=> $b} keys %{${$ref{'DEL'}}{$chr}}){
		    next if (exists ${${$match_ref{$type}}{$chr}}{$bp});
		    $dellen = ${${$ref{'DEL'}}{$chr}}{$bp};
		    my $delend = $bp + $dellen - 1;
		    $hit_bp = $bp;
		    $hit_len = $dellen;
		    if ((abs ($pos - $bp) <= $var_sd) and (abs ($end - $delend) <= $var_sd)){
			$flag = 1;
		    }
		    if (($pos >= $bp) and ($pos <= $delend)){
			if ($end <= $delend){
			    if ($len >= $dellen * $min_overlap_ratio2){
				$flag = 1;
			    }
			}
			else{
			    if (($delend - $pos >= $dellen * $min_overlap_ratio2) and ($delend - $pos >= $len * $min_overlap_ratio2)){
				$flag = 1;
			    }
			}
		    }
		    elsif (($bp >= $pos) and ($bp <= $end)){
			if ($delend <= $end){
			    if ($dellen >= $len * $min_overlap_ratio2){
				$flag = 1;
			    }
			}
			else{
			    if (($end - $bp >= $dellen * $min_overlap_ratio2) and ($end - $bp >= $len * $min_overlap_ratio2)){
				$flag = 1;
			    }
			}
		    }
                    if ($flag == 1){
                        ${${$match_ref{$type}}{$chr}}{$bp} = 1;
                        $hit_len = $dellen;
			last;
                    }
		}
	    }
	    elsif ($type eq 'DUP'){
		my $duplen = 0;
		my $hit_bp = 0;
		foreach my $bp (sort {$a <=> $b} keys %{${$ref{'DUP'}}{$chr}}){
		    next if (exists ${${$match_ref{$type}}{$chr}}{$bp});
		    $duplen = ${${$ref{'DUP'}}{$chr}}{$bp};
		    my $dupend = $bp + $duplen - 1;
		    $hit_bp = $bp;
		    $hit_len = $duplen;
		    if ((abs ($pos - $bp) <= $var_sd) and (abs ($end - $dupend) <= $var_sd)){
			$flag = 1;
			${${$match_ref{$type}}{$chr}}{$bp} = 1;
			last;
		    }
		    if (($pos >= $bp) and ($pos <= $dupend)){
			if ($end <= $dupend){
			    if ($len >= $duplen * $min_overlap_ratio2){
				$flag = 1;
			    }
			}
			else{
			    if (($dupend - $pos >= $duplen * $min_overlap_ratio2) and ($dupend - $pos >= $len * $min_overlap_ratio2)){
				$flag = 1;
			    }
			}
		    }
		    elsif (($bp >= $pos) and ($bp <= $end)){
			if ($dupend <= $end){
			    if ($duplen >= $len * $min_overlap_ratio2){
				$flag = 1;
			    }
			}
			else{
			    if (($end - $bp >= $duplen * $min_overlap_ratio2) and ($end - $bp >= $len * $min_overlap_ratio2)){
				$flag = 1;
			    }
			}
		    }
                    if ($flag == 1){
                        ${${$match_ref{$type}}{$chr}}{$bp} = 1;
                        $hit_len = $duplen;
			last;
                    }
		}
		if ($flag == 0){
                    foreach my $bp (sort {$a <=> $b} keys %{$ref_ins{$chr}}){
                        last if ($bp > $end + $len);
                        next if (exists ${${$match_ref{'INS'}}{$chr}}{$bp});
                        my $inslen = ${$ref_ins{$chr}}{$bp};
                        next if ($inslen < $min_sv_len);
                        my $insend = $bp + $inslen - 1;
                        if (($pos <= $bp) and ($end >= $insend)){
                            if ($inslen >= $len * $min_overlap_ratio2){
                                $flag = 2;
                            }
                        }
                        elsif (($pos >= $bp) and ($end <= $insend)){
                            if ($len >= $inslen * $min_overlap_ratio2){
                                $flag = 2;
                            }
                        }
                        elsif (($pos >= $bp) and ($pos <= $insend)){
                            my $overlap = $insend - $pos + 1;
                            if (($overlap >= $len * $min_overlap_ratio2) and ($overlap >= $inslen * $min_overlap_ratio2)){
                                $flag = 2;
                            }
                        }
                        elsif (($end >= $bp) and ($end <= $insend)){
                            my $overlap = $end - $bp + 1;
                            if (($overlap >= $len * $min_overlap_ratio2) and ($overlap >= $inslen * $min_overlap_ratio2)){
                                $flag = 2;
                            }
                        }
                        if ($flag == 2){
                            ${${$match_ref{'INS'}}{$chr}}{$bp} = 1;
                            my $bp1 = $bp - $inslen;
                            my $bp2 = $bp + $inslen;
                            $hit_len = $inslen;
                            ${${$match_ref{'INS'}}{$chr}}{$bp1} = 1 if (exists ${$ref_ins{$chr}}{$bp1});
                            ${${$match_ref{'INS'}}{$chr}}{$bp2} = 1 if (exists ${$ref_ins{$chr}}{$bp2});
                            last;
                        }
                    }
                }
	    }
	    elsif ($type eq 'INS'){
		my $duplen = 0;
		my $hit_bp = 0;
		foreach my $bp (sort {$a <=> $b} keys %{${$ref{'INS'}}{$chr}}){
		    next if (exists ${${$match_ref{$type}}{$chr}}{$bp});
		    $duplen = ${${$ref{'INS'}}{$chr}}{$bp};
		    $hit_bp = $bp;
		    my $dupend = $bp + $duplen - 1;
		    if (abs ($pos - $bp) <= $ins_sd){
			$flag = 1;
			${${$match_ref{$type}}{$chr}}{$bp} = 1;
		    }
		    last if ($flag > 0);
		}
                if ($flag == 0){
                    foreach my $bp (sort {$a <=> $b} keys %{$ref_dup{$chr}}){
                        last if ($bp > $end + $ins_sd);
                        next if (exists ${${$match_ref{'DUP'}}{$chr}}{$bp});
                        $duplen = ${$ref_dup{$chr}}{$bp};
                        my $dupend = $bp + $duplen - 1;
                        next if ($dupend < $pos - $ins_sd);
                        $hit_bp = $bp;
                        if ((abs ($pos - $bp) <= $ins_sd) or (abs ($pos - $dupend) <= $ins_sd)){
                            $flag = 2;
                            ${${$match_ref{'DUP'}}{$chr}}{$bp} = 1;
                            last;
                        }
                    }
                }
	    }
	    elsif ($type eq 'INV'){
		my $invlen = 0;
		my $hit_bp = 0;
		foreach my $bp (sort {$a <=> $b} keys %{${$ref{'INV'}}{$chr}}){
		    next if (exists ${${$match_ref{$type}}{$chr}}{$bp});
		    $invlen = ${${$ref{'INV'}}{$chr}}{$bp};
		    my $invend = $bp + $invlen - 1;
		    $hit_bp = $bp;
		    $hit_len = $invlen;
		    if ((abs ($pos - $bp) <= $var_sd) and (abs ($end - $invend) <= $var_sd)){
			$flag = 1;
		    }
		    if (($pos >= $bp) and ($pos <= $invend)){
			if ($end <= $invend){
			    if ($len >= $invlen * $min_overlap_ratio2){
				$flag = 1;
			    }
			}
			else{
			    if (($invend - $pos >= $invlen * $min_overlap_ratio2) and ($invend - $pos >= $len * $min_overlap_ratio2)){
				$flag = 1;
			    }
			}
		    }
		    elsif (($bp >= $pos) and ($bp <= $end)){
			if ($invend <= $end){
			    if ($invlen >= $len * $min_overlap_ratio){
				$flag = 1;
			    }
			}
			else{
			    if (($end - $bp >= $invlen * $min_overlap_ratio2) and ($end - $bp >= $len * $min_overlap_ratio2)){
				$flag = 1;
			    }
			}
		    }
                    if ($flag == 1){
                        ${${$match_ref{$type}}{$chr}}{$bp} = 1;
			last;
                    }
		}
	    }
	    if ($flag >= 1){
                $match_num{$type} ++;
                $ID1_MIE_TP{$type} ++ if ($ID1_parent_flag == 1) and ($trio_flag1 == 0);
                $ID2_MIE_TP{$type} ++ if ($ID2_parent_flag == 1) and ($ID2 ne '') and ($trio_flag2 == 0);
                my $recal_flag = 0;
		if ($type eq 'DEL'){
		    if ($sv_size eq 'SS'){
			if ($hit_len > 100){
                            $recal_flag = 1;
			}
		    }
		    elsif ($sv_size eq 'S'){
			if (($hit_len <= 100) or ($hit_len > 1000)){
                            $recal_flag = 1;
			}
		    }
		    elsif ($sv_size eq 'M'){
			if (($hit_len <= 1000) or ($hit_len > 100000)){
                            $recal_flag = 1;
			}
		    }
		    elsif ($sv_size eq 'L'){
			if ($hit_len <= 100000){
                            $recal_flag = 1;
			}
		    }
		}
		elsif ($type eq 'DUP'){
		    if ($sv_size eq 'S'){
			if ($hit_len > 1000){
                            $recal_flag = 1;
			}
		    }
		    elsif ($sv_size eq 'M'){
			if (($hit_len <= 1000) or ($hit_len > 100000)){
                            $recal_flag = 1;
			}
		    }
		    elsif ($sv_size eq 'L'){
			if ($hit_len <= 100000){
                            $recal_flag = 1;
			}
		    }
		}
                if ($recal_flag == 1){
                    $recal_num{$type} ++;
                    $ID1_MIE_recal{$type} ++ if ($ID1_parent_flag == 1) and ($trio_flag1 == 1);
                    $ID2_MIE_recal{$type} ++ if ($ID2_parent_flag == 1) and ($ID2 ne '') and ($trio_flag2 == 1);
                }
                if (($output_vcf == 1) and ($set2 =~ /=/)){
		    if ($hg19_flag == 0){
			print OUT2 "TP ${$vcf{$chr}}{$pos}\n";
		    }
		    else{
			my @line2 = split (/\t/, ${$vcf{$chr}}{$pos});
			$line2[0] = 'chr' . $line2[0];
			my $new_line2 = join ("\t", @line2);
			print OUT2 "TP $new_line2\n";
		    }
                }
		$pre_info{$type} = "$chr=$pos=$len=1";
	    }
            else{
                if (($output_vcf == 1) and ($set2 =~ /=/)){
		    if ($hg19_flag == 0){
			print OUT2 "FP ${$vcf{$chr}}{$pos}\n";
		    }
		    else{
			my @line2 = split (/\t/, ${$vcf{$chr}}{$pos});
			$line2[0] = 'chr' . $line2[0];
			my $new_line2 = join ("\t", @line2);
			print OUT2 "FP $new_line2\n";
		    }
                }
		$pre_info{$type} = "$chr=$pos=$len=0";
	    }
	    if ($flag == 2){
                $uncall_num{$type} ++ if ($type =~ /DEL|DUP/);
                $ID1_MIE_TP2{$type} ++ if ($ID1_parent_flag == 1) and ($trio_flag1 == 0);
                $ID2_MIE_TP2{$type} ++ if ($ID2_parent_flag == 1) and ($ID2 ne '') and ($trio_flag2 == 0);
                $pre_info{$type} = "$chr=$pos=$len=1";
	    }
	}
    }
    
    print "# $set2\n";
    print OUT "# $set2\n" if ($vcf_out ne '');
    foreach my $type (sort keys %call_num){
        $ID1_MIE_TP{$type} = 0 if (!exists $ID1_MIE_TP{$type});
        $ID1_MIE_TP2{$type} = 0 if (!exists $ID1_MIE_TP2{$type});
        $ID1_MIE_recal{$type} = 0 if (!exists $ID1_MIE_recal{$type});
        $ID2_MIE_TP{$type} = 0 if (!exists $ID2_MIE_TP{$type});
        $ID2_MIE_TP2{$type} = 0 if (!exists $ID2_MIE_TP2{$type});
        $ID2_MIE_recal{$type} = 0 if (!exists $ID2_MIE_recal{$type});
        $recal_num{$type} = 0 if (!exists $recal_num{$type});
        $uncall_num{$type} = 0 if (!exists $uncall_num{$type});
        my $recall = 0;
        my $precis = 0;
        if (($ID1_parent_flag == 0) and ($ID2_parent_flag == 0)){
            $recall = int (($match_num{$type} - $recal_num{$type}) / $ref_num{$type} * 10000 + 0.5) / 100 if ($type eq 'DEL');
            $recall = int (($match_num{$type} - $recal_num{$type} + 0.5) / ($ref_num{$type} + $uncall_num{$type}) * 10000) / 100 if ($type eq 'DUP');
            $recall = int ($match_num{$type} / ($ref_num{$type} + $uncall_num{$type}) * 10000 + 0.5) / 100 if ($type eq 'INS');
            $recall = int ($match_num{$type} / ($ref_num{$type}) * 10000 + 0.5) / 100 if ($type eq 'INV');
            $precis = int ($match_num{$type} / $call_num{$type} * 10000 + 0.5) / 100;
        }
        elsif (($ID2 eq '') or ($ID2_parent_flag == 0)){
	    my $match_num = 0;
            $match_num = $match_num{$type} if (exists $match_num{$type});
            $match_num -= $ID1_MIE_TP{$type};
            $recall = int (($match_num - $ID1_MIE_recal{$type}) / $ref_num{$type} * 10000 + 0.5) / 100 if ($type eq 'DEL');
            $recall = int (($match_num - $ID1_MIE_recal{$type}) / ($ref_num{$type} + $ID1_MIE_TP2{$type}) * 10000 + 0.5) / 100 if ($type eq 'DUP');
            $recall = int ($match_num / ($ref_num{$type} + $ID1_MIE_TP2{$type}) * 10000 + 0.5) / 100 if ($type eq 'INS');
            $recall = int ($match_num / $ref_num{$type} * 10000 + 0.5) / 100 if ($type eq 'INV');
            $precis = int ($match_num / $call_num{$type} * 10000 + 0.5) / 100;
        }
        elsif ($ID1_parent_flag == 0){
	    my $match_num = 0;
            $match_num = $match_num{$type} if (exists $match_num{$type});
            $match_num -= $ID2_MIE_TP{$type};
            $recall = int (($match_num - $ID2_MIE_recal{$type}) / $ref_num{$type} * 10000 + 0.5) / 100 if ($type eq 'DEL');
            $recall = int (($match_num - $ID2_MIE_recal{$type}) / ($ref_num{$type} + $ID2_MIE_TP2{$type}) * 10000 + 0.5) / 100 if ($type eq 'DUP');
            $recall = int ($match_num / ($ref_num{$type} + $ID2_MIE_TP2{$type}) * 10000 + 0.5) / 100 if ($type eq 'INS');
            $recall = int ($match_num / $ref_num{$type} * 10000 + 0.5) / 100 if ($type eq 'INV');
            $precis = int ($match_num / $call_num{$type} * 10000 + 0.5) / 100;
        }
	else{
            my $recall_1 = 0;
            my $recall_2 = 0;
            my $precis_1 = 0;
            my $precis_2 = 0;
	    my $match_num = 0;
            $match_num = $match_num{$type} if (exists $match_num{$type});
            my $match_num_1 = $match_num - $ID1_MIE_TP{$type};
            my $match_num_2 = $match_num - $ID2_MIE_TP{$type};
            $precis_1 = int ($match_num_1 / $call_num{$type} * 10000 + 0.5) / 100;
            $precis_2 = int ($match_num_2 / $call_num{$type} * 10000 + 0.5) / 100;
            if ($type eq 'DEL'){
                $recall_1 = int (($match_num_1 - $ID1_MIE_recal{$type}) / $ref_num{$type} * 10000 + 0.5) / 100;
                $recall_2 = int (($match_num_2 - $ID2_MIE_recal{$type}) / $ref_num{$type} * 10000 + 0.5) / 100;
            }
            elsif ($type eq 'DUP'){
                $recall_1 = int (($match_num_1 - $ID1_MIE_recal{$type}) / ($ref_num{$type} + $ID1_MIE_TP2{$type}) * 10000 + 0.5) / 100;
                $recall_2 = int (($match_num_2 - $ID2_MIE_recal{$type}) / ($ref_num{$type} + $ID2_MIE_TP2{$type}) * 10000 + 0.5) / 100;
            }
            elsif ($type eq 'INS'){
                $recall_1 = int ($match_num_1 / ($ref_num{$type} + $ID1_MIE_TP2{$type}) * 10000 + 0.5) / 100;
                $recall_2 = int ($match_num_2 / ($ref_num{$type} + $ID2_MIE_TP2{$type}) * 10000 + 0.5) / 100;
            }
            elsif ($type eq 'INV'){
                $recall_1 = int ($match_num_1 / $ref_num{$type} * 10000 + 0.5) / 100;
                $recall_2 = int ($match_num_2 / $ref_num{$type} * 10000 + 0.5) / 100;
            }
            
            $recall = int (($recall_1 + $recall_2) / 2 * 100 + 0.5) / 100;
            $precis = int (($precis_1 + $precis_2) / 2 * 100 + 0.5) / 100;
        }
        $recall = int ($recall * 10 + 0.5) / 10 if ($recall >= 0.1);
        $precis = int ($precis * 10 + 0.5) / 10 if ($precis >= 0.1);
        my $overlap_num = $call_num{$type};
        if ($vcf_out eq ''){
	    print "$type ($sv_size)\t$recall\t$precis\t$overlap_num\n";
	}
	else{
	    print OUT "$type ($sv_size)\t$recall\t$precis\t$overlap_num\n" if ($vcf_out ne '');
	}
    }
}
close (OUT) if ($vcf_out ne '');
print "\n";
