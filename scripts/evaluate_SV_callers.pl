#!/usr/bin/perl -w
use strict;
use File::Basename;
use FindBin qw/$Bin $Script/;
use Getopt::Long;
use Pod::Usage;


my $script_path = $Bin;

my $ref_sv_simA = "$script_path/../Ref_SV/Sim-A.SV.vcf";

my $ref_sv_NA12878 = "$script_path/../Ref_SV/NA12878_DGV-2016_LR-assembly.vcf";

my $ref_mei = "$script_path/../Simulated_data_in_our_study/Sim-MEI/Sim-MEI.chr17.vcf";
my $ref_vei = "$script_path/../Simulated_data_in_our_study/Sim-VEI/Sim-VEI.chr17.vcf";
my $ref_numt = "$script_path/../Simulated_data_in_our_study/Sim-NUMT/Sim-NUMT.chr17.vcf";

my $gap_bed = "$script_path/../Ref_SV/gap.bed";

my $region_bed = '';

my $var_parent1 = '';
my $var_parent2 = '';

my $min_reads = 0;

my $read_set = 1;

#my @min_reads = (2, 3, 4, 5, 6, 7, 8, 9, 10);
my @min_reads = (2, 3, 4, 5, 6, 7, 8, 9, 10, 12);

my @size = ('A', 'S', 'M', 'L');

my $ref_sv = 'N';

my $ref_sv2 = '';

my $sv_type = 'ALL';

my $target_chr = 'all';

my $include_y = 0;

my $min_sv_len = 50;

my $max_sv_len = 2000000;

my $min_ref_len = 30;

my $max_ref_len = 2000000;

my $var_sd = 125;
#my $ins_sd = 1000;
my $ins_sd = 200;
my $ctx_sd = 1000;
my $mei_sd = 125;

my $min_ctx = 20000;

my $min_overlap_ratio = 0.5;
my $min_overlap = 0.5;

my $min_qual = 90;
my $min_del_qual = 0;

my $eval_gt = 0;
my $eval_bp = 0;

my $ins_eval = 0;

my $min_reads_parent = 2;
my %min_reads_parent;

my $substract_ME = 0;

my $out_TF = 0;

my $help;

GetOptions(
    'ref|r=s' => \$ref_sv,
    'ref2|r2=s' => \$ref_sv2,
    'parent1|p1=s' => \$var_parent1,
    'parent2|p2=s' => \$var_parent2,
    'sv_type|st=s' => \$sv_type,
    'region_bed|rb=s' => \$region_bed,
    'chr|c=s' => \$target_chr,
    'len|l=i' => \$min_sv_len,
    'xlen|xl=i' => \$max_sv_len,
    'ref_len|rl=i' => \$min_ref_len,
    'ref_xlen|rxl=i' => \$max_ref_len,
    'min_read|mr=i' => \$min_reads,
    'read_set|rs=i' => \$read_set,
    'min_ovl|mo=f' => \$min_overlap,
    'min_ins|mins=i' => \$ins_sd,
    'eval_gt|eg' => \$eval_gt,
    'eval_bp|eb' => \$eval_bp,
    'ins|i' => \$ins_eval,
    'in_y|y' => \$include_y,
    'out_tf|of=i' => \$out_TF,
    'sub_me|sm' => \$substract_ME,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

my $var_child = shift @ARGV;

my $var_base = basename ($var_child);
$var_base = $1 if ($var_base =~ /(.+)\./);

=head1 SYNOPSIS

  evaluate_SV_callers.pl <option> [vcf file]
  evaluate_SV_callers.pl <option> -p1 <parent1 vcf file> -p2 <parent2 vcf file> [vcf file] (for integrating trio data)

  Options:
   --ref or -r <STR>        reference SV type (A|N), A: Sim-A, N: NA12878 [default: N]
   --sv_type or -st <STR>   SV type (ALL|DEL|DUP|INS|INV|TRA) if specifying multi type, e.g., DEL,INS [default: ALL]
   --parent1 or -p1 <STR>   SV vcf file of parent1 (optional)
   --parent2 or -p2 <STR>   SV vcf file of parent2 (optional)
   --ref2 or -r2 <STR>	    reference SV vcf file (optional)
   --ref_dir or -rd <STR>   directory containing reference SV vcf files (optional)
   --chr or -c <STR>        target chromosome to be analyzed [all or chr name(s), e.g., 4,5,6,X; default: all]
   --region_bed or -rb <STR> bed file indicating the regions to be analyzed (optional)
   --len or -l <INT>        minimum size (bp) of SV [default: 50]
   --xlen or -xl <INT>      maximum size of SV [default: 2000000]
   --ref_len or -rl <INT>   minimum size of reference SV [default: 30]
   --ref_xlen or -rxl <INT> maximum size of reference SV [default: 2000000]
   --min_read or -mr <INT>  minimum number of reads supporting an SV [default: 0]
   --read_set or -rs <INT>  read set of minimum number of reads (1: 2,3,4,5,6,7,8,9,10,12; 2: 2,3,5,7,9,11,13,15,17,19,30,40) [default: 1]
   --min_ovl or -mo <FLOAT> minimum rate of reciprocal overlap between called non-INS-SVs and reference non-INS-SVs [default: 0.5 for NA12878, 0.8 for Sim-A > 1Kb SVs, 0.6 for Sim-A <= 1 Kb SVs]
   --min_ins or -mins <INT> maximum allowable length between the breakpoints of called INSs and reference INSs [default: 200]
   --eval_gt or -eg         determine genotype accuracy for each SV type with the Sim-A data [default: false]
   --eval_bp or -eb         determine breakpoint and SV size accuracy for each SV type with the Sim-A data [default: false]
   --ins or -i              evaluate accuracy for INSs with called insertion size [default: false]
   --in_y or -y             include chrY [default: false]
   --out_tf or -of <INT>    output TP, FO, or TP/FP calls of vcf file (0: not output, 1: output TP calls, 2: output FP calls, 3: output TP and FP calls with labels [default: 0]
   --sub_me or -sm          substract Mendelian inheritance errors from true positive calls to calculate precision and recall [default: false]
   --help or -h             output help message
   
=cut

if ((($eval_gt == 1) or ($eval_bp == 1)) and ($ref_sv ne 'A')){
    die "eval_gt or eval_bp is only effective when ref is \'A\': \n";
}

my $parent_flag = 0;
$parent_flag = 1 if ($var_parent1 ne '') or ($var_parent2 ne '');

@min_reads = (2, 3, 5, 7, 9, 11, 13, 15, 17, 19, 30, 40) if ($read_set == 2);
@min_reads = (2, 3, 5, 7, 9, 11, 13, 15, 17, 19) if ($read_set == 3);

my $child_vcf_base = basename ($var_child);
if ($child_vcf_base =~ /PBHoney|Sniffles/){
    $min_reads_parent = 2;
}
elsif ($child_vcf_base =~ /RAPTR/){
    $min_reads_parent = 0;
    $min_reads_parent{'DEL'} = 8;
    $min_reads_parent{'DUP'} = 4;
    $min_reads_parent{'INS'} = 4;
}
elsif ($child_vcf_base =~ /SVfinder/){
    $min_reads_parent = 0;
    $min_reads_parent{'DEL'} = 5;
    $min_reads_parent{'INS'} = 4;
    $min_reads_parent{'INV'} = 5;
}

if ($min_overlap == 0){
    $min_overlap_ratio = 0.5 if (uc $ref_sv eq 'N');
}
else{
    $min_overlap_ratio = $min_overlap;
}

@min_reads = (2, 3, 5, 7, 9, 11, 13, 15, 17, 19, 30, 40) if ($read_set == 2);

my $ref_sv_tag = $ref_sv;


if (($ref_sv eq 'A') or ($ref_sv eq 'a')){
    $ref_sv = $ref_sv_simA;
    if ($ref_sv2 ne ''){
	$ref_sv = $ref_sv2;
    }
}
elsif (($ref_sv eq 'N') or ($ref_sv eq 'n')){
    $ref_sv = $ref_sv_NA12878;
    if ($ref_sv2 ne ''){
	$ref_sv = $ref_sv2;
    }
}
elsif ($ref_sv =~ /^[Mm][Ee]/){
    $ref_sv = $ref_mei;
}
elsif ($ref_sv =~ /^[Vv][Ee]/){
    $ref_sv = $ref_vei;
}
elsif ($ref_sv =~ /^[Mm][Tt]|^NU/){
    $ref_sv = $ref_numt;
}

$target_chr = 'all' if ($target_chr eq 'ALL');

my %sv_type;

if ($sv_type ne 'ALL'){
    if ($sv_type !~ /,/){
        $sv_type{$sv_type} = 1;
    }
    else{
        my @type = split (/,/, $sv_type);
        map{$sv_type{$_} = 1} @type;
    }
}

print "< Parameter: Min SV length: $min_sv_len, Allowed BP diff: $mei_sd, Ref-SV: MEI >\n" if ($ref_sv eq $ref_mei);
print "< Parameter: Min SV length: $min_sv_len, Allowed BP diff: $mei_sd, Ref-SV: VEI >\n" if ($ref_sv eq $ref_vei);
print "< Parameter: Min SV length: $min_sv_len, Allowed BP diff: $mei_sd, Ref-SV: NUMT >\n" if ($ref_sv eq $ref_numt);

my %ref;
my %refGT;
my %ref_ins;
my %ref_dup;
my %match_ref;
my %ref_mei_num;
my $ref_del_num = 0;
my $ref_ins_num = 0;
my $ref_inv_num = 0;
my $ref_dup_num = 0;
my $ref_del_s = 0;
my $ref_del_ss = 0;
my $ref_ins_s = 0;
my $ref_dup_s = 0;
my $ref_inv_s = 0;
my $ref_del_m = 0;
my $ref_ins_m = 0;
my $ref_dup_m = 0;
my $ref_inv_m = 0;
my $ref_del_l = 0;
my $ref_ins_l = 0;
my $ref_dup_l = 0;
my $ref_inv_l = 0;

my $ref_alu_num = 0;
my $ref_l1_num = 0;
my $ref_sva_num = 0;
my $ref_hervk_num = 0;

my $ref_vei_num = 0;

my $ref_numt_num = 0;

my %ref_info;

my %gap;

my %region;
my $total_region = 0;

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

if ($region_bed ne ''){
    open (FILE, $region_bed) or die "$region_bed is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^#|^$/);
        my ($chr, $start, $end) = split (/\t/, $line);
        ${$region{$chr}}{$start} = $end;
        $total_region += $end - $start + 1;
    }
    close (FILE);
    print STDERR "Total genomic size to be analyzed: $total_region\n";
}

open (FILE, $ref_sv) or die "$ref_sv is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#[^\dXY]+/);
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    next if ($chr =~ /^#/);
    my $chr2 = $chr;
    $chr =~ s/^chr// if ($chr =~ /^chr/);
    next if ($target_chr ne 'all') and (($chr ne $target_chr) and ($target_chr !~ /,$chr,|,$chr$|^$chr,/));
    next if (($chr eq 'Y')) and ($include_y == 0);
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    my $type2 = $line[2];
    $type2 = 'ALU' if ($type2 eq 'Alu') or ($type2 eq 'alu');
    $type2 = 'LINE1' if ($type2 eq 'L1');
    $type = 'DUP' if ($type2 eq 'tandem');
    $type = 'MEI' if ($type2 eq 'ALU') or ($type2 eq 'LINE1') or ($type2 eq 'HERVK') or ($type2 eq 'SVA');
    next if ($type2 eq 'INS-DUP');
    my $svlen = 0;
    $svlen = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
    next if ($svlen < $min_ref_len) and ($type ne 'INS');
    next if ($svlen > $max_ref_len);

    my $end = $pos + $svlen - 1;
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
    
    if ($total_region > 0){
        my $region_overlap = 0;
        if (exists $region{$chr}){
            foreach my $rpos (sort {$a <=> $b} keys %{$region{$chr}}){
                my $rend = ${$region{$chr}}{$rpos};
                if (($pos >= $rpos) and ($pos <= $rend)){
                    if ($type eq 'INS'){
                        $region_overlap = 1;
                        last;
                    }
                    else{
                        if ($rend - $pos + 1 >= $svlen * 0.5){
                            $region_overlap = 1;
                            last;
                        }
                    }
                }
                elsif (($type ne 'INS') and ($rpos >= $pos) and ($rpos <= $end)){
                    if ($end - $rpos + 1 >= $svlen * 0.5){
                        $region_overlap = 1;
                        last;
                    }
                }
                last if ($rpos > $end);
            }
        }
        if ($region_overlap == 0){
            next;
        }
    }
    
    my $GT = 'HM';
    $GT = 'HT' if ($line[7] =~ /HT/);
    ${${$ref{$type}}{$chr}}{$pos} = $svlen;
    ${${$refGT{$type}}{$chr}}{$pos} = $GT if ($ref_sv eq $ref_sv_simA);
    ${${$ref_info{$type}}{$chr}}{$pos} = $line[7];
    ${${$ref{'INS'}}{$chr}}{$pos} = $svlen if ($type =~ /VEI|MEI|NUMT/);
    if (($type eq 'INS') or ($type eq 'MEI')){
        ${$ref_ins{$chr}}{$pos} = $svlen;
        ${$ref_ins{$chr}}{$pos - $svlen} = $svlen;
    }
    elsif ($type eq 'DUP'){
        ${$ref_dup{$chr}}{$pos} = $svlen;
        ${$ref_dup{$chr}}{$pos - $svlen} = $svlen;
    }
    
    $ref_del_num ++ if ($type eq 'DEL');
    $ref_ins_num ++ if ($type eq 'INS');
    $ref_dup_num ++ if ($type eq 'DUP');
    $ref_inv_num ++ if ($type eq 'INV');
    if ($type eq 'MEI'){
        $ref_alu_num ++ if ($type2 eq 'ALU');
        $ref_l1_num ++ if ($type2 eq 'LINE1');
        $ref_sva_num ++ if ($type2 eq 'SVA');
        $ref_hervk_num ++ if ($type2 eq 'HERVK');
        $ref_ins_num ++;
	$ref_mei_num{$type2} ++;
    }
    elsif ($type eq 'VEI'){
        $ref_vei_num ++;
        $ref_ins_num ++;
	$ref_mei_num{$type} ++;
    }
    elsif ($type eq 'NUMT'){
        $ref_numt_num ++;
        $ref_ins_num ++;
	$ref_mei_num{$type} ++;
    }
    if ($svlen <= 1000){
        if (($svlen <= 100) and ($type eq 'DEL')){
            $ref_del_ss ++;
        }
        $ref_del_s ++ if ($type eq 'DEL') and ($svlen > 100);
        $ref_ins_s ++ if ($type eq 'INS');
        $ref_dup_s ++ if ($type eq 'DUP');
        $ref_inv_s ++ if ($type eq 'INV');
    }
    elsif (($svlen > 1000) and ($svlen <= 100000)){
        $ref_del_m ++ if ($type eq 'DEL');
        $ref_ins_m ++ if ($type eq 'INS');
        $ref_dup_m ++ if ($type eq 'DUP');
        $ref_inv_m ++ if ($type eq 'INV');
    }
    elsif ($svlen > 100000){
        $ref_del_l ++ if ($type eq 'DEL');
        $ref_ins_l ++ if ($type eq 'INS');
        $ref_dup_l ++ if ($type eq 'DUP');
        $ref_inv_l ++ if ($type eq 'INV');
    }
}
close (FILE);

my $out_file = "$var_base.eval.txt";
open (OUT, "> $out_file");
open (OUT, "> $out_file");
if ($ref_sv eq $ref_mei){
    print "Ref-ALU: $ref_alu_num\n";
    print "Ref-L1: $ref_l1_num\n";
    print "Ref-SVA: $ref_sva_num\n";
    print "Ref-HERVK: $ref_hervk_num\n";
    print OUT "Ref-ALU: $ref_alu_num\n";
    print OUT "Ref-L1: $ref_l1_num\n";
    print OUT "Ref-SVA: $ref_sva_num\n";
    print OUT "Ref-HERVK: $ref_hervk_num\n";
}
elsif ($ref_sv eq $ref_vei){
    print "Ref-VEI: $ref_vei_num\n";
    print OUT "Ref-VEI: $ref_vei_num\n";
}
elsif ($ref_sv eq $ref_numt){
    print "Ref-NUMT: $ref_numt_num\n";
    print OUT "Ref-NUMT: $ref_numt_num\n";
}
else{
    print "Ref-DEL: $ref_del_num\ttinny (<= 100 bp): $ref_del_ss\tshort (<= 1.0 kp): $ref_del_s\tmiddle (<= 100 kb): $ref_del_m\tlarge (> 100 kb): $ref_del_l\n";
    print "Ref-INS: $ref_ins_num\tshort (<=1.0 kp): $ref_ins_s\tmiddle (<= 100 kb): $ref_ins_m\tlarge (> 100 kb): $ref_ins_l\n";
    print "Ref-DUP: $ref_dup_num\tshort (<=1.0 kp): $ref_dup_s\tmiddle (<= 100 kb): $ref_dup_m\tlarge (> 100 kb): $ref_dup_l\n";
    print "Ref-INV: $ref_inv_num\tshort (<=1.0 kp): $ref_inv_s\tmiddle (<= 100 kb): $ref_inv_m\tlarge (> 100 kb): $ref_inv_l\n";
    print OUT "Ref-DEL: $ref_del_num\ttinny (<= 100 bp): $ref_del_ss\tshort (<=1.0 kp): $ref_del_s\tmiddle (<= 100 kb): $ref_del_m\tlarge (> 100 kb): $ref_del_l\n";
    print OUT "Ref-INS: $ref_ins_num\tshort (<=1.0 kp): $ref_ins_s\tmiddle (<= 100 kb): $ref_ins_m\tlarge (> 100 kb): $ref_ins_l\n";
    print OUT "Ref-DUP: $ref_dup_num\tshort (<=1.0 kp): $ref_dup_s\tmiddle (<= 100 kb): $ref_dup_m\tlarge (> 100 kb): $ref_dup_l\n";
    print OUT "Ref-INV: $ref_inv_num\tshort (<=1.0 kp): $ref_inv_s\tmiddle (<= 100 kb): $ref_inv_m\tlarge (> 100 kb): $ref_inv_l\n";
}

my %match_del_num;
my %match_ins_num;
my %match_inv_num;
my %match_dup_num;
my %match_mei_num;

my %match_GT;
my %match_GTNA;

my %bp_diff;
my %len_diff;

my %call_del_num;
my %call_ins_num;
my %call_inv_num;
my %call_dup_num;
my %call_mei_num;

my %nocall_dup_num;
my %nocall_ins_num;
my %recal_del_num;
my %recal_dup_num;

my %recal_dup_ins;
my %recal_ins_dup;

my %pre_info;
my %overlap;
$overlap{'DEL'} = 0;
$overlap{'DUP'} = 0;
$overlap{'INV'} = 0;
$overlap{'INS'} = 0;
my $gt_flag = 0;

my %parent1;
my %parent2;
my %mendel_err;
my %TP;

if ($var_parent1 ne ''){
    open (FILE, $var_parent1) or die "$var_parent1 is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        if ($line =~ /^#|^$/){
            next;
        }
        my @line = split (/\t/, $line);
        my $chr = $line[0];
        next if ($target_chr ne 'all') and (($chr ne $target_chr) and ($target_chr !~ /,$chr,|,$chr$|^$chr,/));
	my $chr2 = $chr;
	$chr =~ s/^chr// if ($chr =~ /^chr/);
        next if ($chr !~ /^\d+$|[XY]/);
        next if (($chr eq 'Y')) and ($include_y == 0);
        my $pos = $line[1];
        my $type = '';
        $type = $line[2] if ($line[2] =~ /^DEL$|^DUP$|^INS$|^INV$|^TRA$|^ALU$|^LINE1$|^L1$|^SVA$|^HERVK$|^ERVK$|^VEI$|^NUMT$/i);
        $type = $1 if ($type eq '') and ($line[7] =~ /SVTYPE=(.+?);/);
        $type = uc $type;
        $type = 'TRA' if ($type eq 'CTX');
        my $class = '';
        $class = 'MEI' if ($type =~ /^ALU$|^LINE1$|^L1$|^SVA$|^HERVK$|^ERVK$/i);
        $type = 'INS' if ($type =~ /^ALU$|^LINE1$|^L1$|^SVA$|^HERVK$|^ERVK$|^VEI$|^NUMT$/i);
        next if (!exists $sv_type{$type}) and ($sv_type ne 'ALL');
        my $len = 0;
        $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
        next if ($len < $min_sv_len) and ($type ne 'INS');
        next if ($len > $max_sv_len);
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
        if ($total_region > 0){
            my $region_overlap = 0;
            if (exists $region{$chr}){
                foreach my $rpos (sort {$a <=> $b} keys %{$region{$chr}}){
                    my $rend = ${$region{$chr}}{$rpos};
                    if (($pos >= $rpos) and ($pos <= $rend)){
                        if ($type eq 'INS'){
                            $region_overlap = 1;
                            last;
                        }
                        else{
                            if ($rend - $pos + 1 >= $len * 0.5){
                                $region_overlap = 1;
                                last;
                            }
                        }
                    }
                    elsif (($type ne 'INS') and ($rpos >= $pos) and ($rpos <= $end)){
                        if ($end - $rpos + 1 >= $len * 0.5){
                            $region_overlap = 1;
                            last;
                        }
                    }
                    last if ($rpos > $end);
                }
            }
            if ($region_overlap == 0){
                next;
            }
        }
        ${${$parent1{$type}}{$chr}}{$pos} = $len;
    }
    close (FILE);
}

if ($var_parent2 ne ''){
    open (FILE, $var_parent2) or die "$var_parent2 is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        if ($line =~ /^#|^$/){
            next;
        }
        my @line = split (/\t/, $line);
        my $chr = $line[0];
        next if ($target_chr ne 'all') and (($chr ne $target_chr) and ($target_chr !~ /,$chr,|,$chr$|^$chr,/));
	my $chr2 = $chr;
	$chr =~ s/^chr// if ($chr =~ /^chr/);
        next if ($chr !~ /^\d+$|[XY]/);
        next if (($chr eq 'Y')) and ($include_y == 0);
        my $pos = $line[1];
        my $type = '';
        $type = $line[2] if ($line[2] =~ /^DEL$|^DUP$|^INS$|^INV$|^TRA$|^ALU$|^LINE1$|^L1$|^SVA$|^HERVK$|^ERVK$|^VEI$|^NUMT$/i);
        $type = $1 if ($type eq '') and ($line[7] =~ /SVTYPE=(.+?);/);
        $type = uc $type;
        $type = 'TRA' if ($type eq 'CTX');
        my $class = '';
        $class = 'MEI' if ($type =~ /^ALU$|^LINE1$|^L1$|^SVA$|^HERVK$|^ERVK$/i);
        $type = 'INS' if ($type =~ /^ALU$|^LINE1$|^L1$|^SVA$|^HERVK$|^ERVK$|^VEI$|^NUMT$/i);
        next if (!exists $sv_type{$type}) and ($sv_type ne 'ALL');
        my $len = 0;
        $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
        next if ($len < $min_sv_len) and ($type ne 'INS');
        next if ($len > $max_sv_len);
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
        if ($total_region > 0){
            my $region_overlap = 0;
            if (exists $region{$chr}){
                foreach my $rpos (sort {$a <=> $b} keys %{$region{$chr}}){
                    my $rend = ${$region{$chr}}{$rpos};
                    if (($pos >= $rpos) and ($pos <= $rend)){
                        if ($type eq 'INS'){
                            $region_overlap = 1;
                            last;
                        }
                        else{
                            if ($rend - $pos + 1 >= $len * 0.5){
                                $region_overlap = 1;
                                last;
                            }
                        }
                    }
                    elsif (($type ne 'INS') and ($rpos >= $pos) and ($rpos <= $end)){
                        if ($end - $rpos + 1 >= $len * 0.5){
                            $region_overlap = 1;
                            last;
                        }
                    }
                    last if ($rpos > $end);
                }
            }
            if ($region_overlap == 0){
                next;
            }
        }
        ${${$parent2{$type}}{$chr}}{$pos} = $len;
    }
    close (FILE);
}

open (FILE, $var_child) or die "$var_child is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#|^$/){
        next;
    }
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    next if ($target_chr ne 'all') and (($chr ne $target_chr) and ($target_chr !~ /,$chr,|,$chr$|^$chr,/));
    my $chr2 = $chr;
    $chr =~ s/^chr// if ($chr =~ /^chr/);
    next if ($chr !~ /^\d+$|[XY]/);
    next if (($chr eq 'Y')) and ($include_y == 0);
    my $pos = $line[1];
    my $type = '';
    $type = $line[2] if ($line[2] =~ /^DEL$|^DUP$|^INS$|^INV$|^TRA$|^ALU$|^LINE1$|^L1$|^SVA$|^HERVK$|^ERVK$|^VEI$|^NUMT$/i);
    $type = $1 if ($line[4] =~ /INS:ME:(.+)/);
    $type = $1 if ($type eq '') and ($line[7] =~ /SVTYPE=(.+?);/);
    $type = uc $type;
    $type = 'TRA' if ($type eq 'CTX');
    my $class = '';
    $class = uc $type if ($type =~ /^ALU$|^LINE1$|^L1$|^SVA$|^HERVK$|^ERVK$|^VEI$|^NUMT$/i);
    $class = 'LINE1' if ($class eq 'L1');
    $class = 'HERVK' if ($class eq 'ERVK');
    $type = 'INS' if ($type =~ /^ALU$|^LINE1$|^L1$|^SVA$|^HERVK$|^ERVK$|^VEI$|^NUMT$/i);
    next if (!exists $sv_type{$type}) and ($sv_type ne 'ALL');
    my $len = 0;
    $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
    next if ($len < $min_sv_len) and ($type ne 'INS');
    next if ($len > $max_sv_len);
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
    my $end2 = $pos + $len - 1;
    my $reads = $1 if ($line[7] =~ /READS=(\d+)/);
    next if ($reads < $min_reads);
    my $GT = 'NA';
    if (@line > 9){
	$GT = 'HM' if ($line[7] =~ /GT=1\/1/) or ($line[9] =~ /^1\/1/);
	$GT = 'HT' if ($line[7] =~ /GT=0\/1/) or ($line[7] =~ /GT=1\/0/) or ($line[9] =~ /^1\/0/) or ($line[9] =~ /^0\/1/);
    }
    else{
	$GT = 'HM' if ($line[7] =~ /GT=1\/1/);
	$GT = 'HT' if ($line[7] =~ /GT=0\/1/) or ($line[7] =~ /GT=1\/0/);
    }
    $gt_flag = 1 if ($GT eq 'HT') or ($GT eq 'HM');
    my $gap_overlap = 0;
    foreach my $gstart (sort {$a <=> $b} keys %{$gap{$chr}}){
        my $gend = ${$gap{$chr}}{$gstart};
        if (($pos >= $gstart) and ($pos <= $gend)){
            if ($type eq 'INS'){
                $gap_overlap = 1;
                last;
            }
            else{
                if ($gend - $pos + 1 >= $len * 0.3){
                    $gap_overlap = 1;
                    last;
                }
            }
        }
        elsif (($type ne 'INS') and ($gstart >= $pos) and ($gstart <= $end)){
            if ($end - $gstart + 1 >= $len * 0.3){
                $gap_overlap = 1;
                last;
            }
        }
        last if ($gstart > $end);
    }
    next if ($gap_overlap == 1);
    if ($total_region > 0){
        my $region_overlap = 0;
        if (exists $region{$chr}){
            foreach my $rpos (sort {$a <=> $b} keys %{$region{$chr}}){
                my $rend = ${$region{$chr}}{$rpos};
                if (($pos >= $rpos) and ($pos <= $rend)){
                    if ($type eq 'INS'){
                        $region_overlap = 1;
                        last;
                    }
                    else{
                        if ($rend - $pos + 1 >= $len * 0.5){
                            $region_overlap = 1;
                            last;
                        }
                    }
                }
                elsif (($type ne 'INS') and ($rpos >= $pos) and ($rpos <= $end)){
                    if ($end - $rpos + 1 >= $len * 0.5){
                        $region_overlap = 1;
                        last;
                    }
                }
                last if ($rpos > $end);
            }
        }
        if ($region_overlap == 0){
            next;
        }
    }
    my $size = 'S';
    if ($type ne 'INS'){
        if ($len <= 1000){
            $size = 'S';
            if (($len <= 100) and ($type eq 'DEL')){
                $size = 'SS';
            }
        }
        elsif (($len > 1000) and ($len <= 100000)){
            $size = 'M';
        }
        elsif ($len > 100000){
            $size = 'L';
        }
    }
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
        if ($type eq 'INS'){
            foreach my $minread (sort {$a <=> $b} @min_reads){
                if ($reads >= $minread){
                    $call_ins_num{$minread} ++;
		    ${$call_mei_num{$class}}{$minread} ++ if ($class ne '');
                }
            }
        }
        else{
            foreach my $minread (sort {$a <=> $b} @min_reads){
                if ($reads >= $minread){
                    ${$call_del_num{$size}}{$minread} ++ if ($type eq 'DEL');
                    ${$call_dup_num{$size}}{$minread} ++ if ($type eq 'DUP');
                    ${$call_inv_num{$size}}{$minread} ++ if ($type eq 'INV');
                    ${$call_del_num{'A'}}{$minread} ++ if ($type eq 'DEL');
                    ${$call_dup_num{'A'}}{$minread} ++ if ($type eq 'DUP');
                    ${$call_inv_num{'A'}}{$minread} ++ if ($type eq 'INV');
                }
            }
        }
    }
    my $min_overlap_ratio2 = $min_overlap_ratio;
    if (($ref_sv_tag eq 'A') and ($min_overlap == 0.5)){
        $min_overlap_ratio2 = 0.8  if ($size eq 'M') or ($size eq 'L');
        $min_overlap_ratio2 = 0.6  if ($size eq 'SS') or ($size eq 'S');
    }
    
    my $trio_flag = 0;
    if ($var_parent1 ne ''){
        foreach my $fpos (sort {$a <=> $b} keys %{${$parent1{$type}}{$chr}}){
            last if ($fpos > $end + $ins_sd);
            my $flen = ${${$parent1{$type}}{$chr}}{$fpos};
            my $fend = $fpos + $flen - 1;
            next if ($fend < $pos - $ins_sd);
            if ($type eq 'INS'){
                if (abs ($pos - $fpos) <= $ins_sd){
                    $trio_flag = 1;
                    last;
                }
            }
            else{
                if (($pos <= $fpos) and ($end >= $fend)){
                    if ($flen >= $len * $min_overlap_ratio2){
                        $trio_flag = 1;
                    }
                }
                elsif (($pos >= $fpos) and ($end <= $fend)){
                    if ($len >= $flen * $min_overlap_ratio2){
                        $trio_flag = 1;
                    }
                }
                elsif (($pos >= $fpos) and ($pos <= $fend)){
                    my $overlap = $fend - $pos + 1;
                    if (($overlap >= $len * $min_overlap_ratio2) or ($overlap >= $flen * $min_overlap_ratio2)){
                        $trio_flag = 1;
                    }
                }
                elsif (($end >= $fpos) and ($end <= $fend)){
                    my $overlap = $end - $fpos + 1;
                    if (($overlap >= $len * $min_overlap_ratio2) or ($overlap >= $flen * $min_overlap_ratio2)){
                        $trio_flag = 1;
                    }
                }
            }
            last if ($trio_flag == 1);
        }
    }
    if (($trio_flag == 0) and ($var_parent2 ne '')){
        foreach my $fpos (sort {$a <=> $b} keys %{${$parent2{$type}}{$chr}}){
            last if ($fpos > $end + $ins_sd);
            my $flen = ${${$parent2{$type}}{$chr}}{$fpos};
            my $fend = $fpos + $flen - 1;
            next if ($fend < $pos - $ins_sd);
            if ($type eq 'INS'){
                if (abs ($pos - $fpos) <= $ins_sd){
                    $trio_flag = 1;
                    last;
                }
            }
            else{
                if (($pos <= $fpos) and ($end >= $fend)){
                    if ($flen >= $len * $min_overlap_ratio2){
                        $trio_flag = 1;
                    }
                }
                elsif (($pos >= $fpos) and ($end <= $fend)){
                    if ($len >= $flen * $min_overlap_ratio2){
                        $trio_flag = 1;
                    }
                }
                elsif (($pos >= $fpos) and ($pos <= $fend)){
                    my $overlap = $fend - $pos + 1;
                    if (($overlap >= $len * $min_overlap_ratio2) or ($overlap >= $flen * $min_overlap_ratio2)){
                        $trio_flag = 1;
                    }
                }
                elsif (($end >= $fpos) and ($end <= $fend)){
                    my $overlap = $end - $fpos + 1;
                    if (($overlap >= $len * $min_overlap_ratio2) or ($overlap >= $flen * $min_overlap_ratio2)){
                        $trio_flag = 1;
                    }
                }
            }
            last if ($trio_flag == 1);
        }
    }
    if (($trio_flag == 0) and ($parent_flag == 1)){
        foreach my $minread (sort {$a <=> $b} @min_reads){
            if ($reads >= $minread){
                ${${$mendel_err{$type}}{'A'}}{$minread} ++;
                ${${$mendel_err{$type}}{$size}}{$minread} ++ if ($type eq 'DEL') or ($type eq 'DUP');
            }
        }
        next if ($substract_ME == 1);
    }
    
    my $flag = 0;
    my $flag2 = 0;
    my $hit_bp = 0;
    my $hit_len = 0;
    if ($type eq 'DEL'){
        my $dellen = 0;
        foreach my $bp (sort {$a <=> $b} keys %{${$ref{'DEL'}}{$chr}}){
            last if ($bp > $end + $ins_sd);
            next if (exists ${${$match_ref{$type}}{$chr}}{$bp});
            $dellen = ${${$ref{'DEL'}}{$chr}}{$bp};
            my $delend = $bp + $dellen - 1;
            next if ($delend < $pos - $ins_sd);
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
                last;
            }
        }
    }
    elsif ($type eq 'DUP'){
        my $duplen = 0;
        foreach my $bp (sort {$a <=> $b} keys %{${$ref{'DUP'}}{$chr}}){
            last if ($bp > $end + $ins_sd);
            next if (exists ${${$match_ref{$type}}{$chr}}{$bp});
            $duplen = ${${$ref{'DUP'}}{$chr}}{$bp};
            my $dupend = $bp + $duplen - 1;
            next if ($dupend < $pos - $ins_sd);
            $hit_bp = $bp;
            $hit_len = $duplen;
            if ((abs ($pos - $bp) <= $var_sd) and (abs ($end - $dupend) <= $var_sd)){
                $flag = 1;
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
                next if ($insend < $pos - $ins_sd);
                $hit_len = $inslen;
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
                    ${${$match_ref{'INS'}}{$chr}}{$bp1} = 1 if (exists ${$ref_ins{$chr}}{$bp1});
                    ${${$match_ref{'INS'}}{$chr}}{$bp2} = 1 if (exists ${$ref_ins{$chr}}{$bp2});
                    last;
                }
            }
        }
    }
    elsif ($type eq 'INS'){
        my $duplen = 0;
        foreach my $bp (sort {$a <=> $b} keys %{${$ref{'INS'}}{$chr}}){
            last if ($bp > $end2 + $ins_sd);
            next if (exists ${${$match_ref{$type}}{$chr}}{$bp});
            my $inslen = ${${$ref{'INS'}}{$chr}}{$bp};
	    my $insend = $bp + $inslen - 1;
            next if ($ins_eval == 0) and ($bp < $pos - $ins_sd);
	    next if ($ins_eval == 1) and ($insend < $pos);
            $hit_bp = $bp;
            if (abs ($pos - $bp) <= $ins_sd){
                if (($ins_eval == 0) or ($len <= 1)){
		    $flag = 1;
		}
		elsif ($ins_eval == 1){
		    my $len_rate = int ($inslen / $len * 100 + 0.5) / 100;
		    if (($len_rate <= 2) and ($len_rate >= 0.5)){
			$flag = 1;
		    }
		}
            }
            if ($flag == 1){
                ${${$match_ref{$type}}{$chr}}{$bp} = 1;
                last;
            }
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
                    if (($ins_eval == 0) or ($len <= 1)){
			$flag = 2;
		    }
		    elsif ($ins_eval == 1){
			my $len_rate = int ($duplen / $len * 100 + 0.5) / 100;
			if (($len_rate <= 2) and ($len_rate >= 0.5)){
			    $flag = 2;
			}
		    }
                }
                if ($flag == 2){
		    ${${$match_ref{'DUP'}}{$chr}}{$bp} = 1;
		    last;
		}
            }
        }
    }
    elsif ($type eq 'INV'){
        my $invlen = 0;
        foreach my $bp (sort {$a <=> $b} keys %{${$ref{'INV'}}{$chr}}){
            last if ($bp > $end + $ins_sd);
            next if (exists ${${$match_ref{$type}}{$chr}}{$bp});
            $invlen = ${${$ref{'INV'}}{$chr}}{$bp};
            my $invend = $bp + $invlen - 1;
            next if ($invend < $pos - $ins_sd);
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
                    if ($invlen >= $len * $min_overlap_ratio2){
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
        foreach my $minread (sort {$a <=> $b} @min_reads){
            if ($type eq 'INS'){
                if ($reads >= $minread){
                    $match_ins_num{$minread} ++;
		    ${$match_mei_num{$class}}{$minread} ++ if ($class ne '');
                }
            }
            else{
                if ($reads >= $minread){
                    ${$match_del_num{'A'}}{$minread} ++ if ($type eq 'DEL');
                    ${$match_dup_num{'A'}}{$minread} ++ if ($type eq 'DUP');
                    ${$match_inv_num{'A'}}{$minread} ++ if ($type eq 'INV');
                    ${$match_dup_num{$size}}{$minread} ++ if ($type eq 'DUP');
                    ${$match_del_num{$size}}{$minread} ++ if ($type eq 'DEL');
                    ${$match_inv_num{$size}}{$minread} ++ if ($type eq 'INV');
                    if ($size eq 'SS'){
                        if ($hit_len > 100){
                            ${$recal_del_num{$size}}{$minread} -- if ($type eq 'DEL');
                            ${$recal_del_num{'S'}}{$minread} ++ if ($type eq 'DEL');
                        }
                    }
                    elsif ($size eq 'S'){
                        if ($hit_len <= 100){
                            ${$recal_del_num{$size}}{$minread} -- if ($type eq 'DEL');
                            ${$recal_del_num{'SS'}}{$minread} ++ if ($type eq 'DEL');
                        }
                        elsif ($hit_len > 1000){
                            ${$recal_del_num{$size}}{$minread} -- if ($type eq 'DEL');
                            ${$recal_del_num{'M'}}{$minread} ++ if ($type eq 'DEL');
                            ${$recal_dup_num{$size}}{$minread} -- if ($type eq 'DUP');
                            ${$recal_dup_num{'M'}}{$minread} ++ if ($type eq 'DUP');
                        }
                    }
                    elsif ($size eq 'M'){
                        if ($hit_len <= 1000){
                            ${$recal_del_num{$size}}{$minread} -- if ($type eq 'DEL');
                            ${$recal_del_num{'S'}}{$minread} ++ if ($type eq 'DEL');
                            ${$recal_dup_num{$size}}{$minread} -- if ($type eq 'DUP');
                            ${$recal_dup_num{'S'}}{$minread} ++ if ($type eq 'DUP');
                        }
                        elsif ($hit_len > 100000){
                            ${$recal_del_num{$size}}{$minread} -- if ($type eq 'DEL');
                            ${$recal_del_num{'L'}}{$minread} ++ if ($type eq 'DEL');
                            ${$recal_dup_num{$size}}{$minread} -- if ($type eq 'DUP');
                            ${$recal_dup_num{'L'}}{$minread} ++ if ($type eq 'DUP');
                        }
                    }
                    elsif ($size eq 'L'){
                        if ($hit_len <= 100000){
                            ${$recal_del_num{$size}}{$minread} -- if ($type eq 'DEL');
                            ${$recal_del_num{'M'}}{$minread} ++ if ($type eq 'DEL');
                            ${$recal_dup_num{$size}}{$minread} -- if ($type eq 'DUP');
                            ${$recal_dup_num{'M'}}{$minread} ++ if ($type eq 'DUP');
                        }
                    }
                }
            }
        }
        if (($eval_bp == 1) and ($flag == 1)){
            if ($type eq 'DEL'){
                $size = 'S' if ($size eq 'SS');
                push @{${$bp_diff{$type}}{$size}}, abs ($pos - $hit_bp) if ($hit_bp > 0);
                push @{${$bp_diff{$type}}{'A'}}, abs ($pos - $hit_bp) if ($hit_bp > 0);
                push @{${$len_diff{$type}}{$size}}, abs ($len - $hit_len) if ($hit_len > 0);
                push @{${$len_diff{$type}}{'A'}}, abs ($len - $hit_len) if ($hit_len > 0);
            }
            elsif ($type eq 'DUP'){
                push @{${$bp_diff{$type}}{$size}}, abs ($pos - $hit_bp) if ($hit_bp > 0);
                push @{${$bp_diff{$type}}{'A'}}, abs ($pos - $hit_bp) if ($hit_bp > 0);
                push @{${$len_diff{$type}}{$size}}, abs ($len - $hit_len) if ($hit_len > 0);
                push @{${$len_diff{$type}}{'A'}}, abs ($len - $hit_len) if ($hit_len > 0);
            }
            elsif ($type eq 'INV'){
                push @{${$bp_diff{$type}}{'A'}}, abs ($pos - $hit_bp) if ($hit_bp > 0);
                push @{${$len_diff{$type}}{'A'}}, abs ($len - $hit_len) if ($hit_len > 0);
            }
            elsif ($type eq 'INS'){
                push @{${$bp_diff{$type}}{'A'}}, abs ($pos - $hit_bp) if ($hit_bp > 0);
                push @{${$len_diff{$type}}{'A'}}, abs ($len - $hit_len) if ($hit_len > 10) and ($len > 10);
                my $diff_bp = abs ($pos - $hit_bp);
#print STDERR "$pos\t$diff_bp\n";
            }
        }
        if ($eval_gt == 1){
            my $ref_GT = '';
            $ref_GT = ${${$refGT{$type}}{$chr}}{$hit_bp} if (exists ${${$refGT{$type}}{$chr}}{$hit_bp});
            $ref_GT = ${${$refGT{'INS'}}{$chr}}{$hit_bp} if ($ref_GT eq '') and ($type eq 'DUP');
            $ref_GT = ${${$refGT{'DUP'}}{$chr}}{$hit_bp} if ($ref_GT eq '') and ($type eq 'INS');
            foreach my $minread (sort {$a <=> $b} @min_reads){
                if ($reads >= $minread){
                    ${${$match_GT{$type}}{'A'}}{$minread} ++ if ($GT eq $ref_GT);
                    ${${$match_GTNA{$type}}{'A'}}{$minread} ++ if ($GT eq 'NA');
                    if ($type =~ /DEL|DUP/){
                        ${${$match_GT{$type}}{$size}}{$minread} ++ if ($GT eq $ref_GT);
                        ${${$match_GTNA{$type}}{$size}}{$minread} ++ if ($GT eq 'NA');
                    }
                }
            }
        }
        ${${$TP{$chr}}{$pos}}{$type} = 1;
        $pre_info{$type} = "$chr=$pos=$len=1";
    }
    else{
        $pre_info{$type} = "$chr=$pos=$len=0";
    }
    if ($flag == 2){
        foreach my $minread (sort {$a <=> $b} @min_reads){
            if ($type eq 'INS'){
                if ($reads >= $minread){
                    $nocall_ins_num{$minread} ++;
                }
            }
            elsif ($type eq 'DUP'){
                if ($reads >= $minread){
                    ${$nocall_dup_num{'A'}}{$minread} ++;
                    ${$nocall_dup_num{$size}}{$minread} ++;
                }
            }
        }
        $pre_info{$type} = "$chr=$pos=$len=1";
    }
}
close (FILE);

foreach my $type (keys %mendel_err){
    foreach my $size (keys %{$mendel_err{$type}}){
        foreach my $minread (sort {$a <=> $b} @min_reads){
            ${${$mendel_err{$type}}{$size}}{$minread} = 0 if (!exists ${${$mendel_err{$type}}{$size}}{$minread});
        }
    }
}

my @MEI = ('ALU', 'LINE1', 'SVA', 'HERVK', 'NUMT', 'VEI');

my $ref_prefix = basename ($ref_sv);
$ref_prefix = $1 if ($ref_prefix =~ /(.+?)\./);
print "<< $ref_prefix >>\n\n";
print OUT "\n<< $ref_prefix >>\n\n";

if ((exists $call_del_num{'A'}) and ($ref_del_num > 0)){
    my $type = 'DEL';
    print "## DEL ##\n";
    print "        \t<Number of supporting reads>\n";
    print "        ";
    print OUT "## DEL ##\n";
    print OUT "        \t<Number of supporting reads>\n";
    print OUT "        ";
    foreach my $minread (sort {$a <=> $b} @min_reads){
        print "\t$minread";
        print OUT "\t$minread";
        if (!exists ${$call_del_num{'A'}}{$minread}){
            ${$call_del_num{'A'}}{$minread} = 0;
            ${${$mendel_err{'DEL'}}{'A'}}{$minread} = 0;
        }
        if (!exists ${$call_del_num{'SS'}}{$minread}){
            ${$call_del_num{'SS'}}{$minread} = 0;
            ${${$mendel_err{'DEL'}}{'SS'}}{$minread} = 0;
        }
        if (!exists ${$call_del_num{'S'}}{$minread}){
            ${$call_del_num{'S'}}{$minread} = 0;
            ${${$mendel_err{'DEL'}}{'S'}}{$minread} = 0;
        }
        if (!exists ${$call_del_num{'M'}}{$minread}){
            ${$call_del_num{'M'}}{$minread} = 0;
            ${${$mendel_err{'DEL'}}{'M'}}{$minread} = 0;
        }
        if (!exists ${$call_del_num{'L'}}{$minread}){
            ${$call_del_num{'L'}}{$minread} = 0;
            ${${$mendel_err{'DEL'}}{'L'}}{$minread} = 0;
        }
        if (!exists ${$match_del_num{'A'}}{$minread}){
            ${$match_del_num{'A'}}{$minread} = 0;
        }
        if (!exists ${$match_del_num{'SS'}}{$minread}){
            ${$match_del_num{'SS'}}{$minread} = 0;
        }
        if (!exists ${$match_del_num{'S'}}{$minread}){
            ${$match_del_num{'S'}}{$minread} = 0;
        }
        if (!exists ${$match_del_num{'M'}}{$minread}){
            ${$match_del_num{'M'}}{$minread} = 0;
        }
        if (!exists ${$match_del_num{'L'}}{$minread}){
            ${$match_del_num{'L'}}{$minread} = 0;
        }
    }
    print "\nCall   (A)\t";
    print OUT "\nCall (A)\t";
    foreach my $read (sort {$a <=> $b} keys %{$call_del_num{'A'}}){
        print "${$call_del_num{'A'}}{$read}\t";
        print OUT "${$call_del_num{'A'}}{$read}\t";
    }
    print "\nRecall (A)\t";
    print OUT "\nRecall (A)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_del_num{'A'}}){
        my $recall = int (${$match_del_num{'A'}}{$read} / $ref_del_num * 1000) / 10;
        print "$recall\t";
        print OUT "$recall\t";
    }
    print "\nPrecis (A)\t";
    print OUT "\nPrecis (A)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_del_num{'A'}}){
        my $precis = 0;
        $precis = int (${$match_del_num{'A'}}{$read} / ${$call_del_num{'A'}}{$read} * 1000) / 10 if (${$call_del_num{'A'}}{$read} > 0);
        print "$precis\t";
        print OUT "$precis\t";
    }
    if ($parent_flag == 1){
        print "\nMIER   (A)\t";
        print OUT "\nMIER (A)\t";
        foreach my $read (sort {$a <=> $b} keys %{${$mendel_err{'DEL'}}{'A'}}){
            my $error = ${${$mendel_err{'DEL'}}{'A'}}{$read};
            my $error_rate = 0;
            $error_rate = int ($error / ${$call_del_num{'A'}}{$read} * 1000) / 10 if (${$call_del_num{'A'}}{$read} > 0);
            print "$error_rate\t";
            print OUT "$error_rate\t";
        }
    }
    if (($eval_gt == 1) and ($gt_flag == 1)){
	print "\nGenotying (A)\t";
	print OUT "\nGenotying (A)\t";
	foreach my $read (sort {$a <=> $b} keys %{${$match_GT{$type}}{'A'}}){
            my $matchGT = 0;
	    $matchGT = ${${$match_GT{$type}}{'A'}}{$read} if (exists ${${$match_GT{$type}}{'A'}}{$read});
	    my $match = ${$match_del_num{'A'}}{$read};
            my $NA = 0;
	    $NA = ${${$match_GTNA{$type}}{'A'}}{$read} if (exists ${${$match_GTNA{$type}}{'A'}}{$read});
	    my $precis = 0;
	    $precis = int ($matchGT / ($match - $NA) * 1000) / 10 if ($match - $NA > 0);
	    print "$precis($match $matchGT $NA)\t";
	    print OUT "$precis($match $matchGT $NA)\t";
	}
    }
    print "\nCall   (SS)\t";
    print OUT "\nCall (SS)\t";
    foreach my $read (sort {$a <=> $b} keys %{$call_del_num{'SS'}}){
        print "${$call_del_num{'SS'}}{$read}\t";
        print OUT "${$call_del_num{'SS'}}{$read}\t";
    }
    print "\nRecall (SS)\t";
    print OUT "\nRecall (SS)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_del_num{'SS'}}){
        my $recall = 0;
        my $norecal = 0;
        $norecal = ${$recal_del_num{'SS'}}{$read} if (exists ${$recal_del_num{'SS'}}{$read});
        $recall = int ((${$match_del_num{'SS'}}{$read} + $norecal) / $ref_del_ss * 1000) / 10 if ($ref_del_ss > 0);
        print "$recall\t";
        print OUT "$recall\t";
    }
    print "\nPrecis (SS)\t";
    print OUT "\nPrecis (SS)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_del_num{'SS'}}){
        my $precis = 0;
        $precis = int (${$match_del_num{'SS'}}{$read} / ${$call_del_num{'SS'}}{$read} * 1000) / 10 if (${$call_del_num{'SS'}}{$read} > 0);
        print "$precis\t";
        print OUT "$precis\t";
    }
    if ($parent_flag == 1){
        print "\nMIER   (SS)\t";
        print OUT "\nMIER (SS)\t";
        foreach my $read (sort {$a <=> $b} keys %{${$mendel_err{'DEL'}}{'SS'}}){
            my $error = 0;
            $error = ${${$mendel_err{'DEL'}}{'SS'}}{$read} if (exists ${${$mendel_err{'DEL'}}{'SS'}}{$read});
            my $error_rate = 0;
            $error_rate = int ($error / ${$call_del_num{'SS'}}{$read} * 1000) / 10 if (${$call_del_num{'SS'}}{$read} > 0);
            print "$error_rate\t";
            print OUT "$error_rate\t";
        }
    }
    if (($eval_gt == 1) and ($gt_flag == 1)){
	print "\nGenotying (SS)\t";
	print OUT "\nGenotying (SS)\t";
	foreach my $read (sort {$a <=> $b} keys %{${$match_GT{$type}}{'SS'}}){
            my $matchGT = 0;
	    $matchGT = ${${$match_GT{$type}}{'SS'}}{$read} if (exists ${${$match_GT{$type}}{'SS'}}{$read});
	    my $match = ${$match_del_num{'SS'}}{$read};
            my $NA = 0;
	    $NA = ${${$match_GTNA{$type}}{'SS'}}{$read} if (exists ${${$match_GTNA{$type}}{'SS'}}{$read});
	    my $precis = 0;
	    $precis = int ($matchGT / ($match - $NA) * 1000) / 10 if ($match - $NA > 0);
	    print "$precis($match $matchGT $NA)\t";
	    print OUT "$precis($match $matchGT $NA)\t";
	}
    }
    print "\nCall   (S)\t";
    print OUT "\nCall (S)\t";
    foreach my $read (sort {$a <=> $b} keys %{$call_del_num{'S'}}){
        print "${$call_del_num{'S'}}{$read}\t";
        print OUT "${$call_del_num{'S'}}{$read}\t";
    }
    print "\nRecall (S)\t";
    print OUT "\nRecall (S)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_del_num{'S'}}){
        my $recall = 0;
        my $norecal = 0;
        $norecal = ${$recal_del_num{'S'}}{$read} if (exists ${$recal_del_num{'S'}}{$read});
        $recall = int ((${$match_del_num{'S'}}{$read} + $norecal) / $ref_del_s * 1000) / 10 if ($ref_del_s > 0);
        print "$recall\t";
        print OUT "$recall\t";
    }
    print "\nPrecis (S)\t";
    print OUT "\nPrecis (S)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_del_num{'S'}}){
        my $precis = 0;
        $precis = int (${$match_del_num{'S'}}{$read} / ${$call_del_num{'S'}}{$read} * 1000) / 10 if (${$call_del_num{'S'}}{$read} > 0);
        print "$precis\t";
        print OUT "$precis\t";
    }
    if ($parent_flag == 1){
        print "\nMIER   (S)\t";
        print OUT "\nMIER (S)\t";
        foreach my $read (sort {$a <=> $b} keys %{${$mendel_err{'DEL'}}{'S'}}){
            my $error = ${${$mendel_err{'DEL'}}{'S'}}{$read};
            my $error_rate = 0;
            $error_rate = int ($error / ${$call_del_num{'S'}}{$read} * 1000) / 10 if (${$call_del_num{'S'}}{$read} > 0);
            print "$error_rate\t";
            print OUT "$error_rate\t";
        }
    }
    if (($eval_gt == 1) and ($gt_flag == 1)){
	print "\nGenotying (S)\t";
	print OUT "\nGenotying (S)\t";
	foreach my $read (sort {$a <=> $b} keys %{${$match_GT{$type}}{'S'}}){
            my $matchGT = 0;
	    $matchGT = ${${$match_GT{$type}}{'S'}}{$read} if (exists ${${$match_GT{$type}}{'S'}}{$read});
	    my $match = ${$match_del_num{'S'}}{$read};
            my $NA = 0;
	    $NA = ${${$match_GTNA{$type}}{'S'}}{$read} if (exists ${${$match_GTNA{$type}}{'S'}}{$read});
	    my $precis = 0;
	    $precis = int ($matchGT / ($match - $NA) * 1000) / 10 if ($match - $NA > 0);
	    print "$precis($match $matchGT $NA)\t";
	    print OUT "$precis($match $matchGT $NA)\t";
	}
    }
    print "\nCall   (M)\t";
    print OUT "\nCall (M)\t";
    foreach my $read (sort {$a <=> $b} keys %{$call_del_num{'M'}}){
        print "${$call_del_num{'M'}}{$read}\t";
        print OUT "${$call_del_num{'M'}}{$read}\t";
    }
    print "\nRecall (M)\t";
    print OUT "\nRecall (M)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_del_num{'M'}}){
        my $recall = 0;
        my $norecal = 0;
        $norecal = ${$recal_del_num{'M'}}{$read} if (exists ${$recal_del_num{'M'}}{$read});
        $recall = int ((${$match_del_num{'M'}}{$read} + $norecal) / $ref_del_m * 1000) / 10 if ($ref_del_m > 0);
        print "$recall\t";
        print OUT "$recall\t";
    }
    print "\nPrecis (M)\t";
    print OUT "\nPrecis (M)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_del_num{'M'}}){
        my $precis = 0;
        $precis = int (${$match_del_num{'M'}}{$read} / ${$call_del_num{'M'}}{$read} * 1000) / 10 if (${$call_del_num{'M'}}{$read} > 0);
        print "$precis\t";
        print OUT "$precis\t";
    }
    if ($parent_flag == 1){
        print "\nMIER   (M)\t";
        print OUT "\nMIER (M)\t";
        foreach my $read (sort {$a <=> $b} keys %{${$mendel_err{'DEL'}}{'M'}}){
            my $error = ${${$mendel_err{'DEL'}}{'M'}}{$read};
            my $error_rate = 0;
            $error_rate = int ($error / ${$call_del_num{'M'}}{$read} * 1000) / 10 if (${$call_del_num{'M'}}{$read} > 0);
            print "$error_rate\t";
            print OUT "$error_rate\t";
        }
    }
    if (($eval_gt == 1) and ($gt_flag == 1)){
	print "\nGenotying (M)\t";
	print OUT "\nGenotying (M)\t";
	foreach my $read (sort {$a <=> $b} keys %{${$match_GT{$type}}{'M'}}){
            my $matchGT = 0;
	    $matchGT = ${${$match_GT{$type}}{'M'}}{$read} if (exists ${${$match_GT{$type}}{'M'}}{$read});
	    my $match = ${$match_del_num{'M'}}{$read};
            my $NA = 0;
	    $NA = ${${$match_GTNA{$type}}{'M'}}{$read} if (exists ${${$match_GTNA{$type}}{'M'}}{$read});
	    my $precis = 0;
	    $precis = int ($matchGT / ($match - $NA) * 1000) / 10 if ($match - $NA > 0);
	    print "$precis($match $matchGT $NA)\t";
	    print OUT "$precis($match $matchGT $NA)\t";
	}
    }
    if ($ref_del_l > 0){
        print "\nCall   (L)\t";
        print OUT "\nCall (L)\t";
        foreach my $read (sort {$a <=> $b} keys %{$call_del_num{'L'}}){
            print "${$call_del_num{'L'}}{$read}\t";
            print OUT "${$call_del_num{'L'}}{$read}\t";
        }
        print "\nRecall (L)\t";
        print OUT "\nRecall (L)\t";
        foreach my $read (sort {$a <=> $b} keys %{$match_del_num{'L'}}){
            my $recall = 0;
            my $norecal = 0;
            $norecal = ${$recal_del_num{'L'}}{$read} if (exists ${$recal_del_num{'L'}}{$read});
            $recall = int ((${$match_del_num{'L'}}{$read} + $norecal) / $ref_del_l * 1000) / 10 if ($ref_del_l > 0);
            print "$recall\t";
            print OUT "$recall\t";
        }
        print "\nPrecis (L)\t";
        print OUT "\nPrecis (L)\t";
        foreach my $read (sort {$a <=> $b} keys %{$match_del_num{'L'}}){
            my $precis = 0;
            $precis = int (${$match_del_num{'L'}}{$read} / ${$call_del_num{'L'}}{$read} * 1000) / 10 if (${$call_del_num{'L'}}{$read} > 0);
            print "$precis\t";
            print OUT "$precis\t";
        }
        if ($parent_flag == 1){
            print "\nMIER   (L)\t";
            print OUT "\nMIER (L)\t";
            foreach my $read (sort {$a <=> $b} keys %{${$mendel_err{'DEL'}}{'L'}}){
                my $error = ${${$mendel_err{'DEL'}}{'L'}}{$read};
                my $error_rate = 0;
                $error_rate = int ($error / ${$call_del_num{'L'}}{$read} * 1000) / 10 if (${$call_del_num{'L'}}{$read} > 0);
                print "$error_rate\t";
                print OUT "$error_rate\t";
            }
        }
        if (($eval_gt == 1) and ($gt_flag == 1)){
            print "\nGenotying (L)\t";
            print OUT "\nGenotying (L)\t";
            foreach my $read (sort {$a <=> $b} keys %{${$match_GT{$type}}{'L'}}){
                my $matchGT = 0;
                $matchGT = ${${$match_GT{$type}}{'L'}}{$read} if (exists ${${$match_GT{$type}}{'L'}}{$read});
                my $match = ${$match_del_num{'L'}}{$read};
                my $NA = 0;
                $NA = ${${$match_GTNA{$type}}{'L'}}{$read} if (exists ${${$match_GTNA{$type}}{'L'}}{$read});
                my $precis = 0;
                $precis = int ($matchGT / ($match - $NA) * 1000) / 10 if ($match - $NA > 0);
                print "$precis($match $matchGT $NA)\t";
                print OUT "$precis($match $matchGT $NA)\t";
            }
        }
    }
    if ($eval_bp == 1){
	foreach my $size (@size){
        print "\n";
        print OUT "\n";
	    my $SDbp = 0;
	    my $SDlen = 0;
	    if (exists ${$bp_diff{$type}}{$size}){
		my $sq_total = 0;
		foreach (@{${$bp_diff{$type}}{$size}}){
		    $sq_total += $_ ** 2;
		}
		$SDbp = int (($sq_total / @{${$bp_diff{$type}}{$size}}) ** 0.5 + 0.5) if (@{${$bp_diff{$type}}{$size}} > 0);
	    }
	    if (exists ${$len_diff{$type}}{$size}){
		my $sq_total = 0;
		foreach (@{${$len_diff{$type}}{$size}}){
		    $sq_total += $_ ** 2;
		}
		$SDlen = int (($sq_total / @{${$len_diff{$type}}{$size}}) ** 0.5 + 0.5) if (@{${$len_diff{$type}}{$size}} > 0);
	    }
	    if (exists ${$bp_diff{$type}}{$size}){
		print "[Break Points SD ($size)]: $SDbp\t[Size SD ($size)]: $SDlen\n";
		print OUT "[Break Points SD ($size)]: $SDbp\t[Size SD ($size)]: $SDlen\n";
	    }
	}
    }
    
    print "\n\n";
    print OUT "\n\n";
}

if ((exists $call_dup_num{'A'}) and ($ref_dup_num > 0)){
     my $type = 'DUP';
    print "## DUP ##\n";
    print "        \t<Number of supporting reads>\n";
    print "        ";
    print OUT "## DUP ##\n";
    print OUT "        \t<Number of supporting reads>\n";
    print OUT "        ";
    foreach my $minread (sort {$a <=> $b} @min_reads){
        print "\t$minread";
        print OUT "\t$minread";
        if (!exists ${$call_dup_num{'A'}}{$minread}){
            ${$call_dup_num{'A'}}{$minread} = 0;
            ${${$mendel_err{'DUP'}}{'A'}}{$minread} = 0;
        }
        if (!exists ${$call_dup_num{'S'}}{$minread}){
            ${$call_dup_num{'S'}}{$minread} = 0;
            ${${$mendel_err{'DUP'}}{'S'}}{$minread} = 0;
        }
        if (!exists ${$call_dup_num{'M'}}{$minread}){
            ${$call_dup_num{'M'}}{$minread} = 0;
            ${${$mendel_err{'DUP'}}{'M'}}{$minread} = 0;
        }
        if (!exists ${$call_dup_num{'L'}}{$minread}){
            ${$call_dup_num{'L'}}{$minread} = 0;
            ${${$mendel_err{'DUP'}}{'L'}}{$minread} = 0;
        }
        if (!exists ${$match_dup_num{'A'}}{$minread}){
            ${$match_dup_num{'A'}}{$minread} = 0;
        }
        if (!exists ${$match_dup_num{'S'}}{$minread}){
            ${$match_dup_num{'S'}}{$minread} = 0;
        }
        if (!exists ${$match_dup_num{'M'}}{$minread}){
            ${$match_dup_num{'M'}}{$minread} = 0;
        }
        if (!exists ${$match_dup_num{'L'}}{$minread}){
            ${$match_dup_num{'L'}}{$minread} = 0;
        }
    }
    print "\nCall   (A)\t";
    print OUT "\nCall (A)\t";
    foreach my $read (sort {$a <=> $b} keys %{$call_dup_num{'A'}}){
        print "${$call_dup_num{'A'}}{$read}\t";
        print OUT "${$call_dup_num{'A'}}{$read}\t";
    }
    print "\nRecall (A)\t";
    print OUT "\nRecall (A)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_dup_num{'A'}}){
        my $dup_ins_hit = 0;
        $dup_ins_hit = ${$nocall_dup_num{'A'}}{$read} if (exists ${$nocall_dup_num{'A'}}{$read});
        my $recall = int (${$match_dup_num{'A'}}{$read} / ($ref_dup_num + $dup_ins_hit) * 1000) / 10;
        print "$recall\t";
        print OUT "$recall\t";
    }
    print "\nPrecis (A)\t";
    print OUT "\nPrecis (A)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_dup_num{'A'}}){
        my $precis = 0;
        my $call = ${$call_dup_num{'A'}}{$read};
        $precis = int (${$match_dup_num{'A'}}{$read} / $call * 1000) / 10 if ($call > 0);
        print "$precis\t";
        print OUT "$precis\t";
    }
    if ($parent_flag == 1){
        print "\nMIER   (A)\t";
        print OUT "\nMIER (A)\t";
        foreach my $read (sort {$a <=> $b} keys %{${$mendel_err{'DUP'}}{'A'}}){
            my $error = ${${$mendel_err{'DUP'}}{'A'}}{$read};
            my $error_rate = 0;
            $error_rate = int ($error / ${$call_dup_num{'A'}}{$read} * 1000) / 10 if (${$call_dup_num{'A'}}{$read} > 0);
            print "$error_rate\t";
            print OUT "$error_rate\t";
        }
    }
    if (($eval_gt == 1) and ($gt_flag == 1)){
	print "\nGenotying (A)\t";
	print OUT "\nGenotying (A)\t";
	foreach my $read (sort {$a <=> $b} keys %{${$match_GT{$type}}{'A'}}){
            my $matchGT = 0;
	    $matchGT = ${${$match_GT{$type}}{'A'}}{$read} if (exists ${${$match_GT{$type}}{'A'}}{$read});
	    my $match = ${$match_dup_num{'A'}}{$read};
            my $NA = 0;
	    $NA = ${${$match_GTNA{$type}}{'A'}}{$read} if (exists ${${$match_GTNA{$type}}{'A'}}{$read});
	    my $precis = 0;
	    $precis = int ($matchGT / ($match - $NA) * 1000) / 10 if ($match - $NA > 0);
	    print "$precis($match $matchGT $NA)\t";
	    print OUT "$precis($match $matchGT $NA)\t";
	}
    }
    print "\nCall   (S)\t";
    print OUT "\nCall (S)\t";
    foreach my $read (sort {$a <=> $b} keys %{$call_dup_num{'S'}}){
        print "${$call_dup_num{'S'}}{$read}\t";
        print OUT "${$call_dup_num{'S'}}{$read}\t";
    }
    print "\nRecall (S)\t";
    print OUT "\nRecall (S)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_dup_num{'S'}}){
        my $norecal = 0;
	$norecal = ${$recal_dup_num{'S'}}{$read} if (exists ${$recal_dup_num{'S'}}{$read});
        my $dup_ins_hit = 0;
        $dup_ins_hit = ${$nocall_dup_num{'S'}}{$read} if (exists ${$nocall_dup_num{'S'}}{$read});
        my $recall = int ((${$match_dup_num{'S'}}{$read} + $norecal) / ($ref_dup_s + $dup_ins_hit) * 1000) / 10;
        print "$recall\t";
        print OUT "$recall\t";
    }
    print "\nPrecis (S)\t";
    print OUT "\nPrecis (S)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_dup_num{'S'}}){
        my $precis = 0;
        my $call = ${$call_dup_num{'S'}}{$read};
        $precis = int (${$match_dup_num{'S'}}{$read} / $call * 1000) / 10 if ($call > 0);
        print "$precis\t";
        print OUT "$precis\t";
    }
    if ($parent_flag == 1){
        print "\nMIER   (S)\t";
        print OUT "\nMIER (S)\t";
        foreach my $read (sort {$a <=> $b} keys %{${$mendel_err{'DUP'}}{'S'}}){
            my $error = ${${$mendel_err{'DUP'}}{'S'}}{$read};
            my $error_rate = 0;
            $error_rate = int ($error / ${$call_dup_num{'S'}}{$read} * 1000) / 10 if (${$call_dup_num{'S'}}{$read} > 0);
            print "$error_rate\t";
            print OUT "$error_rate\t";
        }
    }
    if (($eval_gt == 1) and ($gt_flag == 1)){
	print "\nGenotying (S)\t";
	print OUT "\nGenotying (S)\t";
	foreach my $read (sort {$a <=> $b} keys %{${$match_GT{$type}}{'S'}}){
            my $matchGT = 0;
	    $matchGT = ${${$match_GT{$type}}{'S'}}{$read} if (exists ${${$match_GT{$type}}{'S'}}{$read});
	    my $match = ${$match_dup_num{'S'}}{$read};
            my $NA = 0;
	    $NA = ${${$match_GTNA{$type}}{'S'}}{$read} if (exists ${${$match_GTNA{$type}}{'S'}}{$read});
	    my $precis = 0;
	    $precis = int ($matchGT / ($match - $NA) * 1000) / 10 if ($match - $NA > 0);
	    print "$precis($match $matchGT $NA)\t";
	    print OUT "$precis($match $matchGT $NA)\t";
	}
    }
    print "\nCall   (M)\t";
    print OUT "\nCall (M)\t";
    foreach my $read (sort {$a <=> $b} keys %{$call_dup_num{'M'}}){
        print "${$call_dup_num{'M'}}{$read}\t";
        print OUT "${$call_dup_num{'M'}}{$read}\t";
    }
    print "\nRecall (M)\t";
    print OUT "\nRecall (M)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_dup_num{'M'}}){
        my $norecal = 0;
	$norecal = ${$recal_dup_num{'M'}}{$read} if (exists ${$recal_dup_num{'M'}}{$read});
        my $dup_ins_hit = 0;
        $dup_ins_hit = ${$nocall_dup_num{'M'}}{$read} if (exists ${$nocall_dup_num{'M'}}{$read});
        my $recall = int ((${$match_dup_num{'M'}}{$read} + $norecal) / ($ref_dup_m + $dup_ins_hit) * 1000) / 10;
        print "$recall\t";
        print OUT "$recall\t";
    }
    print "\nPrecis (M)\t";
    print OUT "\nPrecis (M)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_dup_num{'M'}}){
        my $precis = 0;
        my $call = ${$call_dup_num{'M'}}{$read};
        $precis = int (${$match_dup_num{'M'}}{$read} / $call * 1000) / 10 if ($call > 0);
        print "$precis\t";
        print OUT "$precis\t";
    }
    if ($parent_flag == 1){
        print "\nMIER   (M)\t";
        print OUT "\nMIER (M)\t";
        foreach my $read (sort {$a <=> $b} keys %{${$mendel_err{'DUP'}}{'M'}}){
            my $error = ${${$mendel_err{'DUP'}}{'M'}}{$read};
            my $error_rate = 0;
            $error_rate = int ($error / ${$call_dup_num{'M'}}{$read} * 1000) / 10 if (${$call_dup_num{'M'}}{$read} > 0);
            print "$error_rate\t";
            print OUT "$error_rate\t";
        }
    }
    if (($eval_gt == 1) and ($gt_flag == 1)){
	print "\nGenotying (M)\t";
	print OUT "\nGenotying (M)\t";
	foreach my $read (sort {$a <=> $b} keys %{${$match_GT{$type}}{'M'}}){
            my $matchGT = 0;
	    $matchGT = ${${$match_GT{$type}}{'M'}}{$read} if (exists ${${$match_GT{$type}}{'M'}}{$read});
	    my $match = ${$match_dup_num{'M'}}{$read};
            my $NA = 0;
	    $NA = ${${$match_GTNA{$type}}{'M'}}{$read} if (exists ${${$match_GTNA{$type}}{'M'}}{$read});
	    my $precis = 0;
	    $precis = int ($matchGT / ($match - $NA) * 1000) / 10 if ($match - $NA > 0);
	    print "$precis($match $matchGT $NA)\t";
	    print OUT "$precis($match $matchGT $NA)\t";
	}
    }
    print "\nCall   (L)\t";
    print OUT "\nCall (L)\t";
    foreach my $read (sort {$a <=> $b} keys %{$call_dup_num{'L'}}){
        print "${$call_dup_num{'L'}}{$read}\t";
        print OUT "${$call_dup_num{'L'}}{$read}\t";
    }
    print "\nRecall (L)\t";
    print OUT "\nRecall (L)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_dup_num{'L'}}){
        my $norecal = 0;
	$norecal = ${$recal_dup_num{'L'}}{$read} if (exists ${$recal_dup_num{'L'}}{$read});
        my $dup_ins_hit = 0;
        $dup_ins_hit = ${$nocall_dup_num{'L'}}{$read} if (exists ${$nocall_dup_num{'L'}}{$read});
        my $recall = int ((${$match_dup_num{'L'}}{$read} + $norecal) / ($ref_dup_l + $dup_ins_hit) * 1000) / 10;
        print "$recall\t";
        print OUT "$recall\t";
    }
    print "\nPrecis (L)\t";
    print OUT "\nPrecis (L)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_dup_num{'L'}}){
        my $precis = 0;
        my $call = ${$call_dup_num{'L'}}{$read};
        $precis = int (${$match_dup_num{'L'}}{$read} / $call * 1000) / 10 if ($call > 0);
        print "$precis\t";
        print OUT "$precis\t";
    }
    if ($parent_flag == 1){
        print "\nMIER   (L)\t";
        print OUT "\nMIER (L)\t";
        foreach my $read (sort {$a <=> $b} keys %{${$mendel_err{'DUP'}}{'L'}}){
            my $error = ${${$mendel_err{'DUP'}}{'L'}}{$read};
            my $error_rate = 0;
            $error_rate = int ($error / ${$call_dup_num{'L'}}{$read} * 1000) / 10 if (${$call_dup_num{'L'}}{$read} > 0);
            print "$error_rate\t";
            print OUT "$error_rate\t";
        }
    }
    if (($eval_gt == 1) and ($gt_flag == 1)){
	print "\nGenotying (L)\t";
	print OUT "\nGenotying (L)\t";
	foreach my $read (sort {$a <=> $b} keys %{${$match_GT{$type}}{'L'}}){
            my $matchGT = 0;
	    $matchGT = ${${$match_GT{$type}}{'L'}}{$read} if (exists ${${$match_GT{$type}}{'L'}}{$read});
	    my $match = ${$match_dup_num{'L'}}{$read};
            my $NA = 0;
	    $NA = ${${$match_GTNA{$type}}{'L'}}{$read} if (exists ${${$match_GTNA{$type}}{'L'}}{$read});
	    my $precis = 0;
	    $precis = int ($matchGT / ($match - $NA) * 1000) / 10 if ($match - $NA > 0);
	    print "$precis($match $matchGT $NA)\t";
	    print OUT "$precis($match $matchGT $NA)\t";
	}
    }
    if ($eval_bp == 1){
        print "\n";
        print OUT "\n";
	foreach my $size (@size){
	    my $SDbp = 0;
	    my $SDlen = 0;
	    if (exists ${$bp_diff{$type}}{$size}){
		my $sq_total = 0;
		foreach (@{${$bp_diff{$type}}{$size}}){
		    $sq_total += $_ ** 2;
		}
		$SDbp = int (($sq_total / @{${$bp_diff{$type}}{$size}}) ** 0.5 + 0.5) if (@{${$bp_diff{$type}}{$size}} > 0);
	    }
	    if (exists ${$len_diff{$type}}{$size}){
		my $sq_total = 0;
		foreach (@{${$len_diff{$type}}{$size}}){
		    $sq_total += $_ ** 2;
		}
		$SDlen = int (($sq_total / @{${$len_diff{$type}}{$size}}) ** 0.5 + 0.5) if (@{${$len_diff{$type}}{$size}} > 0);
	    }
	    if (exists ${$bp_diff{$type}}{$size}){
		print "[Break Points SD ($size)]: $SDbp\t[Size SD ($size)]: $SDlen\n";
		print OUT "[Break Points SD ($size)]: $SDbp\t[Size SD ($size)]: $SDlen\n";
	    }
	}
    }
    print "\n\n";
    print OUT "\n\n";
}

if ((exists $call_ins_num{3}) and ($ref_ins_num > 0)){
     my $type = 'INS';
    print "## INS ##\n";
    print "        \t<Number of supporting reads>\n";
    print "        ";
    print OUT "## INS ##\n";
    print OUT "        \t<Number of supporting reads>\n";
    print OUT "        ";
    foreach my $minread (sort {$a <=> $b} @min_reads){
        print "\t$minread";
        print OUT "\t$minread";
        if (!exists $call_ins_num{$minread}){
            $call_ins_num{$minread} = 0;
            ${${$mendel_err{'INS'}}{'A'}}{$minread} = 0;
        }
        if (!exists $match_ins_num{$minread}){
            $match_ins_num{$minread} = 0;
        }
    }
    print "\nCall   (A)\t";
    print OUT "\nCall (A)\t";
    foreach my $read (sort {$a <=> $b} keys %call_ins_num){
        print "$call_ins_num{$read}\t";
        print OUT "$call_ins_num{$read}\t";
    }
    print "\nRecall (A)\t";
    print OUT "\nRecall (A)\t";
    foreach my $read (sort {$a <=> $b} keys %match_ins_num){
        my $dup_ins_hit = 0;
        $dup_ins_hit = ${$nocall_dup_num{'A'}}{$read} if (exists ${$nocall_dup_num{'A'}}{$read});
        my $ins_dup_hit = 0;
        $ins_dup_hit = $nocall_ins_num{$read} if (exists $nocall_ins_num{$read});
        my $recall = int ($match_ins_num{$read} / ($ref_ins_num + $ins_dup_hit) * 1000) / 10;
        print "$recall\t";
        print OUT "$recall\t";
    }
    print "\nPrecis (A)\t";
    print OUT "\nPrecis (A)\t";
    foreach my $read (sort {$a <=> $b} keys %match_ins_num){
        my $precis = 0;
        my $call = $call_ins_num{$read};
        $precis = int ($match_ins_num{$read} / $call * 1000) / 10 if ($call > 0);
        print "$precis\t";
        print OUT "$precis\t";
    }
    if ($parent_flag == 1){
        print "\nMIER   (A)\t";
        print OUT "\nMIER (A)\t";
        foreach my $read (sort {$a <=> $b} keys %{${$mendel_err{'INS'}}{'A'}}){
            my $error = ${${$mendel_err{'INS'}}{'A'}}{$read};
            my $error_rate = 0;
            $error_rate = int ($error / $call_ins_num{$read} * 1000) / 10 if ($call_ins_num{$read} > 0);
            print "$error_rate\t";
            print OUT "$error_rate\t";
        }
    }
    if (($eval_gt == 1) and ($gt_flag == 1)){
	print "\nGenotying (A)\t";
	print OUT "\nGenotying (A)\t";
	foreach my $read (sort {$a <=> $b} keys %{${$match_GT{$type}}{'A'}}){
            my $matchGT = 0;
	    $matchGT = ${${$match_GT{$type}}{'A'}}{$read} if (exists ${${$match_GT{$type}}{'A'}}{$read});
	    my $match = ${$match_dup_num{'A'}}{$read};
            my $NA = 0;
	    $NA = ${${$match_GTNA{$type}}{'A'}}{$read} if (exists ${${$match_GTNA{$type}}{'A'}}{$read});
	    my $precis = 0;
	    $precis = int ($matchGT / ($match - $NA) * 1000) / 10 if ($match - $NA > 0);
	    print "$precis($match $matchGT $NA)\t";
	    print OUT "$precis($match $matchGT $NA)\t";
	}
    }
    if ($eval_bp == 1){
        print "\n";
        print OUT "\n";
	my $SDbp = 0;
	my $SDlen = 0;
	if (exists ${$bp_diff{$type}}{'A'}){
	    my $sq_total = 0;
	    foreach (@{${$bp_diff{$type}}{'A'}}){
		$sq_total += $_ ** 2;
	    }
	    $SDbp = int (($sq_total / @{${$bp_diff{$type}}{'A'}}) ** 0.5 + 0.5) if (@{${$bp_diff{$type}}{'A'}} > 0);
	}
	if (exists ${$len_diff{$type}}{'A'}){
	    my $sq_total = 0;
	    foreach (@{${$len_diff{$type}}{'A'}}){
		$sq_total += $_ ** 2;
	    }
	    $SDlen = int (($sq_total / @{${$len_diff{$type}}{'A'}}) ** 0.5 + 0.5) if (@{${$len_diff{$type}}{'A'}} > 0);
	}
	if (exists ${$bp_diff{$type}}{'A'}){
	    print "[Break Points SD ('A')]: $SDbp\t[Size SD ('A')]: $SDlen\n";
	    print OUT "[Break Points SD ('A')]: $SDbp\t[Size SD ('A')]: $SDlen\n";
	}
    }
    print "\n\n";
    print OUT "\n\n";
}
if ((scalar keys %call_mei_num > 0) and (($ref_sv eq $ref_mei) or ($ref_sv eq $ref_vei) or ($ref_sv eq $ref_numt))){
    foreach my $mei (@MEI){
        next if (!exists $ref_mei_num{$mei});
	foreach my $minread (sort {$a <=> $b} @min_reads){
	    ${$call_mei_num{$mei}}{$minread} = 0 if (!exists ${$call_mei_num{$mei}}{$minread});
	    ${$match_mei_num{$mei}}{$minread} = 0 if (!exists ${$match_mei_num{$mei}}{$minread});
	}
    }
    foreach my $mei (@MEI){
	next if (!exists $ref_mei_num{$mei});
	print "\n# $mei\n";
	print "Call (A)\t";
	print OUT "Call (A)\t";
	foreach my $minread (sort {$a <=> $b} keys %{$call_mei_num{$mei}}){
	    print "${$call_mei_num{$mei}}{$minread}\t";
	    print OUT "${$call_mei_num{$mei}}{$minread}\t";
	}
	print "\nRecall (A)\t";
	print OUT "\nRecall (A)\t";
	foreach my $read (sort {$a <=> $b} keys %{$match_mei_num{$mei}}){
	    my $match = ${$match_mei_num{$mei}}{$read};
	    my $recall = int ($match / $ref_mei_num{$mei} * 1000) / 10;
	    print "$recall\t";
	    print OUT "$recall\t";
	}
	print "\nPrecis (A)\t";
	print OUT "\nPrecis (A)\t";
	foreach my $read (sort {$a <=> $b} keys %{$match_mei_num{$mei}}){
	    my $match = ${$match_mei_num{$mei}}{$read};
	    my $precis = 0;
	    $precis = int ($match / ${$call_mei_num{$mei}}{$read} * 1000) / 10 if (${$call_mei_num{$mei}}{$read} > 0);
	    print "$precis\t";
	    print OUT "$precis\t";
	}
	if (($eval_gt == 1) and ($gt_flag == 1)){
	    print "\nGenotying (A)\t";
	    print OUT "\nGenotying (A)\t";
	    foreach my $read (sort {$a <=> $b} keys %{${$match_GT{$mei}}{'A'}}){
		my $matchGT = 0;
		$matchGT = ${${$match_GT{mei}}{'A'}}{$read} if (exists ${${$match_GT{$mei}}{'A'}}{$read});
		my $match = ${$match_mei_num{'A'}}{$read};
		my $NA = 0;
		$NA = ${${$match_GTNA{$mei}}{'A'}}{$read} if (exists ${${$match_GTNA{$mei}}{'A'}}{$read});
		my $precis = 0;
		$precis = int ($matchGT / ($match - $NA) * 1000) / 10 if ($match - $NA > 0);
		print "$precis($match $matchGT $NA)\t";
		print OUT "$precis($match $matchGT $NA)\t";
	    }
	}
    }
    print "\n\n";
    print OUT "\n\n";
}

if ((exists $call_inv_num{'A'}) and ($ref_inv_num > 0)){
     my $type = 'INV';
    print "## INV ##\n";
    print "        \t<Number of supporting reads>\n";
    print "        ";
    print OUT "## INV ##\n";
    print OUT "        \t<Number of supporting reads>\n";
    print OUT "        ";
    foreach my $minread (sort {$a <=> $b} @min_reads){
        print "\t$minread";
        print OUT "\t$minread";
        if (!exists ${$call_inv_num{'A'}}{$minread}){
            ${$call_inv_num{'A'}}{$minread} = 0;
            ${${$mendel_err{'INV'}}{'A'}}{$minread} = 0;
        }
        if (!exists ${$call_inv_num{'S'}}{$minread}){
            ${$call_inv_num{'S'}}{$minread} = 0;
        }
        if (!exists ${$call_inv_num{'M'}}{$minread}){
            ${$call_inv_num{'M'}}{$minread} = 0;
        }
        if (!exists ${$call_inv_num{'L'}}{$minread}){
            ${$call_inv_num{'L'}}{$minread} = 0;
        }
        if (!exists ${$match_inv_num{'A'}}{$minread}){
            ${$match_inv_num{'A'}}{$minread} = 0;
        }
        if (!exists ${$match_inv_num{'S'}}{$minread}){
            ${$match_inv_num{'S'}}{$minread} = 0;
        }
        if (!exists ${$match_inv_num{'M'}}{$minread}){
            ${$match_inv_num{'M'}}{$minread} = 0;
        }
        if (!exists ${$match_inv_num{'L'}}{$minread}){
            ${$match_inv_num{'L'}}{$minread} = 0;
        }
    }
    print "\nCall   (A)\t";
    print OUT "\nCall (A)\t";
    foreach my $read (sort {$a <=> $b} keys %{$call_inv_num{'A'}}){
        print "${$call_inv_num{'A'}}{$read}\t";
        print OUT "${$call_inv_num{'A'}}{$read}\t";
    }
    print "\nRecall (A)\t";
    print OUT "\nRecall (A)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_inv_num{'A'}}){
        my $recall = int (${$match_inv_num{'A'}}{$read} / $ref_inv_num * 1000) / 10;
        print "$recall\t";
        print OUT "$recall\t";
    }
    print "\nPrecis (A)\t";
    print OUT "\nPrecis (A)\t";
    foreach my $read (sort {$a <=> $b} keys %{$match_inv_num{'A'}}){
        my $precis = 0;
        $precis = int (${$match_inv_num{'A'}}{$read} / ${$call_inv_num{'A'}}{$read} * 1000) / 10 if (${$call_inv_num{'A'}}{$read} > 0);
        print "$precis\t";
        print OUT "$precis\t";
    }
    if ($parent_flag == 1){
        print "\nMIER   (A)\t";
        print OUT "\nMIER (A)\t";
        foreach my $read (sort {$a <=> $b} keys %{${$mendel_err{'INV'}}{'A'}}){
            my $error = ${${$mendel_err{'INV'}}{'A'}}{$read};
            my $error_rate = 0;
            $error_rate = int ($error / ${$call_inv_num{'A'}}{$read} * 1000) / 10 if (${$call_inv_num{'A'}}{$read} > 0);
            print "$error_rate\t";
            print OUT "$error_rate\t";
        }
    }
    if (($eval_gt == 1) and ($gt_flag == 1)){
	print "\nGenotying (A)\t";
	print OUT "\nGenotying (A)\t";
	foreach my $read (sort {$a <=> $b} keys %{${$match_GT{$type}}{'A'}}){
            my $matchGT = 0;
	    $matchGT = ${${$match_GT{$type}}{'A'}}{$read} if (exists ${${$match_GT{$type}}{'A'}}{$read});
	    my $match = ${$match_inv_num{'A'}}{$read};
            my $NA = 0;
	    $NA = ${${$match_GTNA{$type}}{'A'}}{$read} if (exists ${${$match_GTNA{$type}}{'A'}}{$read});
	    my $precis = 0;
	    $precis = int ($matchGT / ($match - $NA) * 1000) / 10 if ($match - $NA > 0);
	    print "$precis($match $matchGT $NA)\t";
	    print OUT "$precis($match $matchGT $NA)\t";
	}
    }
    if ($eval_bp == 1){
        print "\n";
        print OUT "\n";
	my $SDbp = 0;
	my $SDlen = 0;
	if (exists ${$bp_diff{$type}}{'A'}){
	    my $sq_total = 0;
	    foreach (@{${$bp_diff{$type}}{'A'}}){
		$sq_total += $_ ** 2;
	    }
	    $SDbp = int (($sq_total / @{${$bp_diff{$type}}{'A'}}) ** 0.5 + 0.5) if (@{${$bp_diff{$type}}{'A'}} > 0);
	}
	if (exists ${$len_diff{$type}}{'A'}){
	    my $sq_total = 0;
	    foreach (@{${$len_diff{$type}}{'A'}}){
		$sq_total += $_ ** 2;
	    }
	    $SDlen = int (($sq_total / @{${$len_diff{$type}}{'A'}}) ** 0.5 + 0.5) if (@{${$len_diff{$type}}{'A'}} > 0);
	}
	if (exists ${$bp_diff{$type}}{'A'}){
	    print "[Break Points SD ('A')]: $SDbp\t[Size SD ('A')]: $SDlen\n";
	    print OUT "[Break Points SD ('A')]: $SDbp\t[Size SD ('A')]: $SDlen\n";
	}
    }
    print "\n\n";
    print OUT "\n\n";
}
close (OUT);

if ($out_TF > 0){
    my $out_file = "$var_base.TP.vcf" if ($out_TF == 1);
    $out_file = "$var_base.FP.vcf" if ($out_TF == 2);
    $out_file = "$var_base.TF.vcf" if ($out_TF == 3);
    open (OUT1, "> $out_file");
    open (FILE, $var_child) or die "$var_child is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        if ($line =~ /^#|^$/){
            next;
        }
        my ($chr, $pos) = split (/\t/, $line);
	my $type = $1 if ($line =~ /SVTYPE=(.+?);/);
        $type = 'INS' if ($type =~ /MEI|VEI|NUMT/);
        if (exists ${${$TP{$chr}}{$pos}}{$type}){
            print OUT1 "$line\n" if ($out_TF == 1);
            print OUT1 "TP $line\n" if ($out_TF == 3);
        }
        else{
            print OUT1 "$line\n" if ($out_TF == 2);
            print OUT1 "FP $line\n" if ($out_TF == 3);
        }
    }
    close (FILE);
    close (OUT1);
}
