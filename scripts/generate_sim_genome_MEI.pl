#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw/$Bin $Script/;

# generated simulated diploid genome and vcf file (sim-genome coordinates) , where AluY, LINE1, SVA, and HERVK mobile elements were introduced
# Making a simulated reference genome where user-defined rates of artificial SNPs and 1- to 6-bp indels are randomly introduced.
# The SNP bases and the indel frequency depending on the length are automatically selected, based on a naturally occuring SNP compoition and indel frequency.
# excluding regions of predefined mobile elements (simulated ME insertions into the annotated ME regions are not allowed)
# introducing MEs are from libraries of similar sequences of MEs (blast search for chr1; min_ident: 90%, min_qcov: 10% for L1 and HERVK, 50% for SVA, and 70% for AluY)

# AluY:L1:SVA:HERVK = 15:3.5:1:1 (300:70:20:20 for chr17)

my $script_path = $Bin;

my $reference = '';		# -r
my $output_prefix = 'out';	# -p
my $Alu_fasta = "$Bin/../Simulated_data/Alu_chr1_mi90mc10.fa";            # -af
my $L1_fasta = "$Bin/../Simulated_data/L1_chr1_mi90mc10.fa";              # -lf
my $SVA_fasta = "$Bin/../Simulated_data/SVA_chr1_mi90mc10.fa";            # -sf
my $HERVK_fasta = "$Bin/../Simulated_data/HERVK_chr1_mi90mc10.fa";        # -hf
my $ME_bed = "$Bin/../Simulated_data/hs37_chr17_rmsk_MEI.bed";            # -mb	for excluding regions of ME insertion (simulated ME insertions into the annotated ME regions are not allowed)
my $snp_rate = 0.001;		# -s
my $indel_rate = 0.0002;	# -i
my $Alu_num = 500;              # -an
my $L1_num = 100;               # -ln
my $SVA_num = 40;               # -sn
my $HERVK_num = 40;             # -hn
my $hetero_rate = 2;            # -ht
my $help;			# -h

GetOptions(
    'ref|r=s' => \$reference,
    'pref|p=s' => \$output_prefix,
    'al_fa|af=s' => \$Alu_fasta,
    'l1_fa|lf=s' => \$L1_fasta,
    'sv_fa|sf=s' => \$SVA_fasta,
    'he_fa|hf=s' => \$HERVK_fasta,
    'me_bed|mb=s' => \$ME_bed,
    'snp=f' => \$snp_rate,
    'indel|i=f' => \$indel_rate,
    'alu_num|an=i' => \$Alu_num,
    'l1_num|ln=i' => \$L1_num,
    'sva_num|sn=i' => \$SVA_num,
    'her_num|hn=i' => \$HERVK_num,
    'het_rate|ht=f' => \$hetero_rate,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;
pod2usage(-verbose => 0) unless (-e $reference);


print STDERR "########## ME simulate options ##########\n";
print STDERR "AluY=$Alu_num L1=$L1_num SVA=$SVA_num HERVK=$HERVK_num snp-rate=$snp_rate indel-rate=$indel_rate het-hm-rate: $hetero_rate pref=$output_prefix ref=$reference\n";

=head1 SYNOPSIS

#### generate_sim_genome_MEI.pl [options] -r <reference file> -p <prefix of output files>
  Output:
   prefix.maternal.fa	      fasta file for maternal simulated reference
   prefix.paternal.fa	      fasta file for paternal simulated reference
   prefix.ME.vcf	          vcf file showing the reference positions and ME species of introduced mobile elements
   prefix.snv.vcf	          vcf file showing the reference positions and bases of introduced SNPs and indels

  Options:
   --ref or -r <STR>          reference fasta file
   --pref or -p <STR>         prefix of outputfiles [default: out]
   --al_fa or -af <STR>       fasta file of blast-searched Alu library to be introduced [default: Alu_chr1_mi90mc10.fa]
   --l1_fa or -lf <STR>       fasta file of blast-searched L1 library to be introduced [default: L1_chr1_mi90mc10.fa]
   --sv_fa or -sf <STR>       fasta file of blast-searched SVA library to be introduced [default: SVA_chr1_mi90mc10.fa]
   --he_fa or -hf <STR>       fasta file of blast-searched HERVK library to be introduced [default: HERVK_chr1_mi90mc10.fa]
   --me_bed or -mb <STR>      bed file of mobile elements [i.e., Alu, L1, SVA, and HERVK] defined for the hs37 reference (optional) [default: hs37_chr17_rmsk_MEI.bed]
   --snp or -s <FLOAT>        genomic mutation rate for SNPs [default: 0.01]
   --indel or -i <FLOAT>      genomic mutation rate for indels [default: 0.002]
   --alu_num or -an <INT>     number of introducing Alu [default: 300]
   --l1_num or -ln <INT>      number of introducing LINE1 [default: 70]
   --sva_num or -sn <INT>     number of introducing SVA [default: 20]
   --her_num or -hn <INT>     number of introducing HERVK [default: 20]
   --het_rate or -ht <FLOAT>  hetero to homo ratio of introducing variants [default: 2]
   --help or -h               output help message

=cut

die "reference file is not specified: \n" if ($reference eq '');

my @ME_file;
push @ME_file, $Alu_fasta, $L1_fasta, $SVA_fasta, $HERVK_fasta;

my $genome_size = 0;
my %chr_pat;
my %chr_mat;
my %chr_length;
my %chr_rate;
my @chr_name;

my $seq = '';
my $header = '';

open (FILE, $reference) or die "$reference is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^>(\S+)/){
	if ($seq ne ''){
	    $chr_length{$header} = length $seq;
	    $genome_size += length $seq;
	    if (length $seq > 1000000){
		my $max_1M = int (length ($seq) / 1000000);
		$max_1M ++ if (length ($seq) % 1000000 > 0);
		for (my $i = 1; $i <= $max_1M; $i++){
		    if (length $seq > 1000000){
			my $subseq = substr ($seq, 0, 1000000, '');
			${$chr_pat{$header}}{$i} = $subseq;
			${$chr_mat{$header}}{$i} = $subseq;
		    }
		    else{
			${$chr_pat{$header}}{$i} = $seq;
			${$chr_mat{$header}}{$i} = $seq;
		    }
		}
	    }
	    else{
		${$chr_pat{$header}}{1} = $seq;
		${$chr_mat{$header}}{1} = $seq;
	    }
	}
	$header = $1;
	push @chr_name, $header;
	$seq = '';
    }
    else{
	$seq .= uc $line;
    }
}
$chr_length{$header} = length $seq;
$genome_size += length $seq;
if (length $seq > 1000000){
    my $max_1M = int (length ($seq) / 1000000);
    $max_1M ++ if (length ($seq) % 1000000 > 0);
    for (my $i = 1; $i <= $max_1M; $i++){
	if (length $seq > 1000000){
	    my $subseq = substr ($seq, 0, 1000000, '');
	    ${$chr_pat{$header}}{$i} = $subseq;
	    ${$chr_mat{$header}}{$i} = $subseq;
	}
	else{
	    ${$chr_pat{$header}}{$i} = $seq;
	    ${$chr_mat{$header}}{$i} = $seq;
	}
    }
}
else{
    ${$chr_pat{$header}}{1} = $seq;
    ${$chr_mat{$header}}{1} = $seq;
}
close (FILE);

foreach (keys %chr_length){
    $chr_rate{$_} = int (($chr_length{$_} / $genome_size) * 1000) / 1000;
}


my %ME_list;
my %ME_num;
my %exclude;

foreach my $file (@ME_file){
    $seq = '';
    $header = '';
print STDERR "$file\n";
    open (FILE, $file) or die "$file is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^>(\S+)/){
	    if ($seq ne ''){
		push @{$ME_list{$header}}, $seq;
		$ME_num{$header} ++;
	    }
	    $seq = '';
	    $header = $1;
	    $header = $1 if ($header =~ /(.+)\./);
	    $header =~ s/-// if ($header =~ /-/);
	}
	else{
	    $seq .= uc $line;
	}
    }
    if ($seq ne ''){
	push @{$ME_list{$header}}, $seq;
	$ME_num{$header} ++;
    }
    close (FILE);
}

if ($ME_bed ne ''){
    open (FILE, "$ME_bed") or die "$ME_bed is not found: $!\n";
    while (my $line = <FILE>){
	chomp $line;
	next if ($line =~ /^#/);
	my @line = split (/\s+/, $line);
	my $chr = $line[0];
	my $start = $line[1];
	my $end = $line[2];
	${$exclude{$chr}}{$start} = $end;
    }
    close (FILE);
}

foreach my $tag (keys %ME_list){
    print STDERR "$tag\t", $ME_num{$tag}, "\t";
    print STDERR scalar @{$ME_list{$tag}}, "\n";
}

my %indel_all;
my %indel_allpos;
my @indel_info;
my @indel_orig_info;
my %count_var;
my $count_var = 0;
my $count_snp_2 = 0;
my $count_overlength = 0;

my $hetero_rate2 = $hetero_rate / (1 + $hetero_rate) / 2;

my $snp_no = $snp_rate * $genome_size;
my $snp_ht_no = int ($snp_no * $hetero_rate2);
my $snp_hm_no = $snp_no - $snp_ht_no * 2;
my $indel_no = $indel_rate * $genome_size;
my $indel_ht_no = $indel_no * $hetero_rate2;
my $indel_hm_no = $indel_no - $indel_ht_no * 2;


my $alu_ht_num = int ($Alu_num * $hetero_rate2);
my $alu_hm_num = $Alu_num - $alu_ht_num * 2;
my $l1_ht_num = int ($L1_num * $hetero_rate2);
my $l1_hm_num = $L1_num - $l1_ht_num * 2;
my $sva_ht_num = int ($SVA_num * $hetero_rate2);
my $sva_hm_num = $SVA_num - $sva_ht_num * 2;
my $her_ht_num = int ($HERVK_num * $hetero_rate2);
my $her_hm_num = $HERVK_num - $her_ht_num * 2;

# the expected distribution of indel 1-6 bp is derived from the japanese human genome (Nature Genet 42, 931 2010 Suppl.Fig.13)
my $indel_1b = $indel_no * 0.5 * 0.66;
my $indel_2b = $indel_no * 0.5 * 0.17;
my $indel_3b = $indel_no * 0.5 * 0.07;
my $indel_4b = $indel_no * 0.5 * 0.07;
my $indel_5b = $indel_no * 0.5 * 0.02;
my $indel_6b = $indel_no * 0.5 * 0.01;

my $indel_1b_ht = int ($indel_1b * $hetero_rate2);
my $indel_1b_hm = $indel_1b - $indel_1b_ht * 2;
my $indel_2b_ht = int ($indel_2b * $hetero_rate2);
my $indel_2b_hm = $indel_2b - $indel_2b_ht * 2;
my $indel_3b_ht = int ($indel_3b * $hetero_rate2);
my $indel_3b_hm = $indel_3b - $indel_3b_ht * 2;
my $indel_4b_ht = int ($indel_4b * $hetero_rate2);
my $indel_4b_hm = $indel_4b - $indel_4b_ht * 2;
my $indel_5b_ht = int ($indel_5b * $hetero_rate2);
my $indel_5b_hm = $indel_5b - $indel_5b_ht * 2;
my $indel_6b_ht = int ($indel_6b * $hetero_rate2);
my $indel_6b_hm = $indel_6b - $indel_6b_ht * 2;


my @nuc = ('A', 'C', 'G', 'T');
my @snp_nuc_A = ('G', 'G', 'G', 'G', 'C', 'T');
my @snp_nuc_G = ('A', 'A', 'A', 'A', 'C', 'T');
my @snp_nuc_C = ('T', 'T', 'T', 'T', 'A', 'G');
my @snp_nuc_T = ('C', 'C', 'C', 'C', 'A', 'G');

my @strand = ('Plus', 'Minus');

my $mei_seq = "$output_prefix.mei.fa";

&make_indel ($indel_1b_ht, $indel_1b_hm, $indel_2b_ht, $indel_2b_hm, $indel_3b_ht, $indel_3b_hm, $indel_4b_ht, $indel_4b_hm, $indel_5b_ht, $indel_5b_hm, $indel_6b_ht, $indel_6b_hm);

&make_snp ($snp_ht_no, $snp_hm_no);

&make_me ('AluY', $alu_ht_num, $alu_hm_num);

&make_me ('LINE1', $l1_ht_num, $l1_hm_num);

&make_me ('SVA', $sva_ht_num, $sva_hm_num);

&make_me ('HERVK', $her_ht_num, $her_hm_num);

&subst_genome ();


print STDERR 'total AluY: ', $count_var{'AluY'}, "\n";
print STDERR 'total LINE1: ', $count_var{'LINE1'}, "\n";
print STDERR 'total SVA: ', $count_var{'SVA'}, "\n";
print STDERR 'total HERVK: ', $count_var{'HERVK'}, "\n";
print STDERR 'total SNPs = ', $count_var{'SNP'}, "\n";
print STDERR 'total INDELs = ', $count_var{'INDEL'}, "\n";
print STDERR 'total susbstituted variants: ', $count_var, "\n";
print STDERR 'number of overlength: ', $count_overlength, "\n";

my $out_mat_fasta = $output_prefix . '.' . 'maternal.fa';
my $out_pat_fasta = $output_prefix . '.' . 'paternal.fa';
my $out_me = $output_prefix . '.' . 'ME.vcf';
my $out_snv = $output_prefix . '.' . 'snv.vcf';

open (OUT1, "> $out_me");
open (OUT2, "> $out_snv");
print OUT1 "#AluY:$Alu_num,LINE1:$L1_num,SVA:$SVA_num,SNP-rate:$snp_rate,Indel-rate:$indel_rate,Hetero-rate:$hetero_rate\n";
print OUT1 "#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
print OUT2 "#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
foreach my $chr (@chr_name){
    my $sum_shift_pos_pat = 0;
    my $sum_shift_pos_mat = 0;
    if (exists $indel_all{$chr}){
	foreach my $pos (sort {$a <=> $b} keys %{$indel_all{$chr}}){
	    my ($info, $parent) = split (/=/, ${$indel_all{$chr}}{$pos});
	    my ($tag, $ref, $alt, $strand) = split (/-/, $info);
	    my $tag2 = $tag;
	    my $alt_len = length $alt;
	    $parent = 'Pat/Mat' if ($parent eq 'HM');
	    my $alt_pos = $pos + $sum_shift_pos_pat if ($parent eq 'Pat');
	    $alt_pos = $pos + $sum_shift_pos_mat if ($parent eq 'Mat');
	    $alt_pos = 'mat:' . ($pos + $sum_shift_pos_mat) . ',pat:' . ($pos + $sum_shift_pos_pat) if ($parent eq 'Pat/Mat');
	    if ($tag eq 'SNP'){
		print OUT2 "$chr\t$pos\t$tag2\t$ref\t$alt\t.\tPASS\tSVTYPE=$tag;SVLEN=$alt_len;ALTPOS=$alt_pos;PARENT=$parent\n";
	    }
	    elsif ($tag eq 'I'){
		if ($parent eq 'Pat/Mat'){
		    $sum_shift_pos_pat += length $alt;
		    $sum_shift_pos_mat += length $alt;
		}
		elsif ($parent eq 'Pat'){
		    $sum_shift_pos_pat += length $alt;
		}
		elsif ($parent eq 'Mat'){
		    $sum_shift_pos_mat += length $alt;
		}
		print OUT2 "$chr\t$pos\t$tag2\t$ref\t$alt\t.\tPASS\tSVTYPE=$tag;SVLEN=$alt_len;ALTPOS=$alt_pos;PARENT=$parent\n";
	    }
	    elsif ($tag eq 'D'){
		if ($parent eq 'Pat/Mat'){
		    $sum_shift_pos_pat -= length $alt;
		    $sum_shift_pos_mat -= length $alt;
		}
		elsif ($parent eq 'Pat'){
		    $sum_shift_pos_pat -= length $alt;
		}
		elsif ($parent eq 'Mat'){
		    $sum_shift_pos_mat -= length $alt;
		}
		print OUT2 "$chr\t$pos\t$tag2\t$ref\t$alt\t.\tPASS\tSVTYPE=$tag;SVLEN=$alt_len;ALTPOS=$alt_pos;PARENT=$parent\n";
	    }
	    else{
		if ($parent eq 'Pat/Mat'){
		    $sum_shift_pos_pat += $alt;
		    $sum_shift_pos_mat += $alt;
		}
		elsif ($parent eq 'Pat'){
		    $sum_shift_pos_pat += $alt;
		}
		elsif ($parent eq 'Mat'){
		    $sum_shift_pos_mat += $alt;
		}
		$tag2 = 'Alu' if ($tag =~ /^Alu/);
		$tag2 = 'LINE1' if ($tag =~ /^LINE/);
		$tag2 = 'SVA' if ($tag =~ /^SVA/);
		$tag2 = 'HERVK' if ($tag =~ /^HERVK/);
		$alt_len = $alt;
		$alt = '.';
		print OUT1 "$chr\t$pos\t$tag2\t$ref\t$alt\t.\tPASS\tSVTYPE=MEI;SVLEN=$alt_len;SVDIR=$strand;ALTPOS=$alt_pos;PARENT=$parent\n";
	    }
	}
    }
}
close (OUT1);
close (OUT2);

my $chr_seq = '';
open (OUT1, "> $out_pat_fasta");
foreach my $chrom (@chr_name){
    if (exists $chr_pat{$chrom}){
	foreach my $chr_frac (sort {$a <=> $b} keys %{$chr_pat{$chrom}}){
	    $chr_seq .= ${$chr_pat{$chrom}}{$chr_frac};
	}
	print OUT1 ">$chrom-paternal\n";
	while (length $chr_seq > 0){
	    if (length $chr_seq > 60){
		my $subseq = substr ($chr_seq, 0, 60, '');
		print OUT1 $subseq, "\n";
	    }
	    else{
		print OUT1 $chr_seq, "\n";
		$chr_seq = '';
	    }
	}
    }
}
close (OUT1);

$chr_seq = '';
open (OUT1, "> $out_mat_fasta");
foreach my $chrom (@chr_name){
    if (exists $chr_mat{$chrom}){
	foreach my $chr_frac (sort {$a <=> $b} keys %{$chr_mat{$chrom}}){
	    $chr_seq .= ${$chr_mat{$chrom}}{$chr_frac};
	}
	print OUT1 ">$chrom-maternal\n";
	while (length $chr_seq > 0){
	    if (length $chr_seq > 60){
		my $subseq = substr ($chr_seq, 0, 60, '');
		print OUT1 $subseq, "\n";
	    }
	    else{
		print OUT1 $chr_seq, "\n";
		$chr_seq = '';
	    }
	}
    }
}
close (OUT1);

###########################################################################################################

sub make_me{
    my ($species, $ht_num, $hm_num) = @_;
    foreach my $chrno (sort keys %chr_rate){			# assign the positions of homozygous mobile elements
	for (my $i = 0; $i < int ($hm_num * $chr_rate{$chrno} + 0.5); $i++){
	    my $pos1 = int (rand ($chr_length{$chrno})) + 1;
	    if (exists ${$indel_allpos{$chrno}}{$pos1}){
		$i = $i - 1;
		next;
	    }
	    my $chr_divnum = 0;
	    $chr_divnum = int ($pos1 / 1000000);
	    my $pos_2 = $pos1 % 1000000;
	    $pos_2 = 1000000 if ($pos_2 == 0);
	    $chr_divnum ++ if ($pos1 % 1000000 > 0);
	    my $ref_base = substr (${$chr_pat{$chrno}}{$chr_divnum}, $pos_2 - 1, 1);
	    if ($ref_base =~ /[^ACGT]/){
		$i = $i - 1;
		next;
	    }
	    if (exists $exclude{$chrno}){
		my $hit_flag = 0;
		foreach my $start (sort {$a <=> $b} keys %{$exclude{$chrno}}){
		    next if ($pos1 < $start);
		    last if ($pos1 > $start + 10000);
		    if ($pos1 >= $start){
			if ($pos1 <= ${$exclude{$chrno}}{$start}){
			    $hit_flag = 1;
			    last;
			}
		    }
		}
		if ($hit_flag == 1){
		    $i = $i - 1;
		    next;
		}
	    }
	    ${$indel_all{$chrno}}{$pos1} = "$species=HM";
	    for (my $i = $pos1 - 10; $i <= $pos1 + 10; $i++){
		${$indel_allpos{$chrno}}{$i} = $species;
	    }
	}
    }
    foreach my $chrno (sort keys %chr_rate){			# assign the positions of heterozygous paternal mobile elements
	for (my $i = 0; $i < int ($ht_num * $chr_rate{$chrno} + 0.5); $i++){
	    my $pos1 = int (rand ($chr_length{$chrno})) + 1;
	    if (exists ${$indel_allpos{$chrno}}{$pos1}){
		$i = $i - 1;
		next;
	    }
	    my $chr_divnum = 0;
	    $chr_divnum = int ($pos1 / 1000000);
	    my $pos_2 = $pos1 % 1000000;
	    $pos_2 = 1000000 if ($pos_2 == 0);
	    $chr_divnum ++ if ($pos1 % 1000000 > 0);
	    my $ref_base = substr (${$chr_pat{$chrno}}{$chr_divnum}, $pos_2 - 1, 1);
	    if ($ref_base =~ /[^ACGT]/){
		$i = $i - 1;
		next;
	    }
	    if (exists $exclude{$chrno}){
		my $hit_flag = 0;
		foreach my $start (sort {$a <=> $b} keys %{$exclude{$chrno}}){
		    next if ($pos1 < $start);
		    last if ($pos1 > $start + 10000);
		    if ($pos1 >= $start){
			if ($pos1 <= ${$exclude{$chrno}}{$start}){
			    $hit_flag = 1;
			    last;
			}
		    }
		}
		if ($hit_flag == 1){
		    $i = $i - 1;
		    next;
		}
	    }
	    ${$indel_all{$chrno}}{$pos1} = "$species=Pat";
	    for (my $i = $pos1 - 10; $i <= $pos1 + 10; $i++){
		${$indel_allpos{$chrno}}{$i} = $species;
	    }
	}
    }
    foreach my $chrno (sort keys %chr_rate){			# assign the positions of heterozygous maternal mobile elements
	for (my $i = 0; $i < int ($ht_num * $chr_rate{$chrno} + 0.5); $i++){
	    my $pos1 = int (rand ($chr_length{$chrno})) + 1;
	    if (exists ${$indel_allpos{$chrno}}{$pos1}){
		$i = $i - 1;
		next;
	    }
	    my $chr_divnum = 0;
	    $chr_divnum = int ($pos1 / 1000000);
	    my $pos_2 = $pos1 % 1000000;
	    $pos_2 = 1000000 if ($pos_2 == 0);
	    $chr_divnum ++ if ($pos1 % 1000000 > 0);
	    my $ref_base = substr (${$chr_pat{$chrno}}{$chr_divnum}, $pos_2 - 1, 1);
	    if ($ref_base =~ /[^ACGT]/){
		$i = $i - 1;
		next;
	    }
	    if (exists $exclude{$chrno}){
		my $hit_flag = 0;
		foreach my $start (sort {$a <=> $b} keys %{$exclude{$chrno}}){
		    next if ($pos1 < $start);
		    last if ($pos1 > $start + 10000);
		    if ($pos1 >= $start){
			if ($pos1 <= ${$exclude{$chrno}}{$start}){
			    $hit_flag = 1;
			    last;
			}
		    }
		}
		if ($hit_flag == 1){
		    $i = $i - 1;
		    next;
		}
	    }
	    ${$indel_all{$chrno}}{$pos1} = "$species=Mat";
	    for (my $i = $pos1 - 10; $i <= $pos1 + 10; $i++){
		${$indel_allpos{$chrno}}{$i} = $species;
	    }
	}
    }
}

sub make_indel {
    my @indel_num = @_;
    my @indel_ht_num;
    my @indel_hm_num;
    my $count = 0;
    foreach (@indel_num){
	$count ++;
	if ($count % 2 == 1){
	    push @indel_ht_num, $_;
	}
	elsif ($count % 2 == 0){
	    push @indel_hm_num, $_;
	}
    }
    foreach my $chrno (sort keys %chr_rate){
	for (my $delno = 1; $delno < 7; $delno++){
	    for (my $i = 0; $i < int ($indel_hm_num[$delno - 1] * $chr_rate{$chrno} + 0.5); $i++){
		my $match_flag = 0;
		my $pos = int (rand ($chr_length{$chrno})) + 1;
		for (my $i = 0; $i <= $delno; $i++){
		    if (exists ${$indel_allpos{$chrno}}{$pos + $i}){
			$i = $i - 1;
			$match_flag = 1;
			last;
		    }
		}
		next if ($match_flag == 1);
		${$indel_all{$chrno}}{$pos} = 'D-' . $delno . '=HM';
		for (my $j = -1; $j <= $delno; $j++){ 
		    ${$indel_allpos{$chrno}}{$pos + $j} = 1;
		}
	    }
	    for (my $i = 0; $i < int ($indel_hm_num[$delno - 1] * $chr_rate{$chrno} + 0.5); $i++){
		my $pos = int (rand ($chr_length{$chrno})) + 1;
		my $ins_base = '';
		for (my $j = 0; $j < $delno; $j++){
		    $ins_base .= $nuc[(rand (4))];
		}
		if (exists ${$indel_allpos{$chrno}}{$pos}){
		    $i = $i - 1;
		    next;
		}
		${$indel_all{$chrno}}{$pos} = 'I-' . $ins_base . '=HM';
		for (my $j = -1; $j <= $delno; $j++){ 
		    ${$indel_allpos{$chrno}}{$pos + $j} = 1;
		}
	    }
	}
    }
    foreach my $chrno (sort keys %chr_rate){
	for (my $delno = 1; $delno < 7; $delno++){
	    for (my $i = 0; $i < int ($indel_ht_num[$delno - 1] * $chr_rate{$chrno} + 0.5); $i++){
		my $match_flag = 0;
		my $pos = int (rand ($chr_length{$chrno})) + 1;
		for (my $i = 0; $i <= $delno; $i++){
		    if (exists ${$indel_allpos{$chrno}}{$pos + $i}){
			$i = $i - 1;
			$match_flag = 1;
			last;
		    }
		}
		next if ($match_flag == 1);
		${$indel_all{$chrno}}{$pos} = 'D-' . $delno . '=Pat';
		for (my $j = -1; $j <= $delno; $j++){ 
		    ${$indel_allpos{$chrno}}{$pos + $j} = 1;
		}
	    }
	    for (my $i = 0; $i < int ($indel_ht_num[$delno - 1] * $chr_rate{$chrno} + 0.5); $i++){
		my $pos = int (rand ($chr_length{$chrno})) + 1;
		my $ins_base = '';
		for (my $j = 0; $j < $delno; $j++){
		    $ins_base .= $nuc[(rand (4))];
		}
		if (exists ${$indel_allpos{$chrno}}{$pos}){
		    $i = $i - 1;
		    next;
		}
		${$indel_all{$chrno}}{$pos} = 'I-' . $ins_base . '=Pat';
		for (my $j = -1; $j <= $delno; $j++){ 
		    ${$indel_allpos{$chrno}}{$pos + $j} = 1;
		}
	    }
	}
    }
    foreach my $chrno (sort keys %chr_rate){
	for (my $delno = 1; $delno < 7; $delno++){
	    for (my $i = 0; $i < int ($indel_ht_num[$delno - 1] * $chr_rate{$chrno} + 0.5); $i++){
		my $match_flag = 0;
		my $pos = int (rand ($chr_length{$chrno})) + 1;
		for (my $i = 0; $i <= $delno; $i++){
		    if (exists ${$indel_allpos{$chrno}}{$pos + $i}){
			$i = $i - 1;
			$match_flag = 1;
			last;
		    }
		}
		next if ($match_flag == 1);
		${$indel_all{$chrno}}{$pos} = 'D-' . $delno . '=Mat';
		for (my $j = -1; $j <= $delno; $j++){ 
		    ${$indel_allpos{$chrno}}{$pos + $j} = 1;
		}
	    }
	    for (my $i = 0; $i < int ($indel_ht_num[$delno - 1] * $chr_rate{$chrno} + 0.5); $i++){
		my $pos = int (rand ($chr_length{$chrno})) + 1;
		my $ins_base = '';
		for (my $j = 0; $j < $delno; $j++){
		    $ins_base .= $nuc[(rand (4))];
		}
		if (exists ${$indel_allpos{$chrno}}{$pos}){
		    $i = $i - 1;
		    next;
		}
		${$indel_all{$chrno}}{$pos} = 'I-' . $ins_base . '=Mat';
		for (my $j = -1; $j <= $delno; $j++){ 
		    ${$indel_allpos{$chrno}}{$pos + $j} = 1;
		}
	    }
	}
    }
}


sub make_snp {
    my ($ht_num, $hm_num) = @_;
    foreach my $chrno (sort keys %chr_rate){
	for (my $i = 0; $i < int ($hm_num * $chr_rate{$chrno} + 0.5); $i++){
	    my $pos = int (rand ($chr_length{$chrno})) + 1;
	    if (exists ${$indel_allpos{$chrno}}{$pos}){
		$i = $i - 1;
		next;
	    }
	    ${$indel_all{$chrno}}{$pos} = 'SNP=HM';
	    ${$indel_allpos{$chrno}}{$pos} = 1;
	}
    }
    foreach my $chrno (sort keys %chr_rate){
	for (my $i = 0; $i < int ($ht_num * $chr_rate{$chrno} + 0.5); $i++){
	    my $pos = int (rand ($chr_length{$chrno})) + 1;
	    if (exists ${$indel_allpos{$chrno}}{$pos}){
		$i = $i - 1;
		next;
	    }
	    ${$indel_all{$chrno}}{$pos} = 'SNP=Pat';
	    ${$indel_allpos{$chrno}}{$pos} = 1;
	}
    }
    foreach my $chrno (sort keys %chr_rate){
	for (my $i = 0; $i < int ($ht_num * $chr_rate{$chrno} + 0.5); $i++){
	    my $pos = int (rand ($chr_length{$chrno})) + 1;
	    if (exists ${$indel_allpos{$chrno}}{$pos}){
		$i = $i - 1;
		next;
	    }
	    ${$indel_all{$chrno}}{$pos} = 'SNP=Mat';
	    ${$indel_allpos{$chrno}}{$pos} = 1;
	}
    }
}


sub subst_genome{
    my $pre_count = 0;
    open (OUT1, "> $mei_seq");
    foreach my $chr (keys %indel_all){
	foreach my $pos (sort {$b <=> $a} keys %{$indel_all{$chr}}){
	    my ($tag, $parent) = split (/=/, ${$indel_all{$chr}}{$pos});
	    if ($chr_length{$chr} < $pos){
		delete ${$indel_all{$chr}}{$pos};
		next;
	    }
	    my $chr_divnum = 0;
	    $chr_divnum = int ($pos / 1000000);
	    my $pos_2 = $pos % 1000000;
	    $pos_2 = 1000000 if ($pos_2 == 0);
	    $chr_divnum ++ if ($pos % 1000000 > 0);
	    my $ref_base = substr (${$chr_pat{$chr}}{$chr_divnum}, $pos_2 - 1, 1);
	    if ($ref_base =~ /[^ACGT]/){
		delete ${$indel_all{$chr}}{$pos};
#print STDERR "# $chr\t$pos\t$pos_2\t$chr_divnum\t$ref_base\n";
		next;
	    }
	    $count_var++;
#print STDERR "$chr\t$pos\t$pos_2\t$chr_divnum\t$ref_base\n" if ($pre_count == $count_var);
print STDERR "substituted variants - $count_var\n" if ($count_var % 10000 == 0);
	    if ($tag eq 'SNP'){
		my $snp_base = $snp_nuc_A[rand (6)] if ($ref_base eq 'A');
		$snp_base = $snp_nuc_G[rand (6)] if ($ref_base eq 'G');
		$snp_base = $snp_nuc_C[rand (6)] if ($ref_base eq 'C');
		$snp_base = $snp_nuc_T[rand (6)] if ($ref_base eq 'T');
		if ($parent eq 'HM'){
		    my $ref_base_2 = substr (${$chr_pat{$chr}}{$chr_divnum}, $pos_2 - 1, 1, $snp_base);
		    $ref_base_2 = substr (${$chr_mat{$chr}}{$chr_divnum}, $pos_2 - 1, 1, $snp_base);
		}
		elsif ($parent eq 'Pat'){
		    my $ref_base_2 = substr (${$chr_pat{$chr}}{$chr_divnum}, $pos_2 - 1, 1, $snp_base);
		}
		elsif ($parent eq 'Mat'){
		    my $ref_base_2 = substr (${$chr_mat{$chr}}{$chr_divnum}, $pos_2 - 1, 1, $snp_base);
		}
		${$indel_all{$chr}}{$pos} = "SNP-$ref_base-$snp_base=$parent";
		$count_snp_2++;
		$count_var{$tag} ++;
	    }
	    elsif ($tag =~ /^I-(.+)/){
		my $ins_base = $1;
		if ($parent eq 'HM'){
		    my $ref_base_2 = substr (${$chr_pat{$chr}}{$chr_divnum}, $pos_2, 0, $ins_base);
		    $ref_base_2 = substr (${$chr_mat{$chr}}{$chr_divnum}, $pos_2, 0, $ins_base);
		}
		elsif ($parent eq 'Pat'){
		    my $ref_base_2 = substr (${$chr_pat{$chr}}{$chr_divnum}, $pos_2, 0, $ins_base);
		}
		elsif ($parent eq 'Mat'){
		    my $ref_base_2 = substr (${$chr_mat{$chr}}{$chr_divnum}, $pos_2, 0, $ins_base);
		}
		${$indel_all{$chr}}{$pos} = "I-$ref_base-$ins_base=$parent";
		$count_var{'INDEL'} ++;
	    }
	    elsif ($tag =~ /^D-(\d+)/){
		my $del_len = $1;
		if (($parent eq 'HM') or ($parent eq 'Pat')){
		    if ($pos_2 + $del_len <= length ${$chr_pat{$chr}}{$chr_divnum}){
			my $ref_base_2 = substr (${$chr_pat{$chr}}{$chr_divnum}, $pos_2, $del_len, '');
			${$indel_all{$chr}}{$pos} = "D-$ref_base-$ref_base_2=$parent";
		    }
		    else{
			if (!exists ${$chr_pat{$chr}}{$chr_divnum + 1}){
			    delete ${$indel_all{$chr}}{$pos};
			    $count_var --;
			    next;
			}
			if ($pos_2 + $del_len - length ${$chr_pat{$chr}}{$chr_divnum} > length ${$chr_pat{$chr}}{$chr_divnum + 1}){
			    delete ${$indel_all{$chr}}{$pos};
			    $count_var --;
			    next;
			}
			my $ref_base_2 = substr (${$chr_pat{$chr}}{$chr_divnum}, $pos_2, length (${$chr_pat{$chr}}{$chr_divnum}) - $pos_2, '');
			my $ref_base_3 = substr (${$chr_pat{$chr}}{$chr_divnum + 1}, 0, $del_len - length ($ref_base_2), '');
			my $del_seq = $ref_base_2 . $ref_base_3;
			${$indel_all{$chr}}{$pos} = "D-$ref_base-$del_seq=$parent";
		    }
		}
		if (($parent eq 'HM') or ($parent eq 'Mat')){
		    if ($pos_2 + $del_len <= length ${$chr_mat{$chr}}{$chr_divnum}){
			my $ref_base_2 = substr (${$chr_mat{$chr}}{$chr_divnum}, $pos_2, $del_len, '');
			${$indel_all{$chr}}{$pos} = "D-$ref_base-$ref_base_2=$parent";
		    }
		    else{
			if (!exists ${$chr_mat{$chr}}{$chr_divnum + 1}){
			    delete ${$indel_all{$chr}}{$pos};
			    $count_var --;
			    next;
			}
			if ($pos_2 + $del_len - length ${$chr_mat{$chr}}{$chr_divnum} > length ${$chr_mat{$chr}}{$chr_divnum + 1}){
			    delete ${$indel_all{$chr}}{$pos};
			    $count_var --;
			    next;
			}
			my $ref_base_2 = substr (${$chr_mat{$chr}}{$chr_divnum}, $pos_2, length (${$chr_mat{$chr}}{$chr_divnum}) - $pos_2, '');
			my $ref_base_3 = substr (${$chr_mat{$chr}}{$chr_divnum + 1}, 0, $del_len - length ($ref_base_2), '');
			my $del_seq = $ref_base_2 . $ref_base_3;
			${$indel_all{$chr}}{$pos} = "D-$ref_base-$del_seq=$parent";
		    }
		}
		$count_var{'INDEL'} ++;
	    }
	    else{
		my $me_num = $ME_num{$tag};
		my $me_seq = ${$ME_list{$tag}}[(rand ($me_num))];
		my $me_len = length $me_seq;
		my $strand = $strand[(rand (2))];
		if ($strand eq 'Minus'){
		    $me_seq = reverse $me_seq;
		    $me_seq =~ tr/ACGT/TGCA/;
		}
		if ($parent eq 'HM'){
		    my $ref_base_2 = substr (${$chr_pat{$chr}}{$chr_divnum}, $pos_2, 0, $me_seq);
		    $ref_base_2 = substr (${$chr_mat{$chr}}{$chr_divnum}, $pos_2, 0, $me_seq);
		}
		elsif ($parent eq 'Pat'){
		    my $ref_base_2 = substr (${$chr_pat{$chr}}{$chr_divnum}, $pos_2, 0, $me_seq);
		}
		elsif ($parent eq 'Mat'){
		    my $ref_base_2 = substr (${$chr_mat{$chr}}{$chr_divnum}, $pos_2, 0, $me_seq);
		}
		${$indel_all{$chr}}{$pos} = "$tag-$ref_base-$me_len-$strand=$parent";
		print OUT1 ">$tag-$chr-$pos $strand $parent\n";
		print OUT1 $me_seq, "\n";
		$count_var{$tag} ++;
	    }
	}
	$pre_count = $count_var;
    }
    close (OUT1);
}

