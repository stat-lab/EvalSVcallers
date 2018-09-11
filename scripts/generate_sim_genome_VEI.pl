#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw/$Bin $Script/;

# generated simulated diploid genome and vcf file (sim-genome coordinates) , where viral sequences with human-host were introduced.
# Making a simulated reference genome where user-defined rates of artificial SNPs and 1- to 6-bp indels are randomly introduced.
# The SNP bases and the indel frequency depending on the length are automatically selected, based on a naturally occuring SNP compoition and indel frequency.
# introducing viral sequences are from libraries obtained from NCBI (http://www.ncbi.nih.gov/genomes/viruses/).

my $script_path = $Bin;

my $reference = '';		# -r
my $output_prefix = 'out';	# -p
my $virus_lib = "$Bin/../Simulated_data/human_virus_ncbi.rn.fa";             # -l
my @virus_list = ();            # -vl
my $snp_rate = 0.001;		# -s
my $indel_rate = 0.0002;	# -i
my $virus_num = 100;            # -vn
my $min_len = 500;              # -ml
my $max_len = 10000;            # -xl
my $max_mut_rate = 0.05;        # -mr
my $hetero_rate = 2;            # -ht
my $help;			# -h

GetOptions(
    'ref|r=s' => \$reference,
    'pref|p=s' => \$output_prefix,
    'lib|l=s' => \$virus_lib,
    'virus_list|vl=s{,}' => \@virus_list,
    'snp=f' => \$snp_rate,
    'indel|i=f' => \$indel_rate,
    'virus_num|vn=i' => \$virus_num,
    'min_len|ml=i' => \$min_len,
    'max_len|xl' => \$max_len,
    'mut_rate|mr=f' => \$max_mut_rate,
    'het_rate|ht=f' => \$hetero_rate,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;
pod2usage(-verbose => 0) unless (-e $reference);


print STDERR "########## Virus insertion simulate options ##########\n";
print STDERR "Number-of-virus: $virus_num min-virus-len: $min_len max-virus-len: $max_len max-mut-rate: $max_mut_rate snp-rate=$snp_rate indel-rate=$indel_rate het-hm-rate: $hetero_rate pref=$output_prefix ref=$reference lib=$virus_lib virus-to-be-inserted=@virus_list\n";

=head1 SYNOPSIS

#### generate_sim_genome_VEI.pl [options] -r <reference file> -p <prefix of output files>
  Output:
   prefix.maternal.fa	      fasta file for maternal simulated reference
   prefix.paternal.fa	      fasta file for paternal simulated reference
   prefix.VE.vcf	          vcf file showing the reference positions of introduced viral elements
   prefix.snv.vcf	          vcf file showing the reference positions and bases of introduced SNPs and indels

  Options:
   --ref or -r <STR>          reference fasta file
   --pref or -p <STR>         prefix of outputfiles [default: out]
   --lib or -l <STR>          fasta file of virus sequence library [default: human_virus_ncbi.rn.fa]
   --virus_list or -vl <STR>  list of virus name to be introduced into reference (e.g., papilloma herpes); all virus sequences in the library are candidates of insertions if not specified
   --virus_num or -vn <INT>   number of virus segmanets to be introduced [default: 100]
   --min_len or -ml <INT>     minimum length of virus fragment [default: 500]
   --max_len or -xl <INT>     maximum length of virus fragment [default: 10000]
   --mut_rate or -mr <FLOAT>  maximum mutation rate of virus sequence [default: 0.05, which means that 1~5% mutations are introduced into a virus fragment]
   --snp or -s <FLOAT>        genomic mutation rate for SNPs [default: 0.01]
   --indel or -i <FLOAT>      genomic mutation rate for indels [default: 0.002]
   --het_rate or -ht <FLOAT>  hetero to homo ratio of introducing variants [default: 2]
   --help or -h               output help message

=cut

die "reference file is not specified: \n" if ($reference eq '');
die "virus library file is not specified: \n" if ($virus_lib eq '');


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


my %VE_list;
my %VE_annot;
my %VE_num;
$seq = '';
$header = '';

open (FILE, $virus_lib) or die "$virus_lib is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^>/){
	if ($seq ne ''){
	    if (@virus_list > 0){
		foreach my $id (@virus_list){
		    if ($header =~ /$id/i){
			my $id2 = $id;
			$id2 = 'HSV' if ($id =~ /herpes/i);
			$id2 = 'HSV' if ($id =~ /Epstein/i);
			$id2 = 'HPV' if ($id =~ /papillo/i);
			$id2 = 'HTLV' if ($id =~ /lymphotropic/i);
			$id2 = 'HBV' if ($id =~ /hepatitis/i);
			$id2 = 'HIV' if ($id =~ /immunodeficiency/i);
			push @{$VE_list{$id2}}, $seq;
			push @{$VE_annot{$id2}}, $header;
		    }
		}
	    }
	    else{
		$header = $1 if ($header =~ /(.+?)\|/);
		push @{$VE_list{$header}}, $seq;
		push @{$VE_annot{$header}}, $header;
	    }
	}
	$seq = '';
	$header = $1 if ($line =~ /^>(.+)/);
    }
    else{
	$seq .= $line;
    }
}
if ($seq ne ''){
    if (@virus_list > 0){
	foreach my $id (@virus_list){
	    if ($header =~ /$id/i){
		my $id2 = $id;
		$id2 = 'HSV' if ($id =~ /herpes/i);
		$id2 = 'HSV' if ($id =~ /Epstein/i);
		$id2 = 'HPV' if ($id =~ /papillo/i);
		$id2 = 'HTLV' if ($id =~ /lymphotropic/i);
		$id2 = 'HBV' if ($id =~ /hepatitis/i);
		$id2 = 'HIV' if ($id =~ /immunodeficiency/i);
		push @{$VE_list{$id2}}, $seq;
		push @{$VE_annot{$id2}}, $header;
	    }
	}
    }
    else{
	$header = $1 if ($header =~ /(.+?)\|/);
	push @{$VE_list{$header}}, $seq;
	push @{$VE_annot{$header}}, $header;
    }
}
close (FILE);

if (scalar keys %VE_list == 0){
    print STDERR "No sequence with specified virus list was found\n";
    die;
}
else{
    print STDERR "## Virus species to be introduced ##\n";
    foreach my $id (keys %VE_list){
	print STDERR "$id: ", scalar @{$VE_list{$id}}, "\n";
    }
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

my $ve_num = int ($virus_num / scalar keys %VE_list);

foreach my $id (keys %VE_list){
    my $ht_num = int ($ve_num * $hetero_rate2);
    my $hm_num = $ve_num - $ht_num * 2;
    ${$VE_num{$id}}{'ht'} = $ht_num;
    ${$VE_num{$id}}{'hm'} = $hm_num;
}

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

my $mei_seq = "$output_prefix.vei.fa";

&make_indel ($indel_1b_ht, $indel_1b_hm, $indel_2b_ht, $indel_2b_hm, $indel_3b_ht, $indel_3b_hm, $indel_4b_ht, $indel_4b_hm, $indel_5b_ht, $indel_5b_hm, $indel_6b_ht, $indel_6b_hm);

&make_snp ($snp_ht_no, $snp_hm_no);

foreach my $id (keys %VE_num){
    &make_ve ($id, ${$VE_num{$id}}{'ht'}, ${$VE_num{$id}}{'hm'});
}

&subst_genome ();


my $ve_info;
my $total_ve = 0;
foreach my $id (keys %VE_list){
    print STDERR "total $id: $count_var{$id}\n";
    $ve_info .= "$id:$count_var{$id},";
    $total_ve += $count_var{$id};
}
$ve_info =~ s/,*$//;
print STDERR "total SNPs: $count_var{'SNP'}\n";
print STDERR "total INDELs: $count_var{'INDEL'}\n";

print STDERR "total virus fragments: $total_ve\n";
print STDERR 'total susbstituted variants: ', $count_var, "\n";
print STDERR 'number of overlength: ', $count_overlength, "\n";

my $out_mat_fasta = $output_prefix . '.' . 'maternal.fa';
my $out_pat_fasta = $output_prefix . '.' . 'paternal.fa';
my $out_ve = $output_prefix . '.' . 'VE.vcf';
my $out_snv = $output_prefix . '.' . 'snv.vcf';

open (OUT1, "> $out_ve");
open (OUT2, "> $out_snv");
print OUT1 "#Total-VE:$total_ve($ve_info)SNP-rate:$snp_rate,Indel-rate:$indel_rate,Hetero-rate:$hetero_rate\n";
print OUT1 "#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
print OUT2 "#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
foreach my $chr (@chr_name){
    my $sum_shift_pos_pat = 0;
    my $sum_shift_pos_mat = 0;
    if (exists $indel_all{$chr}){
	foreach my $pos (sort {$a <=> $b} keys %{$indel_all{$chr}}){
	    my ($info, $parent, $annot) = split (/=/, ${$indel_all{$chr}}{$pos});
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
		$tag2 = 'VEI';
		$alt_len = $alt;
		$alt = '.';
		print OUT1 "$chr\t$pos\t$tag2\t$ref\t$alt\t.\tPASS\tSVTYPE=$tag;SVLEN=$alt_len;SVDIR=$strand;ALTPOS=$alt_pos;PARENT=$parent;SOURCE=$annot\n";
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

sub make_ve{
    my ($species, $ht_num, $hm_num) = @_;
    foreach my $chrno (sort keys %chr_rate){			# assign the positions of homozygous mobile elements
	for (my $i = 0; $i < int ($hm_num * $chr_rate{$chrno} + 0.5); $i++){
	    my $pos1 = int (rand ($chr_length{$chrno}));
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
	    ${$indel_all{$chrno}}{$pos1} = "$species=HM";
	    for (my $i = $pos1 - 10; $i <= $pos1 + 10; $i++){
		${$indel_allpos{$chrno}}{$i} = $species;
	    }
	}
    }
    foreach my $chrno (sort keys %chr_rate){			# assign the positions of heterozygous paternal mobile elements
	for (my $i = 0; $i < int ($ht_num * $chr_rate{$chrno} + 0.5); $i++){
	    my $pos1 = int (rand ($chr_length{$chrno}));
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
	    ${$indel_all{$chrno}}{$pos1} = "$species=Pat";
	    for (my $i = $pos1 - 10; $i <= $pos1 + 10; $i++){
		${$indel_allpos{$chrno}}{$i} = $species;
	    }
	}
    }
    foreach my $chrno (sort keys %chr_rate){			# assign the positions of heterozygous maternal mobile elements
	for (my $i = 0; $i < int ($ht_num * $chr_rate{$chrno} + 0.5); $i++){
	    my $pos1 = int (rand ($chr_length{$chrno}));
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
		my $ve_num = scalar @{$VE_list{$tag}};
		my $select_num = int (rand ($ve_num));
		my $ve_seq = ${$VE_list{$tag}}[$select_num];
		my $ve_annot = ${$VE_annot{$tag}}[$select_num];
		my $max_ve_len = $max_len;
		$max_ve_len = length $ve_seq if ($max_len > length $ve_seq);
		my $ve_len = int (rand ($max_ve_len - $min_len)) + $min_len;
		my $last_pos = length ($ve_seq) - $ve_len;
		my $start = int (rand ($last_pos));
		my $ve_frag = substr ($ve_seq, $start, $ve_len);
		
		my $mut_rate = int (rand ($max_mut_rate * 100) + 0.1) / 100;
		my $mut_num = int ($ve_len * $mut_rate);
		if ($mut_num > 0){
		    my %upos;
		    for (my $i = 1; $i <= $mut_num; $i++){
			my $mpos = int (rand ($ve_len));
			next if (exists $upos{$mpos});
			my $refbase = substr ($ve_frag, $mpos, 1);
			my $snp_base = $snp_nuc_A[rand (6)] if ($refbase eq 'A');
			$snp_base = $snp_nuc_G[rand (6)] if ($refbase eq 'G');
			$snp_base = $snp_nuc_C[rand (6)] if ($refbase eq 'C');
			$snp_base = $snp_nuc_T[rand (6)] if ($refbase eq 'T');
			substr ($ve_frag, $mpos, 1, $snp_base);
			$upos{$mpos} = 1;
		    }
		}
		my $strand = $strand[(rand (2))];
		if ($strand eq 'Minus'){
		    $ve_seq = reverse $ve_seq;
		    $ve_seq =~ tr/ACGT/TGCA/;
		}
		if ($parent eq 'HM'){
		    my $ref_base_2 = substr (${$chr_pat{$chr}}{$chr_divnum}, $pos_2, 0, $ve_frag);
		    $ref_base_2 = substr (${$chr_mat{$chr}}{$chr_divnum}, $pos_2, 0, $ve_frag);
		}
		elsif ($parent eq 'Pat'){
		    my $ref_base_2 = substr (${$chr_pat{$chr}}{$chr_divnum}, $pos_2, 0, $ve_frag);
		}
		elsif ($parent eq 'Mat'){
		    my $ref_base_2 = substr (${$chr_mat{$chr}}{$chr_divnum}, $pos_2, 0, $ve_frag);
		}
		${$indel_all{$chr}}{$pos} = "$tag-$ref_base-$ve_len-$strand=$parent=$ve_annot";
		print OUT1 ">$tag-$chr-$pos $strand $parent $ve_annot\n";
		print OUT1 $ve_frag, "\n";
		$count_var{$tag} ++;
	    }
	}
	$pre_count = $count_var;
    }
    close (OUT1);
}

