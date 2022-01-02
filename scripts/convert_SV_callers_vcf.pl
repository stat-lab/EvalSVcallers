#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib $Bin;

# covert an output file from the specified SV callers to the vcf format compatible to our evaluation study

my $tool = '';
my $class = '';
my $min_sv_size = 30;
my $max_sv_size = 20000000;
my $non_human = 0;
my $help;

GetOptions(
    'tool|t=s' => \$tool,
    'class|c=s' => \$class,
    'len|l=i' => \$min_sv_size,
    'xlen|xl=i' => \$max_sv_size,
    'non_human|nh' => \$non_human,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  convert_SV_callers_vcf.pl <option> [output file(s) from a SV caller] > [output vcf file]
  [Example] convert_SV_callers_vcf.pl -t Pindel pindel.DEL.out pindel.DUP.out > Pindel.vcf

  Options:
   --tool or -t <STR>       a tool name (69 tools used in our study) [mandatory]
   --class or -c <STR>      specify when tool is Mobster or PBHoney (MEI|NUMT|VEI for Mobster, NGM for PBHoney-NGM)
   --len or -l <INT>        minimum size (bp) of DEL/DUP/INV to be kept [default: 30]
   --xlen or -xl <INT>      maximum size (bp) of SV to be kept [default: 20000000]
   --non_human or -nh       sample is non-human species [default: false]
   --help or -h             output help message
   
=cut

if (($tool eq 'PBHoney') and ($class eq 'NGM')){
    $tool = 'PBHoney-NGM';
}
elsif ($tool eq 'Mobster'){
    die "--class option not specified or inproperly specified:\n" if ($class !~ /^MEI$|^NUMT$|^VEI$/);
}

my $script_path = $Bin;

my @convert_script = <$script_path/convert_$tool*.pl>;
@convert_script = <$script_path/vcf_convert/convert_$tool*.pl> if (@convert_script == 0);

if (@convert_script == 0){
    die "TOOL name specified is absent in 69 tools we analyzed or does not match the names:\n";
}
elsif (@convert_script > 1){
    die "TOOL name specified is redundant. Specify more strict name:\n";
}

my %vcf;
foreach my $var_file (@ARGV){
    my @output = `$convert_script[0] $var_file $class` if ($tool eq 'Mobster');
    @output = `$convert_script[0] $var_file` if ($tool ne 'Mobster');
    foreach (@output){
        chomp $_;
        next if ($_ =~ /^#|^$/);
        my ($chr, $pos) = split (/\t/, $_);
        next if ($chr !~ /^c*h*r*[\dXY]+$/) and ($non_human == 0);
        my $type = $1 if ($_ =~ /SVTYPE=(.+?);/);
        my $chr02d = $chr;
        $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
        my $len = 0;
        $len = $1 if ($_ =~ /SVLEN=-*(\d+)/);
        next if ($len < $min_sv_size) and ($type ne 'INS');
        next if ($len > $max_sv_size);
        ${${$vcf{$chr02d}}{$pos}}{$type} = $_;
    }
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
        foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
            print "${${$vcf{$chr}}{$pos}}{$type}\n";
        }
    }
}

