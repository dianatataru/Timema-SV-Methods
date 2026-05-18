#!/usr/bin/perl

# filter vcf files based on coverage (SNPs and individuals)
# Usage: perl script.pl --keep KeepSNPs.txt --keepinds KeepInds.txt --out output.vcf input.vcf [input2.vcf ...]

use strict;
use warnings;
use Getopt::Long;

my $keep_file      = '';
my $keep_inds_file = '';
my $out_file       = '';

GetOptions(
    'keep=s'     => \$keep_file,
    'keepinds=s' => \$keep_inds_file,
    'out=s'      => \$out_file,
) or die "Usage: perl script.pl --keep <KeepSNPs.txt> --keepinds <KeepInds.txt> --out <output.vcf> <input.vcf> [input2.vcf ...]\n";

die "Error: --keep file required\n"               unless $keep_file;
die "Error: --keepinds file required\n"           unless $keep_inds_file;
die "Error: --out file required\n"                unless $out_file;
die "Error: at least one input VCF required\n"    unless @ARGV;

## read SNP keep flags
my @keep_snps;
open(my $snp_fh, '<', $keep_file) or die "Failed to read keep file: $keep_file\n";
while (<$snp_fh>) {
    chomp;
    push(@keep_snps, $_);
}
close($snp_fh);

## read individual keep flags
my @keep_inds;
open(my $inds_fh, '<', $keep_inds_file) or die "Failed to read keep inds file: $keep_inds_file\n";
while (<$inds_fh>) {
    chomp;
    push(@keep_inds, $_);
}
close($inds_fh);

## get indices of individuals to keep (0-based, offset by 9 for VCF fixed columns)
## VCF columns 1-9 are: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT -- always kept
my @ind_cols;
for my $i (0..$#keep_inds) {
    if ($keep_inds[$i] == 1) {
        push(@ind_cols, $i + 9);
    }
}
my $kept_inds = scalar @ind_cols;

foreach my $in (@ARGV) {
    my @keep_snps_copy = @keep_snps;
    my $cnt_snps = 0;

    open(my $in_fh,  '<', $in)       or die "Could not read the infile = $in\n";
    open(my $out_fh, '>', $out_file) or die "Could not write the outfile = $out_file\n";

    while (<$in_fh>) {
        chomp;
        my $flag = 0;

        if (m/^\#\#/) {
            ## meta-header lines (##INFO, ##FORMAT etc.), always print as-is
            print $out_fh "$_\n";
            next;
        }
        elsif (m/^\#/) {
            ## column header line (#CHROM POS ID ...), subset to kept individuals
            my @fields  = split(/\t/, $_);
            my @fixed   = @fields[0..8];
            my @samples = @fields[@ind_cols];
            print $out_fh join("\t", @fixed, @samples) . "\n";
            next;
        }
        elsif (m/^\S+/) {
            ## sequence/variant line -- check SNP filter
            $flag = shift(@keep_snps_copy);
            if ($flag == 1) {
                $cnt_snps++;
            }
        }
        else {
            print "Warning, failed to match line:\n$_\n";
            $flag = 0;
        }

        if ($flag == 1) {
            ## subset to kept individual columns
            my @fields  = split(/\t/, $_);
            my @fixed   = @fields[0..8];
            my @samples = @fields[@ind_cols];
            print $out_fh join("\t", @fixed, @samples) . "\n";
        }
    }

    close($in_fh);
    close($out_fh);

    print "Finished filtering $in\nRetained $cnt_snps variable loci and $kept_inds individuals\n";
}
