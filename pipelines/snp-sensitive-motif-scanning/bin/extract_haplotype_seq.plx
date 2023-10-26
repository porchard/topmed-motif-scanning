#!/usr/bin/perl -w
# $Id: extract_haplotype_seq.plx 4582 2011-12-20 18:14:59Z pchines $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);
use GTB::File::VCF;
use GTB::FASTA;
our %Opt;
our $Ofh;   # output filehandle
our $Lfh;   # leftover indels filehandle
our @Bfh;   # alt vars per hap BED files
our $Fasta;
our $EMPTY = q{};
our @Allele = qw(REF ALT);

=head1 NAME

extract_haplotype_seq.plx - extract haplotypes around variants

=head1 SYNOPSIS

Extract all possible haplotypes involing variants from a VCF file, with 29bp of
flanking sequence.  All variants within 29bp of each other will be combined.

  extract_haplotype_seq.plx -all -o out.fa in.vcf

Extract non-reference haplotypes for specified sample (SAMP1), extending the
haplotype to include all variants within 29bp of one another (assumes phased
VCF):

  extract_haplotype_seq.plx -s SAMP1 -o out.fa in.vcf

For complete documentation, run C<extract_haplotype_seq.plx -man>

=head1 DESCRIPTION

Extract haplotypes.

=cut

#------------
# Begin MAIN 
#------------

process_commandline();
$Fasta = GTB::FASTA->new($Opt{genome});
my $vcf = GTB::File::VCF->new(file => $ARGV[0]);
$Ofh = Open($Opt{output}, 'w');
if ($Opt{leftovers}) {
    $Lfh = Open($Opt{leftovers}, 'w');
    print $Lfh $vcf->skipped_lines();
}
if ($Opt{bed}) {
    for my $i (1,2) {
        push @Bfh, Open("$Opt{bed}.$i.bed.gz", 'w');
    }
}
my @block;
while (my $rh = $vcf->next_hash) {
    my $size = length($rh->{ALT}) - length($rh->{REF});
    if ($rh->{REF} !~ /^[ACGT]+$/ || $rh->{ALT} !~ /^[ACGT]+$/
            || $size > $Opt{insertion} || -$size > $Opt{deletion}) {
        if ($Lfh) {
            print $Lfh $vcf->current_line();
        }
        next;
    }
    $rh->{POSEND} = $rh->{POS} - 1 + length($rh->{REF});
    while (@block && ($Opt{single}
                || $block[0]{CHROM} ne $rh->{CHROM}
                || $block[0]{POSEND} + $Opt{flanking} < $rh->{POS}
            )) {
        process_block(\@block);
    }
    push @block, $rh;
}
while (@block) {
    process_block(\@block);
}
close $Ofh or die "Error closing $Opt{output}, $!\n";

#------------
# End MAIN
#------------

sub process_block {
    my ($ra_block) = @_;
    my %written;
    my $ra_haps = $Opt{sample} ? make_haplotypes($ra_block, $Opt{sample})
                               : make_all_haplotypes($ra_block);
    for my $ra_hap (@$ra_haps) {
        my $chrom  = $ra_hap->[0]{var}{CHROM};
        my $rstart = $ra_hap->[0]{var}{POS} - $Opt{flanking};
        my $rcurr = $rstart;
        my $radj  = $rstart;
        my $hname = $chrom;
        my $seq = $EMPTY;
        for my $v (@$ra_hap) {
            my $pos = $v->{var}{POS};
            my $offset = $pos - $radj;
            my $ar = substr($Allele[$v->{allele}],0,1);
            if ($pos > $rcurr) {
                $seq .= lc $Fasta->seq($chrom, $rcurr, $pos - 1);
            }
            elsif ($pos < $rcurr) {
                warn "Variants overlap; please rewrite VCF file to combine"
                    . " $chrom:$pos $v->{var}{REF}/$v->{var}{ALT}"
                    . " with previous variant\nSkipping this variant for now.\n";
                next;
            }
            # OPTION: test ref allele vs ref seq?
            $hname .= "_$offset$ar\@$pos$v->{var}{REF}/$v->{var}{ALT}";
            $seq   .= $v->{var}{$Allele[$v->{allele}]};
            $rcurr  = $v->{var}{POS} + length($v->{var}{REF});
            $radj  += length($v->{var}{REF})
                - length($v->{var}{$Allele[$v->{allele}]});
        }
        $hname .= ":$rstart-" . ($rcurr+$Opt{flanking}-1);
        $seq .= lc $Fasta->seq($chrom, $rcurr, $rcurr+$Opt{flanking}-1);
        if (!$written{ $hname }) {
            write_seq($hname, $seq);
            $written{ $hname } = 1;
        }
    }
    @$ra_block = ();
}

sub make_haplotypes {
    my ($ra_block, $sample) = @_;
    my @haps = ([],[]);
    for my $var (@$ra_block) {
        my $gt = $var->{$sample} || die "No genotypes for sample $sample\n";
        if ($gt =~ /^0(?:\D0)?(:|$)/ && !$Opt{reference}) {
            next;   # only ref alleles, so don't need to handle specially
        }
        if ($gt =~ /^(\d)\|(\d)\b/) {
            for (my $i = 0; $i < @haps; ++$i) {
                my $a = $i % 2 ? $2 : $1;
                if ($a || $Opt{reference}) {
                    push @{ $haps[$i] }, { allele => $a, var => $var };
                }
            }
        }
        elsif ($gt =~ /^(\d)\/(\d)\b/) {
            if ($1 eq $2) {
                for my $h (@haps) {
                    if ($1 || $Opt{reference}) {
                        push @$h, { allele => $1, var => $var };
                    }
                }
            }
            elsif ($Opt{bed}) {
                warn "Can't create more than two haps when using --bed"
                    . " option, so ignoring\nunphased het at"
                    . " $var->{CHROM}:$var->{POS}:$var->{REF}/$var->{ALT}\n";
            }
            else {
                warn "Phase unknown for $sample het at"
                    . " $var->{CHROM}:$var->{POS}:$var->{REF}/$var->{ALT}\n";
                warn "Applying both alleles to make all possible haps\n";
                # copy haps to double, add each allele to one copy
                my @r = $Opt{reference} ? ({allele => 0, var => $var}) : ();
                my @a = ({allele => 1, var => $var});
                my @cp1 = map { [@$_, @r] } @haps;
                my @cp2 = map { [@$_, @a] } @haps;
                @haps = (@cp1, @cp2);
            }
        }
        elsif ($gt !~ /^\./) {
            die "Error parsing genotype '$gt' for $sample at"
                . " $var->{CHROM}:$var->{POS}:$var->{REF}/$var->{ALT}\n";
        }
    }
    if ($Opt{bed}) {
        for my $i (0,1) {
            my $bfh = $Bfh[$i];
            for my $rh (@{ $haps[$i] }) {
                my $var = $rh->{var};
                my $ps = $var->{POS} - 1;
                if (length($var->{REF}) != length($var->{ALT})
                        && substr($var->{REF},0,1) eq substr($var->{ALT},0,1)) {
                    ++$ps;
                }
                my $chr = $var->{CHROM};
                $chr =~ s/^(?:chr)?/chr/;
                print $bfh join("\t", $chr, $ps,
                        $var->{POS} + length($var->{REF}) - 1,
                        "$var->{CHROM}:$var->{POS}_$var->{REF}/$var->{ALT}",
                        ), "\n";
            }
        }
    }
    @haps = grep { scalar(@$_) } @haps;
    return \@haps;
}

sub make_all_haplotypes {
    my ($ra_block) = @_;
    my @haps;
    for (my $i = $Opt{reference} ? 0 : 1; $i < 2 ** scalar(@$ra_block); ++$i) {
        for (my $j = 0; $j < @$ra_block; ++$j) {
            # could also be accomplished by doubling, as with unphased
            my $a = $i & (2 ** $j) ? 1 : 0;
            if ($a || $Opt{reference}) {
                push @{ $haps[$i] }, { allele => $a, var => $ra_block->[$j] };
            }
        }
    }
    @haps = grep { $_ } @haps;
    return \@haps;
}

sub write_seq {
    my ($name, $seq) = @_;
    print $Ofh ">$name\n$seq\n";
}

sub process_commandline {
    %Opt = (genome  => 'hg19',
            flanking=> 29,
            output  => '-',
            );
    GetOptions(\%Opt, qw(all bed=s deletion=i flanking=i genome=s insertion=i
        leftovers=s output=s reference sample|s=s single
        manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) {
        die "extract_haplotype_seq.plx, ", q$Revision: 4582 $, "\n";
    }
    if (!@ARGV) {
        pod2usage("Input VCF file is required");
    }
    if (!$Opt{all} && !$Opt{sample}) {
        pod2usage("Either --sample or --all haplotypes must be specified");
    }
    if ($Opt{all} && $Opt{sample}) {
        pod2usage("--all cannot be combined with --sample");
    }
    if ($Opt{bed} && !$Opt{sample}) {
        pod2usage("--bed can only be used with --sample");
    }
    if ($Opt{flanking} < 1) {
        pod2usage("--flanking must be a positive integer");
    }
    for my $o (qw(deletion insertion)) {
        if (!defined $Opt{$o}) {
            $Opt{$o} = $Opt{flanking} + 1;
        }
        elsif ($Opt{$o} < 0) {
            pod2usage("--$o must be a non-negative value");
        }
    }
}

=head1 OPTIONS

=over 4

=item B<--all>

Write all possible haplotypes.

=item B<--bed> PREFIX

Will write PREFIX.1.bed.gz and PREFIX.2.bed.gz with locations of alternate
alleles in haplotype 1 and 2 from specified sample.  Will not work unless
C<--sample> is provided.

=item B<--deletion> N

Maximum size of deletion that will be included in haplotypes.  Default, one
more than size of flanking region.

=item B<--flanking> 29

Bases of flanking sequence to extract on each side of variant.

=item B<--genome> hg19

=item B<--genome> /data/projects/solexa/reference/fasta/hg19/hg19.mfa.

Genome build or multi-fasta file.

=item B<--insertion> N

Maximum size of insertion that will be included in haplotypes.  Default, one
more than size of flanking region.

=item B<--leftovers> FILE

If supplied, will write any variants that do not meet criteria
(multi-allelic, indels that exceed size constraint, etc.) to the specified
file.

=item B<--output> FILE

Destination file, otherwise written to stdout.

=item B<--reference>

Include reference haplotypes.

=item B<--sample> NAME

Evaluate haplotypes for this sample.  This is required, unless C<--all> is specified.

=item B<--single>

Don't write haplotypes; write single variants only.

=item B<--help|--manual>

Display documentation.  C<--help> gives a brief synopsis, C<--manual> provides
complete documentation.

=back

=head1 AUTHOR

 Peter Chines - pchines@mail.nih.gov
 Narisu Narisu - narisu@mail.nih.gov

=head1 LEGAL

This software/database is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software/database is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use. 

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S.  Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose. 

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
FUSION Research Group" as the citation. 

=cut
