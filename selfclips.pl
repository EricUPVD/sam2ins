#!/usr/bin/perl
use strict;
use warnings;

#getting the infile name from arguments
die "Please provide a clip file (sam2clip output). Command line is : perl selfclips.pl <clip_data file>" if @ARGV < 1;
die "Please provide the selfclips.conf file in the working directory" unless -e "selfclips.conf";
my $clipfile = shift;        #clip_data file
$clipfile =~ /clip_(.+)\.tsv/;
my $outname = $1;

#getting the data from the configuration file
my $tsd_min;            # minimum size of the tsd
my $tsd_max;            # maximum size of the tsd
my $min_clipped_size;   # minimum number of clipped bases for a "robust" clipped read
my $min_robust_clip;    # min nb clipped reads at left and at right of the breakpoint that satisfy to the MIN_CLIPPED_SIZE condition
my $min_for_bp;         # minimun number of clipped reads for a valid bp
my $flank;              # nb of nucleotides added at each part of a false positive breakpoint to create a "no landing site"

open CONF, "selfclips.conf";
while (<CONF>) {
    chomp;
    $tsd_min = (split "=")[1] if $_ =~ /^TSD_MIN/;
    $tsd_max = (split "=")[1] if $_ =~ /^TSD_MAX/;
    $min_clipped_size = (split "=")[1] if $_ =~ /^MIN_CLIPPED_SIZE/;
    $min_robust_clip = (split "=")[1] if $_ =~ /^MIN_ROBUST_CLIP/;
    $min_for_bp = (split "=")[1] if $_ =~ /^MIN_FOR_BP/;
    $flank = (split "=")[1] if $_ =~ /^FLANK/;
    }
close CONF;

die "Please provide a value for TSD_MIN" if $tsd_min !~ /\d+/;
die "Please provide a value for TSD_MAX" if $tsd_max !~ /\d+/;
die "Please provide a value for MIN_CLIPPED_SIZE" if $min_clipped_size !~ /\d+/;
die "Please provide a value for MIN_ROBUST_CLIP" if $min_robust_clip !~ /\d+/;
die "Please provide a value for MIN_FOR_BP" if $min_for_bp !~ /\d+/;
die "Please provide a value for FLANK" if $flank !~ /\d+/;

print "\n_______Analysis of ".$clipfile."__\n";

#storing bedfile data
open IN, $clipfile;
my %MS_clip;
my %SM_clip;

while (<IN>) {
    chomp;
    my ($fmate_chr, $fmate_st, $fmate_end, $fmate_type, $fmate_clip_nb, $pair_type) = (split/\t/)[1,2,3,4,5,11];
    if ($pair_type =~ /^clip/) {
        #fmate is MS
        if ($fmate_type =~ /^MS/) {
            #fmate is clipped at the right end. The "ins" position is therefore the end of the alignment.
            $MS_clip{$fmate_chr}{$fmate_end} .= "$fmate_clip_nb,";
        #fmate is SM
        } elsif ($fmate_type =~ /^SM/) {
            #fmate is clipped at the left end. The "ins" position is therefore the start of the alignment.
            $SM_clip{$fmate_chr}{$fmate_st} .= "$fmate_clip_nb,";
        }
    }
}
close IN;
print "reads data stored, analysis begins\n";


open OUT, ">selfclip_$outname.bed";
foreach my $chr (sort keys %MS_clip) {
    foreach my $ms_pos (sort{$a<=>$b} keys %{$MS_clip{$chr}}) {
        my $nb_ms_clip = 0;
        my $nb_robust_ms_clip = 0;
        my $nb_sm_clip = 0;
        my $nb_robust_sm_clip = 0;
        my @ms_clips = split/,/,$MS_clip{$chr}{$ms_pos};
        foreach my $clipped_nb (@ms_clips) {
            $nb_ms_clip ++;
            $nb_robust_ms_clip ++ if $clipped_nb >= $min_clipped_size;
        }
        if ($nb_robust_ms_clip >= $min_robust_clip and $nb_ms_clip >= $min_for_bp) {
            my $ms_st = $ms_pos - $flank;
            $ms_st = $ms_pos if $ms_pos < $flank;
            my $ms_end = $ms_pos + $flank;
            print OUT "$chr\t$ms_st\t$ms_end\n" ;
        } else {
            my $ref_sm_pos = 0;
            foreach my $sm_pos ($ms_pos - $tsd_max..$ms_pos - $tsd_min) {
                if (exists $SM_clip{$chr}{$sm_pos}) {
                    $ref_sm_pos = $sm_pos;
                    my @sm_clips = split/,/,$SM_clip{$chr}{$sm_pos};
                    foreach my $clipped_nb (@sm_clips) {
                        $nb_sm_clip ++;
                        $nb_robust_sm_clip ++ if $clipped_nb >= $min_clipped_size;
                    }
                }
            }
            if ($nb_ms_clip + $nb_sm_clip >= $min_for_bp and $nb_robust_ms_clip + $nb_robust_sm_clip >= $min_robust_clip) {
                if (($ms_pos - $ref_sm_pos) + 1 >= $flank * 2) {
                    print OUT "$chr\t$ref_sm_pos\t$ms_pos\n";
                } else {
                    my $missing = $flank * 2 - (($ms_pos - $ref_sm_pos)+1);
                    my $flank2 = $missing/2;
                    my $sm_st = int($ref_sm_pos - $flank2);
                    $sm_st = $ref_sm_pos if $ref_sm_pos < $flank2;
                    my $sm_end = int($ms_pos + $flank2);
                    print OUT "$chr\t$sm_st\t$sm_end\n"
                }
            }
        }
    }
}
foreach my $chr (sort keys %SM_clip) {
    foreach my $sm_pos (sort{$a<=>$b} keys %{$SM_clip{$chr}}) {
        my $nb_sm_clip = 0;
        my $nb_robust_sm_clip = 0;
        my @sm_clips = split/,/,$SM_clip{$chr}{$sm_pos};
        foreach my $clipped_nb  (@sm_clips) {
            $nb_sm_clip ++;
            $nb_robust_sm_clip ++ if $clipped_nb >= $min_clipped_size;
        }
        if ($nb_robust_sm_clip >= $min_robust_clip and $nb_sm_clip >= $min_for_bp) {
            my $sm_st = $sm_pos - $flank;
            $sm_st = $sm_pos if $sm_pos < $flank;
            my $sm_end = $sm_pos + $flank;
            print OUT "$chr\t$sm_st\t$sm_end\n"  ;
        }
    }
}
close OUT;
