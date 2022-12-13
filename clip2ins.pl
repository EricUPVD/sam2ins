#!/usr/bin/perl
use strict;
use warnings;

#getting the infile name from arguments
die "Please provide a clip file (sam2clip output). Command line is : perl clip2bp.pl <clip_data file>" if @ARGV < 1;
die "Please provide the clip2ins.conf file in the working directory" unless -e "clip2ins.conf";
my $clipfile = shift;        #clip_data file
$clipfile =~ /clip_(.+)\.tsv/;
my $outname = $1;

#getting the data from the configuration file
my $read_length;
my $max_for_conc;           # max size for an alignment to be considered as concordant (inherited from sam2clip.pl)
my $tsd_min;                # minimum size of the tsd
my $tsd_max;                # maximum size of the tsd
my $dist_from_bp;           # maximum distance from the bp for unclip_disc mates
my $min_clipped_size;       # minimum number of clipped bases for a "robust" clipped read
my $min_robust_clip;        # min nb of ms or sm clipped reads that satisfy to the $min_clipped_size condition
my $min_sum_robust_clip;    # min nb of total clipped reads that satisfy to the $min_clipped_size condition
my $min_left;               # min nb of ms clipped or unclip left reads (no effect if <= $min_robust_clip)
my $min_for_bp;             # min nb of reads to validate a bp
my $min_size_for_blast;     # minimum size of fragment to be blasted against TE database
my $min_evalue;             # minimum evalue for TE identification by blast
my $min_te_pident;          # minimum percent identity for TE identification by blast
my $min_TE_pcent;           # minimun percent of clipped or disc reads that matched to the same TE or TE family
my $num_threads;            # number of threads used by blast for TE identification
my $nls;                    # no landing sites data
my $TEdb;                   # TE database
my $clipfasta;              # fasta sequences extracted from $clipfile
my $path_to_bedtools;       # path to bedtools command
my $path_to_blast;          # path to blast command

open CONF, "clip2ins.conf";
while (<CONF>) {
    chomp;
    $read_length = (split "=")[1] if $_ =~ /^READ_LENGTH/;
    $max_for_conc = (split "=")[1] if $_ =~ /^MAX_FOR_CONC/;
    $tsd_min = (split "=")[1] if $_ =~ /^TSD_MIN/;
    $tsd_max = (split "=")[1] if $_ =~ /^TSD_MAX/;
    $dist_from_bp = (split "=")[1] if $_ =~ /^DIST_FROM_BP/;
    $min_clipped_size = (split "=")[1] if $_ =~ /^MIN_CLIPPED_SIZE/;
    $min_robust_clip = (split "=")[1] if $_ =~ /^MIN_ROBUST_CLIP/;
    $min_sum_robust_clip = (split "=")[1] if $_ =~ /^MIN_SUM_ROBUST_CLIP/;
    $min_left = (split "=")[1] if $_ =~ /^MIN_LEFT_AND_RIGTH/;
    $min_for_bp = (split "=")[1] if $_ =~ /^MIN_FOR_BP/;
    $min_size_for_blast = (split "=")[1] if $_ =~ /^MIN_SIZE_FOR_BLAST/;
    $min_evalue = (split "=")[1] if $_ =~ /^MIN_EVALUE/;
    $min_te_pident = (split "=")[1] if $_ =~ /^MIN_TE_PIDENT/;
    $num_threads = (split "=")[1] if $_ =~ /^NUM_THREADS/;
    $min_TE_pcent = (split "=")[1] if $_ =~ /^MIN_TE_PCENT/;
    $nls = (split "=")[1] if $_ =~ /^NLS/;
    $TEdb = (split "=")[1] if $_ =~ /^TEDB/;
    $clipfasta = (split "=")[1] if $_ =~ /^CLIPFASTA/;
    $path_to_bedtools = (split "=")[1] if $_ =~ /^PATH_TO_BEDTOOLS/;
    $path_to_blast = (split "=")[1] if $_ =~ /^PATH_TO_BLAST/;
    }
close CONF;

die "Please provide a value for READ_LENGTH" if $read_length !~ /\d+/;
die "Please provide a value for MAX_FOR_CONC" if $max_for_conc !~ /\d+/;
die "Please provide a value for TSD_MIN" if $tsd_min !~ /\d+/;
die "Please provide a value for TSD_MAX" if $tsd_max !~ /\d+/;
die "Please provide a value for DIST_FROM_BP" if $dist_from_bp !~ /\d+/;
die "Please provide a value for MIN_CLIPPED_SIZE" if $min_clipped_size !~ /\d+/;
die "Please provide a value for MIN_ROBUST_CLIP" if $min_robust_clip !~ /\d+/;
die "Please provide a value for MIN_SUM_ROBUST_CLIP" if $min_sum_robust_clip !~ /\d+/;
die "Please provide a value for MIN_LEFT_AND_RIGTH" if $min_left !~ /\d+/;
die "Please provide a value for MIN_FOR_BP" if $min_for_bp !~ /\d+/;
die "Please provide a value for MIN_SIZE_FOR_BLAST" if $min_size_for_blast !~ /\d+/;
die "Please provide a value for MIN_EVALUE" if $min_evalue !~ /\d+/;
die "Please provide a value for MIN_TE_PIDENT" if $min_te_pident !~ /\d+/;
die "Please provide a value for NUM_THREADS" if $num_threads !~ /\d+/;
die "Please provide a value for MIN_TE_PCENT" if $min_TE_pcent !~ /\d+/;
die "Please provide a value for NLS" if $nls !~ /\H+/;
die "Please provide a value for TEDB" if $TEdb !~ /\H+/;
die "Please provide a value for PATH_TO_BEDTOOLS" if $path_to_bedtools !~ /\H+/;
die "Please provide a value for PATH_TO_BLAST" if $path_to_blast !~ /\H+/;

my $dist_from_bord = $max_for_conc - $read_length * 2;  # maximum distance form the borders of the TE for disc mates or clipped parts
my $clipfastaidx = "$clipfasta.fai";
my $min_right = $min_left;  # min nb of sm clipped or unclip right reads (no effect if <= $min_robust_clip)

print "\n_______Analysis of ".$clipfile."__\n";

#storing clipfile data
open IN, $clipfile;
my %MS_clip;
my %SM_clip;
my %unclip_disc;
while (<IN>) {
    chomp;
    my ($read_name, $fmate_chr, $fmate_st, $fmate_end, $fmate_type, $fmate_clip_nb, $smate_chr, $smate_start, $smate_end, $smate_type, $smate_clip_nb, $pair_type) = split/\t/;
    if ($pair_type =~ /^clip/) {
        #fmate is MS
        if ($fmate_type =~ /^MS/) {
            #fmate is clipped at the right end. The "ins" position is therefore the end of the alignment.
            $MS_clip{$fmate_chr}{$fmate_end} .= "$_,";
        #fmate is SM
        } elsif ($fmate_type =~ /^SM/) {
            #fmate is clipped at the left end. The "ins" position is therefore the start of the alignment.
            $SM_clip{$fmate_chr}{$fmate_st} .= "$_,";
        }
    } elsif ($pair_type =~ /^unclip/) {
        $unclip_disc{$fmate_chr}{$fmate_st} .= "$_,";
    }
}
close IN;
print "reads data stored\n";

#storing nls data
open IN, $nls;
my %no_landing_sites;
while (<IN>) {
    chomp;
    my ($nls_chr, $nls_st, $nls_end) = split/\t/;
    $no_landing_sites{$nls_chr}{$nls_st} = $_;
}
close IN;
print "nls data stored, analysis begins... ";


open OUT1, ">", "bp_".$outname.".bed";
my $i = 0;
foreach my $chr (sort keys %MS_clip) {

    # looking for ms reads that overlap the size of a tsd with sm reads
    foreach my $ms_pos (sort{$a<=>$b} keys %{$MS_clip{$chr}}) {
        my @ms_clips = split/,/,$MS_clip{$chr}{$ms_pos};
        my $nb_robust_ms_clip = 0;
        my $nb_left = 0;
        foreach my $mc (@ms_clips) {
            $nb_left ++;
            my $clipped_nb = (split/\t/, $mc)[5];
            $nb_robust_ms_clip ++ if $clipped_nb >= $min_clipped_size;
        }

        # a minimum nb of robust clippped reads is necessary
        if ($nb_robust_ms_clip >= $min_robust_clip) {
            my $sm_loop = 0;
            foreach my $sm_pos (($ms_pos - $tsd_max) + 1..($ms_pos - $tsd_min) + 1) {
                my @selected_sm;
                my @selected_ud;
                if (exists $SM_clip{$chr}{$sm_pos}) {
                    $sm_loop ++;
                    my @sm_clips = split/,/,$SM_clip{$chr}{$sm_pos};
                    my $nb_robust_sm_clip = 0;
                    my $nb_right = 0;
                    foreach my $sc (@sm_clips) {
                        $nb_right ++;
                        my $clipped_nb = (split/\t/, $sc)[5];
                        $nb_robust_sm_clip ++ if $clipped_nb >= $min_clipped_size;
                    }

                    # a minimum nb of clippped reads is necessary, particularly of the robust type
                    if ($nb_robust_sm_clip >= $min_robust_clip and $nb_robust_ms_clip + $nb_robust_sm_clip >= $min_sum_robust_clip) {
                        push @selected_sm, $sm_pos;
                        push @selected_sm, $SM_clip{$chr}{$sm_pos};

                        # get unclip_disc mates that exist on both sides of the bp
                        # at the left side
                        foreach my $left_pos ($ms_pos - $dist_from_bp..$ms_pos) {
                            if (exists $unclip_disc{$chr}{$left_pos}) {
                                $nb_left ++;
                                push @selected_ud, "$unclip_disc{$chr}{$left_pos};L";
                            }
                        }
                        # at the right side
                        foreach my $right_pos ($sm_pos..$sm_pos + $dist_from_bp - $read_length) {
                            if (exists $unclip_disc{$chr}{$right_pos}) {
                                $nb_right++;
                                push @selected_ud, "$unclip_disc{$chr}{$right_pos};R";
                            }
                        }

                        # a minimum nb of discordant or clipped reads is necessary on both sides of the bp position
                        if ($nb_left >= $min_left and $nb_right >= $min_right and $nb_left + $nb_right >= $min_for_bp) {

                            ####################################################
                            #   At this stage, for each validated ms_pos       #
                            #   @selected_sm contains :                        #
                            #   - sm_pos                                       #
                            #   - SM_clip data (values from %SM_clip)          #
                            #   @selected_ud contains :                        #
                            #   - ud_data (values from %unclip_disc) + R or L  #
                            ####################################################

                            my $sm_pos = shift @selected_sm;
                            my $tsd_length = ($ms_pos-$sm_pos) +1;

                            # the sm_pos should not lie in NLSs
                            # get nls start positions located before bp_pos
                            my @nls_starts = grep {$_ <= $sm_pos} (sort{$a<=>$b} keys %{$no_landing_sites{$chr}});
                            # get nls start position located immediately before sm_pos
                            my $nls_st = pop @nls_starts;
                            # get nls end position corresponding to nls_st
                            if (defined $nls_st) {  # $sm_pos may be located before the first nls of the chr
                                my $nls_end = (split/\t/,$no_landing_sites{$chr}{$nls_st})[2];
                                next if $nls_end >= $sm_pos;
                            }

                            my $blast_query_count = 0;      # nb of blast queries indicating the nb of candidate mates for matching a TE
                            open BED, ">", "$chr:$sm_pos-$ms_pos.bed";
                            #get ms_clip mates or parts that lie within annotated TEs to identify the transposed one
                            #IMPORTANT : BED starts are zero-based and BED ends are one-based
                            foreach my $msc (@ms_clips) {
                                my $name = (split/\t/,$msc) [0];
                                my ($fmate_chr, $fmate_st, $fmate_end, $fmate_type, $fmate_clipped) = (split/\t/,$msc) [1,2,3,4,5];
                                my ($smate_chr, $smate_st, $smate_end, $smate_type, $smate_clipped) = (split/\t/,$msc) [6,7,8,9,10];# smate is distant from bp if disc, at bp if conc
                                my $pair_type = (split/\t/,$msc) [11];
                                #$ms_clip_count ++;
                                if ($fmate_clipped >= $min_size_for_blast) {
                                    #clipped part of fmate wanted
                                    my $fmate_clip_st;
                                    my $fmate_clip_end;
                                    if ($fmate_type =~ /\+/) {
                                        $fmate_clip_st = ($read_length - $fmate_clipped);
                                        $fmate_clip_end = $read_length;
                                    } else {
                                        $fmate_clip_st = 0;
                                        $fmate_clip_end = $fmate_clipped + 1;
                                    }
                                    print BED "$name/1\t$fmate_clip_st\t$fmate_clip_end\t$name"."_$fmate_type"."BP\n" if $fmate_type =~ /MS1/;
                                    print BED "$name/2\t$fmate_clip_st\t$fmate_clip_end\t$name"."_$fmate_type"."BP\n" if $fmate_type =~ /MS2/;
                                    $blast_query_count++;
                                }
                                if ($pair_type eq "clip_conc" and $smate_type =~ /^MS/) {
                                    #$ms_clip_count ++;
                                    #clipped part of smate wanted
                                    if ($smate_clipped >= $min_size_for_blast) {
                                        my $smate_clip_st;
                                        my $smate_clip_end;
                                        if ($smate_type =~ /\+/) {
                                            $smate_clip_st = ($read_length - $smate_clipped);
                                            $smate_clip_end = $read_length;
                                        } else {
                                            $smate_clip_st = 0;
                                            $smate_clip_end = $smate_clipped + 1;
                                        }
                                        print BED "$name/1\t$smate_clip_st\t$smate_clip_end\t$name"."_$smate_type"."BP\n" if $smate_type =~ /MS1/;
                                        print BED "$name/2\t$smate_clip_st\t$smate_clip_end\t$name"."_$smate_type"."BP\n" if $smate_type =~ /MS2/;
                                        $blast_query_count++;
                                    }
                                } elsif ($pair_type eq "clip_disc") {
                                    if ($smate_type =~ /^UN/) {
                                        #whole smate wanted
                                        print BED "$name/1\t1\t$read_length\t$name"."_$smate_type"."L\n" if $smate_type =~ /UN1/;
                                        print BED "$name/2\t1\t$read_length\t$name"."_$smate_type"."L\n" if $smate_type =~ /UN2/;
                                        $blast_query_count++;
                                    } elsif ($smate_type =~ /^MS/) {
                                        #match part of smate wanted
                                        if ($smate_end - $smate_st >= $min_size_for_blast - 1) {
                                            my $smate_match_st;
                                            my $smate_match_end;
                                            if ($smate_type =~ /\+/) {
                                                $smate_match_st = 0;
                                                $smate_match_end = ($read_length - $smate_clipped);
                                            } else {
                                                $smate_match_st = $smate_clipped;
                                                $smate_match_end = $read_length;
                                            }
                                            print BED "$name/1\t$smate_match_st\t$smate_match_end\t$name"."_$smate_type"."L\n" if $smate_type =~ /MS1/;
                                            print BED "$name/2\t$smate_match_st\t$smate_match_end\t$name"."_$smate_type"."L\n" if $smate_type =~ /MS2/;
                                            $blast_query_count++;
                                        }
                                    } elsif ($smate_type =~ /^SM/) {
                                        #match part of smate wanted
                                        if ($smate_end - $smate_st >= $min_size_for_blast - 1) {
                                            my $smate_match_st;
                                            my $smate_match_end;
                                            if ($smate_type =~ /\+/) {
                                                $smate_match_st = $smate_clipped;
                                                $smate_match_end = $read_length;
                                            } else {
                                                $smate_match_st = 0;
                                                $smate_match_end = ($read_length - $smate_clipped);
                                            }
                                            print BED "$name/1\t$smate_match_st\t$smate_match_end\t$name"."_$smate_type"."L\n" if $smate_type =~ /SM1/;
                                            print BED "$name/2\t$smate_match_st\t$smate_match_end\t$name"."_$smate_type"."L\n" if $smate_type =~ /SM2/;
                                            $blast_query_count++;
                                        }
                                    }
                                }
                            }

                            #get sm_clip mates or parts that lie within annotated TEs to identify the transposed one
                            my $sm_clip_data = shift @selected_sm;
                            foreach my $smc (split/,/,$sm_clip_data) {
                                my $name = (split/\t/,$smc) [0];
                                my ($fmate_chr, $fmate_st, $fmate_end, $fmate_type, $fmate_clipped) = (split/\t/,$smc) [1,2,3,4,5];
                                my ($smate_chr, $smate_st, $smate_end, $smate_type, $smate_clipped) = (split/\t/,$smc) [6,7,8,9,10];# smate is distant from bp if disc, at bp if conc
                                my $pair_type = (split/\t/,$smc) [11];
                                if ($fmate_clipped >= $min_size_for_blast) {
                                    #clipped part of fmate wanted
                                    my $fmate_clip_st;
                                    my $fmate_clip_end;
                                    if ($fmate_type =~ /\+/) {
                                        $fmate_clip_st = 0;
                                        $fmate_clip_end = $fmate_clipped + 1;
                                    } else {
                                        $fmate_clip_st = ($read_length - $fmate_clipped);
                                        $fmate_clip_end = $read_length;
                                    }
                                    print BED "$name/1\t$fmate_clip_st\t$fmate_clip_end\t$name"."_$fmate_type"."BP\n" if $fmate_type =~ /SM1/;
                                    print BED "$name/2\t$fmate_clip_st\t$fmate_clip_end\t$name"."_$fmate_type"."BP\n" if $fmate_type =~ /SM2/;
                                    $blast_query_count++;
                                }
                                if ($pair_type eq "clip_conc" and $smate_type =~ /^SM/) {
                                    #clipped part of smate wanted
                                    if ($smate_clipped >= $min_size_for_blast) {
                                        my $smate_clip_st;
                                        my $smate_clip_end;
                                        if ($smate_type =~ /\+/) {
                                            $smate_clip_st = 0;
                                            $smate_clip_end = $smate_clipped + 1;
                                        } else {
                                            $smate_clip_st = ($read_length - $smate_clipped);
                                            $smate_clip_end = $read_length;
                                        }
                                        print BED "$name/1\t$smate_clip_st\t$smate_clip_end\t$name"."_$smate_type"."BP\n" if $smate_type =~ /SM1/;
                                        print BED "$name/2\t$smate_clip_st\t$smate_clip_end\t$name"."_$smate_type"."BP\n" if $smate_type =~ /SM2/;
                                        $blast_query_count++;
                                    }
                                } elsif ($pair_type eq "clip_disc") {
                                    if ($smate_type =~ /^UN/) {
                                        #whole smate wanted
                                        print BED "$name/1\t0\t$read_length\t$name"."_$smate_type"."R\n" if $smate_type =~ /UN1/;
                                        print BED "$name/2\t0\t$read_length\t$name"."_$smate_type"."R\n" if $smate_type =~ /UN2/;
                                        $blast_query_count++;
                                    } elsif ($smate_type =~ /^SM/) {
                                        #match part of smate wanted
                                        if ($smate_end - $smate_st >= $min_size_for_blast - 1) {
                                            my $smate_match_st;
                                            my $smate_match_end;
                                            if ($smate_type =~ /\+/) {
                                                $smate_match_st = $smate_clipped;
                                                $smate_match_end = $read_length;
                                            } else {
                                                $smate_match_st = 0;
                                                $smate_match_end = ($read_length - $smate_clipped);
                                            }
                                            print BED "$name/1\t$smate_match_st\t$smate_match_end\t$name"."_$smate_type"."R\n" if $smate_type =~ /SM1/;
                                            print BED "$name/2\t$smate_match_st\t$smate_match_end\t$name"."_$smate_type"."R\n" if $smate_type =~ /SM2/;
                                            $blast_query_count++;
                                        }
                                    } elsif ($smate_type =~ /^MS/) {
                                        #match part of smate wanted
                                        if ($smate_end - $smate_st >= $min_size_for_blast - 1) {
                                            my $smate_match_st;
                                            my $smate_match_end;
                                            if ($smate_type =~ /\+/) {
                                                $smate_match_st = 0;
                                                $smate_match_end = ($read_length - $smate_clipped);
                                            } else {
                                                $smate_match_st = $smate_clipped + 0;
                                                $smate_match_end = $read_length;
                                            }
                                            print BED "$name/1\t$smate_match_st\t$smate_match_end\t$name"."_$smate_type"."R\n" if $smate_type =~ /MS1/;
                                            print BED "$name/2\t$smate_match_st\t$smate_match_end\t$name"."_$smate_type"."R\n" if $smate_type =~ /MS2/;
                                            $blast_query_count++;
                                        }
                                    }
                                }
                            }
                            foreach my $sel_ud (@selected_ud) {
                                my ($ud_data, $ud_side) = split/;/,$sel_ud;
                                #get unclip_disc mates that lie within annotated TEs to identify the transposed one
                                foreach my $ud (split/,/,$ud_data) {
                                    my ($ud_name, $ud_smate_chr, $ud_smate_st, $ud_smate_end, $ud_smate_type, $ud_smate_clipped) = (split/\t/,$ud) [0,6,7,8,9,10]; # smate is distant from bp
                                    if ($ud_smate_type =~ /^UN/) {
                                        print BED "$ud_name/1\t0\t$read_length\t$ud_name"."_$ud_smate_type"."R\n" if $ud_side eq "R" and $ud_smate_type =~ /UN1/;
                                        print BED "$ud_name/2\t0\t$read_length\t$ud_name"."_$ud_smate_type"."R\n" if $ud_side eq "R" and $ud_smate_type =~ /UN2/;
                                        print BED "$ud_name/1\t0\t$read_length\t$ud_name"."_$ud_smate_type"."L\n" if $ud_side eq "L" and $ud_smate_type =~ /UN1/;
                                        print BED "$ud_name/2\t0\t$read_length\t$ud_name"."_$ud_smate_type"."L\n" if $ud_side eq "L" and $ud_smate_type =~ /UN2/;
                                        $blast_query_count++;
                                    } elsif ($ud_smate_type =~ /^SM/) {
                                        #match part of smate wanted
                                        if ($ud_smate_end - $ud_smate_st >= $min_size_for_blast - 1) {
                                            my $ud_smate_match_st;
                                            my $ud_smate_match_end;
                                            if ($ud_smate_type =~ /\+/) {
                                                $ud_smate_match_st = $ud_smate_clipped;
                                                $ud_smate_match_end = $read_length;
                                            } else {
                                                $ud_smate_match_st = 0;
                                                $ud_smate_match_end = ($read_length - $ud_smate_clipped);
                                            }
                                            print BED "$ud_name/1\t$ud_smate_match_st\t$ud_smate_match_end\t$ud_name"."_$ud_smate_type"."R\n" if $ud_side eq "R" and $ud_smate_type =~ /SM1/;
                                            print BED "$ud_name/2\t$ud_smate_match_st\t$ud_smate_match_end\t$ud_name"."_$ud_smate_type"."R\n" if $ud_side eq "R" and $ud_smate_type =~ /SM2/;
                                            print BED "$ud_name/1\t$ud_smate_match_st\t$ud_smate_match_end\t$ud_name"."_$ud_smate_type"."L\n" if $ud_side eq "L" and $ud_smate_type =~ /SM1/;
                                            print BED "$ud_name/2\t$ud_smate_match_st\t$ud_smate_match_end\t$ud_name"."_$ud_smate_type"."L\n" if $ud_side eq "L" and $ud_smate_type =~ /SM2/;
                                            $blast_query_count++;
                                        }
                                    } elsif ($ud_smate_type =~ /^MS/) {
                                        #match part of smate wanted
                                        if ($ud_smate_end - $ud_smate_st >= $min_size_for_blast - 1) {
                                            my $ud_smate_match_st;
                                            my $ud_smate_match_end;
                                            if ($ud_smate_type =~ /\+/) {
                                                $ud_smate_match_st = 0;
                                                $ud_smate_match_end = ($read_length - $ud_smate_clipped);
                                            } else {
                                                $ud_smate_match_st = $ud_smate_clipped;
                                                $ud_smate_match_end = $read_length;
                                            }
                                            print BED "$ud_name/1\t$ud_smate_match_st\t$ud_smate_match_end\t$ud_name"."_$ud_smate_type"."R\n" if $ud_side eq "R" and $ud_smate_type =~ /MS1/;
                                            print BED "$ud_name/2\t$ud_smate_match_st\t$ud_smate_match_end\t$ud_name"."_$ud_smate_type"."R\n" if $ud_side eq "R" and $ud_smate_type =~ /MS2/;
                                            print BED "$ud_name/1\t$ud_smate_match_st\t$ud_smate_match_end\t$ud_name"."_$ud_smate_type"."L\n" if $ud_side eq "L" and $ud_smate_type =~ /MS1/;
                                            print BED "$ud_name/2\t$ud_smate_match_st\t$ud_smate_match_end\t$ud_name"."_$ud_smate_type"."L\n" if $ud_side eq "L" and $ud_smate_type =~ /MS2/;
                                            $blast_query_count++;
                                        }
                                    }
                                }
                            }
                            close BED;

                            #Looking for TEs that align (blast) with the disc mates or mates clipped parts
                            #getting the sequence of the read from the clip_$outname.fasta file using bedtools getfasta and blasting it against the TE database
                             system "$path_to_bedtools getfasta -fi $clipfasta -bed $chr:$sm_pos-$ms_pos.bed -fo $chr:$sm_pos-$ms_pos.fasta -name ;
                                    $path_to_blast -task blastn -query $chr:$sm_pos-$ms_pos.fasta -db $TEdb -out $chr:$sm_pos-$ms_pos.bl -evalue $min_evalue -perc_identity $min_te_pident -max_target_seqs 1 -num_threads $num_threads -outfmt '6 qseqid qstart qend qlen sseqid sstart send slen pident evalue'";

                            unlink "$chr:$sm_pos-$ms_pos.bed","$chr:$sm_pos-$ms_pos.fasta";

                            #parsing the blast output
                            my %te_data;
                            open IN, "$chr:$sm_pos-$ms_pos.bl";

                            while (<IN>) {
                                chomp;
                                my ($read_desc, $q_st, $q_end, $q_len, $te_name, $hit_st, $hit_end, $te_length, $pc_ident) = (split/\t/)[0,1,2,3,4,5,6,7,8];
                                next if $q_st != 1 and $q_end != $q_len; # no alignment of only internal part of read allowed, one extremity at least must be engaged
                                my $hpos = $hit_st;
                                $hpos = $hit_end if $hit_st > $hit_end; #hsp coordinates are function of the orientation of the alignment
                                #defining the ext limits for the alignment
                                my $in_limit_left = $dist_from_bord;
                                my $in_limit_right = $te_length - $dist_from_bord;
                                $in_limit_left = $te_length and $in_limit_right = 1 if $te_length < $dist_from_bord;
                                # some alignments should start near each extremity of the TE
                                $te_data{$te_name} .= "extL," if ($hpos <= $in_limit_left);
                                $te_data{$te_name} .= "extR," if ($hpos >= $in_limit_right);
                                my ($read_name, $read_type) = split/_/,$read_desc;
                                if ($read_type =~ /MS.{2}BP/) {
                                    $te_data{$te_name} .= "msBP,";
                                } elsif ($read_type =~ /SM.{2}BP/) {
                                    $te_data{$te_name} .= "smBP,";
                                } elsif ($read_type =~ /[^ext]L/) {
                                    $te_data{$te_name} .= "L,";
                                } elsif ($read_type =~ /[^ext]R/) {
                                    $te_data{$te_name} .= "R,";
                                }
                            }
                            close IN;

                            # final selection
                            my $selected_te;
                            my $test = "NO";
                            foreach my $te_name (keys %te_data) {
                                my @atts = split /,/, $te_data{$te_name};
                                my $smr_score = 0; # nb of sm or R reads
                                my $msl_score = 0; # nb of ms or L reads
                                my $extL_score = 0; # nb reads matching the left of the TE extremity
                                my $extR_score = 0; # nb reads matching the right of the TE extremity
                                for (@atts) {
                                    $msl_score++ if $_ ne "extL" and ($_ =~ /L/ or $_ eq "msBP");
                                    $smr_score++ if $_ ne "extR" and ($_ =~ /R/ or $_ eq "smBP");
                                    $extL_score++ if $_ eq "extL";
                                    $extR_score++ if $_ eq "extR";
                                }
                                if ($msl_score > 0
                                and $smr_score > 0
                                and $msl_score + $smr_score >= $blast_query_count * $min_TE_pcent
                                and $extL_score > 0
                                and $extR_score > 0
                                ){
                                    $test = "YES";
                                    $selected_te .= "$te_name";
                                }
                            }
                            if ($test eq "YES") {
                                print OUT1 "$chr\t$sm_pos\t$ms_pos\t$selected_te\n";
                            } else {
                                unlink "$chr:$sm_pos-$ms_pos.bl";
                            }
                        }
                    }
                    last;
                } # end if (exists $SM_clip{$chr}{$sm_pos})
            } # end sm_pos
        }
    }# end ms_pos
    print "\n$chr done";
}# end chr
print "\n";
close OUT1;
