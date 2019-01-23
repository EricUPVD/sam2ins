#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

#getting the sam and the configuration file
die "Please provide a sam file. Command line is : perl sam2clip.pl <sam file> " if @ARGV < 1;
die "Please provide the sam2clip.conf file in the working directory" unless -e "sam2clip.conf";
my $infile = shift;
my $outname = fileparse($infile, (".sam"));

#getting and checking the parameters from the configuration file
my $read_length;
my $max_for_conc; # max distance for concordant alignment
my $max_edit_dist; # max nb of mismatches allowed for a alignment to be considered
my $min_clipped_nb; # min number of clipped nt for a read to be considered as clipped
open CONF, "sam2clip.conf";
while (<CONF>) {
    chomp;
    $read_length = (split /=/)[1] if $_ =~ /^READ_LENGTH/;
    $max_for_conc = (split /=/)[1] if $_ =~ /^MAX_FOR_CONC/;
    $max_edit_dist = (split /=/)[1] if $_ =~ /^MAX_EDIT_DIST/;
    $min_clipped_nb = (split /=/)[1] if $_ =~ /^MIN_CLIP_NB/;
}
close CONF;
die "Please provide a value for READ_LENGTH" if $read_length !~ /\d+/;
die "Please provide a value for MAX_FOR_CONC" if $max_for_conc !~ /\d+/;
die "Please provide a value for MAX_EDIT_DIST" if $max_edit_dist !~ /\d+/;
die "Please provide a value for MIN_CLIP_NB" if $min_clipped_nb !~ /\d+/;
print "Parameters : OK\n";

my %disc_flag = (81 => "-,+", 97 => "+,-", 65 => "+,+", 113 => "-,-"); # flag information is only necessary for mate1
my %conc_flag = (99 => "+,-", 83 => "-,+");

# The initial sam file may be to big to be stored in a single hash because of memory limits.
# Therefore, one hash is made for each chromosome and the initial sam file re-parsed each time.
# The outputs are gathered into only one tsv file and one sam file;

my @ref_seqs; # Store the names of the ref_seqs;
open OUT1, ">", "clip_".$outname.".tsv";
open OUT2, ">", "clip_".$outname.".sam";

#Parsing the header of the sam file to get the names of the ref_seqs
open IN, $infile;
my $j = 0;
while (<IN>) {
    chomp;
    $j++;
    die "Please provide a sam file with a header" if ($j==1 and $_ !~ /^@/);
    if ($_ =~ /^@/) {
        die "Please provide a sam file sorted by query name" if ($_ =~ /SO:/ and $_ !~ /queryname/);
        print OUT2 "$_\n";
        if ($_ =~ /\SQ/) {
            my $ref_seq_data = (split/\t/)[1];
            my $ref = (split/:/,$ref_seq_data)[1];
            push @ref_seqs, $ref;
        }
    } else {
        last;
    }
}
close IN;
print "Ref Seq names stored\n";

# parsing the sam file chr by chr
foreach my $ref (@ref_seqs) {
    my %sam_data;
    open IN, $infile;
    my $i = 0;
    my $keep = "No";
    while (<IN>) {
        chomp;
        next if ($_ =~ /^@/);
        $i++;
        my ($read_name, $read_chr)  = (split/\t/)[0,2];
        if ($i % 2 == 1) {              # odd lines contain data for mate1
            if ($read_chr eq $ref) {    # get mate1
                $sam_data{$read_name} .= "$_~";
                $keep = "Yes";
            } else {
                $keep = "No";
                next;
            }
        } elsif ($i % 2 == 0) {         # even lines contain data for mate2
            if ($keep eq "Yes") {       # get mate2
                $sam_data{$read_name} .= "$_~";
            } elsif ($keep eq "No") {
                next;
            }
        }
    }
    close IN;
    print "$ref data stored, analysis begins\n";

    foreach my $read_name (sort keys %sam_data) {
        my ($mate1_data, $mate2_data) = split /~/, $sam_data{$read_name};
        my @mate1_datalist = split/\t/,$mate1_data;
        my ($mate1_chr, $mate1_pos, $mate1_flag, $mate1_cigar) = @mate1_datalist[2,3,1,5];
        next if $mate1_cigar eq "*"; # ignoring pairs with unmapped read
        my @mate2_datalist = split/\t/,$mate2_data;
        my ($mate2_chr, $mate2_pos, $mate2_flag, $mate2_cigar) = @mate2_datalist[2,3,1,5];
        next if $mate2_cigar eq "*"; # ignoring pairs with unmapped read
        my $nm_tag = 0;
        while (my $tag = pop @mate1_datalist) {
            if ($tag =~/^NM/) {
                $nm_tag = (split/:/,$tag)[2];
                last;
            }
        }
        next if $nm_tag > $max_edit_dist; #ignoring pairs if one read displays too much mismatches
        while (my $tag = pop @mate2_datalist) {
            if ($tag =~/^NM/) {
                $nm_tag = (split/:/,$tag)[2];
                last;
            }
        }
        next if $nm_tag > $max_edit_dist; #ignoring pairs if one read displays too much mismatches

        #looking for reads belonging to CONCORDANT pairs based on flag
        foreach my $cf (keys %conc_flag) {
            if ($mate1_flag == $cf) {
                my ($mate1_orient, $mate2_orient) = (split/,/,$conc_flag{$cf});

                #looking for clipped reads
                my $clipped_nb1 = 0;
                my $clipped_nb2 = 0;

                #mate 1 is clipped and mate 2 is unclipped
                if ($mate1_cigar =~ /S/ and $mate2_cigar !~ /S/) {
                    my $mate2_type = "UN2";
                    my $mate1_type;
                    #MS type
                    if ($mate1_cigar =~ /^\d{2,}M/ and $mate1_cigar =~ /(\d{1,})S$/ ) {
                        $mate1_type = "MS1";
                        $clipped_nb1 = $1;
                    #SM type
                    } elsif ($mate1_cigar =~ /\d{2,}M$/ and $mate1_cigar =~ /^(\d{1,})S/) {
                        $mate1_type = "SM1";
                        $clipped_nb1 = $1;
                    } else { #some reads have a xxSxxMxxS pattern
                        last;
                    }
                    #SAM output for selected reads
                    print OUT2 "$mate1_data\n$mate2_data\n";
                    #A sam file reports the leftmost position of the alignment. The end is calculated based on the nb of
                    #matching nucleotides which depends on whether the read is clipped or not;
                    my $mate1_end = ($mate1_pos + ($read_length - $clipped_nb1)); # calculation is the same whether it is a MS or SM clip
                    $mate1_data = "$mate1_chr\t$mate1_pos\t$mate1_end\t$mate1_type$mate1_orient\t$clipped_nb1";
                    my $mate2_end = $mate2_pos + $read_length;
                    $mate2_data = "$mate2_chr\t$mate2_pos\t$mate2_end\t$mate2_type$mate2_orient\t$clipped_nb2";
                    #clipped mate is printed first, a min nb of clipped bases being necessary
                    print OUT1 "$read_name\t$mate1_data\t$mate2_data\tclip_conc\n" if $clipped_nb1 >= $min_clipped_nb;

                #mate 2 is clipped and mate 1 is unclipped
                } elsif ($mate2_cigar =~ /S/ and $mate1_cigar !~ /S/) {
                    my $mate1_type = "UN1";
                    my $mate2_type;
                    #MS type
                    if ($mate2_cigar =~ /^\d{2,}M/ and $mate2_cigar =~ /(\d{1,})S$/) {
                        $mate2_type = "MS2";
                        $clipped_nb2 = $1;
                    #SM type
                    } elsif ($mate2_cigar =~ /\d{2,}M$/ and $mate2_cigar =~ /^(\d{1,})S/) {
                        $mate2_type = "SM2";
                        $clipped_nb2 = $1;
                    } else { #some reads have a xxSxxMxxS pattern
                        last;
                    }
                    #SAM output for selected reads
                    print OUT2 "$mate1_data\n$mate2_data\n";
                    #A sam file reports the leftmost position of the alignment. The end is calculated based on the nb of
                    #matching nucleotides which depends on whether the read is clipped or not;
                    my $mate2_end = ($mate2_pos + ($read_length - $clipped_nb2)); # calculation is the same whether it is a MS or SM clip
                    $mate2_data = "$mate2_chr\t$mate2_pos\t$mate2_end\t$mate2_type$mate2_orient\t$clipped_nb2";
                    my $mate1_end = $mate1_pos + $read_length;
                    $mate1_data = "$mate1_chr\t$mate1_pos\t$mate1_end\t$mate1_type$mate1_orient\t$clipped_nb1";
                    #clipped mate is printed first, a min nb of clipped bases being necessary
                    print OUT1 "$read_name\t$mate2_data\t$mate1_data\tclip_conc\n" if $clipped_nb2 >= $min_clipped_nb ;

                #both mate1 and mate2 are clipped : if concordant flag, they are both MS or both SM
                } elsif ($mate1_cigar =~ /S/ and $mate2_cigar =~ /S/) {
                    #both MS type
                    my $mate1_type;
                    my $mate2_type;
                    if ($mate1_cigar =~ /^\d{2,}M/ and $mate1_cigar =~ /(\d{1,})S$/) {
                        $mate1_type = "MS1";
                        $clipped_nb1 = $1;
                        if ($mate2_cigar =~ /^\d{2,}M/ and $mate2_cigar =~ /(\d{1,})S$/) {
                            $mate2_type = "MS2";
                            $clipped_nb2 = $1;
                        }
                    #both SM type
                    } elsif ($mate1_cigar =~ /\d{2,}M$/ and $mate1_cigar =~ /^(\d{1,})S/) {
                        $mate1_type = "SM1";
                        $clipped_nb1 = $1;
                        if ($mate2_cigar =~ /\d{2,}M$/ and $mate2_cigar =~ /^(\d{1,})S/) {
                            $mate2_type = "SM2";
                            $clipped_nb2 = $1;
                        }
                    } else { #some reads have a xxSxxMxxS pattern
                        last;
                    }
                    if ($clipped_nb1 > 0 and $clipped_nb2 > 0) {
                        #SAM output for selected reads
                        print OUT2 "$mate1_data\n$mate2_data\n";
                        #A sam file reports the leftmost position of the alignment. The end is calculated based on the nb of
                        #matching nucleotides which depends on whether the read is clipped or not;
                        my $mate1_end = ($mate1_pos + ($read_length - $clipped_nb1)); # calculation is the same whether it is a MS or SM clip
                        $mate1_data = "$mate1_chr\t$mate1_pos\t$mate1_end\t$mate1_type$mate1_orient\t$clipped_nb1";
                        my $mate2_end = ($mate2_pos + ($read_length - $clipped_nb2)); # calculation is the same whether it is a MS or SM clip
                        $mate2_data = "$mate2_chr\t$mate2_pos\t$mate2_end\t$mate2_type$mate2_orient\t$clipped_nb2";
                        #both are clipped, both are at the bp : only one record is necessary
                        print OUT1 "$read_name\t$mate1_data\t$mate2_data\tclip_conc\n";
                    }
                }
                last;
            }
        }

        #looking for reads belonging to DISCORDANT pairs based on flag
        foreach my $df (keys %disc_flag) {
            if ($mate1_flag == $df) {
                my ($mate1_orient, $mate2_orient) = (split/,/,$disc_flag{$df});

                #looking for reads belonging to CONCORDANT pairs despite of a discordant flag based on size of the template
                my $clipped_nb1 = 0;
                my $clipped_nb2 = 0;
                if ($mate2_chr eq $mate1_chr and abs($mate1_pos - $mate2_pos + $read_length + 1) <= $max_for_conc) {

                    #mate 1 is clipped and mate 2 is unclipped
                    if ($mate1_cigar =~ /S/ and $mate2_cigar !~ /S/) {
                        my $mate2_type = "UN2";
                        my $mate1_type;
                    #MS type
                    if ($mate1_cigar =~ /^\d{2,}M/ and $mate1_cigar =~ /(\d{1,})S$/ ) {
                        $mate1_type = "MS1";
                        $clipped_nb1 = $1;
                    #SM type
                    } elsif ($mate1_cigar =~ /\d{2,}M$/ and $mate1_cigar =~ /^(\d{1,})S/) {
                        $mate1_type = "SM1";
                        $clipped_nb1 = $1;
                    } else { #some reads have a xxSxxMxxS pattern
                        last;
                    }
                    #SAM output for selected reads
                    print OUT2 "$mate1_data\n$mate2_data\n";
                    #A sam file reports the leftmost position of the alignment. The end is calculated based on the nb of
                    #matching nucleotides which depends on whether the read is clipped or not;
                    my $mate1_end = ($mate1_pos + ($read_length - $clipped_nb1)); # calculation is the same whether it is a MS or SM clip
                    $mate1_data = "$mate1_chr\t$mate1_pos\t$mate1_end\t$mate1_type$mate1_orient\t$clipped_nb1";
                    my $mate2_end = $mate2_pos + $read_length;
                    $mate2_data = "$mate2_chr\t$mate2_pos\t$mate2_end\t$mate2_type$mate2_orient\t$clipped_nb2";
                    #clipped mate is printed first, a min nb of clipped bases being necessary
                    print OUT1 "$read_name\t$mate1_data\t$mate2_data\tclip_conc\n" if $clipped_nb1 >= $min_clipped_nb;

                    #mate 2 is clipped and mate 1 is unclipped
                    } elsif ($mate2_cigar =~ /S/ and $mate1_cigar !~ /S/) {
                        my $mate1_type = "UN1";
                        my $mate2_type;
                        #MS type
                        if ($mate2_cigar =~ /^\d{2,}M/ and $mate2_cigar =~ /(\d{1,})S$/) {
                            $mate2_type = "MS2";
                            $clipped_nb2 = $1;
                        #SM type
                        } elsif ($mate2_cigar =~ /\d{2,}M$/ and $mate2_cigar =~ /^(\d{1,})S/) {
                            $mate2_type = "SM2";
                            $clipped_nb2 = $1;
                        } else { #some reads have a xxSxxMxxS pattern
                            last;
                        }
                        #SAM output for selected reads
                        print OUT2 "$mate1_data\n$mate2_data\n";
                        #A sam file reports the leftmost position of the alignment. The end is calculated based on the nb of
                        #matching nucleotides which depends on whether the read is clipped or not;
                        my $mate2_end = ($mate2_pos + ($read_length - $clipped_nb2)); # calculation is the same whether it is a MS or SM clip
                        $mate2_data = "$mate2_chr\t$mate2_pos\t$mate2_end\t$mate2_type$mate2_orient\t$clipped_nb2";
                        my $mate1_end = $mate1_pos + $read_length;
                        $mate1_data = "$mate1_chr\t$mate1_pos\t$mate1_end\t$mate1_type$mate1_orient\t$clipped_nb1";
                        #clipped mate is printed first, a min nb of clipped bases being necessary
                        print OUT1 "$read_name\t$mate2_data\t$mate1_data\tclip_conc\n" if $clipped_nb2 >= $min_clipped_nb;

                    #both mate1 and mate2 are clipped : if concordant flag, they are both MS or both SM
                    } elsif ($mate1_cigar =~ /S/ and $mate2_cigar =~ /S/) {
                        #both MS type
                        my $mate1_type;
                        my $mate2_type;
                        if ($mate1_cigar =~ /^\d{2,}M/ and $mate1_cigar =~ /(\d{1,})S$/) {
                            $mate1_type = "MS1";
                            $clipped_nb1 = $1;
                            if ($mate2_cigar =~ /^\d{2,}M/ and $mate2_cigar =~ /(\d{1,})S$/) {
                                $mate2_type = "MS2";
                                $clipped_nb2 = $1;
                            }
                        #both SM type
                        } elsif ($mate1_cigar =~ /\d{2,}M$/ and $mate1_cigar =~ /^(\d{1,})S/) {
                            $mate1_type = "SM1";
                            $clipped_nb1 = $1;
                            if ($mate2_cigar =~ /\d{2,}M$/ and $mate2_cigar =~ /^(\d{1,})S/) {
                                $mate2_type = "SM2";
                                $clipped_nb2 = $1;
                            }
                        } else { #some reads have a xxSxxMxxS pattern
                            last;
                        }
                        if ($clipped_nb1 > 0 and $clipped_nb2 > 0) {
                            #SAM output for selected reads
                            print OUT2 "$mate1_data\n$mate2_data\n";
                            #A sam file reports the leftmost position of the alignment. The end is calculated based on the nb of
                            #matching nucleotides which depends on whether the read is clipped or not;
                            my $mate1_end = ($mate1_pos + ($read_length - $clipped_nb1)); # calculation is the same whether it is a MS or SM clip
                            $mate1_data = "$mate1_chr\t$mate1_pos\t$mate1_end\t$mate1_type$mate1_orient\t$clipped_nb1";
                            my $mate2_end = ($mate2_pos + ($read_length - $clipped_nb2)); # calculation is the same whether it is a MS or SM clip
                            $mate2_data = "$mate2_chr\t$mate2_pos\t$mate2_end\t$mate2_type$mate2_orient\t$clipped_nb2";
                            #both are clipped, both are at the bp : only one record is necessary
                            print OUT1 "$read_name\t$mate1_data\t$mate2_data\tclip_conc\n";
                        }
                    }
                } else {
                    #looking for reads belonging to REAL DISCORDANT pairs based on size of the template
                    #unclipped reads from discordant pairs to reinforce identification of the TE
                    if ($mate1_cigar eq "$read_length"."M" and $mate2_cigar eq "$read_length"."M" ) {
                        my $mate1_type = "UN1";
                        my $mate2_type = "UN2";
                        #SAM output for selected reads
                        print OUT2 "$mate1_data\n$mate2_data\n";
                        my $mate1_end = $mate1_pos + $read_length;
                        $mate1_data = "$mate1_chr\t$mate1_pos\t$mate1_end\t$mate1_type$mate1_orient\t$clipped_nb1";
                        my $mate2_end = $mate2_pos + $read_length;
                        $mate2_data = "$mate2_chr\t$mate2_pos\t$mate2_end\t$mate2_type$mate2_orient\t$clipped_nb2";
                        # don't know which mate will be near the bp : mate1 and mate2 are then both printed first
                        print OUT1 "$read_name\t$mate1_data\t$mate2_data\tunclip_disc\n";
                        print OUT1 "$read_name\t$mate2_data\t$mate1_data\tunclip_disc\n";

                    #clipped reads from discordant pairs to both reinforce the localization of the insertion site and to identify the TE
                    #mate 1 is clipped and mate 2 is unclipped
                    } elsif ($mate1_cigar =~ /S/ and $mate2_cigar !~ /S/) {
                        my $mate2_type = "UN2";
                        my $mate1_type;
                        #MS type
                        if ($mate1_cigar =~ /^\d{2,}M/ and $mate1_cigar =~ /(\d{1,})S$/ ) {
                            $mate1_type = "MS1";
                            $clipped_nb1 = $1;
                        #SM type
                        } elsif ($mate1_cigar =~ /\d{2,}M$/ and $mate1_cigar =~ /^(\d{1,})S/) {
                            $mate1_type = "SM1";
                            $clipped_nb1 = $1;
                        } else { #some reads have a xxSxxMxxS pattern
                            last;
                        }
                        #SAM output for selected reads
                        print OUT2 "$mate1_data\n$mate2_data\n";
                        #A sam file reports the leftmost position of the alignment. The end is calculated based on the nb of
                        #matching nucleotides which depends on whether the read is clipped or not;
                        my $mate1_end = ($mate1_pos + ($read_length - $clipped_nb1)); # calculation is the same whether it is a MS or SM clip
                        $mate1_data = "$mate1_chr\t$mate1_pos\t$mate1_end\t$mate1_type$mate1_orient\t$clipped_nb1";
                        my $mate2_end = $mate2_pos + $read_length;
                        $mate2_data = "$mate2_chr\t$mate2_pos\t$mate2_end\t$mate2_type$mate2_orient\t$clipped_nb2";
                        #clipped mate is printed first, a min nb of clipped bases being necessary
                        print OUT1 "$read_name\t$mate1_data\t$mate2_data\tclip_disc\n"if $clipped_nb1 >= $min_clipped_nb;
                        #tagged as unclip_disc if a min nb of clipped bases is not reached
                        print OUT1 "$read_name\t$mate1_data\t$mate2_data\tunclip_disc\n"if $clipped_nb1 < $min_clipped_nb;
                        # to detect locally unclipped mate with distant clipped mate
                        print OUT1 "$read_name\t$mate2_data\t$mate1_data\tunclip_disc\n";

                    #mate 2 is clipped and mate 1 is unclipped
                    } elsif ($mate2_cigar =~ /S/ and $mate1_cigar !~ /S/) {
                        my $mate1_type = "UN1";
                        my $mate2_type;
                        #MS type
                        if ($mate2_cigar =~ /^\d{2,}M/ and $mate2_cigar =~ /(\d{1,})S$/) {
                            $mate2_type = "MS2";
                            $clipped_nb2 = $1;
                        #SM type
                        } elsif ($mate2_cigar =~ /\d{2,}M$/ and $mate2_cigar =~ /^(\d{1,})S/) {
                            $mate2_type = "SM2";
                            $clipped_nb2 = $1;
                        } else { #some reads have a xxSxxMxxS pattern
                            last;
                        }
                        #SAM output for selected reads
                        print OUT2 "$mate1_data\n$mate2_data\n";
                        #A sam file reports the leftmost position of the alignment. The end is calculated based on the nb of
                        #matching nucleotides which depends on whether the read is clipped or not;
                        my $mate2_end = ($mate2_pos + ($read_length - $clipped_nb2)); # calculation is the same whether it is a MS or SM clip
                        $mate2_data = "$mate2_chr\t$mate2_pos\t$mate2_end\t$mate2_type$mate2_orient\t$clipped_nb2";
                        my $mate1_end = $mate1_pos + $read_length;
                        $mate1_data = "$mate1_chr\t$mate1_pos\t$mate1_end\t$mate1_type$mate1_orient\t$clipped_nb1";
                        #clipped mate is printed first, a min nb of clipped bases being necessary
                        print OUT1 "$read_name\t$mate2_data\t$mate1_data\tclip_disc\n" if $clipped_nb2 >= $min_clipped_nb;
                        #tagged as unclip_disc if a min nb of clipped bases is not reached
                        print OUT1 "$read_name\t$mate2_data\t$mate1_data\tunclip_disc\n" if $clipped_nb2 < $min_clipped_nb;
                        # to detect locally unclipped mate with distant clipped mate
                        print OUT1 "$read_name\t$mate1_data\t$mate2_data\tunclip_disc\n";

                    #both are clipped : when discordant if one is MS, the other can be MS or SM
                    } elsif ($mate1_cigar =~ /S/ and $mate2_cigar =~ /S/) {
                        my $mate1_type;
                        my $mate2_type;
                        # mate1 MS
                        if ($mate1_cigar =~ /^\d{2,}M/ and $mate1_cigar =~ /(\d{1,})S$/) {
                            $mate1_type = "MS1";
                            $clipped_nb1 = $1;
                            # mate2 MS
                            if ($mate2_cigar =~ /^\d{2,}M/ and $mate2_cigar =~ /(\d{1,})S$/) {
                                $mate2_type = "MS2";
                                $clipped_nb2 = $1;
                            #mate2 SM
                            } elsif ($mate2_cigar =~ /\d{2,}M$/ and $mate2_cigar =~ /^(\d{1,})S/){
                                $mate2_type = "SM2";
                                $clipped_nb2 = $1;
                            }
                        # mate1 SM
                        } elsif ($mate1_cigar =~ /\d{2,}M$/ and $mate1_cigar =~ /^(\d{1,})S/) {
                            $mate1_type = "SM1";
                            $clipped_nb1 = $1;
                            # mate2 SM
                            if ($mate2_cigar =~ /\d{2,}M$/ and $mate2_cigar =~ /^(\d{1,})S/) {
                                $mate2_type = "SM2";
                                $clipped_nb2 = $1;
                            # mate2 MS
                            } elsif ($mate2_cigar =~ /^\d{2,}M/ and $mate2_cigar =~ /(\d{1,})S$/){
                                $mate2_type = "MS2";
                                $clipped_nb2 = $1;
                            }
                        } else { #some reads have a xxSxxMxxS pattern
                            last;
                        }
                        if ($clipped_nb1 > 0 and $clipped_nb2 > 0) {
                            #SAM output for selected reads
                            print OUT2 "$mate1_data\n$mate2_data\n";
                            #A sam file reports the leftmost position of the alignment. The end is calculated based on the nb of
                            #matching nucleotides which depends on whether the read is clipped or not;
                            my $mate1_end = ($mate1_pos + ($read_length - $clipped_nb1)); # calculation is the same whether it is a MS or SM clip
                            $mate1_data = "$mate1_chr\t$mate1_pos\t$mate1_end\t$mate1_type$mate1_orient\t$clipped_nb1";
                            my $mate2_end = ($mate2_pos + ($read_length - $clipped_nb2)); # calculation is the same whether it is a MS or SM clip
                            $mate2_data = "$mate2_chr\t$mate2_pos\t$mate2_end\t$mate2_type$mate2_orient\t$clipped_nb2";
                            #both are clipped, don't know which mate will be near the bp : both mates are printed first
                            print OUT1 "$read_name\t$mate1_data\t$mate2_data\tclip_disc\n";
                            print OUT1 "$read_name\t$mate2_data\t$mate1_data\tclip_disc\n";
                        }
                    }
                }
                last;
            }
        }
    }
}
close OUT1;
close OUT2;
