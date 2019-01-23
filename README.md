# sam2ins

Perl scripts to uncover transposable elements insertions from SAM output

### General guidelines
The starting point is a SAM file produced by bowtie2. This SAM file must contain a header and must be sorted by queryname  
Each perl script (name.pl) is provided with its companion configuration file (name.conf)  
The configuration files contain all the parameters (default parameters are indicated, specific ones should be added) not indicated in the command line for the corresponding script to work  
The only 2 arguments in the command line are the perl script and the main file to parse  
All scripts and configuration files must be present in the same working directory

### Step 1 : sam2clip.pl
Command line  : `perl sam2clip.pl SAM_file`  
Input : SAM file from bowtie2 analysis  
Output : clip data file (TSV format [see description below],T0_sam2clip_chr01.tsv is provided as an example), SAM file and FASTA file describing the clipped reads  
Note : the SAM output is not used in next step but may be useful after step 2 for visual inspection of the insertions in IGV or similar genome browser

The 12 columns clip data file contains :
- col 1 : name of the read pair
- col 2-6 : data describing the alignment of the first mate as follows :
- col 2,3,4 : coordinates of the alignment on the chromosomes
- col 5 : code describing the alignment (see below)) 
- col 6 : number of clipped nucleotides
- col 7-11 : data describing the alignment of the second mate as above
- col 12 : description of the aligned pair (eg : "clip_conc" means : concordant alignment, at least one mate clipped)
Code describing the alignements (col 5 and 10) : eg, SM2+ means that the left part of the read is soft-clipped (S), the right part is aligned (M for Match), the read is number 2 in the initial pair and the read is aligned into its initial orientation (+). UN means unclipped.

### Step 2 : clip2ins.pl
Command line : `clip2ins.pl clip_data_file`  
Input : clip data file from step 1  
Output : BED file describing the insertions (T0_clip2ins.bed is provided as an example)  
Additional necessary software : BLAST, BEDTOOLS  
Additional necessary data : nls.bed (see below), TE database in FASTA format, FASTA file from step 1
 
### Optional step : selfclips.pl
Command line : `perl selfclips.pl clip_data_file`  
Input : clip data file from step 1  
Output : BED file listing the false positive breakpoints (nls.bed)  
Note : A nls.bed file is necessary for clip2ins.pl to work but this step is optional because a nls.bed file is provided
