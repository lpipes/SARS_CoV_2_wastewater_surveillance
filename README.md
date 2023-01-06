# Method for estimating relative proportions of SARS-CoV-2 strains from wastewater samples
We present a method for estimating the relative proportions of SARS-CoV-2 strains from wastewater samples. The method uses an initial step to remove unlikely strains, imputation of missing nucleotides using the global SARS-CoV-2 phylogeny, and an Expectation-Maximization (EM) algorithm for obtaining maximum likelihood estimates of the proportions of different strains in a sample.

View the manuscript here: <a href="https://www.sciencedirect.com/science/article/pii/S266723752200203X">https://www.sciencedirect.com/science/article/pii/S266723752200203X</a>

The method has 2 different components: estimating proportions of SARS-CoV-2 strains and imputation of SARS-CoV-2 reference strains. To estimate proportions of SARS-CoV-2 strains from short-read sequencing data on a pre-built database of 5,254,113 non-redundant strains from November 2, 2022. Download the required database associated files:

Imputed Multiple Sequence Alignment (-i): <a href="http://149.165.153.149/download_gz/seqs_final_out.fasta.gz">seqs_final_out.fasta.gz<a/> [34GB]

Variants file (-v): <a href="http://149.165.153.149/download_variants/variants.txt">variants.txt</a>

MSA reference file (-g): <a href="http://149.165.153.149/download_seq/EPI_ISL_402124.fasta">EPI_ISL_402124.fasta</a>

<a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">bowtie2</a> must also be installed and in your path.

Download a list of the redundant strains in the database (labeled c1, c2, c3, etc.): <a href="http://149.165.153.149/download_redundant/redundant.txt">redundant.txt</a>
 
# Installation
To install
	
	git clone https://github.com/lpipes/SARS_CoV_2_wastewater_surveillance.git
	cd SARS_CoV_2_wastewater_surveillance/eliminate_strains
	make
	cd ../imputation
	make

# Cleaning your reads
To run the method. First, please quality filter the reads. The method requires high quality reads as the variant calling in the method is sensitive to bad calls. I recommend the FASTX_toolkit (https://github.com/agordon/fastx_toolkit) in this way with your reads (quality filtering as well as removing DNA damage at both 5' and 3' ends of reads):
	
	#First quality filter reads
	fastq_quality_trimmer -v -t 35 -i reads.fastq -o reads_trimmed1.fastq -Q33
	#Second remove 15 bases from the end of reads and remove sequences with length less than 65bp
	fastx_trimmer -m 65 -t 15 -i reads_trimmed1.fastq -o reads_trimmed2.fastq
	#Third remove 15 bases from the start of the reads
	fastx_trimmer -f 15 -i reads_trimmed2.fastq -o reads_trimmed3.fastq

Running `-d` will clean your reads in this way. Please make sure fastq_quality_trimmer is installed and in your path.

# Usage
`eliminate_strains` filters unlikely SARS-CoV-2 genomes, prints a mismatch matrix, and also runs the `EM_C_LLR.R` program, which needs to be in your path.

	eliminate_strains [OPTIONS]
	
		-h, --help				
		-i, --infile [REQUIRED,FILE]		MSA FASTA of SARS-CoV-2 reference strains
		-s, --samfile [REQUIRED,FILE]		output sam file to print alignments
		-f, --freq [REQUIRED,decimal]		allele frequency to filter unlikely strains [default: 0.01]
		-o, --outfile [REQUIRED,FILE]		output file to print mismatch matrix for EM algorithm
		-v, --variant_sites [REQUIRED,FILE]	list of variant sites
		-g, --msa-reference [REQUIRED,FILE]	MSA reference index
		-p, --paired				using paired-reads
		-0, --single_end_file [FILE]		single-end reads
		-1, --forward_file [FILE]		if using paired-reads, the forward reads file
		-2, --reverse_file [FILE]		if using paired-reads, the reverse reads file
		-e, --EM-error [decimal]		error rate for EM algorithm
		-d, --clean-my-reads                    Clean reads with fastq_quality_trimmer [must have FASTQ reads]
		-c, --coverage [integer]		number of reads needed to calculate allele freq [default: 50]
		-a, --fasta				reads are in FASTA format [default: FASTQ]
		-l, --llr				Perform the LLR procedure
		-m, --min [decimal]			Minimum strains remaining to invoke iterative procedure [default: 100]
		-x, --max [decimal]			Maximum strains remaining for EM algorithm [default: 10000]
		-b, --print-allele-counts [FILE]	Print allele counts to file
		-t, --cores [decimal]			Number of cores [default: 1]
		-n, --no-read-sam			Don't thread, don't read in sam file to memory
		-r, --print-deletions [FILE]		Print sites with deletions
		-j, --threshold-for-deleted-sites	Threshold to print deleted sites [default: 0.001]

To use the bowtie2 database compatible with the pre-built database use the following bowtie2 database:
	
	eliminate_strains/MN908947.3.fasta

Example command to run eliminate_strains with single-end FASTQ reads (with a maximum of 35,000 strains remaining after elimination):

	eliminate_strains -i seqs_final_out.fasta.gz -s alignment.sam -f 0.01 -o mismatch.txt -v variants.txt -0 reads_trimmed3.fastq -e 0.005 -g EPI_ISL_402124.fasta -x 35000

Example command to run eliminate_strains with paired-end FASTQ reads (with a maximum of 35,000 strains remaining after elimination, and without reading the sam file into memory):

	eliminate_strains -i seqs_final_out.fasta.gz -s alignment.sam -f 0.01 -o mismatch.txt -v variants.txt -p -1 reads_trimmed_forward.fastq -2 reads_trimmed_reverse.fastq -n -e 0.005 -g EPI_ISL_402124.fasta -x 35000

After `eliminate_strains` is run, run `EM_C_LLR.R` on the `-o` mismatch matrix (see below).

To impute and build a new database use `sarscov2_imputation`. You must have an account with GISAID and download the Audacity tree and the masked MSA file. You can run our imputation pipeline in the `imputation_scripts` directory (fill in the X with appropriate download dates from GISAID):
	
	./imputation.sh GISAID-hCoV-19-phylogeny-XXXX-XX-XX.zip mmsa_XXXX-XX-XX.tar.xz

`imputation.sh` requires `sarscov2_imputation` in your path and requires the R `ape` package version 5.6-3 or later (<a href="https://github.com/emmanuelparadis/ape">https://github.com/emmanuelparadis/ape</a>).

	sarscov2_imputation [OPTIONS]
	
		-h, --help						usage: -i [Input MSA FASTA] -o [Output file] -t [Tree file]
		-i, --MSA [INFILE, REQUIRED]				multiple sequence alignment in FASTA format
		-m, --final [OUTFILE, REQUIRED if tree imputation]	final out file
		-o, --outfile [OUTFILE, REQUIRED]			imputed output file in FASTA format
		-t, --tree [OUTFILE, REQUIRED if tree imputation]	rooted phylogenetic tree in Newick format
		-c, --common_allele					impute with the most common allele (no tree required)
		-l, --limit [INT]					limit of leaf nodes within a clade [default:10000]
		-v, --variants [OUTFILE, REQUIRED if tree imputation]	file to print variants
		-r, --redundant [OUTFILE, REQUIRED]			file to print redundant strains
	
`sarscov2_imputation` uses <a href="https://github.com/DavidLeeds/hashmap">David Leeds' hashmap</a>.

# Performance
<img src="https://github.com/lpipes/SARS_CoV_2_wastewater_surveillance/blob/main/single_end_300bp.png?raw=true">
<b>Figure 1</b>. Estimated proportions for simulated 300 bp single-end reads with five replicates for when the sample truly contains 1 (A), 3 (B), 5 (C), or 10 (D) strains out of a total of 1,499,078 non-redundant candidate strains in the database. The red dashed lines indicate the true proportion of each strain. 'Other' indicates the sum of estimated proportions for all strains that are not truly represented in the sample.
<img src="https://github.com/lpipes/SARS_CoV_2_wastewater_surveillance/blob/main/bayarea.png?raw=true">
<b>Figure 2</b>. Estimated proportions of the top 25 strains estimated from wastewater shotgun sequencing data from Crits-Cristoph et al. (2021) and their log-likelihood ratios. Strains with an asterisk (*) are identical with other strains. 

## Notes for EM.R
Usage: EM_C_LLR.R [-[-mismatch|i] [<character>]] [-[-error_rate|e] <double>] [-[-filter|f] [<double>]] [-[-llr|l]] [-[-num_show|n] <integer>] [-[-help|h]] [-[-site_LLR|s]] [-[-variant|v] <character>] [-[-read_count|b] <character>] [-[-ref_file|r] <character>] [-[-deletion_file|d] <character>] [-[-core_num|c] <integer>]
	
    -i|--mismatch  [REQUIRED,FILE]                  Name of the file containing the mismatch matrix. (This should be a TXT file.)
    -e|--error_rate                                 Assumed error rate in the EM algorithm (default = 0.005).
    -f|--filter  [REQUIRED,decimal]                 The allele frequency cut-off we used to remove ‘unlikely’ strains.
    -l|--llr                                        Use this to perform the LLR procedure for each strain.
    -n|--num_show [REQUIRED,decimal]                Maximal number of strains to show in the plot (default = 10).
    -h|--help                                       Show help.
    -s|--site_LLR                                   Perform the LLR procedure for each site.
    -v|--variant [REQUIRED,FILE]                    Name of the file containing variant sites.
    -b|--read_count [REQUIRED if test LLR,FILE]     Name of the file containing allele counts from reads.
    -r|--ref_file [REQUIRED if test LLR,FILE]       MSA FASTA of SARS-CoV-2 reference strains.
    -d|--deletion_file [REQUIRED if test LLR,FILE]  Deletion proportion for each site.
    -c|--core_num [REQUIRED,decimal]                Number of cores to use
  
 (Rely on getopt package.)

The estimated proportion of candidate strains will be in the .csv file whose name begins with "em_output_". Also, a barplot will be in the .pdf file whose name starts with "proportion_plot_". Only those strains with an estimated proportion larger than 1% will appear in this plot. If there are **unidentifiable strains**, groups of unidentifiable strains will be printed to a .txt file whose name begins with "Unidentifiable_Strains_" and their group index will appear in the estimated proportions instead of the originial names. If there is no identifiability issue, this .txt file will not be produced.

The time cost of each step will be printed to the standard output stream. It will also print **Errors**  to the standard output stream. Currently, if the input contains less than **20** reads, the code will report there are not enough reads.

I'm still seeking ways to run SQUAREM in a parallel way.

# Data from manuscript
Simulations used in the manuscript can be downloaded at <a href="https://doi.org/10.5281/zenodo.5838942">https://doi.org/10.5281/zenodo.5838942</a>. The imputed MSA can be downloaded at <a href="https://doi.org/10.5281/zenodo.5838946">https://doi.org/10.5281/zenodo.5838946</a>. Identical strains are contained in the headers of the MSA separated by colons.

<a href="https://zenodo.org/badge/latestdoi/308476526"><img src="https://zenodo.org/badge/308476526.svg" alt="DOI"></a>

# References

Crits-Christoph A, Kantor RS, Olm MR, Whitney ON, Al-Shayeb B, Lou YC, Flamholz A, Kennedy LC, Greenwald H, Hinkle A, et al. Genome sequencing of sewage detects regionally prevalent SARS-CoV-2 variants. MBio. 2021;12(1):e02703–20.

# Citations
Pipes, L., Chen, Z., Afanaseva, S., & Nielsen, R. (2022). Estimating the relative proportions of SARS-CoV-2 haplotypes from wastewater samples. Cell reports methods, 100313.
