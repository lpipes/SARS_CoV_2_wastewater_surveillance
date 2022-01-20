# Method for estimating relative proportions of SARS-CoV-2 strains from wastewater samples
We present a method for estimating the relative proportions of SARS-CoV-2 strains from wastewater samples. The method uses an initial step to remove unlikely strains, imputation of missing nucleotides using the global SARS-CoV-2 phylogeny, and an Expectation-Maximization (EM) algorithm for obtaining maximum likelihood estimates of the proportions of different strains in a sample.

View the preprint here: <a href="https://www.medrxiv.org/content/10.1101/2022.01.13.22269236v2">https://www.medrxiv.org/content/10.1101/2022.01.13.22269236v2</a>

The method has 4 different components: imputation, redundancy removal, filtering unlikely strains, and running the EM algorithm. To estimate proportions of SARS-CoV-2 strains on a pre-built database of 1,499,078 non-redundant strains <a href="https://doi.org/10.5281/zenodo.5838946">download the imputed multiple sequence alignment<a/> and run the `eliminate_strains` program.

# Installation
To install
	
	git clone https://github.com/lpipes/SARS_CoV_2_wastewater_surveillance.git
	cd SARS_CoV_2_wastewater_surveillance/eliminate_strains
	make
	cd ../imputation
	make

# Usage
`eliminate_strains` filters unlikely SARS-CoV-2 genomes, prints a mismatch matrix, and also runs the `EM_C_LLR.R` program.

	eliminate_strains [OPTIONS]
	
	-h, --help				
	-i, --infile [REQUIRED]			MSA FASTA of SARS-CoV-2 reference strains
	-s, --samfile [REQUIRED]		output sam file to print alignments
	-f, --freq [REQUIRED]			allele frequency to filter unlikely strains
	-o, --outfile [REQUIRED]		output file to print mismatch matrix for EM algorithm
	-v, --variant_sites [REQUIRED]		list of variant sites
	-d, --bowtie-db [REQUIRED]		path to Wuhan-Hu-1 bowtie2 database
	-p, --paired				using paired-reads
	-0, --single_end_file			single-end reads
	-1, --forward_file			if using paired-reads, the forward reads file
	-2, --reverse_file			if using paired-reads, the reverse reads file
	-e, --EM-error				error rate for EM algorithm
	-c, --coverage				number of reads needed to calculate allele freq [default: 1]
	
To use the bowtie2 database compatible with the pre-built database use the following bowtie2 database:
	
	eliminate_strains/EPI_ISL_402124.fasta

And the following variant sites file:
	
	eliminate_strains/global_variants.txt

To impute and build a new database use `sarscov2_imputation`.


	sarscov2_imputation [OPTIONS]
	
		-h, --help						usage: -i [Input MSA FASTA] -o [Output file] -t [Tree file]
		-i, --MSA [INFILE, REQUIRED]				multiple sequence alignment in FASTA format
		-m, --final [OUTFILE, REQUIRED if tree imputation]	final out file
		-o, --outfile [OUTFILE, REQUIRED]			imputed output file in FASTA format
		-t, --tree [OUTFILE, REQUIRED if tree imputation]	rooted phylogenetic tree in Newick format
		-c, --common_allele					impute with the most common allele (no tree required)
		-l, --limit [INT]					limit of leaf nodes within a clade [default:10000]
		-v, --variants [OUTFILE, REQUIRED if tree imputation]	file to print variants
	
`sarscov2_imputation` uses <a href="https://github.com/DavidLeeds/hashmap">David Leeds' hashmap</a>.

# Performance
<img src="https://github.com/lpipes/SARS_CoV_2_wastewater_surveillance/blob/main/single_end_300bp.png?raw=true">
<b>Figure 1</b>. Estimated proportions for simulated 300 bp single-end reads with five replicates for when the sample truly contains 1 (A), 3 (B), 5 (C), or 10 (D) strains out of a total of 1,499,078 non-redundant candidate strains in the database. The red dashed lines indicate the true proportion of each strain. 'Other' indicates the sum of estimated proportions for all strains that are not truly represented in the sample.
<img src="https://github.com/lpipes/SARS_CoV_2_wastewater_surveillance/blob/main/bayarea.png?raw=true">
<b>Figure 2</b>. Estimated proportions of the top 25 strains estimated from wastewater shotgun sequencing data from Crits-Cristoph et al. (2021) and their log-likelihood ratios. Strains with an asterisk (*) are identical with other strains. 

## Notes for EM.R
Usage: EM_C_LLR.R [-[-mismatch|i] [<character>]] [-[-error_rate|e] <double>] [-[-filter|f] [<double>]] [-[-llr|l]] [-[-num_show|n] <integer>] [-[-help|h]]
    
    -i|--mismatch      Name of file containing the mismatch matrix. (This should be a TXT file.)
    
    -e|--error_rate    Assumed error rate in the EM algorithm (default = 0.005).
    
    -f|--filter        The allele frequency cut-off we used to remove ‘unlikely’ strains.
    
    -l|--llr           Use this to perform the LLR procedure.
    
    -n|--num_show      Maximal number of strains to show in the plot (default = 10).
    
    -h|--help          Show help.
    
 (Rely on getopt package.)

The estimated proportion of candidate strains will be in the .csv file whose name begins with "em_output_". Also, a barplot will be in the .pdf file whose name starts with "proportion_plot_". Only those strains with an estimated proportion larger than 1% will appear in this plot. If there are **unidentifiable strains**, groups of unidentifiable strains will be printed to a .txt file whose name begins with "Unidentifiable_Strains_" and their group index will appear in the estimated proportions instead of the originial names. If there is no identifiability issue, this .txt file will not be produced.

The time cost of each step will be printed to the standard output stream. It will also print **Errors**  to the standard output stream. Currently, if the input contains less than **20** reads, the code will report there are not enough reads.

I'm still seeking ways to run SQUAREM in a parallel way.

# Data from manuscript
Simulations used in the manuscript can be downloaded at <a href="https://doi.org/10.5281/zenodo.5838942">https://doi.org/10.5281/zenodo.5838942</a>. The imputed MSA can be downloaded at <a href="https://doi.org/10.5281/zenodo.5838946">https://doi.org/10.5281/zenodo.5838946</a>. Identical strains are contained in the headers of the MSA separated by colons.

# References

Crits-Christoph A, Kantor RS, Olm MR, Whitney ON, Al-Shayeb B, Lou YC, Flamholz A, Kennedy LC, Greenwald H, Hinkle A, et al. Genome sequencing of sewage detects regionally prevalent SARS-CoV-2 variants. MBio. 2021;12(1):e02703–20.

# Citations
Pipes L, Chen Z, Afanaseva S, Nielsen R (2022) Estimating the relative proportions of SARS-CoV-2 strains from wastewater samples. <a href="https://www.medrxiv.org/content/10.1101/2022.01.13.22269236v2">medRxiv</a>.
