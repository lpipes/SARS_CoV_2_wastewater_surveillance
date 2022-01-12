# SARS_CoV_2_wastewater_surveillance
We present a method for estimating the relative proportions of SARS-CoV-2 strains from wastewater samples. The method uses an initial step to remove unlikely strains, imputation of missing nucleotides using the global SARS-CoV-2 phylogeny, and an Expectation-Maximization (EM) algorithm for obtaining maximum likelihood estimates of the proportions of different strains in a sample.
# Performance
Estimated proportions for simulated 300 bp single-end reads with five replicates for when the sample truly contains 1 (A), 3 (B), 5 (C), or 10 (D) strains out of a total of 1,499,078 non-redundant candidate strains in the database. The red dashed lines indicate the true proportion of each strain. 'Other' indicates the sum of estimated proportions for all strains that are not truly represented in the sample.
<img src="https://github.com/lpipes/SARS_CoV_2_wastewater_surveillance/blob/main/single_end_300bp.png?raw=true">
Estimated proportions of the top 25 strains estimated from wastewater shotgun sequencing data from Crits-Cristoph et al. (2021) and their log-likelihood ratios. Strains with an asterisk (*) are identical with other strains. EPI\_ISL\_682010* is identical to EPI\_ISL\_682025, EPI\_ISL\_1373628, EPI\_ISL\_1373632, and EPI\_ISL\_1373659. EPI\_ISL\_451226* is identical to EPI\_ISL\_451227 and EPI\_ISL\_455983. EPI\_ISL\_625508* is identical to EPI\_ISL\_625520, EPI\_ISL\_672318, EPI\_ISL\_672449, EPI\_ISL\_739003, EPI\_ISL\_739029, EPI\_ISL\_739135, EPI\_ISL\_739161, EPI\_ISL\_739207, and EPI\_ISL\_739286. EPI\_ISL\_1859609* is identical to EPI\_ISL\_1859762. EPI\_ISL\_510925* is identical to EPI\_ISL\_510926. EPI\_ISL\_426109* is identical to  EPI\_ISL\_486012, EPI\_ISL\_570168, EPI\_ISL\_570172, EPI\_ISL\_576500, and EPI\_ISL\_576501. EPI\_ISL\_1074397* is identical to EPI\_ISL\_2190584. EPI\_ISL\_517805* is identical to EPI\_ISL\_527398 and EPI\_ISL_137362.
<img src="https://github.com/lpipes/SARS_CoV_2_wastewater_surveillance/blob/main/bayarea.png?raw=true">
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
