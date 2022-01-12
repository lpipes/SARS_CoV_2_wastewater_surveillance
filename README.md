# SARS_CoV_2_wastewater_surveillance
<img src="https://github.com/lpipes/SARS_CoV_2_wastewater_surveillance/blob/main/single_end_300bp.png?raw=true">
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
