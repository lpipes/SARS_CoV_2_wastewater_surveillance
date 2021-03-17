# SARS_CoV_2_wastewater_surveillance

## Notes for EM.R
To run EM algorithm, use command: Rscript EM.R [mismatch matrix.txt] [error rate (optional, default = 0.005)]

The estimated proportion of candidate strains will be in the .csv file whose name begins with "em_output_". Also, a barplot will be in the .pdf file whose name starts with "proportion_plot_". Only those strains with an estimated proportion larger than 1% will appear in this plot. If there are **unidentifiable strains**, groups of unidentifiable strains will be printed to a .txt file whose name begins with "Unidentifiable_Strains_" and their group index will appear in the estimated proportions instead of the originial names. If there is no unidentifiability issue, this .txt file will not be produced.

The time cost of each step, the most likely strain when assuming only one strain exists, and the likelihood ratio will be printed to the standard output stream. It will also print **Errors**  to the standard output stream. Currently, if the input contains less than **20** reads, the code will report there are not enough reads.

Currently, the code will also generate a file called Rplot.pdf, which is useless. (I'm trying to fix this.)

I'm still seeking ways to run SQUAREM in a parallel way.
