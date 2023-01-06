#ifndef _GLOBAL_
#define _GLOBAL_

#define FASTA_MAXLINE 30000
#define	MAX_CIGAR 1000
#define MAX_READ_LENGTH 1000

typedef struct Options{
	char fasta[1000];
	char sam[1000];
	double freq;
	char outfile[1000];
	int remove_identical;
	char variant[1000];
	int paired;
	int fasta_format;
	char single_end_file[1000];
	char forward_end_file[1000];
	char reverse_end_file[1000];
	double error;
	int clean_reads;
	int coverage;
	int llr;
	int min_strains;
	int max_strains;
	char print_counts[1000];
	char MSA_reference[1000];
	int number_of_cores;
	int no_read_bam;
	char print_deletions[1000];
	double deletion_threshold;
}Options;

typedef struct resultsStruct{
	char **mismatch;
}resultsStruct;

typedef struct thread_struct{
	int start;
	int end;
	int thread_number;
	int max_sam_length;
	int length_of_MSA;
	int number_of_strains;
	int number_of_strains_remaining;
	resultsStruct *str; 
}thread_struct;
extern char** resize_MSA;
extern char** resize_names_of_strains;
extern int* reference_index;
extern char** sam_results;
#endif
