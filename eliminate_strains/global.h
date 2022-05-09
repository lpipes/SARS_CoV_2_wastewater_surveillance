#ifndef _GLOBAL_
#define _GLOBAL_

#define FASTA_MAXLINE 50000
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
	char bowtie_reference_db[1000];
	double error;
	int coverage;
	int llr;
	char print_counts[1000];
}Options;

typedef struct blob{
	char* key;
	size_t data_len;
	char* data;
}blob;

#endif
