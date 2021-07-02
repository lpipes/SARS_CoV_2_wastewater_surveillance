#ifndef _GLOBAL_
#define _GLOBAL_

#define FASTA_MAXLINE 40000

typedef struct Options{
	char fasta[1000];
	int remove_identical;
	int print_variant;
	char variant[1000];
	char out_MSA[1000];
	char reference[1000];
}Options;

typedef struct blob{
	char* key;
	size_t data_len;
	char* data;
}blob;

#endif
