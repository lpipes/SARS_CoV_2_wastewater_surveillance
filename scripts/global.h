#ifndef _GLOBAL_
#define _GLOBAL_

#define FASTA_MAXLINE 40000
#define	MAX_CIGAR 1000

typedef struct Options{
	char fasta[1000];
	char sam[1000];
	double freq;
	char outfile[1000];
	int remove_identical;
}Options;

#endif
