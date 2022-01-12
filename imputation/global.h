#ifndef _GLOBAL_
#define _GLOBAL_

#define MAX_NODENAME 100
#define FASTA_MAXLINE 50000
#define NUMCAT 1
#define MINBL 0.000001
#define MAXBL 2.0
#define MAXNUMBEROFINDINSPECIES 500 /*maximum number of individuals belonging to a species*/
#define MAXFILENAME 500
#define STATESPACE 20 /*number of categories in approximation of gamma distribution for Ne. Must be at least 4 because some of the memory is used for the nucleotide model*/
#define MAX_NODE_LIST 100000
#define MAX_POLYTOMIES 100
#define EPSILON 0.01
typedef struct node{
	int up[2];
	int down;
	double bl;
	char *name;
	int nd;
	int depth;
	int pos;
	double **likenc;
	double **posteriornc;
}node;

typedef struct Options{
	char outfile[MAXFILENAME];
	char msa[MAXFILENAME];
	char tree[MAXFILENAME];
	int common;
	int limit;
}Options;

#endif
