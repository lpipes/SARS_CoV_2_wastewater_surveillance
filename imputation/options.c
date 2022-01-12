#include "options.h"

static struct option long_options[]=
{
	{"help", no_argument, 0, 'h'},
	{"MSA", required_argument, 0, 'i'},
	{"tree", optional_argument, 0, 't'},
	{"common_allele", optional_argument, 0, 'c'},
	{"outfile", optional_argument, 0, 'o'},
	{"limit", optional_argument, 0, 'l'},
	{0,0,0,0}
};

char usage[] = "\nsarscov2_imputation [OPTIONS]\n\
	\n\
	-h, --help			usage: -i [Input MSA FASTA] -o [Output file] -t [Tree file]\n\
	-i, --MSA [REQUIRED]		multiple sequence alignment in FASTA format\n\
	-o, --outfile [REQUIRED]	imputed output file in FASTA format\n\
	-t, --tree			phylogenetic tree in Newick format\n\
	-c, --common_allele		impute with the most common allele (no tree required)\n\
	-l, --limit			limit of leaf nodes within a clade [default:10000]\n\
	\n";

void print_help_statement(){
	printf("%s", &usage[0]);
	return;
}

void parse_options(int argc, char **argv, Options *opt){
	int option_index, success;
	char c;
	if (argc==1){
		print_help_statement();
		exit(0);
	}
	while(1){
		c=getopt_long(argc,argv,"hci:t:o:l:",long_options, &option_index);
		if (c==-1) break;
		switch(c){
			case 'h':
				print_help_statement();
				exit(0);
				break;
			case 'i':
				success = sscanf(optarg, "%s", opt->msa);
				if (!success)
					fprintf(stderr, "Invalid fasta file\n");
				break;
			case 't':
				success = sscanf(optarg, "%s", opt->tree);
				if (!success)
					fprintf(stderr, "Invalid tree file\n");
				break;
			case 'o':
				success = sscanf(optarg, "%s", opt->outfile);
				if (!success)
					fprintf(stderr, "Invalid output file\n");
				break;
			case 'l':
				success = sscanf(optarg, "%d", &(opt->limit));
				if (!success)
					fprintf(stderr, "Could not read limit number\n");
				break;
			case 'c':
				opt->common=1;
				break;
		}
	}
}
