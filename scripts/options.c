#include "options.h"

static struct option long_options[]=
{
	{"help", no_argument, 0, 'h'},
	{"infile", required_argument, 0, 'i'},
	{"samfile", required_argument, 0, 's'},
	{"freq", required_argument, 0, 'f'},
	{"outfile", required_argument, 0, 'o'},
	{"variant_sites", required_argument, 0, 'v'},
	{"remove_identical", no_argument, 0, 'r'},
	{0,0,0,0}
};

char usage[] = "\neliminate_strains [OPTIONS]\n\
	\n\
	-h, --help				usage: -i [Input MSA FASTA]\n\
	-i, --infile [REQUIRED]			MSA FASTA\n\
	-s, --samfile [REQUIRED]		SAM\n\
	-f, --freq [REQUIRED]			allele frequency to filter\n\
	-o, --outfile [REQUIRED]		outfile\n\
	-v, --variant_sites [REQUIRED]		variant sites file\n\
	-r, --remove_identical			remove identical sequences\n\
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
		c=getopt_long(argc,argv,"hri:s:f:o:v:",long_options, &option_index);
		if (c==-1) break;
		switch(c){
			case 'h':
				print_help_statement();
				exit(0);
				break;
			case 'r':
				opt->remove_identical=1;
				break;
			case 'i':
				success = sscanf(optarg, "%s", opt->fasta);
				if (!success)
					fprintf(stderr, "Invalid fasta file\n");
				break;
			case 's':
				success = sscanf(optarg, "%s", opt->sam);
				if (!success)
					fprintf(stderr, "Invalid sam file\n");
				break;
			case 'f':
				success = sscanf(optarg, "%lf", &(opt->freq));
				if (!success)
					fprintf(stderr, "Invalid freq\n");
				break;
			case 'o':
				success = sscanf(optarg, "%s", opt->outfile);
				if (!success)
					fprintf(stderr, "Invalid out file\n");
				break;
			case 'v':
				success = sscanf(optarg, "%s", opt->variant);
				if (!success)
					fprintf(stderr, "Invalid file\n");
				break;
		}
	}
}
