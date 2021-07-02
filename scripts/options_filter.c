#include "options_filter.h"

static struct option long_options[]=
{
	{"help", no_argument, 0, 'h'},
	{"infile", required_argument, 0, 'i'},
	{"variant_sites", required_argument, 0, 'v'},
	{"remove_identical", no_argument, 0, 'r'},
	{"out_MSA", required_argument, 0, 'o'},
	{"print_variant_sites", no_argument, 0, 'p'},
	{0,0,0,0}
};

char usage[] = "\nremove_redundant [OPTIONS]\n\
	\n\
	-h, --help				usage: -i [Input MSA FASTA]\n\
	-i, --infile [REQUIRED]			MSA FASTA\n\
	-v, --variant_sites [REQUIRED]		variant sites file to print\n\
	-r, --remove_identical			remove identical sequences\n\
	-o, --out_MSA [REQUIRED]		out MSA\n\
	-p, --print_variant_sites		print variant sites\n\
	-f, --reference_pos [REQUIRED]		print reference positions\n\
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
		c=getopt_long(argc,argv,"hpri:v:o:f:",long_options, &option_index);
		if (c==-1) break;
		switch(c){
			case 'h':
				print_help_statement();
				exit(0);
				break;
			case 'r':
				opt->remove_identical=1;
				break;
			case 'p':
				opt->print_variant=1;
				break;
			case 'i':
				success = sscanf(optarg, "%s", opt->fasta);
				if (!success)
					fprintf(stderr, "Invalid fasta file\n");
				break;
			case 'o':
				success = sscanf(optarg, "%s", opt->out_MSA);
				if (!success)
					fprintf(stderr, "Invalid fasta file\n");
				break;
			case 'f':
				success = sscanf(optarg, "%s", opt->reference);
				if (!success)
					fprintf(stderr, "Invalid file\n");
				break;
			case 'v':
				success = sscanf(optarg, "%s", opt->variant);
				if (!success)
					fprintf(stderr, "Invalid file\n");
				break;
		}
	}
}
