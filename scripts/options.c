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
	{"paired", no_argument, 0, 'p'},
	{"single_end", required_argument, 0, '0'},
	{"forward_read", required_argument, 0, '1'},
	{"reverse_read", required_argument, 0, '2'},
	{"bowtie-db", required_argument, 0, 'd'},
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
	-p, --paired				paired-reads\n\
	-0, --single_end_file [REQUIRED]	single-end fasta\n\
	-1, --forward_file [REQUIRED]		forward fasta\n\
	-2, --reverse_file [REQUIRED]		reverse fasta\n\
	-d, --bowtie-db [REQUIRED]		bowtie db\n\
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
		c=getopt_long(argc,argv,"hpri:s:f:o:v:0:1:2:d:",long_options, &option_index);
		if (c==-1) break;
		switch(c){
			case 'h':
				print_help_statement();
				exit(0);
				break;
			case 'd':
				success = sscanf(optarg, "%s", opt->bowtie_reference_db);
				if (!success)
					fprintf(stderr, "Invalid bowtie db\n");
				break;
			case 'r':
				opt->remove_identical=1;
				break;
			case 'p':
				opt->paired=1;
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
			case '0':
				success = sscanf(optarg, "%s", opt->single_end_file);
				if (!success)
					fprintf(stderr, "Invalid fasta file\n");
				break;
			case '1':
				success = sscanf(optarg, "%s", opt->forward_end_file);
				if (!success)
					fprintf(stderr, "Invalid fasta file\n");
				break;
			case '2':
				success = sscanf(optarg, "%s", opt->reverse_end_file);
				if (!success)
					fprintf(stderr, "Invalid fasta file\n");
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
