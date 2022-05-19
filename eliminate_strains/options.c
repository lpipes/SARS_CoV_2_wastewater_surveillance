#include "options.h"

static struct option long_options[]=
{
	{"help", no_argument, 0, 'h'},
	{"infile", required_argument, 0, 'i'},
	{"samfile", required_argument, 0, 's'},
	{"freq", required_argument, 0, 'f'},
	{"outfile", required_argument, 0, 'o'},
	{"variant_sites", required_argument, 0, 'v'},
	{"paired", no_argument, 0, 'p'},
	{"single_end", required_argument, 0, '0'},
	{"forward_read", required_argument, 0, '1'},
	{"reverse_read", required_argument, 0, '2'},
	{"bowtie-db", required_argument, 0, 'd'},
	{"EM-error",required_argument,0, 'e'},
	{"coverage",required_argument,0,'c'},
	{"fasta",no_argument,0,'a'},
	{"llr",no_argument,0,'l'},
	{"min",required_argument,0,'m'},
	{"max",required_argument,0,'x'},
	{"print-allele-counts",required_argument,0,'b'},
	{"cores",required_argument,0,'t'},
	{"msa-reference",required_argument,0,'g'},
	{0,0,0,0}
};

char usage[] = "\neliminate_strains [OPTIONS]\n\
	\n\
	-h, --help				\n\
	-i, --infile [REQUIRED,FILE]		MSA FASTA of SARS-CoV-2 reference strains\n\
	-s, --samfile [REQUIRED,FILE]		output sam file to print alignments\n\
	-f, --freq [REQUIRED,decimal]		allele frequency to filter unlikely strains [default: 0.01]\n\
	-o, --outfile [REQUIRED,FILE]		output file to print mismatch matrix for EM algorithm\n\
	-v, --variant_sites [REQUIRED,FILE]	list of variant sites\n\
	-d, --bowtie-db [REQUIRED,FILE]		path to Wuhan-Hu-1 bowtie2 database\n\
	-p, --paired				using paired-reads\n\
	-0, --single_end_file [FILE]		single-end reads\n\
	-1, --forward_file [FILE]		if using paired-reads, the forward reads file\n\
	-2, --reverse_file [FILE]		if using paired-reads, the reverse reads file\n\
	-e, --EM-error [decimal]		error rate for EM algorithm\n\
	-c, --coverage [integer]		number of reads needed to calculate allele freq [default: 50]\n\
	-a, --fasta				reads are in FASTA format [default: FASTQ]\n\
	-l, --llr				Perform the LLR procedure\n\
	-m, --min [decimal]			Minimum strains remaining to invoke iterative procedure [default: 100]\n\
	-x, --max [decimal]			Maximum strains remaining for EM algorithm [default: 10000]\n\
	-b, --print-allele-counts [FILE]	Print allele counts to file\n\
	-t, --cores [decimal]			Number of cores [default: 1]\n\
	-g, --msa-reference [FILE]		MSA reference index\n\
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
		c=getopt_long(argc,argv,"hplai:s:f:o:v:0:1:2:d:e:t:c:m:x:b:g:",long_options, &option_index);
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
			case 'b':
				success = sscanf(optarg, "%s", opt->print_counts);
				if (!success)
					fprintf(stderr, "Invalid counts file\n");
				break;
			case 'p':
				opt->paired=1;
				break;
			case 'a':
				opt->fasta_format=1;
				break;
			case 'l':
				opt->llr=1;
				break;
			case 'g':
				success = sscanf(optarg, "%s", opt->MSA_reference);
				if (!success)
					fprintf(stderr, "Invalid fasta file\n");
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
			case 't':
				success = sscanf(optarg, "%d", &(opt->number_of_cores));
				if( !success)
					fprintf(stderr, "Invalid number of cores\n");
				break;
			case 'c':
				success = sscanf(optarg, "%d", &(opt->coverage));
				if (!success)
					fprintf(stderr, "Invalid freq\n");
				break;
			case 'm': 
				success = sscanf(optarg, "%d", &(opt->min_strains));
				if (!success)
					fprintf(stderr, "Invalid min strains\n");
				break;
			case 'x': 
				success = sscanf(optarg, "%d", &(opt->max_strains));
				if (!success)
					fprintf(stderr, "Invalid min strains\n");
				break;
			case 'e':
				success = sscanf(optarg, "%lf", &(opt->error));
				if (!success)
					fprintf(stderr, "Invalid error rate\n");
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
