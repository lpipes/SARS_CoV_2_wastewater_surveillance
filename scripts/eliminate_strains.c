#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "global.h"
#include <string.h>
#include <ctype.h>

int setMSALength(FILE* MSA_file){
	char buffer [FASTA_MAXLINE];
	int length = 0;
	int i=0;
	int iter=0;
	while( fgets(buffer,FASTA_MAXLINE,MSA_file) != NULL ){
		if (buffer[0] != '>'){
			for(i=0; buffer[i]!='\0'; i++){
				length++;
			}
		}else if (iter==0){
			iter++;
		}else{
			break;
		}
	}
	return length;
}

void setNumStrains(FILE* MSA_file, int* strain_info){
	char buffer [FASTA_MAXLINE];
	int i=0;
	int numstrains=0;
	int maxname=0;
	while( fgets(buffer,FASTA_MAXLINE,MSA_file) != NULL ){
		if (buffer[0] == '>'){
			int length = strlen(buffer) -1;
			if (length > maxname){
				maxname=length;
			}
			numstrains++;
		}
	}
	strain_info[0]=numstrains;
	strain_info[1]=maxname;
}

void readInMSA(FILE* MSA_file, char** MSA, int** allele_frequency, char** names, int* reference, int length_of_MSA){
	char buffer [FASTA_MAXLINE];
	int i=0;
	int index=-1;
	int found_ref=0;
	while( fgets(buffer,FASTA_MAXLINE,MSA_file) != NULL ){
		if (buffer[0] == '>'){
			index++;
			for(i=1; buffer[i]!='\n'; i++){
				names[index][i-1]=buffer[i];
			}
			names[index][i-1]='\0';
			if (strcmp(names[index],"NC_045512v2")==0){
				found_ref=1;
			}
			if (strcmp(names[index],"hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30|China")==0){
				found_ref=1;
			}
		}else{
			int size = strlen(buffer);
			for(i=0; i<size; i++){
				if (buffer[i]=='A' || buffer[i]=='a'){
					MSA[index][i]='A';
					allele_frequency[i][0]++;
				}else if (buffer[i]=='G' || buffer[i]=='g'){
					MSA[index][i]='G';
					allele_frequency[i][1]++;
				}else if (buffer[i]=='C' || buffer[i]=='c'){
					MSA[index][i]='C';
					allele_frequency[i][2]++;
				}else if (buffer[i]=='T' || buffer[i]=='t'){
					MSA[index][i]='T';
					allele_frequency[i][3]++;
				}else if (buffer[i]=='-'){
					MSA[index][i]='-';
				}else{
					MSA[index][i]=-1;
				}
			}
			if (found_ref==1){
				int placement = 0;
				for(i=0; i<length_of_MSA; i++){
					if ( MSA[index][i] != '-' ){
						reference[placement] = i;
						placement++;
					}
				}
				found_ref=0;
			}
		}
	}
}

void imputeNucMat(int number_of_strains, int length_of_MSA, int** MSA, int* imputation){
	int i,j;
	for( i=0; i<number_of_strains; i++){
		for(j=0; j<length_of_MSA; j++){
			if ( MSA[i][j] == -1 ){
				MSA[i][j] = imputation[j] + 1;
			}
		}
	}
}

void findMaxAllele(int length_of_MSA, int** allele_frequency, int* imputation){
	int i,j;
	for(i=0; i<length_of_MSA; i++){
		int max=0;
		int max_index=0;
		for(j=0; j<4; j++){
			if (allele_frequency[i][j] > max){
				max = allele_frequency[i][j];
				max_index = j;
			}
		}
		imputation[i]=max_index;
	}
}

int removeIdenticalStrains(int number_of_strains, int length_of_MSA, int** MSA, int* identical, char** names_of_strains, int maxname){
	int i,j, k;
	int mismatch=0;
	int index=0;
	int* identical2 = (int *)malloc(number_of_strains*sizeof(int));
	for(i=0; i<number_of_strains; i++){
		identical2[i]=0;
	}
	for(i=0; i<number_of_strains; i++){
		//printf("working on %d\n",i);
		for(j=i+1; j<number_of_strains; j++){
			mismatch=0;
			for(k=0; k<length_of_MSA; k++){
				if ( MSA[i][k] != MSA[j][k] ){
					mismatch++;
				}
			}
			if (mismatch==0){
				printf("Identical seq found!\n");
				int placement=0;
				for(k=0; k<number_of_strains; k++){
					if (identical[k]==-1){
						placement=k;
						break;
					}
				}
				int found=0;
				for(k=0; k<number_of_strains; k++){
					if (identical[k]==j){
						found=1;	
					}
				}
				if (found==0){
					identical[placement]=j;
				}
			}
		}
	}
	int number_of_identical_strains=0;
	for(i=0; i<number_of_strains; i++){
		if (identical[i]==-1){ break; }
		number_of_identical_strains++;
		//printf("%d\n",identical[i]);
	}
	printf("Number of identical strains: %d\n",number_of_identical_strains);
	for(i=0; i<number_of_strains; i++){
		int ident=0;
		for(j=0; j<number_of_identical_strains; j++){
			if (i==identical[j]){
				ident=1;
			}
		}
		if (ident==0){
			identical2[i] = i;
		}else{
			identical2[i] = -1;
		}
	}
	for(i=0; i<number_of_strains; i++){
		identical[i] = -1;
	}
	int placement = 0;
	for(i=0; i<number_of_strains; i++){
		if ( identical2[i] != -1 ){
			identical[placement] = identical2[i];
			placement++;
		}else{
			memset(names_of_strains[i],'\0',maxname);
		}
	}
	free(identical2);
	return number_of_strains-number_of_identical_strains;
}

int calculateAlleleFreq(FILE* sam, double** allele, int* reference, int length_of_MSA, char** MSA, int number_of_strains, char** names_of_strains, double freq_threshold, int maxname, struct timespec tstart, struct timespec tend, int* identical, int number_of_different_strains, int number_of_variant_sites, int* variant_sites){
	int i,j;
	char buffer [FASTA_MAXLINE];
	char *s;
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	while( fgets(buffer,FASTA_MAXLINE,sam) != NULL ){
		if ( buffer[0] != '@'){
			s = strtok(buffer,"\t");
			for(i=0; i<3; i++){
				s = strtok(NULL,"\t");
			}
			int position=0;
			sscanf(s, "%d", &position);
			position--;
			for(i=0; i<6; i++){
				s = strtok(NULL,"\t");
			}
			char* sequence = s;
			int size = strlen(sequence);
			for(i=0; i<size; i++){
				//printf("pos is %d\n",reference[i+position]);
				if (sequence[i]=='A' || sequence[i]=='a' ){
					allele[reference[i+position]][0]++;
				}else if ( sequence[i]=='G' || sequence[i]=='g'){
					allele[reference[i+position]][1]++;
				}else if ( sequence[i]=='C' || sequence[i]=='c'){
					allele[reference[i+position]][2]++;
				}else if (sequence[i]=='T' || sequence[i]=='t'){
					allele[reference[i+position]][3]++;
				}
			}
			//printf("pos %d: %lf %lf %lf %lf\n",663,allele[663][0],allele[663][1],allele[663][2],allele[663][3]);
			//exit(1);
		}
	}
	for(i=0; i<length_of_MSA; i++){
		double total=0;
		for(j=0; j<4; j++){
			total=total+allele[i][j];
		}
		//printf("pos: %d total is %lf\n",i,total);
		for(j=0; j<4; j++){
			allele[i][j] = allele[i][j]/total;
		}
		//printf("%lf %lf %lf %lf\n",allele[i][0],allele[i][1],allele[i][2],allele[i][3]);
	}
	int length_of_reference = 0;
	for (i=0; i<length_of_MSA; i++){
		if ( reference[i]==-1 ){
			break;
		}else{
			length_of_reference++;
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	printf("Eliminating strains...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	int number_remaining=0;
	int next=0;
	for(i=0; i<length_of_reference; i++){
		if ( i==variant_sites[next] ){
			for(j=0; j<number_of_different_strains; j++){
				if ( MSA[identical[j]][reference[i]] != '-' ){
					int base;
					if ( MSA[identical[j]][reference[i]] == 'A' ){
						base=0;
					}else if ( MSA[identical[j]][reference[i]] == 'G' ){
						base=1;
					}else if ( MSA[identical[j]][reference[i]] == 'C' ){
						base=2;
					}else if ( MSA[identical[j]][reference[i]] == 'T' ){
						base=3;
					}
					if ( allele[reference[i]][base] < freq_threshold ){
						//memset(names_of_strains[identical[j]],'\0',maxname);
						names_of_strains[identical[j]][0] = '\0';
					}
				}
			}
			next++;
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	for(i=0; i<number_of_strains; i++){
		if ( names_of_strains[i][0] != '\0' ){
			printf("Remaining strain: %s\n",names_of_strains[i]);
			number_remaining++;
		}
	}
	printf("Number of strains remaining is %d\n",number_remaining);
	return number_remaining;
}

void writeMismatchMatrix( FILE* outfile, FILE* samfile, char** MSA, int* strains_kept, int length_of_MSA, int number_of_strains, int number_of_strains_remaining, char** names_of_strains, int* reference){
	int i,j;
	char buffer [FASTA_MAXLINE];
	char *s;
	fprintf(outfile,"qName\tblockSizes");
	for(i=0; i<number_of_strains_remaining; i++){
		fprintf(outfile,"\t%s",names_of_strains[strains_kept[i]]);
	}
	fprintf(outfile,"\n");
	while( fgets(buffer,FASTA_MAXLINE,samfile) != NULL ){
		if ( buffer[0] != '@'){
			s = strtok(buffer,"\t");
			fprintf(outfile,"%s",s);
			for(i=0; i<3; i++){
				s = strtok(NULL,"\t");
			}
			int position=0;
			sscanf(s, "%d", &position);
			position--;
			for(i=0; i<2; i++){
				s = strtok(NULL,"\t");
			}
			char alignment_size [MAX_CIGAR];
			sscanf(s, "%s", &(alignment_size));
			int size = strlen(alignment_size);
			char number_to_convert [MAX_CIGAR];
			memset(number_to_convert,'\0',MAX_CIGAR);
			int placement = 0;
			int alignment_size_int = 0;
			for(i=0; i<size; i++){
				if (isalpha(alignment_size[i])){
					if (alignment_size[i]=='M'){
						alignment_size_int += atoi(number_to_convert);
						memset(number_to_convert,'\0',MAX_CIGAR);
						placement=0;
					}
				}else{
					number_to_convert[placement]=alignment_size[i];
					placement++;
				}
			}
			fprintf(outfile,"\t%d",alignment_size_int);
			for(i=0; i<4; i++){
				s = strtok(NULL,"\t");
			}
			char* sequence = s;
			size = strlen(sequence);
			for(i=0; i<number_of_strains_remaining; i++){
				int number_of_mismatches = 0;
				for(j=0; j<size; j++){
					if ( sequence[j] == 'A' || sequence[j] == 'a' ){
						if ( MSA[strains_kept[i]][reference[j+position]] != 'A' && MSA[strains_kept[i]][reference[j+position]] != '-'){
							number_of_mismatches++;
						}
					}else if ( sequence[j] == 'G' || sequence[j] == 'g' ){
						if (MSA[strains_kept[i]][reference[j+position]] != 'G' && MSA[strains_kept[i]][reference[j+position]] != '-'){
							number_of_mismatches++;
						}
					}else if ( sequence[j] == 'C' || sequence[j] == 'c' ){
						if (MSA[strains_kept[i]][reference[j+position]] != 'C' && MSA[strains_kept[i]][reference[j+position]] != '-' ){
							number_of_mismatches++;
						}
					}else if ( sequence[j] == 'T' || sequence[j] == 't' ){
						if (MSA[strains_kept[i]][reference[j+position]] != 'T' && MSA[strains_kept[i]][reference[j+position]] != '-'){
							number_of_mismatches++;
						}
					}
				}
				fprintf(outfile,"\t%d",number_of_mismatches);
			}
			fprintf(outfile,"\n");
		}
	}
}

int* readVariantSitesFile(FILE* variant_sites_file, int* number_of_sites){
	char buffer [20];
	int iter=0;
	int* variant_sites;
	int placement=0;
	while( fgets(buffer,20,variant_sites_file) != NULL){
		if (iter==0){
			number_of_sites[0]=atoi(buffer);
			variant_sites = (int*)malloc(number_of_sites[0]*sizeof(int));
			memset(variant_sites,0,number_of_sites[0]);
			iter++;
		}else{
			variant_sites[placement]=atoi(buffer);
			placement++;
		}
	}
	return variant_sites;
}

int main(int argc, char **argv){
	struct timespec tstart={0,0}, tend={0.0};
	Options opt;
	opt.remove_identical=0;
	parse_options(argc, argv, &opt);
	FILE* MSA_file;
	if (( MSA_file = fopen(opt.fasta,"r")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	int length_of_MSA=0;
	length_of_MSA=setMSALength(MSA_file);
	fclose(MSA_file);
	printf("Length of MSA: %d\n",length_of_MSA);
	if (( MSA_file = fopen(opt.fasta,"r")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	int *strain_info = (int*)malloc(2*sizeof(int));
	setNumStrains(MSA_file,strain_info);
	fclose(MSA_file);
	int number_of_strains = strain_info[0];
	int max_name_length = strain_info[1];
	free(strain_info);
	printf("Number of strains: %d\n",number_of_strains);
	//printf("Maxname length: %d\n",max_name_length);
	char** MSA = (char **)malloc(number_of_strains*sizeof(char *));
	int i,j;
	for(i=0; i<number_of_strains; i++){
		MSA[i] = (char *)malloc(length_of_MSA*sizeof(char));
	}
	char** names_of_strains = (char**)malloc(number_of_strains*sizeof(char *));
	for(i=0; i<number_of_strains; i++){
		names_of_strains[i] = (char *)malloc((max_name_length+1)*sizeof(char));
	}
	int** allele_max = (int **)malloc(length_of_MSA*sizeof(int *));
	for(i=0; i<length_of_MSA; i++){
		allele_max[i] = (int *)malloc(4*sizeof(int));
		for(j=0; j<4; j++){
			allele_max[i][j]=0;
		}
	}
	int* reference = (int *)malloc(length_of_MSA*sizeof(int));
	for(i=0; i<length_of_MSA; i++){
		reference[i] = -1;
	}
	printf("Reading in MSA...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	if (( MSA_file = fopen(opt.fasta,"r")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	readInMSA(MSA_file,MSA,allele_max,names_of_strains,reference,length_of_MSA);
	fclose(MSA_file);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	//int* imputation = (int *)malloc(length_of_MSA*sizeof(int));
	//findMaxAllele(length_of_MSA,allele_max,imputation);
	//for(i=0; i<length_of_MSA; i++){
	//	free(allele_max[i]);
	//}
	//free(allele_max);
	//imputeNucMat(number_of_strains,length_of_MSA,MSA,imputation);
	//free(imputation);
	int* identical = (int *)malloc(number_of_strains*sizeof(int));
	int number_of_identical_strains=number_of_strains;
	if (opt.remove_identical==0){
		for(i=0; i<number_of_strains; i++){
			identical[i]=i;
		}
	}
	//if (opt.remove_identical==1){
	//	printf("Finding identical sequences\n");
	//	clock_gettime(CLOCK_MONOTONIC, &tstart);
	//	for(i=0; i<number_of_strains; i++){
	//		identical[i]=-1;
	//	}
	//	number_of_identical_strains=removeIdenticalStrains(number_of_strains,length_of_MSA,MSA,identical,names_of_strains,max_name_length);
	//	clock_gettime(CLOCK_MONOTONIC, &tend);
	//	printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	//}	
	FILE* variant_sites_file;
	int* variant_sites;
	int* number_of_variant_sites_p = (int *)malloc(sizeof(int));
	if ((variant_sites_file = fopen(opt.variant,"r")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	variant_sites=readVariantSitesFile(variant_sites_file,number_of_variant_sites_p);
	fclose(variant_sites_file);
	int number_of_variant_sites = number_of_variant_sites_p[0];
	free(number_of_variant_sites_p);
	printf("Calculating allele frequencies...\n");
	FILE* sam_file;
	if ((sam_file = fopen(opt.sam,"r")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	double** allele_frequency = (double **)malloc(length_of_MSA*sizeof(double *));
	for(i=0; i<length_of_MSA; i++){
		allele_frequency[i] = (double *)malloc(4*sizeof(double));
		for(j=0; j<4; j++){
			allele_frequency[i][j]=0;
		}
	}
	int number_of_strains_remaining=0;
	number_of_strains_remaining=calculateAlleleFreq(sam_file,allele_frequency,reference,length_of_MSA,MSA,number_of_strains,names_of_strains,opt.freq,max_name_length,tstart,tend,identical,number_of_identical_strains,number_of_variant_sites,variant_sites);
	fclose(sam_file);
	for(i=0; i<length_of_MSA; i++){
		free(allele_frequency[i]);
	}
	free(allele_frequency);
	int* strains_kept = (int *)malloc(number_of_strains_remaining*sizeof(int));
	int placement=0;
	for(i=0; i<number_of_strains; i++){
		if (names_of_strains[i][0] != '\0' ){
			strains_kept[placement]=i;
			placement++;
		}
	}
	printf("Creating mismatch matrix...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	FILE *outfile;
	if (( outfile = fopen(opt.outfile,"w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	if ((sam_file = fopen(opt.sam,"r")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	writeMismatchMatrix(outfile,sam_file,MSA,strains_kept,length_of_MSA,number_of_strains,number_of_strains_remaining,names_of_strains,reference);
	fclose(sam_file);
	fclose(outfile);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fseconds\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	free(strains_kept);
	free(reference);
	for(i=0; i<number_of_strains; i++){
		free(names_of_strains[i]);
		free(MSA[i]);
	}
	free(names_of_strains);
	free(MSA);
}
