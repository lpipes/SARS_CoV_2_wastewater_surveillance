#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "global_filter.h"
#include "hashmap.h"
#include <string.h>

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

int readInMSA(FILE* MSA_file, int** MSA, int** allele_frequency, char** names, int* reference, int length_of_MSA){
	char buffer [FASTA_MAXLINE];
	int i=0;
	int index=-1;
	int found_ref=0;
	int ref_index=0;
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
		}else{
			int size = strlen(buffer);
			for(i=0; i<size; i++){
				if (buffer[i]=='A' || buffer[i]=='a'){
					MSA[index][i]=1;
					allele_frequency[i][0]++;
				}else if (buffer[i]=='G' || buffer[i]=='g'){
					MSA[index][i]=2;
					allele_frequency[i][1]++;
				}else if (buffer[i]=='C' || buffer[i]=='c'){
					MSA[index][i]=3;
					allele_frequency[i][2]++;
				}else if (buffer[i]=='T' || buffer[i]=='t'){
					MSA[index][i]=4;
					allele_frequency[i][3]++;
				}else if (buffer[i]=='-'){
					MSA[index][i]=0;
				}else{
					MSA[index][i]=-1;
				}
			}
			if (found_ref==1){
				int placement = 0;
				for(i=0; i<length_of_MSA; i++){
					if ( MSA[index][i] != 0 ){
						reference[placement] = i;
						placement++;
					}
				}
				found_ref=0;
				ref_index=index;
			}
		}
	}
	return ref_index;
}

void imputeNucMat(int number_of_strains, int length_of_MSA, int** MSA, int* imputation){
	int i,j;
	for( i=0; i<number_of_strains; i++){
		for(j=0; j<length_of_MSA; j++){
			if ( MSA[i][j] == -1){
				MSA[i][j] = imputation[j] + 1;
			}
			if (MSA[i][j] == 0 ){
				MSA[i][j] = imputation[j] +1;
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

int removeIdenticalStrains(int number_of_strains, int length_of_MSA, int** MSA, int* identical, char** names_of_strains, int maxname, Options opt, int ref_index, int start_of_ref, int end_of_polyA, int* strains_remaining){
	int i,j, k;
	int mismatch=0;
	int index=0;
	int* identical2 = (int *)malloc(number_of_strains*sizeof(int));
	HASHMAP(char, struct blob) hash;
	//hashmap_init(&hash, hashmap_hash_string, hashmap_compare_string, 2*number_of_strains);
	hashmap_init(&hash, hashmap_hash_string, strcmp);
	char** seq = (char **)malloc(number_of_strains*sizeof(char*));
	for(i=0; i<number_of_strains; i++){
		seq[i] = (char *)malloc((length_of_MSA+1)*sizeof(char));
		memset(seq[i],'\0',length_of_MSA+1);
	}
	struct blob *b;
	b = malloc(sizeof(blob));
	//for(i=length_of_MSA-1; i>=0; i--){
	//	if ( MSA[ref_index][i] != 1 && MSA[ref_index][i] != 0){
	//		end_of_polyA = i;
	//		break;
	//	}
	//}
	//int start_of_ref = 0;
	//for(i=0; i<length_of_MSA; i++){
	//	if ( MSA[ref_index][i] != 0 ){
	//		start_of_ref = i;
	//		break;
	//	}
	//}
	for(i=0; i<number_of_strains; i++){
		for(j=start_of_ref; j<end_of_polyA+1; j++){
			if ( MSA[i][j] == 0 ){
				seq[i][j-start_of_ref] = '-';
			}else if (MSA[i][j] == 1 ){
				seq[i][j-start_of_ref] = 'A';
			}else if (MSA[i][j] == 2 ){
				seq[i][j-start_of_ref] = 'G';
			}else if (MSA[i][j] == 3 ){
				seq[i][j-start_of_ref] = 'C';
			}else if (MSA[i][j] == 4 ){
				seq[i][j-start_of_ref] = 'T';
			}
		}
		seq[i][end_of_polyA+1-start_of_ref]='\0';
		//b->key = seq[i];
		//b->data = names_of_strains[i];
		hashmap_put(&hash, seq[i], names_of_strains[i]);
	}
	hashmap_put(&hash, seq[ref_index], names_of_strains[ref_index]);
	const char *key;
	const char *value;
	FILE* remove_ident_file;
	int remaining_strains = 0;
	if (( remove_ident_file = fopen(opt.out_MSA,"w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	hashmap_foreach(key, value, &hash){
		fprintf(remove_ident_file,">%s\n",value);
		fprintf(remove_ident_file,"%s\n",key);
		remaining_strains++;
	}
	fclose(remove_ident_file);
	printf("remaining strains: %d\n",remaining_strains);
	//strains_remaining = (int *)malloc(sizeof(int)*remaining_strains);
	j=0;
	hashmap_foreach(key,value, &hash){
		for(i=0; i<number_of_strains; i++){
			if (strcmp(value,names_of_strains[i])==0){
				strains_remaining[j]=i;
				j++;
			}
		}
	}
	//exit(1);
	//for(i=0; i<number_of_strains; i++){
	//	identical2[i]=0;
	//}
	/*for(i=0; i<number_of_strains; i++){
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
	}*/
	/*int number_of_identical_strains=0;
	for(i=0; i<number_of_strains; i++){
		if (identical[i]==-1){ break; }
		number_of_identical_strains++;
		//printf("%d\n",identical[i]);
	}
	printf("Number of identical strains: %d\n",number_of_strains-remaining_strains);
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
	free(identical2);*/
	//return number_of_strains-number_of_identical_strains;
	return remaining_strains;
}

int calculateAlleleFreq(FILE* sam, double** allele, int* reference, int length_of_MSA, int** MSA, int number_of_strains, char** names_of_strains, double freq_threshold, int maxname, struct timespec tstart, struct timespec tend, int* identical, int number_of_different_strains){
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
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	printf("Eliminating strains...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	int number_remaining=0;
	for(i=0; i<length_of_MSA; i++){
		for(j=0; j<number_of_different_strains; j++){
			if ( MSA[identical[j]][reference[i]] > 0 ){
				if ( allele[reference[i]][MSA[identical[j]][reference[i]]-1] < freq_threshold ){
					memset(names_of_strains[identical[j]],'\0',maxname);
				}
			}
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

void writeMismatchMatrix( FILE* outfile, FILE* samfile, int** MSA, int* strains_kept, int length_of_MSA, int number_of_strains, int number_of_strains_remaining, char** names_of_strains, int* reference){
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
			for(i=0; i<6; i++){
				s = strtok(NULL,"\t");
			}
			char* sequence = s;
			int size = strlen(sequence);
			for(i=0; i<number_of_strains_remaining; i++){
				int number_of_mismatches = 0;
				for(j=0; j<size; j++){
					if ( sequence[i] == 'A' || sequence[i] == 'a' ){
						if ( MSA[strains_kept[i]][reference[i+position]] != 1 && MSA[strains_kept[i]][reference[i+position]] != 0){
							number_of_mismatches++;
						}
					}else if ( sequence[i] == 'G' || sequence[i] == 'g' ){
						if (MSA[strains_kept[i]][reference[i+position]] != 2 && MSA[strains_kept[i]][reference[i+position]] != 0){
							number_of_mismatches++;
						}
					}else if ( sequence[i] == 'C' || sequence[i] == 'c' ){
						if (MSA[strains_kept[i]][reference[i+position]] != 3 && MSA[strains_kept[i]][reference[i+position]] != 0 ){
							number_of_mismatches++;
						}
					}else if ( sequence[i] == 'T' || sequence[i] == 't' ){
						if (MSA[strains_kept[i]][reference[i+position]] != 4 && MSA[strains_kept[i]][reference[i+position]] != 0){
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

void findVariantSites(int number_of_strains, int length_of_MSA, int** MSA, int* invariant_sites, Options opt, int ref_index, int start_of_ref, int end_of_polyA, int* strains_remaining, int* reference){
	int i,j;
	int number_of_sites=0;
	//int start_of_ref = 0;
	//int end_of_polyA=0;
	//for(i=length_of_MSA-1; i>=0; i--){
	//	if ( MSA[ref_index][i] != 1 && MSA[ref_index][i] != 0){
	//		end_of_polyA = i;
	//		break;
	//	}
	//}
	//for(i=0; i<length_of_MSA; i++){
	//	if ( MSA[ref_index][i] != 0 ){
	//		start_of_ref = i;
	//		break;
	//	}
	//}
	for(i=start_of_ref; i<end_of_polyA+1; i++){
		int base=MSA[0][i];
		int invariant=0;
		for(j=0; j<number_of_strains; j++){
			if ( base != MSA[strains_remaining[j]][i] ){
				invariant++;
			}
		}
		if (invariant==0){
			invariant_sites[i]=i;
		}else{
			number_of_sites++;
		}
	}
	FILE *invariant_file;
	if (( invariant_file = fopen(opt.variant,"w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	fprintf(invariant_file,"%d\n",number_of_sites);
	for(i=start_of_ref; i<end_of_polyA+1; i++){
		if (invariant_sites[i]==-1){
			fprintf(invariant_file,"%d\n",i-start_of_ref);
		}
	}
	fclose(invariant_file);
	/*FILE *reference_file;
	if (( reference_file = fopen(opt.reference,"w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	for(i=0; i<end_of_polyA+1; i++){
		if (reference[i]==-1){ break;}
		if (reference[i] > end_of_polyA){ break;}
		fprintf(reference_file,"%d\n",reference[i]-start_of_ref);
	}
	fclose(reference_file);*/
}

void printNewRefFile(FILE* new_ref_file, int number_of_identical_strains, int* identical, int length_of_MSA, int** MSA, int* invariant_sites, char** names_of_strains){
	int i,j;
	for(i=0; i<number_of_identical_strains; i++){
		fprintf(new_ref_file,">%s\n",names_of_strains[identical[i]]);
		for(j=0; j<length_of_MSA; j++){
			if ( invariant_sites[j]==-1 ){
				if ( MSA[identical[i]][j] == 0 ){
					fprintf(new_ref_file,"-");
				} else if ( MSA[identical[i]][j] == 1 ){
					fprintf(new_ref_file,"A");
				}else if (MSA[identical[i]][j] == 2 ){
					fprintf(new_ref_file,"G");
				}else if (MSA[identical[i]][j] == 3 ){
					fprintf(new_ref_file,"C");
				}else if (MSA[identical[i]][j] == 4 ){
					fprintf(new_ref_file,"T");
				}
			}
		}
		fprintf(new_ref_file,"\n");
	}
}

int findStartOfRef(int** MSA, int length, int ref_index){
	int i;
	int start=0;
	for(i=0; i<length; i++){
		if (MSA[ref_index][i] != 0 && MSA[ref_index][i] != -1){
			start=i;
			break;
		}
	}
	printf("start is %d\n",start);
	return start;
}

int findEndOfRef(int** MSA, int length, int ref_index){
	int i;
	int end_of_polyA=0;
	for(i=length-1; i>=0; i--){
		if(MSA[ref_index][i] !=0 && MSA[ref_index][i] != 1 && MSA[ref_index][i] != -1){
			end_of_polyA=i;
			break;
		}
	}
	printf("end is %d\n",end_of_polyA);
	return end_of_polyA;
}

int main(int argc, char **argv){
	struct timespec tstart={0,0}, tend={0.0};
	Options opt;
	opt.remove_identical=0;
	opt.print_variant=0;
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
	int** MSA = (int **)malloc(number_of_strains*sizeof(int *));
	int i,j;
	for(i=0; i<number_of_strains; i++){
		MSA[i] = (int *)malloc(length_of_MSA*sizeof(int));
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
	int ref_index=readInMSA(MSA_file,MSA,allele_max,names_of_strains,reference,length_of_MSA);
	fclose(MSA_file);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	int start_of_ref = findStartOfRef(MSA,length_of_MSA,ref_index);
	int end_of_ref = findEndOfRef(MSA,length_of_MSA,ref_index);
	FILE *reference_file;
	if (( reference_file = fopen(opt.reference,"w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	for(i=0; i<end_of_ref+1; i++){
		if (reference[i]==-1){ break;}
		if (reference[i] > end_of_ref){ break;}
		fprintf(reference_file,"%d\n",reference[i]-start_of_ref);
	}
	fclose(reference_file);
	int number_of_identical_strains = 0;
	int* identical = (int *)malloc(number_of_strains*sizeof(int));
	int* strains_remaining = (int*)malloc(number_of_strains*sizeof(int));
	for(i=0; i<number_of_strains; i++){
		strains_remaining[i]=-1;
	}
	if (opt.remove_identical==1){
		int* imputation = (int *)malloc(length_of_MSA*sizeof(int));
		findMaxAllele(length_of_MSA,allele_max,imputation);
		for(i=0; i<length_of_MSA; i++){
			free(allele_max[i]);
		}
		free(allele_max);
		imputeNucMat(number_of_strains,length_of_MSA,MSA,imputation);
		free(imputation);
		number_of_identical_strains=number_of_strains;
	//if (opt.remove_identical==0){
	//	for(i=0; i<number_of_strains; i++){
	//		identical[i]=i;
	//	}
	//}
	//if (opt.remove_identical==1){
		printf("Finding identical sequences\n");
		clock_gettime(CLOCK_MONOTONIC, &tstart);
		for(i=0; i<number_of_strains; i++){
			identical[i]=-1;
		}
		number_of_identical_strains=removeIdenticalStrains(number_of_strains,length_of_MSA,MSA,identical,names_of_strains,max_name_length,opt,ref_index,start_of_ref,end_of_ref,strains_remaining);
		clock_gettime(CLOCK_MONOTONIC, &tend);
		printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	}
	int* variant_sites = (int *)malloc(length_of_MSA*sizeof(int));
	if (opt.print_variant==1){
		for(i=0; i<length_of_MSA; i++){
			variant_sites[i]=-1;
		}
		printf("Finding invariant sites\n");
		findVariantSites(number_of_identical_strains,length_of_MSA,MSA,variant_sites,opt,ref_index,start_of_ref,end_of_ref,strains_remaining,reference);
		//free(variant_sites);
	}
	//FILE* new_ref_file;
	//if (( new_ref_file = fopen("out_ref_test.fasta","w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	//printNewRefFile(new_ref_file,number_of_identical_strains,identical,length_of_MSA,MSA,variant_sites,names_of_strains);
	//fclose(new_ref_file);
	free(identical);
	free(variant_sites);
	free(strains_remaining);
	//exit(1);	
	//printf("Calculating allele frequencies...\n");
	//FILE* sam_file;
	//if ((sam_file = fopen(opt.sam,"r")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	//double** allele_frequency = (double **)malloc(length_of_MSA*sizeof(double *));
	//for(i=0; i<length_of_MSA; i++){
	//	allele_frequency[i] = (double *)malloc(4*sizeof(double));
	//	for(j=0; j<4; j++){
	//		allele_frequency[i][j]=0;
	//	}
	//}
	//int number_of_strains_remaining=0;
	//number_of_strains_remaining=calculateAlleleFreq(sam_file,allele_frequency,reference,length_of_MSA,MSA,number_of_strains,names_of_strains,opt.freq,max_name_length,tstart,tend,identical,number_of_identical_strains);
	//fclose(sam_file);
	//for(i=0; i<length_of_MSA; i++){
	//	free(allele_frequency[i]);
	//}
	//free(allele_frequency);
	//int* strains_kept = (int *)malloc(number_of_strains_remaining*sizeof(int));
	//int placement=0;
	//for(i=0; i<number_of_strains; i++){
	//	if (names_of_strains[i][0] != '\0' ){
	//		strains_kept[placement]=i;
	//		placement++;
	//	}
	//}
	//printf("Creating mismatch matrix...\n");
	//clock_gettime(CLOCK_MONOTONIC, &tstart);
	//FILE *outfile;
	//if (( outfile = fopen(opt.outfile,"w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	//if ((sam_file = fopen(opt.sam,"r")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	//writeMismatchMatrix(outfile,sam_file,MSA,strains_kept,length_of_MSA,number_of_strains,number_of_strains_remaining,names_of_strains,reference);
	//fclose(sam_file);
	//fclose(outfile);
	//clock_gettime(CLOCK_MONOTONIC, &tend);
	//printf("Took %.5fseconds\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	//free(strains_kept);
	free(reference);
	for(i=0; i<number_of_strains; i++){
		free(names_of_strains[i]);
		free(MSA[i]);
	}
	free(names_of_strains);
	free(MSA);
}
