#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include "global.h"

int setMSALength(gzFile MSA_file){
	char buffer [FASTA_MAXLINE];
	int length = 0;
	int i=0;
	int iter=0;
	while( gzgets(MSA_file,buffer,FASTA_MAXLINE) != NULL ){
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

void setNumStrains(gzFile MSA_file, int* strain_info){
	char buffer [FASTA_MAXLINE];
	int i=0;
	int numstrains=0;
	int maxname=0;
	while( gzgets(MSA_file,buffer,FASTA_MAXLINE) != NULL ){
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

int readInMSA(gzFile MSA_file, char** MSA, int** allele_frequency, char** names, int length_of_MSA){
	char buffer [FASTA_MAXLINE];
	int i=0;
	int index=-1;
	int found_ref=0;
	int ref_index = 0;
	while( gzgets(MSA_file,buffer,FASTA_MAXLINE) != NULL ){
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
					MSA[index][i]='\0';
				}
			}
			if (found_ref==1){
				ref_index=index;
				/*int placement = 0;
				for(i=0; i<length_of_MSA; i++){
					if ( MSA[index][i] != '-' ){
						reference[placement] = i;
						placement++;
					}
				}*/
				found_ref=0;
			}
		}
	}
	return ref_index;
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

int calculateAlleleFreq_paired(FILE* sam, double** allele, int* reference, int length_of_MSA, char** MSA, int number_of_strains, char** names_of_strains, double freq_threshold, int maxname, struct timespec tstart, struct timespec tend, int* identical, int number_of_different_strains, int number_of_variant_sites, int* variant_sites){
	int i,j;
	char buffer [FASTA_MAXLINE];
	char *s;
	int cigar[MAX_CIGAR];
	char cigar_chars[MAX_CIGAR];
	for(i=0; i<MAX_CIGAR; i++){
		cigar[i]=0;
		cigar_chars[i]='\0';
	}
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	while( fgets(buffer,FASTA_MAXLINE,sam) != NULL ){
		if ( buffer[0] != '@'){
			char* buffer_copy = strdup(buffer);
			s = strtok(buffer,"\t");
			for(i=0; i<3; i++){
				s = strtok(NULL,"\t");
			}
			int position=0;
			sscanf(s, "%d", &position);
			position--;
			s=strtok(NULL,"\t");
			s=strtok(NULL,"\t");
			char *cigar_string;
			cigar_string = strdup(s);
			//printf("%s\n",cigar_string);
			if ( strcmp(cigar_string,"*")!=0 ){
			char *copy = strdup(cigar_string);
			char *res = strtok(cigar_string, "MID");
			int index=0;
			while(res){
				int from = res - cigar_string + strlen(res);
				//printf("%s\n",res);
				int cigar_count = 0;
				sscanf(res, "%d", &cigar_count);
				//printf("%d\n",cigar_count);
				res = strtok ( NULL, "MID");
				int to = res != NULL ? res-cigar_string : strlen(copy);
				//printf("%.*s\n", to-from, copy+from);
				char cigar_char = '\0';
				sscanf(copy+from, "%c", &cigar_char);
				//printf("cigar char: %c\n",cigar_char);
				cigar[index]=cigar_count;
				cigar_chars[index]=cigar_char;
				index++;
			}
			//printf("%s\n",cigar_string);
			free(copy);
			//printf("S: %s\n",buffer_copy);
			s = strtok(buffer_copy,"\t");
			//printf("S: %s\n",s);
			for(i=0; i<9; i++){
				s = strtok(NULL,"\t");
			}
			char* sequence = s;
			int cigar_char_count=index;
			index=0;
			int start=0;
			int start_ref=0;
			for(i=0; i< cigar_char_count; i++){
				//printf("cigar_chars[%d]: %c\n",i,cigar_chars[i]);
			}
			//printf("cigar_char_count: %d\n",cigar_char_count);
			for(i=0;i<cigar_char_count;i++){
				//printf("cigar[%d]=%d\n",i,cigar[i]);
				for(j=0; j<cigar[i]; j++){
					//printf("cigar_chars[%d]: %c\n",i,cigar_chars[i]);
					if ( cigar_chars[i] == 'M' || cigar_chars[i] == 'I'){
						if (sequence[j+start]=='A' || sequence[j+start]=='a'){
							if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
								allele[reference[j+start_ref+position]][0]++;
							}
						}else if ( sequence[j+start]=='G' || sequence[j+start]=='g'){
							if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
								allele[reference[j+start_ref+position]][1]++;
							}
						}else if ( sequence[j+start]=='C' || sequence[j+start]=='c'){
							if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
								allele[reference[j+start_ref+position]][2]++;
							}
						}else if (sequence[j+start]=='T' || sequence[j+start]=='t'){
							if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
								allele[reference[j+start_ref+position]][3]++;
							}
						}
					}
				}
				if (cigar_chars[i] == 'M'){
					start = cigar[i]+start;
					start_ref = cigar[i]+start_ref;
				}
				if (cigar_chars[i] == 'I'){
					start = cigar[i] + start;
				}
				if (cigar_chars[i] == 'D'){
					start_ref = cigar[i] + start_ref;
				}
			}
			}
			free(buffer_copy);
			free(cigar_string);
		}
	}
	for(i=0; i<length_of_MSA; i++){
		double total=0;
		for(j=0; j<4; j++){
			total=total+allele[i][j];
		}
		for(j=0; j<4; j++){
			allele[i][j] = allele[i][j]/total;
		}
	}
	int length_of_reference = 0;
	for (i=0; i<length_of_MSA; i++){
		if ( reference[i]==-1 ){
			break;
		}else{
			length_of_reference++;
		}
	}
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
			//printf("Remaining strain: %s\n",names_of_strains[i]);
			number_remaining++;
		}
	}
	printf("Number remaining: %d\n",number_remaining);
	return number_remaining;
}

int calculateAlleleFreq(FILE* sam, double** allele, int* reference, int length_of_MSA, char** MSA, int number_of_strains, char** names_of_strains, double freq_threshold, int maxname, struct timespec tstart, struct timespec tend, int* identical, int number_of_different_strains, int number_of_variant_sites, int* variant_sites){
	int i,j;
	char buffer [FASTA_MAXLINE];
	char *s;
	int cigar[MAX_CIGAR];
	char cigar_chars[MAX_CIGAR];
	for(i=0; i<MAX_CIGAR; i++){
		cigar[i]=0;
		cigar_chars[i]='\0';
	}
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	while( fgets(buffer,FASTA_MAXLINE,sam) != NULL ){
		if ( buffer[0] != '@'){
			char* buffer_copy = strdup(buffer);
			s = strtok(buffer,"\t");
			for(i=0; i<3; i++){
				s = strtok(NULL,"\t");
			}
			int position=0;
			sscanf(s, "%d", &position);
			position--;
			s=strtok(NULL,"\t");
			s=strtok(NULL,"\t");
			char *cigar_string;
			cigar_string = strdup(s);
			char *copy = strdup(cigar_string);
			char *res = strtok(cigar_string, "MID");
			int index=0;
			while(res){
				int from = res - cigar_string + strlen(res);
				int cigar_count = 0;
				sscanf(res, "%d", &cigar_count);
				res = strtok ( NULL, "MID");
				int to = res != NULL ? res-cigar_string : strlen(copy);
				char cigar_char = '\0';
				sscanf(copy+from, "%c", &cigar_char);
				cigar[index]=cigar_count;
				cigar_chars[index]=cigar_char;
				index++;
			}
			free(copy);
			free(cigar_string);
			s = strtok(buffer_copy,"\t");
			for(i=0; i<9; i++){
				s = strtok(NULL,"\t");
			}
			int cigar_char_count=index;
			int start=0;
			int start_ref=0;
			char* sequence = s;
			for(i=0; i<cigar_char_count; i++){
				for(j=0; j<cigar[i]; j++){
					if( cigar_chars[i] == 'M' ){
						if (sequence[j+start]=='A' || sequence[j+start]=='a' ){
							if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
								allele[reference[j+start_ref+position]][0]++;
							}
						}else if ( sequence[j+start]=='G' || sequence[j+start]=='g'){
							if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
								allele[reference[j+start_ref+position]][1]++;
							}
						}else if ( sequence[j+start]=='C' || sequence[j+start]=='c'){
							if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
								allele[reference[j+start_ref+position]][2]++;
							}
						}else if (sequence[j+start]=='T' || sequence[j+start]=='t'){
							if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
								allele[reference[j+start_ref+position]][3]++;
							}
						}
					}
				}
				if (cigar_chars[i] == 'M'){
					start = cigar[i]+start;
					start_ref = cigar[i]+start_ref;
				}
				if (cigar_chars[i] == 'I'){
					start = cigar[i] + start;
				}
				if (cigar_chars[i] == 'D'){
					start_ref = cigar[i] + start_ref;
				}
			}
			free(buffer_copy);
		}
	}
	for(i=0; i<length_of_MSA; i++){
		double total=0;
		for(j=0; j<4; j++){
			total=total+allele[i][j];
		}
		for(j=0; j<4; j++){
			allele[i][j] = allele[i][j]/total;
		}
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
			//printf("Remaining strain: %s\n",names_of_strains[i]);
			number_remaining++;
		}
	}
	printf("Number of strains remaining is %d\n",number_remaining);
	return number_remaining;
}

int dec2bin( int n){
	int binaryNum[32];
	int i=0;
	while(n>0){
		binaryNum[i] = n%2;
		n = n/2;
		i++;
	}
	if (binaryNum[2]==1){
		return -1;
	}else if ( binaryNum[3]==1){
		return 2;
	}else{
		return binaryNum[6];
	}
}

void writeMismatchMatrix_paired( FILE* outfile, FILE* samfile, char** MSA, int* strains_kept, int length_of_MSA, int number_of_strains, int number_of_strains_remaining, char** names_of_strains, int* reference){
	int i,j,k;
	char buffer [FASTA_MAXLINE];
	char *s;
	int cigar[MAX_CIGAR];
	char cigar_chars[MAX_CIGAR];
	for(i=0; i<MAX_CIGAR; i++){
		cigar[i]=0;
		cigar_chars[i]='\0';
	}
	int* number_of_mismatches = (int *)malloc(number_of_strains_remaining*sizeof(int));
	fprintf(outfile,"qName\tblockSizes");
	for(i=0; i<number_of_strains_remaining; i++){
		fprintf(outfile,"\t%s",names_of_strains[strains_kept[i]]);
	}
	fprintf(outfile,"\n");
	int alignment_size;
	while( fgets(buffer,FASTA_MAXLINE,samfile) != NULL ){
		if ( buffer[0] != '@'){
			char* buffer_copy = strdup(buffer);
			s = strtok(buffer,"\t");
			char* name = strdup(s);
			s = strtok(NULL, "\t");
			int decimal=0;
			sscanf(s, "%d", &decimal);
			decimal=dec2bin(decimal);
			if (decimal==1){
				fprintf(outfile,"%s",name);
			}
			if (decimal==2){
				fprintf(outfile,"%s",name);
			}
			free(name);
			for(i=0; i<2; i++){
				s = strtok(NULL,"\t");
			}
			int position=0;
			sscanf(s, "%d", &position);
			position--;
			s=strtok(NULL,"\t");
			s=strtok(NULL,"\t");
			char *cigar_string;
			cigar_string = strdup(s);
			char *copy = strdup(cigar_string);
			char *res = strtok(cigar_string, "MID");
			int index=0;
			while(res){
				int from = res - cigar_string + strlen(res);
				int cigar_count = 0;
				sscanf(res, "%d", &cigar_count);
				res = strtok ( NULL, "MID");
				int to = res != NULL ? res-cigar_string : strlen(copy);
				char cigar_char = '\0';
				sscanf(copy+from, "%c", &cigar_char);
				cigar[index]=cigar_count;
				cigar_chars[index]=cigar_char;
				index++;
			}
			free(copy);
			free(cigar_string);
			s = strtok(buffer_copy,"\t");
			for(i=0; i<9; i++){
				s = strtok(NULL,"\t");
			}
			char* sequence = s;
			int cigar_char_count=index;
			index=0;
			int start=0;
			int start_ref=0;
			if (decimal==1 || decimal==2){
				alignment_size=0;
			}
			if ( decimal==1 || decimal==2){
				for(i=0; i<number_of_strains_remaining; i++){
					number_of_mismatches[i] = 0;
				}
			}
			if (decimal != -1){
			for(i=0; i<number_of_strains_remaining; i++){
				start=0;
				start_ref=0;
					for(j=0;j<cigar_char_count;j++){
						for(k=0; k<cigar[j]; k++){
							if ( cigar_chars[j] == 'M' ){
								if ( sequence[k+start] == 'A' || sequence[k+start] == 'a'){
									if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
										if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'A' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
											number_of_mismatches[i]++;
										}
									}
								}else if ( sequence[k+start] == 'G' || sequence[k+start] == 'g'){
									if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
										if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'G' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
											number_of_mismatches[i]++;
										}
									}
								}else if ( sequence[k+start] == 'C' || sequence[k+start] == 'c'){
									if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
										if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'C' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
											number_of_mismatches[i]++;
										
									}}
								}else if ( sequence[k+start] == 'T' || sequence[k+start] == 't'){
									if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
										if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'T' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
											number_of_mismatches[i]++;
										}
									}
								}
							}
							if ( cigar_chars[j] == 'D' ){
								if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'){
										number_of_mismatches[i]++;
									}
								}
							}
							if ( cigar_chars[j] == 'I' ){
								if ( sequence[k+start] == 'A' || sequence[k+start] == 'a'){
									if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
										if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'A'){
											number_of_mismatches[i]++;
										}
									}
								}
								if ( sequence[k+start] == 'G' || sequence[k+start] == 'g'){
									if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
										if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'G'){
											number_of_mismatches[i]++;
										}
									}
								}
								if ( sequence[k+start] == 'C' || sequence[k+start] == 'c'){
									if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
										if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'C'){
											number_of_mismatches[i]++;
										}
									}
								}
								if ( sequence[k+start] == 'T' || sequence[k+start] == 't'){
									if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
										if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'T'){
											number_of_mismatches[i]++;
										}
									}
								}
							}
						}
						if (cigar_chars[j] == 'M'){
							start = cigar[j]+start;
							start_ref = cigar[j]+start_ref;
							//if (i==0){
							//	alignment_size = alignment_size + cigar[j];
							//}
						}
						if (cigar_chars[j] == 'I'){
							start = cigar[j] + start;
						}
						if (cigar_chars[j] == 'D'){
							start_ref = cigar[j] + start_ref;
						}
						if (i==0){
							alignment_size = alignment_size + cigar[j];
						}
					}
			}
			}
			if (decimal==0 || decimal==2){
				fprintf(outfile,"\t%d",alignment_size);
				for(i=0; i<number_of_strains_remaining; i++){
					fprintf(outfile,"\t%d",number_of_mismatches[i]);
				}
				fprintf(outfile,"\n");
			}
			free(buffer_copy);
		}
	}
	free(number_of_mismatches);
}

void writeMismatchMatrix( FILE* outfile, FILE* samfile, char** MSA, int* strains_kept, int length_of_MSA, int number_of_strains, int number_of_strains_remaining, char** names_of_strains, int* reference){
	int i,j,k;
	char buffer [FASTA_MAXLINE];
	char *s;
	int cigar[MAX_CIGAR];
	char cigar_chars[MAX_CIGAR];
	for(i=0; i<MAX_CIGAR; i++){
		cigar[i]=0;
		cigar_chars[i]='\0';
	}
	int* number_of_mismatches = (int *)malloc(number_of_strains_remaining*sizeof(int));
	fprintf(outfile,"qName\tblockSizes");
	for(i=0; i<number_of_strains_remaining; i++){
		fprintf(outfile,"\t%s",names_of_strains[strains_kept[i]]);
	}
	fprintf(outfile,"\n");
	int alignment_size;
	while( fgets(buffer,FASTA_MAXLINE,samfile) != NULL ){
		if ( buffer[0] != '@'){
			char* buffer_copy = strdup(buffer);
			s = strtok(buffer,"\t");
			fprintf(outfile,"%s",s);
			s = strtok(NULL, "\t");
			int decimal=0;
			sscanf(s, "%d", &decimal);
			decimal=dec2bin(decimal);
			for(i=0; i<2; i++){
				s = strtok(NULL,"\t");
			}
			int position=0;
			sscanf(s, "%d", &position);
			position--;
			for(i=0; i<2; i++){
				s = strtok(NULL,"\t");
			}
			char *cigar_string;
			cigar_string = strdup(s);
			char *copy = strdup(cigar_string);
			char *res = strtok(cigar_string, "MID");
			int index=0;
			while(res){
				int from = res - cigar_string + strlen(res);
				int cigar_count = 0;
				sscanf(res, "%d", &cigar_count);
				res = strtok ( NULL, "MID");
				int to = res != NULL ? res-cigar_string : strlen(copy);
				char cigar_char = '\0';
				sscanf(copy+from, "%c", &cigar_char);
				cigar[index]=cigar_count;
				cigar_chars[index]=cigar_char;
				index++;
			}
			free(copy);
			free(cigar_string);
			s = strtok(buffer_copy,"\t");
			for(i=0; i<9; i++){
				s = strtok(NULL,"\t");
			}
			char* sequence = s;
			int cigar_char_count=index;
			int start=0;
			int start_ref=0;
			alignment_size=0;
			for(i=0; i<number_of_strains_remaining; i++){
				number_of_mismatches[i] = 0;
			}
			free(buffer_copy);
			/*char alignment_size [MAX_CIGAR];
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
			size = strlen(sequence);*/
			for(i=0; i<number_of_strains_remaining; i++){
				start=0;
				start_ref=0;
				for(j=0; j<cigar_char_count; j++){
					for(k=0; k<cigar[j]; k++){
						if (cigar_chars[j] == 'M' ){
							if ( sequence[k+start] == 'A' || sequence[k+start] == 'a' ){
								if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'A' && MSA[strains_kept[i]][reference[k+position+start_ref]] != '-' ){
									if ( reference[k+start_ref+position] < length_of_MSA ){
										number_of_mismatches[i]++;
									}
								}
							}else if ( sequence[k+start] == 'G' || sequence[k+start] == 'g' ){
								if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'G' && MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'){
									if ( reference[k+start_ref+position] < length_of_MSA ){
										number_of_mismatches[i]++;
									}
								}
							}else if ( sequence[k+start] == 'C' || sequence[k+start] == 'c' ){
								if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'C' && MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'){
									if ( reference[k+start_ref+position] < length_of_MSA ){
										number_of_mismatches[i]++;
									}
								}
							}else if ( sequence[k+start] == 'T' || sequence[k+start] == 't' ){
								if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'T' && MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'){
									if ( reference[k+start_ref+position] < length_of_MSA ){
										number_of_mismatches[i]++;
									}
								}
							}
						}
					}
					if (cigar_chars[j] == 'M'){
						start = cigar[j]+start;
						start_ref = cigar[j]+start_ref;
						if (i==0){
							alignment_size = alignment_size + cigar[j];
						}
					}
					if (cigar_chars[j] == 'I'){
						start = cigar[j] + start;
					}
					if (cigar_chars[j] == 'D'){
						start_ref = cigar[j] + start_ref;
					}
				}
			}
			fprintf(outfile,"\t%d",alignment_size);
			for(i=0; i<number_of_strains_remaining; i++){
				fprintf(outfile,"\t%d",number_of_mismatches[i]);
			}
			fprintf(outfile,"\n");
		}
	}
	free(number_of_mismatches);
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

void readReferencePositionsFile(FILE* file, int* reference){
	char buffer[20];
	int placement=0;
	while( fgets(buffer,20,file) != NULL){
		reference[placement]=atoi(buffer);
		placement++;
	}
}

int findEndOfPolyA(char **MSA, int length_of_MSA, int ref_index, int* reference){
	int i;
	int length_of_ref = 0;
	int end_of_polyA=reference[length_of_ref-1];
	for(i=0; i<length_of_MSA; i++){
		if ( reference[i] == -1 ){
			length_of_ref=i;
			break;
		}
	}
	for(i=length_of_ref-2; i>=0; i--){
		if ( MSA[ref_index][reference[i]] != 'A' ){
			end_of_polyA = reference[i];
			break;
		}
	}
	return end_of_polyA;
}

int main(int argc, char **argv){
	struct timespec tstart={0,0}, tend={0.0};
	Options opt;
	opt.remove_identical=0;
	opt.paired=0;
	opt.error=0.005;
	parse_options(argc, argv, &opt);
	char* buffer = (char*)malloc(FASTA_MAXLINE*sizeof(char));
	memset(buffer,'\0',FASTA_MAXLINE);
	if (opt.paired==1){
		sprintf(buffer,"bowtie2 --all -f -x %s -1 %s -2 %s -S %s",opt.bowtie_reference_db,opt.forward_end_file,opt.reverse_end_file,opt.sam);
		//printf("%s\n",buffer);
		system(buffer);
	}else{
		sprintf(buffer,"bowtie2 --all -f -x %s -U %s -S %s",opt.bowtie_reference_db,opt.single_end_file,opt.sam);
		system(buffer);
	}
	free(buffer);
	gzFile MSA_file = Z_NULL;
	if (( MSA_file = gzopen(opt.fasta,"r")) == Z_NULL ) fprintf(stderr, "File could not be opened.\n");
	int length_of_MSA=0;
	length_of_MSA=setMSALength(MSA_file);
	gzclose(MSA_file);
	printf("Length of MSA: %d\n",length_of_MSA);
	if (( MSA_file = gzopen(opt.fasta,"r")) == Z_NULL ) fprintf(stderr, "File could not be opened.\n");
	int *strain_info = (int*)malloc(2*sizeof(int));
	setNumStrains(MSA_file,strain_info);
	gzclose(MSA_file);
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
	FILE* reference_file;
	if (( reference_file = fopen(opt.reference,"r")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	readReferencePositionsFile(reference_file,reference);
	fclose(reference_file);
	printf("Reading in MSA...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	if (( MSA_file = gzopen(opt.fasta,"r")) == Z_NULL ) fprintf(stderr, "File could not be opened.\n");
	int ref_index = readInMSA(MSA_file,MSA,allele_max,names_of_strains,length_of_MSA);
	gzclose(MSA_file);
	for(i=0; i<length_of_MSA; i++){
		free(allele_max[i]);
	}
	free(allele_max);
	//int end_of_polyA = findEndOfPolyA(MSA,length_of_MSA,ref_index,reference);
	//length_of_MSA = end_of_polyA;
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
	if ( opt.paired==1 ){
		number_of_strains_remaining=calculateAlleleFreq_paired(sam_file,allele_frequency,reference,length_of_MSA,MSA,number_of_strains,names_of_strains,opt.freq,max_name_length,tstart,tend,identical,number_of_identical_strains,number_of_variant_sites,variant_sites);
	}else{
		number_of_strains_remaining=calculateAlleleFreq(sam_file,allele_frequency,reference,length_of_MSA,MSA,number_of_strains,names_of_strains,opt.freq,max_name_length,tstart,tend,identical,number_of_identical_strains,number_of_variant_sites,variant_sites);
	}
	fclose(sam_file);
	for(i=0; i<length_of_MSA; i++){
		free(allele_frequency[i]);
	}
	free(allele_frequency);
	free(identical);
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
	if (opt.paired==1){
		writeMismatchMatrix_paired(outfile,sam_file,MSA,strains_kept,length_of_MSA,number_of_strains,number_of_strains_remaining,names_of_strains,reference);
	}else{
		writeMismatchMatrix(outfile,sam_file,MSA,strains_kept,length_of_MSA,number_of_strains,number_of_strains_remaining,names_of_strains,reference);
	}
	fclose(sam_file);
	fclose(outfile);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fseconds\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	free(strains_kept);
	free(reference);
	free(variant_sites);
	for(i=0; i<number_of_strains; i++){
		free(names_of_strains[i]);
		free(MSA[i]);
	}
	free(names_of_strains);
	free(MSA);
	buffer = (char*)malloc(FASTA_MAXLINE*sizeof(char));
	memset(buffer,'\0',FASTA_MAXLINE);
	sprintf(buffer,"Rscript EM.R %s %lf",opt.outfile,opt.error);
	system(buffer);
	free(buffer);
}
