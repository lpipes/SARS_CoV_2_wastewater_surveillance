#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include <assert.h>
#include <pthread.h>
#include "global.h"
#include "options.h"

char** resize_MSA;
char** resize_names_of_strains;
int* reference_index;
char** sam_results;
int setMSALength(gzFile MSA_file){
	char buffer [FASTA_MAXLINE];
	int length = 0;
	int i=0;
	int iter=0;
	while( gzgets(MSA_file,buffer,FASTA_MAXLINE) != NULL ){
		if (buffer[0] != '>'){
			for(i=0; buffer[i]!='\n'; i++){
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

int readInMSA(gzFile MSA_file, char** MSA, char** names, int length_of_MSA){
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
					//allele_frequency[i][0]++;
				}else if (buffer[i]=='G' || buffer[i]=='g'){
					MSA[index][i]='G';
					//allele_frequency[i][1]++;
				}else if (buffer[i]=='C' || buffer[i]=='c'){
					MSA[index][i]='C';
					//allele_frequency[i][2]++;
				}else if (buffer[i]=='T' || buffer[i]=='t'){
					MSA[index][i]='T';
					//allele_frequency[i][3]++;
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

int calculateAlleleFreq_paired(FILE* sam, double** allele, int length_of_MSA, char** MSA, int number_of_strains, char** names_of_strains, double freq_threshold, int maxname, struct timespec tstart, struct timespec tend, int number_of_variant_sites, int* variant_sites, int coverage, int min_strains_remaining, int max_strains_remaining, char print_counts[], int* max_sam_length){
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
	int first_in_pair=0;
	int second_in_pair=0;
	int first_seq_length = 0;
	int first_seq_cigar[MAX_CIGAR];
	char first_seq_cigar_chars[MAX_CIGAR];
	int first_seq_cigar_char_count=0;
	char* first_seq = (char*)malloc(MAX_READ_LENGTH*sizeof(char));
	memset(first_seq,'\0',MAX_READ_LENGTH);
	int first_end_pos = 0;
	int first_seq_start_pos = 0;
	int visited[MAX_READ_LENGTH];
	//memset(visited,-1,MAX_READ_LENGTH);
	for(i=0; i<MAX_READ_LENGTH; i++){
		visited[i]=-1;
	}
	while( fgets(buffer,FASTA_MAXLINE,sam) != NULL ){
		if ( buffer[0] != '@'){
			max_sam_length[1]++;
			char* buffer_copy = strdup(buffer);
			int length_of_sam = strlen(buffer);
			if ( length_of_sam > max_sam_length[0] ){
				max_sam_length[0]=length_of_sam;
			}
			s = strtok(buffer,"\t");
			//char* name = strdup(s);
			s = strtok(NULL,"\t");
			int decimal=0;
			sscanf(s, "%d", &decimal);
			decimal=dec2bin(decimal);
			if (decimal==1){
				first_in_pair = 1;
			}else if (decimal==0){
				second_in_pair = 1;
			}
			//free(name);
			for(i=0; i<2; i++){
				s = strtok(NULL, "\t");
			}
			int position=0;
			sscanf(s, "%d", &position);
			position--;
			if( first_in_pair == 1 ){
				first_seq_start_pos = position;
			}
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
				if (first_in_pair==1){
					first_seq_cigar[index]=cigar_count;
				}
				cigar_chars[index]=cigar_char;
				if (first_in_pair==1){
					first_seq_cigar_chars[index]=cigar_char;
				}
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
			if (first_in_pair==1){
				strcpy(first_seq,sequence);
				first_seq_length = strlen(first_seq);
			}
			int cigar_char_count=index;
			if (first_in_pair==1){
				first_seq_cigar_char_count=index;
			}
			double alignment_score = 0;
			int matches=0;
			int ns = 0;
			int is = 0;
			for( i=0; i<cigar_char_count; i++){
				if ( cigar_chars[i] == 'M' ){
					matches = matches + cigar[i];
				}
				if ( cigar_chars[i] == 'I'){
					ns = ns + cigar[i];
				}
				if (cigar_chars[i] == 'N' ){
					is = is + cigar[i];
				}
			}
			if ( cigar_char_count > 0 ){
				alignment_score = matches/(matches+ns+is);
			}
			index=0;
			int start=0;
			int start_ref=0;
			int k=0;
			int l=0;
			int visited_place=0;
			if (second_in_pair==1 /*&& alignment_score > 0.8*/){
				int start_ref1 =0;
				int start1 = 0;
				for(i=0; i<first_seq_cigar_char_count; i++){
					for(j=0; j<first_seq_cigar[i]; j++){
						int position_in_MSA = reference_index[start_ref1 + first_seq_start_pos + j];
						if ( position_in_MSA != -1 && position_in_MSA >= position){
						int start2=0;
						int start_ref2=0;
						for(k=0; k<cigar_char_count; k++){
							for(l=0; l<cigar[k]; l++){
								if ( position_in_MSA !=-1 && position_in_MSA == reference_index[start_ref2 + position + l] && cigar_chars[k] == 'M' && first_seq_cigar_chars[i] == 'M'){
									if ( first_seq[start1 + j] != sequence[start2 + l] && first_seq_start_pos + j <= variant_sites[number_of_variant_sites-1]){
										if ( first_seq[start1+j]=='A' || first_seq[start1+j]=='a' ){
											allele[position_in_MSA][0]--;
										}else if ( first_seq[start1+j]=='G' || first_seq[start1+j]=='g' ){
											allele[position_in_MSA][1]--;
										}else if ( first_seq[start1+j]=='C' || first_seq[start1+j]=='c' ){
											allele[position_in_MSA][2]--;
										}else if ( first_seq[start1+j]=='T' || first_seq[start1+j]=='t' ){
											allele[position_in_MSA][3]--;
										}
									}
									visited[visited_place]=position_in_MSA;
									visited_place++;
								}
							}
							if ( cigar_chars[k] == 'M' ){
								start2 = start2+cigar[k];
								start_ref2 = start_ref2 + cigar[k];
							}
							if ( cigar_chars[k] == 'I' ){
								start2 = start2 + cigar[k];
							}
							if ( cigar_chars[k] == 'D' ){
								start_ref2 = start_ref2 + cigar[k];
							}
						}
						}
					}
					if (first_seq_cigar_chars[i] == 'M'){
						start1 = start1+first_seq_cigar[i];
						start_ref1 = start_ref1 + first_seq_cigar[i];
					}
					if (first_seq_cigar_chars[i] == 'I'){
						start1 = start1+first_seq_cigar[i];
					}
					if (first_seq_cigar_chars[i] == 'D'){
						start_ref1 = first_seq_cigar[i] + start_ref1;
					}
				}
			}
			//if ( alignment_score > 0.8){
			for(i=0;i<cigar_char_count;i++){
				//printf("cigar[%d]=%d\n",i,cigar[i]);
				for(j=0; j<cigar[i]; j++){
					//printf("cigar_chars[%d]: %c\n",i,cigar_chars[i]);
					int skip=0;
					int position_in_MSA = reference_index[j+start_ref+position];
					for(k=0; k<visited_place; k++){
						if ( position_in_MSA != -1 && visited[k]==position_in_MSA ){
							skip=1;
						}
					}
					if ( cigar_chars[i] == 'M' /*|| cigar_chars[i] == 'I'*/){
						if (sequence[j+start]=='A' || sequence[j+start]=='a'){
							//if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
							//	allele[reference[j+start_ref+position]][0]++;
							//}
							if ( position_in_MSA != -1 && position_in_MSA < length_of_MSA && skip==0){
								allele[position_in_MSA][0]++;
							}
						}else if ( sequence[j+start]=='G' || sequence[j+start]=='g'){
							//if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
							//	allele[reference[j+start_ref+position]][1]++;
							//}
							if ( position_in_MSA != -1 && position_in_MSA<length_of_MSA && skip==0){
								allele[position_in_MSA][1]++;
							}
						}else if ( sequence[j+start]=='C' || sequence[j+start]=='c'){
							//if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
							//	allele[reference[j+start_ref+position]][2]++;
							//}
							if ( position_in_MSA != -1 && position_in_MSA< length_of_MSA && skip==0){
								allele[position_in_MSA][2]++;
							}
						}else if (sequence[j+start]=='T' || sequence[j+start]=='t'){
							//if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
							//	allele[reference[j+start_ref+position]][3]++;
							//}
							if ( position_in_MSA != -1 && position_in_MSA< length_of_MSA && skip==0){
								allele[position_in_MSA][3]++;
							}
						}
						first_end_pos = position_in_MSA;
					}
					if ( cigar_chars[i] == "I" ){
						first_end_pos = position_in_MSA;
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
			//}
			free(buffer_copy);
			free(cigar_string);
			if ( second_in_pair == 1 ){
				memset(first_seq,'\0',MAX_READ_LENGTH);
				first_seq_cigar_char_count=0;
			}
			first_in_pair = 0;
			second_in_pair = 0;
		}
	}
	free(first_seq);
	int covered=0;
	for(i=0; i<length_of_MSA; i++){
		double total=0;
		for(j=0;j<4; j++){
			total=total+allele[i][j];
		}
		if (total>0){
			covered++;
		}
	}
	int* covered_sites = (int *)malloc(covered*sizeof(int));
	for(i=0; i<covered; i++){
		covered_sites[i]=-1;
	}
	int k=0;
	if ( print_counts[0] != '\0' ){
		FILE* allele_counts_file;
		if (( allele_counts_file = fopen(print_counts,"w")) == (FILE *) NULL ) fprintf(stderr, "Allele Counts file could not be opened.\n");	
		fprintf(allele_counts_file,"position\tA\tG\tC\tT\n");
		for(i=0; i<length_of_MSA; i++){
			fprintf(allele_counts_file,"%d\t%lf\t%lf\t%lf\t%lf\n",i,allele[i][0],allele[i][1],allele[i][2],allele[i][3]);	
		}
		fclose(allele_counts_file);
	}
	for(i=0; i<length_of_MSA; i++){
		double total=0;
		for(j=0; j<4; j++){
			total=total+allele[i][j];
		}
		if (total>= coverage){
			covered_sites[k]=i;
			k++;
		}
		for(j=0; j<4; j++){
			allele[i][j] = allele[i][j]/total;
		}
	}
	printf("Number of sites not covered: %d\n",length_of_MSA-k);
	for(i=0; i<number_of_variant_sites; i++){
		int found=0;
		for(j=0; j<covered; j++){
			if (variant_sites[i] == covered_sites[j]){
				found=1;
			}
		}
		if (found==0){
			variant_sites[i]=-1;
		}
	}
	free(covered_sites);
	int temp_num_var_sites=0;
	for(i=0; i<number_of_variant_sites; i++){
		if ( variant_sites[i] >= 0 && variant_sites[i] < length_of_MSA-1){
			temp_num_var_sites++;
		}
	}
	int* variant_sites_updated = (int*)malloc(temp_num_var_sites*sizeof(int));
	k=0;
	for(i=0; i< number_of_variant_sites; i++){
		if ( variant_sites[i] != -1 && variant_sites[i] < length_of_MSA-1){
			variant_sites_updated[k]=variant_sites[i];
			k++;
		}
	}
	free(variant_sites);
	int** bad_bases = (int **)malloc(length_of_MSA*sizeof(int *));
	for(i=0; i<length_of_MSA; i++){
		bad_bases[i] = (int *)malloc(4*sizeof(int));
		for(j=0; j<4; j++){
			bad_bases[i][j]=0;
		}
	}
	for(i=0; i<length_of_MSA; i++){
		for(j=0; j<4; j++){
			if ( allele[i][j] < freq_threshold ){
				bad_bases[i][j] = 1;
			}
		}
	}
	int* bad_bases_count = (int *)malloc(length_of_MSA*sizeof(int));
	for(i=0; i<length_of_MSA; i++){
		bad_bases_count[i] = 0;
	}
	for(i=0; i<length_of_MSA; i++){
		for(j=0; j<4; j++){
			bad_bases_count[i] += bad_bases[i][j];
		}
	}
	char** bad_base_char = (char**)malloc(length_of_MSA*sizeof(char*));
	for(i=0; i<length_of_MSA; i++){
		bad_base_char[i] = (char *)malloc(4*sizeof(char));
		for(j=0; j<4; j++){
			bad_base_char[i][j] = '\0';
		}
	}

	for(i=0; i<length_of_MSA; i++){
		j=0;
		//for(j=0; j<4-bad_bases_count[i]; j++){
			for(k=0; k<4; k++){
				if (bad_bases[i][k]==0){
					if ( k==0 ){
						bad_base_char[i][j] = 'A';
						j++;
					}else if (k==1){
						bad_base_char[i][j] = 'G';
						j++;
					}else if (k==2){
						bad_base_char[i][j] = 'C';
						j++;
					}else if (k==3){
						bad_base_char[i][j] = 'T';
						j++;
					}
					bad_bases[i][k]=1;
				}
			}
		//}
	}
	for(i=0; i<length_of_MSA; i++){
		free(bad_bases[i]);
	}
	free(bad_bases);
	//int length_of_reference = 0;
	//for (i=0; i<length_of_MSA; i++){
	//	if ( reference[i]==-1 ){
	//		break;
	//	}else{
	//		length_of_reference++;
	//	}
	//}
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	int number_remaining=number_of_strains;
	int base;
	int count;
	int var_count=0;
	int number_removed=0;
	int run_loop=1;
	int number_of_iterations=1;
	int* incompat_counter = (int*)malloc(number_of_strains*sizeof(int));
	for(i=0; i<number_of_strains; i++){
		incompat_counter[i]=0;
	}
	while(run_loop==1){
		printf("iteration %d\n",number_of_iterations);
		for(i=0; i<number_of_strains; i++){
			incompat_counter[i]=0;
		}
		number_remaining=number_of_strains;
	for(i=0; i<temp_num_var_sites; i++){
		number_removed=0;
		for(j=0; j<number_of_strains; j++){
				/*if ( MSA[identical[j]][reference[i]] != '-' ){
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
				}*/
			//if ( names_of_strains[j][0] != '\0'){
			if (incompat_counter[j] < number_of_iterations){
				count=0;
				for(k=0; k<4-bad_bases_count[variant_sites_updated[i]]; k++){
					if ( MSA[j][variant_sites_updated[i]] != bad_base_char[variant_sites_updated[i]][k] ){
						count++;
					}
				}
				if ( count==(4-bad_bases_count[variant_sites_updated[i]]) ){
					/*if ( number_remaining ==1 ){
						printf("site %d\n",variant_sites_updated[i]);
						printf("last one: %s\n",names_of_strains[j]);
						printf("MSA[%d][%d]: %c\n",j,variant_sites_updated[i],MSA[j][variant_sites_updated[i]]);
						printf("bad_base_char[%d][0]: %c\n",variant_sites_updated[i],bad_base_char[variant_sites_updated[i]][0]);
					}*/
					//printf("removing %d: %s at site %d\n",j,names_of_strains[j],variant_sites_updated[i]);
					//if ( strcmp(names_of_strains[j],"EPI_ISL_4510987")==0){
					//	printf("removing %d: %s at site %d\n",j,names_of_strains[j],variant_sites_updated[i]);
					//}
					//names_of_strains[j][0] = '\0';
					//if ( incompat_counter[j]==number_of_iterations){
						//number_remaining--;
						//number_removed++;
						incompat_counter[j]++;
					//}
				}	
			}
		}
	}
	for(i=0; i<number_of_strains; i++){
		if ( incompat_counter[i] == number_of_iterations ){
			number_remaining--;
			number_removed++;
		}
	}
	if (number_remaining >= min_strains_remaining && number_remaining < max_strains_remaining){
		printf("exiting loop. %d remaining\n",number_remaining);
		run_loop=0;
	}else if ( number_remaining >= max_strains_remaining ){
		printf("%d strains remaining. exiting...\n",number_remaining);
		exit(1);	
	}else{
		printf("there are %d remaining... \n",number_remaining);
		number_of_iterations++;
	}
	}
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	for(i=0; i<number_of_strains; i++){
		if ( incompat_counter[i] == number_of_iterations ){
			names_of_strains[i][0]='\0';
		}
	}
	for(i=0; i<length_of_MSA; i++){
		free(bad_base_char[i]);
	}
	free(bad_base_char);
	free(bad_bases_count);
	free(variant_sites_updated);
	number_remaining=0;
	for(i=0; i<number_of_strains; i++){
		if ( names_of_strains[i][0] != '\0' ){
			//printf("Remaining strain: %s\n",names_of_strains[i]);
			number_remaining++;
		}
	}
	printf("Number remaining: %d\n",number_remaining);
	free(incompat_counter);
	return number_remaining;
}

int calculateAlleleFreq(FILE* sam, double** allele, int length_of_MSA, char** MSA, int number_of_strains, char** names_of_strains, double freq_threshold, int maxname, struct timespec tstart, struct timespec tend, int number_of_variant_sites, int* variant_sites, int coverage, char print_counts[]){
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
							//if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
							//	allele[reference[j+start_ref+position]][0]++;
							//}
							if ( j+start_ref+position < length_of_MSA ){
								allele[j+start_ref+position][0]++;
							}
						}else if ( sequence[j+start]=='G' || sequence[j+start]=='g'){
							//if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
							//	allele[reference[j+start_ref+position]][1]++;
							//}
							if ( j+start_ref+position < length_of_MSA ){
								allele[j+start_ref+position][1]++;
							}
						}else if ( sequence[j+start]=='C' || sequence[j+start]=='c'){
							//if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
							//	allele[reference[j+start_ref+position]][2]++;
							//}
							if ( j+start_ref+position < length_of_MSA ){
								allele[j+start_ref+position][2]++;
							}
						}else if (sequence[j+start]=='T' || sequence[j+start]=='t'){
							//if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
							//	allele[reference[j+start_ref+position]][3]++;
							//}
							if ( j+start_ref+position < length_of_MSA ){
								allele[j+start_ref+position][3]++;
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
	if ( print_counts[0] != '\0' ){
		FILE* allele_counts_file;
		if (( allele_counts_file = fopen(print_counts,"w")) == (FILE *) NULL ) fprintf(stderr, "Allele Counts file could not be opened.\n");	
		fprintf(allele_counts_file,"position\tA\tG\tC\tT\n");
		for(i=0; i<length_of_MSA; i++){
			fprintf(allele_counts_file,"%d\t%lf\t%lf\t%lf\t%lf\n",i,allele[i][0],allele[i][1],allele[i][2],allele[i][3]);	
		}
		fclose(allele_counts_file);
	}
	int covered=0;
	for(i=0; i<length_of_MSA; i++){
		double total=0;
		for(j=0; j<4; j++){
			total=total+allele[i][j];
		}
		if ( total>0 ){
			covered++;
		}
	}
	int* covered_sites = (int *)malloc(covered*sizeof(int));
	for(i=0; i<covered; i++){
		covered_sites[i] = -1;
	}
	int k=0;
	for(i=0; i<length_of_MSA; i++){
		double total=0;
		for(j=0; j<4; j++){
			total=total+allele[i][j];
		}
		if (total >= coverage){
			covered_sites[k]=i;
			k++;
		}
		for(j=0; j<4; j++){
			allele[i][j] = allele[i][j]/total;
		}
	}
	printf("Number of sites not covered: %d\n",length_of_MSA-k);
	for(i=0; i<number_of_variant_sites; i++){
		int found=0;
		for(j=0; j<covered; j++){
			if (variant_sites[i] == covered_sites[j]){
				found=1;
			}
		}
		if (found==0){
			variant_sites[i]=-1;
		}
	}
	free(covered_sites);
	int temp_num_var_sites = 0;
	for(i=0; i<number_of_variant_sites; i++){
		if ( variant_sites[i] >= 0 && variant_sites[i] < length_of_MSA-1){
			temp_num_var_sites++;
		}
	}
	int* variant_sites_updated = (int*)malloc(temp_num_var_sites*sizeof(int));
	k=0;
	for(i=0; i< number_of_variant_sites; i++){
		if ( variant_sites[i] != -1 && variant_sites[i] < length_of_MSA-1){
			variant_sites_updated[k]=variant_sites[i];
			k++;
		}
	}
	free(variant_sites);
	//int length_of_reference = 0;
	//for (i=0; i<length_of_MSA; i++){
	//	if ( reference[i]==-1 ){
	//		break;
	//	}else{
	//		length_of_reference++;
	//	}
	//}
	int** bad_bases = (int **)malloc(length_of_MSA*sizeof(int *));
	for(i=0; i<length_of_MSA; i++){
		bad_bases[i] = (int *)malloc(4*sizeof(int));
		for(j=0; j<4; j++){
			bad_bases[i][j]=0;
		}
	}
	for(i=0; i<length_of_MSA; i++){
		for(j=0; j<4; j++){
			if ( allele[i][j] < freq_threshold ){
				bad_bases[i][j] = 1;
			}
		}
	}
	int* bad_bases_count = (int *)malloc(length_of_MSA*sizeof(int));
	for(i=0; i<length_of_MSA; i++){
		bad_bases_count[i] = 0;
	}
	for(i=0; i<length_of_MSA; i++){
		for(j=0; j<4; j++){
			bad_bases_count[i] += bad_bases[i][j];
		}
	}
	char** bad_base_char = (char**)malloc(length_of_MSA*sizeof(char*));
	for(i=0; i<length_of_MSA; i++){
		bad_base_char[i] = (char *)malloc(4*sizeof(char));
		for(j=0; j<4; j++){
			bad_base_char[i][j] = '\0';
		}
	}

	for(i=0; i<length_of_MSA; i++){
		j=0;
		//for(j=0; j<4-bad_bases_count[i]; j++){
			for(k=0; k<4; k++){
				if (bad_bases[i][k]==0){
					if ( k==0 ){
						bad_base_char[i][j] = 'A';
						j++;
					}else if (k==1){
						bad_base_char[i][j] = 'G';
						j++;
					}else if (k==2){
						bad_base_char[i][j] = 'C';
						j++;
					}else if (k==3){
						bad_base_char[i][j] = 'T';
						j++;
					}
					bad_bases[i][k]=1;
				}
			}
		//}
	}
	for(i=0; i<length_of_MSA; i++){
		free(bad_bases[i]);
	}
	free(bad_bases);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	printf("Eliminating strains...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	int number_remaining=0;
	int base;
	int count;
	int var_count=0;
	/*int* nums_of_strains = (int*)malloc(number_of_different_strains*sizeof(int));
	for(i=0; i<number_of_different_strains; i++){
		nums_of_strains[i]=-1;
	}
	int place=0;*/
	for(i=0; i<temp_num_var_sites; i++){
			//if ( variant_sites[i] == covered_sites[next] ){
		for(j=0; j<number_of_strains; j++){
			if ( names_of_strains[j][0] != '\0'){
				/*if ( MSA[j][variant_sites_updated[i]] == 'A' ){
					base=0;
				}else if ( MSA[j][variant_sites_updated[i]] == 'G' ){
					base=1;
				}else if ( MSA[j][variant_sites_updated[i]] == 'C' ){
					base=2;
				}else if ( MSA[j][variant_sites_updated[i]] == 'T' ){
					base=3;
				}
				if ( allele[variant_sites_updated[i]][base] < freq_threshold ){
					//memset(names_of_strains[identical[j]],'\0',maxname);
					names_of_strains[j][0] = '\0';
				}*/
			count=0;
			for(k=0; k<4-bad_bases_count[variant_sites_updated[i]]; k++){
				if ( MSA[j][variant_sites_updated[i]] != bad_base_char[variant_sites_updated[i]][k] ){
					count++;
				}
			}
			if ( count==(4-bad_bases_count[variant_sites_updated[i]]) ){
				names_of_strains[j][0] = '\0';
			}
			}
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	for(i=0; i<length_of_MSA; i++){
		free(bad_base_char[i]);
	}
	free(bad_base_char);
	free(bad_bases_count);
	free(variant_sites_updated);
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
		//unpaired
		return 2;
	}else{
		//if 1 this is first pair, if 0 this is second pair
		return binaryNum[6];
	}
}

void *writeMismatchMatrix_paired(void *ptr){
	int i,j,k;
	struct thread_struct *tstr = (thread_struct *) ptr;
	resultsStruct *results=tstr->str;
	int max_sam_length=tstr->max_sam_length;
	int length_of_MSA=tstr->length_of_MSA;
	int number_of_strains=tstr->number_of_strains;
	int number_of_strains_remaining=tstr->number_of_strains_remaining;
	int thread_number=tstr->thread_number;
	char buffer [FASTA_MAXLINE];
	char *s;
	int cigar[MAX_CIGAR];
	char cigar_chars[MAX_CIGAR];
	for(i=0; i<MAX_CIGAR; i++){
		cigar[i]=0;
		cigar_chars[i]='\0';
	}
	int* number_of_mismatches = (int *)malloc(number_of_strains_remaining*sizeof(int));
	int alignment_size;
	int first_in_pair=0;
	int second_in_pair=0;
	int first_seq_cigar[MAX_CIGAR];
	char first_seq_cigar_chars[MAX_CIGAR];
	int first_seq_cigar_char_count=0;
	char* first_seq = (char*)malloc(MAX_READ_LENGTH*sizeof(char));
	memset(first_seq,'\0',MAX_READ_LENGTH);
	int first_seq_start_pos = 0;
	int visited[MAX_READ_LENGTH];
	int visited_place=0;
	//memset(visited,-1,MAX_READ_LENGTH);
	for(i=0; i<MAX_READ_LENGTH; i++){
		visited[i]=-1;
	}
	int line_number=0;
	//while( fgets(buffer,FASTA_MAXLINE,samfile) != NULL ){
	int line_count;
	char* resultsPath = (char*)malloc((max_sam_length+300000)*sizeof(char));
	int index_mismatch=0;
	char* context=NULL;
	for(line_count=tstr->start; line_count < tstr->end; line_count++){
		line_number++;
		char* buffer_copy = strdup(sam_results[line_count]);
		s = strtok_r(sam_results[line_count],"\t",&context);
		char* name = strdup(s);
		s = strtok_r(NULL, "\t", &context);
		int decimal=0;
		sscanf(s, "%d", &decimal);
		decimal=dec2bin(decimal);
		if (decimal==1){
			for(i=0; i<max_sam_length+300000; i++){
				resultsPath[i]='\0';
			}
			strcpy(resultsPath,name);
			strcat(resultsPath,"\t");
			first_in_pair=1;
		}else if (decimal==0){
			second_in_pair=1;
		}
		if (decimal==2){
			for(i=0; i<max_sam_length+300000; i++){
				resultsPath[i]='\0';
			}
			strcpy(resultsPath,name);
			strcat(resultsPath,"\t");
		}
		free(name);
		for(i=0; i<2; i++){
			s = strtok_r(NULL,"\t",&context);
		}
		int position=0;
		sscanf(s, "%d", &position);
		position--;
		if( first_in_pair == 1 ){
			first_seq_start_pos = position;
		}
		s=strtok_r(NULL,"\t",&context);
		s=strtok_r(NULL,"\t",&context);
		char *cigar_string;
		cigar_string = strdup(s);
		char *copy = strdup(cigar_string);
		char *res = strtok_r(cigar_string, "MID", &context);
		int index=0;
		while(res){
			int from = res - cigar_string + strlen(res);
			int cigar_count = 0;
			sscanf(res, "%d", &cigar_count);
			res = strtok_r ( NULL, "MID", &context);
			int to = res != NULL ? res-cigar_string : strlen(copy);
			char cigar_char = '\0';
			sscanf(copy+from, "%c", &cigar_char);
			cigar[index]=cigar_count;
			if ( first_in_pair==1){
				first_seq_cigar[index]=cigar_count;
			}
			cigar_chars[index]=cigar_char;
			if (first_in_pair==1){
				first_seq_cigar_chars[index]=cigar_char;
			}
			index++;
		}
		free(copy);
		free(cigar_string);
		s = strtok_r(buffer_copy,"\t",&context);
		for(i=0; i<9; i++){
			s = strtok_r(NULL,"\t",&context);
		}
		char* sequence = s;
		if(first_in_pair==1){
			strcpy(first_seq,sequence);
			first_seq_cigar_char_count=index;
		}
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
		int l=0;
		int m=0;
		if (decimal != -1){
			visited_place=0;
			if (second_in_pair==1){
				int start_ref1=0;
				int start1=0;
				for(i=0; i<number_of_strains_remaining; i++){
					visited_place=0;
					start1=0;
					start_ref1=0;
					for(j=0; j<first_seq_cigar_char_count; j++){
						for(k=0;k<first_seq_cigar[j]; k++){
							int start2=0;
							int start_ref2=0;
							if (start_ref1 + first_seq_start_pos + k >= position){
								for(l=0; l<cigar_char_count; l++){
									int position_in_MSA = reference_index[start_ref1+first_seq_start_pos + k];
									for(m=0; m<cigar[l]; m++){
										if ( start_ref1 + first_seq_start_pos + k == start_ref2 + position + m && first_seq_cigar_chars[j]=='M'){
											if ( cigar_chars[l] == 'M'){
											if ( first_seq[start1 + k] != sequence[start2 + m] ){
												if ( first_seq[start1+k]=='A' || first_seq[start1+k]=='a'){
													if ( resize_MSA[i][position_in_MSA] != 'A' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0' ){
														if ( position_in_MSA < length_of_MSA){
															number_of_mismatches[i]--;
															assert(number_of_mismatches[i] >= 0);
														}
													}
												}else if ( first_seq[start1+k]=='G' || first_seq[start1+k]=='g'){
													if ( resize_MSA[i][position_in_MSA] != 'G' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0' ){
														if ( position_in_MSA < length_of_MSA){
															number_of_mismatches[i]--;
															assert(number_of_mismatches[i] >= 0);
														}
													}
												}else if ( first_seq[start1+k]=='C' || first_seq[start1+k]=='c'){
													if ( resize_MSA[i][position_in_MSA] != 'C' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0' ){
														if ( position_in_MSA < length_of_MSA){
															number_of_mismatches[i]--;
															assert(number_of_mismatches[i] >= 0);
														}
													}
												}else if ( first_seq[start1+k]=='T' || first_seq[start1+k]=='t'){
													if ( resize_MSA[i][position_in_MSA] != 'T' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0' ){
														if ( position_in_MSA < length_of_MSA){
															number_of_mismatches[i]--;
															assert(number_of_mismatches[i] >= 0);
														}
													}
												}
											}
											}
												/*if (cigar_chars[l]== 'D'){
													if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0'){
														if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
															number_of_mismatches[i]--;
														}
													}
												}*/
												/*if (cigar_chars[l] == 'I'){
													if (first_seq[k+start1] == 'A' || first_seq[k+start1] == 'a'){
														if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'A' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
															if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
																number_of_mismatches[i]--;
															}
														}
													}else if ( first_seq[start1+k]=='G' || first_seq[start1+k]=='g'){
														if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'G' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
															if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
																number_of_mismatches[i]--;
															}
														}
													}else if ( first_seq[start1+k]=='C' || first_seq[start1+k]=='c'){
														if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'C' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
															if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
																number_of_mismatches[i]--;
															}
														}
													}else if ( first_seq[start1+k]=='T' || first_seq[start1+k]=='t'){
														if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'T' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
															if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
																number_of_mismatches[i]--;
															}
														}
													}
												}*/
												if (cigar_chars[l]=='M'){
												visited[visited_place]=start_ref1 + first_seq_start_pos + k;
												visited_place++;
												}
											}
										}
									//}
									if ( cigar_chars[l]=='M'){
										start2 = start2+cigar[l];
										start_ref2 = start_ref2 + cigar[l];
									}
									if ( cigar_chars[l] == 'I'){
										start2 = start2 + cigar[l];
									}
									if ( cigar_chars[l] == 'D'){
										start_ref2 = start_ref2 + cigar[l];
									}
								}
								}
								//if ( first_seq_cigar_chars[j] != 'I' ){
								//	start_ref1 = start_ref1+first_seq_cigar[j];
								//	start1 = start1 + first_seq_cigar[j];
								//}
							}
							if (first_seq_cigar_chars[j] == 'M'){
								start1 = start1+first_seq_cigar[j];
								start_ref1 = start_ref1 + first_seq_cigar[j];
							}
							if (first_seq_cigar_chars[j] == 'I'){
								start1 = start1+first_seq_cigar[j];
							}
							if (first_seq_cigar_chars[j] == 'D'){
								start_ref1 = first_seq_cigar[j] + start_ref1;
							}
						}
					}
				}
				for(i=0; i<number_of_strains_remaining; i++){
				start=0;
				start_ref=0;
					for(j=0;j<cigar_char_count;j++){
						for(k=0; k<cigar[j]; k++){
							int skip=0;
							int position_in_MSA = reference_index[k+position+start_ref]; 
							for(l=0; l<visited_place; l++){
								if ( position_in_MSA == visited[l] ){
									skip=1;
								}
							}
							if ( cigar_chars[j] == 'M' ){
								if ( sequence[k+start] == 'A' || sequence[k+start] == 'a'){
									//if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'A' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
									//		number_of_mismatches[i]++;
									//	}
									//}
									if ( resize_MSA[i][position_in_MSA] != 'A' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0'){
										if ( position_in_MSA < length_of_MSA && skip==0){
											number_of_mismatches[i]++;
										}
									}
								}else if ( sequence[k+start] == 'G' || sequence[k+start] == 'g'){
									//if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'G' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
									//		number_of_mismatches[i]++;
									//	}
									//}
									if (resize_MSA[i][position_in_MSA] != 'G' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0'){
										if ( position_in_MSA < length_of_MSA && skip==0){
											number_of_mismatches[i]++;
										}
									}
								}else if ( sequence[k+start] == 'C' || sequence[k+start] == 'c'){
									//if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'C' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
									//		number_of_mismatches[i]++;
									//	
									//}}
									if (resize_MSA[i][position_in_MSA] != 'C' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0'){
										if ( position_in_MSA < length_of_MSA && skip==0){
											number_of_mismatches[i]++;
										}
									}
								}else if ( sequence[k+start] == 'T' || sequence[k+start] == 't'){
									//if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'T' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
									//		number_of_mismatches[i]++;
									//	}
									//}
									if (resize_MSA[i][position_in_MSA] != 'T' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0'){
										if ( position_in_MSA < length_of_MSA && skip==0){
											number_of_mismatches[i]++;
										}
									}
								}
							}
							if ( cigar_chars[j] == 'D' ){
								//if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'){
								//		number_of_mismatches[i]++;
								//	}
								//}
								if (resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0'){
									if ( position_in_MSA < length_of_MSA && skip==0){
									//	number_of_mismatches[i]++;
									}
								}
							}
							if ( cigar_chars[j] == 'I' ){
								if ( sequence[k+start] == 'A' || sequence[k+start] == 'a'){
									//if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'A'){
									//		number_of_mismatches[i]++;
									//	}
									//}
									if ( resize_MSA[i][position_in_MSA] != 'A' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0'){
										if ( position_in_MSA < length_of_MSA && skip==0){
											number_of_mismatches[i]++;
										}
									}
								}
								if ( sequence[k+start] == 'G' || sequence[k+start] == 'g'){
									//if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'G'){
									//		number_of_mismatches[i]++;
									//	}
									//}
									if (resize_MSA[i][position_in_MSA] != 'G' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0'){
										if ( position_in_MSA < length_of_MSA && skip==0){
											number_of_mismatches[i]++;
										}
									}
								}
								if ( sequence[k+start] == 'C' || sequence[k+start] == 'c'){
									//if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'C'){
									//		number_of_mismatches[i]++;
									//	}
									//}
									if (resize_MSA[i][position_in_MSA] != 'C' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0'){
										if ( position_in_MSA < length_of_MSA && skip==0){
											number_of_mismatches[i]++;
										}
									}
								}
								if ( sequence[k+start] == 'T' || sequence[k+start] == 't'){
									//if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'T'){
									//		number_of_mismatches[i]++;
									//	}
									//}
									if (resize_MSA[i][position_in_MSA] != 'T' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0'){
										if ( position_in_MSA < length_of_MSA && skip==0){
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
						//if (i==0){
						//	alignment_size = alignment_size + cigar[j];
						//}
					}
			}
			}
			if (decimal==0 || decimal==2){
				if ( visited_place > 0 ){
					alignment_size = alignment_size-visited_place;
				}
				//fprintf(outfile,"\t%d",alignment_size);
				char *num = NULL;
				asprintf(&num,"%d",alignment_size);
				strcat(resultsPath,num);
				free(num);
				for(i=0; i<number_of_strains_remaining; i++){
					//fprintf(outfile,"\t%d",number_of_mismatches[i]);
					char *num2 = NULL;
					asprintf(&num2,"\t%d",number_of_mismatches[i]);
					strcat(resultsPath,num2);
					free(num2);
				}
				//fprintf(outfile,"\n");
				strcpy(results->mismatch[index_mismatch],resultsPath);
				index_mismatch++;
			}
			free(buffer_copy);
			first_in_pair=0;
			second_in_pair=0;
	}
	free(number_of_mismatches);
}

void writeMismatchMatrix( FILE* outfile, FILE* samfile, char** MSA, int* strains_kept, int length_of_MSA, int number_of_strains, int number_of_strains_remaining, char** names_of_strains){
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
								//if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'A' && MSA[strains_kept[i]][reference[k+position+start_ref]] != '-' ){
								if ( MSA[strains_kept[i]][k+position+start_ref] != 'A' && MSA[strains_kept[i]][k+position+start_ref] != '-' && MSA[strains_kept[i]][k+position+start_ref] != '\0'){
									//if ( reference[k+start_ref+position] < length_of_MSA ){
									if ( k+start_ref+position < length_of_MSA ){
										number_of_mismatches[i]++;
									}
								}
							}else if ( sequence[k+start] == 'G' || sequence[k+start] == 'g' ){
								//if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'G' && MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'){
								if (MSA[strains_kept[i]][k+position+start_ref] != 'G' && MSA[strains_kept[i]][k+position+start_ref] != '-' && MSA[strains_kept[i]][k+position+start_ref] != '\0'){
									if ( k+start_ref+position < length_of_MSA ){
										number_of_mismatches[i]++;
									}
								}
							}else if ( sequence[k+start] == 'C' || sequence[k+start] == 'c' ){
								//if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'C' && MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'){
								if (MSA[strains_kept[i]][k+position+start_ref] != 'C' && MSA[strains_kept[i]][k+position+start_ref] != '-' && MSA[strains_kept[i]][k+position+start_ref] != '\0'){
									//if ( reference[k+start_ref+position] < length_of_MSA ){
									if ( k+start_ref+position < length_of_MSA ){
										number_of_mismatches[i]++;
									}
								}
							}else if ( sequence[k+start] == 'T' || sequence[k+start] == 't' ){
								//if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'T' && MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'){
								if (MSA[strains_kept[i]][k+position+start_ref] != 'T' && MSA[strains_kept[i]][k+position+start_ref] != '-' && MSA[strains_kept[i]][k+position+start_ref] != '\0'){
									//if ( reference[k+start_ref+position] < length_of_MSA ){
									if ( k+start_ref+position < length_of_MSA ){
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
int process_problematic_sites(int* problematic_sites){
	FILE* file;
	if (( file = fopen("problematic_sites_sarsCov2.vcf","r")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	char buffer[1000];
	char name[30];
	int position;
	char ch1[1];
	char ch2[1];
	char ch3[1];
	char ch4[1];
	char s1[10];
	char s2[30];
	int i=0;
	while( fgets(buffer,1000,file) != NULL){
		if ( buffer[0] != "#" ){
			sscanf(buffer,"%s\t%d\t%c\t%c\t%c\t%c\t%s\t%s",&name,&position,&ch1,&ch2,&ch3,&ch4,&s1,&s2);	
			problematic_sites[i]=position;
			i++;	
		}
	}	
	fclose(file);
	return i;
}
void reallocate_memory(int number_of_strains_remaining, int* strains_kept, char** MSA, char** names_of_strains, int maxname, int number_of_total_strains){
	int i,j;
	int counter=0;
	for(i=0; i<number_of_strains_remaining; i++){
		strcpy(resize_names_of_strains[counter],names_of_strains[strains_kept[i]]);
		strcpy(resize_MSA[counter],MSA[strains_kept[i]]);
		counter++;
	}
	for(i=0; i<number_of_total_strains; i++){
		free(names_of_strains[i]);
		free(MSA[i]);
	}
	free(names_of_strains);
	free(MSA);
}
void readInSamFile(FILE* sam_file){
	char buffer[FASTA_MAXLINE];
	int i,j,k;
	i=0;
	while( fgets(buffer,FASTA_MAXLINE,sam_file) != NULL ){
		if ( buffer[0] != '@' ){
			strcpy(sam_results[i],buffer);
			i++;
		}
	}
}
void adjust_start_end(int start, int end, int* return_arr, int max_sam_length){
	char* s;
	char* buffer_copy = (char*)malloc((max_sam_length+1)*sizeof(char));
	strcpy(buffer_copy,sam_results[start]);
	s = strtok(buffer_copy,"\t");
	s = strtok(NULL, "\t");
	int decimal=0;
	sscanf(s, "%d", &decimal);
	decimal=dec2bin(decimal);
	if (decimal==0){ //second in pair
		return_arr[0]=start-1;
	}else{
		return_arr[0]=start;
	}
	strcpy(buffer_copy,sam_results[end]);
	s = strtok(buffer_copy,"\t");
	s = strtok(NULL, "\t");
	decimal=0;
	sscanf(s, "%d", &decimal);
	decimal=dec2bin(decimal);
	if (decimal==1){
		return_arr[1]=end+1;
	}else{
		return_arr[1]=end;
	}
	free(buffer_copy);
}
int main(int argc, char **argv){
	struct timespec tstart={0,0}, tend={0.0};
	Options opt;
	opt.remove_identical=0;
	opt.paired=0;
	opt.error=0.005;
	opt.coverage=50;
	opt.fasta_format=0;
	opt.freq=0.01;
	opt.llr=0;
	opt.min_strains=500;
	opt.max_strains=10000;
	opt.number_of_cores=1;
	memset(opt.print_counts,'\0',1000);
	memset(opt.MSA_reference,'\0',1000);
	parse_options(argc, argv, &opt);
	char* buffer = (char*)malloc(FASTA_MAXLINE*sizeof(char));
	memset(buffer,'\0',FASTA_MAXLINE);
	if (opt.paired==1 && opt.fasta_format==1){
		sprintf(buffer,"bowtie2 --all -f -x %s -1 %s -2 %s -S %s",opt.bowtie_reference_db,opt.forward_end_file,opt.reverse_end_file,opt.sam);
		system(buffer);
	}else if (opt.paired==0 && opt.fasta_format==1){
		sprintf(buffer,"bowtie2 --all -f -x %s -U %s -S %s",opt.bowtie_reference_db,opt.single_end_file,opt.sam);
		system(buffer);
	}else if (opt.paired==1 && opt.fasta_format==0){
		sprintf(buffer,"bowtie2 --all -x %s -1 %s -2 %s -S %s",opt.bowtie_reference_db,opt.forward_end_file,opt.reverse_end_file,opt.sam);
		system(buffer);
	}else if (opt.paired==0 && opt.fasta_format==0){
		sprintf(buffer,"bowtie2 --all -x %s -U %s -S %s",opt.bowtie_reference_db,opt.single_end_file,opt.sam);
		system(buffer);
	}
	free(buffer);
	gzFile MSA_file = Z_NULL;
	if (( MSA_file = gzopen(opt.fasta,"r")) == Z_NULL ) fprintf(stderr, "File could not be opened.\n");
	int length_of_MSA=0;
	length_of_MSA=setMSALength(MSA_file);
	gzclose(MSA_file);
	printf("Length of MSA: %d\n",length_of_MSA);
	int i,j,k,l;
	int* problematic_sites = (int*)malloc(1000*sizeof(int));
	for(i=0; i<1000; i++){
		problematic_sites[i]=-1;
	}
	int number_of_problematic_sites=process_problematic_sites(problematic_sites);
	reference_index = (int*)malloc(30000*sizeof(int));
	for(i=0; i<30000; i++){
		reference_index[i]=0;
	}
	align_references(number_of_problematic_sites,problematic_sites,opt.MSA_reference);
	if (( MSA_file = gzopen(opt.fasta,"r")) == Z_NULL ) fprintf(stderr, "File could not be opened.\n");
	int *strain_info = (int*)malloc(2*sizeof(int));
	setNumStrains(MSA_file,strain_info);
	gzclose(MSA_file);
	//strain_info[0]=3503116;
	//strain_info[1]=17;
	int number_of_strains = strain_info[0];
	int max_name_length = strain_info[1];
	free(strain_info);
	printf("Number of strains: %d\n",number_of_strains);
	printf("Maxname length: %d\n",max_name_length);
	char** MSA = (char **)malloc(number_of_strains*sizeof(char *));
	for(i=0; i<number_of_strains; i++){
		MSA[i] = (char *)malloc((length_of_MSA+1)*sizeof(char));
	}
	char** names_of_strains = (char**)malloc(number_of_strains*sizeof(char *));
	for(i=0; i<number_of_strains; i++){
		names_of_strains[i] = (char *)malloc((max_name_length+1)*sizeof(char));
	}
	/*int* reference = (int *)malloc(length_of_MSA*sizeof(int));
	for(i=0; i<length_of_MSA; i++){
		reference[i] = -1;
	}
	FILE* reference_file;
	if (( reference_file = fopen(opt.reference,"r")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	readReferencePositionsFile(reference_file,reference);
	fclose(reference_file);*/
	printf("Reading in MSA...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	if (( MSA_file = gzopen(opt.fasta,"r")) == Z_NULL ) fprintf(stderr, "File could not be opened.\n");
	int ref_index = readInMSA(MSA_file,MSA,names_of_strains,length_of_MSA);
	gzclose(MSA_file);
	//for(i=0; i<length_of_MSA; i++){
	//	free(allele_max[i]);
	//}
//	free(allele_max);
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
	/*int* identical = (int *)malloc(number_of_strains*sizeof(int));
	int number_of_identical_strains=number_of_strains;
	if (opt.remove_identical==0){
		for(i=0; i<number_of_strains; i++){
			identical[i]=i;
		}
	}*/
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
	printf("Number of variant sites: %d\n",number_of_variant_sites);
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
	int *max_sam_length = (int*)malloc(2*sizeof(int));
	max_sam_length[0]=0;
	max_sam_length[1]=0;
	if ( opt.paired==1 ){
		number_of_strains_remaining=calculateAlleleFreq_paired(sam_file,allele_frequency,length_of_MSA,MSA,number_of_strains,names_of_strains,opt.freq,max_name_length,tstart,tend,number_of_variant_sites,variant_sites,opt.coverage,opt.min_strains,opt.max_strains,opt.print_counts,max_sam_length);
	}else{
		number_of_strains_remaining=calculateAlleleFreq(sam_file,allele_frequency,length_of_MSA,MSA,number_of_strains,names_of_strains,opt.freq,max_name_length,tstart,tend,number_of_variant_sites,variant_sites,opt.coverage,opt.print_counts);
	}
	fclose(sam_file);
	for(i=0; i<length_of_MSA; i++){
		free(allele_frequency[i]);
	}
	free(allele_frequency);
	//free(identical);
	int* strains_kept = (int *)malloc(number_of_strains_remaining*sizeof(int));
	int placement=0;
	for(i=0; i<number_of_strains; i++){
		if (names_of_strains[i][0] != '\0' ){
			strains_kept[placement]=i;
			placement++;
		}
	}
	resize_names_of_strains = (char**)malloc(number_of_strains_remaining*sizeof(char*));
	resize_MSA = (char**)malloc(number_of_strains_remaining*sizeof(char*));
	for(i=0; i<number_of_strains_remaining; i++){
		resize_names_of_strains[i] = (char*)malloc((max_name_length+1)*sizeof(char));
		resize_MSA[i] = (char*)malloc((length_of_MSA+1)*sizeof(char));
		for(j=0; j<length_of_MSA+1; j++){
			resize_MSA[i][j]='\0';
		}
		for(j=0; j<max_name_length+1; j++){
			resize_names_of_strains[i][j]='\0';
		}
	}

	reallocate_memory(number_of_strains_remaining,strains_kept,MSA,names_of_strains,max_name_length,number_of_strains);
	printf("Creating mismatch matrix...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	FILE *outfile;
	if (( outfile = fopen(opt.outfile,"w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	if ((sam_file = fopen(opt.sam,"r")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	if (opt.paired==1){
		sam_results = (char**)malloc(max_sam_length[1]*sizeof(char*));
		for(i=0; i<max_sam_length[1]; i++){
			sam_results[i]=(char*)malloc((max_sam_length[0]+1)*sizeof(char));
			for(j=0; j<max_sam_length[0]+1; j++){
				sam_results[i][j]='\0';
			}
		}
		readInSamFile(sam_file);
		fclose(sam_file);
		fprintf(outfile,"qName\tblockSizes");
		for(i=0; i<number_of_strains_remaining; i++){
			fprintf(outfile,"\t%s",resize_names_of_strains[i]);
		}
		fprintf(outfile,"\n");
		pthread_t threads[opt.number_of_cores];//array of our threads
		thread_struct tstr[opt.number_of_cores];//array of stuct that contains input and output for each thread
		for(i=0; i<opt.number_of_cores; i++){
			tstr[i].start = 0;
			tstr[i].end = 0;
			tstr[i].str = malloc(sizeof(struct resultsStruct));
			tstr[i].number_of_strains_remaining=number_of_strains_remaining;
			tstr[i].number_of_strains=number_of_strains;
			tstr[i].length_of_MSA=length_of_MSA;
		}
		int divideFile, start, end;
		divideFile = max_sam_length[1]/opt.number_of_cores;
		j=0;
		for(i=0; i<opt.number_of_cores; i++){
			start=j;
			end=j+divideFile;
			if ( i==opt.number_of_cores-1){
				end=max_sam_length[1]-1;
			}
			int* adjust_ends = (int*)malloc(2*sizeof(int));
			if ( tstr[i].start == 1 ){
				start++;
			}
			adjust_ends[0]=0;
			adjust_ends[1]=0;
			adjust_start_end(start,end,adjust_ends,max_sam_length[0]);
			if ( start != adjust_ends[0] ){
				tstr[i-1].end = tstr[i-1].end-1;
			}
			if (end != adjust_ends[1] ){
				tstr[i+1].start = tstr[i+1].start+1;
			}
			if ( i==opt.number_of_cores-1){
				adjust_ends[1]=max_sam_length[1];
			}
			tstr[i].start = adjust_ends[0];
			tstr[i].end =adjust_ends[1];
			tstr[i].thread_number=i;
			tstr[i].max_sam_length=max_sam_length[0];
			j=j+divideFile;
			tstr[i].str->mismatch =(char**) malloc((end-start)*(sizeof(char *)));
			for(k=0; k<end-start; k++){
				tstr[i].str->mismatch[k]=(char*)malloc((max_sam_length[0]+300000)*sizeof(char));
				for(l=0; l<max_sam_length[0]+300000; l++){
					tstr[i].str->mismatch[k][l]='\0';
				}
			}
			free(adjust_ends);
		}
		for(i=0; i<opt.number_of_cores; i++){
			printf("thread: %d start: %d end: %d\n",tstr[i].thread_number,tstr[i].start,tstr[i].end);
		}
		for(i=0; i<opt.number_of_cores; i++){
			pthread_create(&threads[i], NULL, writeMismatchMatrix_paired, &tstr[i]);
		}
		for(i=0; i<opt.number_of_cores; i++){
			pthread_join(threads[i], NULL);
		}
		for(i=0; i<opt.number_of_cores; i++){
			for ( j=0; j<(tstr[i].end-tstr[i].start); j++){
				if ( tstr[i].str->mismatch[j][0]=='\0' ){
					break;
				}
				fprintf(outfile,"%s\n",tstr[i].str->mismatch[j]);
			}
		}
			//writeMismatchMatrix_paired(outfile,sam_file,resize_MSA,length_of_MSA,number_of_strains,number_of_strains_remaining,resize_names_of_strains,reference_index);
	}else{
		writeMismatchMatrix(outfile,sam_file,MSA,strains_kept,length_of_MSA,number_of_strains,number_of_strains_remaining,names_of_strains);
	}
	fclose(outfile);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fseconds\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	//free(reference);
	//free(variant_sites);
	//for(i=0; i<number_of_strains; i++){
		//free(names_of_strains[i]);
	//	free(MSA[i]);
	//}
	//free(names_of_strains);
	//free(MSA);
	if ( number_of_strains_remaining > 0 ){
		free(strains_kept);
	}else{
		printf("No strains remaining exiting...\n");
		exit(1);
	}
	buffer = (char*)malloc(FASTA_MAXLINE*sizeof(char));
	memset(buffer,'\0',FASTA_MAXLINE);
	if (opt.llr==1){
		sprintf(buffer,"Rscript EM_C_LLR_updated.R -i %s -f %lf -e %lf -l -s -v %s -r %s -b %s",opt.outfile,opt.freq,opt.error,opt.variant,opt.fasta,opt.print_counts);
	}else{
		sprintf(buffer,"Rscript EM_C_LLR_updated.R -i %s -f %lf -e %lf -s -v %s -r %s -b %s",opt.outfile,opt.freq,opt.error,opt.variant,opt.fasta,opt.print_counts);
	}
	system(buffer);
	free(buffer);
}
