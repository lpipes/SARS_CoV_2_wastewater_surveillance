#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include <assert.h>
#include <math.h>
#include "hashmap.h"
#include "math.h"
#include "global.h"
#include "options.h"

int tip=0;
int comma=0;
int open=0;
double minVariance = 99999999999999999;
int minVarNode = -1;
HASHMAP(char, struct specmap) sm;

void printtree(node* tree, int root, int numspec){
	int i;
	for (i=0; i<2*numspec-1; i++){
		if (tree[i].up[0] != -1){
			printf("Node %i: up: (%i, %i) down: %i, nd: %i, ",i,tree[i].up[0],tree[i].up[1],tree[i].down,tree[i].nd);
		}else{
			printf("(Node %i): up: (%i, %i) down: %i, nd: %i, name: %s",i,tree[i].up[0],tree[i].up[1],tree[i].down,tree[i].nd,tree[i].name);
		}
		if (i != root){
			printf(" (bl: %lf)\n",tree[i].bl);
		}else{
			printf("\n");
		}
	}
}
void readNodeList(FILE *infile, char** nodeIDs){
	char buffer[FASTA_MAXLINE];
	char *s;
	int size;
	char nodename[MAX_NODENAME];
	int i;
	int index=0;
	while ( fgets(buffer,FASTA_MAXLINE,infile) != NULL){
		s = strtok(buffer,"\n");
		size = strlen(s);
		for(i=0; i<size; i++){
			nodename[i]=buffer[i];
		}
		nodename[size]='\0';
		strcpy(nodeIDs[index],nodename);
		index++;
	}
}
void linknodes(node* tree, int i,int j,int node) /*linking i down to node and j down to node*/
{
	tree[node].up[0]=j;
	tree[node].up[1]=i;
	tree[i].down=node;
	tree[j].down=node;
}

/*subfunction needed by ‘getclade*/
int specsearch(node* tree, FILE *treefile, char **nodeIDs,int numspec,/*struct specmap* sm,*/ double average){
	char ch;
	int i=0;
	char specname[MAX_NODENAME];
	ch = fgetc(treefile);
	if ((ch!=')')&&(ch!='(')&&(ch!=',')&&(ch!=' ')&&(ch!='\t')&&(ch!='\n')&&(ch!=EOF)){
		ungetc(ch, treefile);
		while ((ch=fgetc(treefile))!=':'&&(i<MAX_NODENAME)){
			if ( isalpha(ch) || isdigit(ch) || ch=='.' || ch=='_' || ch=='/' || ch=='-'){
				specname[i]=ch;	
				i++;
			}
		}
		specname[i]='\0';
		/*if (i==0){
			ungetc(ch, treefile);
			return 0;
		}*/
		//This is sooo slow... make faster
		if ( i !=0 ){
		int index = hashmap_get(&sm,specname);
		//int tmp1=0;
		//for(tmp1=0;tmp1<numspec;tmp1++){
		//	if ( !strcmp(nodeIDs[tmp1],specname) ){
		//		tip = tmp1+1;
		//		break;
		//	}
		//}
		tip = index + 1;
		strcpy(tree[tip+numspec-2].name,specname);
		tree[tip+numspec-2].up[0]=-1;
		tree[tip+numspec-2].up[1]=-1;
		fscanf(treefile,"%le",&tree[tip+numspec-2].bl);
		if ( tree[tip+numspec-2].bl == 0 ){
			tree[tip+numspec-2].bl = EPSILON/average;
			//tree[tip+numspec-2].bl = 1.0000000000000001e-09;
		}else{
			tree[tip+numspec-2].bl = tree[tip+numspec-2].bl/average;
			//tree[tip+numspec-2].bl = tree[tip+numspec-2].bl/29892;
		}
		}
		return 1;
	}else{
		ungetc(ch, treefile);
		return 0;
	}
}
/*subfunction needed by ‘getclade*/
int getnodenumb(FILE *treefile){
	char c;
	int i,j=0;
	fpos_t position;
	i=0;
	fgetpos(treefile, &position);
	do{
		c=fgetc(treefile);
		if (c==',')
			i++;
		if (c=='(')
			j=j-1;
		if (c==')')
			j++;
	}
	while ((j<0)&&(c!=EOF));
	fsetpos(treefile,&position);
	return (i+comma+1);
}
/*some old code for reading a Newick tree*/
int getclade(node* tree, FILE *treefile, char** nodeIDs, int numspec, /*struct specmap* sm,*/ double average){
	int n1, n2, n3;
	char ch;
	n1=-1;
	n2=-1;
	n3=-1;
	do{
		if (specsearch(tree,treefile,nodeIDs,numspec,/*&sm,*/average)==1){
			/*tip++;*/
			/*if ( open > 1 ){
				printf("found polytomy\n");
				int new_num = getnodenumb(treefile);
				printf("new_num parent %d\n",new_num);
				printf("child %d\n",tip+numspec-2);
				//if ( tree[new_num].up[0] == -1){
				//	new_num = tree[new_num].down;
				//}
				tree[tip+numspec-2].down = new_num;
				tree[tip+numspec-2].up[0] = -1;
				tree[tip+numspec-2].up[1] = -1;
				tree[new_num].up[0] = tip+numspec-2;
			}*/
			return tip+(numspec-1);
		}
		ch = fgetc(treefile);
		if (ch==','){
			comma++;
			//open++;
		}
		if (ch==')'){
			if ((ch=fgetc(treefile))!=':'){
				ungetc(ch,treefile);
			}else{
				do{
					ch=(fgetc(treefile));
				}while((ch=='\n')||(ch==' '));
				ungetc(ch,treefile);
				fscanf(treefile,"%le",&tree[n3-1].bl);
				if ( tree[n3-1].bl == 0 ){
					tree[n3-1].bl = EPSILON/average;
				}else{
					tree[n3-1].bl = tree[n3-1].bl/average;
				}
			}
			//if ((ch=fgetc(treefile))==';'){	
				return n3;
			//}else{
			//	ungetc(ch,treefile);
			//}
		}
		if (ch=='('){
			//open=0;
			n3=getnodenumb(treefile);
			n1=getclade(tree,treefile,nodeIDs,numspec,/*&sm,*/average);
			n2=getclade(tree,treefile,nodeIDs,numspec,/*&sm,*/average);
			if ( n1!=-1 && n2!=-1 && n3!=-1){linknodes(tree,n1-1,n2-1,n3-1);}
			if ( n1 ==-1 || n2==-1 || n3 ==-1 ){ ungetc(ch,treefile); }
		}
	}
	while (ch!=';');
}
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
void readInMSA( gzFile MSA_file, int** MSA, int length_of_MSA){
	char buffer [length_of_MSA+1];
	int seq_idx=0;
	int i=0;
	int index=-1;
	while( gzgets(MSA_file,buffer,length_of_MSA+1) != NULL ){
		if ( buffer[0] != '>'){
			int size = strlen(buffer);
			for(i=0; i<size; i++){
				if (buffer[i]=='A' || buffer[i]=='a'){
					MSA[index][seq_idx]=0;
				}else if ( buffer[i]=='C' || buffer[i]=='c'){
					MSA[index][seq_idx]=1;
				}else if (buffer[i]=='G' || buffer[i]=='g'){
					MSA[index][seq_idx]=2;
				}else if (buffer[i]=='T' || buffer[i]=='t'){
					MSA[index][seq_idx]=3;
				}else{
					MSA[index][seq_idx]=4;
				}
				seq_idx++;
			}
		} else {
			index++;
			seq_idx=0;
		}
	}
}
void readInMSA_for_partition(node* tree, gzFile MSA_file, int** MSA, int* leaf_nodes, int length_of_MSA, int max_name, int size_of_partition){
	char buffer [length_of_MSA+1];
	int i=0;
	int j=0;
	int index=-1;
	int choose = 0;
	char name[max_name];
	memset(name,'\0',max_name);
	while( gzgets(MSA_file,buffer,length_of_MSA+1) != NULL){
		if ( buffer[0] == '>' && j < size_of_partition){
			for(i=1; buffer[i]!='\n'; i++){
				name[i-1]=buffer[i];
			}
			name[i-1]='\0';
			if ( strcmp(name,tree[leaf_nodes[j]].name)==0 ){
				choose=1;
				j++;
				index++;
			}	
		}else if ( choose==1 ){
			int size = strlen(buffer);
			for(i=0; i<size; i++){
				if (buffer[i]=='A' || buffer[i]=='a'){
					MSA[index][i]=0;
				}else if ( buffer[i]=='C' || buffer[i]=='c'){
					MSA[index][i]=1;
				}else if (buffer[i]=='G' || buffer[i]=='g'){
					MSA[index][i]=2;
				}else if (buffer[i]=='T' || buffer[i]=='t'){
					MSA[index][i]=3;
				}else{
					MSA[index][i]=4;
				}
			}
			choose=0;
			if ( j == size_of_partition){ break; }
		}
	}
}
int readInImputed(char** nodeIDs, char** MSA, int* reference, int max_name, int length_of_MSA,int start, int end, int ref_length, Options opt, int numspec){
	char buffer[FASTA_MAXLINE];
	//char *buffer = NULL;
	//char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int i=0;
	int j=0;
	int k=0;
	int index=-1;
	char name[max_name];
	char sequence[length_of_MSA];
	int ref_index = -1;
	printf("creating hash...\n");
	HASHMAP(char, struct blob) hash;
	hashmap_init(&hash, hashmap_hash_string, strcmp);
	//gzFile MSA_file = Z_NULL;
	//if ((MSA_file=gzopen(opt.outfile,"r"))==Z_NULL ){ puts("Cannot open MSA file!"); exit(-1);}
	FILE* MSA_file;
	MSA_file = fopen(opt.outfile,"r");
	if (MSA_file == 0){
		perror("ERROR OPENING FILE!");
		exit(-1);
	}
	//if (NULL==(MSA_file=fopen(opt.outfile,"r"))){ puts("Cannot open file!"); exit(-1);}
	//while( ( read=getline(&buffer, &len, MSA_file)) != -1){
	while( fgets(buffer,FASTA_MAXLINE,MSA_file) != NULL){
		if ( buffer[0] == '>'){
			index++;
			for(i=1; buffer[i]!='\n'; i++){
				name[i-1]=buffer[i];
			}
			name[i-1]='\0';
			strcpy(nodeIDs[index],name);
			if (strcmp(name,"EPI_ISL_402124")==0){
				ref_index=index;
			}
		}else{
			int size = strlen(buffer);
			j=0;
			for(i=start; i<end; i++){
				if ( i==reference[j] ){
					sequence[j]=buffer[i];
					j++;
				}
			}
			sequence[j]='\0';
			strcpy(MSA[index],sequence);
			//hashmap_put(&hash,MSA[index],nodeIDs[index]);
			struct blob *c;
			c=hashmap_get(&hash, MSA[index]);
			if (c==NULL){
				struct blob *b;
				b = malloc(sizeof(*b));
				b->key = MSA[index];
				b->data = nodeIDs[index];
				b->data_len = strlen(nodeIDs[index]);
				hashmap_put(&hash, b->key, b);
			}else{
				char x = ':';
				for(k=0; k<MAX_STRAINS*max_name; k++){
					if ( c->data[k] == '\0' ){
						break;
					}		
				}
				c->data[k] = ':';
				strcat(c->data,nodeIDs[index]);
				hashmap_put(&hash, c->key, c);
			}
		}
	}
	fclose(MSA_file);
	//hashmap_put(&hash, MSA[ref_index], nodeIDs[index]);
	const char *key;
	const char *value;
	FILE* remove_ident_file;
	int remaining_strains=0;
	struct blob *d;
	printf("printing final out file\n");
	if (( remove_ident_file = fopen(opt.out_MSA,"w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	hashmap_foreach(key, d, &hash){
		fprintf(remove_ident_file,">%s\n",d->data);
		fprintf(remove_ident_file,"%s\n",key);
		remaining_strains++;
	}
	fclose(remove_ident_file);
	//printf("remaining strains: %d\n",remaining_strains);
	//strains_remaining=malloc(remaining_strains*sizeof(int));
	//for(i=0; i<remaining_strains; i++){
	//	strains_remaining[i]=-1;
	//}
	//j=0;
	//hashmap_foreach(key,value, &hash){
	//	for(i=0; i<numspec; i++){
	//		if (strcmp(value,nodeIDs[i])==0){
	//			strains_remaining[j]=i;
	//			j++;
	//		}
	//	}
	//}
	return remaining_strains;
}
int readInNodes(gzFile MSA_file, char** nodeIDs, int length_of_MSA, /*struct blob *sm,*/ int* reference, int* start_end){
	char buffer[length_of_MSA+1];
	int i=0;
	int index=-1;
	int save=0;
	int ref_length=0;
	int iter=0;
	//HASHMAP(char, struct specmap) sm;
	//hashmap_init(&sm, hashmap_hash_string, strcmp);
	while( gzgets(MSA_file,buffer,length_of_MSA+1) != NULL ){
		if ( buffer[0] == '>'){
			index++;
			for(i=1; buffer[i]!='\n'; i++){
				nodeIDs[index][i-1]=buffer[i];
			}
			nodeIDs[index][i-1]='\0';
			if ( strcmp(nodeIDs[index],"EPI_ISL_402124")==0){
				save=1;
			}
			hashmap_put(&sm,nodeIDs[index],index);
		}else if ( save==1){
			int iter = 0;
			for(i=0; i<buffer[i]!='\n'; i++){
				if( buffer[i] == 'A' || buffer[i] == 'G' || buffer[i] == 'C' || buffer[i] == 'T' || buffer[i] =='a' || buffer[i] =='c' || buffer[i]=='t' || buffer[i]=='g'){
					start_end[0] = i;
					break;
				}
			}
			//printf("start is %d\n",start_end[0]);
			for(i=length_of_MSA-1; i>=0; i--){
				if ( buffer[i]=='G' || buffer[i] == 'C' || buffer[i] == 'T' || buffer[i]=='g' || buffer[i]=='c' || buffer[i]=='t'){
					start_end[1]=i;
					break;
				}
			}
			//printf("end is %d\n",start_end[1]);
			for(i=start_end[0]; i<start_end[1]+1; i++){
				if ( buffer[i]!= '-' ){
					iter++;
				}
			}
			ref_length = iter;
			iter=0;
			for(i=start_end[0]; i<start_end[1]+1; i++){
				if ( buffer[i]!= '-' ){
					reference[iter]=i;
					iter++;
				}
			}
			save=0;
		}
	}
	return ref_length;
}
void maketransitionmatrixnc(int n, double t, double LRVECnc[4][4], double RRVECnc[4][4], double RRVALnc[4], double PMATnc[2][4][5]){
	int i, j, k;
	double EXPOS[4];
	for (k=0; k<4; k++){
		EXPOS[k] = exp(t*RRVALnc[k]);
		if ( EXPOS[k] < 0.000000 ){
			printf("EXPOS[k]=%e\n",EXPOS[k]);
		}
		assert(EXPOS[k]>=0.000000);
	}

	for (i=0; i<4; i++)
  {
    for (j=0; j<4; j++)
    {
      PMATnc[n][i][j] = 0.0;
      for (k=0; k<4; k++){
		  PMATnc[n][i][j] =  PMATnc[n][i][j] + RRVECnc[k][j]*LRVECnc[i][k]*EXPOS[k];
	  }
	}
  }
}
void makeconnc(node* tree, int node, int root,int numspec, int numbase,double lambda,int** seq, double* UFCnc, double LRVECnc[4][4], double RRVECnc[4][4], double RRVALnc[4], double PMATnc[2][4][5],  int number_of_leaves)

{
  int i, j, child, seqn, site;
  double L, max;
	int l=0;
  child = tree[node].up[0];


  if (tree[child].up[0]==-1){
    maketransitionmatrixnc(0, lambda*tree[child].bl, LRVECnc, RRVECnc, RRVALnc, PMATnc);
    //seqn=child-numspec+1;
    //int this_node=-1;
    //int this_node = hashmap_get(&node_list_coding,child);
    /*for(l=0; l<number_of_leaves; l++){
	if ( child == node_list_coding[l]){
		this_node = l;
		break;
	}
    }*/
    int this_node = tree[child].pos;
    seqn=this_node;
    for (site=0; site<numbase; site++){
      for (i=0; i<4; i++)
        tree[node].likenc[site][i] = PMATnc[0][i][seq[seqn][site]];
       if (site==0 && node ==root) printf("node 10 (b=%i): %lf %lf %lf %lf\n",seq[seqn][site],PMATnc[0][0][seq[seqn][site]],PMATnc[0][1][seq[seqn][site]],PMATnc[0][2][seq[seqn][site]],PMATnc[0][3][seq[seqn][site]]);
    }
  }
  else {
    makeconnc(tree,child,root,numspec,numbase,lambda,seq,UFCnc,LRVECnc, RRVECnc, RRVALnc, PMATnc, number_of_leaves);
    maketransitionmatrixnc(0, lambda*tree[child].bl, LRVECnc, RRVECnc, RRVALnc, PMATnc);
    for (site=0; site<numbase; site++)
    {
      for (i=0; i<4; i++){
        tree[node].likenc[site][i]=0.0;
        for (j=0; j<4; j++)
          tree[node].likenc[site][i] += PMATnc[0][i][j]*tree[child].likenc[site][j];
      }
    }
  }
  child = tree[node].up[1];
  if (tree[child].up[1]==-1){
    maketransitionmatrixnc(0,lambda*tree[child].bl,LRVECnc, RRVECnc, RRVALnc, PMATnc);
    //seqn=child-numspec+1;
	//int this_node = -1;
	//int this_node = hashmap_get(&node_list_coding,child);
	/*for(l=0; l<number_of_leaves; l++){
		if ( child == node_list_coding[l] ){
			this_node = l;
			break;
		}
	}*/
	int this_node = tree[child].pos;
	seqn = this_node;
    for (site=0; site<numbase; site++)
      for (i=0; i<4; i++)
        tree[node].likenc[site][i] = tree[node].likenc[site][i]*PMATnc[0][i][seq[seqn][site]];
  }
  else {
    makeconnc(tree,child, root,numspec,numbase,lambda,seq,UFCnc,LRVECnc, RRVECnc, RRVALnc, PMATnc, number_of_leaves);
    maketransitionmatrixnc(0,lambda*tree[child].bl,LRVECnc, RRVECnc, RRVALnc, PMATnc);
    for (site=0; site<numbase; site++)
    {
      max=0.0;
      for (i=0; i<4; i++){
        L=0.0;
        for (j=0; j<4; j++)
          L += PMATnc[0][i][j]*tree[child].likenc[site][j];
        //underflow control
        if ((tree[node].likenc[site][i] = tree[node].likenc[site][i]*L)>max)
          max = tree[node].likenc[site][i];
      }
      //might be worth with some code here dealing with the case of max=0+epsilon

	if (max<0.00000000001){
		printf("Warning, max = %lf\n",max);
	}

      for (i=0; i<4; i++)
        tree[node].likenc[site][i]=tree[node].likenc[site][i]/max;
      UFCnc[site] = UFCnc[site] + log(max);
    }
  }
}
void inittransitionmatrixnc(double pi[4], double par[10], double LRVECnc[4][4], double RRVECnc[4][4], double RRVALnc[4])

{
  int i, j;
  double sum, piT, RIVAL[4], RIVEC[4][4],  A[4][4], workspace[8];
  for (i=0; i<8; i++){
	  workspace[i]=0;
  }
  A[0][1]=pi[1]*par[4];
  A[0][2]=pi[2]*par[5];
  A[0][3]=pi[3]*par[6];
  A[1][0]=pi[0]*par[4];
  A[1][2]=pi[2]*par[7];
  A[1][3]=pi[3]*par[8];
  A[2][0]=pi[0]*par[5];
  A[2][1]=pi[1]*par[7];
  A[2][3]=pi[3]; //unscaled rate of GT = 1.0
  A[3][0]=pi[0]*par[6];
  A[3][1]=pi[1]*par[8];
  A[3][2]=pi[2]; //unscaled rate of GT = 1.0

  for (i=0; i<4; i++)
  {
    A[i][i]=0.0;
    sum=0.0;
    for (j=0; j<4; j++)
      sum = sum + A[i][j];
    A[i][i] = -sum;
  }
  if (eigen(1, A[0], 4, RRVALnc, RIVAL, RRVECnc[0], RIVEC[0], workspace) != 0)
  {
    printf("Transitions matrix did not converge or contained non-real values!\n");
    exit(-1);
  }

  for (i=0; i<4; i++)
    for (j=0; j<4; j++){
		LRVECnc[i][j] = RRVECnc[i][j];
	}
  if (matinv(RRVECnc[0],4, 4, workspace) != 0)
    printf("Could not invert matrix!\nResults may not be reliable!\n");
}
double getlike_gamma(double parameters[], node* tree, int numspec, int numbase, int root, int** MSA, double LRVECnc[4][4], double RRVECnc[4][4], double RRVALnc[4],double PMATnc[2][4][5], int number_of_leaves){
  /*GTR + gamma model
    par[1]: inverse function for piA;
    par[2]: inverse function for piC;
    par[3]: inverse function for piG;
    par[4]: AC;
    par[5]: AG;
    par[6]: AT;
    par[7]: CG;
    par[8]: CT;
    par[9]: alpha;
    unscaled rate of G<->T defined to be 1.0;
    */

double stand, L, loclike, **locloglike, max, pi[4], gampar[2], d, like = 0.0;
  int i, j, k;
  stand = 1.0+parameters[1]+parameters[2]+parameters[3];
  pi[0]=parameters[1]/stand;
  pi[1]=parameters[2]/stand;
  pi[2]=parameters[3]/stand;
  pi[3]=1.0-pi[0]-pi[1]-pi[2];
  gampar[0]=gampar[1]=parameters[9]; //We are setting alpha=beta to keep a constant mean to avoid identifiability issues.  This is not the same as a standard gammma.
  double* UFCnc = malloc((numbase)*(sizeof(double)));
  double* statevector = malloc(NUMCAT*(sizeof(double)));
  locloglike = malloc(numbase*(sizeof(double *)));
  for (i=0; i<numbase; i++){
	  locloglike[i] = malloc(NUMCAT*(sizeof(double)));
  } 
	definegammaquantiles(NUMCAT, gampar, statevector);
  //FIX THIS IF USING MULTIPLE CATEGORIES!!!!
  statevector[0]=1.0;
  inittransitionmatrixnc(pi, parameters, LRVECnc, RRVECnc, RRVALnc);
  for (j=0; j<NUMCAT; j++){
    for (i=0; i<numbase; i++)
      UFCnc[i]=0.0;
    makeconnc(tree,root,root,numspec,numbase,statevector[j],MSA,UFCnc,LRVECnc, RRVECnc, RRVALnc, PMATnc, number_of_leaves);
    for (i=0; i<numbase; i++){
      L=0.0;
      for (k=0;k<4;k++)
        L += tree[root].likenc[i][k]*pi[k];
      if (L>0.0) locloglike[i][j] = log(L) + UFCnc[i];
    }
  }
  for (i=0; i<numbase; i++){
    loclike=0.0;
    max = -100000000000.0;
    for (j=0; j<NUMCAT; j++){
      if (locloglike[i][j]>max) //underflow protection
        max=locloglike[i][j];
    }
    for (j=0; j<NUMCAT; j++){
      d=locloglike[i][j]-max;
      if (d>-100)
        loclike += exp(d);
    }
    like = like + log(loclike) + max;
  }
  free(statevector);
  for (i=0; i<numbase; i++)
    free(locloglike[i]);
  free(locloglike);
  free(UFCnc);
  return -like + (double)numbase*log((double)NUMCAT);
}
void clearGlobals(double parameters[10], double LRVECnc[4][4], double RRVECnc[4][4], double RRVALnc[4], double PMATnc[2][4][5]){
	int i,j, k;
	for(i=4;i<4;i++){
		RRVALnc[i]=0;
		for(i=4;i<4;i++){
			LRVECnc[i][j]=0;
			RRVECnc[i][j]=0;
		}
	}
	for (i=0;i<4;i++){
		PMATnc[0][i][4]=1.0;//missing data
		PMATnc[1][i][4]=1.0;//missing data
	}
	parameters[0]=0.0;
	for(i=1;i<10;i++){
		parameters[i]=1.0;
	}
}
void makeposterior_nc(node* tree, int node, int numbase, int numspec, int** seq, double** templike_nc, double LRVECnc[4][4], double RRVECnc[4][4], double RRVALnc[4], double PMATnc[2][4][5], int number_of_leaves,int parent)
{
  int i,j, s, /*parent,*/ otherb, child1, child2, b;
  double bl, max;
	int l=0;

  child1 = tree[node].up[0];
  child2 = tree[node].up[1];
  //parent = tree[node].down;
  bl = tree[node].bl;
  maketransitionmatrixnc(0, bl, LRVECnc, RRVECnc, RRVALnc, PMATnc);
  if ((otherb = tree[parent].up[0])==node)
    otherb = tree[parent].up[1];
  maketransitionmatrixnc(1, tree[otherb].bl, LRVECnc, RRVECnc, RRVALnc, PMATnc);
  for (s=0; s<numbase; s++){
    if (tree[otherb].up[0]>-1){
      for (i=0; i<4; i++){
        templike_nc[s][i]=0;
        for (j=0; j<4; j++)
          templike_nc[s][i] += tree[otherb].likenc[s][j]*PMATnc[1][i][j];
        templike_nc[s][i]=templike_nc[s][i]*tree[parent].posteriornc[s][i];
      }
    }
    else{
      //b=seq[otherb-numspec+1][s];
      //int this_node=-1;
      //int this_node = hashmap_get(&leaf_nodes_coding,otherb);
	/*for(l=0; l<number_of_leaves; l++){
		if ( otherb == leaf_nodes_coding[l] ){
			this_node = l;
			break;	
		}	
	}*/
	   int this_node = tree[otherb].pos;
      b=seq[this_node][s];
      for (i=0; i<4; i++)
        templike_nc[s][i] = PMATnc[1][i][b]*tree[parent].posteriornc[s][i];
    }
    for (i=0; i<4; i++){
      tree[node].posteriornc[s][i]=0.0;
      max=0.0;
      for (j=0; j<4; j++){
        if ((tree[node].posteriornc[s][i] = tree[node].posteriornc[s][i] + PMATnc[0][i][j]*templike_nc[s][j])>max)
          max=tree[node].posteriornc[s][i];//more underflow protection
      }
    }
    for (i=0; i<4; i++)
      tree[node].posteriornc[s][i]=tree[node].posteriornc[s][i]/max;
  }
  if (tree[child1].up[0]>-1)
    makeposterior_nc(tree,child1,numbase,numspec,seq,templike_nc,LRVECnc,RRVECnc,RRVALnc,PMATnc,number_of_leaves,node);
  if (tree[child2].up[0]>-1)
    makeposterior_nc(tree,child2,numbase,numspec,seq,templike_nc,LRVECnc,RRVECnc,RRVALnc,PMATnc,number_of_leaves,node);
}
//THIS FUNCTION IS IMPLEMENTED WITHOUT GAMMA CATEGORIES
void getposterior_nc(double parameters[10], node* tree,  int numspec, int numbase, int root, int** MSA, double LRVECnc[4][4], double RRVECnc[4][4], double RRVALnc[4], double PMATnc[2][4][5], int number_of_children, int* node_list, int number_of_leaves){
	int i, j, s, k, parent, b, notdonebefore;
	int l=0;
	double p, sum, pi[4], stand, **templike;
	//if ( tree[root].nd == 1){
	//	root = tree[root].down;
	//}
	getlike_gamma(parameters,tree,numspec,numbase,root,MSA,LRVECnc,RRVECnc,RRVALnc,PMATnc,number_of_leaves); //need to call likelihood again
	double ** templike_nc = malloc(numbase*(sizeof(double *)));
	for (i=0; i<numbase; i++){
		templike_nc[i]=malloc(4*(sizeof(double)));
  	}
	stand = 1.0+parameters[1]+parameters[2]+parameters[3];
	pi[0]=parameters[1]/stand;
	pi[1]=parameters[2]/stand;
	pi[2]=parameters[3]/stand;
	pi[3]=1.0-pi[0]-pi[1]-pi[2];
	for (s=0; s<numbase; s++){
		for (i=0; i<4; i++){
			tree[root].posteriornc[s][i] = 1.0;
		}
	}
	if (tree[tree[root].up[0]].up[0]>-1) makeposterior_nc(tree,tree[root].up[0],numbase,numspec,MSA,templike_nc,LRVECnc,RRVECnc,RRVALnc,PMATnc,number_of_leaves,root);
	if (tree[tree[root].up[1]].up[0]>-1) makeposterior_nc(tree,tree[root].up[1],numbase,numspec,MSA,templike_nc,LRVECnc,RRVECnc,RRVALnc,PMATnc,number_of_leaves,root);
	for (j=0; j<number_of_children; j++){
		if (tree[node_list[j]].up[0]>-1){
			for (s=0; s<numbase; s++){
				sum = 0.0;
				for (i=0; i<4; i++){
					sum = sum + (tree[node_list[j]].posteriornc[s][i]=tree[node_list[j]].likenc[s][i]*tree[node_list[j]].posteriornc[s][i]*pi[i]);
				}
				for (i=0; i<4; i++){
					tree[node_list[j]].posteriornc[s][i] = tree[node_list[j]].posteriornc[s][i]/sum;
				}
			}
		}else {
			for (s=0; s<numbase; s++) {
				notdonebefore=1;
				//int this_node=-1;
				//int this_node = hashmap_get(&leaf_nodes_coding,node_list[j]);
				/*for(l=0; l<number_of_leaves; l++){
					if ( node_list[j] == leaf_nodes_coding[l] ){
						this_node = l;
						break;
					}
				}*/
				//b=MSA[node_list[j]-numspec+1][s];
				int this_node = tree[node_list[j]].pos;
				b=MSA[this_node][s];
				if (b==4){
					if (notdonebefore==1) {
						maketransitionmatrixnc(0, tree[node_list[j]].bl,LRVECnc,RRVECnc,RRVALnc,PMATnc);
						notdonebefore=0;
					}
					parent=tree[node_list[j]].down;
					sum = 0.0;
					for (i=0; i<4; i++){
						tree[node_list[j]].posteriornc[s][i]=0.0;
						for (k=0; k<4; k++){
							tree[node_list[j]].posteriornc[s][i] += tree[parent].posteriornc[s][k]*PMATnc[0][i][k];
						}
						sum = sum + (tree[node_list[j]].posteriornc[s][i]=tree[node_list[j]].posteriornc[s][i]*pi[i]);
					}
					for (i=0; i<4; i++){
						tree[node_list[j]].posteriornc[s][i] = tree[node_list[j]].posteriornc[s][i]/sum;
					}
				}else{
					for (i=0; i<4; i++){
						if (i==b){
							tree[node_list[j]].posteriornc[s][i]=1.0;
						}else{
							tree[node_list[j]].posteriornc[s][i]=0.0;
          					}
        				}
      				}
			}
		}
	}
	for (i=0; i<numbase; i++){
		free(templike_nc[i]);
	}
	free(templike_nc);
}
void impute(node* tree, int node, int** MSA, int length_of_MSA, int numspec, int* leaf_nodes_coding, int number_of_leaves){
	int i,j, k;
	double minimum;
	int index;
	int this_node = -1;
	for(i=0; i<number_of_leaves; i++){
		if ( node == leaf_nodes_coding[i] ){
			this_node = i;
			break;
		}
	}
		for(j=0; j<length_of_MSA; j++){
			if (MSA[this_node][j]==4){
				minimum=1.0-tree[node].posteriornc[j][0];
				index=0;
				for(k=0; k<4; k++){
					if ( minimum > 1.0-tree[node].posteriornc[j][k]){
						minimum=1.0-tree[node].posteriornc[j][k];
						index=k;
					}
				}	
				MSA[this_node][j]=index;
			}
		}
}
void print_imputed_MSA(node* tree, int** MSA, int length_of_MSA, FILE* outfile, int* node_list_coding, int number_of_leaves){
	int i,j,k;
	for(i=0; i<number_of_leaves; i++){
		fprintf(outfile,">%s\n",tree[node_list_coding[i]].name);
		for(j=0; j<length_of_MSA; j++){
			char base;
			if (MSA[i][j]==0){
				base='A';
			}else if (MSA[i][j]==1){
				base='C';
			}else if (MSA[i][j]==2){
				base='G';
			}else if (MSA[i][j]==3){
				base='T';
			}else{
				printf(" UH OH\n");
			}
			fprintf(outfile,"%c",base);
		}
		fprintf(outfile,"\n");
	}
}
void calculate_max_allele(int** MSA, int number_of_strains, int length_of_MSA, int** allele_frequency, int* imputation){
	int i,j;
	for(i=0; i<number_of_strains; i++){
		for(j=0; j<length_of_MSA; j++){
			if ( MSA[i][j] == 0 ){
				allele_frequency[j][0]++;
			}else if ( MSA[i][j] == 1){
				allele_frequency[j][1]++;
			}else if (MSA[i][j] == 2){
				allele_frequency[j][2]++;
			}else if (MSA[i][j] == 3){
				allele_frequency[j][3]++;
			}
		}
	}
	for(i=0; i<length_of_MSA; i++){
		int max=0;
		int max_index=0;
		for(j=0; j<4; j++){
			if(allele_frequency[i][j] > max){
				max = allele_frequency[i][j];
				max_index = j;
			}
		}
		imputation[i]=max_index;
	}
}
void print_max_allele(FILE* outfile, int** MSA, int number_of_strains, int length_of_MSA, int* imputation, char** nodeIDs){
	int i,j;
	for(i=0; i<number_of_strains; i++){
		fprintf(outfile,">%s\n",nodeIDs[i]);
		for(j=0; j<length_of_MSA; j++){
			if ( MSA[i][j] == 4 ){
				if ( imputation[j] == 0 ){
					fprintf(outfile,"A");
				}else if (imputation[j] == 1){
					fprintf(outfile,"C");
				}else if(imputation[j] == 2){
					fprintf(outfile,"G");
				}else if (imputation[j] ==3){
					fprintf(outfile,"T");
				}
			}else{
				if (MSA[i][j] == 0 ){
					fprintf(outfile,"A");
				}else if (MSA[i][j] == 1){
					fprintf(outfile,"C");
				}else if (MSA[i][j] == 2){
					fprintf(outfile,"G");
				}else if (MSA[i][j] == 3){
					fprintf(outfile,"T");
				}
			}
		}
		fprintf(outfile,"\n");
	}
}
int get_number_descendants(node *tree, int node){
	if (tree[node].up[0]==-1) return (tree[node].nd=1);
	else return (tree[node].nd=(get_number_descendants(tree,tree[node].up[0])+get_number_descendants(tree,tree[node].up[1])));
}
void assignDepth(node* tree, int node0, int node1, int depth){
	if( node0 != -1 && node1 != -1){
		tree[node0].depth = depth;
		tree[node1].depth = depth;
		assignDepth(tree,tree[node0].up[0], tree[node0].up[1], depth+1);
		assignDepth(tree, tree[node1].up[0], tree[node1].up[1], depth+1);
	}
}
void find_descendants(node* tree, int node, int depth, int* descendants, int* node_list){
	int child0 = tree[node].up[0];
	int child1 = tree[node].up[1];
	int i;
	if ( child0 != -1 && child1 != -1 ){
		if ( tree[node].depth == depth ){
			//printf("node %d: %d\n",node,tree[node].nd);
			//printf("checking depth %d...\n",tree[node].depth);
			for(i=0; i<MAX_NODE_LIST; i++){
				if ( descendants[i]==-1){
					descendants[i]=tree[node].nd;
					node_list[i]=node;
					break;
				}
			}
			//checklimit(tree,node,tree[node].depth);
			//return tree[node].depth;
		}
		find_descendants(tree,child0,depth,descendants,node_list);
		find_descendants(tree,child1,depth,descendants,node_list);
	}
	//return tree[node].depth;
}
int find_max_depth(node* tree, int numspec){
	int i;
	int max=0;
	for(i=numspec; i<2*numspec-1; i++){
		if ( tree[i].depth > max ){
			max = tree[i].depth;
		}
	}
	return max;
}
void get_children(node* tree, int node, int* node_list){
	int child0 = tree[node].up[0];
	int child1 = tree[node].up[1];
	int i;
	if ( child0 == -1 && child1 == -1){
		for(i=0; i<MAX_NODE_LIST; i++){
			if (node_list[i]==node){
				break;
			}
			if (node_list[i]==-1){
				node_list[i]=node;
				break;
			}
		}
		return;
	}else{
		for(i=0; i<MAX_NODE_LIST; i++){
			if (node_list[i]==node){
				break;
			}
			if (node_list[i]==-1){
				node_list[i]=node;
				break;
			}
		}
		if (child0 != -1){
			get_children(tree,child0,node_list);
		}
		if (child1 != -1){
			get_children(tree,child1,node_list);
		}
	}
}
int number_of_nodes_at_depth(node* tree, int depth, int numspec){
	int i,j;
	int max=0;
	int max2=0;
	for(i=0; i<=depth; i++){
		max=0;
		for(j=0; j<2*numspec-1; j++){
			if ( tree[j].depth == i ){
				max++;
			}
		}
		if (max > max2){
			max2 = max;
		}
	}
	return max2;
}
void find_nodes_at_depth(node* tree, int node, int* node_list, int depth, int max){
	int child0 = tree[node].up[0];
	int child1 = tree[node].up[1];
	int i;
	if ( tree[node].depth == depth ){
		for(i=0; i<max; i++){
			if(node_list[i] == -1){
				node_list[i] = node;
				break;
			}
		}
	}
	if (child0 != -1 && child1 != -1){
		find_nodes_at_depth(tree,child0,node_list,depth,max);
		find_nodes_at_depth(tree,child1,node_list,depth,max);
	}
}
void find_in_node_list(node* tree, int node, int target_node, int* node_list, int iter){
	int child0=-1;
	int child1=-1;
	if ( node != -1 ){
		child0 = tree[node].up[0];
		child1 = tree[node].up[1];
	}
	if ( node == target_node){
		node_list[iter]=-1;
	}
	if ( child0 != -1 && child1 != -1 && node != -1){
		find_in_node_list(tree,child0,target_node,node_list,iter);
		find_in_node_list(tree,child1,target_node,node_list,iter);
	}
}
void findMinVariance(node* tree, int node, int size){
	int child0 = tree[node].up[0];
	int child1 = tree[node].up[1];
	int parent = tree[node].down;
	int i=0;
	if (node==-1){ return; }
	if (parent == -1){
		findMinVariance(tree,child0,size);
		findMinVariance(tree,child1,size);
		return;
	}
	if (child0 == -1 ){ return; }
	if (child1 == -1 ){ return; }
	int num_children0 = tree[child0].nd;
	int num_children1 = tree[child1].nd;
	int num_ancestors = size - num_children1 - num_children0;
	double mean = (double)(num_children0 + num_children1 + num_ancestors )/3;
	double variance = (double)((num_ancestors-mean)*(num_ancestors-mean) + (num_children0-mean)*(num_children0-mean) + (num_children1-mean)*(num_children1-mean))/3;
	if (minVariance > variance){
		minVariance = variance;
		minVarNode = node;
	}
	findMinVariance(tree,child0,size);
	findMinVariance(tree,child1,size);
	return;
}
int findStartNode(node* tree,int node){
	if ( node != -1){
		if (tree[node].down==-1){
			return node;
		}
		findStartNode(tree,tree[node].down);
	}
}
void recurseMinVar(node* tree, int node, int limit, int* node_list){
	if ( tree[node].nd > limit){
		minVariance = 99999999999999999;
		minVarNode = -1;
		findMinVariance(tree,node,tree[node].nd);
		if ( tree[ tree[ minVarNode ].down ].down != -1 ){
			
		if ( tree[tree[ tree[ minVarNode].down ].down].up[0] == tree[ minVarNode].down && tree[tree [minVarNode].down].up[0] != minVarNode){
			tree[tree[ tree[ minVarNode].down ].down].up[0] = tree[tree [minVarNode].down].up[0];
			tree [ tree[tree [minVarNode].down].up[0] ].down = tree[ tree[ minVarNode].down ].down;
		}else if ( tree[tree[ tree[ minVarNode].down ].down].up[1] == tree[ minVarNode].down && tree[tree [minVarNode].down].up[1] != minVarNode){
			tree[tree[ tree[ minVarNode].down ].down].up[1] = tree[tree [minVarNode].down].up[1];
			tree [ tree[tree [minVarNode].down].up[1] ].down = tree[ tree[ minVarNode].down ].down;
		}else if ( tree[tree[ tree[ minVarNode].down ].down].up[0] == tree[ minVarNode].down && tree[tree [minVarNode].down].up[1] != minVarNode){
			tree[tree[ tree[ minVarNode].down ].down].up[0] = tree[tree [minVarNode].down].up[1];
			tree[ tree[tree [minVarNode].down].up[1] ].down = tree[ tree[ minVarNode].down ].down;
		}else if ( tree[tree[ tree[ minVarNode].down ].down].up[1] == tree[ minVarNode].down && tree[tree [minVarNode].down].up[0] != minVarNode){
			tree[tree[ tree[ minVarNode].down ].down].up[1] = tree[tree [minVarNode].down].up[0];
			tree[ tree[tree [minVarNode].down].up[0] ].down = tree[ tree[ minVarNode].down ].down;
		}
		}else{
			if ( tree[ tree[minVarNode].down].up[1] == minVarNode ){
				tree[tree[tree[minVarNode].down].up[0]].down = -1;
				int i;
				int dont_add=0;
				for(i=0; i<MAX_NODE_LIST; i++){
					if ( node_list[i]==-1){
						break;
					}
				}
				int number_of_nodes = i;
				for(i=0; i<number_of_nodes; i++){
					if (node_list[i]==tree[tree[minVarNode].down].up[0]){
						dont_add=1;
					}
				}
				if ( dont_add==0 && tree[ tree[tree[minVarNode].down].up[0] ].nd < limit ){
					node_list[i] = tree[tree[minVarNode].down].up[0];
				}else{
					recurseMinVar(tree,tree[tree[minVarNode].down].up[0],limit,node_list);
				}
			}	
			if ( tree[ tree[minVarNode].down].up[0] == minVarNode ){
				tree[tree[tree[minVarNode].down].up[1]].down = -1;
				int i;
				for(i=0; i<MAX_NODE_LIST; i++){
					if ( node_list[i]==-1){
						break;
					}
				}
				int number_of_nodes = i;
				int dont_add=0;
				for(i=0; i<number_of_nodes; i++){
					if (node_list[i] == tree[tree[minVarNode].down].up[1]){
						dont_add=1;
					}
				}
				if ( dont_add==0 && tree[ tree[tree[minVarNode].down].up[1] ].nd < limit ){
					node_list[i] = tree[tree[minVarNode].down].up[1];
				}else{
					recurseMinVar(tree,tree[tree[minVarNode].down].up[1],limit,node_list);
				}
			}
		}
		int local_min_var = minVarNode;
		int start = findStartNode(tree,local_min_var);
		tree[start].nd = tree[start].nd - tree[local_min_var].nd;
		get_number_descendants(tree,start);
		tree[ tree[local_min_var].up[0] ].down = -1;
		tree[ tree[local_min_var].up[1] ].down = -1;
		int i;
		for(i=0; i<MAX_NODE_LIST; i++){
			if ( node_list[i]==-1){
				break;
			}
		}
		int number_of_nodes = i;
		int dont_add=0;
		for ( i=0; i<number_of_nodes; i++){
			if ( node_list[i] == start ){
				dont_add = 1;
			}
		}
		if ( dont_add == 0 && tree[start].nd < limit ){
			node_list[number_of_nodes] = start;
			number_of_nodes++;
		}
		dont_add = 0;
		for(i=0; i<number_of_nodes; i++){
			if ( node_list[i] == tree[local_min_var].up[0] ){
				dont_add = 1;
			}
		}
		if ( dont_add == 0 && tree[tree[local_min_var].up[0]].nd< limit ){
			node_list[number_of_nodes] = tree[local_min_var].up[0];
			number_of_nodes++;
		}
		dont_add = 0;
		for(i=0; i<number_of_nodes; i++){
			if ( node_list[i] == tree[local_min_var].up[1] ){
				dont_add = 1;
			}
		}
		if (dont_add == 0 && tree[tree[local_min_var].up[1]].nd < limit){
			node_list[number_of_nodes] = tree[local_min_var].up[1];
			number_of_nodes++;
		}
		//node_list[i]=start;
		//node_list[i+1]=tree[minVarNode].up[0];
		//node_list[i+2]=tree[minVarNode].up[1];
		if ( tree[start].nd > limit){
			recurseMinVar(tree,start,limit,node_list);
		}
		if ( tree[tree[local_min_var].up[0]].nd > limit){
			recurseMinVar(tree,tree[local_min_var].up[0],limit,node_list);
		}
		if ( tree[tree[local_min_var].up[1]].nd > limit){
			recurseMinVar(tree,tree[local_min_var].up[1],limit,node_list);
		}
	}
}
void swap(int* xp, int* yp){
	int temp = *xp;
	*xp = *yp;
	*yp = temp;
}
void selectionSort(int* arr, int n){
	int i, j, min_idx;
	for(i=0; i<n-1; i++){
		min_idx=i;
		for(j=i+1; j<n; j++){
			if ( arr[j] < arr[min_idx]){
				min_idx = j;
			}
		}
		swap(&arr[min_idx], &arr[i]);
	}
}
void getSequenceLengths(gzFile file, int* lengths){
	char buffer [FASTA_MAXLINE];
	int length = 0;
	int i;
	int iter=0;
	while( gzgets(file,buffer,FASTA_MAXLINE) != NULL ){
		if (buffer[0] != '>'){
			length=0;
			for(i=0; buffer[i]!='\n';i++){
				if ( buffer[i] != '-' ){
					length++;
				}
			}
			lengths[iter]=length;
			iter++;
		}
	}	
}
double calculateAvgSeqLength(int numspec, int* sequence_lengths){
	int i;
	double average = 0;
	long double sum = 0;
	for(i=0; i<numspec; i++){
		sum = sequence_lengths[i] + sum;
	}
	average = sum/numspec;
	return average;
}
void findVariantSites(char** MSA, int ref_length, int numspec, Options opt, int* invariant_sites, int* reference){
	int i=0;
	int j=0;
	int number_of_sites=0;
	for(i=0; i<ref_length; i++){
		char base = MSA[0][reference[i]];
		int invariant=0;
		for(j=1; j<numspec; j++){
			if (base != MSA[j][reference[i]] ){
				invariant++;
			}
		}
		if (invariant==0){
			invariant_sites[i]=reference[i];
		}
	}
	for(i=0; i<ref_length; i++){
		if ( reference[i] == -1 ){ break; }
		if (invariant_sites[i] == -1){
			number_of_sites++;
		}
	}
	FILE *invariant_file;
	if (( invariant_file = fopen(opt.variant,"w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
	fprintf(invariant_file,"%d\n",number_of_sites);
	for(i=0; i<ref_length; i++){
		if ( reference[i] == -1 ){ break; }
		if (invariant_sites[i] == -1 ){
			fprintf(invariant_file,"%d\n",i);
		}
	}
	fclose(invariant_file);
}
int main(int argc, char **argv){
	Options opt;
	opt.common=0;
	opt.limit=10000;
	parse_options(argc,argv,&opt);
	struct timespec tstart={0,0}, tend={0.0};
	//FILE* file;
	//char buffer[FASTA_MAXLINE];
	//int c;
	//file=fopen("hello.fasta","r");
	//if (file){
	//	while ((c=getc(file)) != EOF )
	//		putchar(c);
	//	fclose(file);
	//}
	//ssize_t read;
	//char *line = NULL;
	//size_t len=0;
	//while((read=getline(&line, &len, file))!= -1){
	//	printf("%s",line);
	//}
	//if ((file=fopen("hello.fasta","r"))==NULL ){ puts("Cannot open MSA file!"); exit(-1);}
	//while( fgets(file,buffer,FASTA_MAXLINE) != NULL){
	//	int countt=0;
	//}
	//printf("finished");
	//fclose(file);
	//exit(1);
	gzFile MSA_file = Z_NULL;
	if ((MSA_file=gzopen(opt.msa,"r"))==Z_NULL ){ puts("Cannot open MSA file!"); exit(-1);}
	int length_of_MSA=0;
	length_of_MSA=setMSALength(MSA_file);
	gzclose(MSA_file);
	printf("Length of MSA: %d\n",length_of_MSA);
	int* strain_info = (int*)malloc(2*sizeof(int));
	if ((MSA_file=gzopen(opt.msa,"r"))==Z_NULL ){ puts("Cannot open MSA file!"); exit(-1);}
	setNumStrains(MSA_file,strain_info);
	gzclose(MSA_file);
	int numspec = strain_info[0];
	int max_name_length = strain_info[1];
	free(strain_info);
	printf("Number of strains: %d\n",numspec);
	int i,j, k, l;
	char **nodeIDs = (char**)malloc(numspec*sizeof(char*));
	for(i=0; i<numspec; i++){
		nodeIDs[i] = malloc((MAX_STRAINS+max_name_length)*sizeof(char));
		for(j=0; j<(MAX_STRAINS+max_name_length); j++){
			nodeIDs[i][j]='\0';
		}
	}
	printf("Reading in MSA...  ");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	//struct hashmap specmap;
	//hashmap_init(&specmap, hashmap_hash_string, hashmap_compare_string, numspec);
	hashmap_init(&sm, hashmap_hash_string, strcmp);
	int* reference = (int*)malloc(length_of_MSA*sizeof(int));
	for(i=0; i<length_of_MSA; i++){
		reference[i]=-1;
	}
	int* start_end = (int*)malloc(2*sizeof(int));
	start_end[0] = -1;
	start_end[1] = -1;
	if ((MSA_file=gzopen(opt.msa,"r"))==Z_NULL ){ puts("Cannot open MSA file!"); exit(-1);}
	int ref_length=readInNodes(MSA_file,nodeIDs,length_of_MSA,/*&sm,*/reference,start_end);
	gzclose(MSA_file);
	int start_of_MSA = start_end[0];
	int end_of_MSA = start_end[1];
	free(start_end);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	printf("Start of MSA: %d\n",start_of_MSA);
	printf("End of MSA: %d\n",end_of_MSA);
	if (opt.common==1){
		printf("Imputing with the most common allele...  ");
		clock_gettime(CLOCK_MONOTONIC, &tstart);
		int** allele_max = (int **)malloc(length_of_MSA*sizeof(int *));
		for(i=0; i<length_of_MSA; i++){
			allele_max[i] = (int *)malloc(4*sizeof(int));
			//memset(allele_max[i],0,4);
			for(j=0; j<4; j++){
				allele_max[i][j]=0;
			}
		}
		int** MSA = (int**)malloc(numspec*sizeof(int*));
		for(i=0; i<numspec; i++){
			MSA[i] = malloc(length_of_MSA*sizeof(int));
			//memset(MSA[i],'\0',length_of_MSA);
			for(j=0; j<length_of_MSA; j++){
				MSA[i][j] = '\0';
			}
		}
		if ((MSA_file=gzopen(opt.msa,"r"))==Z_NULL ){ puts("Cannot open MSA file!"); exit(-1);}
		readInMSA(MSA_file,MSA,length_of_MSA);
		gzclose(MSA_file);
		int* imputation = (int*)malloc(length_of_MSA*sizeof(int));
		calculate_max_allele(MSA,numspec,length_of_MSA,allele_max,imputation);
		FILE *allele_max_out_file;
		if (NULL==(allele_max_out_file=fopen(opt.outfile,"w"))){ puts("cannot open file!"); exit(-1); }
		print_max_allele(allele_max_out_file,MSA,numspec,length_of_MSA,imputation,nodeIDs);
		fclose(allele_max_out_file);
		free(imputation);
		//for(i=0; i<numspec; i++){
		//	free(nodeIDs[i]);
		//}
		//free(nodeIDs);
		for(i=0; i<length_of_MSA; i++){
			free(allele_max[i]);
		}
		free(allele_max);
		for(i=0; i<numspec; i++){
			free(MSA[i]);
		}
		free(MSA);
		clock_gettime(CLOCK_MONOTONIC, &tend);
		printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	}else{
		printf("Calculating average sequence length...\n");
		int* sequence_lengths = (int*)malloc(numspec*sizeof(int));
		for(i=0; i<numspec; i++){
			sequence_lengths[i] = -1;
		}
		if ((MSA_file=gzopen(opt.msa,"r"))==Z_NULL ){ puts("Cannot open MSA file!"); exit(-1);}
		getSequenceLengths(MSA_file,sequence_lengths);
		gzclose(MSA_file);
		double average = calculateAvgSeqLength(numspec,sequence_lengths);
		printf("Average sequence length is %lf\n",average);
		printf("Imputing using the phylogenetic tree...\n");
		struct node* tree;
		tree=malloc((2*numspec-1)*(sizeof(struct node)));
		for(i=0; i<(2*numspec-1); i++){
			tree[i].up[0] = -1;
			tree[i].up[1] = -1;
			tree[i].down = -1;
			tree[i].nd = 0;
			tree[i].bl =-1;
			tree[i].depth = 0;
			tree[i].pos = -1;
			if ( i >= numspec-1){
				tree[i].name = (char*)malloc(max_name_length*sizeof(char));
				memset(tree[i].name,'\0',max_name_length);
			}
		}
		/*struct multinode* multinode_tree;
		multinode_tree=malloc((2*numspec-1)*(sizeof(struct multinode)));
		for(i=0; i<(2*numspec-1); i++){
			multinode_tree[i].down = -1;
			for(j=0; j<MAX_POLYTOMIES; j++){
				multinode_tree[i].up[j] = -1;
			}
			multinode_tree[i].name = (char*)malloc(max_name_length*sizeof(char));
			memset(multinode_tree[i].name,'\0',max_name_length);
		}*/
		printf("Reading in tree file...  \n");
		clock_gettime(CLOCK_MONOTONIC, &tstart);
		FILE *treefile;
		/*if (NULL==(treefile=fopen(opt.tree,"r"))){ puts("Cannot open tree!"); exit(-1);}
		int c=0;
		int para=0;
		while( (c = fgetc(treefile))!= EOF){
			if ( (char)c == '(' ){
				para++;
			}
		}
		fclose(treefile);
		printf("para: %d\n",para);*/
		if (NULL==(treefile=fopen(opt.tree,"r"))){ puts("Cannot open tree!"); exit(-1);}
		int root=getclade(tree,treefile,nodeIDs,numspec,/*&sm,*/average)-1;
		fclose(treefile);
		clock_gettime(CLOCK_MONOTONIC, &tend);
		printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		//printtree(tree,root,numspec);
		//for(i=0; i<numspec; i++){
		//	free(nodeIDs[i]);
		//}
		//free(nodeIDs);
		//hashmap_destroy(&sm);
		tree[root].down=-1;
		get_number_descendants(tree,root);
		assignDepth(tree, tree[root].up[0], tree[root].up[1], 1);
		int* node_list = (int*)malloc(MAX_NODE_LIST*sizeof(int));	
		for(i=0; i<MAX_NODE_LIST; i++){
			node_list[i]=-1;
		}
		//minVariance = (tree[tree[root].up[0]].nd*tree[tree[root].up[0]].nd + tree[tree[root].up[1]].nd*tree[tree[root].up[1]].nd)/3;
		clock_gettime(CLOCK_MONOTONIC, &tstart);
		recurseMinVar(tree,root,opt.limit,node_list);
		for(i=0; i<MAX_NODE_LIST; i++){
			if(node_list[i]==-1){
				break;
			}
		}
		if ( i==0 ){
			node_list[0]=root;
			i=1;
		}
		int max_iter=i;
		printf("Splitting tree into %d partitions...\n",max_iter);
		clock_gettime(CLOCK_MONOTONIC, &tend);
		printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		clock_gettime(CLOCK_MONOTONIC, &tstart);
		/*for(i=0; i<max_iter; i++){
			int* node_list_children = (int*)malloc(MAX_NODE_LIST*sizeof(int));
			for(j=0; j<MAX_NODE_LIST; j++){
				node_list_children[j]=-1;
			}
			get_children(tree,node_list[i],node_list_children);
			if ( node_list_children[0] == -1){
				node_list_children[0] = node_list[i];
			}
			for(j=0; j<MAX_NODE_LIST; j++){
				if (node_list_children[j] == 3117130 ){
					printf("partition %d\n",i);
				}
			}
		}
		exit(1);*/
		if ( max_iter > 0 ){
		for(i=0; i<max_iter; i++){
			printf("Imputing partition %d size %d\n",i,tree[node_list[i]].nd);
			int* node_list_children = (int*)malloc(MAX_NODE_LIST*sizeof(int));
			for(j=0; j<MAX_NODE_LIST; j++){
				node_list_children[j]=-1;
			}
			get_children(tree,node_list[i],node_list_children);
			if ( node_list_children[0] == -1){
				node_list_children[0] = node_list[i];
			}
			for(j=0; j<MAX_NODE_LIST; j++){
				if ( node_list_children[j] == -1){
					break;
				}else{
					tree[node_list_children[j]].likenc = malloc(length_of_MSA*sizeof(double*));
					tree[node_list_children[j]].posteriornc = malloc(length_of_MSA*sizeof(double*));
					for(k=0; k<length_of_MSA; k++){
						tree[node_list_children[j]].likenc[k]=malloc(4*sizeof(double));
						tree[node_list_children[j]].posteriornc[k]=malloc(4*sizeof(double));
						for(l=0; l<4; l++){
							tree[node_list_children[j]].likenc[k][l]=0;
							tree[node_list_children[j]].posteriornc[k][l]=0;
						}
					}
				}
			}
			int number_of_children = j;
			int number_of_leaves = 0;
			for(j=0; j<number_of_children; j++){
				if ( tree[node_list_children[j]].up[0] == -1){
					number_of_leaves++;
				}
			}
			int* leaf_nodes = (int*)malloc(number_of_leaves*sizeof(int));
			k=0;
			for(j=0; j<number_of_children; j++){
				if ( tree[node_list_children[j]].up[0] == -1){
					leaf_nodes[k]=node_list_children[j];
					k++;	
				}
			}
			selectionSort(leaf_nodes, k);
			int** MSA = (int**)malloc(k*sizeof(int*));
			for(j=0; j<k; j++){
				MSA[j]= malloc(length_of_MSA*sizeof(int));
				memset(MSA[j],0,length_of_MSA);
			}
			if ((MSA_file=gzopen(opt.msa,"r"))==Z_NULL ){ puts("Cannot open MSA file!"); exit(-1);}
			readInMSA_for_partition(tree,MSA_file,MSA,leaf_nodes,length_of_MSA,max_name_length,k);
			gzclose(MSA_file);
			//struct hashmap leaf_nodes_coding;
			//hashmap_init(&leaf_nodes_coding, hashmap_hash_string, hashmap_compare_string, k);
			for(j=0; j<k; j++){
				//hashmap_put(&leaf_nodes_coding,leaf_nodes[j],j);
				tree[leaf_nodes[j]].pos = j;
			}
			double LRVECnc[4][4], RRVECnc[4][4], RRVALnc[4], PMATnc[2][4][5];
			double parameters[10] = {0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
			clearGlobals(parameters,LRVECnc,RRVECnc,RRVALnc,PMATnc);
			getposterior_nc(parameters,tree,numspec,length_of_MSA,node_list[i],MSA,LRVECnc,RRVECnc,RRVALnc,PMATnc,number_of_children,node_list_children,k);
			for(j=0; j<number_of_children; j++){
				if ( tree[node_list_children[j]].up[0] == -1){
					impute(tree,node_list_children[j],MSA,length_of_MSA,numspec,leaf_nodes,k);
				}
			}
			if ( i==0 ){
				FILE *outfile;
				if (NULL==(outfile=fopen(opt.outfile,"w"))){ puts("Cannot open file!"); exit(-1);}
				print_imputed_MSA(tree,MSA,length_of_MSA,outfile,leaf_nodes,k);
				fclose(outfile);
			}else{
				FILE *outfile;
				if (NULL==(outfile=fopen(opt.outfile,"a"))){ puts("Cannot open file!"); exit(-1);}
				print_imputed_MSA(tree,MSA,length_of_MSA,outfile,leaf_nodes,k);
				fclose(outfile);
			}
			free(leaf_nodes);
			for(j=0; j<k; j++){
				free(MSA[j]);
			}
			free(MSA);
			for(j=0; j<MAX_NODE_LIST; j++){
				if ( node_list_children[j] == -1){
					break;
				}else{
					for(k=0; k<length_of_MSA; k++){
						free(tree[node_list_children[j]].likenc[k]);
						free(tree[node_list_children[j]].posteriornc[k]);
					}
					free(tree[node_list_children[j]].likenc);
					free(tree[node_list_children[j]].posteriornc);
					if ( tree[node_list_children[j]].up[0] == -1 && tree[node_list_children[j]].up[1] == -1){
						free(tree[node_list_children[j]].name);
					}
				}
			}
			free(node_list_children);
		}
		}
		clock_gettime(CLOCK_MONOTONIC, &tend);
		printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		clock_gettime(CLOCK_MONOTONIC, &tstart);
		free(node_list);
		clock_gettime(CLOCK_MONOTONIC, &tend);
		printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	}
	printf("Allocating memory to read in imputed file\n");
	char** MSAi = (char**)malloc(numspec*sizeof(char*));
	for(i=0; i<numspec; i++){
		MSAi[i] = (char*)malloc(length_of_MSA*sizeof(char));
	}
	//FILE* Imputed = NULL;
	//if ((Imputed=fopen(opt.outfile,"r"))==NULL ){ puts("Cannot open MSA file!"); exit(-1);}
	//FILE *outfile;
        //if (NULL==(outfile=fopen(opt.outfile,"r"))){ puts("Cannot open file!"); exit(-1);}
	//char buffer[FASTA_MAXLINE];
	//while( fgets(outfile,buffer,FASTA_MAXLINE) != NULL){
	//	int count=0;
	//}
	//printf("finished\n");
	//fclose(outfile);
	//exit(1);
	/* Clear nodeIDs */
	for(i=0; i<numspec; i++){
		for(j=0; j<(MAX_STRAINS+max_name_length); j++){
			nodeIDs[i][j]='\0';
		}
	}
	int num_remaining=readInImputed(nodeIDs,MSAi,reference,max_name_length,length_of_MSA,start_of_MSA,end_of_MSA,ref_length,opt,numspec);
	int* variant_sites = (int*)malloc(ref_length*sizeof(int));
	for(i=0; i<ref_length; i++){
		variant_sites[i]=-1;
	}
	findVariantSites(MSAi,ref_length,numspec,opt,variant_sites,reference);
	free(variant_sites);
	free(reference);
	for(i=0; i<numspec; i++){
		free(MSAi[i]);
		free(nodeIDs[i]);
	}
	free(MSAi);
	free(nodeIDs);
}
