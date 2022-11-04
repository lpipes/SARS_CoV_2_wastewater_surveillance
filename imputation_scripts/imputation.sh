#!/bin/bash
TREEZIP=$1
MSA_MASKED=$2
echo "$TREEZIP"
echo "$MSA_MASKED"
unzip $TREEZIP
tar xf $MSA_MASKED
MSA_MASKED_NAME=`basename $MSA_MASKED .tar.xz`
TREE_NAME=`basename $TREEZIP .zip`
echo "$MSA_MASKED_NAME"
echo "$TREE_NAME"
arrIN=(${MSA_MASKED_NAME//_/ })
FASTA_FILE=${arrIN[1]}
echo "$FASTA_FILE"
./select_SARS_CoV_2.pl ${TREE_NAME}/metadata.csv ${MSA_MASKED_NAME}/${FASTA_FILE}_masked.fa
rm ${MSA_MASKED_NAME}/${FASTA_FILE}_masked.fa
gzip seqs.fasta &
#if missing > 0
if [ -f "missing.txt" ]
then
	echo "Pruning missing strains..."
	Rscript --vanilla cleantree.R ${TREE_NAME}/global.tree missing.txt
	mv global_out.tree global.tree
	mv global.tree ${TREE_NAME}
fi
perl -i -pe 's/node\_[0-9]+//g' ${TREE_NAME}/global.tree
./sarcov2_imputation -i seqs.fasta.gz -t ${TREE_NAME}/global.tree -l 50000 -v variants.txt -o imputed.fasta -m final_imputed.fasta -r redundant.txt
rm imputed.fasta
gzip final_imputed.fasta
rm -rf ${TREE_NAME}
rm -rf ${MSA_MASKED_NAME}
#final database file is final_imputed.fasta.gz
#creating bowtie2 database
bowtie2-build -f EPI_ISL_402124.fasta EPI_ISL_402124.fasta
#transfer following bowtie2 indexes to virtual machine
# EPI_ISL_402124.fasta.1.bt2
# EPI_ISL_402124.fasta.2.bt2
# EPI_ISL_402124.fasta.3.bt2
# EPI_ISL_402124.fasta.4.bt2
# EPI_ISL_402124.fasta.rev.1.bt2
# EPI_ISL_402124.fasta.rev.2.bt2
