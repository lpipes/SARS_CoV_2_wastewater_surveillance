#!/bin/bash
FILE=$1
CASE=$2
echo $FILE
echo $CASE
SAM="test_input/test_input/reads_case${CASE}.sam"
MSA="filtered_virus_MSA.fasta"
TIMES="times/times_${CASE}.txt"
EMOUT="times/mismatches_output_${CASE}.txt"
LINE="============"
start=$SECONDS
bowtie2 -a --threads 1 -x wuhCor1 -f -U $FILE -S $SAM 
duration=$(( SECONDS - start ))
echo "${LINE} Time for alignment: ${duration} ${LINE}" >> $TIMES 

perl sam2mismatches.pl $MSA $SAM 
duration=$(( SECONDS - start ))
echo "${LINE} Time for mismatches: ${duration} ${LINE}" >> $TIMES

perl elimination_sam2mismatches.pl $MSA $SAM $EMOUT 
duration=$(( SECONDS - start ))
echo "${LINE} Time for filter: ${duration} ${LINE}" >> $TIMES
