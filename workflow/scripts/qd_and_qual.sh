#!/bin/bash

INP=$1
QD=$2
QUAL=$3
LOG=$4

(bcftools query -f '%QUAL\t%INFO/QD\n' $INP | 
awk -v qd_file=${QD}-unsrt -v qual_file=${QUAL}-unsrt  '
	BEGIN {OFS="\t"}

	{
		qual[int($1)]++; 
		qd[int($2)]++;
	}

	END {
		for(i in qual) print i, qual[i] > qual_file;
    	for(i in qd) print i, qd[i] > qd_file;
	}
') 2> $LOG;

(sort -grk1 ${QUAL}-unsrt > ${QUAL} && rm ${QUAL}-unsrt) 2>>$LOG;
(sort -grk1 ${QD}-unsrt > ${QD} && rm ${QD}-unsrt) 2>>$LOG
