#!/bin/bash

set -e

echo "python /usr/local/bin/biasFilter.py $1 --mapq=1 --baseq=1 --tempFolder=/var/spool/cwl/ $2 $3 $4 /var/spool/cwl/filtered.vcf"
python /usr/local/bin/biasFilter.py $1 --mapq=1 --baseq=1 --tempFolder=${TMPDIR} $2 $3 $4 ${PWD}/filtered.vcf
mkdir -p ${PWD}/filtered_qcSummary

