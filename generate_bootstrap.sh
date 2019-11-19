#!/bin/bash -l

OUTDIR=/path/to/output
POP_OR_IND=$1
DATE=`date +%m%d%y`
RUN_NAME=msmc_${POP_OR_IND}_${DATE}
nchr=`wc -l SCAFFOLDS.txt` # number of chromosomes in reference genome

#input for the bootstrapping
BS_INPUT=`for s in `cat SCAFFOLDS.txt`; do find ${OUTDIR}/input/ -maxdepth 1 -name "*${POP_OR_IND}.${s}*.txt"; done`

#output from the bootstrapping 
BS_OUTPUT=${OUTDIR}/bootstrap/${POP_OR_IND}.${RUN_NAME}.bootstrap

multihetsep_bootstrap.py -n 50 -s 20000000 --chunks_per_chromosome 10 --nr_chromosomes $nchr $BS_OUTPUT $BS_INPUT

cd ${OUTDIR}/bootstrap

ls -d *bootstrap_* > ${OUTDIR}/${POP_OR_IND}.${RUN_NAME}.bs_file_list.txt
