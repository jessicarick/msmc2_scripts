#!/bin/bash

# all parameters come from the msmc_param control file
# make edits there before using this script!
scriptdir=$(dirname "$0")
source ${scriptdir}/msmc_params.sh

IND=$1

#input for the bootstrapping
BS_INPUT=`for s in `cat SCAFFOLDS.txt`; do find ${OUTDIR}/input/ -maxdepth 1 -name "*${POP_OR_IND}.${s}*.txt"; done`

#output from the bootstrapping 
BS_OUTPUT=${OUTDIR}/bootstrap/${POP_OR_IND}.${RUN_NAME}.bootstrap

echo "generating bootstraps for ${POP_OR_IND}.${RUN_NAME}"
multihetsep_bootstrap.py -n 50 -s 20000000 --chunks_per_chromosome 10 --nr_chromosomes $nchr $BS_OUTPUT $BS_INPUT

cd ${OUTDIR}/bootstrap
ls -d *bootstrap_* > ${OUTDIR}/${POP_OR_IND}.${RUN_NAME}.bs_file_list.txt

##### run msmc ####
MSMC_BS=$(for x in `cat ${OUTDIR}/${POP_OR_IND}.${RUN_NAME}.bs_file_list.txt`; do find $BASEDIR/bootstrap/$x -maxdepth 2 -name "bootstrap_multihetsep*.txt"; done)

#MSMC_BS=`ls ${OUTDIR}/bootstrap/${FOLDER}/*.txt`
MSMC_OUTPUT=${OUTDIR}/bootstrap/msmc_output.${RUN_NAME}

echo "running msmc2 on bootstraps for ${POP_OR_IND}.${RUN_NAME}"
msmc2 -t 20 -p $P_PAR -o $MSMC_OUTPUT -I 0,1 $MSMC_BS

echo "done with msmc bootstraps for ${POP_OR_IND}.${RUN_NAME}"
