#!/bin/sh

# to run this script:
# ./msmc_3_runMSMC.sh

# all parameters come from the msmc_param control file
# make edits there before using this script!
source msmc_param

if [ $NR_IND == 1 ]; then

	find ${OUTDIR}/input/msmc_input.${POP_OR_IND}.*.txt -size 0 -delete
	ls ${OUTDIR}/input/msmc_input.${POP_OR_IND}.*.txt | grep -v $sex_chr > ${OUTDIR}/input/SCAFS_INPUT_${POP_OR_IND}
else
	for i in `cat ${POP_OR_IND}_IND`
       do echo $i
       IND=$i

       for s in `cat SCAFFOLDS.txt`
               do echo $s
               ls ${OUTDIR}/input/msmc_input.${IND}.${s}.txt >> ${OUTDIR}/input/SCAFS_INPUT_${POP_OR_IND}
      done
  done
fi

module load msmc2

### Report settings/parameters:
date
echo "Script: msmc_3_onepop.sh"
echo "Run name: $RUN_NAME"
echo "SNP calling method: $METHOD"
echo "Period setting: $P_PAR"
echo "Nr of individuals (1 or 2+): $NR_IND"
echo "Population or individuals ID: $POP_OR_IND"
echo "Individual: "
echo "Scaffolds: SCAFS_INPUT_${POP_OR_IND}"
echo "Iterations: 100"


if [ $NR_IND == 1 ]
        then
        echo "Running MSMC for one individual"
        MSMC_INPUT=`cat ${OUTDIR}/input/SCAFS_INPUT_${POP_OR_IND}`
        MSMC_OUTPUT=${OUTDIR}/output/msmc_output.${RUN_NAME}

        if [ -f "${OUTDIR}/input/SCAFS_INPUT_${POP_OR_IND}" ]
                then
                        echo "MSMC_INPUTS: SCAFS_INPUT_${POP_OR_IND}_noLG9"
                        echo "MSMC_OUTPUT: $MSMC_OUTPUT"
                else
                        echo "MSMC_INPUT does not exist! Exiting now"
                        exit 1
        fi

        msmc2 -t 16 -p $P_PAR -i 100 -o $MSMC_OUTPUT -I 0,1 $MSMC_INPUT

        mv $MSMC_OUTPUT*loop.txt ${OUTDIR}/output/log_and_loop/
        mv $MSMC_OUTPUT*log ${OUTDIR}/output/log_and_loop/
else
        echo "Running MSMC for $NR_IND individuals"
        MSMC_INPUT=`cat ${OUTDIR}/input/SCAFS_INPUT_${POP_OR_IND}`
        MSMC_OUTPUT=${OUTDIR}/output/msmc_output.${POP_OR_IND}.${RUN_NAME}
       	n=$(expr ${NR_IND} - 2)
       	INDEX=$(for num in `seq 0 ${n}`; do echo -n "${num},"; done; echo ${NR_IND})
        
        msmc2 -t 16 -p $P_PAR -i 100 -o ${MSMC_OUTPUT} -I `echo $INDEX` $MSMC_INPUT

        mv $MSMC_OUTPUT*loop.txt ${OUTDIR}/output/log_and_loop/
        mv $MSMC_OUTPUT*log ${OUTDIR}/output/log_and_loop/
fi


echo "done running msmc2 for ${POP_OR_IND}"
date