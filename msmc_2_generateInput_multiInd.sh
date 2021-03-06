#!/bin/bash

# all parameters come from the msmc_param control file
# make edits there before using this script!
source msmc_param

### Variables:
IND=`cat $1`
POP=$2

for s in `cat SCAFFOLDS.txt`
        do SCAFFOLD=$s

        MSMC_INPUT=${OUTDIR}/input/msmc_input.${POP}.${SCAFFOLD}.txt

        printf "\n \n \n \n"
        date
        echo "Script: msmc_2_generateInput_multiInd"
        echo "Individuals: ${IND}"
        echo "Population: $POP"
        echo "Scaffold: ${SCAFFOLD}"
        echo "Method: ${METHOD}"
        echo "MSMC input file: ${MSMC_INPUT}"

        for ind in $IND
                do INDMASK=`ls ${OUTDIR}/mask/ind_mask.${ind}.${SCAFFOLD}.bed.gz`
                echo "--mask=$INDMASK " >> ${OUTDIR}/mask/${POP}.mask_file.$SCAFFOLD
                INDVCF=`ls ${OUTDIR}/vcf/${ind}.${SCAFFOLD}.msmc.vcf.gz`
                echo $INDVCF >> ${OUTDIR}/vcf/${POP}.vcf_file.${SCAFFOLD}
        done

### Generate MSMC input files:
        if [ $METHOD == samtools ]
                then
                MASK_GENOME=${OUTDIR}/mask/${prefix}_${SCAFFOLD}.mask.${k}.50.bed.gz

                echo "MAPPABILITY MASK: ${MASK_GENOME}"
                echo "Creating MSMC input file WITH individual mask (samtools)"
                #${MSMCTOOLS}/generate_multihetsep.py --negative_mask=$MASK_REPEATS --mask=$MASK_INDIV $VCF > $MSMC_INPUT # with repeat mask
                generate_multihetsep.py `cat ${OUTDIR}/mask/${POP}.mask_file.${SCAFFOLD}` --mask=$MASK_GENOME `cat ${OUTDIR}/vcf/${POP}.vcf_file.${SCAFFOLD}` > ${MSMC_INPUT} # without repeat mask

        elif [ $METHOD == gatk ]
                then
                echo "Creating MSMC input file WITHOUT individual mask (gatk)"
                MASK_GENOME=`ls ${OUTDIR}/mask/${prefix}_${SCAFFOLD}.mask.${k}.50.bed.gz`
                #msmc-tools/generate_multihetsep.py --negative_mask=$MASK_REPEATS $VCF > $MSMC_INPUT # with repeat mask
                generate_multihetsep.py --mask=$MASK_GENOME `cat ${OUTDIR}/vcf/${POP}.vcf_file.${SCAFFOLD}` > $MSMC_INPUT # without repeat mask
        fi

done

echo "Done with script."
date

####