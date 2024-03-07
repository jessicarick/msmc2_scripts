#!/bin/bash

# all parameters come from the msmc_param control file
# make edits there before using this script!
scriptdir=$(dirname "$0")
source ${scriptdir}/msmc_params.sh

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

                # NOTE THAT THIS WAS CHANGED 10 FEB 2024
                # AND HAS NOT YET BEEN TESTED TO MAKE SURE IT WORKS
                                
                echo "Creating individual mask. Note that your input VCF should include ALL sites (variant & invariant)."
                MASK_INDIV=${OUTDIR}/mask/ind_mask.${IND}.${SCAFFOLD}.${METHOD}.bed.gz
                                
                VCF_OUT=${VCF}.parsed.vcf
                vcfAllSiteParser.py $SCAFFOLD $MASK_INDIV $VCF_OUT
                                
                echo "Creating MSMC input file with new individual mask"
                                
                generate_multihetsep.py --mask=$MASK_INDIV --mask=$MASK_GENOME $VCF > $MSMC_INPUT # with new repeat mask
        fi

done

echo "Done with script."
date

####
