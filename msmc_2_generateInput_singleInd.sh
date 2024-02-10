#!/bin/bash

# all parameters come from the msmc_param control file
# make edits there before using this script!
scriptdir=$(dirname "$0")
source ${scriptdir}/msmc_params.sh

### Variables:
IND=$1

module list

for s in `cat SCAFFOLDS.txt`
        do echo "working on scaffold $s"
        SCAFFOLD=$s
        VCF=`ls ${OUTDIR}/vcf/${IND}.${SCAFFOLD}.${METHOD}.vcf.gz`

        #MASK_REPEATS=repeats.bed.gz # Needs to be gzipped
        MSMC_INPUT=${OUTDIR}/input/msmc_input.${IND}.${SCAFFOLD}.txt

        printf "\n \n \n \n"
        date
        echo "Script: msmc_2_generateInput_singleInd"
        echo "Individual: ${IND}"
        echo "Scaffold: ${SCAFFOLD}"
        echo "Phasing: ${PHASING}"
        echo "Method: ${METHOD}"
        echo "MSMC input file: ${MSMC_INPUT}"
        echo "VCF: ${VCF}"

        if [ -f "$VCF" ]
                then
                        echo "VCF exists, starting creation of input for MSMC2!"

### Generate MSMC input files:
                if [ $METHOD == samtools ]
                        then
                                MASK_INDIV=${OUTDIR}/mask/ind_mask.${IND}.${SCAFFOLD}.${METHOD}.bed.gz # store indiv.mask file path
                                MASK_GENOME=${OUTDIR}/mask/prefix_${SCAFFOLD}.mask.${k}.50.bed.gz
                                echo "MASK: ${MASK_INDIV}"
                                echo "MAPPABILITY MASK: ${MASK_GENOME}"
                                echo "Creating MSMC input file WITH individual mask (samtools)"

                                generate_multihetsep.py --mask=${MASK_INDIV} --mask=${MASK_GENOME} $VCF > $MSMC_INPUT # without repeat mask

                elif [ $METHOD != samtools ]
                        then
                                # NOTE THAT THIS WAS CHANGED 10 FEB 2024
                                # AND HAS NOT YET BEEN TESTED TO MAKE SURE IT WORKS
                                
                                echo "Creating individual mask. Note that your input VCF should include ALL sites (variant & invariant)."
                                MASK_INDIV=${OUTDIR}/mask/ind_mask.${IND}.${SCAFFOLD}.${METHOD}.bed.gz
                                MASK_GENOME=${OUTDIR}/mask/prefix_${SCAFFOLD}.mask.${k}.50.bed.gz
                                
                                VCF_OUT=${VCF}.parsed.vcf
                                vcfAllSiteParser.py $SCAFFOLD $MASK_INDIV $VCF_OUT
                                
                                echo "Creating MSMC input file with new individual mask"
                                
                                generate_multihetsep.py --mask=$MASK_INDIV --mask=$MASK_GENOME $VCF > $MSMC_INPUT # with new repeat mask
                fi

        else
                echo "VCF does not exist, moving on to next scaffold"
        fi

        echo "Done with ${SCAFFOLD}; moving on to next scaffold"
done;

echo "Done with script."
date

####
