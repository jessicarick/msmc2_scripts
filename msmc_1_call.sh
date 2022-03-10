#!/bin/sh

# all parameters come from the msmc_param control file
# make edits there before using this script!
scriptdir=$(dirname "$0")
source ${scriptdir}/msmc_params.sh

## VARIABLES:
IND=$1
BAMFILE=${BAMDIR}/aln_${IND}.sorted.bam

printf "\n \n \n \n"
date
echo "Current script: msmc_1_call.sh"
echo "Individual: $IND"
echo "Bamfile: $BAMFILE"

if [ -f "${BAMFILE}.bai" ]
        then
                echo "Bamfile index already exists, moving on!"
        else
                echo "Bamfile index does not exist, creating index"
                samtools index $BAMFILE ${BAMFILE}.bai
fi


for s in `cat SCAFFOLDS.txt`; \
        do echo "Handling scaffold $s"; \

        ### Calculate mean coverage (to be used as input for bamCaller.py):
        MEANCOV=`samtools depth -r $s $BAMFILE | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` # calculate mean coverage
        echo ${IND}.${s} $MEANCOV >> ${OUTDIR}/coverage_samtoolsDepth_${IND}.txt # save mean coverage in separate file
        echo "Mean coverage for this individual, scaffold ${s}: $MEANCOV"

        ### Generate a single-sample VCF and a mask-file:
        MASK_IND=${OUTDIR}/mask/ind_mask.${IND}.${s}.${METHOD}.bed.gz # Individual mask file to be created
        VCF=${OUTDIR}/vcf/${IND}.${s}.${METHOD}.vcf # VCF file to be created

        #If genome isn't indexed, add:
        #samtools faidx $GENOME
        if [ "$METHOD" == "samtools" ]; then
                echo "starting samtools alignment"
                bcftools mpileup -Ou -r ${s} --threads 16 -f $GENOME $BAMFILE | bcftools call -c --threads 16 -V indels | bamCaller.py $MEANCOV $MASK_IND > ${VCF}

                ## Only DP > 9:
                #samtools mpileup -q 20 -Q 20 -C 50 -u -r $SCAFFOLD -f $GENOME $BAMFILE | bcftools call -c -V indels | bcftools view -i 'INFO/DP>9' | $MSMCTOOLS/bamCaller.py $MEANCOV $MASK_IND | gzip -c > $VCF.gz

                # -q = min. mapping qual; -Q = min. base qual; -C = coefficient for downgrading mapping qual for reads w/ excessive mismatches; -u = generate uncompressed VCF/BCF; -r = only specified region; -f = fasta.
                # bcftools: "call" will call SNPs/indels; "-V indels" skips indels; -c = consensus caller.
        fi
        echo "done with scaffold ${s}; moving on to next scaffold"
done

### Report:
echo "Filtered VCF and mask created for all scaffolds for ${IND}."
echo "Done with script."
date
