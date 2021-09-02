#!/bin/sh

source msmc_params.sh

IND=$1
BAMFILE=${BAMDIR}/aln_${IND}.sorted.bam

module load gcc
module load miniconda3

source activate new_env # needs to have whatshap installed via conda

echo "working with individual $IND"

if [ "$PHASING" == "whatshap" ]; then
  for s in `cat $SCAFFOLDS`
        do echo "working with scaffold $s"
        if [ -f $OUTDIR/vcf/${IND}.${s}.${prefix}.minDP10.${phasing}.${method}.vcf.gz ]; then
                echo "phased VCF already exists; moving onto next scaffold"
        else
                echo "phased VCF does not exist; phasing VCF for scaffold $s"
                sed -i 's/^ //g' ${OUTDIR}/vcf/${IND}.${s}.${prefix}.minDP10.${method}.vcf

                whatshap phase --reference $GENOME --ignore-read-groups -o ${OUTDIR}/vcf/${IND}.${s}.${prefix}.minDP10.${phasing}.${method}.vcf.gz ${OUTDIR}/vcf/${IND}.${s}.${prefix}.minDP10.${method}.vcf $BAMFILE
                whatshap stats --tsv=$OUTDIR/stats/${IND}.${s}.${prefix}.minDP10.${phasing}.stats.tsv $OUTDIR/vcf/${IND}.${s}.${prefix}.minDP10.${phasing}.${method}.vcf.gz 
        fi
  done
else
  echo "Phasing parameter is not whatshap; exiting now"
fi

echo "finished with individual $IND"
date
