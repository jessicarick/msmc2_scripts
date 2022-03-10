module load gcc
module load msmc2
module load samtools
module load bcftools
module load vcftools
module load bwa
module load python/3.6.3

MSMCTOOLS=/path/to/msmc-tools-master # folder with msmc-tools binaries
PATH=$PATH:$MSMCTOOLS

OUTDIR=/path/to/output # main directory for output files

# make directories for intermediate files-- will fail if these don't exist
mkdir -p ${OUTDIR}/vcf
mkdir -p ${OUTDIR}/mask
mkdir -p ${OUTDIR}/input

# for msmc_1_call.sh
GENOME=/path/to/reference.fa # reference genome fasta
prefix=prefix # prefix of genome masks
BAMDIR=/path/to/bamfiles/ # directory with bamfiles
k=35

# for msmc_3_generateInput.sh
NR_IND=1 # number of individuals in analysis
POP_OR_IND=POPNAME # name of individual or population being analyzed for script 3
DATE=`date +%m%d%y`
RUN_NAME=msmc_${POP_OR_IND}_${DATE}
P_PAR=1*2+25*1+1*2+1*3 

nchr=`wc -l SCAFFOLDS.txt` # number of chromosomes in reference genome
sex_chr=LG9 # name of sex chromosome to omit in analyses

METHOD="samtools"
PHASING="unphased"
