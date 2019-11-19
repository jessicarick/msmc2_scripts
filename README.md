# msmc2_scripts
MSMC Tutorial and Scripts for running MSMC2.

## MSMC2 Workflow and Code, v1.2 (Nov 2019)

*This guide is intended to get you started on MSMC analyses, and to
provide scripts for you to work off of. These scripts may not be perfect
for your analyses, however, so please also read the manuals for the
different programs being used (e.g. guides to
[PSMC](https://github.com/lh3/psmc)/[MSMC](https://github.com/stschiff/msmc/blob/master/guide.md)/[MSMC2](https://github.com/stschiff/msmc2))
so that you understand and can make good choices about parameter
settings and the workflow you’re using.*

To start with, you need (1) a reference genome, (2) whole genome
resequencing data, (3) [the MSMC2
program](https://github.com/stschiff/msmc2) (already installed as a
module on Teton), and (4) the [MSMC helper
tools](https://github.com/stschiff/msmc-tools). The more complete your
genome is, the better your plots will look (in my experience).
Theoretically you could do this with *de novo* WGS assembly, but I don’t
think it’s recommended, as the analyses are performed per chromosome. If
your genome is in a bunch of pieces, one way to fix this for MSMC is to
concatenate groups of the scaffolds together, separated by a bunch of
’N’s. Here’s how you could do that:

``` sh
echo '> Genome_name' > genome_new.fasta
cat genome.fa | sed 's/^>.*$/\n\(N\)\{100\}/g' | tr -d '\n' | 
    fold -w 80 >> genome_new.fasta
```

Another (probably better) thing to do could be to use
[Chromosomer](https://github.com/gtamazian/chromosomer) or a similar
program to assemble your draft genome to a more-complete published
genome. If you’re interested in how the quality of your reference genome
may affect your inferences, here’s a neat paper using simulations to
test this: [Patton et al. 2019,
MBE](https://doi.org/10.1093/molbev/msz191).

Once you have a good genome, you’ll also need the following scripts:

1.  `run_snpable.sh` – this generates a "mappability mask" for your
    reference genome using Heng Li’s [SNPable Regions
    program](http://lh3lh3.users.sourceforge.net/snpable.shtml)

2.  `submit_1.txt` and `msmc_1_call.sh` – this generates `vcf` and
    `mask` files for each individual and each chromosome. The submission
    script loops this script over all individuals and all chromosomes.

3.  `submit_2.txt` and `msmc_2_generateInput_singleInd.sh` (if working
    with one indv) – this uses the `vcf` and `mask` files to generate a
    msmc input file when you will be analyzing each individual
    separately.

4.  `submit_2_multi.txt` and `msmc_2_generateInput_multiInd.sh` (if
    multiple ind per pop), which will generate input files for all
    individuals and all chromosomes for runs where you intend to combine
    multiple individuals.

5.  `msmc_3_runMSMC_onepop.sh` (if all ind belong to same pop)

6.  `msmc_3_runMSMC_twopop.sh` (if looking at cross-coalescence)

7.  bootstrap scripts: `generate_bootstrap.sh`, `msmc_bootstrap.sh`,
    `submit_bootstrap.sh`

These scripts are appended to this guide and can be found on Github at  
<https://github.com/jessicarick/msmc2_scripts>.

## Alignment

The MSMC process starts with sorted `.bam` files, so if you haven’t
already, you’ll need to align your `fastq` reads to your reference
genome. You can do this in whatever way you fancy, such as with a `bwa`
or `bowtie` pipeline.

## Step 0 - Create Mappability Mask

One more preparation step before you begin (which can also be run
concurrently with `Step 1` below) is to create a mappability mask for
your genome. To do this, you need the scripts for
[SNPable](http://lh3lh3.users.sourceforge.net/snpable.shtml), which
you’ll find linked on Heng Li’s website. You’ll also need the
`makeMappabilityMask.py` script from the
[msmc-tools-master](https://github.com/stschiff/msmc-tools) set of
scripts.

This process involves first extracting all *k*-mer subsequences from the
genome as read sequences, then aligning them back to the genome to get
an estimate of how “mappable" different regions are. This is done using
the `run_snpable2.sh` script, which requires `samtools` and `bwa`, in
addition to the `splitfa`, `gen_raw_mask.pl`, and `gen_mask` scripts
from `SNPable` and `msmc-tools-master`.

``` sh
# running on command line
./run_snpable2.sh prefix /path/to/reference.fa

# submitting as a SLURM job
sbatch run_snpable2_slurm.sh prefix /path/to/reference.fa
```

## Step 1 - Call Variants

From here, you will need to edit the first couple of lines of the
`msmc_1_call.sh` script for your data. The `BASEDIR` is the folder where
you want to be working; `GENOME` is the (full) path to the reference
genome that you aligned your reads to; `BAMDIR` is the directory
containing all of your sorted bamfiles; `BAMFILE` indicates the syntax
for how your bamfiles are named (e.g. `${id}.bowtie2.sorted.bam`).

The script uses the `samtools-bcftools-vcftools` pipeline for calling
variants, but could be modified if you’d prefer a different pipeline.
You’ll also need a file called `SCAFFOLDS.txt` with all of the scaffold
names in your reference genome and a file called `INDS.txt` with all of
your individual IDs. If you don’t already have the list of scaffolds,
you can create one using this simple one-liner:

``` sh
grep '^>' reference.fa | sed 's/>//' > SCAFFOLDS.txt
```

Once you have this file and you’ve edited both the `msmc_1_call.sh` and
`submit_1.txt` scripts to correspond to your files, you can run this
first step using:

``` sh
source submit_1.txt
```

**NOTE:** if you are planning to phase your data, the time to do this is
after generating VCF files, and before generating input in the next
step. Alternatively, you can generate VCF files however you’d like and
then jump into the MSMC process with `Step 2`.

### Step 1.5 - Phasing

If you are planning to phase your data, this is the time\! I won’t be
explaining how to do it, but common phasing programs to use are
[Beagle](https://faculty.washington.edu/browning/beagle/beagle.html)
([Browning &
Browning 2007](https://www.sciencedirect.com/science/article/pii/S0002929707638828)),
[fastPHASE](http://stephenslab.uchicago.edu/software.html#fastphase)
([Scheet &
Stephens 2006](https://www.sciencedirect.com/science/article/pii/S000292970763701X)),
and
[SHAPEIT2](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)
([Delaneau et al. 2013](https://doi.org/10.1016/j.ajhg.2013.09.002)).

## Step 2 - Generate Input

Now that you have `vcf` files and `mask` files for each individual and
each chromosome/scaffold, you can use these to generate the input for
MSMC2 using `msmc_2_generateInput_singleInd.sh`. This script is for
analyzing a single genome/individual separately; if you have multiple
individuals per population/species that you want to combine in analyses,
then follow the instructions in “MSMC on Multiple Individuals" below.
**NOTE**: I would recommend doing a first-pass analysis with all
individuals separate to make sure that they’re all reasonably similar,
and then combining multiple individuals in a later run for greater
statistical power.

The workhorse of this script is the line

``` sh
generate_multihetsep.py --mask=${MASK_INDIV} --mask=${MASK_GENOME} $VCF > $MSMC_INPUT
```

The output of this script, the MSMC input, has a list of heterozygous
sites in the genome, and the distance between these heterozygous sites.

## Step 3 - Run MSMC

Now that we have input files for all chromosomes for an individual, we
can run MSMC on these. You will need to edit the
`msmc_3_runMSMC_slurm.sh` or `msmc_3_runMSMC.sh` script to work with
your data and the parameters you’d like to use for running the analysis.
The things you might want to change are `P_PAR` (the timesteps used for
binning the data for anlaysis) and `-i` (number of iterations). I
recommend going with the default `P_PAR` to start, and then potentially
changing it later depending on how your data look. The scripts for
running MSMC can be called as follows:

``` sh
# for a single individual
sbatch msmc_3_runMSMC_slurm.sh 1 indID
# or without using slurm
./msmc_3_runMSMC.sh 1 indID
```

The outputs will be stored as `*.loop.text`, `*.log`, and `*.final.txt`.
This last file can be imported in to R to plot the results\! You’ll want
to have the mutation rate and generation time ready so that you can
scale Ne and time to individuals and years, respectively, like so:

``` R
### R Script for plotting MSMC Results ###
data <- msmc_output.final.txt

mu <- [mutation rate]
gen <- [generation time]

time <- (data$left_time_boundary/mu*gen)
pop.size <- (1/data$lambda)/(2*mu)

plot(time, pop.size, type="s", 
    xlab="log Years before present", 
    ylab="Effective Population Size", 
    log="x")
##########################################
```

From here, you may need to adjust the time steps and iterations for your
MSMC run. There are more details on what these mean, and why you might
want to adjust them, in the [MSMC
guide](https://github.com/stschiff/msmc/blob/master/guide.md).

### MSMC on Multiple Genomes

To combine multiple individuals, the steps generally follow what is
detailed above. However, in step 2, you will instead use the script
`msmc_2_generateInput_multiInd.sh`, which will combine all of the
individuals from a given population. For this script, which is submitted
using `submit_2_multi.txt`, you will need a file with the names of the
individuals you want combined in your MSMC run.

In step 3, you will then need to change `-I`, which specifies which
haplotypes to analyze. For one individual, this is `0,1`, for two
individuals `0,1,2,3`, and so on. This will be done automatically in the
“multiInd" script, which can simply be run as follows:

``` sh
# for multiple individuals in a population
sbatch msmc_3_runMSMC_slurm.sh numInd popID
# or without using slurm
./msmc_3_runMSMC.sh numInd popID
```

### Cross-Coalescence

If you’re interested in looking at the split between two populations or
species, then you’ll want to do things slightly differently, and will
need to actually make three MSMC runs. You should generate one combined
input file with all of the individuals you want to include in analyses.
Then, you’ll run MSMC three times, using different indices: once for
each of the two populations separately, and once for the two populations
combined. The command for running MSMC on two individuals from each of
two populations would look something like this:

``` sh
msmc2 -I 0,1,2,3 -o within_pop1 `cat INPUT_LIST.txt`
msmc2 -I 4,5,6,7 -o within_pop2 `cat INPUT_LIST.txt`
msmc2 -P 0,0,0,0,1,1,1,1 -o between_pop1-2 `cat INPUT_LIST.txt`
```

You then need to combine all three of these runs into a combined msmc
output file using `combineCrossCoal.py` (from `msmc-tools-master`):

``` sh
python combineCrossCoal.py \
    between_pop1-2.final.txt \
    within_pop1.final.txt \
    within_pop2.final.txt > combined_pop1-2.final.txt
```

This output can then be plotted in `R`, using the same transformations
as described above.

## Step 4 - Bootstraps

Once your MSMC plots are looking reasonable, then you can create
bootstraps to determine your confidence in your findings. This requires
randomly subsampling your data and running MSMC on those subsampled
data. To do this, you’ll use the `generate_bootstrap.sh` and
`msmc_bootstrap.sh` scripts. The first will generate the bootstrap data,
and the second will actually run the MSMC program on those subsampled
data.

## SCRIPTS

*Scripts originally written by Jelmer Poelstra, and modified by Rachel
Williams (July 2017) and subsequently by J. Rick (July 2018).*

### run\_snpable.sh

-----

``` sh
#!/bin/sh

#########################################
## change SBATCH commands as necessary ##
#########################################

#SBATCH --job-name SLURM_snpable
#SBATCH --account=passerinagenome
#SBATCH -o stdout_snpable
#SBATCH -e stderr_snpable
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrick@uwyo.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=2-00:00:00
#SBATCH --mem=0

date

echo "loading modules"

module load bwa
module load samtools
msmctools = /project/WagnerLab/genomic_resources/msmc_scripts/msmc-tools-master
PATH=$PATH:/project/passerinagenome/jrick/msmc/snpable/scripts
PATH=$PATH:msmctools

echo "modules loaded, read to go!"

ref=/path/to/ref.fa
k=35
out=ref_prefix
basedir=/path/to/folder

mkdir ${basedir}/snpable/
cd ${basedir}/snpable/

echo "Starting extraction of overlapping ${k}-mer subsequences"

splitfa $ref $k | split -l 20000000
cat x* >> ${out}_split.$k

echo "Aligning ${k}-mer reads to the genome with BWA, then converting to sam file"

# the genome needs to be indexed prior to this step-- if it has not already been indexed, run:
echo "indexing $ref"
bwa index $ref

echo "aligning reads to genome with BWA and converting to sam"
bwa aln -t 8 -R 1000000 -O 3 -E 3 ${ref} ${out}_split.${k} > ${out}_split.${k}.sai
bwa samse -f ${out}_split.${k}.sam $ref ${out}_split.${k}.sai ${out}_split.${k}

echo "reads aligned, starting to generate rawMask"
gen_raw_mask.pl ${out}_split.${k}.sam > ${out}_rawMask.${k}.fa

echo "raw mask created as ${out}_rawMask.35.fa, now generating final mask with stringency 50%"
gen_mask -l ${k} -r 0.5 ${out}_rawMask.${k}.fa > ${out}_mask.${k}.50.fa

##################################################################################
## note: it may be easiest to just edit the mappability script by hand          ##
## and then run it, rather than going through this process inside this script   ##
##################################################################################

echo "final mask saved as ${out}_mask.${k}.50.fa; now creating mappability mask"
grep 'with open' ${msmctools}/makeMappabilityMask.py | 
	sed 's/.*/with open\("${out}\_mask\.${k}.50.fa", "r"\) as f:/' 
    > ${msmctools}/${ref}.makeMappability.py
grep 'mask =' ${msmctools}/${ref}.makeMappability.py | 
	sed 's/.*/mask = MaskGenerator\("${basedir}\/mask\/${ref}\_\{\}\.mask\.35\.50\.bed\.gz".format(chr), chr)/'
    > ${basedir}/${ref}.makeMappability.py
python ${ref}.makeMappability.py

echo "done!"
date
```

### submit\_1.txt

-----

``` sh
for i in `cat INDS.txt`; \
do echo $i; \
IND=$i; \
sbatch --account=passerinagenome \
--job-name=aln_${i} \
--mail-type=END \
--output=outs/stdout_${i}_msmc1 \
--error=outs/stderr_${i}_msmc1 \
--nodes=1 \
--ntasks-per-node=8 \
--time=0-04:00:00 \
msmc_1_call.sh `echo unphased` $i; \
done
```

### msmc\_1\_call.sh

-----

``` sh
#!/bin/sh -l

## VARIABLES: 
## these are read in from submit_1.txt
PHASING=$1
IND=$2 

date

BASEDIR=/path/to/msmcfolder
GENOME=/path/to/reference.fa # reference genome fasta
BAMDIR=/path/to/bamdir
MSMCTOOLS=/project/WagnerLab/genomic_resources/msmc_scripts/msmc-tools-master
BAMFILE=`ls ${BAMDIR}/*_${IND}_*.bowtie2.sorted.bam`

module load samtools
module load bcftools
module load python/3.4.3
PATH=$PATH:$MSMCTOOLS

echo "modules loaded:"
module list

printf "\n \n \n \n"
date
echo "Current script: msmc_1_call.sh"
echo "Individual: $IND"
echo "Bamfile: $BAMFILE"
echo "Loaded modules:"
module list

if [ -f $GENOME ]
        then
                echo "Genome does not need to be unzipped, moving on!"
        else
                echo "Genome needs to be unzipped, unzipping now"
                gunzip -c ${GENOME}.gz > ${GENOME}
fi

if [ -f ${BAMFILE}.csi ]
        then
                echo "Bamfile index already exists, moving on!"
        else
                echo "Bamfile index does not exist, indexing bamfile"
                samtools index -c $BAMFILE
fi

echo "indexing genome"
samtools faidx $GENOME

### Calculate mean coverage (to be used as input for bamCaller.py):
echo 'ind.scaffold meancov' > ${BASEDIR}/coverage/coverage_samtoolsDepth_${IND}.txt
for SCAFFOLD in `cat SCAFFOLDS.txt`;
        do echo "working with $SCAFFOLD"
        MEANCOV=`samtools depth -r ${SCAFFOLD} ${BAMFILE} | 
        	awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | 
            tr ',' '.'` 
        echo "${IND}.${SCAFFOLD} $MEANCOV" >> 
        	${BASEDIR}/coverage/coverage_samtoolsDepth_${IND}.txt
        echo "Mean coverage for this individual, scaffold ${SCAFFOLD}: $MEANCOV"

### Generate a single-sample VCF and a mask-file:
	mkdir ${BASEDIR}/mask
	mkdir ${BASEDIR}/vcf
        
        MASK_IND=${BASEDIR}/mask/ind_mask.${IND}.${SCAFFOLD}.bed.gz # mask file to create
        VCF=${BASEDIR}/vcf/${IND}.${SCAFFOLD}.msmc.vcf # VCF file to create

        echo "starting samtools mpileup"

        samtools mpileup -q 20 -Q 20 -C 50 -u -r ${SCAFFOLD} -f $GENOME ${BAMFILE} | 
        	bcftools call -c -V indels | 
            bamCaller.py $MEANCOV ${MASK_IND} > ${VCF}
        sed 's/^\s#/#/' ${VCF} | gzip -c > ${VCF}.gz

## these values can be changed to filter differently
## -q = min. mapping qual
## -Q = min. base qual
## -C = coefficient for downgrading mapping qual for reads w/ excessive mismatches
## -u = generate uncompressed VCF/BCF
## -r = only specified region
## -f = fasta.
# bcftools: "call" will call SNPs/indels
## "-V indels" skips indels
## -c = consensus caller.

### Report:
        echo "Filtered VCF and mask created for $IND, scaffold $SCAFFOLD."
        done

echo "Done with script."

date
```

### submit\_2.txt

-----

``` sh
for i in `cat INDS.txt`; \
do echo $i; \
sbatch --account=passerinagenome \
--output=outs/stdout_${i}_msmc2 \
--error=outs/stderr_${i}_msmc2 \
--nodes=1 \
--ntasks-per-node=1 \
--time=0-4:00:00 \
--mail-type=END \
msmc_2_generateInput_singleInd.sh `echo unphased` $i; \
done
```

### msmc\_2\_generateInput\_singleInd.sh

-----

``` sh
#!/bin/bash -l

date

### Variables, read from submit_2.txt:
PHASING=$1
IND=$2

METHOD=samtools

BASEDIR=/path/to/msmcfolder
MSMCTOOLS=/project/WagnerLab/genomic_resources/msmc_scripts/msmc-tools-master

module load samtools
module load python/3.4.3
PATH=$PATH:$MSMCTOOLS

module list

for SCAFFOLD in `cat SCAFFOLDS.txt`; \
        do echo $SCAFFOLD; \
        VCF=`ls ${BASEDIR}/vcf/${IND}.${SCAFFOLD}.msmc.vcf.gz`; \

        MSMC_INPUT=${BASEDIR}/input/msmc_input.${IND}.${SCAFFOLD}.txt

        printf "\n \n \n \n"
        date
        echo "Script: msmc_2_generateInput_singleInd"
        echo "Individual: ${IND}"
        echo "Scaffold: ${SCAFFOLD}"
        echo "Phasing: ${PHASING}"
        echo "Method: ${METHOD}"
        echo "MSMC input file: ${MSMC_INPUT}"
        echo "VCF: ${VCF}"

        if [ -f $VCF ]
        	then
            echo "VCF exists, moving on!"
			### Generate MSMC input files:
            if [ $METHOD == samtools ]
                  then
                  MASK_INDIV=${BASEDIR}/mask/ind_mask.${IND}.${SCAFFOLD}.bed.gz
                  MASK_GENOME=${BASEDIR}/mask/pamoena_${SCAFFOLD}.mask.35.50.bed.gz
                  echo "MASK: ${MASK_INDIV}"
                  echo "MAPPABILITY MASK: ${MASK_GENOME}"
                  echo "Creating MSMC input file WITH individual mask (samtools)"
                  generate_multihetsep.py --mask=${MASK_INDIV} --mask=${MASK_GENOME} $VCF >
                   	$MSMC_INPUT

            elif [ $METHOD == gatk ]
                  then
                  echo "Creating MSMC input file WITHOUT individual mask (gatk)"
                  ${MSMCTOOLS}/generate_multihetsep.py $VCF > $MSMC_INPUT
            fi

        else
                echo "VCF does not exist, moving on to next scaffold"
        fi

   echo "Done with $SCAFFOLD; moving on to next scaffold"
done

echo "Done with script."
date
exit 0

####
```

### msmc\_3\_runMSMC\_onepop.sh

-----

``` sh
#!/bin/sh

#########################################
## change SBATCH commands as necessary ##
#########################################

#SBATCH --job-name=msmc3
#SBATCH --account=passerinagenome
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --output=outs/stdout_msmc3
#SBATCH --error=outs/stderr_msmc3
#SBATCH --time=0-01:00:00
#SBATCH --mail-type=ALL

BASEDIR=/path/to/msmcfolder

### Variables:
#Had to put loop in script because it wouldn't carry the values over
for i in `cat INDS.txt`
        do echo $i
        IND=$i

P_PAR=1*2+25*1+1*2+1*3
NR_IND=1
POP_OR_IND=$IND
RUN_NAME=msmc_${IND}
METHOD="samtools + mask"

find ${basedir}/input/msmc_input.${IND}.*.txt -size 0 -delete
ls ${basedir}/input/msmc_input.${IND}.*.txt > SCAFS_INPUT_${IND}

module load msmc2

### Report settings/parameters:
date
echo "Script: msmc_3_runMSMC_onepop.sh"
echo "Run name: $RUN_NAME"
echo "SNP calling method: $METHOD"
echo "Period setting: $P_PAR"
echo "Nr of individuals (1 or 2+): $NR_IND"
echo "Population or individuals ID: $POP_OR_IND"
echo "Individual: $IND"
echo "Scaffolds: `cat SCAFS_INPUT_${IND}`"
echo "Iterations: 100"


if [ $NR_IND == 1 ]
        then
        echo "Running MSMC for one individual"
        MSMC_INPUT=`cat SCAFS_INPUT_${IND}`
        MSMC_OUTPUT=${BASEDIR}/output/msmc_output.${RUN_NAME}

        if [ -f $MSMC_INPUT ]
                then
                        echo "MSMC_INPUTS: $MSMC_INPUT"
                        echo "MSMC_OUTPUT: $MSMC_OUTPUT"
                else
                        echo "MSMC_INPUT does not exist! Exiting now"
                        exit 1
        fi

        msmc2 -t 8 -p $P_PAR -i 100 -o $MSMC_OUTPUT -I 0,1 $MSMC_INPUT

        mv $MSMC_OUTPUT*loop.txt ${BASEDIR}/output/log_and_loop/
        mv $MSMC_OUTPUT*log ${BASEDIR}/output/log_and_loop/
else
        echo "Running MSMC for $NR_IND individuals"
        MSMC_INPUT=`cat SCAFS_INPUT`
        MSMC_OUTPUT=${BASEDIR}/output/msmc_output.${POP_OR_IND}.${RUN_NAME}

        msmc2 -t 8 -p $P_PAR -i 100 -o ${MSMC_OUTPUT} -I 0,1 $MSMC_INPUT
fi

echo "DONE_WITH_SCRIPT"
date
```

### generate\_bootstrap.sh

-----

``` sh

#!/bin/bash -l

for i in `cat INDS.txt`; \
do echo $i;\
IND=`echo ${i[@]}`; \
SCAFS=`cat SCAFFOLDS.txt`; \
done

BASEDIR=/path/to/msmcfolder
RUN_NAME=msmc_bootstrap
POP_OR_IND=$i
NUM_CHR=`cat SCAFFOLDS.txt | wc -l`
NUM_BS=50 # number of bootstraps to perform

PATH=$PATH:/project/WagnerLab/genomic_resources/msmc_scripts/msmc-tools-master

#input for the bootstrapping
BS_INPUT=`for c in $SCAFS; 
	do find ${BASEDIR}/input -maxdepth 1 -name "*${POP_OR_IND}.${c}*.txt"; done`

#output from the bootstrapping 
BS_OUTPUT=${BASEDIR}/bootstrap/${POP_OR_IND}.${RUN_NAME}.bootstrap

multihetsep_bootstrap.py -n $NUM_BS \
    -s 20000000 \
    --chunks_per_chromosome 10 \
    --nr_chromosomes \
    $NUM_CHR \
    $BS_OUTPUT \ 
    $BS_INPUT

cd ${BASEDIR}/bootstrap

ls -d *${POP_OR_IND}*bootstrap_* > ${BASEDIR}/bootstrap/${POP_OR_IND}_bs_file_list.txt
```

### submit\_bootstrap.sh

-----

``` sh
BASEDIR=/path/to/msmcfolder
for i in `cat INDS`; do echo $i; \
for q in `cat ${BASEDIR}/bootstrap/${i}_bs_file_list.txt`; \
do echo $q;
sbatch --account=WagnerLab \
--nodes=1 --ntasks-per-node=4 \
--time=0-01:00:00 \
--output=/project/WagnerLab/jrick/msmc_Sept2017/outs/stdout_bootstrap_$q \
--error=/project/WagnerLab/jrick/msmc_Sept2017/outs/stderr_bootstrap_$q \
msmc_bootstrap.sh $q $i; \
done; done
```

### msmc\_bootstrap.sh

-----

``` sh

#!/bin/bash -l

FOLDER=$1
POP_OR_IND=$2

BASEDIR=/path/to/msmcfolder
RUN_NAME=msmc_bootstrap
P_PAR=1*2+25*1+1*2+1*3 ## this needs to match the original MSMC run

BS_LIST=${BASEDIR}/bootstrap/${POP_OR_IND}_bs_file_list.txt

module load msmc2

##### run msmc alone #### 
MSMC_BS=`ls ${BASEDIR}/bootstrap/${FOLDER}/*.txt`
MSMC_OUTPUT=${BASEDIR}/bootstrap/msmc_output.${RUN_NAME}.${FOLDER}

msmc2 -t 20 -p $P_PAR -o $MSMC_OUTPUT -I 0,1 $MSMC_BS
```
