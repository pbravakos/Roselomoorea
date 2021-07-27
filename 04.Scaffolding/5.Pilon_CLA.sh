#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=126675
# #SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Pilon"
#SBATCH --output=Pilon_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END
	
	This script takes as input one argument. 
	For Strain 01 that would be: 
	sbatch ${ScriptName} Strain01
	This script runs from the master folder of Pilon and then we change directory to the folder of each specific Strain folder.
	
	NOTE:
	Here we use the fastq files after Prinseq even if further downstream filtering has been done (e.g. contamination removal) because we want as much information as possible to correct our contigs!!

END
}

if [[ -n $SLURM_JOB_ID ]];  then
    ScriptName=$(scontrol show job $SLURM_JOBID | awk '/Command=/ {print $1}' | awk -F '[ =]' '{print $2}' | grep -Eo "[^/]+$")
    # Some job specific info
    echo "Job ID is = " $SLURM_JOBID
    echo "SLURM cluster name = " $SLURM_CLUSTER_NAME
    echo "SLURM partition = " $SLURM_JOB_PARTITION
    echo "SLURM node list = " $SLURM_JOB_NODELIST
    echo "SLURM num of nodes = " $SLURM_JOB_NUM_NODES
    echo "SLURM number of tasks = " $SLURM_NTASKS
    echo "SLURM memory per node = " $SLURM_MEM_PER_NODE
    echo "SLURM memory per cpu = " $SLURM_MEM_PER_CPU
    echo "working directory = " $SLURM_SUBMIT_DIR
    echo "=================================================="
    echo "SBATCΗ job started " `date`
    echo "=================================================="
    echo
else 
    ScriptName=${0##*/}
    generalInfo
    exit 1
fi

# Check that an argument has been given in the correct form.
if [[ $# -ne 1 ]] || [[ ! $1 =~ ^Strain[0-9]{2}$ ]]; then
   echo "Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9!" >&2
   generalInfo >&2
   exit 1
fi


# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
StrainCode="HK7V7BBXY"
OutputDir=${HOME}/Titlos_ktisis/Pilon/${StrainX}

PilonDir="${HOME}/Software/Pilon"
SealerDir="${HOME}/Titlos_ktisis/Abyss/${StrainX}"
PrinseqDir="${HOME}/Titlos_ktisis/Prinseq/${StrainX}"
ContigSealer=${StrX}_CLA_Sealer_scaffold.fa

PESorted=${StrX}_PE_Pilon_sorted.bam
SESorted=${StrX}_SE_Pilon_sorted.bam

PE1=${StrX}_${StrainCode}_prinseq_good_R_1.fastq
PE2=${StrX}_${StrainCode}_prinseq_good_R_2.fastq
SE=${StrX}_${StrainCode}_prinseq_good_singletons.fastq


# Parameters to Pilon which can be changed
FlankBases=1 # Controls how much of the well-aligned reads will be used; this many bases at each end of the good reads will be ignored (default 10).
GapMarginBases=2000   # Closed gaps must be within this number of bases of true size to be closed (100000)
Kmer=77    # Kmer size used by internal assembler (default 47)
#--------------------------------------------------------------------------------------------------------------------------

MemMB=$((${SLURM_MEM_PER_NODE:-$((SLURM_MEM_PER_CPU*SLURM_NTASKS))}))

export LC_ALL=en_US.UTF-8
export PATH="${HOME}/Software/bwa-mem2:$PATH"
export PATH="${HOME}/Software/samtools-1.10:$PATH"

[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}
cd ${OutputDir}

## Create soft link to contig fasta files. This is done because contig fasta files need to be in the same folder as the bwa.sh file for the pipeline to work!
#ln -s ${SealerDir}/${ContigSealer}

#echo
## Print Bwa version
#( bwa 3>&1 1>&2- 2>&3- ) | head -n 3
#echo

#echo
## Print samtools version
#( samtools 3>&1 1>&2- 2>&3- ) | head -n 3
#echo

## Create bwa index for Spades
#bwa-mem2 index ${ContigSealer}
## Align reads to contigs with bwa mem and bwa bwasw. bwasw seems not to work well with paired end reads. For this reason we use bwa mem.
##Turn sam to bam and finally sort by coordinates
#bwa-mem2 mem -t $SLURM_NTASKS ${ContigSealer} ${PrinseqDir}/${PE1} ${PrinseqDir}/${PE2} |\
#samtools view -hu -@ $SLURM_NTASKS - |\
#samtools sort -@ $SLURM_NTASKS - -o ${PESorted}
##Index the sorted bam file
#samtools index ${PESorted}

##Repeat the same pipeline for the single end reads.
##1st
#bwa-mem2 mem -t $SLURM_NTASKS ${ContigSealer} ${PrinseqDir}/${SE} |\
#samtools view -hu -@ $SLURM_NTASKS - |\
#samtools sort -@ $SLURM_NTASKS - -o ${SESorted}
##Index the sorted bam file
#samtools index ${SESorted}


##Remove files that are not needed for downstream analysis 
#rm *.fasta.*
echo
echo
echo
echo

java -Xmx${MemMB}M -jar ${PilonDir}/pilon-1.24.jar --genome ${ContigSealer} \
							--frags ${PESorted} \
							--unpaired ${SESorted} \s
							--output ${StrX}_Pilon_CLA \
							--outdir ${OutputDir} \
							--K $Kmer \
							--flank $FlankBases \
							--gapmargin $GapMarginBases \
							--fix "gaps","local","amb","breaks" \
							--iupac  \
							--verbose \
							--debug
							# --tracks
							#--threads $SLURM_NTASKS \
							
## Remove files not needed anymore
#rm *.bam* *.fa.*

## Rename the Pilon fasta scaffolds output to something better!
#mv ${StrX}_Pilon_.fasta ${StrX}_Pilon_CLA.fasta


##-------------------------------------------------------------------------------
## Finally create a job dependency to move this job's output to working directory.
#cd $SLURM_SUBMIT_DIR
#MoveOutput=${StrainX}_job_${SLURM_JOBID}_mv_output.sh       

## Create the bash file, to move the SLURM output.
#cat > ${MoveOutput} << EOF
##!/bin/bash
##SBATCH --partition=fast
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --mem-per-cpu=100
##SBATCH --job-name="mv"
##SBATCH --output=/dev/null

#mv ${SLURM_JOB_NAME}_job_${SLURM_JOB_ID}.out ${OutputDir}

#EOF

#echo
## Start a job dependency to move the Sdtout file to the output directory.
#sbatch --dependency=afterany:"$SLURM_JOB_ID" ${MoveOutput}

#rm ${MoveOutput}


echo "==============================="
echo "SLRUM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
elapsed=$SECONDS
echo Time taken: `printf '%dd %dh:%dm:%ds\n' $((elapsed/86400)) $((elapsed%86400/3600)) $((elapsed%3600/60)) $((elapsed%60))`

exit 0
