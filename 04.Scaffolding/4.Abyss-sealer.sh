#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=126675
#SBATCH --job-name="AbyssSealer"
#SBATCH --output=AbyssSealer_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END

	This script takes as input one argument. It runs from the master folder of Abyss Sealer and then we change directory to the directory of each specific Strain.	
	For Strain 01 the correct usage is: 
	$0 Strain01

	NOTE:
	We are using the pre-filtered reads here! These are the ones right after fastp and before prinseq. This is why trim-quality is required as an option.

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
   echo "Error: Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9!" >&2
   generalInfo >&2
   exit 1
fi


# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
StrainCode="HK7V7BBXY"

OutputDir=${HOME}/Titlos_ktisis/Abyss/${StrainX}

AbyssDir=${HOME}/Software/abyss-2.3.1
KrakenDir=${HOME}/Titlos_ktisis/Kraken/${StrainX}
ReadDir=${HOME}/Titlos_ktisis/Fastp/${StrainX}
Read1=${StrX}_${StrainCode}_R1_fastp.fastq
Read2=${StrX}_${StrainCode}_R2_fastp.fastq
ContigKraken=${StrX}_CLA_Kraken.fasta

#---------------------------------------------------------------------------------------------------------------------------------

export LC_ALL=en_US.UTF-8

[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}
cd ${OutputDir}


${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=25 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k25.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=39 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k39.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=49 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k49.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=59 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k59.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=69 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k69.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=81 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k81.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=91 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k91.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=101 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k101.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=117 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k117.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Sealer/abyss-sealer -vv -k25 -k39 -k49 -k59 -k69 -k81 -k91 -k101 -k117 -o ${StrX}_CLA_Sealer -S ${KrakenDir}/${ContigKraken} -i k25.bloom -i k39.bloom -i k49.bloom -i k59.bloom -i k69.bloom -i k81.bloom -i k91.bloom -i k101.bloom -i k117.bloom --max-gap-length=1000 --max-branches=3000 --max-paths=500 --fix-errors --threads=$SLURM_NTASKS --trim-quality=16 


#rm *.bloom


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


echo "============================================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
elapsed=$SECONDS
echo Time taken: `printf '%dd %dh:%dm:%ds\n' $((elapsed/86400)) $((elapsed%86400/3600)) $((elapsed%3600/60)) $((elapsed%60))`

exit 0
