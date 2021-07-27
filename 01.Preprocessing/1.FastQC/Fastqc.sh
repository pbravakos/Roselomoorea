#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=6400
# #SBATCH --mem=128000
#SBATCH --job-name="FastQC"
#SBATCH --output=FastQC_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END

	This script takes as input one argument. 
	For Strain 01 that would be: 
	script.sh Strain01
	This script runs from the master folder of FastQC and then we change directory to the folder of each specific Strain folder.
	We assume that in the master folder there are already folders named Strain01 Strain02 etc.

	Here we take the raw fastq files, taken from the illumina machine, and check their quality.

	We have already created soft links of each raw fastq paired reads inside the correct directory

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
if [[ $# -ne 1 ]] || [[ ! $1 =~ Strain[0-9]{2} ]]; then
    echo "Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9'!" >&2
    generalInfo >&2
    exit 1
fi


echo

echo
# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}

OutputDir=${HOME}/Titlos_ktisis/FastQC/${StrainX}
RawReadsDir=${HOME}/Titlos_ktisis/Fastq/${StrainX}
StrainCode="HK7V7BBXY"
Pair1Fastq=${StrX}_R1_${StrainCode}.fastq.gz
Pair2Fastq=${StrX}_R2_${StrainCode}.fastq.gz

export PATH="${HOME}/Software/FastQC:$PATH"
export LC_ALL=en_US.UTF-8

[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}

cd ${OutputDir}
echo
mkdir tmp
echo
fastqc --version
echo


fastqc --threads $SLURM_NTASKS --noextract --dir tmp -o ${OutputDir} ${RawReadsDir}/${Pair1Fastq} ${RawReadsDir}/${Pair2Fastq} 

<< ////
	
	NOTE:
	-t --threads    Specifies the number of files which can be processed
                        simultaneously.  Each thread will be allocated 250MB of
                        memory so you shouldn't run more threads than your
                        available memory will cope with, and not more than
                        6 threads on a 32 bit machine

////


#-------------------------------------------------------------------------------
# Finally create a job dependency to move this job's output to working directory.
cd $SLURM_SUBMIT_DIR
MoveOutput=${StrainX}_job_${SLURM_JOBID}_mv_output.sh       

# Create the bash file, to move the SLURM output.
cat > ${MoveOutput} << EOF
#!/bin/bash
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100
#SBATCH --job-name="mv"
#SBATCH --output=/dev/null

mv ${SLURM_JOB_NAME}_job_${SLURM_JOB_ID}.out ${OutputDir}

EOF

echo
# Start a job dependency to move the Sdtout file to the output directory.
sbatch --dependency=afterany:"$SLURM_JOB_ID" ${MoveOutput}

rm ${MoveOutput}


echo
echo "==========================================================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
elapsed=$SECONDS
echo Time taken: `printf '%dd %dh:%dm:%ds\n' $((elapsed/86400)) $((elapsed%86400/3600)) $((elapsed%3600/60)) $((elapsed%60))`

exit 0
