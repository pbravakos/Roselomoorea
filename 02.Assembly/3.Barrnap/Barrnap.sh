#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=126675
#SBATCH --job-name="Barrnap"
#SBATCH --output=Barrnap_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END
	
	This script takes as input one argument. 
	For Strain 01 that would be: 
	sbatch ${ScriptName} Strain01
	This script runs from the master folder of Barrnap and then we change directory to the folder of each specific Strain folder.
	As input we use here the Spades output, contig fasta files.

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
   echo "Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9!" >&2
   generalInfo >&2
   exit 1
fi


# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
StrainCode="HK7V7BBXY"
Iteration=Iter2  # This parameter can be changed every time we want to have more than one Spades runs saved on separate directories. 

OutputDir=${HOME}/Titlos_ktisis/Barrnap/${StrainX}

BarrnapDir=${HOME}/Software/barrnap/bin

InputFasta=${HOME}/Titlos_ktisis/Spades/${StrainX}/${Iteration}/scaffolds.fasta


export PATH="${HOME}/Software/barrnap-0.9/bin:$PATH"
export LC_ALL=en_US.UTF-8


echo
echo "Barrnap version:"
barrnap --version
echo
echo


[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}

cd ${OutputDir}
barrnap --threads $SLURM_NTASKS --kingdom 'bac' --outseq ${StrX}_${StrainCode}_ribosomal_seqs_ONLY.fasta $InputFasta


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




<< ////
	FURTHER STEPS NEEDED:
	Blast search the barrnap results. 
	Download the most relevant genomes and use them as references from now on.
	Save the Downloaded NCBI genomes in the References directory on the server!
	Reference sequences can be checked with Bandage if blasted against the scaffolds from Spades output. From these blast results and by coloring the blast hits in Bandage we can see whether the Reference genome is actually a reference i.e. if it has blast hits with most of the scaffolds or not.
////

echo "==============================="
echo "SLRUM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
elapsed=$SECONDS
echo Time taken: `printf '%dd %dh:%dm:%ds\n' $((elapsed/86400)) $((elapsed%86400/3600)) $((elapsed%3600/60)) $((elapsed%60))`

exit 0
