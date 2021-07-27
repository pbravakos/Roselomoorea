#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=126675
#SBATCH --job-name="Spades"
#SBATCH --output=Spades_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END

	This script takes as input one argument. It runs from the master folder of Spades and then we change directory to the directory of each specific Strain.	
	For Strain 01 the correct usage is: 
	sbatch ${ScriptName} Strain01

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
   echo "Error: Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9!" >&2
   generalInfo >&2
   exit 1
fi


# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
StrainCode="HK7V7BBXY"
Iteration=Iter1  # This parameter can be changed every time we want to have more than one Spades runs saved on separate directories. 

SpadesDir=${HOME}/Software/SPAdes-3.15.2-Linux/bin
OutputDir=${HOME}/Titlos_ktisis/Spades/${StrainX}/${Iteration}

Yaml=${StrainX}_reads.yaml
TempDir=tmp_${StrainX}

kmers=23,31,43,55,71,83,91,103,111,123 # list of k-mer sizes (must be odd and less than 128) [default: 'auto']

# Reads for the yaml file.
FastqDir=${HOME}/Titlos_ktisis/BBnorm/${StrainX}
Read1=${StrX}_${StrainCode}_R1_bbnorm.fastq
Read2=${StrX}_${StrainCode}_R2_bbnorm.fastq
#Merged=${StrX}_${StrainCode}_merged.fastq
Singletons=${StrX}_${StrainCode}_prinseq_good_singletons.fastq

MemGB=$((${SLURM_MEM_PER_NODE:-$((SLURM_MEM_PER_CPU*SLURM_NTASKS))}/1024)) # Memory in Gb. SPAdes terminates if it reaches this limit. The default value is 250 Gb.


export LC_ALL=en_US.UTF-8

#--------------------------SANITY CHECKS -------------------------------------------------------------------#

[[ -z ${MemGB} ]] && { echo "Memory requirement for Spades is not set">&2; exit 1; }


if [[ ! -s ${FastqDir}/${Read1} || ! -s ${FastqDir}/${Read2} ]]
then
    echo "Read files ${Read1} and/or ${Read2} cannot be found in ${FastqDir}" >&2
    generalInfo >&2
    exit 1 
fi

if [[  ! -s ${FastqDir}/${Singletons} ]]
then
    echo "Read files ${Singletons} cannot be found in ${FastqDir}" >&2
    generalInfo >&2
    exit 1 
fi

#-------------------------------------------------------------------------------------------------------------

cleanup()
{
  echo
  echo "Caught Exit Signal. Cleaning up."
  echo
  rm -rf ${TempDir} ${Yaml} input_dataset.yaml params.txt dataset.info run_spades.yaml \
  run_spades.sh K* before_rr.fasta first_pe_contigs.fasta misc pipeline_state 
  # assembly_graph.fastg
  # rm -rf scaffolds.paths contigs.paths
  echo "Done cleanup ... quitting."
  echo
}

echo
echo Memory used by Spades is ${MemGB} GB
echo


echo
echo "Spades version:"
python3 ${SpadesDir}/spades.py -v
echo
echo

[[ ! -d ${OutputDir} ]] && mkdir ${OutputDir}
cd ${OutputDir}


# Check whether there is already a previous run of Spades and if yes start from the checkpoint.
Checkpoint=$(find ${OutputDir} -name checkpoint.dat)
if [[ -z ${Checkpoint} ]]
then
# Create the ${Yaml} yaml file.
cat > ${Yaml} << EOF
- "left reads":
  - "${FastqDir}/${Read1}"
  "orientation": "fr"
  "right reads":
  - "${FastqDir}/${Read2}"
  "single reads":
  - "${FastqDir}/${Singletons}"
  "type": "paired-end"
EOF

    python3 ${SpadesDir}/spades.py 	--only-assembler \
    					-k ${kmers} \
    					--isolate \
    					--dataset ${Yaml} \
    					-t ${SLURM_NTASKS} \
    					--memory ${MemGB} \
    					-o ${OutputDir} \
    					--checkpoints "last" \
    					--tmp-dir ${TempDir}
    					# --careful \ Not used because it is not compatible with the --isolate option
    					#  \
    # We want to perform the cleanup only when the program runs to completion. Otherwise we want to be able to continue from the checkpoint.
    # Thus we run the cleanup only with exit code 0. No cleanup will be performed if the exit code is not 0.
    [[ $? -eq 0 ]] && cleanup
else
    echo 'Continue Spades'
    echo
    # Run the program
    python3 ${SpadesDir}/spades.py --continue -o ${OutputDir}
    [[ $? -eq 0 ]] && cleanup
fi



# NEXT STEPS:
# Open fastg output with Bandage to check for contamination, and either continue with scaffolding, or start the optional decontamination pipeline.



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
echo "SBATCH job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
elapsed=$SECONDS
echo Time taken: `printf '%dd %dh:%dm:%ds\n' $((elapsed/86400)) $((elapsed%86400/3600)) $((elapsed%3600/60)) $((elapsed%60))`

exit 0

