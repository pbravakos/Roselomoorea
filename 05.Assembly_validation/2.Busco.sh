#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=126675
# #SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Busco"
#SBATCH --output=Busco_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END
	
	This script takes as input one argument. 
	For Strain 01 that would be: 
	./$0 Strain01
	This script runs from the master folder of Busco and then we change directory to the folder of each specific Strain folder.
	
	NOTE:
	Download relative busco database from:
	https://busco-data.ezlab.org/v4/data/lineages/ 
	
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
if [[ $# -ne 1  || ! $1 =~ ^Strain[0-9]{2}$ ]]; then
   echo "Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9!" >&2
   generalInfo >&2
   exit 1
fi


# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
StrainCode="HK7V7BBXY"
OutputDir=${HOME}/Titlos_ktisis/Busco

ContigPilonDir=${HOME}/Titlos_ktisis/Pilon/${StrainX}
ContigSpadesDir=${HOME}/Titlos_ktisis/Spades/${StrainX}/Iter2
BuscoDir=${HOME}/Software/busco/bin
BuscoGammaDir=${HOME}/Software/busco/bacillales_odb10

BuscoPilonOut=${StrX}_Pilon_Busco
BuscoSpadesOut=${StrX}_Spades_Busco
ContigSpades=scaffolds.fasta
ContigPilon=${StrX}_Pilon_CLA.fasta

# AUGUSTUS_CONFIG_PATH=${HOME}/Software/augustus-3.3.3/config
declare -x AUGUSTUS_CONFIG_PATH="${HOME}/Software/augustus-3.3.3/config"
export BUSCO_CONFIG_FILE="${HOME}/Software/busco/config/config.ini"
export AUGUSTUS_CONFIG_PATH
export PATH="${HOME}/Software/augustus-3.3.3/scripts:$PATH"
export PATH="${HOME}/Software/augustus-3.3.3/bin:$PATH"
export PATH="${HOME}/Software/hmmer-3.3.2/bin:$PATH"
export PATH="${HOME}/Software/ncbi-blast-2.10.1+/bin:$PATH"
export PATH="${HOME}/Software/busco/scripts:$PATH"
export LC_ALL=en_US.UTF-8

#-----------------------------------------------------------------------------------------------------

[[ ! -d ${OutputDir} ]] && mkdir ${OutputDir}
cd ${OutputDir}

## Run busco 2 times. One for the Pilon final scaffolds and one for the Spades assembler output scaffolds.
#ln -s ${ContigPilonDir}/${ContigPilon}

#python ${BuscoDir}/busco --out ${BuscoPilonOut} \
#				--in ${ContigPilonDir} \
#				--lineage ${BuscoGammaDir} \
#				--mode genome \
#				--cpu $SLURM_NTASKS \
#				--long \
#				--force
#				# --tmp ./tmp
#				# --restart  # Continue a run that had already partially completed.
#				# --config ${HOME}/Software/busco/config/config.ini \

# Run again!

ln -s ${ContigSpadesDir}/${ContigSpades}
python ${BuscoDir}/busco --out ${BuscoSpadesOut} \
				--in ${ContigSpades} \
				--lineage ${BuscoGammaDir} \
				--mode genome \
				--cpu $SLURM_NTASKS \
				--long \
				--force
				# --tmp ./tmp
				# --restart  # Continue a run that had already partially completed.
				# --config ${HOME}/Software/busco/config/config.ini \




# Next steps
# Compare the two outcomes.



#{
## Finally create a job dependency to move this job's output to working directory.
#cd $SLURM_SUBMIT_DIR
#MoveOutput=${OSDStation}_mv_output.sh       

## Create the bash file, to move the SLURM output.
#cat > ${MoveOutput} << EOF
##!/bin/bash
##SBATCH --partition=fast
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --mem-per-cpu=100
##SBATCH --job-name="mv"

#mv ${SLURM_JOB_NAME}_job_${SLURM_JOB_ID}.out ${OutputDir}

#EOF

#echo
## Start a job dependency to move the Sdtout file to the output directory.
#sbatch --dependency=afterany:"$SLURM_JOB_ID" ${MoveOutput}

#rm ${MoveOutput} 
#}

echo "======================================================================"
echo "SLRUM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
elapsed=$SECONDS
echo Time taken: `printf '%dd %dh:%dm:%ds\n' $((elapsed/86400)) $((elapsed%86400/3600)) $((elapsed%3600/60)) $((elapsed%60))`

exit 0
