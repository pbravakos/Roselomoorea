#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=126675
# #SBATCH --mem-per-cpu=6400
#SBATCH --job-name="KaijuDB"
#SBATCH --output=KaijuDB_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END
	
	This script takes as input one argument. 
	For Strain 01 that would be: 
	sbatch ${ScriptName} Strain01
        
        Create a Refseq Databse for Kaiju.
        This task requires a lot of memory.

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

# Select where to build the database.
OutputDir=${HOME}/Software/kaiju/kaijudb



export LC_ALL=en_US.UTF-8
export PATH="${HOME}/Software/kaiju/bin:$PATH"

kaiju -h


# Print the Prokka version we are going to use in the analysis.
prokka --version
echo
which prokka
echo
#-------------------------------------------------------------------------------------------------
	[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}
	cd ${OutputDir}
	kaiju-makedb -s refseq -t $SLURM_NTASKS -v
	
	#rm -r *.bwt *.sa genome
	# Kaiju only needs the files refseq/kaiju_db_refseq.fmi, nodes.dmp, and names.dmp
	# The remaining files can be deleted.
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
