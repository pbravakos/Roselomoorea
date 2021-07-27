#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=126675
# #SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Blastp"
#SBATCH --output=Blastp_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END

	This script takes as input one argument.
	For Strain 01 the correct usage is: 
	sbatch ${ScriptName} Strain01

	We are running Blastp on the amino acid fasta file that was created by Interproscan.

	NOTE
	We are going to use these results to run Blast2GO.

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
if [[ $# -ne 1 || ! $1 =~ ^Strain[0-9]{2}$ ]]; then
   echo "Error: Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9!" >&2
   generalInfo >&2
   exit 1
fi


# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
StrainCode="HK7V7BBXY"
Assembler="Pilon"  # Choose either "Pilon" or "Spades"
Prefix=${StrX}_${Assembler}

OutputDir=${HOME}/Titlos_ktisis/BlastP/${Prefix}

AnvioDir=${HOME}/Titlos_ktisis/Anvio/Str28_${Assembler}

# Blast parameters that can be changed
OutFile=${Prefix}_blastp_Anvio.xml
QueryFasta=${StrX}_InterProScan_amino-acid-sequences.fasta
DataBase="nr"
Task=blastp
WordSize=7  # Has to be less than 8 for blastp.
Eval=1e-15
OutFormat="16" 
           
MaxSeqs=10  # Maximum number of aligned sequences to keep. Not applicable for outfmt <= 4 
MaxDesc=20
MaxAln=15



# Export the Blastn binary to the PATH.
export PATH="${HOME}/Software/ncbi-blast-2.10.1+/bin:$PATH"

# Export the local database directory.
export BLASTDB="/mnt/dbs/ncbi/nr"


#--------------------------SANITY CHECKS -------------------------------------------------------------------#
# Check for the existence of the CLA fasta output. 
if [[ ! -s ${AnvioDir}/${QueryFasta} ]]
then
    echo "${QueryFasta} cannot be found in ${CLADir}. Please run CLA first." >&2
    echo >&2
    generalInfo >&2
    exit 1
fi

#--------------------------------------------------------------------------------------------------------------------

# Print the blastp version that is going to be used for the subsequent analyses.
echo
blastp -version
echo

[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}
cd ${OutputDir}

# Create a link of the InterproScan output protein fasta file.
if [[ ! -s ${QueryFasta} ]]; then
    ln -s ${AnvioDir}/${QueryFasta}
fi

blastp -query ${QueryFasta} \
	-db "$DataBase" \
	-out ${OutFile} \
	-task $Task \
	-word_size $WordSize \
	-evalue $Eval \
	-outfmt "$OutFormat" \
	-num_threads $SLURM_NTASKS \
	-max_target_seqs $MaxSeqs


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
