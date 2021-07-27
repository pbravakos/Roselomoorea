#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Kraken"
#SBATCH --output=Kraken_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END
	
	This script takes as input one argument. 
	For Strain 01 that would be: 
	./$0 Strain01
	This script runs from the master folder of Kraken and then we change directory to the folder of each specific Strain folder.

	NOTE:
	We are supposed to run this script after CLA, in order to filter out any scaffolds that is considered to be contamination.

	IMPORTANT:
	Please, read the report from the Kraken output and evaluate it. This is very important, because based on this report we try to filter out the scaffolds. Based on the report evaluation, make any adjustments to this script e.g. by changing the search pattern in the grep command.
	IMPORTANT:
	Here we take into account only the classified scaffolds by Kraken. Please check the Kraken output (StdOut) to check that ALL scaffolds have been classified. In case there are unclassified scaffolds, these should probably be included in the final output.

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

OutputDir=${HOME}/Titlos_ktisis/Kraken/${StrainX}

KrakenDir=/mnt/big/Metagenomics/kraken-2.0.7/installation
KrakenDB=/mnt/dbs/kraken2-db
KrakenClassified=Kraken_Classified_CLA_scaffolds.fasta
KrakenOutput=Kraken_CLA_scaffolds.out
KrakenReport=Kraken_CLA_scaffolds.report

CLADir=${HOME}/Titlos_ktisis/CLA/${StrainX}/CLA-Results

#-----------------------------------------------------------------------------------------------------

export LC_ALL=en_US.UTF-8

[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}
cd ${OutputDir}

# Create a link of the CLA output file.
ln -s ${CLADir}/final_Scaffolds.fa

# Run Kraken
${KrakenDir}/kraken2 --db ${KrakenDB} \
			--classified-out ${KrakenClassified} \
			--output ${KrakenOutput} \
			--use-names \
			--threads $SLURM_NTASKS \
			final_Scaffolds.fa

# Get the Kraken report. 
${KrakenDir}/kraken2 --report ${KrakenReport} --threads $SLURM_NTASKS --db ${KrakenDB} $KrakenClassified

# Extract the fasta headers that have been classified as Bacillus by Kraken.
grep -P -o "[Sa-z_0-9]*\t[A-z]*[Bb]acil" ${KrakenOutput} | grep -P -o "S[a-z_0-9]*" > Bacillus_fasta.headers

# Create the final fasta output. (We do this in order to make sure that it is empty before going to the while loop in the next step.)
cat /dev/null > ${StrX}_CLA_Kraken.fasta

# Extract (using awk) the Bacillus fasta from the CLA output into a new file. We also delete with sed all the empty lines (for a reason we get an empty line in the end of the file!)
while read line; do
    awk -v pattern="$line" 'BEGIN {RS=">"} $0 ~ pattern {print ">"$0}' final_Scaffolds.fa | sed -r '/^\s*$/d' >> ${StrX}_CLA_Kraken.fasta
 done < Bacillus_fasta.headers


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
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
elapsed=$SECONDS
echo Time taken: `printf '%dd %dh:%dm:%ds\n' $((elapsed/86400)) $((elapsed%86400/3600)) $((elapsed%3600/60)) $((elapsed%60))`

exit 0
