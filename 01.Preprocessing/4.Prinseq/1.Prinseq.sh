#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6400
# #SBATCH --mem=128000
#SBATCH --job-name="Prinseq_before_trimming"
#SBATCH --output=Prinseq_before_trimming_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END

	This script takes as input one argument. 
	For Strain 01 that would be: 
	script.sh Strain01
	This script runs from the master folder of Prinseq and then we change directory to the folder of each specific Strain folder.

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


# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
PrinseqDir=${HOME}/Software/prinseq-lite-0.20.4
OutputDir=${HOME}/Titlos_ktisis/Prinseq/${StrainX}

FastpDir=${HOME}/Titlos_ktisis/Fastp/${StrainX}

StrainCode="HK7V7BBXY"
FastpR1=${StrX}_${StrainCode}_R1_fastp.fastq
FastpR2=${StrX}_${StrainCode}_R2_fastp.fastq

PrinseqGoodOut=${StrX}_${StrainCode}_prinseq_good_R
PrinseqBadOut=${StrX}_${StrainCode}_prinseq_bad_R
PrinseqLog=${StrX}_${StrainCode}_prinseq_trimming.log

# Prinseq parameters that can be changed
OutFormat=3   # 1 (FASTA only), 2 (FASTA and QUAL), 3 (FASTQ), 4 (FASTQ and FASTA), or 5 (FASTQ, FASTA and QUAL) 
MinLen=41  # Filter sequence shorter than $MinLen
MinQualMin=17  # Filter sequence with quality score mean below $MinQualMin
NsMaxN=1  # Filter sequence with more than $NsMaxN Ns
TrimQualRight=16  # Trim sequence by quality score from the 3'-end with this threshold score
TrimQualWindow=10   # The sliding window size used to calculate quality score by type
TrimQualStep=5  # Step size used to move the sliding window
TrimQualType=min  # Type of quality score calculation to use. Allowed options are min, mean, max and sum
TrimQualRule=lt  # Rule to use to compare quality score to calculated value. Allowed options are lt (less than), gt (greater than) and et (equal to)

LCMethod=entropy # "dust" or "entropy"
LCTreshold=40  # between 0 and 100. The dust method uses this as maximum allowed score and the entropy method as minimum allowed value.
Derep=1 # 1 (exact duplicate), 2 (5' duplicate), 3 (3' duplicate), 4 (reverse complement exact duplicate), 5 (reverse complement 5'/3' duplicate)

 
#export PATH="${HOME}/Software/prinseq-lite-0.20.4:$PATH"
export LC_ALL=en_US.UTF-8

. /etc/profile.d/modules.sh
module load perl/5.32.0
#==================================================================================================

# Check that Fastp files have been generated in the previous step of this pipeline.
if [[ ! -s ${FastpDir}/${FastpR1} || ! -s ${FastpDir}/${FastpR2} ]]
then
    echo "The previous step of this pipeline (Fastp) has not produced the expected files ${FastpR1} and ${FastpR2}. " >&2
    generalInfo >&2
    exit 1
fi


[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}

cd ${OutputDir}

ln -s ${FastpDir}/${FastpR1}
ln -s ${FastpDir}/${FastpR2}


#perl ${PrinseqDir}/prinseq-lite.pl	-fastq ${FastpR1} \
#					-fastq2 ${FastpR2} \
#					-graph_data ${StrX}_${StrainCode}_before_trimming.gd \
#					-out_good null \
#					-out_bad null \
#					-log ${StrX}_${StrainCode}_before_trimming.log

#perl ${PrinseqDir}/prinseq-graphs-noPCA.pl	-i ${StrX}_${StrainCode}_before_trimming.gd \
#						-html_all


perl ${PrinseqDir}/prinseq-lite.pl	-fastq ${FastpR1} \
					-fastq2 ${FastpR2} \
					-out_format $OutFormat \
					-min_len $MinLen \
					-min_qual_mean $MinQualMin \
					-ns_max_n $NsMaxN \
					-noniupac \
					-trim_qual_right $TrimQualRight \
					-trim_qual_type $TrimQualType \
					-trim_qual_rule $TrimQualRule \
					-trim_qual_window $TrimQualWindow \
					-trim_qual_step $TrimQualStep \
					-out_good ${PrinseqGoodOut} \
					-out_bad ${PrinseqBadOut} \
					-log ${PrinseqLog} \
					-lc_method ${LCMethod} \
					-lc_threshold ${LCTreshold}


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

echo
echo "==============================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
elapsed=$SECONDS
echo Time taken: `printf '%dd %dh:%dm:%ds\n' $((elapsed/86400)) $((elapsed%86400/3600)) $((elapsed%3600/60)) $((elapsed%60))`

exit 0
