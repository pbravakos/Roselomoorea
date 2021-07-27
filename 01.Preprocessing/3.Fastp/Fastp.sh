#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=16
# #SBATCH --mem=128000
# #SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Fastp"
#SBATCH --output=Fastp_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END

	This script takes as input one argument. 
	For Strain 01 that would be: 
	script.sh Strain01
	This script runs from the master folder of Fastp and then we change directory to the folder of each specific Strain folder.

	Attention:
	Fastp can take up to 16 threads no more. If more than 16 threads are inserted, it runs normally but uses only the 16 threads.
	
	NOTE:
	Fastp cannot produce singleton reads, meaning, it filters out either both reads or none of them. In order to be able to produce singleton reads we have to repeat 

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
OutputDir=${HOME}/Titlos_ktisis/Fastp/${StrainX}

ReadsDir=${HOME}/Titlos_ktisis/BBsplit/${StrainX}
StrainCode="HK7V7BBXY"
Read1=${StrainX}_${StrainCode}_R1_NOphix.fastq
Read2=${StrainX}_${StrainCode}_R2_NOphix.fastq


# Fastp parameters that can be changed
TrimFront1=6  # Trimm that number of bases in the front of the read 1
TrimFront2=6  # Trimm that number of bases in the front of the read 2
TrimTail1=1  # Trimm that number of bases in the tali of the read 1
TrimTail2=1  # Trimm that number of bases in the tail of the read 2
CutMeanQuality=17  # When triiming by a mean, this is the mean phred quality score that will be kept in the sliding window
QualifyBaseQuality=15  # The phred quality score that determines whether a base will be characterized as qualified or not.
UnqualPercLim=20   # Percentage of reads that are allowed to be unqualified in a read passing the filters.
NBaseLimit=1   # if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
LengthRequired=40  #  reads shorter than length_required will be discarded, default is 15. (int [=15])
Compression=4 # compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4. (int [=4]) 
OverLenReq=30  # the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default. (int [=30])
OverDiffLim=4  # the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default. (int [=5])
OveDiffPercLim=10  # the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%.
PolyGMinLen=10 # the minimum length to detect polyG in the read tail. 10 by default. (int [=10])
OverRepreSampling=20 # one in (--overrepresentation_sampling) reads will be computed for overrepresentation analysis (1~10000), smaller is slower, default is 20. (int [=20])
ComplTres=30 # the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).
# Below are the Truseq adapters
Adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"  # The Illumina adapter1
Adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  # The Illumina adapter2



export PATH="${HOME}/Software/fastp:$PATH"
export PATH="/home1/${USER}/Software/pigz-2.4/bin:$PATH"
export LC_ALL=en_US.UTF-8


[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}

cd ${OutputDir}

echo
echo "Fastp Version:"
fastp --version
echo


# Preliminary output
HtmlPrelimOut=${StrainX}_preliminary_Fastp.html
JsonPrelimOut=${StrainX}_preliminary_Fastp.json
Read1PrelimOut=${StrainX}_R1_fastp_preliminary.fastq.gz
Read2PrelimOut=${StrainX}_R2_fastp_preliminary.fastq.gz

if [[ ! -s ${HtmlPrelimOut} ]]; then
    echo "Now we will just get the statistics without any actuall filtering"
    # First run fastq without any actual filtering, just to check on the quality and select the apropriate filters.
    fastp --in1 ${ReadsDir}/${Read1} \
          --in2 ${ReadsDir}/${Read2} \
          --out1 ${Read1PrelimOut} \
          --out2 ${Read2PrelimOut} \
          --compression 1 \
          --disable_quality_filtering \
          --disable_adapter_trimming \
          --html ${HtmlPrelimOut} \
          --json ${JsonPrelimOut} \
          --report_title "${StrainX} preliminary report" \
          --thread $SLURM_NTASKS \
          --verbose
    rm ${Read1PrelimOut} ${Read2PrelimOut} ${JsonPrelimOut}
    
    # Run this script once more, to filter the reads.
    sbatch --dependency=afterok:"$SLURM_JOB_ID" $SLURM_SUBMIT_DIR/${ScriptName}
else
   echo "Now we will do the filtering and trimming!"
   # Run again, this time with the filters active.   
    fastp --in1 ${ReadsDir}/${Read1}\
          --in2 ${ReadsDir}/${Read2} \
          --out1 ${StrX}_${StrainCode}_R1_fastp.fastq \
          --out2 ${StrX}_${StrainCode}_R2_fastp.fastq \
          --compression ${Compression} \
          --dont_overwrite \
          --adapter_sequence ${Adapter1} \
          --adapter_sequence_r2 ${Adapter2} \
          --html ${StrX}_${StrainCode}_Fastp.html \
          --json ${StrX}_${StrainCode}_Fastp.json \
          --report_title "Initial Quality report for ${StrainX} ${StrainCode}" \
          --thread ${SLURM_NTASKS} \
          --verbose \
          --trim_front1 ${TrimFront1} \
          --trim_front2 ${TrimFront2} \
          --trim_tail1 ${TrimTail1} \
          --trim_tail2 ${TrimTail2} \
          --disable_trim_poly_g \
          --cut_tail \
          --cut_mean_quality ${CutMeanQuality} \
          --qualified_quality_phred ${QualifyBaseQuality} \
          --unqualified_percent_limit ${UnqualPercLim} \
          --n_base_limit ${NBaseLimit} \
          --length_required ${LengthRequired} \
          --correction \
          --overrepresentation_analysis
fi








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
