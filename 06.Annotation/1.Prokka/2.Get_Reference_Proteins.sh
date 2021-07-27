#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
# #SBATCH --mem=126675
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Seqret"
#SBATCH --output=Seqret_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END
	
	To run this script:
	sbatch ${ScriptName}
        
        Here we will download the gff3 file along with the fasta genome sequence file from BSGatlas and try 
        to convert these two files to genbank with emboss seqret software.
        The new genbank file will be used as input in Prokka.
        BSGatlas provides access to the annotation of Bacillus subtilis 168.
        This strain seems to be the most relevant to our reference strain Rossellomorea vietnamensis strain 151. 
	
	NOTE:
	To download EMBOSS we followed the instructions found here:
	http://emboss.open-bio.org/html/adm/ch01s01.html
	1.1.2. Downloading by Anonymous FTP
	
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

OutputDir=${HOME}/Titlos_ktisis/References/${StrainX}

EmbossDir=/home1/pbravakos/Software/EMBOSS-6.6.0/bin

# The correct option should be found from the download section in BSGatlas webpage.
BSGatlasVersion=BSGatlas_v1.1


export LC_ALL=en_US.UTF-8


# Print the Prokka version we are going to use in the analysis.
echo "Emboss seqret version: "
${EmbossDir}/seqret -version
echo
#-------------------------------------------------------------------------------------------------
[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}
cd ${OutputDir}


wget -N https://rth.dk/resources/bsgatlas/annotation/${BSGatlasVersion}.gff

if [[ ! -s ${BSGatlasVersion}.fasta ]]
then
    wget https://rth.dk/resources/bsgatlas/annotation/genome.fna
    mv genome.fna ${BSGatlasVersion}.fasta
fi

if [[ -s ${BSGatlasVersion}.gff && -s ${BSGatlasVersion}.fasta ]]
then
    ${EmbossDir}/seqret -sequence ${BSGatlasVersion}.fasta \
    			-feature \
    			-fformat gff \
    			-fopenfile ${BSGatlasVersion}.gff \
    			-osformat genbank \
    			-osname_outseq ${BSGatlasVersion} \
    			-ofdirectory_outseq gbk_file \
    			-auto

    [[ $? -eq 0 ]] && echo "Seqret run succesfully!! Files can be found in ${OutputDir}"
    echo
    echo
    # Rename the newly created genbank file to an extension that Prokka undestands.			
    if [[ -s ${BSGatlasVersion}.genbank ]]
    then
        mv ${BSGatlasVersion}.genbank ${BSGatlasVersion}.gb
    fi
        
else
    echo "${BSGatlasVersion}.gff and/or ${BSGatlasVersion}.fasta could not be downloaded. Please check how to download these files."
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


echo "==============================="
echo "SLRUM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
elapsed=$SECONDS
echo Time taken: `printf '%dd %dh:%dm:%ds\n' $((elapsed/86400)) $((elapsed%86400/3600)) $((elapsed%3600/60)) $((elapsed%60))`

exit 0
