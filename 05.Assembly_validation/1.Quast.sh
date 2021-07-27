#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=126675
# #SBATCH --mem-per-cpu=6400
# Memory per node specification is in MB. It is optional.
#SBATCH --job-name="Quast"
#SBATCH --output=Quast_job_%j.out
# #SBATCH --error=Sderr_Quast_job_%j.txt
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END



generalInfo () {
    cat <<END
	
	This script takes as input one argument. 
	For Strain 01 that would be: 
	bash ${ScriptName} Strain01
	This script runs from the master folder of Pilon and then we change directory to the folder of each specific Strain folder.
	
	NOTE:
	Here we use the fastq files after Prinseq even if further downstream filtering has been done (e.g. contamination removal) because we want as much information as possible to correct our contigs!!
	IMPORTANT:
	Quast overloads the system, and there is currently an admin script that stops any job that exceeds a specific cpu load.
	In order to avoid this from happening, we can do one of the following options:
	1) Run the script with a reduced number of cpus e.g. 14 instead of the full 20 and allocate the whole memory to prohibit anyone else using the node.
	2) Run cpulimit https://github.com/opsengine/cpulimit. This software finds processes that have a cpu load above a specific treshold and reduce the load to a defined value.

	NOTE:
	To create the blast database we did the following:
	wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
	wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz


	makeblastdb -in nt.fasta -out nt -parse_seqids -dbtype nucl -input_type fasta -blastdb_version 5 -max_file_sz 300GB
	The -parse_seqids option is required to keep the original sequence identifiers. Otherwise makeblastdb will generate its own identifiers

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
StrNum=${StrainX/Strain/}
StrainCode="HK7V7BBXY"
OutputDir=${HOME}/Titlos_ktisis/Quast

QuastDir=${HOME}/Software/quast-5.1.0rc1
ReadDir=${HOME}/Titlos_ktisis/Prinseq/${StrainX}
SpadesDir=${HOME}/Titlos_ktisis/Spades/${StrainX}/Iter2
# CeleraDir=${HOME}/Pseudomonas/Celera/${StrainX}/9-terminator
PilonDir=${HOME}/Titlos_ktisis/Pilon/${StrainX}
RefDir=${HOME}/Titlos_ktisis/References/${StrainX}
#ContigDir=${HOME}/Pseudomonas/Seqtk/${StrainX}


Paired1=${StrX}_${StrainCode}_prinseq_good_R_1.fastq
Paired2=${StrX}_${StrainCode}_prinseq_good_R_2.fastq
Single=${StrX}_${StrainCode}_prinseq_good_singletons.fastq

RefGenome=Bacillus_vietnamensis_strain_151-6.fasta
RefFeatures=Rossellomorea_vietnamensis_strain_151-6.gff3
#RefGenome=`find ${RefDir} -name '*.fasta'`
#RefFeatures=`find ${RefDir} -name '*.gff3'`


export LC_ALL=en_US.UTF-8
export PATH="${HOME}/Software/gms2_linux_64:$PATH"
export PATH="${HOME}/Software/barrnap-0.9/bin:$PATH"
export BUSCO_CONFIG_FILE="${HOME}/Software_OLD/busco/config/config.ini"
declare -x AUGUSTUS_CONFIG_PATH="${HOME}/Software/augustus-3.3.3/config"
export PATH="${HOME}/Software/cpulimit-0.3.2/src:$PATH"
export PATH="${HOME}/Software/bwa-0.7.17:$PATH"
export PATH="${HOME}/Software/bwa-mem2:$PATH"


#------------------------------------------------------------------------------------------------------------

# Change directory to the Strain specific folder

echo "			This is the QUAST analysis for ${StrainX}"


[[ ! -d ${OutputDir} ]] && mkdir ${OutputDir}

# Run Quast!!
python3 ${QuastDir}/quast.py -r ${RefDir}/${RefGenome} \
				--min-contig 200 \
				-o ${StrX}_results \
				--threads $SLURM_NTASKS \
				--split-scaffolds \
				--features ${RefDir}/${RefFeatures} \
				--labels ${StrX}_CLA,${StrX}_Spades \
				--k-mer-stats \
				--k-mer-size 101 \
				--gene-finding \
				--conserved-genes-finding \
				--ambiguity-usage all \
				--rna-finding \
				--pe1 ${ReadDir}/${Paired1} \
				--pe2 ${ReadDir}/${Paired2} \
				--single ${ReadDir}/${Single} \
				 ${PilonDir}/${StrX}_Pilon_CLA.fasta ${SpadesDir}/scaffolds.fasta
				 
				 

## First, check whether the current number of reserved cpus is above a specific treshold.
#if (( SLURM_NTASKS > $(sinfo --Node --long | grep "${SLURM_JOB_NODELIST} " | awk '{print $5}')-6 ))
#then
#    # Quast overloads the system and it is being stopped automatically by the system.
#    # For this reason we run the "cpulimit" software to stop the overloading of the cpus.              
#    while sleep 120        # execute this loop every that many seconds
#    do 
#        # get the processes that use more than 125% of CPU
#        OverloadPIDs=$(top -b -n1 -H | tail -n +8 | awk -v user=${USER::6} '$2~user && $9>125 {print $1}')
#        for i in ${OverloadPIDs}
#        do
#            # Only run cpulimit if this process is not already being controlled.
#            if [[ -z $(ps aux | grep cpulimit | grep "\-p ${i}") ]]
#            then
#                # limit the processes to 105%
#                cpulimit -p ${i} -l 105 -z &
#            fi  
#        done
#        # If there are no quast (or ${SLURM_JOB_NAME}) related process then exit this loop.
#        if (( $(ps ux --user ${USER} | grep -i quast | wc -l) <= 1 ))
#        then
#            echo "No more processes related to quast!"
#            echo "Exiting the loop!"
#            break
##        # If cpulimit is currently controlling more than ${SLURM_NTASKS} cpus, then exit this loop.
##        elif  (( $(pgrep cpulimit | wc -l) >= SLURM_NTASKS ))
##        then
##            echo "Cpulimit is already running in many processes!"
##            echo "Exiting the loop!"
##            break
#        fi  
#    done
#fi
#QuastPIDs=$(ps aux | grep -i quast | awk '{print $2}')

#wait ${QuastPIDs} 2> /dev/null
#kill $(pgrep cpulimit)

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
