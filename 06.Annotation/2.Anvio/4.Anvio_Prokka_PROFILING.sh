#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="3rd_Anvio_Prokka"
#SBATCH --output=3rd_Anvio_Prokka_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END

	This is the 3rd part of the Anvio Prokka pipeline. 
	Check the tutorial:
	http://merenlab.org/2016/06/22/anvio-tutorial-v2/ "Anvi'o User Tutorial for Metagenomic Workflow"

	This script takes as input one argument. 
	For Strain 01 that would be: 
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

# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
StrainCode="HK7V7BBXY"
Assembler="Spades"  # Choose either "Pilon" or "Spades"
Prefix=${StrX}_${Assembler}

OutputDir=${HOME}/Titlos_ktisis/Anvio/${Prefix}
CntgDB=${Prefix}_contigs.db


if [[ ${Assembler} == "Pilon" ]]
then
    ContigDir=${HOME}/Titlos_ktisis/${Assembler}/${StrainX}
    Contig=${StrX}_${Assembler}_CLA.fasta

elif [[ ${Assembler} == "Spades" ]]
then
    ContigDir=${HOME}/Titlos_ktisis/${Assembler}/${StrainX}/Iter2
    Contig=scaffolds.fasta
else
    echo "Choose either 'Spades' or 'Pilon' for the 'Assembler'" >&2
    echo >&2
    generalInfo >&2
    exit 1
fi


ReadDir=${HOME}/Titlos_ktisis/BBnorm/${StrainX}

Paired1=${StrX}_${StrainCode}_R1_bbnorm.fastq
Paired2=${StrX}_${StrainCode}_R2_bbnorm.fastq
Single=${StrX}_${StrainCode}_prinseq_good_singletons.fastq

Bam_PE_File=${Prefix}_PE_sorted.bam
BAM_SE_File=${Prefix}_SE_sorted.bam


ProfileDirPE=${Prefix}_PE_BamProfile
ProfileDirSE=${Prefix}_SE_BamProfile



# Anvio command options that can be altered
ProfileDBMinCngLength=1000 # Default 2500, do not go below 1000


. /etc/profile.d/modules.sh
. "${HOME}/Software/miniconda3/etc/profile.d/conda.sh"
conda activate anvio-7
module load R/3.6.1
export LC_ALL=en_US.UTF-8
export PATH="${HOME}/Software/Diamond:$PATH"
export PATH="${HOME}/Software/gms2_linux_64:$PATH"
export PATH="${HOME}/Software/ncbi-blast-2.10.1+/bin:$PATH"
export PATH="${HOME}/Software/bwa-0.7.17:$PATH"
export PATH="${HOME}/Software/bwa-mem2:$PATH"


#--------------------------SANITY CHECKS -------------------------------------------------------------------#
# Check that all the reads can be found
if [[ ! -s ${ReadDir}/$Paired1 ]] && [[ ! -s ${ReadDir}/$Paired2 ]] && [[ ! -s ${ReadDir}/$Single ]]
then
    echo "Some or all of the fastq files are missing. Please check the folder ${ReadDir}" >&2
    echo  >&2
    generalInfo >&2
    exit 1
fi

# Check that the contigs fasta file exists in the working folder
if [[ ! -s ${ContigDir}/${Contig} ]]
then
    echo "The contigs fasta file is not present in the working directory." >&2
    echo  >&2
    generalInfo >&2
    exit 1
fi


#------------------------------------------------------------------------------------------------------------
# Create output directory.
[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}
cd ${OutputDir}

# We also need the same fasta at the working directory for BWA to work!
ln -s ${ContigDir}/${Contig}

echo
echo "				The analysis begins!!" 
echo

# ATTENTION!
# BWA needs the contig file in the same directory where it is running!!

echo '		BWA MAPPING BEGINS'
echo
echo "BWA version is 0.7.17" 
echo
echo
# Create bwa index for ${Assembler}
bwa index ${Contig}
# Align reads to contigs with bwa mem.
# Turn sam to bam and finally sort by coordinates
bwa mem -t $SLURM_NTASKS ${Contig} ${ReadDir}/${Paired1} ${ReadDir}/${Paired2} |\
samtools view -hu -@ $SLURM_NTASKS - |\
samtools sort -@ $SLURM_NTASKS - -o ${Bam_PE_File}
# Index the sorted bam file
samtools index ${Bam_PE_File}

bwa mem -t $SLURM_NTASKS ${Contig} ${ReadDir}/${Single} |\
samtools view -hu -@ $SLURM_NTASKS - |\
samtools sort -@ $SLURM_NTASKS - -o ${BAM_SE_File}
samtools index ${BAM_SE_File}



echo 
echo '-------------------------BWA MAPPING COMPLETE----------------------------------------'
echo

echo
echo '		Profile BAM files'
anvi-profile --contigs-db ${CntgDB} --num-threads $SLURM_NTASKS --input-file ${Bam_PE_File} \
             --output-dir ${ProfileDirPE} --overwrite-output-destinations \
             --sample-name ${Prefix}_PE_Profile --min-contig-length ${ProfileDBMinCngLength} \
             --cluster-contigs --profile-SCVs

anvi-profile --contigs-db ${CntgDB} --num-threads $SLURM_NTASKS --input-file ${BAM_SE_File} \
             --output-dir ${ProfileDirSE} --overwrite-output-destinations \
             --sample-name ${Prefix}_SE_Profile --min-contig-length ${ProfileDBMinCngLength} \
             --cluster-contigs --profile-SCVs



echo
echo '--------------------------------BAM FILE PROFILING ANALYSIS COMPLETED-----------------------'
echo
echo
# touch ${Prefix}_ProfileDB_Anvio_decription.txt
echo "${Prefix} Profile database created from the ${Assembler} output scaffold fasta files. Gene finding was done with Prodigal (V2.6.3: February, 2016) and annotation in Prokka v1.14.6, using a genus specific Database, after downloading all Rossellomorea complete genomes in genbank format from Refseq. Also genes from GenemarkS2 (http://exon.gatech.edu/GeneMark/genemarks2.cgi) gene caller have been added, only when these gene calls were not found by Prodigal. There are two different profiles. One corresponds to the paired end reads, one to the single end reads. Scaffolds were created only from one Illumina run!" >  ${Prefix}_ProfileDB_Anvio_decription.txt
echo
echo '		Merge BAM PROFILES'
echo

anvi-merge --contigs-db ${CntgDB} --output-dir ${Prefix}_Merged_Profiles \
           --sample-name ${Prefix}_Merged_Profiles --overwrite-output-destinations \
           --description ${Prefix}_ProfileDB_Anvio_decription.txt \
           ${ProfileDirPE}/PROFILE.db ${ProfileDirSE}/PROFILE.db

echo
echo '------------------------MERGING PROFILES COMPLETED----------------------------------------'
echo
echo "			Create a new Collection uniting all the bins into one!"

# In our case probably it is meaningful to unite all bins in one new collection and examine them all together.
anvi-script-add-default-collection -p ${Prefix}_Merged_Profiles/PROFILE.db -c ${Prefix}_contigs.db

# Subsequently, in order to view the new collection use the following command:
# ONLY ON GUI ENVIRONMENTS! anvi-interactive -c {CntgDB} -p ${Prefix}_Merged_Profiles/PROFILE.db --gene-mode --collection-name "DEFAULT" --bin-id EVERYTHING

echo 
echo "-------------------------------DEFAULT COLLECTION CREATED!---------------------------------"
echo

#rm ${Prefix}_ProfileDB_Anvio_decription.txt *.bam* *.fasta.*
#echo

<< ////
	NEXT STEPS:
	Download to a local pc the following files and folders (with all their contents): $CntgDB $ProfileDirMeReads
	Start the interactive interface on a gui enabled pc by typing:
	anvi-interactive -p ${Prefix}_Merged_Profiles/PROFILE.db -c ${CntgDB} --title ${Prefix}_Anvio_results --taxonomic-level t_species

	The next commands are not actually needed! Everything can be done from the interactive interface of the above command!
	Instead of contigs one can also explore the genes AND their annotations:
	anvi-interactive -c{CntgDB} -p ${Prefix}_Merged_Profiles/PROFILE.db --gene-mode \
                         --collection-name "CONCOCT" --bin-id Bin_1 #Bin_1 can be substituted with Bin_2 and so on
	To check the available collections and the number of bins:
	anvi-show-collections-and-bins -p ${Prefix}_Merged_Profiles/PROFILE.db
	anvi-interactive -p ${Prefix}_Merged_Profiles/PROFILE.db -c {CntgDB} --list-collections
////

echo "======================================================================"
echo "SLRUM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
elapsed=$SECONDS
echo Time taken: `printf '%dd %dh:%dm:%ds\n' $((elapsed/86400)) $((elapsed%86400/3600)) $((elapsed%3600/60)) $((elapsed%60))`

exit 0
