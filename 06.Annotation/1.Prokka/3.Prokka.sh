#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=126675
# #SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Prokka"
#SBATCH --output=Prokka_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END
	
	This script takes as input one argument. 
	For Strain 01 that would be: 
	sbatch ${ScriptName} Strain01
		
	This script runs from the master folder of Prokka and then we change directory to the folder of each specific Strain folder.

	IMPORTANT!!
	tbl2asn is a dependency for Prokka and has to be updated every year! 
	Check before running Prokka that you have the latest version on the path! 
	Download latest version from here:
	https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/
	
	NOTE:
	Check all the Prokka dependencies here:
	https://vicbioinformatics.com/software.prokka.shtml
	
#	wget http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
#	wget http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
#	wget http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/RF01690.fa.gz # Family: Bacillaceae-1 (RF01690)
#	gunzip RF01690.fa.gz
#	gunzip Rfam.cm.gz
#	/home1/pbravakos/Software/infernal-1.1.4-linux-intel-gcc/binaries/cmpress Rfam.cmscan
	

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
#StrainCode="HK7V7BBXY"

# IMPORTANT
# Choose the correct scaffold assembler!
Assembler="Spades"  # Choose either "Pilon" or "Spades"
OutputDir=${HOME}/Titlos_ktisis/Prokka


#Refdir=${HOME}/Titlos_ktisis/References/${StrainX}
#RefGb=BSGatlas.gb

# Parameters that could be changed
# Genus=Pseudomonas
Kingdom=Bacteria
Gram=pos   # -/neg +/pos (default '')
Gcode=11
ContigLen=200
Eval=1e-9
Genus="Rossellomorea"


if [[ ${Assembler} == "Pilon" ]]
then
    ContigDir=${HOME}/Titlos_ktisis/${Assembler}/${StrainX}       
    Contig=${StrX}_Pilon_CLA.fasta

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


if [[ ! -s ${ContigDir}/${Contig} ]]
then
    echo "Contig ${Contig} cannot be found in ${ContigDir}" >&2
    echo "Please place ${Contig} in the correct directory first! " >&2
    echo >&2
    generalInfo >&2
    exit 1
fi

#if [[ ! -s ${Refdir}/${RefGb} ]]
#then
#    echo "Reference ${RefGb} cannot be found in ${Refdir}" >&2
#    echo "Please place ${RefGb} in the correct directory first! " >&2
#    echo >&2
#    generalInfo >&2
#    exit 1
#fi

export LC_ALL=en_US.UTF-8
export PATH="${HOME}/Software/Prodigal/bin:$PATH"
export PATH="${HOME}/Software/ncbi-blast-2.10.1+/bin:$PATH"
export PATH="${HOME}/Software/hmmer-2.3.2/src:$PATH"
export PATH="${HOME}/Software/hmmer-3.3.2/bin:$PATH"
export PATH="${HOME}/Software/aragorn1.2.38:$PATH"
export PATH="${HOME}/Software/tbl2asn:$PATH"
export PATH="${HOME}/Software/infernal-1.1.4-linux-intel-gcc/binaries:$PATH"
export PATH="${HOME}/Software/infernal-1.1.4-linux-intel-gcc/easel/miniapps:$PATH"
export PATH="${HOME}/Software/barrnap-0.9/bin:$PATH"
export PATH="${HOME}/Software/minced-0.4.2:$PATH"
export PATH="${HOME}/Software/rnammer-1.2/bin:$PATH"
export PATH="${HOME}/Software/signalp-5.0b/bin:$PATH"
export PATH="${HOME}/Software/prokka/bin:$PATH"


# Print the Prokka version we are going to use in the analysis.
prokka --version
echo
which prokka
echo
#-------------------------------------------------------------------------------------------------
[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}
cd ${OutputDir}

# Create a soft link of the contigs fasta file
ln -s ${ContigDir}/${Contig}
## Create a soft link of the reference genbank file
#ln -s ${Refdir}/${RefGb}

prokka --force \
	--outdir ${StrX}_Prokka_${Assembler} \
	--prefix ${StrX}_Prokka_${Assembler} \
	--genus ${Genus} \
	--kingdom ${Kingdom} \
	--gcode ${Gcode} \
	--gram ${Gram} \
	--evalue ${Eval} \
	--mincontiglen ${ContigLen} \
	--cpus $SLURM_NTASKS \
	--rfam \
	${Contig}
#   	--proteins ${RefGb} \	


rm ${Contig}

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
