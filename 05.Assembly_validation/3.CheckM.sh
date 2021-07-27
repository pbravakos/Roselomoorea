#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=126675
# #SBATCH --mem-per-cpu=20000
#SBATCH --job-name="CheckM"
#SBATCH --output=CheckM_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END
	
	This script takes as input one argument. 
	For Strain 01 that would be: 
	sbatch ${ScriptName} Strain01
	This script runs from the master folder of CheckM and then we change directory to the folder of each specific Strain folder.
	
	In order to run this bash file for different strains changes have to take place:
	a) In the initial parameters given in the start of this file like the fastq read files (and associated changes of directory structure in each command).
	b) In the taxon_set command parameters for the taxonomy analysis (i.e. change the genus or species). 
	c) Someone could also change (for better or worse!) the individual parameters of each given command.
	
	NOTE:
	We have already run the command 
	checkm data setRoot /home1/pbravakos/Software/CheckM-1.1.3/Database
	
	NOTE:
	To check the available species, genera and families run
	checkm taxon_list

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
OutputDir=${HOME}/Titlos_ktisis/CheckM/${StrX}_Spades
Assembler=Spades  # Options are: "Pilon" or "Spades"

#CheckMDir=${HOME}/Software/CheckM-1.1.3/bin

# Functions are not exported by default to be made available in subshells.
. "${HOME}/Software/miniconda3/etc/profile.d/conda.sh"
conda activate
export LC_ALL=en_US.UTF-8

ReadDir=${HOME}/Titlos_ktisis/BBnorm/${StrainX}

# IMPORTANT
# Choose the correct scaffold directory!
#ContigDir=${HOME}/Titlos_ktisis/${Assembler}/${StrainX}       # Choose this for Pilon!
ContigDir=${HOME}/Titlos_ktisis/${Assembler}/${StrainX}/Iter2  # Choose this for Spades!

# IMPORTANT
# Choose the correct scaffold from bellow!
ContigAssembler=scaffolds.fasta
#ContigAssembler=${StrX}_Pilon_CLA.fasta



StrXAssembler=${StrX}_${Assembler}


Paired1=${StrX}_${StrainCode}_R1_bbnorm.fastq
Paired2=${StrX}_${StrainCode}_R2_bbnorm.fastq
Single=${StrX}_${StrainCode}_prinseq_good_singletons.fastq

Bam_PE_File=${StrX}_PE_${Assembler}_sorted.bam
BAM_SE_File=${StrX}_SE_${Assembler}_sorted.bam



family='Bacillaceae'
genus='Bacillus'
# species=''

BinFolder=${StrX}_bin_folder
TempFolder=${StrX}_temp
PlotFolder=${StrX}_plot_folder
TaxFolder=${StrX}_taxonomy_output
LinFolder=${StrX}_lineage_output
TreeFolder=${StrX}_output_tree


CovFile=${StrX}_${Assembler}.coverage
LinMarkers=${StrX}_${Assembler}_lineage.markers
TaxMarkers=${StrX}_${Assembler}_taxonomy.markers
TetraFile=${StrX}_${Assembler}.tetra



#-------------------------------------------------------------------------------------------------
[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}
cd ${OutputDir}

mkdir ${TempFolder}
mkdir ${BinFolder}
cd ${BinFolder}
ln -s ${ContigDir}/${ContigAssembler} ${StrXAssembler}.fna
cd ..
# We also need the same fasta at the working directory for BWA to work!
ln -s ${ContigDir}/${ContigAssembler}


echo '		CHECKM PIPELINE FOR GENOME SCAFFOLDS'
echo

echo 
echo '		BWA MAPPING BEGINS'
echo 
# Create bwa index for fasta assemblies
bwa index ${ContigAssembler}
# Align reads to contigs with bwa mem and bwa bwasw.
# Turn sam to bam and finally sort by coordinates
bwa mem -t $SLURM_NTASKS ${ContigAssembler} ${ReadDir}/${Paired1} ${ReadDir}/${Paired2} |\
samtools view -hu -@ $SLURM_NTASKS - |\
samtools sort -@ $SLURM_NTASKS - -o ${Bam_PE_File}
# Index the sorted bam file
samtools index ${Bam_PE_File}


bwa mem -t $SLURM_NTASKS ${ContigAssembler} ${ReadDir}/${Single} |\
samtools view -hu -@ $SLURM_NTASKS - |\
samtools sort -@ $SLURM_NTASKS - -o ${BAM_SE_File}
samtools index ${BAM_SE_File}



echo 
echo '!!!!!!!!!!!!!!!BWA MAPPING COMPLETE!!!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		GC PLOT ANALYSIS BEGINS'
echo 

# ATTENTION!! IMPORTANT!!!
# matplotilib version HAS to be less than 2.2 for this to run. Except for some plots (like for tetranuclotides) for which the version HAS to be 1.3.1!!!!!!


checkm gc_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 8 --gc_window_size 2000 --gc_bin_width 0.01 ${BinFolder} ${PlotFolder} 50

echo 
echo '!!!!!!!!!!!!!!!ANALYSIS OF GC PLOT COMPLETE!!!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		Nx ANALYSIS BEGINS'
echo 

checkm nx_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 3.5 --step_size 0.01 ${BinFolder} ${PlotFolder}

echo ''
echo '!!!!!!!!!!!!Nx ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!'
echo ''
echo ''
#echo '		Sequence Length ANALYSIS BEGINS'
#echo ''

#checkm len_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 6.5 ${BinFolder} ${PlotFolder}

#echo 
#echo '!!!!!!!!!!!!Length ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		LENGTH HISTOGRAM ANALYSIS BEGINS'
echo 

checkm len_hist --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 6.5 ${BinFolder} ${PlotFolder}

echo ''
echo '!!!!!!!!!!!!LENGTH HISTOGRAM ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!!'
echo ''
echo ''
echo '		COVERAGE ANALYSIS BEGINS'
echo ''

checkm coverage --all_reads --threads $SLURM_NTASKS ${BinFolder} ${CovFile} ${Bam_PE_File} ${BAM_SE_File}

echo 
echo '!!!!!!!!!!!!COVERAGE ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		TETRANUCLEOTIDE ANALYSIS BEGINS'
echo 

checkm tetra -t 1 ${BinFolder}/${StrXAssembler}.fna ${TetraFile}

echo ''
echo '!!!!!!!!!!!!!!!!TETRA ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!'
echo ''
echo '		TREE ANALYSIS BEGINS'
echo ''

checkm tree --ali --nt --threads $SLURM_NTASKS --pplacer_threads $SLURM_NTASKS --tmpdir ${TempFolder} ${BinFolder} ${TreeFolder}
echo ''
checkm tree_qa --out_format 1 --file ${StrXAssembler}_tree_format1.tsv --tab_table --tmpdir ${TempFolder} ${TreeFolder}
echo ''
checkm tree_qa --out_format 2 --file ${StrXAssembler}_tree_format2.tsv --tab_table --tmpdir ${TempFolder} ${TreeFolder}
echo ''
# Format 3 is nwk tree but tips have only img genome ids.
checkm tree_qa --out_format 3 --file ${StrXAssembler}_tree_format3.nwk --tmpdir ${TempFolder} ${TreeFolder}
echo ''
# Format 3 is nwk tree but tips of the tree have img genome ids and genome names
checkm tree_qa --out_format 4 --file ${StrXAssembler}_tree_format4.nwk --tmpdir ${TempFolder} ${TreeFolder}
echo ''
checkm tree_qa --out_format 5 --file ${StrXAssembler}_tree_format5.msa --tmpdir ${TempFolder} ${TreeFolder}

echo 
echo '!!!!!!!!!!!!!TREE ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		LINEAGE ANALYSIS BEGINS'

checkm lineage_set --tmpdir ${TempFolder}  --unique 10 --multi 10 ${TreeFolder} ${LinMarkers} 

checkm analyze --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${BinFolder} ${LinFolder}

# AAI stands for Amino Acid Identity in the checkm qa manual. Look here https://github.com/Ecogenomics/CheckM/wiki/Genome-Quality-Commands
checkm qa --out_format 1 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_lineage_qa_format1.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

checkm qa --out_format 2 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_lineage_qa_format2.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

checkm qa --out_format 3 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_lineage_qa_format3.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

checkm qa --out_format 4 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_lineage_qa_format4.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

checkm qa --out_format 5 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_lineage_qa_format5.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

checkm qa --out_format 6 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_lineage_qa_format6.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

checkm qa --out_format 7 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_lineage_qa_format7.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

checkm qa --out_format 8 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_lineage_qa_format8.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

checkm qa --out_format 9 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_lineage_qa_format9.fasta --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

echo 
echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!LINEAGE ANALYSIS COMPLETE!!!!!!!!!!!!!!!'
echo 
echo 
echo '		TAXONOMY ANALYSIS BEGINS'
echo 

# The list of available taxonomic-specific marker sets which can be inserted as input in the taxon_set command can be viewed with the taxon_list command.

checkm taxon_set --tmpdir ${TempFolder} genus ${genus} ${TaxMarkers}
# A more specialized command depending on the specific species under analysis would be to specify a species name instead of a genus e.g.:
# NOT NEDED! checkm taxon_set --tmpdir ${TempFolder} species ${species} ${TaxMarkers}

checkm analyze --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${BinFolder} ${TaxFolder}

checkm qa --out_format 1 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_taxonomy_qa_format1.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

checkm qa --out_format 2 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_taxonomy_qa_format2.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

checkm qa --out_format 3 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_taxonomy_qa_format3.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

checkm qa --out_format 4 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_taxonomy_qa_format4.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

checkm qa --out_format 5 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_taxonomy_qa_format5.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

checkm qa --out_format 6 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_taxonomy_qa_format6.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

checkm qa --out_format 7 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_taxonomy_qa_format7.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

checkm qa --out_format 8 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_taxonomy_qa_format8.tsv --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

checkm qa --out_format 9 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXAssembler}_taxonomy_qa_format9.fasta --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

echo 
echo '!!!!!!!!!!!!!!!!!!!TAXONOMY ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 
echo 
#echo '		PLOTS FOR completeness AND contamination ANALYSIS BEGINS'
#echo 

#checkm bin_qa_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --row_height 0.3 ${TaxFolder} ${BinFolder} ${PlotFolder}
#cd ${PlotFolder}
#mv bin_qa_plot.pdf ${StrXAssembler}_taxonomy_Contamination.pdf
#cd ..
#echo 
#checkm bin_qa_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --row_height 0.3 ${LinFolder} ${BinFolder} ${PlotFolder}
#cd ${PlotFolder}
#mv bin_qa_plot.pdf ${StrXAssembler}_lineage_Contamination.pdf
#cd ..
#echo ''
#echo '!!!!!!!!!!!!!!!PLOTS FOR completeness AND contamination ANALYSIS COMPLETE!!!!!!!!!!!!!'
echo 
echo 
echo '		Coding plot analysis BEGINS'
echo 

checkm coding_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 3.5 --cd_window_size 5000 --cd_bin_width 0.01 ${TaxFolder} ${BinFolder} ${PlotFolder} 50

echo 
echo '!!!!!!!!!!!!!CODING PLOT ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!'
echo 
echo 
#echo '		PLOT OF GC and COVERAGE BEGINS'

#checkm par_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 6.5 ${TaxFolder} ${BinFolder} ${PlotFolder} ${CovFile}

#echo 
#echo '!!!!!!!!!!!!!!!!!!!PLOT OF GC and COVERAGE COMPLETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		GC, CD, and TD distribution plots ANALYSIS BEGINS'
echo
checkm dist_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 8 --gc_window_size 5000 --td_window_size 5000 --cd_window_size 5000 ${TaxFolder} ${BinFolder} ${PlotFolder} ${TetraFile} 50
echo
echo '!!!!!!!!!!!!!!!!!!!!!GC, CD, and TD distribution plots ANALYSIS COMPLETE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo
echo ''

echo 
echo '		PLOTS ANALYSIS BEGINS'
echo 
echo 
## ATTENTION!!!! NOT WORKING!!!
## Tetra pca needs matplotlib 1.3.1 to run!!! Otherwise it doesn't run at all!!!!!!!!!!!!
#checkm tetra_pca --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 6.5 ${BinFolder} ${PlotFolder} ${TetraFile}

echo 
echo 

# ATTENTION!!!! NOT WORKING!!!
# Tetra plot needs matplotlib 1.3.1 to run!!! Otherwise it doesn't run at all!!!!!!!!!!!!

checkm tetra_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 3.5 --td_window_size 1000 --td_bin_width 0.01 ${TaxFolder} ${BinFolder} ${PlotFolder} ${TetraFile} 50
echo
echo 
echo 
echo "GC bias plot"
echo
echo

# ATTENTION!!
# NOT WORKING!! Probably needs matplotlib 1.3.1 to run!!!!!!
checkm gc_bias_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 3.5 --window_size 1000 --all_reads --threads $SLURM_NTASKS ${BinFolder} ${PlotFolder} ${Bam_PE_File}

echo 
echo '!!!!!!!!!!!!!!!!PLOTS ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		MAPPED READS % ANALYSIS BEGINS'

checkm profile --file ${StrXAssembler}_mapped_reads_percentage.tsv --tab_table ${CovFile} 

echo '!!!!!!!!!!!!!!!!!!!!MAPPED READS % ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!!!!!!!'


#rm -r ${TempFolder} *.bam* *fasta.*



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
