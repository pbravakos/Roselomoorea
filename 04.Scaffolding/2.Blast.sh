#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=126675
# #SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Blastn"
#SBATCH --output=Blastn_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END

	This script takes as input one argument.
	For Strain 01 the correct usage is: 
	sbatch ${ScriptName} Strain01

	We are running Blast on the scaffolds that we were created by CLA (which was after the Spades Assembly!).

	IMPORTANT!
	Based on the results from this BLAST search we are going to select the reference genomes for the downstream analysis! In order to achieve that, we may have to rerun CLA, and then repeat this BLAST step again! It is like an endless cycle, but when we are quite sure about the reference we selected, we can escape the cycle!

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

OutputDir=${HOME}/Titlos_ktisis/BlastN/${StrainX}
CLADir=${HOME}/Titlos_ktisis/CLA/${StrainX}/CLA-Results

# Blast parameters that can be changed
OutFile=${StrX}_${StrainCode}_blastn_CLA
QueryMultiFasta=${StrX}_${StrainCode}_CLA_scaffolds.fasta
DataBase=nt
Task=blastn
WordSize=15
Eval=1e-15
OutFormat="6 qseqid sacc length pident qlen qcovs evalue bitscore staxids sscinames"  # also stitle same as sscinames
           
MaxSeqs=10  # Maximum number of aligned sequences to keep. Not applicable for outfmt <= 4 
MaxDesc=20
MaxAln=15



# Export the Blastn binary to the PATH.
export PATH="${HOME}/Software/ncbi-blast-2.11.0+/bin:$PATH"

# Export the local database directory.
export BLASTDB="/mnt/dbs/ncbi/nt"


#--------------------------SANITY CHECKS -------------------------------------------------------------------#
# Check for the existence of the CLA fasta output. 
if [[ ! -s ${CLADir}/final_Scaffolds.fa ]]
then
    echo "final_Scaffolds.fa cannot be found in ${CLADir}. Please run CLA first." >&2
    generalInfo >&2
    exit 1
fi

#--------------------------------------------------------------------------------------------------------------------

# Print the blastn version that is going to be used for the subsequent analyses.
echo
blastn -version
echo

[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}
cd ${OutputDir}

# Create a link of the CLA output scafold fasta file.
if [[ ! -e ${QueryMultiFasta} ]]; then
    ln -s ${CLADir}/final_Scaffolds.fa ${QueryMultiFasta}
fi

# Next, we will split the multifasta input file, in order to have the blast output separated for each scaffold. Otherwise it is too difficult to read and interpret the results!
# Split multifasta scaffolds file to multiple single fasta files and print those files with a new name.
awk -v a="$StrX" '/^>/{s=a"_CLA_Scaffold_"++d".fasta"} {print > s}' ${QueryMultiFasta}

# Find how many sequences exist in the multifasta file, and assign this number to a new parameter.
numSeqs=$(grep -c ">" ${QueryMultiFasta})

# Check that $numSeqs is greater than 1. Otherwise exit the script.
if [[ $numSeqs -lt 1 ]]; then
    echo "${QueryMultiFasta} seems to not have any fasta sequences. Please check the file."  >&2
    exit 1
fi

# Do a blast search for each of the single fasta files and save the output to separate files.
for i in $(seq 1 $numSeqs)
do
    QueryFasta=${StrX}_CLA_Scaffold_${i}.fasta
    blastn -query ${QueryFasta} \
	-db $DataBase \
	-out  ${OutFile}_Scaffold_${i}.tsv \
	-task $Task \
	-word_size $WordSize \
	-evalue $Eval  \
	-outfmt "$OutFormat" \
	-num_threads $SLURM_NTASKS \
	-max_target_seqs $MaxSeqs  
	# This is for the computer readable format. Check also the Note-dicussion at the end of the current script.
    rm ${QueryFasta}
done
   
   
#blastn -query ${QueryFasta} \
#	-db $DataBase \
#	-out ${OutFile}_Scaffold_${i}.html \
#	-task $Task \
#	-word_size $WordSize \
#	-evalue $Eval \
#	-remote \
#	-num_alignments $MaxAln \
#	-html

## We repeat the blast search but this time the output is a single tsv file!
#OutFile=${StrX}_blastn_CLA.tsv



# NOTE:
# We could run the two blastn searches (the first is considered the html output and the second the tsv output) with srun, each on a single core (in total two cores!) and gain some time, but we choose not to do it that way because it could mean that one blast search could end first and then one reserved core would not be used till the second blast search ends.


## Next, we apply some sorting for easier interpretation of the results. This step is used only when the output is not html but tsv!
#sort -t$'\t' -r -n -k5 $OutFile | awk 'BEGIN{FS="\t"; OFS="\t";print "Query Seq-id","Subject accession","Subject Title","Expect value","Alignment length","Percentage of identical matches","Query Coverage Per Subject","Query sequence length"} { print $0 }' > ${StrX}_blastn_CLA_final.tsv

#-num_threads $SLURM_NTASKS \


<< ////

	NOTE:
	The invocation using the parameter "-max_target_seqs 1" simply returns the first good hit found in the database, not the best hit as one would assume. Worse yet, the output produced depends on the order in which the sequences occur in the database. For the same query, different results will be returned by BLAST when using different versions of the database even if all versions contain the same best hit for this database sequence. Even ordering the database in a different way would cause BLAST to return a different "top hit" when setting the max_target_seqs parameter to 1.
	The confusion is further compounded by the fact that in the online BLAST portal, the max_target_seqs parameter behaves in the expected way – the best (rather than first) N hits are returned.
	
	Reference: Misunderstood parameter of NCBI BLAST impacts the correctness of bioinformatics workflows 
	           https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty833/5106166 
	

	Τhis issue is not specific to -max_target_seqs (used with the computer readable output format). With human readable output formats, the maximum of -num_descriptions and -num_alignments is used in exactly the same way during the BLAST search.
	
	Reference: https://blastedbio.blogspot.com/2018/11/blast-max-alignment-limits-repartee-one.html


	We examined the new example and it became clear that the demonstrated behavior was a bug, resulting from an overly aggressive optimization, introduced in 2012 for BLASTN and MegaBLAST (DNADNA alignments). This bug has been fixed in the BLAST+ 2.8.1 release, due out in December 2018. The aberrant behavior seems to occur only in alignments with an extremely large number of gaps, which is the case in the example provided by Shah and collaborators.
	It is also important to be clear that the BLAST website at blast.ncbi.nlm.nih.gov uses the same software libraries as the BLAST+ executables and should produce the same results when the same
parameters are used 
	
	Reference: Reply to the paper: Misunderstood parameters of NCBI BLAST impacts the correctness of bioinformatics workflows 
	           https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty1026/5259186

////


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
