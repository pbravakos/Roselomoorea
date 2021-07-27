#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=126675
#SBATCH --job-name="2nd_Anvio_Prokka"
#SBATCH --output=2nd_Anvio_Prokka_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END

	This script takes as input one argument. 
	For Strain 01 that would be: 
	sbatch ${ScriptName} Strain01

	This is the Anvio pipeline. PART TWO

	IMPORTANT!
	This script uses the busco bacillales single copy genes. It should be changed if we do not have bacillales!!
	We can use the program busco_hmms_for_anvio (https://github.com/guyleonard/busco_hmms_for_anvio) to download 
       and extract the bacillales_odb10.2020-03-06 hmm busco dataset of single copy genes.
       This dataset can be used as input to the anvi-run-hmms command
       We must edit the following files:
       "kind.txt" and change it to bacteria, if we have downloaded the bacterial dataset
	"noise_cutoff_terms.txt" if we want a value other than '1e-25"

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
   echo "Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9'!" >&2
   echo >&2
   generalInfo >&2
   exit 1
fi


# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
Assembler="Spades"  # Choose either "Pilon" or "Spades"

OutputDir=${HOME}/Titlos_ktisis/Anvio/${StrX}_${Assembler}

CntgDB=${StrX}_${Assembler}_contigs.db
InterproOutputDir=${OutputDir}/${StrX}_INTEPRO

InterProDir=${HOME}/Software/interproscan-5.52-86.0
BuscoDBDir=${HOME}/Software/busco/bacillales_anvio    # This is a Database specifically for gammaproteobacteria

KaijuDB=${HOME}/Software/kaiju/kaijudb

anviGetSeq='anvi-get-sequences-for-gene-calls'
anviImpTax='anvi-import-taxonomy-for-genes'

# InterproScan options that can change
Input=${StrX}_InterProScan_amino-acid-sequences.fasta   # Optional, path to fasta file that should be loaded on Master startup. 
                                                        # Alternatively, in CONVERT mode, the InterProScan 5 XML file to convert. 
OutputFileBase=${StrX}_${Assembler}_interpro-output   # Optional, base output filename (relative or absolute path). 
                                                   	# Note that this option, the --output-dir (-d) option and the --outfile (-o) option are mutually exclusive.  
                                                   	# The appropriate file extension for the output format(s) will be appended automatically. 
                                                   	# By default the input file path/name will be used. 


. /etc/profile.d/modules.sh
. "${HOME}/Software/miniconda3/etc/profile.d/conda.sh"
conda activate anvio-7
module load R/3.6.1
export LC_ALL=en_US.UTF-8
export PATH="${HOME}/Software/Diamond:$PATH"
export PATH="${HOME}/Software/gms2_linux_64:$PATH"
export PATH="${HOME}/Software/ncbi-blast-2.10.1+/bin:$PATH"
export PATH="${HOME}/Software/kaiju/bin:$PATH"

#--------------------------SANITY CHECKS -------------------------------------------------------------------#
# Check that Kaiju Database is installed.
if [[ ! -e ${KaijuDB}/refseq/kaiju_db_refseq.fmi ]] && [[ ! -e ${KaijuDB}/names.dmp ]] && [[ ! -e ${KaijuDB}/nodes.dmp ]]
then
    echo >&2
    echo "Kaiju Database was not installed properly. Please check the database installation again." >&2
    exit 1
fi

# Check whether interproscan executable exists or not
if [[ ! -r ${InterProDir}/interproscan.sh ]]; then
    echo >&2
    echo "Interproscan executable is missing!" >&2
    exit 1
fi

# Check that the Busco Directory can be found
if [[ ! -d $BuscoDBDir ]]; then
    echo "Busco Directory is missing! It was expected to be found here: $BuscoDBDir" >&2
    exit 1
fi


#------------------------------------------------------------------------------------------------------------
# Create output directory.
[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}
cd ${OutputDir}


# For safety reasons copy the amino acid fasta file to a new one!
if [[ -s ${StrX}_amino_acid_seqs.faa ]]
then
    cp ${StrX}_amino_acid_seqs.faa ${StrX}_InterProScan_amino-acid-sequences.fasta
else   
    echo >&2 
    echo "${StrX}_amino_acid_seqs.faa proteins file cannot be found" >&2
    echo "Please run the first part of this analysis to create the file" >&2
    exit 1;
fi


echo
echo
echo '		Start Interproscan analysis'
echo
# Instructions can be found here http://merenlab.org/2016/06/18/importing-functions/

<< ////

	IMPORTANT!
	The "cpu" option just overrides both of the following values in the INTERPROSCAN properties file:
	#Number of embedded workers at start time
	number.of.embedded.workers=1
	#Maximum number of embedded workers
	maxnumber.of.embedded.workers=35
	https://github.com/ebi-pf-team/interproscan/issues/41
	
	
	
	Available analyses:
                      TIGRFAM (15.0) : TIGRFAMs are protein families based on hidden Markov models (HMMs).
                      Phobius (1.01) : A combined transmembrane topology and signal peptide predictor.
                         SFLD (4) : SFLD is a database of protein families based on hidden Markov models (HMMs).
                  SUPERFAMILY (1.75) : SUPERFAMILY is a database of structural and functional annotations for all proteins and genomes.
                      PANTHER (15.0) : The PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification System is a unique resource that classifies genes by their functions, using publish
ed scientific experimental evidence and evolutionary relationships to predict function even in the absence of direct experimental evidence.
                       Gene3D (4.3.0) : Structural assignment for whole genes and genomes using the CATH domain structure database.
                        Hamap (2020_05) : High-quality Automated and Manual Annotation of Microbial Proteomes.
                        Coils (2.2.1) : Prediction of coiled coil regions in proteins.
              ProSiteProfiles (2021_01) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them.
                        SMART (7.1) : SMART allows the identification and analysis of domain architectures based on hidden Markov models (HMMs).
                          CDD (3.18) : CDD predicts protein domains and families based on a collection of well-annotated multiple sequence alignment models.
                       PRINTS (42.0) : A compendium of protein fingerprints - a fingerprint is a group of conserved motifs used to characterise a protein family.
                        PIRSR (2021_02) : PIRSR is a database of protein families based on hidden Markov models (HMMs) and Site Rules.
              ProSitePatterns (2021_01) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them.
                         Pfam (33.1) : A large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs).
                   MobiDBLite (2.0) : Prediction of intrinsically disordered regions in proteins.
                        PIRSF (3.10) : The PIRSF concept is used as a guiding principle to provide comprehensive and non-overlapping clustering of UniProtKB sequences into a hierarchical order to reflect their evolutionary relationships.

Deactivated analyses:
        SignalP_GRAM_NEGATIVE (4.1) : Analysis SignalP_GRAM_NEGATIVE is deactivated, because the resources expected at the following paths do not exist: bin/signalp/4.1/signalp
                        TMHMM (2.0c) : Analysis TMHMM is deactivated, because the resources expected at the following paths do not exist: bin/tmhmm/2.0c/decodeanhmm, data/tmhmm/2.0c/TMHMM2.0c.model
                  SignalP_EUK (4.1) : Analysis SignalP_EUK is deactivated, because the resources expected at the following paths do not exist: bin/signalp/4.1/signalp
        SignalP_GRAM_POSITIVE (4.1) : Analysis SignalP_GRAM_POSITIVE is deactivated, because the resources expected at the following paths do not exist: bin/signalp/4.1/signalp
        
        
        The following packages will be UPDATED:                                                                                                                                                                    
                                                                                                                                                                                                           
  ca-certificates                      2020.12.5-ha878542_0 --> 2021.5.30-ha878542_0                                                                                                                       
  certifi                          2020.12.5-py36h5fab9bb_1 --> 2021.5.30-py36h5fab9bb_0                                                                                                                   
  openssl                                 1.1.1i-h7f98852_0 --> 1.1.1k-h7f98852_0
	
////




echo
${InterProDir}/interproscan.sh --version
echo
echo
# Run interproscan for functional annotation
${InterProDir}/interproscan.sh --input ${Input} \
                               --output-file-base ${OutputFileBase} \
				--cpu $SLURM_NTASKS \
				--formats TSV,GFF3,JSON,XML \
				--goterms \
				--iprlookup \
				--pathways \
				--seqtype "p" \
				--disable-precalc \
				--verbose \
				--applications CDD,Coils,Gene3D,Hamap,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,PIRSR,ProSitePatterns,ProSiteProfiles,SFLD,SMART,SUPERFAMILY,TIGRFAM,Phobius

echo
# Move the Interpro output to a new directory, in case we delete the files in the future, to be able to save this directory at least.
[[ ! -d ${InterproOutputDir} ]] && mkdir -p ${InterproOutputDir}
mv ${OutputFileBase}* ${InterproOutputDir}
cp ${InterproOutputDir}/${OutputFileBase}.tsv .


# From the above list we have excluded SignalP and TMHMM, which we are going to run downstream in our analysis.
# The ProDom (2006.1) analysis has been deprecated
# HTML format is deprecated and will be removed in the second quarter of 2021. Instead, you can choose to use a JSON output and generate a graphical output in PNG, PDF, etc


echo
echo "---------------------------INTERPROSCAN ANALYSIS COMPLETED!!-------------------------------------------"
echo
echo 

echo
echo "Anvio version:"
echo
anvi-self-test --version
echo
echo


echo "			Import Pfam database search results into Anvio"
echo

 
 # The following command will search the Pfam database and add a source named "Pfam" into Anvio
 anvi-run-pfams --contigs-db ${CntgDB} --num-threads $SLURM_NTASKS

echo
echo
echo "-----------------------------------PFAM ANALYSIS COMPLETED!----------------------------------------------------"
echo
echo

echo '		Import NCBI COGs annotations into anvio'
# Instructions can be found here http://merenlab.org/2016/06/22/anvio-tutorial-v2/ and http://merenlab.org/2016/10/25/cog-annotation/
# NCBI’s COG database is quite outdated and not maintained but still awesome. The release is of December 2014! 
echo

# Prerequisite for this command to run is to have completed the anvi-setup-ncbi-cogs first!
anvi-run-ncbi-cogs --contigs-db ${CntgDB} --num-threads $SLURM_NTASKS --search-with blastp

echo
echo '----------------------------COG ANNOTATION COMPLETED---------------------------------------'
echo
echo
echo
echo "		Import Kaiju taxonomy into contigs Database"

# Instructions can be found here http://merenlab.org/2016/06/18/importing-taxonomy/
echo
# Fist get the nucleotide gene calls
${anviGetSeq} -c ${CntgDB} -o ${StrX}_gene_calls.fna
# File ${StrX}_gene_calls.fna is nucleotide sequences and for this reason Anvio exports every gene sequence
# so the good news is that we can use this file without any further analysis!!!
echo
echo '---------------------Nucleotide gene calls fasta file was created------------------------------'
echo
echo
echo "			Kaiju analysis starts now!!!"
echo
echo
# Instructions on how to run kaiju can be found here https://github.com/bioinformatics-centre/kaiju/blob/master/README.md
# IMPORTANT!!
# The following commands need the Refseq database for kaiju.
# We have downloaded the whole NCBI Refseq. This procedure takes a lot of time and space!!!

# Run the Kaiju classifier
kaiju -t ${KaijuDB}/nodes.dmp \
      			-f ${KaijuDB}/refseq/kaiju_db_refseq.fmi \
      			-i ${StrX}_gene_calls.fna \
      			-o ${StrX}_${Assembler}_gene_calls_RefSeq.tab \
      			-z $SLURM_NTASKS \
      			-v

# Next, add taxon names 
kaiju-addTaxonNames -t ${KaijuDB}/nodes.dmp \
              		-n ${KaijuDB}/names.dmp \
              		-i ${StrX}_${Assembler}_gene_calls_RefSeq.tab \
              		-o ${StrX}_${Assembler}_gene_calls_RefSeq.names \
              		-r superkingdom,phylum,order,class,family,genus,species

# At this point its not a bad idea to make a copy of the contigs database –just in case
cp ${CntgDB} ${CntgDB}.bak

# Finally, run the anvio parser for Kaiju
${anviImpTax} --input-files ${StrX}_${Assembler}_gene_calls_RefSeq.names \
                       --contigs-db ${CntgDB} \
                       --parser kaiju \
                       --just-do-it

echo
echo '------------------------KAIJU TAXONOMOMY ANALYSIS COMPLETED------------------------------------'
echo
echo
echo
echo '		Parse HMM hits to the contigs database'
echo
# Check this tutorial:
# http://merenlab.org/2015/06/25/screening-cultivars/ "Removing contaminants from cultivars with anvi'o"


anvi-run-hmms --num-threads $SLURM_NTASKS -c ${CntgDB} -H ${BuscoDBDir}
echo
echo
echo
# The following command will use the default anvio bacterial single copy gene database.
anvi-run-hmms --contigs-db ${CntgDB} --num-threads $SLURM_NTASKS

echo
echo '-------------------------------HMM HITS ANALYSIS COMPLETE-------------------------------------'
echo
echo
echo
echo '		How many genomes?'
echo

anvi-script-gen_stats_for_single_copy_genes.py ${CntgDB}
anvi-script-gen_stats_for_single_copy_genes.R ${CntgDB}.hits ${CntgDB}.genes --output_prefix=${StrX}_single_copy_genes
# The result of the above command is a pdf for later viewing!
echo
echo
echo "--------------------------SINGLE COPY GENES ANALYSIS IS COMPLETE---------------------------------"
echo
echo
echo
echo '		Contigs Statistics'
echo

# Shows simple stats of the contigs database to assess your assembly output, and also estimate the number of bacterial genomes to recover
anvi-display-contigs-stats --report-as-text --output-file ${StrX}_contig_stats.tab ${CntgDB}

echo
echo
echo '--------------------------------CONTIGS STATISTICS ANALYSIS COMPLETE----------------------------------'
echo
echo
echo
echo "			Export PANTHER annotations for online analysis"
echo
<< ////
	NOTE
	The sequence identifiers (gene calls) have been mapped to PANTHER HMM IDs by InterproScan and can be used in the gene list analysis tools.
	Panther Generic Mapping File
	For IDs from organisms other than the 82 organisms in the PANTHER database, user-generated data containing mappings between those IDs and their corresponding PANTHER IDs can be used. The file must be tab-delimited and must contain the following columns: the first column can contain a list of unique IDs from the user; the second column should be the corresponding PANTHER family or subfamily ID (e.g., PTHR10078 or PTHR10078:SF6), and is used to look up the association with GO and PANTHER terms (molecular function, biological process and pathway).
	from: https://www.nature.com/articles/nprot.2013.092

////


# We export the Panther annotation from Interproscan results and we also make the Panther IDs unique.
awk -F'[\t]' '$4 ~ /PANTHER/ {print $1"\t"$5}' ${OutputFileBase}.tsv | sort -k2 -u > ${StrX}_${Assembler}_Panther.tab


# FURTHER STEPS NEEDED!

echo
echo
echo "----------------------------------------------PANTHER FILE READY FOR UPLOAD!-------------------------------------------------------" 
echo
echo
echo " 			Export PIRSF annotations for online analysis"
echo

# We first export the PIRSF annotation from Interproscan results and we also make the PIRSF IDs unique.
awk '$4 ~ /PIRSF/{print $5}' ${OutputFileBase}.tsv | sort -u > ${StrX}_${Assembler}_PIRSF_InterproScan_IDS.txt

echo
echo "----------------------------------------------PIRSF FILE READY FOR UPLOAD!-------------------------------------------------------" 
echo
<< ////
	FURTHER STEPS NEEDED!

	FOR PANTHER
	Now it is time to go online to http://www.pantherdb.org/ and upload the file "${StrX}_${Assembler}_Panther.tab" as Panther Generic Mapping File. Please, do not forget to change the Reference List from Human to Pseudomonas aeruginosa.

	FOR PIRSF
	Now it is time to go online to https://pir.georgetown.edu/pirwww/search/batch_sf.shtml and copy and paste the IDs found in the file "${StrX}_${Assembler}_PIRSF_InterproScan_IDS.txt". Then, click on the "Save Results As: Table", download the output and rename the file to "${StrX}_${Assembler}_PIRSF.tsv". Upload "${StrX}_${Assembler}_PIRSF.tsv" to the correct folder on the server, and we are ready to go!

	NOTE:
	Commands if we want to make contamination graphs based on the hmm profiles:
	anvi-script-gen_stats_for_single_copy_genes.py ${StrX}_Pilon_contigs.db # This will create the files ${StrX}_Pilon_contigs.db.hits and ${StrX}_Pilon_contigs.db.genes 
	anvi-script-gen_stats_for_single_copy_genes.R ${StrX}_Pilon_contigs.db.hits ${StrX}_Pilon_contigs.db.genes


////
echo
echo
# Remove files no longer needed!
# rm ${StrX}_InterProScan_columns_12456_filtered.tsv ${StrX}_InterProScanVariousDatabasesImportable.tsv ${StrX}_InterProScan_columns_1213_filtered.tsv AnvioInterProAccessionsImportable.tsv ${StrX}_InterProScan_columns_14_filtered.tsv AnvioGoTermsImportable.tsv ${StrX}_InterProScan_columns_15_filtered.tsv AnvioInterProPathwaysImportable.tsv ko00001.keg KO_Orthology_ko00001.txt KeggOrthology_Table1.txt ${StrX}_KeggAnnotations-AnviImportable.txt


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
