#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=126675
# #SBATCH --mem-per-cpu=6400
#SBATCH --job-name="5th_Anvio_Prokka"
#SBATCH --output=5th_Anvio_Prokka_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END

generalInfo () {
    cat <<END

	This is the 5th part of the Anvio Prokka pipeline. 
	
	This script takes as input one argument. 
	For Strain 01 that would be: 
	sbatch ${ScriptName} Strain01

	IMPORTANT BEFORE RUNNING THIS SCRIPT:
	For gfftools:
	We downloaded gfftools with the command:
	git clone https://github.com/ihh/gfftools.git
	
	For emapper:
	We download emapper with the command:
	git clone https://github.com/eggnogdb/eggnog-mapper.git
	We created a conda environment with python3.7 and downloaded all requirements.
	We downloaded the databases in the directory ${HOME}/Software/eggnog-mapper/Database
	with the commands:
	download_eggnog_data.py -H -d 91061 -y -P -M --data_dir ${HOME}/Software/eggnog-mapper/Database
	create_dbs.py -y -m diamond --dbname bacteria --taxa Bacteria

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
    generalInfo >&2
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


#AntiSmashDir=${HOME}/Software/antismash-4.2.0
PfamDir=${HOME}/Software/Pfam32.0
OutputAntiSmashDir=${HOME}/Titlos_ktisis/Anvio/${Prefix}/Antismash
OutputEmapperDir=${HOME}/Titlos_ktisis/Anvio/${Prefix}/Emapper
AnvioDir=${HOME}/Titlos_ktisis/Anvio/${Prefix}

Genbank2GffDir=${HOME}/Software/gfftools
EmapperDir=${HOME}/Software/eggnog-mapper


EggNOGParseKEGGDir=${HOME}/Software/GhostKoalaParser  # I have slightly modified the KEGG-to-anvio executable.

AnvioGffONLYProkka=${StrX}_ONLY_PROKKA_amino_acid_seqs.gff

# Emapper parameters that can be changed
TaxScope=gproNOG # This is to limit the search on a specific database e.g. gamma Proteobacteria.
MaxSeqLen=5000
Eval=0.001


. "${HOME}/Software/miniconda3/etc/profile.d/conda.sh"
export LC_ALL=en_US.UTF-8


# Create output directory.
[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}
cd ${OutputDir}

if [[ ! -d ${OutputAntiSmashDir} ]]; then
    mkdir -p ${OutputAntiSmashDir}
fi

if [[ ! -d ${OutputEmapperDir} ]]; then
    mkdir -p ${OutputEmapperDir}
fi


#--------------------------SANITY CHECKS -------------------------------------------------------------------#
# Check that the Contigs file already exists in the working folder.
if [[ ! -s ${ContigDir}/$Contig ]]; then
    echo "The contigs fasta file $Contig is missng from ${ContigDir}. " >&2
    exit 1
fi

if [[ ! -s ko00001.keg ]]; then 
    wget -N -O ko00001.keg "https://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=htext"
fi

# Check that KEGG-to-anvio exists.
if [[ ! -s ${EggNOGParseKEGGDir}/KEGG-to-anvio ]]; then 
    echo "kegg-to-anvio has not been installed properly.  Check the directory ${EggNOGParseKEGGDir}." >&2
    exit 1
fi

# Check that the gff annotation file with only Prokka genes exists
if [[ ! -s $AnvioGffONLYProkka ]]; then
    echo "The Anvio gff file containing only prokka gene calls is missing!" >&2
    exit 1
fi


# Check that gfftools have been installed.
if [[ ! -s ${Genbank2GffDir}/genbank2gff.pl ]]; then
    echo "Gfftools have not been installed. Check the folder ${Genbank2GffDir}" >&2
    exit 1
fi


echo
#echo "		Antismash Analysis begins"
#echo

#conda activate antismash

## Create a soft link of the contigs fasta file (This file should already exist, from previous parts of the Anvio pipeline!)
#[[ ! -s ${ContigDir}/${Contig} ]] && ln -s ${ContigDir}/${Contig}

## IMPORTANT
## We are going to use the gff file which has the gene caller ids, only from Prokka CDS! We have already created this file in the first part of the pipeline.

#echo
#echo " 				Start Antismash Analysis!!!"
#echo

#echo "Antismash verstion"
#antismash --version
#echo
#echo

#antismash 	--cpus $SLURM_NTASKS \
#		--verbose \
#		--debug \
#		--taxon bacteria \
#		--minlength 200 \
#		--genefinding-gff3 ${AnvioGffONLYProkka} \
#		--output-dir ${OutputAntiSmashDir} \
#		--output-basename ${Prefix}_Antismash \
#		--html-title ${Prefix}_Antismash \
#		--html-description "${Prefix} Antismash analysis" \
#		${Contig}
#		
#		# Options not used		 
#		# --executable-paths 
#		# --cassis  # CASSIS cluster prediction only works for fungal sequences.
##		--fullhmmer \
##		--clusterhmmer \
##		--tigrfam \
##		--smcog-trees \
##		--cb-general \
##		--cb-subcluster \
##		--cb-knownclusters \
##		--asf \
##		--pfam2go \
##		--rre \
##		--cc-mibig \
#				 
#				 
#conda deactivate				 
#echo
#echo
#echo "---------------------------------ANTISMASH ANALYSIS COMPLETED!-------------------------------------------------"
#echo
#echo
#echo
#echo "				Import Antismash SmCOG Analysis into Anvio"
#echo

#conda activate anvio-7


#<< ////
#	Secondary metabolism gene families (smCOGs) analysis attempts to allocate each gene in the detected gene clusters to a secondary 
#	metabolism-specific gene family using profile hidden Markov models specific for the conserved sequence region characteristic of this family.
#////

#${Genbank2GffDir}/genbank2gff.pl ${Prefix}_Antismash.gbk > ${StrX}_Antismash.gff

#sed -E -i "s/[\']//g;s/\(Score//g;s/\)\;/;/g;s/gene=/=/g;s/NADH\://g" ${StrX}_Antismash.gff

#printf '%s\t%s\t%s\t%s\t%s\n' 'gene_callers_id' 'source' 'accession' 'function' 'e_value' > AnvioSMCOGsImportable.tsv

#awk -F'[\t=:;]' '$3 ~ /CDS/ && $12 ~ /smCOG/ {print $10"\t"$12"\t"$13"\t"$14"\t"0}' ${StrX}_Antismash.gff >> AnvioSMCOGsImportable.tsv

#mv AnvioSMCOGsImportable.tsv ${AnvioDir}

#cd ${AnvioDir}

#anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioSMCOGsImportable.tsv

#conda deactivate

#echo
#echo "----------------------------------ANTISMASH SmCOG RESULTS PARSED INTO ANVIO!!--------------------------------------------"
echo
echo
echo
echo "		Eggnog Mapper Analysis begins!!"
echo


conda activate emapper
export EGGNOG_DATA_DIR=${HOME}/Software/eggnog-mapper/Database
export PATH="${HOME}/Software/Diamond:$PATH"
export PATH="${HOME}/Software/ncbi-blast-2.10.1+/bin:$PATH"


export PATH="${HOME}/Software/hmmer-3.3.2/bin:$PATH"

echo
echo "empapper version:"		
emapper.py --version
echo
[[ ! -d ${OutputEmapperDir}/temp ]] && mkdir -p ${OutputEmapperDir}/temp


## DIAMOND SEARCH IS NOT USED BECAUSE WE USE HMMER SEARCH INSTEAD!
### diamond search
##emapper.py	-i ${StrX}_InterProScan_amino-acid-sequences.fasta \
##		--itype proteins \
##		--evalue 0.001 \
##		-m diamond \
##		--dmnd_db ${HOME}/Software/eggnog-mapper/Database/bacteria.dmnd \
##		--sensmode ultra-sensitive \
##		--matrix BLOSUM62 \
##		--seed_ortholog_evalue 0.001 \
##		--tax_scope Bacilli \
##		--tax_scope_mode inner_narrowest \
##		--go_evidence all \
##		--output ${Prefix}_diamond \
##		--output_dir ${OutputEmapperDir} \
##		--temp_dir ${OutputEmapperDir}/temp \
##		--cpu $SLURM_NTASKS 
#		

# hmmer search
emapper.py 	-i ${StrX}_InterProScan_amino-acid-sequences.fasta \
		--itype proteins \
		--evalue 0.001 \
		-m hmmer \
		--dbtype hmmdb \
		--database 91061 \
		--clean_overlaps hmmsearch_clans \
		--seed_ortholog_evalue 0.001 \
		--tax_scope Bacilli \
		--tax_scope_mode inner_narrowest \
		--go_evidence all \
		--output ${Prefix}_hmmer \
		--output_dir ${OutputEmapperDir} \
		--temp_dir ${OutputEmapperDir}/temp \
		--cpu $SLURM_NTASKS 

## NOT USED
##		--override
##		--data_dir ${HOME}/Software/eggnog-mapper/Database \

conda deactivate

#echo
#echo "------------------------------EGGNOG MAPPER ANALYSIS COMPLETED!---------------------------------"
echo
echo "			Import Eggnog Annotations into Anvio"
echo

conda activate anvio-7

# Create from the downloaded ko00001.keg file a tab delimited file 
# where the first column corresponds to the broadest classification, 
# the fifth corresponds to the gene itself

{
kegfile="ko00001.keg"
while read -r prefix content
do
    case "$prefix" in A) col1="$content";; 
                      B) col2="$content" ;; 
                      C) col3="$content";; 
                      D) echo -e "$col1\t$col2\t$col3\t$content";;
    esac 
done < <(sed '/^[#!+]/d;s/<[^>]*>//g;s/^./& /' < "$kegfile") > KO_Orthology_ko00001.txt
}

echo


printf '%s\t%s\t%s\t%s\t%s\n' 'gene_callers_id' 'source' 'accession' 'function' 'e_value' > AnvioEGGNOGImportable.tsv


grep -v '^#' ${OutputEmapperDir}/${Prefix}_hmmer.emapper.annotations | awk  -F"\t" '{print $1"\t""EGGNOG""\t"$2"\t"$8"; "$9"\t"$3}' | sed -E 's/\b[0-9]+\.//g' >> AnvioEGGNOGImportable.tsv 

anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioEGGNOGImportable.tsv 



grep -v '^#' ${OutputEmapperDir}/${Prefix}_hmmer.emapper.annotations | awk -F"\t" '{print "genecall_"$1"\t"$12}' | sed 's/ko://g;s/\t-/\t/g' > ${StrX}_EggNOGKEGGParserInput.txt

# Fist format the KEGG results into a format that can be parsed into anvio
${EggNOGParseKEGGDir}/KEGG-to-anvio --KeggDB KO_Orthology_ko00001.txt \
					-i ${StrX}_EggNOGKEGGParserInput.txt \
					-o ${StrX}_EggNOGKEGGAnnotations-AnviImportable.tsv


# Second import the results to anvio
anvi-import-functions --contigs-db ${CntgDB} --input-files ${StrX}_EggNOGKEGGAnnotations-AnviImportable.tsv

rm KO_Orthology_ko00001.txt AnvioEGGNOGImportable.tsv ${StrX}_EggNOGKEGGParserInput.txt


conda deactivate
echo
echo "----------------------------------EGGNOG ANNOTATIONS PARSED INTO ANVIO!------------------------------------------------------------"
echo
# Depracated!
#anvi-script-run-eggnog-mapper --contigs-db ${CntgDB} \
#				--annotation ${StrX}_diamond.emapper.annotations \
#				--use-version 1.0.3\
#				--num-threads $SLURM_NTASKS

<< ////
	For a manual download of the eggonog database please do the following:
	1) Download all HMM models for your target taxonomic level. 
	For instance, if Bacteria, download them from http://eggnogdb.embl.de/download/eggnog_4.5/data/bactNOG/
	2) Build a HMMER database using hmmpress. 
	For this, you will need to concatenate all models in a single file (i.e. cat bactNOG_hmm/*.hmm > bactDB.hmmer), 
	and run hmmpress bactDB.hmmer.
	3) Use hmmscan to query your sequences against the database.For instance: hmmscan bactDB.hmmer MyQueryFasta.fa
	
	More specifically:
	The following commands should run only once!
	Download the gammaproteobacteria database, concatenate the hmm files and run hmmpress!	
	wget http://eggnogdb.embl.de/download/eggnog_4.5/data/gproNOG/gproNOG.hmm.tar.gz
	tar xvf gproNOG.hmm.tar.gz
	cd gproNOG_hmm
	cat *.hmm > gproNOG_DB.hmmer
	hmmpress gproNOG_DB.hmmer
	hmmscan -o eggnog_${StrainX}_hmmscan.output \
					--tblout eggnog_${StrainX}_hmmscan_per_sequence.table \
					--domtblout eggnog_${StrainX}_hmmscan_per_domain.table \
					--pfamtblout eggnog_${StrainX}_hmmscan_pfam-like.table \
					--noali \
					-E 0.0001 \
					--domE 0.0001 --cpu $SLURM_NTASKS \
					${GProteoNOGDB}/gproNOG_DB.hmmer \
					${StrX}_InterProScan_amino-acid-sequences.fasta




wget http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/91061/91061_hmms.tar.gz
tar xvf 91061_hmms.tar.gz
cd 91061
cat *.hmm > 91061.hmm
/home1/pbravakos/Software/hmmer-3.3.2/bin/hmmpress 91061.hmm
/home1/pbravakos/Software/hmmer-3.3.2/bin/hmmscan -o eggnog_Strain28_hmmscan.output \
						--tblout eggnog_Strain28_hmmscan_per_sequence.table \
						--domtblout eggnog_Strain28_hmmscan_per_domain.table \
						--pfamtblout eggnog_Strain28_hmmscan_pfam-like.table \
						--noali \
						-E 0.0001 \
						--domE 0.0001 \
						--cpu 1 \
						91061.hmm \
						/home1/pbravakos/Titlos_ktisis/Anvio/Str28_Pilon/Str28_InterProScan_amino-acid-sequences.fasta

////

echo
echo "======================================================================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
elapsed=$SECONDS
echo Time taken: `printf '%dd %dh:%dm:%ds\n' $((elapsed/86400)) $((elapsed%86400/3600)) $((elapsed%3600/60)) $((elapsed%60))`

exit 0
