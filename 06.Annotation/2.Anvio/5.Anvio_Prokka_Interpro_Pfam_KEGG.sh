#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="4th_Anvio_Prokka"
#SBATCH --output=4th_Anvio_Prokka_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


generalInfo () {
    cat <<END

	This is the Anvio pipeline. PART FOUR
	Here we try to parse into anvio, most of the databases from the output of Interpro, but also run independently TMHMM, signalP, LipoP.

	This script takes as input one argument. 
	For Strain 01 that would be: 
	sbatch ${ScriptName} Strain01

	IMPORTANT!
	Instead of running GhostKOALA which is mostly for metagenomes we prefered to run BlastKOALA, although not mentioned in the tutorial, but the procedure is exactly the same!
	We have run the online BlastKOALA https://www.kegg.jp/blastkoala/, downloaded the annotation file, renamed it to "${StrX}_${Assembler}_user_ko.txt" and uploaded it to the Anvio folder on the server.
	IMPORTANT!
	We have already run the Panther Generic Mapping online, at http://www.pantherdb.org/, downloaded pantherGeneList.txt, renamed it to ${StrX}_${Assembler}_Panther_Bacilus_subtilis.tsv.
	We also have uploaded the PIRSF annotations on the server with the name ${StrX}_PIRSF.tsv
	
	IMPORTANT!
	Download MetaCyc database from 
	Got to https://biocyc.org/group?id=:ALL-PATHWAYS&orgid=META# 
	Choose export to spreadsheet file
	Then choose "frame IDs" and rename this file to ${MetaCycDBIDs}
	Repeat the same process but the second time
	Choose "Common names" and rename this file to ${MetaCycDBAnnot}

END
}

echo



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


InterproTSVOut=${Prefix}_interpro-output.tsv

TMHMMDir=${HOME}/Software/tmhmm-2.0c/bin
LipoPDir=${HOME}/Software/LipoP1.0a
MetaCycDBDir=${HOME}/Titlos_ktisis/MetaCyc
MetaCycDBIDs=MetaCyc_IDs.txt
MetaCycDBAnnot=MetaCyc_annotations.txt


# SignalP parameters that can change
Gram=gram+ #options: euk, gram+, gram-
Format=short #options: 'short', 'long', 'summary' or 'all'



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
export PATH="${HOME}/Software/signalp/bin:$PATH"


#------------------------------------------------------------------------------------------------------------
# Create output directory.
[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}
cd ${OutputDir}


#--------------------------SANITY CHECKS -------------------------------------------------------------------#

# Check that the interpo output tsv file exists
if [[ ! -s ${OutputDir}/$InterproTSVOut ]]; then
   echo " The Interpo output file $InterproTSVOut is missing from directory ${OutputDir}" >&2
   echo >&2
   generalInfo >&2
fi

# Check whether the programs needed are installed
if [[ ! -r "${TMHMMDir}/tmhmm"  && ! -r "${SignalPDir}/signalp" && ! -r "${LipoPDir}/LipoP" ]]; then
    echo "Missing executables (TMHMM, SignalP or LipoP)!" >&2
    echo >&2
    generalInfo >&2
    exit 1
fi

## Check that we have uploaded the files for Panther and PIRSF annotations.
#if [[ ! -s ${OutputDir}/${Prefix}_pantherGeneList.txt ]] &&  [[ ! -s ${OutputDir}/${Prefix}_PIRSF.tsv ]]; then 
#    echo "Please upload the missing file(s), ${Prefix}_pantherGeneList.txt or ${Prefix}_PIRSF.tsv!" >&2
#    exit 1
#fi

# Check whether we have uploaded the MetaCyc Database, from the website https://biocyc.org/group?id=:ALL-PATHWAYS&orgid=META# 
if [[ ! -s $MetaCycDBDir/$MetaCycDBIDs && ! -s $MetaCycDBDir/$MetaCycDBAnnot ]]; then 
    echo "MetaCyc database files have not been uploaded in the correct folder $MetaCycDBDir. Please upload them!" >&2
    echo >&2
    generalInfo >&2
    exit 1
fi

# Check whether the files we need exist and if not actually download them!
if [[ ! -s br08901.keg ]]; then
    wget -N -O br08901.keg "https://www.genome.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext"
fi

if [[ ! -s ko01000.keg ]]; then
    wget -N -O ko01000.keg "https://www.genome.jp/kegg-bin/download_htext?htext=ko01000.keg&format=htext"
fi

if [[ ! -e descriptions.pl ]]; then 
    wget "http://smart.embl-heidelberg.de/smart/descriptions.pl"
fi

if [[ ! -e scop.annotation.1.73.txt ]]; then
    wget "http://supfam.org/SUPERFAMILY/function/scop.annotation.1.73.txt"
fi

#if [[ ! -e $SLURM_SUBMIT_DIR/cath-b-newest-names.gz ]]; then
#    wget "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/daily-release/newest/cath-b-newest-names.gz"
#fi

if [[ ! -e cath-b-newest-names.gz ]]; then
    wget "http://download.cathdb.info/cath/releases/daily-release/newest/cath-b-newest-names.gz"
fi

if [[ ! -e ReactomePathways.txt ]]; then
    wget https://reactome.org/download/current/ReactomePathways.txt
fi

# Next we apply a simple grep, with search patterns that normally should exist in the files of interest.
if  ! grep --quiet -E 'SM[]0-9]{5}' descriptions.pl  &&  ! grep --quiet -P "\tsf\t" scop.annotation.1.73.txt  &&  ! zgrep --quiet "Mainly Alpha" cath-b-newest-names.gz && ! grep -P --quiet "R-ATH-[0-9]+\t" ReactomePathways.txt; then
    echo "Some or all of the annotation files needed for the analyses (SMART, SUPERFAMILY, Gene3D, KEGG and Reactome) are probably problematic. Please check the corresponding files and/or the websites" >&2
    exit 1
fi

#--------------------------SANITY CHECKS END! -------------------------------------------------------------------#



echo
echo
echo '		Import Interproscan annotations into anvio'
echo
echo
<< ////
	Next we will import functions from the Interpro tsv file manually to the Contig Database!
	Because Anvio does not keep track of the start and stop codons related to the annotated regions inside the genes, it is expected that some of the annotations will refer more than once to the same gene, even after the filtering we apply bellow. Just keep in mind that because of this, and because the graphical representations of genes in the Anvio interface can show only one annotation per source per gene, it is expected that the search results for a specific term will return annotations that cannot be found at the graphical representation! Anvio in theory chooses to represent only the gene annotation with the best e-value, in case of a conflict. 
	Check also here:
	https://github.com/merenlab/anvio/issues/999
	
	HOW we extract the tsv interpo output info
	Fist we export the columns we are interested in, delete lines with empty fields, remove dashes (-) from the e-value column. Also replace spaces with @@ (I stil do not understand why i have to do this, but otherwise the awk filtering applied at the end does not work!). We also remove Panther, PIRSF, SMART, SUPERFAMILY, Gene3D because we are going to parse these annotations later on. We also remove all fields (lines) that do not have an e-value asigned to the annotation. We sort by the evalues, having the lowest e-value first and then we apply the awk filtering which removes any duplicates, that have the same gene caller id and source. So as a result of the filtering we keep only the lines with the lowest e-values. Finally we append all the results that did not have an e-value, after filtering out duplicate records, and replacing the e-value with a random high number e.g 5000, in order to help anvio always select to present the sources that do have an assigned e-value, in case of a conflict.
	One of the Databases searched by Interpro is Pfam. In the current version of Interpro the Pfam database is version 33.1. On the other hand Anvio has its own function to search Pfam and uses the lastest database which is version 32. We want to have both results, so we have to rename the source "Pfam"" from Interpro to something else, because otherwise one result will overwrite the other one. 
	NOTE
	Because the filtering we apply is extensive, we also have to keep the interpo output files, in case we really want to get all the information associated to a specific annotation.
////

awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$4,$5,$6,$9}' ${InterproTSVOut} | \
sed -E 's/-$//g;s/ /@@/g'| \
awk '$3!="" && $4!="" && $2!="Gene3D" && $2 !="Pfam" && $2 !="SMART" && $5!=""' | \
sort -k5g -t$'\t'| \
awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} !a[$1$2]++' > ${StrX}_InterProScanVariousDatabasesImportable.tsv


awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$4,$5,$6,$9}' ${InterproTSVOut} | \
sed -E 's/-$//g;s/ /@@/g'| \
awk '$3!="" && $4!="" && $2!="Gene3D" && $2 !="Pfam" && $2 !="SMART" && $5==""' | \
sort -g -k5,1 -t$'\t' | \
awk 'BEGIN{FS="\t";OFS="\t"} !a[$1$4]++ {print $1,$2,$3,$4,"5000"}' >> ${StrX}_InterProScanVariousDatabasesImportable.tsv

sed -i 's/@@/ /g' ${StrX}_InterProScanVariousDatabasesImportable.tsv
#;s/Pfam/PfamInterPro/g  # Removed from the above sed command!

anvi-import-functions --contigs-db ${CntgDB} --input-files ${StrX}_InterProScanVariousDatabasesImportable.tsv


# Second, import Interpo annotations, by Interpro!
# We export the columns of interest from the interpo output and create a new file with the evalue column stripped of any dash (-), then we sort by the scientific numbers of the evalues and save to a new file.
awk 'BEGIN{FS="\t";OFS="\t"} $12!="-" && $13!="-" {print $1,"Interpro",$12,$13,$9}' ${InterproTSVOut} | \
sed -E '/\t-$/d' | \
sort -k5g -t$'\t' > InterpoAnnotationsFiltered.tsv

# We append to the previously created file the missing annotations, which have no evalue associated with them.
awk 'BEGIN{FS="\t";OFS="\t"} $5 ~ /-/ && $12!="-" && $13!="-" {print $1,"Interpro",$12,$13,$9}' ${InterproTSVOut} | \
sed -E 's/\-//g' >> InterpoAnnotationsFiltered.tsv

# And we create a file ready to be imported into Anvio. 
# Please notice the array filtering which leaves us with only one interpo accession per gene caller id. This filtering is aggressive, but necessary (my opinion!!).
awk 'BEGIN{FS="\t"; OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} !a[$1]++' InterpoAnnotationsFiltered.tsv  > AnvioInterProAccessionsImportable.tsv
echo
# Finally, we import the Interpo (by Interpo!) annotations into Anvio!
anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioInterProAccessionsImportable.tsv
echo

# Third, remove all duplicates again on the initial file for importing GoTerm annotations
awk -F'\t' '!a[$1$14]++' ${InterproTSVOut} > ${StrX}_InterProScan_columns_14_filtered.tsv

# Import GoTerms
printf '%s\t%s\t%s\t%s\t%s\n' 'gene_callers_id' 'source' 'accession' 'function' 'e_value' > AnvioGoTermsImportable.tsv

# Export the columns we are interested in, delete lines with empty fields, remove dashes (-) from the e-value column and append to file.
awk 'BEGIN{FS="\t";OFS="\t"} $14!="-" && $9!="-" {print $1,"Interpo_GoTerms","",$14,$9}' ${StrX}_InterProScan_columns_14_filtered.tsv | \
awk '$4!=""' | \
sed -e 's/-$//g' >> AnvioGoTermsImportable.tsv
echo
anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioGoTermsImportable.tsv
echo

# Finally, remove all duplicates again on the initial file for importing Pathways annotations
awk -F'\t' '!a[$1$15]++' ${InterproTSVOut} > ${StrX}_InterProScan_columns_15_filtered.tsv

# Import Pathways found by the interposcan search.
printf '%s\t%s\t%s\t%s\t%s\n' 'gene_callers_id' 'source' 'accession' 'function' 'e_value' > AnvioInterProPathwaysImportable.tsv

# Export the columns we are interested in, delete lines with empty fields, remove dashes (-) from the e-value column and append to file.
awk 'BEGIN{FS="\t";OFS="\t"} $15!="-" && $9!="-" {print $1,"PathwayInterPro","",$15,$9}' ${StrX}_InterProScan_columns_15_filtered.tsv | \
awk '$4!=""' | \
sed -e 's/-$//g' >> AnvioInterProPathwaysImportable.tsv
echo
anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioInterProPathwaysImportable.tsv
echo

rm InterpoAnnotationsFiltered.tsv

## The next step is not needed since we parse the interproscan output in the previous awk commands!!!.
# NOT NEEDED! anvi-import-functions --contigs-db ${CntgDB} --input-files ${InterproTSVOut} --parser interproscan 

echo
echo '--------------------------------INTERPROSCAN ANNOTATION IMPORTED INTO ANVIO-------------------------------------'
echo
echo
echo "				Signal Peptide Analysis Begins"
echo
echo
signalp -version	
echo
echo
echo
signalp -fasta ${StrX}_InterProScan_amino-acid-sequences.fasta \
	-format ${Format} \
	-org ${Gram} \
	-mature \
	-gff3 \
	-verbose \
	-prefix ${Prefix}_signalp

echo
echo "-------------------------------SIGNALP ANALYSIS COMPLETED!!--------------------------------------"
echo
echo
echo "		Import SignalP predictions into Anvio"
echo

printf '%s\t%s\t%s\t%s\t%s\n' 'gene_callers_id' 'source' 'accession' 'function' 'e_value' > AnvioSignalPImportable.tsv



awk -F'\t' '$3 ~ /signal/ {print $1"\t"$2"\t""-""\t"$3"\t"0}' ${Prefix}_signalp.gff3 >> AnvioSignalPImportable.tsv

anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioSignalPImportable.tsv

echo
echo "-------------------------SIGNALP PREDICTIONS WERE IMPORTED INTO ANVIO!!----------------------------------"
echo
echo
echo "				TMHMM Analysis Begins!"
echo

# We will run TMHMM twice.
# First, we will search for Transmembrane Helices on all the amino acids.
${TMHMMDir}/tmhmm ${StrX}_InterProScan_amino-acid-sequences.fasta > ${Prefix}_TMHMM_Anvio_aa.txt


# Next, we are going to run TMHMM on the fasta output of SignalP we got in the previous steps. 
# The advantage of using the processed FASTA entries instead of the entire sequences is that we get rid of the 
# false-positive transmembrane helix predictions that TMHMM often makes for SPs

${TMHMMDir}/tmhmm ${Prefix}_signalp_mature.fasta > ${Prefix}_TMHMM_signalp.txt

echo
echo "-------------------------------TMHMM ANALYSIS COMPLETED!!--------------------------------------"
echo
echo
echo "		Import TMHMM predictions into Anvio"
echo

# For the more general search we are going to keep only the results that actually gave some TMHs results. (Notice the "Number of predicted TMHs:  [1-9]" part)
egrep -o "[0-9]{1,4} Number of predicted TMHs:  [1-9]{1,3}" ${Prefix}_TMHMM_Anvio_aa.txt |sed 's/  / /g' | awk '{print$1"\t"$2,$3,$4,$5,$6}'| awk -F"\t" '{print$1"\t""TMHMM2.0""\t""-""\t"$2"\t""0"}' > AnvioTMHMM_temporary.tsv

# For the more specific results, from the SignalP input, we are going to keep everything, either we got a result or not (Notice the "Number of predicted TMHs:  [0-9]" part). This is to make sure that we always get the more accurate results, and avoid any false positives from the more general search.
egrep -o "[0-9]{1,4} Number of predicted TMHs:  [0-9]{1,3}" ${Prefix}_TMHMM_signalp.txt |sed 's/  / /g' | awk '{print$1"\t"$2,$3,$4,$5,$6}'| awk -F"\t" '{print$1"\t""TMHMM2.0""\t""-""\t"$2"\t""0"}' > AnvioTMHMMsignalP_temporary.tsv


printf '%s\t%s\t%s\t%s\t%s\n' 'gene_callers_id' 'source' 'accession' 'function' 'e_value' > AnvioTMHMMImportable.tsv

# Finally, we want to merge the two results. Because the results from the signalP input are more reliable, we want to keep these and discard the results from the more generals search, when there is conflict.
# For this reason we place AnvioTMHMMsignalP_temporary.tsv before AnvioTMHMM_temporary.tsv so that the correct, newer, lines are preferred when eliminating duplicates.
TAB=`echo -e "\t"` # sort for some reason cannot accept "\t" directly, so we have to pass it to an argument.
sort -t"$TAB" -k1,1n -us AnvioTMHMMsignalP_temporary.tsv AnvioTMHMM_temporary.tsv >> AnvioTMHMMImportable.tsv


anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioTMHMMImportable.tsv

rm ${Prefix}_TMHMM_Anvio_aa.txt ${Prefix}_TMHMM_signalp.txt AnvioTMHMM_temporary.tsv AnvioTMHMMsignalP_temporary.tsv

echo
echo "-------------------------TMHMM PREDICTIONS WERE IMPORTED INTO ANVIO!!----------------------------------"
echo
echo
echo "					Import LipoP predictions into Anvio"
echo

# We are going to import the lipoproteins annotation.

# Next we run the command, format the output according to our needs, and output it in a format Anvio can accept.
${LipoPDir}/LipoP ${StrX}_InterProScan_amino-acid-sequences.fasta -short -noplot > ${StrX}_LipoP_Anvio_aa.txt

grep "cleavage=" ${StrX}_LipoP_Anvio_aa.txt | \
sed -E 's/^# //g;s/^([0-9]*) /\1\t/g;s/ score=[0-9]*\.[0-9]* margin=[0-9]*\.[0-9]* cleavage=/, cleavage AA:/g;s/ Pos/, Pos/g;s/SpII/lipoprotein signal peptide/g;s/SpI/signal peptidase I/g' | \
awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} {print $1,"LipoP","",$2,"0"}' > AnvioLipoPImportable.tsv


# Finally, we parse these results into Anvio.
anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioLipoPImportable.tsv

rm  ${StrX}_LipoP_Anvio_aa.txt

echo
echo "-------------------------LipoP PREDICTIONS WERE IMPORTED INTO ANVIO!!----------------------------------"
echo
echo
echo
##echo "		Import Panther Annotations into Anvio"
##echo 

### We have run the Panther Generic Mapping online, at http://www.pantherdb.org/, downloaded pantherGeneList.txt, renamed it to ${StrX}_pantherGeneList.txt

### First we export the Panther annotations but this time we change the order of the columns, and also export e-values. 
### In the first column we place the panther annotations, and in the second we place the gene caller ids.
##awk -F'[\t]' '$4 ~ /PANTHER/ {print $5"\t"$1"\t"$9}' ${InterproTSVOut} > ${StrX}_Panther_INVERSE.tab

##<< ////
##	Next, we append to each line of the "${StrX}_Panther_INVERSE.tab" - using awk arrays which will hold the common column into memory and print the matching values of the first file to the second - the correct Panther annotation from the downloaded gene list file "${StrX}_pantherGeneList.txt", and thus correlate the annotations to the gene caller ids we have from Anvio. 
##	We also reverse sort based on the Panther IDs, and remove lines which have the same gene caller ID and e-value. This is because we want to keep only the subfamilies of the Panther IDs, which contain more information than the families. Interpro for a reason ouputs both families and subfamilies with the same e-value for each gene caller. 
##	We could apply this filtering one step earlier, when we created the "${StrX}_Panther_INVERSE.tab" file, but we choose to do it bellow, because it is easier to understand, why we are applying this kind of filtering, because of the extra information by the Panther annotations.
##	Finally we save each column of the "${StrX}_pantherGeneList.txt" to a new file, which can be imported later into Anvio. 
##////

### Fist we will import the Panther Protein annotations into Anvio.
##awk -F"\t" 'NR==FNR{a[$1]=$7;next}{print $0,"\t"a[$1]}' ${StrX}_pantherGeneList.txt ${StrX}_Panther_INVERSE.tab | \
##sort -r | awk '!a[$2$3]++' | \
##awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} {print $2,"PANTHER Protein",$1,$4,$3}' > AnvioPANTHERProteinClassAnnotationsImportable.tsv

##anvi-import-functions --contigs-db ${CntgDB}  --input-files AnvioPANTHERProteinClassAnnotationsImportable.tsv


### Next we will import Go Terms. We will use some pattern matching and subtitution, in order to get the desired outcome.
### First Go Molecular Functions
##awk 'BEGIN{FS="\t";OFS="\t"} {print $1, $4}' ${StrX}_pantherGeneList.txt | \
##sed 's/[()]//g' | \
##sed -E 's/([a-zA-Z, -]+)(GO:[0-9]{7})/\2\t\1/g;s/;/\t/g' | \
##awk  'BEGIN{FS="\t";OFS="\t"} {print $1, "F:"$2";"" F:"$4";"" F:"$6";"" F:"$8" F:"$10" F:"$12" F:"$14,"F:"$3";"" F:"$5";"" F:"$7";"" F:"$9" F:"$11" F:"$13" F:"$15}' | sed -E 's/([ ]?F:[;]?){2,}//g;s/ (F:)+\t/\t/g;s/ (F:)+$//g' > PANTHER_GO-Slim_Molecular_Function.tsv

##awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$2;b[$1]=$3;next} {print $1,$2,$3,a[$1],b[$1]}' PANTHER_GO-Slim_Molecular_Function.tsv ${StrX}_Panther_INVERSE.tab | \
##awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} $4!="" {print $2,"PANTHER GO-Slim Function",$4,$5,$3}' | awk '!a[$1$3]++' > AnvioPANTHERFunctionalAnnotationsImportable.tsv

##anvi-import-functions --contigs-db ${CntgDB}  --input-files AnvioPANTHERFunctionalAnnotationsImportable.tsv


### Second, Go Biological Processes
##awk 'BEGIN{FS="\t";OFS="\t"} {print $1, $5}' ${StrX}_pantherGeneList.txt | \
##sed 's/[()]//g' | \
##sed -E 's/([a-zA-Z, -]+)(GO:[0-9]{7})/\2\t\1/g;s/;/\t/g' | \
##awk  'BEGIN{FS="\t";OFS="\t"} {print $1, "P:"$2";"" P:"$4";"" P:"$6";"" P:"$8" P:"$10" P:"$12" P:"$14,"P:"$3";"" P:"$5";"" P:"$7";"" P:"$9" P:"$11" P:"$13" P:"$15}' | \
##sed -E 's/([ ]?P:[;]?){2,}//g;s/ (P:)+\t/\t/g;s/ (P:)+$//g' > PANTHER_GO-Slim_Biological_Process.tsv

##awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$2;b[$1]=$3;next} {print $1,$2,$3,a[$1],b[$1]}' PANTHER_GO-Slim_Biological_Process.tsv ${StrX}_Panther_INVERSE.tab | \
##awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} $4!="" {print $2,"PANTHER GO-Slim Process",$4,$5,$3}' | awk '!a[$1$3]++' > AnvioPANTHERProcessesAnnotationsImportable.tsv

##anvi-import-functions --contigs-db ${CntgDB}  --input-files AnvioPANTHERProcessesAnnotationsImportable.tsv


### And finally, Go Cellular components.
##awk 'BEGIN{FS="\t";OFS="\t"} {print $1, $6}' ${StrX}_pantherGeneList.txt | \
##sed 's/[()]//g' | \
##sed -E 's/([a-zA-Z, -]+)(GO:[0-9]{7})/\2\t\1/g;s/;/\t/g' | \
##awk  'BEGIN{FS="\t";OFS="\t"} {print $1, "C:"$2";"" C:"$4";"" C:"$6";"" C:"$8" C:"$10" C:"$12" C:"$14,"C:"$3";"" C:"$5";"" C:"$7";"" C:"$9" C:"$11" C:"$13" C:"$15}' | \
##sed -E 's/([ ]?C:[;]?){2,}//g;s/ (C:)+\t/\t/g;s/ (C:)+$//g' > PANTHER_GO-Slim_Cellular_Component.tsv

##awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$2;b[$1]=$3;next} {print $1,$2,$3,a[$1],b[$1]}' PANTHER_GO-Slim_Cellular_Component.tsv ${StrX}_Panther_INVERSE.tab | \
##awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} $4!="" {print $2,"PANTHER GO-Slim Cellular",$4,$5,$3}' | awk '!a[$1$3]++' > AnvioPANTHERCellularAnnotationsImportable.tsv

##anvi-import-functions --contigs-db ${CntgDB}  --input-files AnvioPANTHERCellularAnnotationsImportable.tsv

##echo
##rm ${StrX}_Panther_INVERSE.tab PANTHER_GO-Slim_Molecular_Function.tsv PANTHER_GO-Slim_Biological_Process.tsv PANTHER_GO-Slim_Cellular_Component.tsv

##echo
##echo "-----------------------------------------------PANTHER ANNOTATIONS WERE IMPORTED INTO ANVIO!---------------------------------------------------------------------------------"
##echo
##echo
##echo "				Import PIRSF Annotations into Anvio"
##echo

### We have already downloaded the file from https://pir.georgetown.edu/pirwww/search/batch_sf.shtml, renamed it to ${StrX}_PIRSF.tsv and uploaded it to the server.

### First we export again the PIRSF annotations but please note the order of the columns. 
### In the first column we place the PIRSF annotations, in the second we place the gene caller ids and lastly we place the e-values.
##awk 'BEGIN{FS="\t";OFS="\t"} $4 ~ /PIRSF/ {print $5,$1,$9}' ${InterproTSVOut} > ${StrX}_PIRSF_Interpro_INVERSE.tsv

### From the file we downloaded from the internet, we remove the first line which is a header, export only the columns of interest, remove some disturbing empty spaces, rename the accession field to match the accession in the file "${StrX}_PIRSF_Interpro_INVERSE.tsv" and save to a new file.
##awk 'BEGIN{FS="\t";OFS="\t"} {if (NR!=1) {print $1,$3}}' ${StrX}_PIRSF.tsv | sed -E 's/ \t/\t/g;s/ $//g; s/^SF/PIRSF/g' > ${StrX}_PIRSF_Anotations.tsv

### Next, we compare the first columns of the two files, and when there is a match, we append the annotation of the first file to the second. Finally, we save to a new file.
##awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$2;next}{print $1,$2,$3,a[$1]}' ${StrX}_PIRSF_Anotations.tsv ${StrX}_PIRSF_Interpro_INVERSE.tsv > PIRSF_temporary_Anvio.tsv

### Create the file to be parsed into Anvio with the correct headers, and parse PIRSF functional annotation
##awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} {print $2,"PIRSF",$1,$4,$3}' PIRSF_temporary_Anvio.tsv | awk '!a[$1$4]++' | awk '!a[$1$3]++' > AnvioPIRSFFunctionalAnnotationsImportable.tsv


### Finally, parse functional annotations into Anvio!
##anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioPIRSFFunctionalAnnotationsImportable.tsv
##echo
##rm ${StrX}_PIRSF_Interpro_INVERSE.tsv ${StrX}_PIRSF_Anotations.tsv PIRSF_temporary_Anvio.tsv

##echo
##echo
##echo "-----------------------------------------------------PIRSF ANNOTATIONS WERE IMPORTED INTO ANVIO--------------------------------------------------------------------------------"
echo
echo
echo "				Import SMART Annotations into Anvio"
echo 

# First download the SMART annotation file from their website!


# From the file we downloaded from the internet, we remove the first two lines which are headers, export only the columns of interest and save it to a new file.
awk  'BEGIN{FS="\t";OFS="\t"} {if (NR > 2) {print $2,$3"; "$4}}'  descriptions.pl | sed 's/  / /g' > SMART_Annotations.tsv

# Now, we export again the SMART annotations but please note the order of the columns. 
# In the first column we place the SMART annotations, in the second we place the gene caller ids and lastly we place the e-values.
awk 'BEGIN{FS="\t";OFS="\t"} $4 ~ /SMART/ {print $5,$1,$9}' ${InterproTSVOut}  > ${StrX}_SMART_Interpro_INVERSE.tsv

# Next, we compare the first columns of the two files, and when there is a match, we append the annotation of the first file to the second. Finally, we save to a new file.
awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$2;next}{print $1,$2,$3,a[$1]}' SMART_Annotations.tsv ${StrX}_SMART_Interpro_INVERSE.tsv > SMART_temporary_Anvio.tsv

# Create the file to be parsed into Anvio with the correct headers, and functional annotations
awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} {print $2,"SMART",$1,$4,$3}' SMART_temporary_Anvio.tsv | awk '!a[$1$3]++' > AnvioSMARTFunctionalAnnotationsImportable.tsv

# Finally, parse functional annotations into Anvio!
anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioSMARTFunctionalAnnotationsImportable.tsv
echo
rm SMART_Annotations.tsv ${StrX}_SMART_Interpro_INVERSE.tsv SMART_temporary_Anvio.tsv

echo
echo
echo "-----------------------------------------------------SMART ANNOTATIONS WERE IMPORTED INTO ANVIO-----------------------------------------------------------------	---------------"
##echo
##echo
##echo "				Import SUPERFAMILY Annotations into Anvio"
##echo

### First download the SUPERFAMILY annotation file from their website!


### From the file we downloaded from the internet we export only the columns of interest, rename the accession field to match the interpo accesions and save it to a new file.
##awk  'BEGIN{FS="\t";OFS="\t"} {print $2,$6}' scop.annotation.1.73.txt | sed 's/^/SSF/g' > SUPERFAMILY_Annotations.tsv

### Now, we export the SUPERFAMILY annotations from the interpro results but please note the order of the columns. 
### In the first column we place the SUPERFAMILY annotations, in the second we place the gene caller ids and lastly we place the e-values.
##awk 'BEGIN{FS="\t";OFS="\t"} $4 ~ /SUPERFAMILY/ {print $5,$1,$9}' ${InterproTSVOut} > ${StrX}_SUPERFAMILY_Interpro_INVERSE.tsv

### Next, we compare the first columns of the two files, and when there is a match, we append the annotation of the first file to the second. We also remove some disturbing trailing white spaces! Finally, we save to a new file.
##awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$2;next}{print $1"\t"$2"\t"$3"\t"a[$1]}' SUPERFAMILY_Annotations.tsv ${StrX}_SUPERFAMILY_Interpro_INVERSE.tsv | sed 's/\s$//g' > SUPERFAMILY_temporary_Anvio.tsv

### Create the file to be parsed into Anvio with the correct headers, and functional annotations
##awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} {print $2,"SUPERFAMILY",$1,$4,$3}' SUPERFAMILY_temporary_Anvio.tsv | awk '!a[$1$3]++' > AnvioSUPERFAMILYFunctionalAnnotationsImportable.tsv

### Finally, parse functional annotations into Anvio!
##anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioSUPERFAMILYFunctionalAnnotationsImportable.tsv
##echo
##rm SUPERFAMILY_Annotations.tsv ${StrX}_SUPERFAMILY_Interpro_INVERSE.tsv SUPERFAMILY_temporary_Anvio.tsv

##echo
##echo
###echo "-----------------------------------------------------SUPERFAMILY ANNOTATIONS WERE IMPORTED INTO ANVIO--------------------------------------------------------------------------------"
echo
echo
echo "				Import Gene3D Annotations into Anvio"
echo

# Fist, we download the newest annotation from the cath database.
# Then we extract the contents to our working directory
gunzip cath-b-newest-names.gz

# From the file we downloaded from the internet replace the first space with a tab, rename the accession field to match the interpo accesions and save it to a new file.
sed 's/^\([0-9.]*\) /\1\t/g; s/^/G3DSA:/g' cath-b-newest-names > Gene3D_Annotations.tsv


awk 'BEGIN{FS="\t";OFS="\t"} $4 ~ /Gene3D/ {print $5,$1,$9}' ${InterproTSVOut} > ${StrX}_Gene3D_Interpro_INVERSE.tsv


awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$2;next}{print $1"\t"$2"\t"$3"\t"a[$1]}' Gene3D_Annotations.tsv ${StrX}_Gene3D_Interpro_INVERSE.tsv > Gene3D_temporary_Anvio.tsv


awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} {print $2,"Gene3D",$1,$4,$3}' Gene3D_temporary_Anvio.tsv | awk '!a[$1$4]++' > AnvioGene3DFunctionalAnnotationsImportable.tsv


# Finally, parse functional annotations into Anvio!
anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioGene3DFunctionalAnnotationsImportable.tsv
echo
rm Gene3D_Annotations.tsv ${StrX}_Gene3D_Interpro_INVERSE.tsv Gene3D_temporary_Anvio.tsv
echo
echo
echo "-----------------------------------------------------Gene3D ANNOTATIONS WERE IMPORTED INTO ANVIO--------------------------------------------------------------------------------"
echo
echo
echo "				Parse the Reactome annotations found by Interpo "
echo

# Concatenate columns 2 and 3 of the downloaded Reactome annotations file.
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2"; "$3}' ReactomePathways.txt > react.txt

# The next block of code extracts the Reactome annotation done by InterproScan and the Reactome gene annotation from the Reactome website and combines it to a new file which holds information about gene caller ids, accessions as well as gene annotations.
{
# We will export only the first 4 columns
for i in {1..4}
do
    # First we export the Reactome annotations that were found by InterproScan.
    awk 'BEGIN{FS="\t"; OFS="\t"} $15 ~ /Reactome:/ {print $1,$9,$15}' ${InterproTSVOut} | \
    sed -E 's/KEGG: [0-9]{5}[0-9.+]+\|//g;s/MetaCyc: [A-Z0-9\-]+\|//g;s/\|Reactome: /\t/g;s/Reactome://g;s/\t R\-/\tR\-/g' | \
    awk -v var=$((2+i)) 'BEGIN{FS="\t"; OFS="\t"} {print $var,$1,$2}' > ReactomeDeleteMe
    
    # Next we combine the annotations from the Reactome website with the ones from Interpro.
    awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$2;next} {print $2,$1,$3,a[$1]}' react.txt ReactomeDeleteMe > Reactome.temp
    
    # We delete lines that have no evalues and those that do not have an accesion number, and sort based on the evalue column.
    sed -E '/\t-\t/d' Reactome.temp | awk  'BEGIN{FS="\t"; OFS="\t"} $2 !=""' | sort -k3g -t$'\t' > Reactome_sorted
    # Next we append the annotations that do not have evalues, at the end of the file (We delete once more the lines without an accesion number).
    awk '$3=="-"' | awk  'BEGIN{FS="\t"; OFS="\t"} $2 !=""' Reactome.temp |sort -n >> Reactome_sorted
    
    # We filter out only the columns we need.
    awk 'BEGIN{FS="\t";OFS="\t"} !a[$1]++ {print $1,$2,$4}' Reactome_sorted > Reactome${i}_filtered.txt
    echo
    # Finally, we join the columns together, based on the gene callers ids.
    if [[ i -eq 2 ]]; then
        #join  Reactome1_filtered.txt  Reactome2_filtered.txt -t $'\t' >  ReactomeJoined1_2.tsv
        awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$2;a[$2]=$3;next} {if ($2 == a[$1]) print $1,$2,$3,$2,$3; else if ($3 == a[$2]) print $1,$2,$3,$2,$3; else print $1,$2,$3,$2"; "a[$1],$3"; "a[$2]}' Reactome1_filtered.txt Reactome2_filtered.txt | sed -E 's/; $//g' > ReactomeJoined1_2.tsv
    fi
    echo
    if [[ i -gt 2 ]]; then
        #join  ReactomeJoined$((i-2))_$((i-1)).tsv  Reactome${i}_filtered.txt -t $'\t' >  ReactomeJoined$((i-1))_${i}.tsv
        awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$2;a[$2]=$3;next} {if ($2 == a[$1]) print $0; else if ($3 == a[$2]) print $0; else print $1,$2,$3,$4"; "a[$1],$5"; "a[$2]}' Reactome${i}_filtered.txt ReactomeJoined$((i-2))_$((i-1)).tsv | sed -E 's/; $//g' >  ReactomeJoined$((i-1))_${i}.tsv
    fi
    echo
done
}

# Finally, create the anvio importable file containing all reactome accessions with their annotations and the corresponding gene callers ids.
awk 'BEGIN{FS="\t"; OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} {print $1,"Reactome",$4,$5,"0"}' ReactomeJoined3_4.tsv | sed -E 's/ \t;//g;s/\t;/;/g;s/(; ;)+//g;s/ ; /;/g;s/; $//g' > AnvioReactomeImportable.tsv

# Import the Reactome annotations into Anvio!
anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioReactomeImportable.tsv

# Remove files that are no longer needed.
rm Reactome?_filtered.txt ReactomeJoined?_?.tsv ReactomeDeleteMe Reactome.temp Reactome_sorted react.txt

echo
echo  "-----------------------------------------------------Reactome ANNOTATIONS WERE IMPORTED INTO ANVIO--------------------------------------------------------------------------------"
echo
echo
echo
echo "					Parse MetaCyc annotations into Anvio"
echo
echo

# Instructions to download the database
# From the website https://biocyc.org/group?id=:ALL-PATHWAYS&orgid=META# choose export to Spreadsheet files...
# From the popup window choose:
# For the MetaCyc annotations: common names
# For the Metacyc IDs: frame IDs

# Upload the two files "Metacyc_IDs.txt" and "Metacyc_annotations.txt" to the server!

# Create local copies of the MetaCyc database, in order to work with them
cp ${MetaCycDBDir}/${MetaCycDBAnnot} .
ln -s ${MetaCycDBDir}/${MetaCycDBIDs}

sed -i -E 's/"//g;s/<[Ii/subpSUB]+>//g;s/&//g;s/;//g' ${MetaCycDBAnnot}


# Paste the two files together and delete the first line which is a header.
paste ${MetaCycDBIDs} ${MetaCycDBAnnot}  | sed '1d' > MetaCyc.anot

# The next block extracts 8 MetaCyc columns from the ${InterproTSVOut} and combines it with the appropriate annotation, from the MetaCyc annotation file we downloaded.
{
for i in {1..8}
do
    awk 'BEGIN{FS="\t"; OFS="\t"} $15 ~ /MetaCyc:/ {print $1,$9,$15}' ${InterproTSVOut} | \
    sed -E 's/KEGG: [0-9]{5}[0-9.+]+\|//g;s/Reactome: \|?[A-Z0-9\-]+\|?//g;s/\|?MetaCyc: /\t/g;s/\|//g;s/\t\t+/\t/g' | \
    awk -v var=$((2+i)) 'BEGIN{FS="\t"; OFS="\t"} {print $var,$1,$2}' > MetaCycDeleteMe

    awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$2;next} {print $2,$1,$3,a[$1]}' MetaCyc.anot MetaCycDeleteMe > MetaCyc.temp

    sed -E '/\t-\t/d' MetaCyc.temp | awk  'BEGIN{FS="\t"; OFS="\t"} $2 !=""' | sort -k3g -t$'\t' > MetaCyc_sorted
    
    awk '$3=="-"' MetaCyc.temp | awk  'BEGIN{FS="\t"; OFS="\t"} $2 !=""' | sort -n >> MetaCyc_sorted
    
     # We filter out only the columns we need.
    awk 'BEGIN{FS="\t";OFS="\t"} !a[$1]++ {print $1,$2,$4}' MetaCyc_sorted > MetaCyc${i}_filtered.txt
    if [[ i -eq 2 ]]; then
        #join MetaCyc$((i-1))_filtered.txt MetaCyc${i}_filtered.txt -t $'\t' > MetaCycJoined1_2.tsv
        awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$2;a[$2]=$3;next} {if ($2 == a[$1]) print $1,$2,$3,$2,$3; else if ($3 == a[$2]) print $1,$2,$3,$2,$3; else print $1,$2,$3,$2"; "a[$1],$3"; "a[$2]}' MetaCyc$((i-1))_filtered.txt MetaCyc${i}_filtered.txt | sed -E 's/; $//g' > MetaCycJoined1_2.tsv
    fi

    if [[ i -gt 2 ]]; then
        #join MetaCycJoined$((i-2))_$((i-1)).tsv MetaCyc${i}_filtered.txt -t $'\t' > MetaCycJoined$((i-1))_${i}.tsv
        awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$2;a[$2]=$3;next} {if ($2 == a[$1]) print $0; else if ($3 == a[$2]) print $0; else print $1,$2,$3,$4"; "a[$1],$5"; "a[$2]}' MetaCyc${i}_filtered.txt MetaCycJoined$((i-2))_$((i-1)).tsv | sed -E 's/; $//g' > MetaCycJoined$((i-1))_${i}.tsv
    fi
done
}

# Next we prepare the Anvio importable tsv file with the correct headers and the correct columns selected by awk.
awk 'BEGIN{FS="\t"; OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} {print $1,"MetaCyc",$4,$5,"0"}' MetaCycJoined7_8.tsv | sed -E 's/ \t;//g;s/\t;/;/g;s/(; ;)+//g;s/ ; /;/g;s/; $//g' > AnvioMetaCycImportable.tsv


# Finally, we import the MetaCyc annotations into Anvio.
anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioMetaCycImportable.tsv

rm MetaCyc?_filtered.txt MetaCycJoined?_?.tsv MetaCycDeleteMe MetaCyc.temp MetaCyc_sorted MetaCyc.anot MetaCyc_IDs.txt MetaCyc_annotations.txt 

echo
echo
echo  "-----------------------------------------------------MetaCyc ANNOTATIONS WERE IMPORTED INTO ANVIO--------------------------------------------------------------------------------"
##echo
##echo
##echo "					Parse Kegg Interpo Pathways into Anvio"
##echo

# NOTE:
# KEGG PATHWAYS ARE NOT PROVIDED ANYMORE BY INTERPROSCAN!!!


### The downloaded files can be found in these web sites:
### https://www.genome.jp/kegg-bin/get_htext?br08901.keg  for the Pathway Maps
### https://www.genome.jp/kegg-bin/get_htext?ko01000.keg for the Enzymes


### First we start with the Kegg Pathway annotation file. First we flatten it, and then reformat it to our needs.
##{
##kegfile="br08901.keg"
##while read -r prefix content
##do
##    case "$prefix" in A) col1="$content";; 
##                      B) col2="$content" ;; 
##                      C) echo -e "$col1\t$col2\t$content";; 
##    esac 
##done < <(sed '/^[#!+]/d;s/<[^>]*>//g;s/^./& /' < "$kegfile") | sed -E 's/  /\t/g' | awk 'BEGIN{FS="\t";OFS="\t"} {print $3,$1", "$2", "$4}' > KO_Pathway_br08901.txt
##}


### Next we reformat the Kegg Enzymes annotations file. 
##{
##kegfile="ko01000.keg"
##while read -r prefix content
##do
##    case "$prefix" in A) col1="$content";; 
##                      B) col2="$content" ;; 
##                      C) col3="$content" ;; 
##                      D) echo -e "$col1\t$col2\t$col3\t$content";; 
##    esac 
##done < <(sed '/^[#!+]/d;s/<[^>]*>//g;s/^./& /' < "$kegfile") | sed -E 's/^[0-9]\. //g;s/([0-9.]+)  /\1\t/g' | awk 'BEGIN{FS="\t";OFS="\t"} {print $6,$1", "$3", "$5", "$7}'  > KO_Enzymes_ko01000.txt
##}

### Next block of code will extract all useful information from the two files (Kegg Enzymes and Pathways) and join the fields that have more than one annotation (from Interpro) into the same field. Notice that because of the array filtering we apply at the first step (i.e. awk 'NF==4' ) we cannot use the join command as we did in previous implementations, because we have varying number of gene caller ids in each iteration of the for loop. Instead we join with awk with the help of arrays (for the nth time in this script!).
##{
##for i in {1..4}
##do
##    awk 'BEGIN{FS="\t"; OFS="\t"} $15 ~ /KEGG:/ {print $1,$9,$15}' ${InterproTSVOut} | \
##    sed -E 's/\|?MetaCyc: [A-Z0-9\-]+\|?//g;s/Reactome: \|?[A-Z0-9\-]+\|?//g;s/\|KEGG:/\t/g;s/KEGG://g;s/\|//g' | \
##    awk -v var=$((2+i)) 'BEGIN{FS="\t"; OFS="\t"} {print $var,$1,$2}' | \
##    sed -E 's/\+/\t/g;s/ //g' | awk 'NF==4' > KeggPathwaysInterproDeleteMe

##    awk 'BEGIN{FS="\t"; OFS="\t"} NR==FNR{a[$1]=$2;next} {print $2,$3,$1,$4,a[$1]}' KO_Pathway_br08901.txt KeggPathwaysInterproDeleteMe > KeggPathwaysInterpro.temp
##    

##    awk 'BEGIN{FS="\t"; OFS="\t"} NR==FNR{a[$1]=$2;next} {print $2,$3"+"$1,$4,"PATHWAY: "$5" ENZYME: "a[$1]}' KO_Enzymes_ko01000.txt KeggPathwaysInterpro.temp > KeggPathwaysEnzymesInterpro.temp

##    sed -E '/\t-\t/d' KeggPathwaysEnzymesInterpro.temp | sort -k3g -t$'\t' > KeggPathwaysEnzymesInterpro_sorted
##    #sed -E '/\t-\t/d' MetaCyc.temp | awk  'BEGIN{FS="\t"; OFS="\t"} $2 !=""' | sort -k3g -t$'\t' > MetaCyc_sorted
##    
##    awk '$3=="-"' KeggPathwaysEnzymesInterpro.temp |sort -n >> KeggPathwaysEnzymesInterpro_sorted

##    awk 'BEGIN{FS="\t"; OFS="\t"} !a[$1]++ {print $1,$2,$4}' KeggPathwaysEnzymesInterpro_sorted  > KeggPathwaysEnzymes${i}_filtered.txt

##    if [[ i -eq 2 ]]; then
##        awk 'BEGIN{FS="\t"; OFS="\t"} NR==FNR{a[$1]=$2;b[$1]=$3;next} {print $1,$2"; "a[$1],$3"; "b[$1]}' KeggPathwaysEnzymes2_filtered.txt KeggPathwaysEnzymes1_filtered.txt > KeggPathwaysEnzymesJoined1_2.tsv
##    fi
##    echo
##    if [[ i -gt 2 ]]; then
##        awk 'BEGIN{FS="\t"; OFS="\t"} NR==FNR{a[$1]=$2;b[$1]=$3;next} {print $1,$2"; "a[$1],$3"; "b[$1]}'  KeggPathwaysEnzymes${i}_filtered.txt KeggPathwaysEnzymesJoined$((i-2))_$((i-1)).tsv > KeggPathwaysEnzymesJoined$((i-1))_${i}.tsv
##    fi
##done
##}

### Next, we remove some of the remnant semicolons, to make the results more readable.
##sed -i 's/; ;//g;s/ ;/;/g' KeggPathwaysEnzymesJoined3_4.tsv

### Next, we create the tab delimited file to be imported into Anvio, with the correct headers etc.
##awk 'BEGIN{FS="\t"; OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} {print $1,"Kegg Pathways + Enzymes",$2,$3,"0"}' KeggPathwaysEnzymesJoined3_4.tsv  > AnvioKeggPathwaysEnzymesImportable.tsv

#### Finally, import the Annotations into Anvio!
###anvi-import-functions --contigs-db ${StrX}_Pilon_contigs.db  --input-files AnvioKeggPathwaysEnzymesImportable.tsv

#### Remove any files not needed anymore.
###rm KeggPathwaysEnzymes?_filtered.txt KeggPathwaysEnzymesJoined?_?.tsv KeggPathwaysInterproDeleteMe KeggPathwaysInterpro.temp KeggPathwaysEnzymesInterpro.temp KeggPathwaysEnzymesInterpro_sorted 

##echo
##echo
##echo "-----------------------------------KEGG INTERPRO ANNOTATIONS WERE PARSED INTO ANVIO---------------------------------------------------"


<< NOTE

https://github.com/ebi-pf-team/interproscan/issues/202

No EC codes in KEGG pathway identifier #202


a-ag commented on Apr 20 • 
Hi,

I am running interproscan to get pathway data on protein sequences, but I don't see any EC codes attached with KEGG pathway identifier. I remember that Interproscan provided EC codes in an Interproscan4 version. Previously, in the interproscan json output, the pathwayXRefs for KEGG mentioned the id as: 00240+2.7.7.7 and now it only has 00240.

Is this something that interproscan does not provide anymore or is there an argument/flag that needs to be passed for EC codes to be available?

Cheers!

 
Contributor
gsn7 commented on Apr 29
You are right the pathwayXRefs for KEGG mentioned the id such as: 00240+2.7.7.7. The change to 00240 was a bug. Unfortunately, we have discontinued providing KEGG pathways due to licencing restrictions. So, latest releases, which could have fixed this problem don't have KEGG pathways.
	
NOTE

echo
echo
echo "======================================================================"
echo "SLRUM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
elapsed=$SECONDS
echo Time taken: `printf '%dd %dh:%dm:%ds\n' $((elapsed/86400)) $((elapsed%86400/3600)) $((elapsed%3600/60)) $((elapsed%60))`

exit 0
