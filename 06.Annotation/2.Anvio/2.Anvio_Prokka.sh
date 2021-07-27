#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="1st_Anvio_Prokka"
#SBATCH --output=1st_Anvio_Prokka_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END


#set -euo pipefail # Check https://dev.to/pencillr/newbie-at-bash-scripting-heres-some-advice-j1a

generalInfo () {
    cat <<END

	This script takes as input one argument. 
	For Strain 01 that would be: 
	sbatch ${ScriptName} Strain01
	This script runs from the master folder of Anvio and then we change directory to the folder of each specific Strain folder.
	We assume that in the master folder there are already folders named Strain01 Strain02 etc.

	This is the Anvio pipeline. PART ONE
	Check the tutorial http://merenlab.org/2016/06/18/importing-functions/ "Importing functions into contigs database"

	For this script to run we should already have finished the Prokka annotation, because we use the output of the Prokka pipeline as input to the contigs database.
	We should also have done the genecalling with GeneMarkS2, online at http://exon.gatech.edu/GeneMark/genemarks2.cgi.

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

echo
# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
Assembler="Spades"  # Choose either "Pilon" or "Spades"
OutputDir=${HOME}/Titlos_ktisis/Anvio/${StrX}_${Assembler}

CntgDB=${StrX}_${Assembler}_contigs.db

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

ProkkaDir=${HOME}/Titlos_ktisis/Prokka/${StrX}_Prokka_${Assembler}
GeneMarkS2Dir=${HOME}/Titlos_ktisis/GeneMarkS2/${StrX}_${Assembler}
#GfftoGenbankDir=${HOME}/Software/Gff_to_genbank
GfftoGenbankDir=${HOME}/Software/bcbb/gff/Scripts/gff

ProkkaGff=${StrX}_Prokka_${Assembler}.gff
GeneMarkGff=${StrX}_${Assembler}_gms2.gff

# Variables for Anvio command options that can be altered
ContigDBSplitLength=20000 # This is mostly for better visualization. Not to be changed.
ContigDBKmerSize=4 # Tetranucleotides is the default and the most used one. Not to be changed.

# We substitute the command anvi-get-sequences-for-gene-calls with a parameter, but this is certainly not necessary!
anviGetSeq='anvi-get-sequences-for-gene-calls'

echo


. "${HOME}/Software/miniconda3/etc/profile.d/conda.sh"
conda activate anvio-7
export LC_ALL=en_US.UTF-8
export PATH="${HOME}/Software/Diamond:$PATH"
export PATH="${HOME}/Software/gms2_linux_64:$PATH"



#--------------------------SANITY CHECKS -------------------------------------------------------------------#
# Check that the files we have access to the files needed to run this pipeline
if [[ ! -s $ContigDir/$Contig || ! -s $ProkkaDir/${ProkkaGff} || ! -s ${GeneMarkS2Dir}/${GeneMarkGff} ]]; then
    echo "Some or all of the files necessary to run this pipeline (Contigs, Prokka annotation, GeneMarkS2 annotation) are missing!" >&2
    echo >&2
    generalInfo >&2
    exit 1
fi

if [[ ! -e ${GfftoGenbankDir}/gff_to_genbank.py ]]; then
    echo "Software to get genbank sequences is missing in ${GfftoGenbankDir}/gff_to_genbank.py" >&2
    echo >&2
    generalInfo >&2
    exit 1
fi

echo
echo "Anvio version:"
echo
anvi-self-test --version
echo
echo
echo "			Anvio Analysis begins!!"
echo

# Create output directory.
[[ ! -d ${OutputDir} ]] && mkdir -p ${OutputDir}
cd ${OutputDir}

echo
# Create a soft link of the contigs fasta file in the working directory
ln -s "${ContigDir}"/"${Contig}"
echo
echo
echo
echo '		Export Prokka gene calls to an importable format for Anvio '
echo

# Get the gene calls and annotations from a gff3 Prokka output in a format that anvio can import.

# We have to subtract one (1) from the start position column. Please check here http://merenlab.org/2016/06/22/anvio-tutorial-v2/#external-gene-calls

# First, export CDS gene calls found by prodigal, and save it to a temporary file. We will add the gene caller ids at the end, after we concatenate all the temporary files with the gene calls into one file.
awk 'BEGIN{OFS="\t"} $3 ~ /CDS/{print $1,($4-1),$5,$7,$8,"1","Prodigal","2.6.3"}' "${ProkkaDir}"/"${ProkkaGff}" > "${StrX}"_Prokka_Anvio_temporary_gene_calls.tsv

# Second, we export Aragorn gene calls and change call_type to "2" (to prevent Anvio from exporting these genes as amino acids!) and append it to the same temporary file we save the CDS annotation.
awk 'BEGIN{OFS="\t"}  $2 ~ /Aragorn:001002/{print $1,($4-1),$5,$7,"0","2","Aragorn","1.2"}' "${ProkkaDir}"/"${ProkkaGff}" >> "${StrX}"_Prokka_Anvio_temporary_gene_calls.tsv

# Third, we export barnap gene calls and change call_type to "2" (to prevent Anvio from exporting these genes as amino acids!) and append it to the same temporary file we save the CDS and Aragorn annotation.
awk 'BEGIN{OFS="\t"} $2 ~ /barrnap:0.9/{print $1,($4-1),$5,$7,"0","2","Barnap","0.9"}' "${ProkkaDir}"/"${ProkkaGff}" >> "${StrX}"_Prokka_Anvio_temporary_gene_calls.tsv

# Fourth, we export Infernal gene calls and change call_type to "2" (to prevent Anvio from exporting these genes as amino acids!) and append it to the same temporary file we save the CDS, Aragorn and barnap annotation.
awk 'BEGIN{OFS="\t"} $2 ~ /Infernal:001001/{print $1,($4-1),$5,$7,"0","2","Infernal","1.1"}' "${ProkkaDir}"/"${ProkkaGff}" >> "${StrX}"_Prokka_Anvio_temporary_gene_calls.tsv

# NOTE
# We do NOT export SignalP from the Prokka gff file, as we will search and import these later on

#------------------------------------------------------------------------------------------------------------------------------------

# Next we will try to export gene calls from GenemarkS2, which we have already run online http://exon.gatech.edu/GeneMark/genemarks2.cgi
# We will export only gene calls that have not been found by Prokka.
# We MUST append the GeneMarkS2 gene calls last, after all the other gene calls have already been added to our temporary file!

# Prepare the two files (one from Prokka, one from GenemarkS2) for comparison
awk '$3 ~ /CDS/  {print $1,$4,$5,$7,$8}' "${ProkkaDir}"/"${ProkkaGff}" > "${StrX}"_Prokka_gff.out
awk '$3 ~ /CDS/  {print $1,$4,$5,$7,$8}' "${GeneMarkS2Dir}"/"${GeneMarkGff}" > "${StrX}"_GeneMark_gff.out

# Check if GeneMarkS2 gene calls can be found in Prokka gene calls. If not, then write these gene calls to a new file.
{
while read -r p; do
  if ! grep -q "$p" "${StrX}"_Prokka_gff.out; then
    echo "$p" >> GeneMarkS2_genes_NOT_in_Prokka.tsv
  fi
done < "${StrX}"_GeneMark_gff.out
}

# Prepare the temporary file (i.e. every column exept for the gene call ids). Please don't forget that we have to subtract 1 (one) from the start position!
awk 'OFS="\t" {print $1,($2-1),$3,$4,$5,"1","GeneMarkS2","3"}' GeneMarkS2_genes_NOT_in_Prokka.tsv >> "${StrX}"_Prokka_Anvio_temporary_gene_calls.tsv

#------------------------------------------------------------------------------------------------------------------------------------

# Create the necessary headers for the output file
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' 'gene_callers_id' 'contig' 'start' 'stop' 'direction' 'partial' 'call_type' 'source' 'version' > "${StrX}"_Prokka_Anvio_gene_calls.tsv

# Finally, we append all the gene calls from the temporary file we created earlier to the final file we are going to use to parse all the Prokka annotations into Anvio. 
# IMPORTANT
# We have to subsitute all the +,- with f,r and print it with the first column being the number of the row (which corresponds to the gene caller id!)
awk 'OFS="\t" {sub("+","f",$4);sub("-","r",$4)} {print NR,$0}' "${StrX}"_Prokka_Anvio_temporary_gene_calls.tsv >> "${StrX}"_Prokka_Anvio_gene_calls.tsv

echo "----------------------PROKKA GENE CALLS WERE EXPORTED!!!!---------------------------------------------------------"
echo
echo
echo '			Create Contigs DATABASE'
echo
# Create the Contig Database description file
echo "${StrX} contigs database created from the ${Assembler} output scaffold fasta files, after filtering out contigs marked either as contamination or plasmids by blast search. The exact pipeline that created these scaffolds is described elswhere. Gene finding was done inside Prokka v1.14.6 with Prodigal 2.6.3. Annotation was Completed with Prokka, using the genus Rossellomorea command in Prokka (after creating the appropriate database). We have also included gene calls from GeneMarkS2, but only those that have not already been found by Prokka (prodigal)" >  "${StrX}"_ContigDB_Anvio_decription.txt

# Finally, create the Contigs Database with the Prokka gene calls.
anvi-gen-contigs-database -f "${Contig}" \
			  --project-name "${StrX}"_${Assembler} \
			  --output-db-path "${CntgDB}" \
                          --description "${StrX}"_ContigDB_Anvio_decription.txt \
			  --split-length "${ContigDBSplitLength}"\
                          --kmer-size "${ContigDBKmerSize}" \
			  --external-gene-calls "${StrX}"_Prokka_Anvio_gene_calls.tsv \
			  --ignore-internal-stop-codons # We ignore stop codons in order to avoid to be forced to a stop by Anvio. But we should not forget to check the stdout file!
echo
echo
echo '---------------------------------CONTIGS DATABASE CREATED!!!---------------------------------------'
echo
echo
echo '		Parse Prokka Functional Annotations into Anvio'
echo

'true' << ////
	Now, we will export functional annotations from the ${ProkkaGff} file output of the Prokka pipeline to an importable tsv file
	
	IMPORTANT no1
	Gene caller ids should be the same in this file as the ones assigned by Anvio to the contigs database. 
	For this reason we use a similar awk command (i.e. the same .gff file, filtered by CDS/Aragorn/Barnap/Infernal) as the one used 
	to parse the gene calls into Anvio in the previous steps.
	
	IMPORTANT no2
	We MUST follow the same order of parsing these annotations as was the order that we parsed the gene calls.
	i.e. First MUST be parsed CDS, second MUST be parsed Aragorn, third MUST be parsed Barnap and last MUST be parsed Infernal!!
////

# Fisrt, export  CDS annotations to an importable file format, and save these annotations to a temporary file
awk -F'[\t=;,]' '$3 ~ /CDS/ {print $0}' "${ProkkaDir}"/"${ProkkaGff}" | \
# We also remove non useful info (with sed) to have finally only the gene accession and function information
sed -E 's/ID=.+sequence://g; s/ID=.+motif://g; s/;locus_tag=.+;product=/\t/g' | \
awk -F"\t" '{print "Prokka""\t"$9"\t"$10}' > AnvioProkkaTEMPFunctionalAnnotationsImportable.tsv

# At this point, we want to have the number of CDS annotations produced by Prokka Prodigal, and save it to a variable for later use.
ProkkaCDSNum=$(wc -l AnvioProkkaTEMPFunctionalAnnotationsImportable.tsv | cut -d' ' -f1)
echo
echo "Number of CDS annotations produced ONLY by Prokka is ${ProkkaCDSNum}"
echo
echo

# Second, we follow the same procedure to export the Aragorn annotation, and append to the same temporary file we saved the CDS annotations.
awk -F'[\t=;]' '$2 ~ /Aragorn:001002/ {print $0}' "${ProkkaDir}"/"${ProkkaGff}" | \
sed -E 's/ID=.+locus_tag=//g; s/;product=/\t/g' | \
awk -F"\t" '{print "Prokka:Aragorn""\t"$9"\t"$10}' >> AnvioProkkaTEMPFunctionalAnnotationsImportable.tsv

# Third, we repeat the same procedure to export the Barnap annotation, and append to the same temporary file we saved the CDS and Aragorn annotations.
awk -F'[\t=;]' '$2 ~ /barrnap:0.9/ {print $0}' "${ProkkaDir}"/"${ProkkaGff}" | \
sed -E 's/ID=.+locus_tag=//g; s/;product=/\t/g' | \
awk -F"\t" '{print "Prokka:Barnap""\t"$9"\t"$10}' >> AnvioProkkaTEMPFunctionalAnnotationsImportable.tsv

# Fourth, we repeat the same procedure to export the Infernal annotation, and append to the same temporary file we saved the CDS, Aragorn and Barnap annotations.
awk -F'[\t=;]' '$2 ~ /Infernal:001001/ {print $0}' "${ProkkaDir}"/"${ProkkaGff}" | \
sed -E 's/ID=.+locus_tag=//g; s/;product=/\t/g' | \
awk -F"\t" '{print "Prokka:Infernal""\t"$9"\t"$10}' >> AnvioProkkaTEMPFunctionalAnnotationsImportable.tsv

#------------------------------------------------------------------------------------------------------------------------------------------

# Create the file to be parsed into Anvio with the correct headers, in order to append later the Prokka functional annotation
printf '%s\t%s\t%s\t%s\t%s\n' 'gene_callers_id' 'source' 'accession' 'function' 'e_value' > AnvioProkkaFunctionalAnnotationsImportable.tsv

# Finally parse the temporary file we created to Anvio, with the correct gene calls, and (after that!!) apply any further filtration, like removing hypothetical proteins etc. Finally we fix the e-value to 0 (zero).
awk '{print NR"\t"$0}' AnvioProkkaTEMPFunctionalAnnotationsImportable.tsv | \
awk -F"\t" '$4!="hypothetical protein" {print $1"\t"$2"\t"$3"\t"$4"\t""0"}' >> AnvioProkkaFunctionalAnnotationsImportable.tsv


# Finally, parse functional annotations into Anvio!
anvi-import-functions --contigs-db "${CntgDB}" --input-files AnvioProkkaFunctionalAnnotationsImportable.tsv

echo
echo '-------------------------PROKKA FUNCTIONAL ANNOTATIONS WERE PARSED INTO ANVIO!-------------------------------'
echo
echo
echo '				Export Amino Acid fasta file with NO empty records. Also export genes to  a gff file'
echo
echo

# We are going to use the amino acid fasta file from the Anvio Prokka pipeline as input for all subsequent analyses
# For this reason, we want to have an amino acid sequence fasta file without any empty fasta records. These empty records could potentially cause problems in the analysis downstream.

# Export amino acids from the the contigs database into a fasta file
$anviGetSeq --contigs-db "${CntgDB}" --get-aa-sequences --output-file "${StrX}"_amino_acid_seqs_with_empty_fields.fasta

# Remove empty fasta records with awk.
# With RS, we declare that the records are separated by >. If the record contains a second line $2 (i.e. is not an empty fasta record), print the record (add a > in front).
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' "${StrX}"_amino_acid_seqs_with_empty_fields.fasta > "${StrX}"_amino_acid_seqs.faa
echo

# Export gff3 from anvi contig database
${anviGetSeq} -c "${StrX}"_${Assembler}_contigs.db -o "${StrX}"_Anvio.gff --export-gff3

echo
echo
echo "--------------------------------------AMINO ACID FASTA FILE READY-------------------------------------------------------------"
echo
echo "		Create new fasta file WITHOUT the GenemarkS2 gene calls, ONLY PROKKA CDS gene calls!. Also create a gff file with ONLY PROKKA CDS annotations"
echo

[[ -e patterns.txt ]] && rm patterns.txt
# Next, we create the file "patterns.txt" which will be used as input in grep. Patterns.txt contains one gene caller id per line, for all genes found by Prokka.
for i in $(seq 1 "$ProkkaCDSNum"); do echo "^>$i\b" >> patterns.txt; done

# Finally, we linearize the FASTA with awk, pipe to grep to filter for items of interest named in patterns.txt, then pipe to tr to delinearize:
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "${StrX}"_amino_acid_seqs.faa | \
grep -Ef patterns.txt - | \
tr "\t" "\n" >  "${StrX}"_ONLY_PROKKA_amino_acid_seqs.faa
echo

# Also, create a gff file with ONLY Prokka CDS gene calls.
awk -v var="$ProkkaCDSNum" -F"[\t=]" '(NR==1) || ($10 <= var )' "${StrX}"_Anvio.gff > "${StrX}"_ONLY_PROKKA_amino_acid_seqs.gff
# IMPORTANT!
# Only this gff file can actually be used because otherwise we get all the RNAs and non coding RNAs as CDS, which is a huge mistake!!! 

echo
echo
echo '--------------------------AMINO ACID FASTA SEQUENCES and GFF file ONLY FROM PROKKA CDS ANNOTATIONS ARE READY!----------------------------------'
echo
echo 
echo "				Create Genbank File only from Prokka gene calls"
echo

# We want to use this file for uploading to the Rast servers, and later parse the Rast annotations into Anvio, but besides this it is not a bad idea to have a genbank file.

if [[ ${Assembler}  == "Pilon" ]]
then
    sed 's/Scaffold_//g; s/_length//g; s/_pilon//g' "${Contig}" > "${StrX}"_contigs_for_gff_to_genbank.fasta
    sed 's/Scaffold_//g; s/_length//g; s/_pilon//g' "${StrX}"_ONLY_PROKKA_amino_acid_seqs.gff > "${StrX}"_ONLY_PROKKA_aa_for_gff_to_genbank.gff
elif [[ ${Assembler}  == "Spades" ]]
then
    sed -E 's/NODE_//g; s/_length//g; s/_cov_[0-9]+\.[0-9]+//g' "${Contig}" > "${StrX}"_contigs_for_gff_to_genbank.fasta
    sed -E 's/NODE_//g; s/_length//g; s/_cov_[0-9]+\.[0-9]+//g' "${StrX}"_ONLY_PROKKA_amino_acid_seqs.gff > "${StrX}"_ONLY_PROKKA_aa_for_gff_to_genbank.gff
fi

# NOTE:
# gff_to_genbank.py does not run with newer versions of Biopython
# We have to downgrade downgrade Biopython to v1.77 by creating a new python3.6 environment and then running pip3 install biopython==1.77
conda deactivate
conda activate biopythonOLD


"${GfftoGenbankDir}"/gff_to_genbank.py "${StrX}"_ONLY_PROKKA_aa_for_gff_to_genbank.gff "${StrX}"_contigs_for_gff_to_genbank.fasta

mv "${StrX}"_ONLY_PROKKA_aa_for_gff_to_genbank.gb "${StrX}"_ONLY_PROKKA_aa.gbk

rm "${StrX}"_contigs_for_gff_to_genbank.fasta "${StrX}"_ONLY_PROKKA_aa_for_gff_to_genbank.gff


# Now we have to go to the Rast server http://rast.nmpdr.org/ and upload the ${StrX}_ONLY_PROKKA_aa.gbk. When we have the results we have to upload them to the working directory in order to import these annotations into Anvio!



echo
echo
echo '		Export amino acid sequences from Anvio for KEGG functional annotation'
echo
# Instructions can be found here http://merenlab.org/2018/01/17/importing-ghostkoala-annotations/
echo


# BlastKOALA does have a gene number limit of 5000 genes. 
# If your anvi’o database has more than that number of genes, you can split the file into multiple FASTA files, 
# and concatenate the BlastKOALA outputs you will be generating during this step.

TotalProteins=$(grep -c ">" "${StrX}"_amino_acid_seqs.faa)

if (( TotalProteins > 5000 ))
then 
    if (( TotalProteins/2 < 5000 ))
    then
         # Check that the output does not already exist. If it does delete the files.
         [[ -e FirstHalf.txt ]] && rm FirstHalf.txt
         [[ -e SecondHalf.txt ]] && rm SecondHalf.txt
         
         Half=$(( TotalProteins/2))
         for i in $(seq 1 "$Half"); do echo "^>$i\b" >> FirstHalf.txt; done
         for i in $(seq "$(( Half+1 ))" "$TotalProteins"); do echo "^>$i\b" >> SecondHalf.txt; done

        # We linearize the FASTA with awk, pipe to grep to filter for items of interest named in FirstHalf.txt, then pipe to tr to delinearize:
        awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "${StrX}"_amino_acid_seqs.faa | \
        grep -Ef FirstHalf.txt - | \
        tr "\t" "\n" >  "${StrX}"_FirstHalf_amino_acid_seqs.faa
        
        # We linearize the FASTA with awk, pipe to grep to filter for items of interest named in SecondHalf.txt, then pipe to tr to delinearize:
        awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "${StrX}"_amino_acid_seqs.faa | \
        grep -Ef SecondHalf.txt - | \
        tr "\t" "\n" >  "${StrX}"_SecondHalf_amino_acid_seqs.faa
        
        # The fasta sequence headers start with a number (e.g. 0 or 1 etc) and the next command adds 'genecall' as a prefix of each header
        sed 's/>/>genecall_/g' "${StrX}"_FirstHalf_amino_acid_seqs.faa > "${StrX}"_FirstHalf_Koala_${Assembler}.fasta
        
        sed 's/>/>genecall_/g' "${StrX}"_SecondHalf_amino_acid_seqs.faa > "${StrX}"_SecondHalf_Koala_${Assembler}.fasta
        
        rm FirstHalf.txt SecondHalf.txt "${StrX}"_FirstHalf_amino_acid_seqs.faa "${StrX}"_SecondHalf_amino_acid_seqs.faa
    else
        echo "Check file ${StrX}_amino_acid_seqs.faa it seems to have more than 10,000 proteins! Something is probably wrong here!"  >&2
    fi
else 
    # The fasta sequence headers start with a number (e.g. 0 or 1 etc) and the next command adds 'genecall' as a prefix of each header
    sed 's/>/>genecall_/g' "${StrX}"_amino_acid_seqs.faa > "${StrX}"_Koala_${Assembler}.fasta
fi




echo
echo '--------------------------------BlastKOALA/KEEG AA FASTA FILE READY FOR UPLOAD!-----------------------------------------------------'


## Remove files not needed downstream in our analysis.
rm patterns.txt AnvioProkkaTEMPFunctionalAnnotationsImportable.tsv "${StrX}"_amino_acid_seqs_with_empty_fields.fasta

# ATTENTION, FURTHER STEPS NEEDED!
# "${StrX}"_Koala_${Assembler}.fasta (Or "${StrX}"_FirstHalf_Koala_${Assembler}.fasta AND "${StrX}"_SecondHalf_Koala_${Assembler}.fasta) file is now ready to be submitted to BlastKOALA for annotation.
# To do this go to the BlastKOALA webserver https://www.kegg.jp/blastkoala/, 
# and click the “Choose File” button underneath the section that says Upload query amino acid sequences in FASTA format. 
# From the menu you will upload the ${StrX}_KEGG_protein-sequences.fasta, and will be asked to provide an email address.

# Also don't forget to run the Rast annotation online at http://rast.nmpdr.org/, selecting not Rasttk, but Rast while keeping the gene calls from the genbank file!


	
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
