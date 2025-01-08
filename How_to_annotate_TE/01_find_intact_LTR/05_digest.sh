#!/bin/bash 
#Program: digest
#SBATCH -J digest
#SBATCH -n 1
#SBATCH --mem=20GB
#SBATCH --array=1-1%1

INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list)
Name=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $1}')
Genome=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $2}')
Fasta=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $2}'| awk 'BEGIN{FS="/"}{print $8}')
GenomeID=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $3}')

mkdir ${Name}
cd ${Name}
for PasslistGFF3 in $(ls -1 /nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/03_retriever/${Name}/${Fasta}*.list.gff3); do
	if [ -e "/nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/03_retriever/${Name}/${Fasta}.mod.pass.list.gff3" ]; then
		
		## header1
		less $PasslistGFF3 | grep '^##sequence-region' > header1.txt 
		
		## header2
		less  $PasslistGFF3  | grep "^${GenomeID}" | cut -f 1 | less | sort | uniq | awk '{print "#"$1}' > header2.txt
		
		## header3
		echo "##gff-version 3" > header3.txt
		
		##content
		less $PasslistGFF3  |\
		grep -v "^##sequence-region" |\
		grep -v "##gff-version 3" |\
		grep -v "##date" |\
		grep -v "##ltr_identity: Sequence identity (0-1) between the left and right LTR region." |\
		grep -v "##tsd: target site duplication." |\
		grep -v "##seqid source sequence_ontology start end score strand phase attributes" |\
		awk 'BEGIN{FS=";"; OFS="\t"}{split($1,a,"\t") ; \
		if ($1 ~/ID=repeat_region_/ ) print a[1],a[2],a[3],a[4],a[5],a[6],"?",a[8],a[9]; \
		else if (a[3] ~/LTR_retrotransposon/) print a[1],a[2]"_"a[3],"LTR_retrotransposon",a[4],a[5],a[6],"?",a[8],$2";"$3; \
		else if ($1 == "###") print "###"; \
		else print a[1],a[2],a[3],a[4],a[5],a[6],"?",a[8],$2}' \
		> content.txt
		
		## conbine
		cat header3.txt header1.txt header2.txt content.txt > combine.txt
		
		## print
		awk '{OFS="\t"}{$1=$1}1' combine.txt > ${Name}.harvestformat.gff3
		
		
	elif [ -e "/nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/03_retriever/${Name}/${Fasta}.pass.list.gff3" ]; then
		
		## header1
		less $PasslistGFF3 | grep '^##sequence-region' > header1.txt 
		
		## header2
		less  $PasslistGFF3  | grep "^${GenomeID}" | cut -f 1 | less | sort | uniq | awk '{print "#"$1}' > header2.txt
		
		## header3
		echo "##gff-version 3" > header3.txt
		
		##content
		less $PasslistGFF3  |\
		grep -v "^##sequence-region" |\
		grep -v "##gff-version 3" |\
		grep -v "##date" |\
		grep -v "##ltr_identity: Sequence identity (0-1) between the left and right LTR region." |\
		grep -v "##tsd: target site duplication." |\
		grep -v "##seqid source sequence_ontology start end score strand phase attributes" |\
		awk 'BEGIN{FS=";"; OFS="\t"}{split($1,a,"\t") ; \
		if ($1 ~/ID=repeat_region_/ ) print a[1],a[2],a[3],a[4],a[5],a[6],"?",a[8],a[9]; \
		else if (a[3] ~/LTR_retrotransposon/) print a[1],a[2]"_"a[3],"LTR_retrotransposon",a[4],a[5],a[6],"?",a[8],$2";"$3; \
		else if ($1 == "###") print "###"; \
		else print a[1],a[2],a[3],a[4],a[5],a[6],"?",a[8],$2}' \
		> content.txt
		
		## conbine
		cat header3.txt header1.txt header2.txt content.txt > combine.txt
		
		## print
		awk '{OFS="\t"}{$1=$1}1' combine.txt > ${Name}.harvestformat.gff3
		
	else
		echo "ERROR: no pass.list.gff3 for intact Gypsy/Copia"
	fi
	
	rm header1.txt
	rm header2.txt
	rm header3.txt
	rm content.txt
	rm combine.txt
	
done


#2 build index

module load genometools/1.6.1
module load  hmmer/3.3

gt suffixerator \
        -db ${Genome} \
        -indexname ${Name}.gtindex \
        -tis -suf -lcp -des -ssp -sds -dna


#3 sort input gff3
## input.gff3: gff3 format returned by LTRharvest.
## input: C1892.scaffold.fa.pass.list.harvestformat.gff3

Gff3=$(echo $Name".harvestformat.gff3")

gt gff3 -sort ${Gff3} > ${Gff3}.sorted.gff3

#4 run digest
gt ltrdigest -pptlen 10 30 \
         -outfileprefix ltr_WG \
         -hmms /nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/05_digest/TE.pfam.retriever_Biodirect.hmm \
         -seqfile ${Genome} \
         -aaout \
         -matchdescstart ${Gff3}.sorted.gff3 \
         > ${Name}.ltrdigest.gff3

echo "DONE"
