#!/bin/bash
#Program: blast
#SBATCH -J 20H
#SBATCH -n1 -c 20
#SBATCH --mem=100GB
#SBATCH --array=1-6%6

module purge

# bach
#module load blast/2.10.1
#module load bedtools/2.29.2

# holst
module load blast/2.13.0
module load bedtools2/2.30.0


INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p Helitron_genome20.missing20230331.list)
#INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p Helitron_genome20.list)
#INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p TIR_Agenome.list2)
#INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p LTR_Agenome.list)
#INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p TIR_Agenome.list)
#INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p TE8_Cgenome_LTR_noII.list)
#INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p TIR_Cgenome_missing.list)
TE=$(echo $INPUT | awk 'BEGIN{FS=";"}{print $1}')
TE_inList=$(echo $INPUT | awk 'BEGIN{FS=";"}{print $3}')
Genome=$(echo $INPUT | awk 'BEGIN{FS=";"}{print $2}')

################################################

#-------------- preparing input

# 1. softlink subject sequence
#mkdir 00_strucFasta
#ln -s OLD_FASTA $Fasta
#00_strucFasta/strucTE.Helitron.fa -> /nfs/project1/Repeat/10_integrate/06_clusterTE/03_mcl/32genomes_Helitron_filter/03_mcl/ACgenome12.struc_Helitron.fa


# 2. get query homo bed
#TE="Helitrom"
#Genome="II"

GenomeFasta=$(less ~/name32.list | grep -w $Genome | cut -f1)
Name=$(less ~/name32.list | grep -w $Genome | cut -f7 -d "/")
Gff=$(echo "/nfs/project1/Repeat/10_integrate/01_CACTA/"$Name"/"$Name".EDTA_scanCACTA.gff3")
Bed=$(echo "/nfs/project1/Repeat/10_integrate/07_homo80/32genomes_Helitron_filter/"$Name"/"$Name".EDTA_scanCACTA."$TE".bed")
Class=$(less /nfs/project1/Repeat/10_integrate/07_homo80/32genomes_Helitron_filter/EDTAgff_class.list | grep -w $TE | cut -f2)

mkdir $Name
cd $Name

# 3. ############################################################################

# 3. filter TE and add ID name - LTRCO, LTRGY and LTRRT
#less ${Gff} | awk -v TE=${TE} -v Class=${Class} 'BEGIN{FS="\t";OFS="\t"}\
#{if ($3 == Class ) {split($9,x,"ID=");split(x[2],y,";");split(y[1],z,"_");\
#if(z[2]!="homo"){ID=TE;NUM=z[2]}\
#else if(z[2]=="homo"){ID="TE_homo";NUM=z[3]};\
#split($9,b,"Name=");split(b[2],NAME,";");\
#split($9,c,"Classification=");split(c[2],CLASS,";");\
#print $1,$4-1,$5,ID"_"NUM";"NAME[1]";"CLASS[1],$6,$7}}'|\
#sort -k1,1 -k2,2n -k3,3n  > ${Bed}

# 3. filter TE and add ID name - TIR
less ${Gff} | awk -v TE=${TE} -v Class=${Class} 'BEGIN{FS="\t";OFS="\t"}\
{if ($3 == Class ) {split($9,x,"ID=");split(x[2],y,";");split(y[1],z,"_");\
if(z[2]!="homo"){ID=TE;NUM=z[3]}\
else if(z[2]=="homo"){ID="TE_homo";NUM=z[3]};\
split($9,b,"Name=");split(b[2],NAME,";");\
split($9,c,"Classification=");split(c[2],CLASS,";");\
print $1,$4-1,$5,ID"_"NUM";"NAME[1]";"CLASS[1],$6,$7}}'|\
sort -k1,1 -k2,2n -k3,3n  > ${Bed}

# 3. ############################################################################


# 4. make a dictionary and replace fasta
Dict=$(echo "/nfs/project1/Repeat/10_integrate/06_clusterTE/32genomes_Helitron_filter/02_strucTEfa/Dict/strucTE_fasta."$Genome"_"$TE".dict")
less ${Bed} | awk 'BEGIN{FS="\t";OFS="\t"}{split($4,a,";");\
print $1":"$2+1"-"$3"|"a[3], a[1]"::"$1":"$2+1"-"$3}' > $Dict

### NOTE: Need a fresh Fasta file before redo!!!!
Fasta=$(echo "/nfs/project1/Repeat/10_integrate/06_clusterTE/32genomes_Helitron_filter/03_mcl/PolyBgenome20.struc_"$TE_inList".fa")
less $Fasta | \
awk -v Genome=$Genome -v TE=$TE 'BEGIN{FS="\t";OFS="\t"}\
{if(NR==FNR){d[">"$1]=">"$2} else {ID=">"Genome"_";if($0 ~ ID){$0=d[$0]}print}}' $Dict - > ${TE}.tmp


echo "Step3 DONE -- filter TE and add ID name"

# 5. TE_Name of homo and struc list
for i in $( less $Bed | cut -f4 | cut -f2 -d ";" | cut -f1,2 -d "_"| sort | uniq ); do 
	echo $i  
	less $Bed | grep $i | cut -f4 | cut -f1 -d ";" | grep -v "homo"|sort|uniq| tr "\n" ";"
	echo ""
	less $Bed | grep $i | cut -f4 | cut -f1 -d ";" | grep "homo" |sort|uniq| tr "\n" ";"
	echo "" 
done | paste - - - > $Name.$TE.TE_Name.tmp

echo "Step4 DONE -- TE_Name of homo and struc list"


# 6. 80-80-80 rules

less $Name.$TE.TE_Name.tmp | \
awk 'BEGIN{FS="\t";OFS="\t"}{split($2,a,";");\
if($2==""){$2="NOstrucTE"}\
else{for(i in a)if($3==""){$3="NOhomoTE"}} \
print}'| sort | uniq > $Name.$TE.TE_Name.NoID.tmp
# ------------------------------------------------------------------------------- Blastn start for each structural TE

# group1 list: no structural TE --> skip blastn
rm -f $Name.$TE.NOstrucTE.list
for EDTA_TE in $(less $Name.$TE.TE_Name.NoID.tmp | awk 'BEGIN{FS="\t"}\
{if ($2 == "NOstrucTE"){print $3}}' | tr ';' '\n'| awk 'NF>0' | uniq)
do
	echo $EDTA_TE >> $Name.$TE.NOstrucTE.list
done

# group2 and group3 input
less $Name.$TE.TE_Name.NoID.tmp | \
awk 'BEGIN{FS="\t";OFS="\t"}\
{if ($2 != "NOstrucTE" || $3 == "NOhomoTE"){print}}' | \
awk 'NF>0' | sort| uniq > $Name.$TE.TE_Name.NoID.strucTE_NOhomoTE.tmp

mv $Name.$TE.TE_Name.80rule.tmp  $Name.$TE.TE_Name.80rule.tmp_2
rm -f $Name.$TE.step5_error.list
#rm -f $Name.$TE.TE_Name.80rule.tmp $Name.$TE.step5_error.list
for strucTE in $(less $Name.$TE.TE_Name.NoID.strucTE_NOhomoTE.tmp | \
awk 'BEGIN{FS="\t"}{if ($2 != "NOstrucTE"){print $2}}' |cut -f2|tr ';' '\n'|awk 'NF>0'|sort|uniq)
do
	Query=$(echo $Genome"."$strucTE".query.fa.tmp")
	Subject=$(echo $Genome"."$strucTE".subject.fa.tmp")
	
	less $Name.$TE.TE_Name.NoID.strucTE_NOhomoTE.tmp | \
	awk 'BEGIN{FS="\t"}{if ($2 != "NOstrucTE"){print}}' |\
	awk -v strucTE=$strucTE 'BEGIN{FS="\t";OFS="\t"}\
	{split($2,a,";");for(i in a)\
	if(a[i] == strucTE){print}}' > $Name.$TE.TE_Name.NoID.$strucTE.blast.tmp
	
	# make Query (a struc-TE)
	Group=$(less $Name.$TE.TE_Name.NoID.$strucTE.blast.tmp | \
	awk -v strucTE=$strucTE 'BEGIN{FS="\t"}\
	{split($2,a,";");\
	{for(i in a)\
	if(a[i] == strucTE && $3=="NOhomoTE"){print "NOhomoTE"}\
	else if (a[i] == strucTE && $3!="NOhomoTE"){print "YES"}}}'|sort | uniq)
	
	# group2: prepare query for blastn
	if [[ "$Group" = "YES" ]]; then
		less ${TE}.tmp | \
		grep "^>" | sed 's/^>//' | grep $Genome"_" | \
		awk -v strucTE=$strucTE 'BEGIN{FS=":"}\
		{if ($1 == strucTE) print}' > $strucTE.query.id.list.tmp
		
		head -n1 $strucTE.query.id.list.tmp > $strucTE.query.id.list.tmp.tmp
		mv $strucTE.query.id.list.tmp.tmp $strucTE.query.id.list.tmp
		
		perl /nfs/project1/lab_script/group_tools/extract_sequence.pl ${TE}.tmp $strucTE.query.id.list.tmp > $Query
		rm $strucTE.query.id.list.tmp
		
		# make Subject (homo-TEs)
		rm -f $Genome.$strucTE.subject.bed.tmp
		
		for EDTA_TE in $(less $Name.$TE.TE_Name.NoID.strucTE_NOhomoTE.tmp | \
		awk -v strucTE=$strucTE 'BEGIN{FS="\t"}\
		{split($2,a,";");if($3 !="NOhomoTE"){for(i in a)if(a[i] == strucTE){print $3}}}' | \
		tr ';' '\n' | awk 'NF>0' | sort | uniq)
		do
			less ${Bed} | \
			awk -v EDTA_TE=$EDTA_TE 'BEGIN{FS="\t";OFS="\t"}{split($4,a,";");if (a[1] == EDTA_TE) {print}}' | \
			awk 'BEGIN{FS="\t";OFS="\t"}{if($3-$2 >= 80) print}' >> $Genome.$strucTE.subject.bed.tmp
		done
		bedtools getfasta -s -name -fi ${GenomeFasta} -bed $Genome.$strucTE.subject.bed.tmp -fo $Genome.$strucTE.subject.bed.fa.tmp
		less $Genome.$strucTE.subject.bed.fa.tmp | grep "^>" | sed 's/^>//' > $strucTE.subject.id.list.tmp
		perl /nfs/project1/lab_script/group_tools/extract_sequence.pl $Genome.$strucTE.subject.bed.fa.tmp $strucTE.subject.id.list.tmp > $Subject
		rm $Genome.$strucTE.subject.bed.fa.tmp $strucTE.subject.id.list.tmp
		
		# make blast db
		makeblastdb -dbtype nucl -in ${Subject} -out ${Subject}
		
		# Blastn
		blastn -query ${Query} -db ${Subject} \
			-perc_identity 80 \
			-task blastn \
			-num_threads 20 \
			-dust no \
			-max_hsps 1 \
			-soft_masking false \
			-outfmt "7 qseqid sseqid qlen slen length pident" | \
			grep -v "^#" | awk 'BEGIN{FS="\t";OFS="\t"}{if(($5/$4)*100 >= 80) print}' > $Name.$TE.$strucTE.TE_Name.80rule.blast.out.tmp
				
		# revise TE_Name tmp
		less $Name.$TE.$strucTE.TE_Name.80rule.blast.out.tmp | awk 'BEGIN{FS="\t";OFS="\t"}\
		{split($2,EDTA_TEname,";");split(EDTA_TEname[2],EDTA_fam,"_");split($1,strucTE,":");\
		print EDTA_fam[1]"_"EDTA_fam[2],strucTE[1]";"}'| sort | uniq > $strucTE.col_12.tmp
			
		less $Name.$TE.$strucTE.TE_Name.80rule.blast.out.tmp | awk 'BEGIN{FS="\t";OFS="\t"}\
		{split($2,EDTA_TEname,";");split(EDTA_TEname[2],EDTA_fam,"_");\
		print EDTA_TEname[1]}'| tr "\n" ";" > $strucTE.col_3.tmp ; echo "" >> $strucTE.col_3.tmp
		
		paste $strucTE.col_12.tmp $strucTE.col_3.tmp >> $Name.$TE.TE_Name.80rule.tmp
		rm $Query $Subject
		rm ${Subject}.nin ${Subject}.nhr ${Subject}.nsq ${Subject}.ndb ${Subject}.not ${Subject}.nto ${Subject}.ntf  ${Subject}.njs
		rm $Genome.$strucTE.subject.bed.tmp
		rm $Name.$TE.$strucTE.TE_Name.80rule.blast.out.tmp
		rm $Name.$TE.TE_Name.NoID.$strucTE.blast.tmp
		rm $strucTE.col_12.tmp $strucTE.col_3.tmp
		
	# group3: no homo-TE --> skip blastn
	elif [[ "$Group" = "NOhomoTE" ]]; then
		for EDTA_FAM in $(less $Name.$TE.TE_Name.NoID.$strucTE.blast.tmp | awk -v strucTE=$strucTE 'BEGIN{FS="\t"}\
		{split($2,a,";");for(i in a)if (a[i] == strucTE && $3 == "NOhomoTE"){print $1}\
		else print ""}' |sort|uniq)
		do
			echo -e "$EDTA_FAM\t"$strucTE";" > $strucTE.col_12.tmp
			Semi=";"; echo $Semi > $strucTE.col_3.tmp
		done
		
		paste $strucTE.col_12.tmp $strucTE.col_3.tmp >> $Name.$TE.TE_Name.80rule.tmp
		rm $Name.$TE.TE_Name.NoID.$strucTE.blast.tmp
		rm $strucTE.col_12.tmp $strucTE.col_3.tmp
		
	else
		echo ""
		echo "ERROR: Perhaps no $strucTE in TE family -- didnt enter 5 genomes comparison."
		echo ""
		echo $strucTE >> $Name.$TE.step5_error.list
	fi
	
done

echo "Step5 DONE -- 80-80-80 rules"
# ------------------------------------------------------------------------------- END structural TE


# 11. Summary
#       Input_1: Number of all TEs in TE annotation
#       Input_2: Number of struc-TEs in TE annotation
#       Input_3: Number of homo-TEs in TE annotation
#       Input_4: Number of struc-TE in TE family
#       Input_5: Number of TE family

Input_1=$(less $Bed | \
cut -f4 | cut -f1 -d ";" | sort | uniq |wc -l)

Input_2=$(less $Bed | \
grep -v homo | cut -f4 | cut -f1 -d ";" | sort | uniq |wc -l)

Input_3=$(less $Bed | \
grep homo | cut -f4 | cut -f1 -d ";" | sort | uniq |wc -l)

Input_4=$(less $Fasta | \
grep "^>" | grep $Genome"_" | sort | uniq | wc -l)

Input_5=$(less /nfs/project1/Repeat/10_integrate/06_clusterTE/CC/04-parse_id/TE_family.list | \
grep $TE | grep $Genome"_" | cut -f1 -d "_" |sort | uniq |wc -l)


#	Sum_1: struc-TE input
#	Sum_2: homo-TE input (with struc-TE in EDTA_TE_group)
#	Sum_3: Number of TE family
#	Sum_4: Number of all TEs in TE family
#	Sum_5: Number of struc-TE in TE family
#	Sum_6: Number of homo-TE pass 80 rule in TE family
#	Sum_7: Duplicated homo-TE

Sum_1=$(less $Name.$TE.TE_Name.tmp | \
cut -f2 | tr ";" "\n" | awk 'NF>0' | sort | uniq|  grep -v homo|wc -l)

Sum_2=$(less $Name.$TE.TE_Name.tmp | \
awk 'BEGIN{FS="\t"}{if ($2 !="")print $3}' | \
tr ";" "\n" | awk 'NF>0' | sort | grep homo|wc -l)

Sum_3=$(less $Name.$TE.TE_family.list | \
cut -f1 | sort| uniq | grep -v NA | wc -l)

Sum_4=$(less $Name.$TE.TE_family.list | \
cut -f2 | tr ";" "\n" | awk 'NF>0' | sort |uniq|wc -l)

Sum_5=$(less $Name.$TE.TE_family.list | \
cut -f2 | tr ";" "\n" | awk 'NF>0' | sort |grep -v  homo|wc -l)

Sum_6=$(less $Name.$TE.TE_family.list| \
cut -f2 | tr ";" "\n" | awk 'NF>0' | sort |  grep homo| uniq|wc -l)

Sum_7=$(less $Name.$TE.TE_family.list| \
cut -f2 | tr ";" "\n" | awk 'NF>0' | sort |  grep homo| uniq -d |wc -l)

echo $TE > $Name.$TE.TE_family.blast.sum
echo ${Input_1} >> $Name.$TE.TE_family.blast.sum
echo ${Input_2} >> $Name.$TE.TE_family.blast.sum
echo ${Input_3} >> $Name.$TE.TE_family.blast.sum
echo ${Input_4} >> $Name.$TE.TE_family.blast.sum
echo ${Input_5} >> $Name.$TE.TE_family.blast.sum
echo ${Sum_1} >> $Name.$TE.TE_family.blast.sum
echo ${Sum_2} >> $Name.$TE.TE_family.blast.sum
echo ${Sum_3} >> $Name.$TE.TE_family.blast.sum
echo ${Sum_4} >> $Name.$TE.TE_family.blast.sum
echo ${Sum_5} >> $Name.$TE.TE_family.blast.sum
echo ${Sum_6} >> $Name.$TE.TE_family.blast.sum
echo ${Sum_7} >> $Name.$TE.TE_family.blast.sum

rm $Genome.$TE.TE_family.dict.tmp
rm $Name.$TE.TE_Name.NoID.tmp
#rm $Name.$TE.TE_Name.tmp
rm $Name.$TE.TE_Name.TE_family_dup.tmp
rm $Name.$TE.NOstrucTE.list
#rm ${TE}.tmp

echo "Step10 DONE -- Summary"

