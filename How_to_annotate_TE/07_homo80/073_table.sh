# USAGE: sh $0
# input: All stat_tables in current directory 
#	stat_table is the output of "/nfs/project1/Repeat/10_integrate/00_stat_TEclass_gff.sh"

cd stat/

# create each table 
for Genome in $(less ../genome.list ); do 
Name=$(less /homes/kochiaying/name32.list | grep -w $Genome | cut -f7 -d "/")
rm -f $Name.table.tmp  
echo -e $Genome"\t"$Genome"\t"$Genome"\t"$Genome > $Name.table.tmp
less  $(ls -1 ./ | grep "stat_"${Name}".table" )  >> $Name.table.tmp 
done

# paste all tables
#paste $(ls -1 *.table.tmp | tr "\n" " " ) | sed 's/genome/genome\t100%\t/g' | \
paste $(ls -1 *.table.tmp | tr "\n" " " ) | \
awk 'BEGIN{FS=OFS="\t"}\
{if($1!="Class-I:SINE"&&$1!="Class-II:DNA"&&$1!="Class-I:LINE:CR1"&&$1!="Class-I:LINE:L1"&&$1!="Class-I:LINE:unknown"&&$1!="Class-I:SINE:tRNA-CR1"&&$1!="Class-I:SINE:unknown"&&$1!="Class-II:TIR:Mu"&&$1!="Class-II:TIR:unknown"&&$1!="Class-II:DNA:Kolobok-H"&&$1!="Class-II:DNA:Maverick"&&$1!="Class-II:DNA:TRIM"&&$1!="Class-II:TIR:Tc1_Mariner"&&$1!="Class-II:TIR:Merlin"&&$1!="Class-II:TIR:hAT"&&$1!="Class-II:DNA:PIF_Harbinger"&&$1!="Class-II:DNA:Penelope"&&$1!="pararetrovirus"&&$1!="Class-II:DNA:unknown"&&$1!="Unknown")\
print}' > TEtable.txt

rm *.table.tmp
