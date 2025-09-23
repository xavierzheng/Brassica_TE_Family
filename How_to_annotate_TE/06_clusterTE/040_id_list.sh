#cd 04-parse_id

#rm -f TE_family.list2

#for i in Helitron; do
#less /nfs/project1/Repeat/10_integrate/06_clusterTE/CC/03-mcl/I_20/dump.out.all5.struc_${i}.blast_out.I20 | \
#awk '{FS="\t";OFS="\t"}\
#{for(i=1;i<=NF;i++)\
#{split($i,a,":");split(a[1],b,"/");split(a[3],c,"_");\
#print b[2]sprintf("%04d",NR)"_"c[1]"_"i,$i}}' >> TE_family.list2
#done

#for i in DTA DTH DTM DTT ; do
#less /nfs/project1/Repeat/10_integrate/06_clusterTE/CC/03-mcl/I_15/dump.out.all5.struc_${i}.blast_out.I15 | \
#awk '{FS="\t";OFS="\t"}\
#{for(i=1;i<=NF;i++)\
#{split($i,a,":");split(a[1],b,"/");split(a[3],c,"_");\
#print b[2]sprintf("%04d",NR)"_"c[1]"_"i,$i}}' >> TE_family.list2
#done

#for i in CACTA LTRCO LTRGY LTRRT ; do
#less /nfs/project1/Repeat/10_integrate/06_clusterTE/CC/03-mcl/I_15/dump.out.${i}.TEvsGENE.blast_out.I15 | \
#awk '{FS="\t";OFS="\t"}\
#{for(i=1;i<=NF;i++)\
#{split($i,a,":");split(a[1],b,"-");split(a[3],c,"_");split(b[3],d,"_");\
#print d[2]sprintf("%04d",NR)"_"c[1]"_"d[3],$i}}' >> TE_family.list2
#done
#cd ../

#-------------------------

cd 04-parse_id

rm -f TE_family.list

for i in Helitron DTA DTH DTM DTT ; do 

less /nfs/project1/Repeat/10_integrate/06_clusterTE/CC/03-mcl/I_15/dump.out.all5.struc_${i}.blast_out.I15 | \
awk '{FS="\t";OFS="\t"}\
{for(i=1;i<=NF;i++)\
{split($i,a,":");split(a[1],b,"/");split(a[3],c,"_");\
print b[2]sprintf("%04d",NR)"_"c[1]"_"i,$i}}' >> TE_family.list

done

for i in CACTA LTRCO LTRGY LTRRT ; do 

less /nfs/project1/Repeat/10_integrate/06_clusterTE/CC/03-mcl/I_15/dump.out.${i}.TEvsGENE.blast_out.I15 | \
awk '{FS="\t";OFS="\t"}\
{for(i=1;i<=NF;i++)\
{split($i,a,":");split(a[1],b,"-");split(a[3],c,"_");split(b[3],d,"_");\
print d[2]sprintf("%04d",NR)"_"c[1]"_"d[3],$i}}' >> TE_family.list

done
cd ../
