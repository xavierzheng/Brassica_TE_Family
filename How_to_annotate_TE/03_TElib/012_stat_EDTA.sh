for ID in $(less ~/name30.list | cut -f2); do Name=$(less ~/name30.list|grep -w $ID|cut -f7 -d "/");Fasta=$(less ~/name30.list|grep -w $ID|cut -f1|cut -f8 -d "/");Genome=$(less ~/name30.list|grep -w $ID|cut -f1); Gff=$(echo $Name/$Fasta.mod.EDTA.TEanno.gff3);less $Gff | grep -v "^#" | sed 's/Classification=MITE\//Classification=DNA\//;s/Classification=DNA\//Classification=TIR\//;s/Classification=TIR\/Helitron/Classification=DNA\/Helitron/;s/Classification=TIR\/PIF_Harbinger/Classification=DNA\/PIF_Harbinger/;s/Classification=Unknown;/Classification=unknown;/;s/Classification=Maverick/Classification=DNA\/Maverick/;s/Classification=Penelope/Classification=DNA\/Penelope/;s/Classification=TIR\/EnSpm_CACTA/Classification=TIR\/DTC/' > $Name.gff.tmp;sh /nfs/project1/Repeat/10_integrate/00_stat_TEclass_gff.sh $Name.gff.tmp  $Genome.fai > stat/stat_$Name.table 
done

sh /nfs/project1/Repeat/10_integrate/002_stat_TEtable.sh


