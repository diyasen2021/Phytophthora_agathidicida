awk '{for(i=1; i<= $3; i=i+10000) if(i < $3 - 10000) print $1"\t"i"\t"(i+10000); else print $1"\t"i"\t"$3}' length.txt > bins10000.bed # create bins from chromosome length file
bedtools nuc -fi ../rawdata/3770/P_agathidicida_3770_56996300_v2.fna -bed bins10000.bed > bedtoolsout.bed #calculate GC for bins
cut -f 1,2,3,5 bedtoolsout.bed > gc10000.bed # final file
