# Kmer count
# kmer 21, 16 threads, 64G of memory, counting kmers between 1 and 10000x
#kmc -k21 -t16 -m64 -ci1 -cs10000 @samplefiles sample_counts tmp # count kmers
#kmc_tools transform 3770_counts histogram sample_k21.hist -cx10000 # makes histogram in *k21.hist

# Run smudgeplot
# conda activate smudgeplot
for file in *hist
do
sample=`ls $file | sed s'/_k21.hist//'`
L=$(smudgeplot.py cutoff $sample\_k21.hist L)
U=$(smudgeplot.py cutoff $sample\_k21.hist U)
kmc_tools transform $sample\_counts -ci"$L" -cx"$U" dump -s $sample\_L"$L"_U"$U".dump
smudgeplot.py hetkmers -o $sample_L"$L"_U"$U" < $sample\_L"$L"_U"$U".dump
smudgeplot.py plot -o $sample -t "$sample" -q 0.99 $sample\_kmer_pairs_coverages.tsv
done

