#Count number of k-mers
jellyfish histo masked_kmer_reads.jf | awk '{print $1*$2}' | awk '{s+=$1} END {print s}'

#masking k-mers results in less k-mers overall, because a larger k-mer size is used for the mask
#However, the histogram of the k-mers shows, that the number of unique k-mers increases
#whereas the number of non-unique k-mers decreases
#==> mask catches more variation (Probably)