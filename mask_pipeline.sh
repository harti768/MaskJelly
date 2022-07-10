READS="test-reads.fa"
SEGMENTS="k41_path_segments.fasta"
REF="test-reference.fa"
VAR="test-variants.vcf"

#Create segment file
echo "Create segmented genome"
PanGenie-graph -r $REF -v $VAR -k 41 -o k41

#Spaced k-mers for reads
echo "Reads: Count k-mers"
jellyfish count -m 41 -s 3G -p 126 -c 7 -C -L 1 -t 10 --if $SEGMENTS $READS -o .raw_kmers.jf
jellyfish dump .raw_kmers.jf > mer_counts.txt
./mask_mers 
jellyfish count -m 31 -s 3G -C -t 10 masked_mer_counts.txt -o masked_kmer_reads.jf

#Spaced k-mers for pangenome
echo "Genome: Count k-mers"
jellyfish count -m 41 -s 3G -C -t 10 $SEGMENTS -o .raw_kmers.jf
jellyfish dump .raw_kmers.jf > mer_counts.txt
./mask_mers 
jellyfish count -m 31 -s 3G -C -t 10 masked_mer_counts.txt -o masked_kmer_genome.jf

#Clean
rm mer_counts.txt masked_mer_counts.txt .raw_kmers.jf