# MaskJelly

A lightweight tool for adding spaced seeds to a list of k-mers produced by Jellyfish or similar tools.
To install `maskjelly',
```sh
git clone https://github.com/harti768/MaskJelly.git
cd MaskJelly/build; make
```

## Requirements

MaskJelly by itself does not require any other libraries.
For the creation of an input file, [Jellyfish](https://github.com/gmarcais/Jellyfish) is required.

## Input Files
The list of k-mers needs to be in .fasta format, where every odd line specifies the abundacy of the k-mer in the following line. For example:
```sh
>10
AAAAAAAAAAA
>3
ACTGTGGTGTG
```

The spaced seed must be specified in a separate file and needs to consist of 1 (care) and 0 (don't care) positions.

### Prepare Data

* Create a list of k-mers with Jellfish:

        jellyfish count -m 31 -s 100M -C -o kmers.jf reads.fa
        jellyfish dump reads.fa > kmers.fa
* Create a mask:

        echo "1101101101101101011011011011011" > mask.txt


## MaskJelly Examples

* Apply mask to list

        maskjelly -i kmers.fa -m mask.txt > masked_kmers.fa

* Specifiy output file

        maskjelly -i kmers.fa -m mask.txt -o masked_kmers.fa

* Specify number of threads

        maskjelly -i kmers.fa -m mask.txt -o masked_kmers.fa -t 4

* Full pipeline

        jellyfish count -m 31 -s 100M -C -o kmers.jf reads.fa
        jellyfish dump reads.fa > kmers.fa
        maskjelly -i kmers.fa -m mask.txt -t 4 > masked_kmers.fa
        jellyfish count -m 21 -s 100M -C -o masked_kmers.jf masked_kmers.fa
