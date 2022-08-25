# MaskJelly

A lightweight tool for counting spaced seed k-mers in combination with Jellyfish.
To install `maskjelly`,
```sh
git clone https://github.com/hhaentze/MaskJelly.git
cd MaskJelly/build; make
```

## Requirements

MaskJelly by itself does not require any other libraries.
For the creation of a compressed k-mer-count file, [Jellyfish](https://github.com/gmarcais/Jellyfish) is recommended.


## Input Files
MaskJelly can mask the k-mers either directly from a fasta file (recommended) or, if only .jf files are given, by using the output of `jellyfish dump`.

For the second case the list of k-mers needs to be in .fasta format, where every odd line specifies the abundacy of the k-mer in the following line. For example:
```sh
>10
AAAAAAAAAAA
>3
ACTGTGGTGTG
```

The mask / spaced seed must be specified in a separate file and needs to consist of 1 (care) and 0 (don't care) positions.


        echo "1101101101101101011011011011011" > mask.txt


## MaskJelly Examples

* Apply mask to a fasta file and create jellyfish file

        maskjelly -i <reads.fa> -m <mask.txt> -f | jellyfish count -m 21 -s 100M -C /dev/stdin

* Specifiy output file (not recommended, as file can get quite large)

        maskjelly -i <reads.fa> -m <mask.txt> -o masked_kmers.fa

* Apply mask to a jellyfish file

        jellyfish dump <kmers.jf> | maskjelly -m <mask.txt> | jellyfish count -m 21 -s 100M -C /dev/std


## Arguments
* -m :    Mask consisting of care (1) and don't care (0) positions. Needs to be saved as file. (mandatory)
* -i :    Input file
* -l :    Limit highest abundance of k-mers. Higher values will not be put out. (Only works for output of 'jellyfish dump')
* -o :    Output file
* -t :    Number of threads
* -f :    (flag) Input is fasta file