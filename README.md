
# What is Yanagi-count?

It is a kmer-based alignment tool with internal support for aligning over transcripts segments created by [Yanagi](https://github.com/HCBravoLab/yanagi). It is a modified version of [RapMap](https://github.com/COMBINE-lab/RapMap/tree/develop-salmon) with the segments support. Rapmap is a testing ground for ideas in quasi-mapping and selective alignment. Currently, RapMap is a stand-alone quasi-mapper that can be used with other tools.  It is also being used as part of [Salmon](https://github.com/COMBINE-lab/salmon) and [Sailfish](https://github.com/kingsfordgroup/sailfish). Yanagi-count, at this point, it is somewhat experimental. In this readme we will use the name RapMap to refer to the aligner itself while using the name Yanagi-count as an alias to RapMap+segments.
# Building Yanagi-count

Same steps to compile RapMap. To build RapMap, you need a C++14 compliant compiler (g++ >= 4.9 and clang >= 3.4) and CMake (>= 3.9).  RapMap is built with the following steps (assuming that `path_to_rapmap` is the toplevel directory where you have cloned this repository):

```
[path_to_rapmap] > mkdir build && cd build
[path_to_rapmap/build] > cmake ..
[path_to_rapmap/build] > make
[path_to_rapmap/build] > make install
[path_to_rapmap/build] > cd ../bin
[path_to_rapmap/bin] > ./rapmap -h
```

This should output the standard help message for rapmap.

# Using Yanagi-count

To use RapMap to map reads, you first have to index your reference transcriptome.  Once the index is created, it can be used to map many different sets of reads.  Assuming that your reference transcriptome segments is in the file `segs.fa` and its meta file is `segs.fa.meta`, you can produce the index as follows:

```
> rapmap quasiindex -i segs_index -t segs.fa --segments segs.fa.meta --keepDuplicates
```

if you want to make use of a minimum perfect hash when indexing (which will lower the memory requirement during mapping), you can add run options `-p -x 4`.
The `-p` option enables the minimum perfect hash and `-x 4` tells RapMap to use up to 4 threads when building the perfect hash (you can specify as many or as few threads as you wish).

The index itself will record whether it was built with the aid of minimum perfect hashing or not, so no extra information concerning this need be provided when mapping.  For the purposes of this example, we'll assume that we wish to map paired-end reads with the first mates in the file `r1.fq.gz` and the second mates in the file `r2.fq.gz`.  We can perform the mapping like so:

```
> rapmap quasimap -i segs_index -s -1 r1.fq.gz -2 r2.fq.gz -o output -t 8 --consensusSlack 0.5 --noOrphans --hardFilter
```

This will tell RapMap to map the paired-end reads using 8 threads, and to write the segment-pair counts into `output.tsv` file.  The `-s` flag tells RapMap to use selective alignment to generate better mappings and to validate the alignment _score_ of hits. 
The command contain optional options that suites the segments situation more rather than RapMap's defaults, specifically the options `--consensusSlack 0.5 --noOrphans --hardFilter`. The `--noOrphans` requires both ends to be mapped otherwise the entire read is discarded. That option is prefered to simplify processing the resulting segment-pair counts. The `--hardFilter` option restrict counting towards the mappings of the highest score.

# Caveats

Yanagi-count the same as RapMap is experimental, and the code, at this point. See [RapMap](https://github.com/COMBINE-lab/RapMap/tree/develop-salmon) for more details.
# License 

See [RapMap](https://github.com/COMBINE-lab/RapMap/tree/develop-salmon) for more details.
