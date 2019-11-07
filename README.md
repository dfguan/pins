# Pins 

Scaffolding tool based on linked reads, Hi-C reads and linkage disequilibrium information. 

## Overview

**pins** is a scaffolding toolkit, it contains three basic programs, namely, pin\_ld a scaffolder based on linkage disequalibrium inforamtion, pin\_hic a scaffolder using Hi-C data, and pin\_10x a scaffolder based on linked reads. It applies a dual selection and local optimal strategy to bridge two contigs and output a SAT file for each iteration, the SAT format is the extension of GFA format which is able to record the scaffolding process, and can also be useful for furthure genomic analysis. **please change the text**


## Dependencies

1. zlib
2. htslib 



## Installation
Run the following commands to intall pins:

```
git clone https://github.com/dfguan/pins.git
cd pins/src && make

```

## Usage

### scaffolding with linked reads
Given linked reads alignment ***bams*** and a draft assembly ***ref***, you can follow the instruction to build scaffolds from contigs.

```
./bin/pin_10x link $bam1 $bam2 $bam3 ... > lk.matrix  # this will calcuate linkage between any pairs of contigs.
./bin/pin_10x build -c $ref.fai lk.matrix > lk.sat # this will generate scaffolding paths. 
./bin/pin_10x getc -c $ref lk.sat > lk.fa # this will generate scaffolds by a given SAT file.
```

### scaffolding with Hi-C reads
Given linked reads alignment ***bams*** and a draft assembly ***ref***, you can follow the instruction to build scaffolds from contigs.

```
./bin/pin_hic link $bam1 $bam2 $bam3 ... > hc.matrix  # this will calcuate linkage between any pairs of contigs.
./bin/pin_hic build -c $ref.fai hc.matrix > hc.sat # this will generate scaffolding paths. 
./bin/pin_hic getc -c $ref hc.sat > hc.fa # this will generate scaffolds by a given SAT file.
```
If you want to build scaffols with Hi-C in *N* rounds, try to run the following command. The final assembly will be named as **scaffols_final.fa**.

```
./bin/pin_hic_it -i $N -c $ref.fai -x $ref $bam1 $bam2 $bam3 ... 
```


## Limitation


## FAQ



## Contact

Wellcome to use and distribute the package. Please use the github webpage to report an issue or email dfguan9@gmail.com with any advice. 
