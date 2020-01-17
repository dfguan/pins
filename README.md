# Pins (拼接)

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
### Scaffolding with linked reads
#### Linked-read preprocessing
Given a list **10xlist** of 10x read files (suppose in fastq.gz format, paired files in a line) and the assembly **asm**, use the following code to get read files alignments. 

```
bwa index $asm
while read -r r1 r2
do
	prefix=`basename $r1 .fastq.gz`
	10x_trim -c -p $prefix $r1 $r2 # generate trimmed read files $prefix_{1,2}.fq.gz
	bwa mem -t12 $asm $prefix_1.fq.gz $prefix_2.fq.gz | samtools view -b - > $prefix.bam
done < $10xlist
```
#### Linked-read scaffolding
With linked-read alignment files **bams** and the draft assembly ***asm***, you can run the following code build scaffolds.

```
samtools faidx $asm
./bin/pin_10x link $bam1 $bam2 $bam3 ... > link.matrix  # this will calcuate the numbers of shared barcode between pairs of contigs.
./bin/pin_10x build -c $asm.fai link.matrix > scaffolds.sat # this will generate scaffolding paths. 
./bin/pin_10x gets -c $asm scaffolds.sat > scaffolds.fa # this will generate scaffolds by a given SAT file.
```


### Scaffolding with Hi-C reads
##### Hi-C Read preprocessing
Given a list **hiclist** of Hi-C read files (suppose in fastq.gz format, paired files in a line) and the assembly **asm**, use the following code to generate Hi-C alignment files. 

```
bwa index $asm
while read -r r1 r2
do
	prefix=`basename $r1 .fastq.gz`
	bwa mem -SP -B10 -t12 $asm $prefix_1.fq.gz $prefix_2.fq.gz | samtools view -b - > $prefix.bam
done < $hiclist
```

##### Hi-C scaffolding


Given Hi-C reads alignment **bams**, a draft assembly **asm** and a output directory **outdir**, if you want to build scaffols with Hi-C in **N** (default: 3) rounds, please try the following commands. The final assembly will be named as **scaffols_final.fa**.

```
samtools faidx $asm 
./bin/pin_hic_it -i $N -c $asm.fai -x $ref -O $outdir $bam1 $bam2 $bam3 ... 
```

Or you want to build scaffolds step by step:
##### Step 1. contact matrix calculation
From a draft assembly：

```
samtools faidx $asm
./bin/pin_hic link $bam1 $bam2 $bam3 ... > link.matrix  # this will calcuate contact numbers between any pairs of contigs.
```

From a **sat** file：

```
./bin/pin_hic link -s $sat $bam1 $bam2 $bam3 ... > link.matrix  # this will calcuate contact numbers between any pairs of contigs.
```

##### Step 2. Scaffolding graph construction
From a draft assembly:

```
/bin/pin_hic build -w100 -k3 -c $asm.fai link.matrix > scaffolds.sat # this will generate scaffolding paths. 
```

From a **sat** file:

```
/bin/pin_hic build -w100 -k3 -s $sat link.matrix > scaffolds.sat # this will generate scaffolding paths. 
```

##### Step 3. Mis-join detection
Given a **sat** file:

```
./bin/pin_hic break $sat $bam1 $bam2 $bam3 ... > scaffs.bk.sat
./bin/pin_hic gets -c $asm scaffs.bk.sat > scaffols_final.fa # get scaffold sequences.
```

A scaffolding pipeline of 3 iterations:

```
samtools faidx $asm
for i in `seq 1 3`
do
	if [ $i -eq 1 ]
	then 
		./bin/pin_hic link $bam1 $bam2 $bam3 ... > links_$i.matrix
		./bin/pin_hic build -w100 -k3 -c $asm.fai links_$i.matrix > scaffolds_$i.sat
	else
		./bin/pin_hic link -s scaffolds_$pi.sat $bam1 $bam2 $bam3 ... > links_$i.matrix
		./bin/pin_hic build -w100 -k3 -s scaffolds_$pi.sat links_$i.matrix > scaffolds_$i.sat 
	fi
	pi=i
done
./bin/pin_hic break -s scaffolds_$i.sat $bam1 $bam2 $bam3 ... > scaffolds_bk.sat 
./bin/pin_hic gets -c $asm scaffs.bk.sat > scaffols_final.fa 
```


### Scaffolding with linkage disequilibrium information
shall be updated soon...


### Output format: SAT (need to be updated)
SAT format is extended from the [GFA 1.0 format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md).

| Tag | Col 1 | Col 2 | Col 3 | Col 4 | Col 5 | Col 6 | Col 7 | Comment | 
|---|---|---|---|---|---|---|---|---| 
| Header | H | VN:Z:1.0 | 
| Sequence | S | SID | LEN | SEQ | 
| Link | L | SID1 | ORI1 | SID2 | ORI2 | CIGAR | wt:f:x | 
| Path | P | PID | LEN | CONTIGS | 
| Scaffold set | C | SCFID | PIDs | 
| Current scaffold set | A | SCFID | 


## Limitation


## FAQ



## Contact

Wellcome to use and distribute the package. Please use the github webpage to report an issue or email dfguan9@gmail.com with any advice. 
