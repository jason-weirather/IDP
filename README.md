IDP
===

IDP is an gene Isoform Detection and Prediction tool from Second Generation Sequencing and PacBio sequencing developed by Kin Fai Au. It offers very reliable gene isoform identification with high sensitivity.

This is a fork of the original IDP. I trained as a postdoc under the author Kin Fai Au, and have occasionally helped people running this software. The purpose of this fork is to make a docker image to aid users who may have trouble setting up their environment. The official distribution is available here:

https://www.healthcare.uiowa.edu/labs/au/IDP/

## How to run the complete IDP pipeline with docker
### (from long and short reads with a reference genome to isoforms)
The pipeline is fairly involved but the `IDP/example/Data` folder provides a sufficient example to explain the pipline.

The first step is to perform error correction on long reads using long and short reads combined.

We have included both LoRDEC and LSC that can accomplish this.  I recommend LoRDEC for speed and comparable performance on larger datasets.

#### 0. Setup the example

First lets get the example.  You should just clone this git repository and use the example data from it. You can clone this repository and we'll work in the example directory.

```bash
$ git clone https://github.com/jason-weirather/IDP.git
$ cd IDP/example
$ gunzip data/*.gz
$ ls -lht data
```
>total 400M
>-rw-r--r-- 1 jasonweirather 161M Sep 17 20:29 lr.fa
>-rw-r--r-- 1 jasonweirather  31M Sep 17 20:21 gene_est.refFlat.txt
>-rw-r--r-- 1 jasonweirather 752K Sep 17 20:17 reference.gpd
>-rw-r--r-- 1 jasonweirather 768K Sep 17 20:16 exp_len.txt
>-rw-r--r-- 1 jasonweirather  61M Sep 17 20:13 chr20.fa
>-rw-r--r-- 1 jasonweirather 148M Sep 17 20:06 sr.fa

#### 1. Correct errors in long reads using short reads

Because of legacy code in the LSC-IDP pipeline we need the long reads in a very specific type of PacBio read name format.  We can convert any fastq or fasta file with that with an included utility in `seq-tools rename_to_pacbio`.  I included this step in case you are using another read type like Iso-seq or Oxford Nanopore reads.

```bash
$ docker run -v $(pwd):/home vacation/idp \
             seq-tools rename_to_pacbio --fasta data/lr.fa \
             --output_table name_conv.txt > lr.pbname.fa
```
`name_conv.txt` contains the new name of each read in case you need to recover that later.

Lets correct that long read file with the short read file using LoRDEC. LSC is also available, but LoRDEC offers far faster speeds and is less resource intensive and has comparable error correction.

```bash
$ docker run -v $(pwd):/home vacation/idp \
             lordec-correct -i lr.pbname.fa -2 data/sr.fa -k 19 -s 3 \
             -o corrected.lr.fa
```

You should get a big progress output from LoRDEC
>[Graph: build branching nodes           ]  100  %   elapsed:   0 min 1  sec    estimated remaining:   0 min 0  sec   cpu:  188.5 %   mem: [ 73,  73, 355] MB 
>Graph has 50402 branching nodes.

The program will look frozen on this step, but its processing, and after minutes (perhaps 10 minutes on a dual core computer) of processing the corrected reads should finish.  You should see the file accumulating in `corrected.lr.fa`.

#### 2. Align the corrected long reads

You could let IDP do this for you, but I caution against it. Its a slow process and the aligners can crash sometimes, so its better to just sort this out now and not deal with it in the IDP run.

First we need a gmap index to use the gmap aligner

```bash
$ docker run -v $(pwd):/home vacation/idp \
             gmap_build -d ./gmapindex data/chr20.fa
```

>Byte-coding: 62807695 values < 255, 217826 exceptions >= 255 (0.3%)
>Writing file ./gmapindex/gmapindex.salcpchilddc...done
>Found 3535448 exceptions

Building a gmap index is a slow and memory intesive process on a full genome.  This is just two chromosomes.
Now lets make the alignment.

We will only take the best alignment since IDP will only use the best alignment.

```bash
$ docker run -v $(pwd):/home vacation/idp \
             gmap --ordered -D ./ -d gmapindex -t 2 -f 1 -n 1 \
             corrected.lr.fa > corrected.lr.psl
```

>Processed 1513 queries in 408.50 seconds (3.70 queries/sec)
>Removed existing memory for shmid 32769
>Removed existing memory for shmid 0

Again, this is a slow process. If I am doing some serious work I take some time to split the file and parallelize it on a cluster. Even this example will take a few minutes to run. (we dont need to wait for it to finish to start the next step.)

#### 3. Align the short reads

I will use `hisat2` to align reads but `runSpliceMap` is included if you want a more classic apporach to the pipeline.

For speed and stability I recommend `hisat2` but it will require an additional processing step on our part.

```bash
$ mkdir hisat2
$ docker run -v $(pwd):/home vacation/idp \
             hisat2-build data/chr20.fa hisat2/hisat2index
```
>Total time for call to driver() for forward index: 00:01:27

Next we actually align the short reads

```bash
docker run -v $(pwd):/home vacation/idp \
             hisat2 -x hisat2/hisat2index -U data/sr.fa -f \
| docker run -i vacation/idp \
             samtools view -Sb - \
| docker run -i -v $(pwd):/home vacation/idp \
             samtools sort -T midsort - -o sr.sorted.bam
```

>1600000 reads; of these:
>  1600000 (100.00%) were unpaired; of these:
>    48 (0.00%) aligned 0 times
>    1591463 (99.47%) aligned exactly 1 time
>    8489 (0.53%) aligned >1 times
>100.00% overall alignment rate

Looks good!

Unfortunately IDP needs a different format than the garden variety bam.

To accomodate this we will need to conver the bam into a SpliceMap format sam,
and also create a junction file like SpliceMap does. We use helper scripts for
this part.

First we make the SAM file

```bash
$ docker run -v $(pwd):/home vacation/idp \
             seq-tools sam_to_splicemap_like_sam sr.sorted.bam \
             -o sr.splicemap-like.sam
```

>processed 1610000 lines. At: chr20:62907177

Next we make the junction file

```bash
$ docker run -v $(pwd):/home vacation/idp \
             seq-tools sam_to_splicemap_junction_bed sr.sorted.bam \
             -r data/chr20.fa -o sr.splicemap-like.junctions.bed
```

>WARNING skipping non-canonical splice (5/187219)
>finished reading sam

#### 4. Run IDP

The psl option is the most convenient way to run IDP since it allows you to do
your own alignment ahead of time as we have done here. 

To make this easier the IDP/examples folder contains a configuration file
that points to the folders we've generated in this example

On a normal run you will create your own configuration file to describe the run.

Now to actually run IDP.  This configuration file has been set to use the files created in this example.

In this example we are using an RPKM absolute and fraction cutoff rather than an FDR.  The FDR does not execute well in small datasets or nonmodel organisms.

```bash
$ docker run -v $(pwd):/home vacation/idp \
             runIDP.py run.cfg 0
```

Lets see if we have an output .. I didn't wait for all the long reads to correct and align here.

```bash
$ wc -l output/*
```

>   56 output/isoform.exp
>   72 output/isoform.gpd
>   65 output/isoform_detection.gpd
>   20 output/isoform_prediction.gpd
>  213 total
