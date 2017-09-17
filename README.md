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

#### 1. Correct errors in long reads using short reads

Lets list whats in our example data directory
```bash
$ docker run -w /Source/IDP/example/Data \
             -h $(pwd):/Source/IDP/example/Data \
             idptest ls -lht
```

>total 209M
-rw-r--r-- 1 root root  44M Sep 17 01:00 lr.fa
-rw-r--r-- 1 root root  48M Sep 17 01:00 sr.fa
-rw-r--r-- 1 root root  26M Sep 17 01:00 uniqueness.chr20.bedGraph
-rw-r--r-- 1 root root  61M Sep 17 01:00 chr20.fa
-rw-r--r-- 1 root root  31M Sep 17 01:00 chr20.gene_est.refFlat.txt
-rw-r--r-- 1 root root 752K Sep 17 01:00 chr20.gpd

There are inside the docker.  We don't have direct access to them, but we can mount them when we use the docker command.

If we were running a real datset we would need the data to be in directories contained in the working directory and to use the `-h` and `-w` to point the docker to the working we are in, which is a bit simpler than the setup for the example.

Lets correct that long read file with the short read file using LoRDEC

```bash
$ docker run -w /Source/IDP/example/Data \
             -h $(pwd):/Source/IDP/example/Data \
             -v $(pwd)/IDP_workflow:/Source/IDP/example/Data/IDP_workflow \
             idptest \
             lordec-correct -i lr.fa -2 sr.fa -k 19 -s 3 \
             -o IDP_workflow/corrected.lr.fa
```

You should get a big progress output from LoRDEC

>[Graph: build branching nodes           ]  100  %   elapsed:   0 min 1  sec    estimated remaining:   0 min 0  sec   cpu:  192.9 %   mem: [ 54,  54, 333] MB 
Graph has 44863 branching nodes.
graph created

The program will look frozen on this step, but its processing, and after minutes (perhaps 10 minutes on a dual core computer) of processing the corrected reads should finish.  You should see the file accumulating in `IDP_workflow` in your working directory.

Lets look at our folder in our working directory `IDP_workflow`

```bash
$ ls -lht IDP_workflow
```

>total 45M
-rw-r--r-- 1 jasonweirather 45M Sep 16 19:55 corrected.lr.fa

Nice there is our corrected short read file.

#### 2. Align the corrected long reads

You could let IDP do this for you, but I caution against it. Its a slow process and the aligners can crash sometimes, so its better to just sort this out now and not deal with it in the IDP run.

First we need a gmap index to use the gmap aligner

```bash
docker run -w /Source/IDP/example/Data \
             -h $(pwd):/Source/IDP/example/Data \
             -v $(pwd)/IDP_workflow:/Source/IDP/example/Data/IDP_workflow \
             idptest \
             gmap_build -d IDP_workflow/gindex -C chr20.fa
```

Building a gmap index is a slow and memory intesive process on a full genome.  This is just part of one chromosome. Lets look and see what it made.

```bash
$ du -hs IDP_workflow/*
```

>45M	IDP_workflow/corrected.lr.fa
1.8M	IDP_workflow/corrected.lr.psl
809M	IDP_workflow/gindex
23M	IDP_workflow/hisat2index.1.ht2
15M	IDP_workflow/hisat2index.2.ht2
4.0K	IDP_workflow/hisat2index.3.ht2
15M	IDP_workflow/hisat2index.4.ht2
25M	IDP_workflow/hisat2index.5.ht2
15M	IDP_workflow/hisat2index.6.ht2
4.0K	IDP_workflow/hisat2index.7.ht2
4.0K	IDP_workflow/hisat2index.8.ht2

Yes it did just make an 809Mb index out of a 60Mb genome file. Lets align to it and get a psl output.

```bash
$ docker run -w /Source/IDP/example/Data \
             -h $(pwd):/Source/IDP/example/Data \
             -v $(pwd)/IDP_workflow:/Source/IDP/example/Data/IDP_workflow \
             idptest \
             gmap -D IDP_workflow -d gindex -t 2 -f 1 \
             IDP_workflow/corrected.lr.fa > IDP_workflow/corrected.lr.psl
```

Again, this is a slow process. If I am doing some serious work I take some time to split the file and parallelize it on a cluster. Even this example will take a few minutes to run. (we dont need to wait for it to finish to start the next step.

#### 3. Align the short reads

I will use `hisat2` to align reads but `runSpliceMap` is included if you want a more classic apporach to the pipeline.

For speed and stability I recommend `hisat2` but it will require an additional processing step on our part.

```bash
$ docker run -w /Source/IDP/example/Data \
             -h $(pwd):/Source/IDP/example/Data \
             -v $(pwd)/IDP_workflow:/Source/IDP/example/Data/IDP_workflow \
             idptest \
             hisat2-build chr20.fa IDP_workflow/hisat2index
```

Lets check the output

```bash
$ ls -lht IDP_workflow/
```

>-rw-r--r--  1 jasonweirather  25M Sep 16 21:08 hisat2index.5.ht2
-rw-r--r--  1 jasonweirather  15M Sep 16 21:08 hisat2index.6.ht2
drwxr-xr-x 15 jasonweirather  510 Sep 16 21:08 gindex
-rw-r--r--  1 jasonweirather  23M Sep 16 21:08 hisat2index.1.ht2
-rw-r--r--  1 jasonweirather  15M Sep 16 21:08 hisat2index.2.ht2
-rw-r--r--  1 jasonweirather   71 Sep 16 21:07 hisat2index.3.ht2
-rw-r--r--  1 jasonweirather  15M Sep 16 21:07 hisat2index.4.ht2
-rw-r--r--  1 jasonweirather   12 Sep 16 21:07 hisat2index.7.ht2
-rw-r--r--  1 jasonweirather    8 Sep 16 21:07 hisat2index.8.ht2
-rw-r--r--  1 jasonweirather  426 Sep 16 21:07 gindex.coords
-rw-r--r--  1 jasonweirather    9 Sep 16 21:07 gindex.sources
-rw-r--r--  1 jasonweirather 1.8M Sep 16 21:05 corrected.lr.psl
-rw-r--r--  1 jasonweirather  45M Sep 16 19:55 corrected.lr.fa
So far so good.  those corrected long reads are still running. :/

Next we actually align the short reads

```bash
docker run -w /Source/IDP/example/Data \
             -h $(pwd):/Source/IDP/example/Data \
             -v $(pwd)/IDP_workflow:/Source/IDP/example/Data/IDP_workflow \
             idptest \
             hisat2 -x IDP_workflow/hisat2index -U sr.fa -f \
| docker run -i idptest samtools view -Sb - \
| docker run -i \
             -w /Source/IDP/example/Data \
             -h $(pwd):/Source/IDP/example/Data \
             -v $(pwd)/IDP_workflow:/Source/IDP/example/Data/IDP_workflow \
             idptest \
             samtools sort -T IDP_workflow/midsort \
             - -o IDP_workflow/sr.sorted.bam
```

>525000 reads; of these:
  525000 (100.00%) were unpaired; of these:
    54 (0.01%) aligned 0 times
    523807 (99.77%) aligned exactly 1 time
    1139 (0.22%) aligned >1 times
99.99% overall alignment rate

Looks good!

Unfortunately IDP needs a different format than the garden variety bam.

To accomodate this we will need to conver the bam into a SpliceMap format sam,
and also create a junction file like SpliceMap does. We use helper scripts for
this part.

First we make the SAM file

```bash
$ docker run -w /Source/IDP/example/Data \
             -h $(pwd):/Source/IDP/example/Data \
             -v $(pwd)/IDP_workflow:/Source/IDP/example/Data/IDP_workflow \
             idptest \
             seq-tools sam_to_splicemap_like_sam IDP_workflow/sr.sorted.bam \
             -o IDP_workflow/sr.splicemap-like.sam
```

> processed 526000 lines. At: chr20:62906952

Next we make the junction file

```bash
$ docker run -w /Source/IDP/example/Data \
             -h $(pwd):/Source/IDP/example/Data \
             -v $(pwd)/IDP_workflow:/Source/IDP/example/Data/IDP_workflow \
             idptest \
             seq-tools sam_to_splicemap_junction_bed \
             IDP_workflow/sr.sorted.bam \
             -r chr20.fa \
             -o IDP_workflow/sr.splicemap-like.junctions.bed
```

>reading reference genome
finished reading reference genome
reading through sam file
WARNING skipping non-canonical splice (11/31760)
finished reading sam

#### 4. Run IDP

The psl option is the most convenient way to run IDP since it allows you to do
your own alignment ahead of time as we have done here. 

To make this easier the IDP/examples folder contains a configuration file
that points to the folders we've generated in this example

The configuration file is pretty key so lets move it where you can see it, in IDP_workflow
```bash
$ docker run -w /Source/IDP/example/Data \
             -h $(pwd):/Source/IDP/example/Data \
             -v $(pwd)/IDP_workflow:/Source/IDP/example/Data/IDP_workflow \
             idptest \
             cp ../run_LRpsl.cfg IDP_workflow/run_LRpsl.cfg
```

On a normal run you will create your own configuration file to describe the run.

Now to actually run IDP.  This configuration file has been set to put the outputs into `IDP_workflow/output_psl`.

```bash
$ docker run -w /Source/IDP/example/Data \
             -h $(pwd):/Source/IDP/example/Data \
             -v $(pwd)/IDP_workflow:/Source/IDP/example/Data/IDP_workflow \
             idptest \
             runIDP.py IDP_workflow/run_LRpsl.cfg 0
```
