---
layout: main
title: DADA2 example workflow
categories: [amplicon, tutorial]
permalink: amplicon/dada2_workflow_ex
---

{% include _amplicon_ex_workflow_toc.html %}

{% include _side_tab_amplicon.html %}

<br>
<center><img src="{{ site.url }}/images/under_construction.jpeg"></center>
<center><h3>UNDER CONSTRUCTION</h3></center>
<br>
<br>

[DADA2](https://benjjneb.github.io/dada2/index.html){:target="_blank"} is a relatively new processing workflow for recovering single-nucleotide resolved Amplicon Sequence Variants (ASVs) from amplicon data – if you're unfamiliar with ASVs, you can read more about ASVs vs OTUs in the [opening caveats](/amplicon/dada2_workflow_ex#opening-caveats) section that follows. Developed and maintained by [@bejcal](https://twitter.com/bejcal){:target="_blank"} et al., DADA2 leverages sequencing quality and abundance information to a greater extent than previously developed tools in order to generate an error model based on your actual data. It then uses this error model to do its best to infer the true biological sequences that generated your data. The original paper can be found [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/){:target="_blank"}, and the DADA2 R package site is [here](https://benjjneb.github.io/dada2/index.html){:target="_blank"}. The DADA2 team already has a great tutorial available [here](https://benjjneb.github.io/dada2/tutorial.html){:target="_blank"}, and I learned my way around DADA2 by following that and reading through the [manual](https://www.bioconductor.org/packages/3.3/bioc/manuals/dada2/man/dada2.pdf){:target="_blank"}. While the overall workflow presented here is the same, I've added 2 main things (in addition to all of my rambling of course): 1) example code for how to generate the standard output files of amplicon processing from the final R objects that DADA2 creates (i.e. a fasta of ASVs, a count table, and a taxonomy table); and 2) a [section](/amplicon/dada2_workflow_ex#16s-and-18s-mixed-together) for people working with 16S and 18S sequences mixed together. 
<br>
<br>

---
---
<br>

## Opening caveats
There are many ways to process amplicon data. Some of the most widely used tools/pipelines include [mothur](https://www.mothur.org/){:target="_blank"}, [usearch](https://drive5.com/usearch/){:target="_blank"}, [vsearch](https://github.com/torognes/vsearch){:target="_blank"}, [Minimum Entropy Decomposition](http://merenlab.org/2014/11/04/med/){:target="_blank"}, [DADA2](https://benjjneb.github.io/dada2/index.html){:target="_blank"}, and [qiime2](https://qiime2.org/){:target="_blank"} (which employs other tools within it). If you are looking solely at a broad level, you will likely get very similar results regardless of which tool you use so long as you make similar decisions when processing your sequences (e.g. decisions about things like minimum abundance filtering). But there is a movement in the community away from the traditional OTU approach and on to single-nucleotide-resolving methods that generate what are called ASVs (amplicon sequence variants). And the reasoning for this is pretty sound, as recently laid out very nicely by [Callahan et al. here](https://www.nature.com/articles/ismej2017119){:target="_blank"}, but the intricacies of the differences may seem a little nebulous at first if you're not used to thinking about these things yet. If you are new to this, know that most of the experts on these things would recommend using a method that resolves ASVs. It may be the case that you'd like a broader level of resolution than that, but it is still best to first generate ASVs and then you can always "back out" your resolution with clustering at some threshold or binning sequences by taxonomy or whatever. This is to say that clustering OTUs isn't inherently a problem, it's how those OTUs are generated that could bring you into the fray. 

Keep in mind here that just as was the case with the [usearch/vsearch example workflow](/amplicon/workflow_ex){:target="_blank"}, none of this is meant to be authoritative. This is simply one example of one workflow. When working with your own data you should of course never follow any pipeline blindly.  
<br>

---
<br>

# Tools used here
So here we'll be using [DADA2](https://benjjneb.github.io/dada2/index.html){:target="_blank"} as our main processing tool. Additionally we'll be using [Brian Bushnell's](https://twitter.com/BBToolsBio){:target="_blank"} very handy [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/){:target="_blank"}, specifically [bbduk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/){:target="_blank"}, to remove our primers and perform quality trimming. DADA2 is available as an [R](https://www.r-project.org/){:target="_blank"} package, with installation instructions provided by the developers [here](https://benjjneb.github.io/dada2/dada-installation.html){:target="_blank"}. If your experience is like mine, it shouldn't give you any trouble if installing on your computer, but you may run into some issues when trying to install on a server where you don't have authorization to do whatever you'd like (a huge thanks to [@phantomBugs](https://twitter.com/phantomBugs){:target="_blank"} for all his help when I was bugging him about that 🙂). If you'd like to use [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/){:target="_blank"}, you can find installation instructions for it from the developers [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/) and another example installation [here](/bash/installing_tools#bbtools){:target="_blank"}.
 
And since dada2 is an R package, you also of course need to have a working installation of R on your computer. If you'd like more info on this, check out the [R basics](/R/basics){:target="_blank"} section before moving forward. 
<br>
<br>

---
<br>

# The data
For a quick overview of the example data we'll be using and where it came from, we are going to work with a subset of the dataset [published here](https://www.frontiersin.org/articles/10.3389/fmicb.2015.01470/full){:target="_blank"}. We were exploring an underwater mountain ~3 km down at the bottom of the Pacific Ocean that serves as a low-temperature (~5-10°C) hydrothermal venting site. This amplicon dataset was generated from DNA extracted from crushed basalts collected from across the mountain with the goal being to begin characterizing the microbial communities of these deep-sea rocks. No one had ever been here before, so as is often the purpose of marker-gene sequencing, this was just a broad-level community survey. The sequencing was done on the Illumina MiSeq platform with 2x300 paired-end sequencing using primers targeting the V4 region of the 16S rRNA gene. There are 20 samples total: 4 extraction "blanks" (nothing added to DNA extraction kit), 2 bottom-water samples, 13 rocks, and one biofilm scraped off of a rock. None of these details are important for you to remember, it's just to give some overview if you care.  

In the following figure, overlain on the map are the rock sample collection locations, and the panes on the right show examples of the 3 distinct types of rocks collected: 1) basalts with highly altered, thick outer rinds (>1 cm); 2) basalts that were smooth, glassy, thin exteriors (~1-2 mm); and 3) one calcified carbonate.

<center><img src="{{ site.url }}/images/dorado.png"></center>

<br>
You can download the required dataset and files by copying and pasting the following commands into your terminal. For speed purposes we're only going to work with about 10% of the full dataset. Altogether the uncompressed size of the working directory is < 300MB.

```
cd ~
curl -O https://AstrobioMike.github.io/tutorial_files/dada2_amplicon_ex_workflow.tar.gz
tar -xzvf dada2_amplicon_ex_workflow.tar.gz
rm dada2_amplicon_ex_workflow.tar.gz
cd dada2_amplicon_ex_workflow/
```

Now, let's get started!
<br>
<br>

---
<br>

# Preprocessing
It's good to try to keep a bird's-eye view of what's going on. So here is an overview of the main processing steps we'll be performing with [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/){:target="_blank"} and [DADA2](https://benjjneb.github.io/dada2/index.html){:target="_blank"}. Don't worry if anything seems unclear right now, we will discuss each at each step.

||Command|What we're doing|
|:--:|:--------:|----------|
|1|`bbduk.sh`/`filterAndTrim()`|removing primers and quality trimming/filtering|
|2|`learnErrors()`|generate an error model of our data|
|3|`derepFastq`|dereplicate sequences|
|4|`dada()`|infer ASVs|
|4|`mergePairs()`|merge forward and reverse reads to further refine ASVs|
|5|`makeSequenceTable()`|generate a count table|
|6|`removeBimeraDenovo()`|screening for and remove chimeras|
|7|`assignTaxonomy()`|assign taxonomy|

And at the end of this we'll do some R magic to generate regular [flat files](/bash/basics#whats-a-plain-text-file){:target="_blank"} for the standard desired outputs: a fasta file of our ASVs, a count table, and a taxonomy table.  

In our working directory there are 20 samples with forward (R1) and reverse (R2) reads with per-base-call quality information, so 40 fastq files (.fq). I typically like to have a file with all the sample names to use for various things throughout, so here's making that file based on how these sample names are formatted:

```bash
ls *_R1.fq | cut -f1 -d "_" > samples
```

If you're not comfortable with that line, consider running through the [bash basics](/bash/basics){:target="_blank"} and/or [six glorious commands](/bash/six_commands){:target="_blank"} pages before going further here 🙂  

## Removing primers and filtering by length
To start, we need to remove the primers from all of these (the primers used for this run are in the "primers.fa" file in our working directory), and here we're going to use [bbduk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/){:target="_blank"} to do that, but we'll need to run it on each sample individually. So we're going to use a wonderful little bash loop to do that.  

First, here's what the command would look like on an individual sample: 

```bash
bbduk.sh in=B1_sub_R1.fq in2=B1_sub_R2.fq \
out=B1_sub_R1_trimmed.fq.gz out2=B1_sub_R2_trimmed.fq.gz \
literal=GTGCCAGCMGCCGCGGTAA,GGACTACHVGGGTWTCTAAT k=10 ordered=t mink=2 \
ktrim=l rcomp=f minlength=220 maxlength=280 tbo tpe
```

There's a lot of info about bbduk at the [guide](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) for it, and it has an extensive help menu you can access by running it with no arguments, `bbduk.sh`, but here's what we've specified above. The in/out in2/out2 are for the forward and reverse reads. We've provided the forward and reverse primers to the `literal` flag, this will create all of the versions of these primers that are represented by the degenerate bases (e.g. M,V, etc.) and search for all possibilities. `k=10` is setting the kmer size we want to search for – the default is 23 for the version I was using ("Last modified April 12, 2018"), but our primers are 19 and 21 bps long, so we need to use a smaller word size. `ordered=t` just keeps the reads in the same order as we gave them to the program. `mink=2` specifies the smallest word size it will check against either edge of a read. `ktrim=l` tells the program to trim everything to the left of the identified primers. `rcomp=f` states *not* to look for the reverse complement – since the reads (300 bps) are longer than the targeted fragment (~290 bps) in this case, both primers are typically on both reads, so if we look for the reverse complements too we would be trimming to the left of them and dumping the whole read). Min and max length were based roughly on 10% smaller and bigger than would be expected after trimming the primers (so this sort of thing is based on your amplicon size and what you're looking to do). `tbo`, which trims primers based on overlap, and `tpe`, which trims forward and reverse reads to the same length, were helpful here because, again, in this case our reads are longer than our amplicons.

Back to our loop, here is how we can run it on all our samples at once, and since we have a lot of samples here, I'm [redirecting](http://localhost:4000/bash/basics#pipes-and-redirectors) the "stderr" (what's printing the stats for each sample) to a file so we can more easily view and keep track of if we're losing a ton of sequences by having that information stored somwhere – instead of just plastered to the terminal window. If you ran the above command on the single sample, first do `rm B1_sub_R*.gz`. 

```bash
for sample in $(cat samples); \
do echo "On sample: $sample"; bbduk.sh in="$sample"_sub_R1.fq in2="$sample"_sub_R2.fq \
out="$sample"_sub_R1_trimmed.fq.gz out2="$sample"_sub_R2_trimmed.fq.gz \
literal=GTGCCAGCMGCCGCGGTAA,GGACTACHVGGGTWTCTAAT k=10 ordered=t mink=2 \
ktrim=l rcomp=f minlength=220 maxlength=280 tbo tpe; \
done 2> bbduk_primer_trimming_stats.txt
```

I'm not going to break down the loop here as we have other fish to fry, but if this looks confusing to you, then check out the pages on bash [basics](/bash/basics){:target="_blank"} and [loops](/bash/loops_to_help){:target="_blank"}. While odd-looking at first, little command-line loops are extremely powerful, and trust me, you can learn to leverage that power more quickly than you'd think!  

Here's a little one-liner to look at what fraction of reads were thrown away due to missing or imperfect primers:

```bash
paste samples <(grep "Total Removed:" bbduk_primer_trimming_stats.txt | cut -f2)
```

<center><img src="{{ site.url }}/images/dada2_primer_trim.png"></center>
<br>
And this is showing us the fraction of reads lost in each sample, which doesn't really go over 4% in this case. With primers removed, we're now ready to jump into R and start using DADA2. 

# Processing with DADA2 in R
Just as with the bash component above, this portion assumes you have some baseline experience with R already. If you aren't familiar with R at all yet it's probably a good iea to run through the [R basics page](/R/basics){:target="_blank"} first. But even if you don't have any experience with R yet, you'll still be able to follow along here and run everything if you'd like. A full R script of everything done here is available in the "R_working_dir" subdirectory called "dada2_example.R" that can be opened in RStudio if you prefer to follow along with that rather than copying and pasting commands from here. Either way, this part really isn't about the R code (right now), it's more about the processing.  

## Setting up our working environment
To get started let's open up RStudio and take care of a few things. If you need to install the DADA2 package, visit the [installation page](https://benjjneb.github.io/dada2/dada-installation.html){:target="_blank"} provided by [@bejcal](https://twitter.com/bejcal){:target="_blank"}. This may require you needing to update your version of R and/or RStudio as well, which can be tricky sometimes. I've found [this page](https://www.linkedin.com/pulse/3-methods-update-r-rstudio-windows-mac-woratana-ngarmtrakulchol/){:target="_blank"} put together by [Woratana Ngarmtrakulchol](https://www.linkedin.com/in/woratana-ngarmtrakulchol-0079b766/){:target="_blank"} to be a lifesaver for me on more than one occasion 🙂

```R
library(dada2)
packageVersion("dada2") # 1.8.0

setwd("~/dada2_amplicon_ex_workflow")

list.files() # make sure what we think is here is actually here

  # setting a few variables we're going to use
samples <- scan("samples", what="character")

  # one holding the file names of all the forward reads, and one with the reverse
forward_reads <- sort(list.files(pattern="_R1_trimmed.fq.gz"))
reverse_reads <- sort(list.files(pattern="_R2_trimmed.fq.gz"))

  # and variables holding file names for the forward and reverse filtered reads we're going to generate with the next function:
filtered_forward_reads <- paste0(samples, "_sub_R1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_sub_R2_filtered.fq.gz")
```

## Quality trimming/filtering
We did a filtering step above with bbduk (where we eliminated reads that had imperfect or missing primers and those that were shorter than 220 bps or longer than 280), but in DADA2 we'll implement a trimming step as well (where we trim reads down based on some quality threshold rather than throwing the read away). Since we're potentially shortening reads further, we're also going to include another minimum-length filtering component. In DADA2, this is done with the `filterAndTrim()` function:  

```R
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, maxEE=c(2,2), rm.phix=TRUE, multithread=TRUE, minLen=175)
```

Here, the first and third arguments ("forward_reads" and "reverse_reads") are our input files, which are our primer-removed output fastq files from bbduk. The second and fourth are the variables holding the file names of the output forward and reverse seqs. And then we have a few parameters explicitly specified. `maxEE` is the quality filtering threshold being applied, and in this case we are saying we want to throw the read away if it is likely to have more than 2 erroneous base calls (we are specifying for both the forward and reverse reads separately). `rm.phix` removes any reads that match the PhiX bacteriophage genome, which is typically added to Illumina sequencing runs for quality monitoring. `multithread` will run the program in parallel if set to `TRUE` or if a number is given specifying how many cores you'd like run. And `minLen` is setting the minimum length reads we want to keep after trimming. The trimming occurring is coming from a default setting, `truncQ`, which is set to 2 unless we specify otherwise, meaning it trims all bases after the first quality score of 2 it comes across in a read. There is also an addition filtering default parameter that is removing any sequences containing any Ns, `maxN`.  

As mentioned, the output read files were named in those variables we made above ("filtered_forward_reads" and "filtered_reverse_reads"), so those were written to files already – which we can see if we run `list.files()` in R, or by checking in our your working directory in the terminal:

<center><img src="{{ site.url }}/images/dada2_filtered_head.png"></center>
<br>
But we also generated an object in R called filtered_out. And that's a matrix holding how many reads went in and how many reads made it out:

```R
class(filtered_out) # matrix
dim(filtered_out) # 20 2

filtered_out

#                            reads.in reads.out
# B1_sub_R1_trimmed.fq.gz        1637      1412
# B2_sub_R1_trimmed.fq.gz         589       472
# B3_sub_R1_trimmed.fq.gz         505       431
# B4_sub_R1_trimmed.fq.gz         508       447
# BW1_sub_R1_trimmed.fq.gz       2304      1988
# BW2_sub_R1_trimmed.fq.gz       6151      5320
# R10_sub_R1_trimmed.fq.gz      11792     10140
# R11_sub_R1_trimmed.fq.gz       9210      7925
# R11BF_sub_R1_trimmed.fq.gz     9273      8202
# R12_sub_R1_trimmed.fq.gz      16057     13836
# R1A_sub_R1_trimmed.fq.gz      12453     10497
# R1B_sub_R1_trimmed.fq.gz      16438     14066
# R2_sub_R1_trimmed.fq.gz       17670     14998
# R3_sub_R1_trimmed.fq.gz       17950     15210
# R4_sub_R1_trimmed.fq.gz       19100     16313
# R5_sub_R1_trimmed.fq.gz       18745     16137
# R6_sub_R1_trimmed.fq.gz       15183     12890
# R7_sub_R1_trimmed.fq.gz        8181      6936
# R8_sub_R1_trimmed.fq.gz       12622     10840
# R9_sub_R1_trimmed.fq.gz        8968      7662
```

DADA2 also holds a handy quality plotting function to let you visual how your reads are doing, `plotQualityProfile()`. By running that on our variable that holds all of our forward and reverse filtered reads, we can easily generate plots for all of them or a subset of them if wanted:

```R
plotQualityProfile(filtered_forward_reads)
plotQualityProfile(filtered_reverse_reads)
plotQualityProfile(filtered_reverse_reads[17:20])
```

<center><img src="{{ site.url }}/images/dada2_first_filter.png"></center>
<br>
The forward reads look pretty great, but the tail end of the reverse reads consistently drops in quality below 30 by about the 200th base and to a median of about 20–25 by the end (chemistry gets tired 😞). In [Phred](https://en.wikipedia.org/wiki/Phred_quality_score){:target="_blank"} talk the difference between a quality score of 40 and a quality score of 20 is an expected error rate of 1 in 10,000 to 1 in 100. So since we have full overlap with these primers and the sequencing performed (515f-806r, 2x300), we can be pretty conservative and trim these up a bit more. I'm going to cut the forward reads at 250 and the reverse reads at 200 in this case: 

```R
  # note that while we're making a new R object, we are providing the same 
  # output files for the reads given by our variables in here, so we are 
  # overwriting the actual output files
filtered_out2 <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, maxEE=c(2,2), rm.phix=TRUE, multithread=TRUE, minLen=175, truncLen=c(250,200))

plotQualityProfile(filtered_reverse_reads[17:20])
```  

<center><img src="{{ site.url }}/images/dada2_second_filter.png"></center>
<br>
Now we're lookin' good.

## Generate an error model of our data
Next

## Dereplication


## Inferring ASVs


## Merge forward and reverse reads


## Generate a count table


## Assigning taxonomy


# Extracting the standard goods from R


# 16S and 18S mixed together?
So I spent a decent amount of time with Josh trying to figure out how we could filter the dada2 mergePairs() objects to figure out which failed-to-merge sequences failed because they were likely not overlapping (and therefore anticipated to be 18S which we wanted to keep) vs those that just failed to merge because of sequencing error (things we wanted to still throw away). And you can get in quite the rabbit hole about this looking at nmismatches and nmatches and such, but ultimately that seemed to be something that would be very dataset dependent (I couldn't mentally tease out how much varying diversity would affect what were good cutoffs even for the same primers). So ultimately I bailed on trying to come up with a general way to do it that way, and decided to parse the reads before putting them into dada2. So here's how I went about this: 1) blast all reads to pr2 database (magicblast can take paired fastq files, which is pretty sweet, so I used that); 2) filter the blast output based on %ID and % of query sequence aligned (of both reads); 3) split the starting fastq files into 16S and 18S; 4) process both separately in dada2 and merge them at the end :) 




