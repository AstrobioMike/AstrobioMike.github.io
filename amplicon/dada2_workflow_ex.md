---
layout: main
title: DADA2 example workflow
categories: [amplicon, tutorial]
permalink: amplicon/dada2_workflow_ex
---

{% include _amplicon_ex_workflow_toc.html %}

{% include _side_tab_amplicon.html %}

[DADA2](https://benjjneb.github.io/dada2/index.html){:target="_blank"} is a relatively new processing workflow for recovering single-nucleotide resolved Amplicon Sequence Variants (ASVs) from amplicon data – if you're unfamiliar with ASVs, you can read more about ASVs vs OTUs in the [opening caveats](/amplicon/dada2_workflow_ex#opening-caveats) section that follows. Developed and maintained by [@bejcal](https://twitter.com/bejcal){:target="_blank"} et al., DADA2 leverages sequencing quality and abundance information to a greater extent than previously developed tools and generates an error model based on your actual data. It then uses this error model to do its best to infer the true biological sequences that generated your data. The original paper can be found [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/){:target="_blank"}, and the DADA2 R package site is [here](https://benjjneb.github.io/dada2/index.html){:target="_blank"}. The DADA2 team already has a great tutorial available [here](https://benjjneb.github.io/dada2/tutorial.html){:target="_blank"}, and I learned my way around DADA2 by following that and reading through the [manual](https://www.bioconductor.org/packages/3.3/bioc/manuals/dada2/man/dada2.pdf){:target="_blank"}. While the overall workflow presented here is the same, I've added 2 main things (in addition to all of my rambling of course): 1) example code for how to generate the standard output files of amplicon processing from the final R objects that DADA2 creates (i.e. a fasta of ASVs, a count table, and a taxonomy table); and 2) a [section](/amplicon/dada2_workflow_ex#16s-and-18s-mixed-together) for people working with 16S and 18S sequences mixed together. 
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
You can download the required dataset and files by copying and pasting the following commands into your terminal. For speed purposes we're only going to work with about 10% of the full dataset. Altogether the uncompressed size of the working directory is ~300MB.

```
cd ~
curl -L -o dada2_amplicon_ex_workflow.tar.gz https://ndownloader.figshare.com/files/11331515
tar -xzvf dada2_amplicon_ex_workflow.tar.gz
rm dada2_amplicon_ex_workflow.tar.gz
cd dada2_amplicon_ex_workflow/
```

Now, let's get started!
<br>
<br>

---
<br>

# Processing overview
It's good to try to keep a bird's-eye view of what's going on. So here is an overview of the main processing steps we'll be performing with [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/){:target="_blank"} and [DADA2](https://benjjneb.github.io/dada2/index.html){:target="_blank"}. Don't worry if anything seems unclear right now, we will discuss each at each step.

||Command|What we're doing|
|:--:|:--------:|----------|
|1|`bbduk.sh`/`filterAndTrim()`|remove primers and quality trim/filter|
|2|`learnErrors()`|generate an error model of our data|
|3|`derepFastq`|dereplicate sequences|
|4|`dada()`|infer ASVs on both forward and reverse reads independently|
|5|`mergePairs()`|merge forward and reverse reads to further refine ASVs|
|6|`makeSequenceTable()`|generate a count table|
|7|`removeBimeraDenovo()`|screen for and remove chimeras|
|8|`assignTaxonomy()`|assign taxonomy|

And at the end of this we'll do some R magic to generate regular [flat files](/bash/basics#whats-a-plain-text-file){:target="_blank"} for the standard desired outputs: a fasta file of our ASVs, a count table, and a taxonomy table.  

In our working directory there are 20 samples with forward (R1) and reverse (R2) reads with per-base-call quality information, so 40 fastq files (.fq). I typically like to have a file with all the sample names to use for various things throughout, so here's making that file based on how these sample names are formatted:

```bash
ls *_R1.fq | cut -f1 -d "_" > samples
```

If you're not comfortable with that line, consider running through the [bash basics](/bash/basics){:target="_blank"} and/or [six glorious commands](/bash/six_commands){:target="_blank"} pages before going further here 🙂  

# Removing primers
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

DADA2 also holds a handy quality plotting function to let you visualize how your reads are doing, `plotQualityProfile()`. By running that on our variable that holds all of our forward and reverse filtered reads, we can easily generate plots for all of them or a subset of them if wanted:

```R
plotQualityProfile(filtered_forward_reads)
plotQualityProfile(filtered_reverse_reads)
plotQualityProfile(filtered_reverse_reads[17:20])
```

The forward reads look pretty great, but the tail end of the reverse reads consistently drops in quality below 30 by about the 200th base and to a median of about 20–25 by the end (chemistry gets tired 😞). Here's the output of last 4 samples' reverse reads:

<center><img src="{{ site.url }}/images/dada2_first_filter.png"></center>

In [Phred](https://en.wikipedia.org/wiki/Phred_quality_score){:target="_blank"} talk the difference between a quality score of 40 and a quality score of 20 is an expected error rate of 1 in 10,000 to 1 in 100. In this case, since we have full overlap with these primers and the sequencing performed (515f-806r, 2x300), we can be pretty conservative and trim these up a bit more. I'm going to cut the forward reads at 250 and the reverse reads at 200 – roughly where both sets reach a median quality of 30 – and then see how things look: 

```R
filtered_out2 <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, maxEE=c(2,2), rm.phix=TRUE, multithread=TRUE, minLen=175, truncLen=c(250,200))
  # note that while we're making a new R object, we are providing the same output files 
  # for the reads given by our variables in here, so we are overwriting the actual output files

plotQualityProfile(filtered_reverse_reads[17:20])
```  

<center><img src="{{ site.url }}/images/dada2_second_filter.png"></center>

Now we're lookin' good.

## Generate an error model of our data
Next up is generating our error model by learning the specific error-signature of our dataset. I mentioned above what expected error rates correspond to a couple [Phred](https://en.wikipedia.org/wiki/Phred_quality_score){:target="_blank"} quality scores, but each sequencing run, even when all goes well, will have its own subtle variations to this. This step tries to assess that for both the forward and reverse reads.
It can be one of the more computationally intensive steps of the workflow, for this dataset on my laptop (2013 MacBook Pro) these each took about 5 minutes. 
```R
err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)
```

The developers have incorporated a plotting function to visualize how well the estimated error rates match up with the observed:

```R
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)
```

The forward and reverse didn't look too different, here's the output of the reverse:

<center><img src="{{ site.url }}/images/dada2_err_plot.png"></center>
<br>
[@bejcal](https://twitter.com/bejcal){:target="_blank"} goes into how to assess this a bit [here](https://benjjneb.github.io/dada2/tutorial.html#learn-the-error-rates){:target="_blank"}. The red line is what is expected based on the quality score, the black line represents the estimate, and the black dots represent the observed. This is one of those cases where this isn't a binary thing like "yes, things are good" or "no, they're not". I imagine over time and seeing outputs like this for multiple datasets you get a better feeling of what to expect and what should be more cause for alarm (as was the case for me with interpreting and making decisions based on quality-score plots like those [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/){:target="_blank"} produces). But generally speaking, you want the observed (black dots) to track well with the estimated (black line). [@bejcal](https://twitter.com/bejcal){:target="_blank"} notes [here](https://benjjneb.github.io/dada2/tutorial.html#learn-the-error-rates){:target="_blank"} that you can try to improve this by increasing the number of bases the function is using (default 100 million).

## Dereplication
Dereplication is a common step in many amplicon processing workflows. Instead of keeping 100 identical sequences and doing all downstream processing to all 100, you can keep/process one of them, and just attach the number 100 to it. When DADA2 dereplicates sequences, it also generates a new quality-score profile of each unique sequence based on the average quality scores of each base of all of the sequences that were replicates of it. 

```R
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples
```

## Inferring ASVs
Here's where DADA2 gets to do what it was born to do, that is to do its best to infer true biological sequences. It does this by incorporating the consensus quality profiles and abundances of each unique sequence, and then figuring out if it's more likely to be of biological origin or more likely to be spurious. You can read more about the details of this in [the paper](https://www.nature.com/articles/nmeth.3869#methods){:target="_blank"} of course or looking through [the site](https://benjjneb.github.io/dada2/index.html){:target="_blank"}. This step can be run on individual samples, which is the least computationally intensive manner, or on all samples together, which increases the function's ability to resolve low-abundance ASVs. Imagine Sample A has 10,000 copies of sequence Z, and Sample B has 1 copy of sequence Z. Sequence Z would likely be filtered out of Sample B even though it was a "true" singleton among perhaps thousands of spurious singletons we needed to remove. Because running all samples together on large datasets can become impractical very quickly, the developers also added a way to try to combine the best of both worlds they refer to as pseudo-pooling, which is demonstrated very nicely [here](https://benjjneb.github.io/dada2/pseudo.html#Pseudo-pooling){:target="_blank"}. This basically provides a way to tell Sample B from the above example that sequence Z is legit. But it's noted at the end of the [pseudo-pooling page](https://benjjneb.github.io/dada2/pseudo.html#Pseudo-pooling){:target="_blank"} that this is not always the best way to go, and it may depend on your experimental design which is likely more appropriate for your data – as usual. There are no one-size-fits-all solutions in bioinformatics! But that's exactly what makes it so damn fun 🙂

```R
dada_forward <- dada(derep_forward, err=err_forward_reads, multithread=TRUE, pool="pseudo")
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, multithread=TRUE, pool="pseudo")
```

## Merge forward and reverse reads
Now DADA2 merges the forward and reverse ASVs to reconstruct our full target amplicon requiring the overlapping region to be identical between the two. By default it requires that at least 12 bps overlap, but in our case the overlap should be much greater. If you remember above we trimmed the forward reads to 250 and the reverse to 200, and our primers were 515f–806r. After cutting off the primers we're expecting a typical amplicon size of around 260 bases, so our typicaly overlap should be up around 190. That's estimated based on *E. coli* 16S rRNA gene positions and very back-of-the-envelope-esque of course, so to allow for true biological variation and such I'm going ot set the minimum overlap for this dataset for 170. I'm also setting the trimOverhang option to `TRUE` in case any of our reads go passed their opposite primers (which I wouldn't expect based on our trimming, but is possible due to the region and sequencing method).

```R
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse, derep_reverse, trimOverhang=TRUE, minOverlap=170)

  # this object holds a lot of information that may be the first place you'd want to look if you want to start poking under the hood
class(merged_amplicons) # list
length(merged_amplicons) # 20, one for each of our samples
names(merged_amplicons) # the names() function gives us the name of each element of the list 

class(merged_amplicons$B1) # each element of the list is a dataframe that can be accessed and manipulated like any ordinary dataframe

names(merged_amplicons$B1) # the names() function on a dataframe gives you the column names
# "sequence"  "abundance" "forward"   "reverse"   "nmatch"    "nmismatch" "nindel"    "prefer"    "accept"
```

## Generate a count table
Now we can generate a count table with the `makeSequenceTable()` function:

```R
seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab) # matrix
dim(seqtab) # 20 2567
```

We can see from the dimensions of the "seqtab" matrix that we have 2,567 ASVs in this case. But it's not very friendly to look at in this form because the actual sequences are our rownames. But we'll make a more traditional count table in a couple steps.

## Chimera identification
DADA2 identifies likely chimeras by aligning each sequence with those that were recovered in greater abundance and then seeing if there are any sequences that can be made exactly by mixing left and right portions of two of the more-abundant ones. These are then removed:

```R
seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=T, verbose=T) # Identified 21 bimeras out of 2567 input sequences.

  # though we only lost 21 sequences, we don't know if they held a lot in terms of abundance, this is one quick way to look at that
sum(seqtab.nochim)/sum(seqtab) # 0.9925714
```

## Overview of counts throughout
The developers' [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html){:target="_blank"} provides an example of nice, quick way to pull out how many reads were dropped at various points of the pipeline. This can serve as a jumping off point if you're left with too few to help point you towards where you should start digging. 

```R
  # set a little function
getN <- function(x) sum(getUniques(x))

  # making a little table
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1], filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN), dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN), nonchim=rowSums(seqtab.nochim), total_perc_reads_lost=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

summary_tab

#       dada2_input filtered dada_f dada_r merged nonchim total_perc_reads_lost
# B1           1637     1412   1491   1498   1460    1460                  89.2
# B2            589      472    535    534    533     533                  90.5
# B3            505      431    465    466    465     465                  92.1
# B4            508      447    453    457    446     446                  87.8
# BW1          2304     1989   2109   2122   2076    2076                  90.1
# BW2          6151     5320   5335   5433   4859    4859                  79.0
# R10         11792    10140  10266  10465   9574    9403                  79.7
# R11BF        9210     7928   7648   7879   7120    6981                  75.8
# R11          9273     8204   8159   8276   7743    7506                  80.9
# R12         16057    13839  12921  13449  11159   11091                  69.1
# R1A         12453    10503  10110  10419   9041    9017                  72.4
# R1B         16438    14069  13513  13964  11699   11653                  70.9
# R2          17670    15001  14715  15177  13054   12995                  73.5
# R3          17950    15214  14864  15333  13145   13082                  72.9
# R4          19100    16321  16703  16940  15278   15212                  79.6
# R5          18745    16146  15502  16080  13544   13455                  71.8
# R6          15183    12891  12618  12973  11034   11023                  72.6
# R7           8181     6939   6782   6982   5931    5919                  72.4
# R8          12622    10843  10882  11062  10174   10051                  79.6
# R9           8968     7664   7649   7825   7146    7099                  79.2
```

## Assigning taxonomy
DADA2 incorporates a function that assigns taxonomy using the [RDP's kmer-based method](https://rdp.cme.msu.edu/classifier/classifier.jsp){:target="_blank"}, original paper [here](http://www.ncbi.nlm.nih.gov/pubmed/17586664){:target="_blank"}. There are some DADA2-formatted databases available [here](https://benjjneb.github.io/dada2/training.html){:target="_blank"}, which is where the silva one came from that is in our working directory called here, but you can use whatever database you'd like following the formatting specified at the bottom of that page. This step took maybe 10 minutes on my laptop. 

```R
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=T, tryRC=T)
```

# Extracting the standard goods from R
The typical standard outputs from amplicon processing are a fasta file, a count table, and a taxonomy table. So here's how you can generate those files from your DADA2 objects in R:

```R
  # giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

  # making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

  # count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.txt", sep="\t", quote=F)

  # tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.txt", sep="\t", quote=F)
```

And now if we look back at our terminal, we can see the fruits of our labor are no longer confined to the R universe:

<center><img src="{{ site.url }}/images/dada2_outfiles.png"></center>
<br>
You can find examples and corresponding code of some of the things you might want to do with these files in the [analysis section](/amplicon/workflow_ex#analysis-in-r){:target="_blank"} of the [usearch/vsearch tutorial](/amplicon/workflow_ex){:target="_blank"}.

# 16S and 18S mixed together?
So if you're using primers [like these that capture 16S and 18S pretty well](https://onlinelibrary.wiley.com/doi/abs/10.1111/1462-2920.13023){:target="_blank"} and then writing cool papers [like this one](https://www.nature.com/articles/nmicrobiol20165){:target="_blank"}, then you might be wondering if/how you can do this while taking advantage of all the awesomeness that is DADA2. The 18S amplified fragments are typically too long to have any overlap occur with the 2x250bp sequencing that is often performed with these "V4V5" primers, so one method that has been employed was to just work with those that didn't successfully merge. So starting from within DADA2, I considered playing with the `mergePairs()` step and messing with the "justConcatenate" flag. The problem is that you also need to account for those that fail the merge just because of poor quality. So my good buddy [Josh](https://twitter.com/KlingJoshua){:target="_blank"} and I spent a few hours trying to figure out how we could filter the merged objects to figure out which failed-to-merge sequences failed because they were likely not overlapping (and therefore anticipated to be 18S which we wanted to keep) vs those that just failed to merge because of sequencing error (things we wanted to still throw away). And you can end up in quite the rabbit hole about this looking at nmismatches and nmatches and such, but ultimately it seemed like even if we found optimal parameters for one dataset, it would likely be very different for the next one. 

So I decided to bail on that approach and wanted to see if there was a reasonable way to parse the reads into 16S and 18S before putting them into DADA2. Here's what I came up with that so far has worked well:

||Command|What we're doing|
|:--:|:--------:|----------|
|1|`magicblast`|blast all reads to the [PR2 database](https://figshare.com/articles/PR2_rRNA_gene_database/3803709){:target="_blank"}, magicblast is built for this 🙂|
|2|`Filtering magicblast output`|based on % ID and % of query sequence aligned (of both reads)|
|3|`Splitting 16S/18S reads`|based on the magicblast filtering|
|4|`Processing both in DADA2`|processing independently and merging at the end|

## 16S/18S example data
For an example of this process, we're going to work with a couple of samples from [the paper](https://www.nature.com/articles/nmicrobiol20165){:target="_blank"} I mentioned above by [David Needham](https://twitter.com/animalkewls){:target="_blank"} and [Jed Fuhrman](https://twitter.com/JedFuhrman){:target="_blank"}. If you'd like to follow along, you can download a small directory with the code and data with these commands:

<br>
<br>
<center><img src="{{ site.url }}/images/under_construction.jpeg"></center>
<center><h3>UNDER CONSTRUCTION</h3></center>
<br>
<br>


```bash
for sample in $(ls *_1.fastq.gz | cut -f1 -d "_"); do echo $sample; done > samples

for sample in $(cat samples); do bbduk.sh in="$sample"_1.fastq.gz in2="$sample"_2.fastq.gz out="$sample"_1_trimmed.fq.gz out2="$sample"_2_trimmed.fq.gz literal=GTGCCAGCMGCCGCGGTAA,CCGYCAATTYMTTTRAGTTT k=10 ordered=t mink=2 ktrim=l rcomp=t minlength=220; done

makeblastdb -in pr2_seqs_with_tax_headers.fa -dbtype nucl -parse_seqids -out pr2_magicblast_db

for sample in $(cat samples); do echo "$sample"; magicblast -db pr2_magicblast_db -query "$sample"_1_trimmed.fq.gz -query_mate "$sample"_2_trimmed.fq.gz -infmt fastq -out "$sample"_mblast_out.txt -outfmt tabular -num_threads 2 -splice F -no_unaligned; done

for sample in $(cat samples); do echo $sample; cut -f1,2,3,7,8,16 "$sample"_mblast_out.txt | sed '1d' | sed '1d' | sed 's/# Fields: //' | tr " " "_" | awk -v OFS='\t' 'NR==1 {$7="%_query_aln"; print $0} NR>1 { print $0, ($5-$4)/$6*100 }' > "$sample"_mblast_out_mod.txt; done
```

## Magicblast

## Filtering magicblast output

## Splitting starting fastq files into 16S/18S

## Processing both in DADA2 and merging at the end


