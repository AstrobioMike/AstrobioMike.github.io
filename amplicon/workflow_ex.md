---
layout: main
title: Example marker-gene workflow
categories: [amplicon, tutorial]
permalink: amplicon/workflow_ex
---

{% include _amplicon_ex_workflow_toc.html %}

{% include _side_tab_amplicon.html %}

<div class="warning">
<h2>WARNING!</h2>
Don't get too lost in the weeds when working with marker-gene data. More often than not the types of detailed questions that will land you in that situation can't be answered with this type of data anyway!</div>

## Opening caveats
Most often a marker-gene analysis is the microbial ecologist's first tool in a vast toolkit. It is primarily used as a broad survey of community structure. As the warning notes above, it is easy to get caught spinning your wheels about a sublte component in your processing pipeline that ultimately has a negligible impact compared to the noise we are working through. What I mean by this is, generally speaking, tag data is not the appropriate tool to answer really meticulous questions. It is a tool for comparing baseline *proxies* of metrics about microbial communities. It is a tool of exploration and hypothesis generation, not hypothesis confirmation.  

There are a lot of things to keep in mind regarding what tag sequencing means, what it doesn't mean, what it can tell you, and what it can't. And at some point I'll give my two cents on all of these things in the [Amplicon Thoughts]({{ site.url}}/amplicon/thoughts) post should anyone be interested. But for now, let's go through a dataset.

There are many ways to process amplicon data (if you need a quick primer on some relevant terminology, visit the [Amplicon main page]({{ site.url }}/amplicon/)). Some of the most widely used tools/pipelines include [mothur](https://www.mothur.org/), [usearch](https://drive5.com/usearch/), [Minimum Entropy Decomposition](http://merenlab.org/2014/11/04/med/), [DADA2](https://benjjneb.github.io/dada2/index.html), and [qiime2](https://qiime2.org/) (which employs other tools within it). As usual, there is no one-size-fits-all, there is no 'best'. And actually in my experience if you make similar decisions when processing your sequences (decisions about things like minimum abundance filtering or clustering thresholds), you most often get very similar results regardless of which of these you use (refreshingly). Here I'll be using [usearch](https://drive5.com/usearch/) mostly because of it's ease of deployment; there is no installation required and there are no dependencies. So if you want to actively follow along with this module you can just download it and you're ready to rock. But please keep in mind that nothing here is meant to be authoritative. This is simply one example of one workflow. When working with your own data you should never follow any pipeline blindly.
<br>
<br>

---
<br>

# Tools used here

This module is going to walkthrough one possible workflow for an amplicon dataset from start to finish. This will entail processing the raw sequences with [usearch](https://drive5.com/usearch/), and analyzing the output and making some visualizations with [R](https://www.r-project.org/) using some great packages like [*vegan*](https://github.com/vegandevs/vegan), [*phyloseq*](http://joey711.github.io/phyloseq/), and [*breakaway*](https://github.com/adw96/breakaway). 

There is a free, 'light-weight' version of usearch available that you will need to download if you wish to follow along here. To get the free version, go to [https://www.drive5.com/usearch/download.html](https://www.drive5.com/usearch/download.html), and fill out a (very) short form in order to have a download link sent to you. This usually happens virtually instantly in my experience. At the time this page is being put together I'm using usearch v10.0.240 on Mac OSX. As mentioned, we're using usearch here because there is no installation required. So once you receive the email and download the file, you only need to put it somewhere in your [PATH](http://localhost:4000/bash/modifying_your_path) and make sure it is 'executable'. Assuming you're working on a Mac and downloaded the same version as noted above, and it was downloaded into your home downloads directory, this should get the job done:


```bash
mv ~/Downloads/usearch10.0.240_i86osx32 /usr/local/bin/usearch
chmod +x /usr/local/bin/usearch
```

The `mv` command above moves the file from the downloads directory into your local bin (which is in your PATH so you'll be able to call usearch from anywhere), and renames it to 'usearch' rather than the longer name with the version info. And the `chmod` command makes it executable, which just means you can run it as a program from the command line. After this is done you should be able to call usearch from any terminal window, and you can test this with `usearch --version`: 

<center><img src="{{ site.url }}/images/usearch_version.png"></center>

<br>
If you wish to follow along with the analysis portion in [R](https://www.r-project.org/), then you also need to have a working installation of R on your computer. If you're unsure if you already have R, you can check by typing `R` in a terminal window. If this launches R rather than giving an error message, you should be good to go (enter `q()` to exit the R environment). If you do not have R you can download it from here for Mac: [https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/). And if you have a relatively newer Mac you will also need to install XQuartz which you can get from here: [https://www.xquartz.org/](https://www.xquartz.org/). 

Lastly, I highly, highly, highly recommend installing RStudio if you don't already have it and use it. RStudio is an interface for R that not only makes everything you will do in R easier and more organized, but it's also invaluable for reproducibility of your analyses as it makes it second-nature to generate and save R scripts of everything you've done. You can download an RStudio installer from [here](https://www.rstudio.com/products/rstudio/download/#download).
<br>
<br>

---
<br>

# The data

We are going to work with a subset of the dataset published [here](https://www.frontiersin.org/articles/10.3389/fmicb.2015.01470/full). We were exploring an underwater mountain ~3 km down at the bottom of the Pacific Ocean that serves as a low-temperature (~5-10°C) hydrothermal venting site. This amplicon dataset was generated from DNA extracted from crushed basalts collected from across the mountain in order to begin characterizing the microbial communities of these deep-sea rocks. It was generated via Illumina MiSeq 2x300 paired-end sequencing using primers targeting the V4 region of the 16S rRNA gene. There are 20 samples total: 4 extraction blanks (nothing added to DNA extraction kit), 2 bottom water samples, 13 rocks, and one biofilm scraped off of a rock. 

<center><img src="{{ site.url }}/images/dorado.png"></center>

<br>
You can download the required dataset and files by copying and pasting the following commands into your terminal. For speed purposes we're only going to work with about 10% of the full dataset. Altogether the uncompressed size of the working directory is < 300MB.

```
cd ~
curl -O https://AstrobioMike.github.io/tutorial_files/amplicon_example_workflow.tar.gz
tar -xvf amplicon_example_workflow.tar.gz
rm amplicon_example_workflow.tar.gz
cd amplicon_example_workflow/
```

Now, let's get started!
<br>
<br>

---
<br>

# Processing

It's good to try to keep a bird's-eye view of what's going on. So here is an overview of the main processing steps we'll be performing with usearch. Don't worry if anything seems unclear right now, we will discuss each at each step.

||Command|What we're doing|
|:--:|:--------:|----------|
|1|`-fastq_mergepairs`|merge forward and reverse reads together|
|2|`-fastx_truncate`|cut off forward and reverse primers|
|3|`-fastq_filter`|quality filter sequences|
|4|`-fastx_uniques`|dereplicate sequences|
|5|`-cluster_otus`/`-unoise3`|cluster sequences into OTUs and/or generate ASVs|
|6|`-otutab`|generate a count table|
|7|`-sintax`|assign taxonomy to OTUs and/or ASVs|

<h4>1) Merging forward and reverse reads</h4>

Depending on the sequencing facility, when paired-end sequencing was performed you may receive one forward and one reverse fastq file with all samples mixed together, or you may get things already demultiplexed (separate files for each individual sample). In this case our samples have already been demultiplexed and the barcodes have been trimmed off. This means for each of our 20 samples we have a forward (R1) and reverse (R2) reads file:

<center><img src="{{ site.url }}/images/ls_all_sample_reads.png"></center>

<br>
The format of the read file names is a letter or 2 specifying which sample it is, followed by a 'sub' (just because these are subsampled from the full file sizes), and then an R1 or R2. All of the .fq files are fastq files (sequences with quality score information), and the .fa files are simply sequences with no quality information. For some basic information about these formats visit the [Amplicon main page]({{ site.url }}/amplicon/).

So our first step is to merge these forward and reverse reads for each sample. This particular sequencing run was 2x300 paired-end over the V4 region, which is roughly ~300 base pairs (bps). So we have virtually full overlap for all of our sequences. The benefit of full-overlap as compared to sequencing a longer region with less overlap is improved base-calling accuracy. The usearch command `-fastq_mergepairs` as run below does a few things: 1) merges the paired reads for each sample into a single sequence; 2) appends a sample name to the merged sequence header to identify which sample it originated from; 3) combines and merged sequences from all samples into one file; 4) calculates new quality scores for each position of the merged sequences (as detailed [here](https://www.drive5.com/usearch/manual/merge_pair.html) if interested). We could list all of the files for each sample, but instead we're going to use a nice little *bash* shortcut called a wildcard. If you're unfamiliar with these, shoot over to the [Wildcards section](http://localhost:4000/bash/basics/#wildcards) of [*bash* basics](/bash/basics) when you have some time to learn about them. For now, you just need to know the `*` represents anything any number of times, so as we use it here it will grab all of the files we need, one at a time.

```
usearch -fastq_mergepairs *_R1.fq -fastqout all_samples_merged.fq -relabel @
```

So here the magic of the `*` wild card is grabbing all forward read files for every file that ends in .fq. We didn't need to specify the reverse reads because usearch finds them automatically for us. We specify the output file with the `-fastqout` argument. And the last part, `-relabel @`, tells the program to append the file name the merged sequence originated from to its header. Along with the command finishing we get some overall summary statistics printed to the terminal, such as how many total sequences we had and what proportion of the paired reads successfully merged. Checks like this are needed in order to have some idea of what's going on with your sequences as you move them along. Here almost 90% of them successfully merged, which is pretty good. 

<h4>2) Cutting off the primers</h4>
Due to the degeneracy of primers they tend to introduce non-biological variation. As such it's important to cut them off, particularly if you are going to perform any sort of single-nucleotide resolution clustering. Here we have our primers in a fasta file that looks like this:

<center><img src="{{ site.url }}/images/primers.png"></center>

<br>
First we want to make sure the primers are where they are supposed to be (that is at the front and end of our sequences). The way this is done in usearch is with the command `-search_oligodb`. It is rather time-consuming, however, so as long as our samples were all from the same sequencing run, we can randomly subset a portion of them to check:

```
usearch -fastx_subsample all_samples_merged.fq -sample_size 5000 -fastqout all_sub_for_primer_check.fq
usearch -search_oligodb all_sub_for_primer_check.fq -db primers.fa -strand both -userout primer_hits.txt -userfields query+qlo+qhi+qstrand
```

Here the first command is randomly subsampling 5,000 sequences from the total and creating a file we name "all_sub_for_primer_check.fq". And then the second command is searching for the primers in that subset of sequences, and generating an output file called "primer_hits.txt" with the results. This file has 4 columns per line: the sequence ID; the start position of a hit by a primer; the final position of a hit by a primer; and then a plus or minus for whether the primer was in the same orientation or if its reverse complement was the match. Looking through that file we consistently see that there is a hit to the first 19 bps of a sequence and the last 20. Here's just a shot of the `head` and `tail`:

<center><img src="{{ site.url }}/images/primer_hits.png"></center>

<br>
This tells us things are as we expect them to be, and we can proceed to trim off the first 19 and the last 20 bps of each sequence to remove our primers.

```
usearch -fastx_truncate all_samples_merged.fq -stripleft 19 -stripright 20 -fastqout all_merged_no_primers.fq
```

<h4>3) Quality filtering</h4>
Quality filtering is a critical in reducing the abundance and impact of spurious sequences. There is an intrinsic error rate to all sequencing technology (and polymerases) that will consistently be generating some portion of sequences that vary slightly from their true biological origin, and this can substantially inflate metrics such as richness and diversity. Quality filtering is one of the steps in place to mitigate that problem. In usearch it is done with the `-fastq_filter` command, which uses calculated expected errors to determine if sequences are filtered out or not (the details can be found [here](https://www.drive5.com/usearch/manual/exp_errs.html)).

```
usearch -fastq_filter all_merged_no_primers.fq -fastq_maxee 1 -fastaout QCd_merged.fa
```

The output tells us very few reads were discarded due to poor quality. And note that our output file from this is now a fasta file, as after quality filtering we no longer need to keep track of the quality scores for each sequence anymore. 

<h4>4) Dereplication</h4>
The dereplication step collapses all identical sequences to one and simply keeps track of how many there were. This serves a couple of purposes. It can save a ton of time in further processing steps because, for instance, you can then just process once a sequence that may appear 10,000 times – rather than processing 10,000 copies of the same exact sequence 10,000 times. And also, further attempts to mitigate the presence of spurious sequences involves merging low-abundance sequences with high-abundance sequences if they are similar enough, which requires knowing the abundance of individual sequences. In usearch this dereplication step is done with the `-fastx_uniques` command:

```
usearch -fastx_uniques QCd_merged.fa -sizeout -relabel Uniq -fastaout unique_seqs.fa
```

We can see from this output that of ~188,000 total sequences, ~95,000 of them were unique, and ~84,000 of those appear only one time (singletons). Having relatively so many singletons is common with amplicon data and is a direct consequence of intrinsic sequencing error rates. It is likely the majority of these are simply 1 or 2 bps diverged from a true biological sequence that is present in greater abundance. This is why we try to merge low-abundance sequences with high-abundance ones that are very similar – this is also one of the reasons why clustering sequences into arbitrary OTUs took off as well as it did, to mitigate sequencing error inflating true biological signals). 

<h4>5) Clustering OTUs and/or generating ASVs</h4>
OTUs refer to 'operational taxonomic units' that are made by grouping sequences into clusters based on similarity. Then for each cluster, a sequence from that cluster is selected to be the representative sequence of that cluster (typically one from the 'center' of that cluster in 'sequence space', or the most abundant from the cluster). In usearch, if you specify 97% OTUs for instance, the program attempts to build OTUs (clusters) such that none of the representative sequences for the OTUs are more similar than 97% to any other.

ASVs on the otherhand refer to 'amplicon sequence variants', which seems to be emerging as the consensus name for sequences derived using single-nucleotide resolution. As the name suggests, single-nucleotide resolution means you can delineate sequences even if they vary by only one base pair. This most often means similarity levels greater than 99% that as such would be lost with any form of OTU clustering. ASVs are believed to represent true biological sequences (with sequencing error for the most part weeded out via abundance-based merging and filtering), and have the advantage of not throwing away resolution by merging similar sequences together. There are several tools available that apply different approaches to do this such as [MED](http://merenlab.org/2014/11/04/med/), [DADA2](https://benjjneb.github.io/dada2/index.html), and usearch's `-unoise3` command. 

For now we will just move forward with ASVs, but feel free to experiment. When I process a new tag dataset I usual run both analyses alongside each other. There can be instances where greater resolution isn't helpful. As usual it depends on the dataset and on what you are trying to find. (Though I'll note there are some people that are pretty adamant that ASVs should be the only way to go moving forward.)

So, to generate our ASVs:

```
usearch -unoise3 unique_seqs.fa -zotus ASVs.fa -tabbedout unoise3.txt
```

By default, usearch's `-unoise3` command filters out any that are less abundant than 8, and removes sequences it suspects are chimeric, details of which can be found [here](https://www.biorxiv.org/content/biorxiv/early/2016/09/09/074252.full.pdf). From the output to the terminal we can see we have ~1700 ASVs. In theory, these are true biological sequences recovered from our 20 samples. 

<h4>6) Generating a count table</h4>
The count table is what tells us how many times each sequence appears in each sample. It is the end-product of all this processing that we can't wait to get into R. The way the count table is generated in usearch is a bit different than most other approaches. Here, now that usearch has identified what it believes to be true biological sequences, the `-otutab` command is used to attempt to map all of our merged sequences to the ASVs. By default it does this by requiring a merged sequence to be >97% similar to an ASV, and if a sequence crosses that threshold for more than one ASV it is counted for the one it is most similar to only. Incorporating another 97% threshold here may seem a bit confusing at first, but usearch's thinking is that more often than not, the 'true' ASV sequences recovered were more than likely the source of the majority of sequences that are within 3% of them, as the source of this 3% variation is believed to be the result of sequencing and pcr errors. Additionally, it is still fundamentally different to have your ASV units built without similarity clustering, as opposed to 97% OTUs. That said, when I ran through this I decided to bump up the required % similarity to 98.5% when mapping our merged sequences to our ASVs to generate our count table. The downside to this is you will lose some data (as less sequences will map to an ASV), but if you have a decent amount of sequences per sample then you can afford to lose some for the sake of being a bit more conservative. 

Due to a minor glitch in usearch where it doesn't seem to register the headers of the ASV sequences it gave them, we need to run this quick `sed` command before the `-otutab` program:

```
sed -i.tmp 's/Zotu/ASV_/' ASVs.fa
usearch -otutab all_merged_no_primers.fq -zotus ASVs.fa -otutabout ASV_counts.txt -id 0.985
```
Note that the input file is our previous merged sequences file after we trimmed off the primers. The output to the screen after it finishes shows us we mapped ~70% of our total sequences. And the "ASV_counts.txt" file is our glorious count table. And if we take a peek at it with `less -S ASV_counts.txt` we can see it is organized with samples as column names, ASV sequence IDs as row names, and their counts per sample in the cells:

<center><img src="{{ site.url }}/images/less_ASV_counts.png"></center>

<br>

<h4>7) Assigning taxonomy to our ASV sequences</h4>
The final step in our sequence processing is to assign taxonomy to our ASV sequences. There are multiple ways to do this as well. I usually run a few different ways and put them all in a table next to each other. Then I can be more somewhat more confident when there is some consistency. Something to keep in mind is taxonomy-assigning software is generally built to be fast. To do this, most rely on kmer frequencies to classify to a reference, rather than alignments. If you begin analyzing your data and some particular OTUs or ASVs emerge as being imporant to the overall story it's a good idea to take that sequence and BLAST it to try to get more information about it. 

Here we will use usearch's `sintax` program with the RDP training set reference fasta available from [here](https://www.drive5.com/usearch/manual/sintax_downloads.html).

```
usearch -sintax ASVs.fa -db rdp_16s_v16_sp.fa -tabbedout ASV_tax_raw.txt -strand both -sintax_cutoff 0.5
```

The ouput taxonomy file, "ASV_tax.txt", contains 4 columns: the ASV ID; taxonomic lineage with confidence scores for each tier; the orientation of the match to the reference; and then an abridged taxonomic lineage including only tiers with confidence scores above the 0.5 cutoff we specified in the command. Information on how the confidence scores are calculated can be found [here](https://www.drive5.com/usearch/manual/tax_conf.html).

And that's that. Now we're just going to run a few housekeeping commands to get things into a more R-friendly format and copy the relevant files into the R subdirectory:

```
sed -i.tmp 's/#OTU ID//' ASV_counts.txt
cp ASV_counts.txt R_working_dir/
bash convert_usearch_tax.sh
cp ASV_tax.txt R_working_dir/
```

And now we're ready to move onto analysis!
<br>
<br>

---
<br>

# Analysis in R
To get started let's open up RStudio, set our current working directory to where we just left our new files, and install the packages we'll need:

```R
setwd("~/amplicon_example_workflow/R_working_dir/")

install.packages("phyloseq")
install.packages("vegan")
install.packages("breakaway")
```

If any of those didn't succeed, for whichever caused a problem run the corresponding command below to try installing it a different way:

```R
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")

install.packages("devtools")
devtools::install_github("vegandevs/vegan")

install.packages("devtools")
devtools::install_github("adw96/breakaway")
```

After all 3 have installed successfully, load them up:

```R
library("phyloseq")
library("vegan")
library("breakaway")
```

And now it's time to read in our data.
 
<br>
<br>

<center><img src="{{ site.url }}/images/under_construction.jpeg"></center>
<center><h3>UNDER CONSTRUCTION</h3></center>
