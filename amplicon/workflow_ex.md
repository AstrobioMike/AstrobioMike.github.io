---
layout: main
title: Example marker-gene workflow
categories: [amplicon, tutorial]
permalink: amplicon/workflow_ex
---

{% include _amplicon_ex_workflow_toc.html %}

{% include _side_tab_amplicon.html %}

This module is going to walkthrough one possible workflow for an amplicon dataset from start to finish (if you need a quick primer on some relevant terminology, visit the [amplicon main page]({{ site.url }}/amplicon/)). This will entail processing the raw sequences with [usearch](https://drive5.com/usearch/), and analyzing the output and making some visualizations with [R](https://www.r-project.org/) using some great packages like [*vegan*](https://github.com/vegandevs/vegan) and [*phyloseq*](http://joey711.github.io/phyloseq/).
<br>
<br>

---
<br>

## Opening caveats
There are many ways to process amplicon data. Some of the most widely used tools/pipelines include [mothur](https://www.mothur.org/), [usearch](https://drive5.com/usearch/), [Minimum Entropy Decomposition](http://merenlab.org/2014/11/04/med/), [DADA2](https://benjjneb.github.io/dada2/index.html), and [qiime2](https://qiime2.org/) (which employs other tools within it). As usual, there is no one-size-fits-all; there is no best. And actually in my experience if you make similar decisions when processing your sequences (decisions about things like minimum abundance filtering or clustering thresholds), you will get very similar results regardless of which of these you use (refreshingly). Here I'll be using [usearch](https://drive5.com/usearch/) mostly because of it's ease of deployment; there is no installation required and there are no dependencies. So if you want to actively follow along with this module you can just download it and you're ready to rock. Keep in mind that nothing here is meant to be authoritative. This is simply one example of one workflow. When working with your own data you should of course never follow any pipeline blindly.
<br>
<br>

---
<br>

# Tools used here
There is a free, 'light-weight' version of usearch available that you will need to download if you wish to follow along here. To get the free version, go to [https://www.drive5.com/usearch/download.html](https://www.drive5.com/usearch/download.html), and fill out a (very) short form in order to have a download link sent to you. This usually happens virtually instantly in my experience. At the time this page is being put together I'm using usearch v10.0.240 on Mac OSX. As mentioned, we're using usearch here because there is no installation required. So once you receive the email and download the file, you only need to put it somewhere in your [PATH](http://localhost:4000/bash/modifying_your_path) and make sure it is 'executable'. Assuming you're working on a Mac and you downloaded the same version as noted above into your default download directory, this should get the job done (the sudo command will require you to enter your login password for your computer):


```bash
sudo mv ~/Downloads/usearch10.0.240_i86osx32 /usr/local/bin/usearch
chmod +x /usr/local/bin/usearch
```

The `mv` command above moves the file from the downloads directory into your local bin (which is in your PATH so you'll be able to call usearch from anywhere), and renames it to 'usearch' rather than the longer name with the version info. The `sudo` part of that just tells your computer that you have the authority to modify things in that `/usr/local/bin/` directory, which is why you need to enter your password. If you were working on a computer where you didn't have authority to do this, we'd need to put the usearch program somewhere else and add it to our [PATH]({{ site.url }}/bash/modifying_your_path). And the `chmod` command makes it executable, which just means you can run it as a program from the command line. After this is done you should be able to call usearch from any terminal window, and you can test this with `usearch --version`: 

<center><img src="{{ site.url }}/images/usearch_version.png"></center>

<br>
If you wish to follow along with the analysis portion in [R](https://www.r-project.org/), then you also need to have a working installation of R on your computer. If you're unsure if you already have R, you can check by typing `R` in a terminal window. If this launches R rather than giving an error message, you should be good to go (enter `q()` to exit the R environment). If you do not have R you can download it from here for Mac: [https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/). And if you have a relatively newer Mac you will also need to install XQuartz which you can get from here: [https://www.xquartz.org/](https://www.xquartz.org/). 

Lastly, I highly, highly, highly recommend installing RStudio if you don't already have it and use it. RStudio is an interface for R that not only makes everything you will do in R easier and more organized, but it's also invaluable for reproducibility of your analyses as it makes it second-nature to generate and save R scripts of everything you've done. You can download an RStudio installer from [here](https://www.rstudio.com/products/rstudio/download/#download).
<br>
<br>

---
<br>

# The data

For a quick overview of the example data we'll be using and where it came from, we are going to work with a subset of the dataset published [here](https://www.frontiersin.org/articles/10.3389/fmicb.2015.01470/full). We were exploring an underwater mountain ~3 km down at the bottom of the Pacific Ocean that serves as a low-temperature (~5-10°C) hydrothermal venting site. This amplicon dataset was generated from DNA extracted from crushed basalts collected from across the mountain with the goal being to begin characterizing the microbial communities of these deep-sea rocks. No one had never been here before, so as is often the purpose of marker-gene sequencing like this, this was going to be our first, broad survey. The sequencing was done on the Illumina MiSeq platform with 2x300 paired-end sequencing using primers targeting the V4 region of the 16S rRNA gene. There are 20 samples total: 4 extraction blanks (nothing added to DNA extraction kit), 2 bottom water samples, 13 rocks, and one biofilm scraped off of a rock. None of these details are important for you to remember, it's just to give some overview.  

In the following figure, overlain on the map are the rock sample collection locations, and the panes on the right show examples of the 3 distinct types of rocks collected: 1) basalts with highly altered, thick outer rinds (>1 cm); 2) basalts that were smooth, glassy, thin exteriors (~1-2 mm); and 3) one calcified carbonate.

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

## Merging forward and reverse reads

Depending on the sequencing facility, when paired-end sequencing was performed you may receive one forward and one reverse fastq file with all samples mixed together, or you may get things already demultiplexed (separate files for each individual sample). In this case our samples have already been demultiplexed and the barcodes have been trimmed off. This means for each of our 20 samples we have a forward (R1) and reverse (R2) reads file:

<center><img src="{{ site.url }}/images/ls_all_sample_reads.png"></center>

<br>
The format of the read file names is a letter or 2 specifying which sample it is, followed by a 'sub' (just because these are subsampled from the full file sizes), and then an R1 or R2. All of the .fq files are fastq files (sequences with quality score information), and the .fa files are simply sequences with no quality information. For some basic information about these formats visit the [Amplicon main page]({{ site.url }}/amplicon/).

So our first step is to merge these forward and reverse reads for each sample. This particular sequencing run was 2x300 paired-end over the V4 region, which is roughly ~300 base pairs (bps). So we have virtually full overlap for all of our sequences. The benefit of full-overlap as compared to sequencing a longer region with less overlap is improved base-calling accuracy. The usearch command `-fastq_mergepairs` as run below does a few things: 1) merges the paired reads for each sample into a single sequence; 2) appends a sample name to the merged sequence header to identify which sample it originated from; 3) combines and merged sequences from all samples into one file; 4) calculates new quality scores for each position of the merged sequences (as detailed [here](https://www.drive5.com/usearch/manual/merge_pair.html) if interested). We could list all of the files for each sample, but instead we're going to use a nice little *bash* shortcut called a wildcard. If you're unfamiliar with these, shoot over to the [wildcards section]({{ site.url }}/bash/basics/#wildcards) of [*bash* basics](/bash/basics) when you have some time to learn about them. For now, you just need to know the `*` represents anything any number of times, so as we use it here it will grab all of the files we need, one at a time.

```
usearch -fastq_mergepairs *_R1.fq -fastqout all_samples_merged.fq -relabel @
```

So here the magic of the `*` wild card is grabbing all forward read files for every file that ends in .fq. We didn't need to specify the reverse reads because usearch finds them automatically for us. We specify the output file with the `-fastqout` argument. And the last part, `-relabel @`, tells the program to append the file name the merged sequence originated from to its header. Along with the command finishing we get some overall summary statistics printed to the terminal, such as how many total sequences we had and what proportion of the paired reads successfully merged. Checks like this are needed in order to have some idea of what's going on with your sequences as you move them along. Here almost 90% of them successfully merged, which is pretty good.

## Cutting off primers
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

## Quality filtering
Quality filtering is a critical in reducing the abundance and impact of spurious sequences. There is an intrinsic error rate to all sequencing technology (and polymerases) that will consistently be generating some portion of sequences that vary slightly from their true biological origin, and this can substantially inflate metrics such as richness and diversity. Quality filtering is one of the steps in place to mitigate that problem. In usearch it is done with the `-fastq_filter` command, which uses calculated expected errors to determine if sequences are filtered out or not (the details can be found [here](https://www.drive5.com/usearch/manual/exp_errs.html)).

```
usearch -fastq_filter all_merged_no_primers.fq -fastq_maxee 1 -fastaout QCd_merged.fa
```

The output tells us very few reads were discarded due to poor quality. And note that our output file from this is now a fasta file, as after quality filtering we no longer need to keep track of the quality scores for each sequence anymore. 

## Dereplication
The dereplication step collapses all identical sequences to one and simply keeps track of how many there were. This serves a couple of purposes. It can save a ton of time in further processing steps because, for instance, you can then just process once a sequence that may appear 10,000 times – rather than processing 10,000 copies of the same exact sequence 10,000 times. And also, further attempts to mitigate the presence of spurious sequences involves merging low-abundance sequences with high-abundance sequences if they are similar enough, which requires knowing the abundance of individual sequences. In usearch this dereplication step is done with the `-fastx_uniques` command:

```
usearch -fastx_uniques QCd_merged.fa -sizeout -relabel Uniq -fastaout unique_seqs.fa
```

We can see from this output that of ~188,000 total sequences, ~95,000 of them were unique, and ~84,000 of those appear only one time (singletons). Having relatively so many singletons is common with amplicon data and is a direct consequence of intrinsic sequencing error rates. It is likely the majority of these are simply 1 or 2 bps diverged from a true biological sequence that is present in greater abundance. This is why we try to merge low-abundance sequences with high-abundance ones that are very similar – this is also one of the reasons why clustering sequences into arbitrary OTUs took off as well as it did, to mitigate sequencing error inflating true biological signals). 

## Clustering OTUs and/or generating ASVs
OTUs refer to 'operational taxonomic units' that are made by grouping sequences into clusters based on similarity. Then for each cluster, a sequence from that cluster is selected to be the representative sequence of that cluster (typically one from the 'center' of that cluster in 'sequence space', or the most abundant from the cluster). In usearch, if you specify 97% OTUs for instance, the program attempts to build OTUs (clusters) such that none of the representative sequences for the OTUs are more similar than 97% to any other.

ASVs on the otherhand refer to 'amplicon sequence variants', which seems to be emerging as the consensus name for sequences derived using single-nucleotide resolution. As the name suggests, single-nucleotide resolution means you can delineate sequences even if they vary by only one base pair. This most often means similarity levels greater than 99% that as such would be lost with any form of OTU clustering. ASVs are believed to represent true biological sequences (with sequencing error for the most part weeded out via abundance-based merging and filtering), and have the advantage of not throwing away resolution by merging similar sequences together. There are several tools available that apply different approaches to do this such as [MED](http://merenlab.org/2014/11/04/med/), [DADA2](https://benjjneb.github.io/dada2/index.html), and usearch's `-unoise3` command. 

For now we will just move forward with ASVs, but feel free to experiment. When I process a new tag dataset I usual run both analyses alongside each other. There can be instances where greater resolution isn't helpful. As usual it depends on the dataset and on what you are trying to find. (Though I'll note there are some people that are pretty adamant that ASVs should be the only way to go moving forward.)

So, to generate our ASVs:

```
usearch -unoise3 unique_seqs.fa -zotus ASVs.fa -tabbedout unoise3.txt
```

By default, usearch's `-unoise3` command filters out any that are less abundant than 8, and removes sequences it suspects are chimeric, details of which can be found [here](https://www.biorxiv.org/content/biorxiv/early/2016/09/09/074252.full.pdf). From the output to the terminal we can see we have ~1700 ASVs. In theory, these are true biological sequences recovered from our 20 samples.

## Generating a count table
The count table is what tells us how many times each sequence appears in each sample. It is the end-product of all this processing that we can't wait to get into R. The way the count table is generated in usearch is a bit different than most other approaches. Here, now that usearch has identified what it believes to be true biological sequences, the `-otutab` command is used to attempt to map all of our merged sequences to the ASVs. By default it does this by requiring a merged sequence to be >97% similar to an ASV, and if a sequence crosses that threshold for more than one ASV it is counted for the one it is most similar to only. Incorporating another 97% threshold here may seem a bit confusing at first, but usearch's thinking is that more often than not, the 'true' ASV sequences recovered were more than likely the source of the majority of sequences that are within 3% of them, as the source of this 3% variation is believed to be the result of sequencing and pcr errors. Additionally, it is still fundamentally different to have your ASV units built without similarity clustering, as opposed to 97% OTUs. That said, when I ran through this I decided to bump up the required % similarity to 98.5% when mapping our merged sequences to our ASVs to generate our count table. The downside to this is you will lose some data (as less sequences will map to an ASV), but if you have a decent amount of sequences per sample then you can afford to lose some for the sake of being a bit more conservative. 

Due to a minor glitch in usearch where it doesn't seem to register the headers of the ASV sequences it gave them, we need to run this quick `sed` command before the `-otutab` program:

```
sed -i.tmp 's/Zotu/ASV_/' ASVs.fa
usearch -otutab all_merged_no_primers.fq -zotus ASVs.fa -otutabout ASV_counts.txt -id 0.985
```
Note that the input file is our previous merged sequences file after we trimmed off the primers. The output to the screen after it finishes shows us we mapped ~70% of our total sequences. And the "ASV_counts.txt" file is our glorious count table. And if we take a peek at it with `less -S ASV_counts.txt` we can see it is organized with samples as column names, ASV sequence IDs as row names, and their counts per sample in the cells:

<center><img src="{{ site.url }}/images/less_ASV_counts.png"></center>

<br>

## Assigning taxonomy
The final step in our sequence processing is to assign taxonomy to our ASV sequences. There are multiple ways to do this as well. I usually run a few different ways and put them all in a table next to each other, and then I can be more somewhat more confident when there is some consistency. But something important to keep in mind is that taxonomy-assigning software is generally built to be fast. To do this, most rely on kmer frequencies to classify to a reference, rather than alignments. If you begin analyzing your data and some particular OTUs or ASVs emerge as being imporant to the overall story, it's a good idea to take those sequences and BLAST them to try to get more information about it. 

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
This portion assumes you already have some baseline experience with R – meaning we won't be breaking down any syntax or going over any basics of R here (that's presented elsewhere). But even if you don't have any experience with R yet, you'll still be able to follow along here and run everything, you'll just be on your own as far as figuring out how the code works. A full R script of everything done here is available in the "R_working_dir" subdirectory called "amplicon_example_analysis.R" that can be opened in R if you prefer to follow along with that rather than copying and pasting commands from here. Either way, this isn't about the R code (right now).

## Setting up our working environment
To get started let's open up RStudio. This portion assumes you have some baseline experience with R already. If you don't, you'll still be able to follow along and piece things together set our current working directory to where we just left our new files, and install the packages we'll need:

```R
setwd("~/amplicon_example_workflow/R_working_dir/")

install.packages("phyloseq")
install.packages("vegan")
install.packages("DESeq2")
install.packages("ggplot2")
```

If any of those didn't succeed, for whichever caused a problem run the corresponding command below to try installing it a different way:

```R
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
biocLite("DESeq2")

install.packages("devtools")
devtools::install_github("vegandevs/vegan")
devtools::install_github("tidyverse/ggplot2")
```

After all 3 have installed successfully, load them up:

```R
library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
```

If you have a problem loading any of these libraries, close and then restart R and try loading the library again. If still no luck, try installing the package again and loading the library. 9 times out of 10, restarting R will solve the problem, so hopefully you're not special :) 

 
## Reading in our data
Now we're going to read in our counts table, our taxonomy table, and a table with information about each sample including temperature of venting fluids at the site of collection ("temp"), the "type" of sample it is, and a qualitative characteristic ("char") of the rocks.

```R
count_tab <- read.table("ASV_counts.txt", header=T, row.names=1, check.names=F)

tax_tab <- as.matrix(read.table("ASV_tax.txt", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))

sample_info_tab <- read.table("sample_info.txt", header=T, row.names=1, check.names=F)
```

## Treatment of "blanks" 
First we want to deal with our "blanks" (those samples labeled B1-B4). In this case these refer to extraction blanks – meaning when doing the DNA extractions, we ran some extractions alongside the normal ones where we didn't add anything to the kit's starting lysis tubes. These are processed the same as the real samples and sent for sequencing as well. You may or may not have these with your data, if you didn't do it this time, be sure to next time you do an amplicon sequencing run :) 

Sequences that show up in the blanks can be the result of over-amplification of contaminant DNA introduced from things like the DNA extraction or sequencing kit reagents, **or** they can be cross-contamination from real samples into the blanks. As such, strictly throwing away any ASV that shows up at all in the blanks (i.e. disregarding that ASV in all samples) is not wise. Imagine an ASV that contains 90,000 out of 100,000 sequences in a true sample, and has 6 sequences out of 20,000 in a blank. In such a case it is almost certain those sequences came to be in the blank due to cross-contamination from the real sample. I've seen similar cases to this in every dataset I've looked at. 

As with most things, there is no one way or best way to do this. What I do (and did in the [paper](https://www.frontiersin.org/articles/10.3389/fmicb.2015.01470/full) these data come from) is compare normalized totals of the ASVs (or OTUs) in the blanks and in the true samples, and use an arbitrary threshold to decide if it is more likely to be contamination from reagents or more likely to have ended up in the blanks due to cross-contamination from a real sample. The arbitrary threshold I apply is that if the sample-normalized count is more than an order of magnitude greater than the blank-normalized count for a given ASV/OTU, then it is kept – if it is not, then it is thrown out and presumed to have been contamination. Here's what this looks like in this case:

```R
  # first we need to get a sum for each ASV across all 4 blanks and all 16 samples
blank_ASV_counts <- rowSums(count_tab[,1:4])
sample_ASV_counts <- rowSums(count_tab[,5:20])

  # now we normalize them, here by dividing the samples' total by 4 – as there are 4x as many samples (16) as there are blanks (4)
norm_sample_ASV_counts <- sample_ASV_counts/4

  # here we're getting which ASVs are deemed likely contaminants based on the threshold noted above:
blank_ASVs <- names(blank_ASV_counts[blank_ASV_counts * 10 > norm_sample_ASV_counts])
length(blank_ASVs) # this approach identifies 54 out of 1720 that are likely to have orginated from contamination

  # looking at the percentage of reads retained for each sample after removing these presumed contaminant ASVs shows that the blanks lost almost all of their sequences, while the samples, other than one of the bottom water samples, lost less than 1% of their sequences, as would be hoped
colSums(count_tab[!rownames(count_tab) %in% blank_ASVs, ]) / colSums(count_tab) * 100
        B1         B2         B3         B4        BW1        BW2        R10      R11BF        R11        R12        R1A        R1B         R2         R3         R4         R5         R6         R7         R8         R9 
 0.6298111  0.4008016  0.9153318  2.4038462 65.4151243 99.5216583 99.7387794 99.8166819 99.4470365 99.5028249 99.7501644 99.7146351 99.6506941 99.4810831 99.6432515 99.6186047 99.8327013 99.0976210 99.7623009 99.2768761
 
  # now that we've used our extraction blanks to remove ASVs that were likely due to contamination, we're going to trim down our count table by removing those sequences and the blank samples from further analysis
filt_count_tab <- count_tab[!rownames(count_tab) %in% blank_ASVs, -c(1:4)]
  # make a filtered sample info table
filt_sample_info_tab<-sample_info_tab[-c(1:4), ]

  # and let's add some colors to the sample info table that are specific to sample types and characteristics that we can use when plotting things
    # we'll color the water samples blue: 
filt_sample_info_tab$color[filt_sample_info_tab$char == "water"] <- "blue"
    # the biofilm sample a darkgreen:
filt_sample_info_tab$color[filt_sample_info_tab$char == "biofilm"] <- "darkgreen"
    # the basalts with highly altered, thick outer rinds (>1 cm) brown ("chocolate4" is the best brown I can find...):
filt_sample_info_tab$color[filt_sample_info_tab$char == "altered"] <- "chocolate4"
    # the basalts with smooth, glassy, thin exteriors black:
filt_sample_info_tab$color[filt_sample_info_tab$char == "glassy"] <- "black"
    # and the calcified carbonate sample an ugly yellow:
filt_sample_info_tab$color[filt_sample_info_tab$char == "carbonate"] <- "darkkhaki"
```

## Beta diversity
Beta diversity involves calculating metrics such as distances or dissimilarities based on pairwise comparisons of samples. Typically the first thing I do when I get a new dataset into R (whether it's tag data, gene expression data, methylation levels, or pretty much anything) is generate some exploratory visualizations like ordinations and hierarchical clusterings. These give you a quick overview of how your samples relate to each other and can be a way to check for problems like batch effects. 

We're going to use Euclidean distances to generate some exploratory visualizations of our samples. Since differences in sampling depths between samples can influence distance/dissimilarity metrics, we first need to somehow normalize across our samples. 

<h4><center>Normalizing for sampling depth</center></h4>
Common ways to do this involve either subsampling each sample down the the lowest sample's depth, or turning counts into proportions of the total for each sample. However, both of these approaches are generally shunned by people I trust when it comes to such topics (i.e., statisticians). For example, in their 2014 PLOS Computational Biology paper, ["Waste not, want not: why rarefying microbiome data is inadmissible"](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531), McMurdie and Holmes argue that a better method of normalizing across samples is to use a variance stabilizing transformation – which fortunately we can do with the [*DESeq2*](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) package.

```R
  # first we need to make a DESeq2 object
deseq_counts <- DESeqDataSetFromMatrix(filt_count_tab, colData = filt_sample_info_tab, design = ~type) # we have to include the colData and design arguments because they are required, as they are needed for further downstream processing by DESeq2, but for our purposes they don't matter
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

  # pulling out our transformed table 
vst_trans_count_tab <- assay(deseq_counts_vst)

  # and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))
```

<h4><center>Hierarchical clustering</center></h4>
Now that we have our Euclidean distance matrix, let's make and plot a hierarchical clustering of our samples.

```R
euc_clust <- hclust(euc_dist, method="ward.D2")

# plot(euc_clust) # hclust objects like this can be plotted as is, but i like to change them to dendrograms for two reasons:
  # 1) it's easier to color the dendrogram plot by groups
  # 2) if you want you can rotate clusters with the rotate() function of the dendextend package

euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- filt_sample_info_tab$color[order.dendrogram(euc_dend)]
labels_colors(euc_dend) <- dend_cols

plot(euc_dend, ylab="VST Euc. dist.")
```

<center><img src="{{ site.url }}/images/hclust.png"></center>

<br>
So from our first peek, the broadest clusters separate the biofilm and water samples from the rocks. The next tier splits the rocks into two groups, with samples R8–R11 separate from the others (R1-R6, and R12). As we can see from the coloring we added, R8-R11 were all of the glassier type of basalt with thin (~1-2 mm), smooth exteriors, while the rest (R1-R6, and R12) had more highly altered, thick (>1 cm) outer rinds (excluding the oddball carbonate which isn't a basalt, R7). This is starting to suggest that level of alteration of the basalt may be correlated with community structure. If we look at the map figure again, we can also see that level of alteration also co-varies with whether samples were collected from the northern or southern end of the outcrop as all of the more highly altered basalts were collected from the northern end.

<center><img src="{{ site.url }}/images/dorado.png"></center>

<br>

<h4><center>Ordination</center></h4>
Generally speaking, ordinations provide visualizations of sample-relatedness based on dimension reduction – this is where the 'multidimensional scaling' term (MDS) fits in. The 'dimensions' could be, for instance, whatever you measured in each sample. Principle coordinates analysis (PCoA) is a type of multidimensional scaling is a type of multidimensional scaling that operates on dissimilarities or distances. Here we're going to generate and plot our PCoA with *phyloseq*, because it is very convenient for doing such things. But because we're still doing beta diversity here, we want to use our transformed table. So we're going to make a new *phyloseq* object with our *DESeq2*-transformed table and generate the PCoA from that.

```R
  # making new phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

  # generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples


plot_ordination(vst_physeq, vst_pcoa, color="char") + 
  labs(col="type") + geom_point(size=1) + 
  geom_text(aes(label=rownames(filt_sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(filt_sample_info_tab$color[order(filt_sample_info_tab$char)])) + 
  theme(legend.position="none")
```

<center><img src="{{ site.url }}/images/pcoa.png"></center>

<br>
This is just providing us with a different initial overview of how our samples relate to each other. Inferences that are consistent with the hierarchical clustering above can be considered a bit more robust if the same general trend emerges from both approaches. It's important to remember that these are **exploratory visualizations** and do not say anything statistically about our samples. But our initial exploration here shows us the rock microbial communities seem to be more similar to each other than they are to the water samples (which is what I wanted to see, so this was a good start to me), and focusing on just the basalts (brown and black labels), these visualizations both suggest their communities may correlate with level of exterior alteration.

## Alpha diversity
Alpha diversity entails using summary metrics that describe individual samples, and it is a very tricky thing when working with amplicon data. There are a lot of tools from macro-ecology that have been co-opted into the microbial ecology world that just simply do not work the same. If and when I use any alpha diversity metrics, I mostly consider them useful for relative comparisons of samples from the same experiment. And I am absolutely going to ask some experts on the subject (ahem, @SherlockPHolmes and @AmyDWallis) to take a look at this part to be certain I'm not leading anyone astray here. But here are a few things I've done before (unless I'm later informed this is wrong – then these were always meant to be examples of what *not* to do).

<h4><center>Rarefaction curves</center></h4>
It is not okay to use rarefaction curves to estimate total richness of a sample, or to extrapolate anything from them really, but I believe they can still be useful depending on the data. Let's generate the plot and then we'll see why with this example. We'll be using the `rarecurve()` function from the package [*vegan*](https://github.com/vegandevs/vegan) here. Note that *vegan* expects rows to be samples and observations (our ASVs here) to be columns, which is why we transpose the table in the command with `t()`.

```R
rarecurve(t(filt_count_tab), step=100, col=filt_sample_info_tab$color, lwd=2, ylab="ASVs")
abline(v=(min(rowSums(t(filt_count_tab)))))
```
<center><img src="{{ site.url }}/images/rarefaction.png"></center>

In this plot, samples are colored the same way as above, and the black vertical line represents the sampling depth of the sample with the least amount of sequences (a bottom water sample, BW1, in this case). This view suggests that the rock samples (all of them) are more diverse and have a greater richness than the water samples or the biofilm sample (based on where they all cross the vertical line of lowest sampling depth). And again, just focusing on the brown and black lines for the two types of basalts we have, they seem to show similar trends within their respective groups that suggest the more highly altered basalts (brown lines) may host more diverse microbial communities than the glassier basalts (black lines). 

<h4><center>Richness and diversity estimates</center></h4>
Here we're going to plot Chao1 richness esimates and Shannon diversity values. Chao1 is a richness estimator, "richness" being the total number of distinct units in your sample, "distinct units" being whatever you happen to be measuring (ASVs in our case here). And Shannon's diversity index is a metric of diversity. The term diversity includes "richness" (the total number of your distinct units) and "evenness" (the relative proportions of all of your distinct units).

We are going to go back to using the [*phyloseq*]() package for this to use the function `plot_richness()` – which the developers kindly provide some examples of [here](https://joey711.github.io/phyloseq/plot_richness-examples.html). 

```R
  # first we need to create a phyloseq object using our un-transformed count table
count_tab_phy <- otu_table(filt_count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)
sample_info_tab_phy <- sample_data(filt_sample_info_tab)

ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

  # and now we can call the plot_richness() function on our phyloseq object
plot_richness(ASV_physeq, color="char", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(filt_sample_info_tab$color[order(filt_sample_info_tab$char)])) +
  theme(legend.title = element_blank())
```

<center><img src="{{ site.url }}/images/plot_richness.png"></center>

Before saying anything about this I'd like to stress again that these are *not* interpretable as "real" numbers of anything (due to the nature of amplicon data), but they can still be useful as relative metrics of comparison. For example, we again see from this that the more highly altered basalts seem to host communities that are more diverse and have a higher richness than the other rocks, and that the water and biofilm samples are less diverse than the rocks. 

And just for another quick example of why *phyloseq* is pretty awesome, let's plot grouping sample types while still coloring by the characteristics column:

```R
plot_richness(ASV_physeq, x="type", color="char", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(filt_sample_info_tab$color[order(filt_sample_info_tab$char)])) +
  theme(legend.title = element_blank())
```

<center><img src="{{ site.url }}/images/plot_richness_by_type.png"></center>

## Taxonomic summaries
Don't forget that the taxonomy called here was done rapidly and by default has to sacrifice some specificity for speed. For the sequences that become important in your story, you should absolutely pull them out and BLAST them, and possibly make phylogenetic trees to get a more robust idea of who they are most closely related to. 

Let's start by making some broad summarization figures. *Phyloseq* is also very useful for parsing things down by taxonomy now that we've got all that information in there already.



## Permutational ANOVA

As we mentioned earlier, we have some information about our samples as well in a table. There are many ways to incorporate this information, but one of the first I typically go to is a permutational ANOVA test to if any of the sample information is predictive of community structure. One way to do this is with the `adonis` function from the *vegan* package. If we look at the help menu for the function, we see this description:

```R
?adonis
```

<center><img src="{{ site.url }}/images/adonis_desc.png"></center>

 
<br>
<br>

---
<br>
# So what now?
Well, now is when you do the science part. Here's where your questions and the experimental design start to guide how you go further. In this paper for instance we ended up incorporating other seafloor basalt studies to identify what looked to be conserved taxa that showed up everywhere (like a sulfur-oxidizing gammaproteobacterium, *Thioprofundum lithotrophicum*, and we also identified that there seems to be a basalt-hosted Thaumarchaeota distinct from those present in the bottom water samples we analyzed (*Nitrosopumilus* sp.). To see more of how this dataset ended up, check out the results and discussion of the [paper](https://www.frontiersin.org/articles/10.3389/fmicb.2015.01470/full). 
<br>
<br>

<center><img src="{{ site.url }}/images/under_construction.jpeg"></center>
<center><h3>UNDER CONSTRUCTION</h3></center>
