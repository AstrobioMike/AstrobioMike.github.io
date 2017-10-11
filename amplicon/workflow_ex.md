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
There are many ways to process amplicon data (if you need a quick primer on some relevant terminology, visit the [Amplicon main page]({{ site.url }}/amplicon/)). Some of the most widely used tools/pipelines include [mothur](https://www.mothur.org/), [usearch](https://drive5.com/usearch/), [Minimum Entropy Decomposition](http://merenlab.org/2014/11/04/med/), [DADA2](https://benjjneb.github.io/dada2/index.html), and [qiime2](https://qiime2.org/) (which employs other tools within it). As usual, there is no one-size-fits-all, there is no 'best'. And actually in my experience if you make similar decisions when processing your sequences (decisions about things like minimum abundance filtering or clustering thresholds), you most often get very similar results regardless of which of these you use (refreshingly). Here I'll be using [usearch](https://drive5.com/usearch/) mostly because of it's ease of deployment; there is no installation required and there are no dependencies. So if you want to actively follow along with this module you can just download it and you're ready to rock. But please keep in mind that nothing here is meant to be authoritative. This is simply one example of one workflow. When working with your own data you should never follow any pipeline blindly.  

Most often a marker-gene analysis is the microbial ecologist's first tool in a vast toolkit. It is primarily used as a broad survey of community structure. As the warning notes above, it is easy to get caught spinning your wheels about a sublte component in your processing pipeline that ultimately has a negligible impact compared to the noise we are working through. What I mean by this is, generally speaking, tag data is not the appropriate tool to answer really meticulous questions. It is a tool for comparing baseline *proxies* of metrics about microbial communities. It is a tool of exploration and hypothesis generation, not hypothesis confirmation.  

There are a lot of things to keep in mind regarding what tag sequencing means, what it doesn't mean, what it can tell you, and what it can't. And at some point I'll give my two cents on all of these things in the [Amplicon Thoughts]({{ site.url}}/amplicon/thoughts) post should anyone be interested. But for now, let's go through a dataset.
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

This directory contains: 1) forward and reverse fastq files for each of the 20 samples; 2) a primers fasta file; 3) the [RDP](https://rdp.cme.msu.edu/) reference fasta we will use for assigning taxonomy; 4) a text file with all usearch commands we will use; and 5) an R subdirectory containing a text file of some sample information and a .R file with all commands we will use in R.  

Now, let's get started!
<br>
<br>

---
<br>

# Processing

It's good to try to keep a bird's-eye view of what's going on. So here is an overview of the main processing steps we'll be performing with usearch. Don't worry if anything seems unclear right now, we will discuss each at each step.

||Command|What we're doing|
|:--:|:--------:|----------|
|1|fastq_mergepairs|merge forward and reverse reads together|
|2|fastx_truncate|cut off forward and reverse primers|
|3|fastq_filter|quality filter sequences|
|4|fastx_uniques|dereplicate sequences|
|5|cluster_otus/unoise3|cluster sequences into OTUs and/or generate ASVs|
|6|otutab|generate a counts table|
|7|sintax|assign taxonomy to OTUs and/or ASVs|

<h4>Merging forward and reverse reads</h4>

The way 
