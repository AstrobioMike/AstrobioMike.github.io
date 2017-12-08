---
layout: main
title: Quality filtering, error correction, and depth normalization
categories: [genomics]
permalink: /genomics/where_to_start
---

{% include _genomics_where_to_start_toc.html %}

{% include _side_tab_genomics.html %}

There are many things you can do to your sequence data before you being any type of analysis, and they can have a tremendous impact on what you're capable of doing with your data, especially if you will be doing any sort of assembly. As with most things in the bioinformatics world, there is no one-size-fits-all SOP for processing genomic sequencing data. Even when we are talking about individual genome sequencing, lots of things can add up to different datasets requring different processing steps or parameters to allow us to pull out the most (and more importantly most accurate) information that we can – things like the intricacies of the genome itself or the quality of the run when it was sequenced, for example. Here we'll look at some of the things you can do to your data in the initial processing steps. Help for installing the following tools can be found [here](/bash/installing_tools), and usage of some are demonstrated in the [*de novo* genome assembly and initial probing page](/genomics/de_novo_assembly).   
<br>

---
---
<br>

# Quality filtering
Typically you will get your sequencing data back from the sequencing facility in fastq formatted files. The defined fastq format is 4 lines per sequence: 1) the sequence identifier (header), preceded by a "@" character; 2) the sequence; 3) a "+" character and possibly the header information repeated; and 4) the quality score information for each individual basecall. With Illumina sequencing, the quality score information is a measure of how confident the software was when it called that particular base position whatever base it did. This isn't a perfect system, as there are still confounding factors like polymerase error and other systematic errors that won't show up in the quality score information, but nonetheless performing some quality-based filtering is essential.  

There is a very handy and user-friendly tool called [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) that is useful for getting an overview of your data's quality. Fastqc can help spot some commonly occurring problems, and can help guide the decisions you make when quality filtering.  

The tool I probably use most often for quality filtering due to its flexibility would be [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). Trimmomatic will help you do things like remove any remaining adapters that may be muddying your data and truncate sequences based on specific quality filtering thresholds. Another quality filtering program I use quite a bit is within the [illumina-utils](https://github.com/merenlab/illumina-utils) collection of scripts provided by [merenlab.org](http://merenlab.org/).  

Whatever you choose to use, quality filtering needs to happen, and it's a good place to start whenever you get a new dataset.  
<br>

---
<br>

# Read error correction
Probably not all that surprisingly, read error correction in the past has made incredible differences in the resulting assemblies and analyses I've been able to pull out of datasets. When I'm comparing assemblies (genome or metagenome) done under different parameters or with different programs, I usually throw one or two in that were assembled without an error correction step, and they are consistently worse (when considering generic assembly stats, but more importantly when considering recovered genomes).  

I'm sure there are many tools that perform this task – and please shoot me a message if you know of and prefer some others so I can add them here – but I haven't yet ventured further than the error correction available with the [SPAdes assembler](http://cab.spbu.ru/software/spades/). Overall SPAdes has given me excellent results with reconstructing genomes from axenic or very clean enrichment cultures, but it can become a bit too memory intensive with some more complicated samples like diverse metagenomes. In those cases, I still run my error correction step through SPAdes with the `--only-error-correction` flag set, and then I take the error corrected reads to another assembler.  

I imagine there may be some scenarios where error correction would hurt more than help (because all things seem to happen with data), and there may be some particular analysis you want to run where error correction might muddy the signal you're looking for, but barring any unusual context, I would advocate to pretty much always run a program that is designed to correct errors – especially if you are going to be assembling.  
<br>

---
<br>

# Depth normalization
I had hoped to return to the impact of sequencing depth on assembly quality to better flesh things out, but I just haven't been able to yet. All I can say for sure right now is that: 1) I was working with paired-end Illumina data from very clean enrichment cultures; 2) the depth of sequencing attained was insanely high on some of them (>1000 or even 2000X); and 3) normalizing by kmer depth with the [*bbnorm* script](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbnorm-guide/) offered within the [bbtools suite of programs](https://jgi.doe.gov/data-and-tools/bbtools/) substantially improved the quality of genomes I was able to recover as compared to when the depth was sky high.  

So while I haven't systematically looked at this with more than one dataset yet, and it probably wouldn't be an issue for most as it's not common to sequence to such crazy depths, but it's something to keep in mind and take a peek at when you're trying to get the best assembly you can.


