---
layout: main
title: Amplicon and metagenomics overview 
permalink: /misc/amplicon_and_metagen
---  

{% include _amp_and_meta_toc.html %}

This page presents a broad-level overview of amplicon sequencing and metagenomics as approaches for microbial ecology and outlines some generic workflows for each listing some relevant tools. Both of these methods are just two more tools of investigation and should be thought of as steps in the process of science rather than end-points (like all tools of science). They are most often applied for exploration and hypothesis generation. 


> "Why are you doing 16S sequencing? That doesn't tell you anything about function."  
> 
> "Why are you measuring nitrogen-fixation rates via acetylene reduction? That doesn't tell you anything about the biochemistry doing it."  
> 
> **Don't assess the utility of a tool based on something it's not supposed to do anyway.**

<br>

___
<br>

**Amplicon sequencing**  
Amplicon sequencing of marker-genes (e.g. 16S, 18S, ITS) involves using specific primers that target a specific gene or gene fragment. It is one of the first tools in the microbial ecologist's toolkit. It is most often a broad-level survey of community composition used to generate hypotheses based on differences between samples (***differences based on recovered gene-copy numbers, not abundance of organisms***).

**Metagenomics**  
Shotgun metagenomic sequencing provides a way to access *all* the DNA of a mixed community. It uses random rather than targeted primers and therefore suffers much less from pcr bias. (It still suffers from other things such as cell-lysis rates which can be dependent on organism cell walls and the extraction method used.)

<br>

___
<br>

# Some capabilities of each
<br>

## Amplicon
* **Useful for:**

    * one *metric* of community composition
        * can say something about relative abundance of gene copies recovered (recovered gene copies ≠ counts of organisms)
        * gives a snapshot of, e.g., "16S gene-fragment copies recovered"
    * can track changes in community structure (as interpreted by recovered gene copy numbers) in response to a treatment and/or across environmental gradients/time, etc.
    * can provide strong support for further investigation, particularly single-nucleotide resolution methods
        * e.g. *Trichodesmium–Alteromonas* story ([starting paper here](https://www.nature.com/articles/ismej201749){:target="_blank"})  
<br>

* **_Not_ useful for:**
    * abundance of organisms (or relative abundance of organisms)
        * recovered gene copies ≠ counts of organisms
            * gene-copy number varies per genome/organism (16S sequence can vary per genome)
            * pcr bias (small scale) -> under/over representation based on primer-binding efficiency
            * pcr bias (large scale) -> "universal" primers, only looking for what we know and don't even catch all of that 
    * function
        * Even if you can highly resolve the taxonomy of something from an amplicon sequence, it is still only one fragment of one gene, and extrapolating to a genome's functional potential is specious due to the highly variable nature of most microbes' accessory genomes. (hypothesis generation is ok, but state it as such)
        * There is no strain-level resolving capability for a single gene. Strain-level resolution needs to be done at the genome level. This does not eliminate the value of using a single-nucleotide method like dada2, but the devil is in the interpretational details.
        
## Metagenomics
* **Useful for:**
    * functional potential
    * insights into the "unculturables" 
    * much better for "relative" abundance (still not true abundance)
        * still some caveats, like cell-lysis efficiencies  
<br>
* **_Not_ useful for:**
    * "activity"
        * neither is transcriptomics or proteomics for that matter – Life is complicated 🙂
<br>

---
<br>
**So, with all that said, should we expect relative abundance of amplicon sequencing to match up with metagenomics?**

No. They are very different methods with each asking different questions and providing different information.

___
<br>
# General workflows

## Amplicon overview

<center><a href="{{ site.url }}/images/amplicon_overview.png"><img src="{{ site.url }}/images/amplicon_overview.png"></a></center>

<p align="right"><a href="https://ndownloader.figshare.com/files/12732065">PDF download</a></p>

<br>

#### A Note on OTUs vs ASVs  

All sequencing technologies make mistakes, and (to a much lesser extent), polymerases make mistakes as well. These mistakes artificially increase the number of unique sequences in a sample, a lot. Clustering similar sequences together (generating OTUs) emerged as one way to mitigate error and summarize data, though often at the cost of resolution. The field is moving towards using solely ASVs, and there is pretty good reasoning for this. This [paper](https://www.nature.com/articles/ismej2017119) nicely lays out the case for that, and the following points attempt to summarize it:  

* OTUs (operational taxonomic units)
    1. cluster sequences into groups based on percent similarity
    2. choose a representative sequence for that group
        * closed reference
            * **\+** can compare across studies
            * **\-** reference biased and constrained
        * de novo
            * **\+** can capture novel diversity
            * **\-** not comparable across studies
            * **\-** diversity of sample affects what OTUs are generated

* ASVs (amplicon sequence variants)
    1. attempt to identify the original biological sequences taking into account error
        * **\+** enables single-nucleotide resolution
        * **\+** can compare across studies
        * **\+** can capture novel diversity

<br>
## Metagenomics overview

<center><a href="{{ site.url }}/images/metagenomics_overview.png"><img src="{{ site.url }}/images/metagenomics_overview.png"></a></center>

<p align="right"><a href="https://ndownloader.figshare.com/files/12732062">PDF download</a></p>


<br>

___
<br>

# Some example tutorials

<h3><b>Amplicon</b></h3>

* **dada2** (ASVs)
    * [Developer tutorial](https://benjjneb.github.io/dada2/tutorial.html)
    * Mike's [dada2 tutorial](https://astrobiomike.github.io/amplicon/dada2_workflow_ex), with a separate section for dealing with 18S mixed in with 16S
* **usearch/vsearch** (ASVs and OTUs)
    * [usearch](https://www.drive5.com/usearch/) is not entirely free, but it has some very useful tools and approaches (e.g. good calculation of hybrid quality scores after merging overlapping reads). There is a lightweight free version that can still do many things.
    * [vsearch](https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline) is completely free and open, and was made in response to usearch not being completely free. It does not have all of the capabilities of usearch however.
    * Mike's [usearch/vsearch tutorial](https://astrobiomike.github.io/amplicon/workflow_ex)
* **mothur** (OTUs only currently)
    * [Developer tutorial](https://www.mothur.org/wiki/MiSeq_SOP)

* **qiime2** 
    * qiime provides an environment that employs other processing tools (like those above) and also provides convenient visualization capabilities and an excellent infrastructure for tracking everything you've done
    * they have [extensive documentation](https://docs.qiime2.org/2018.6/) and an [amplicon tutorial here](https://docs.qiime2.org/2018.6/tutorials/moving-pictures/)

<br>

---
<br>

<h3><b>Metagenomics</b></h3>

As you might guess, this is not as straightforward as the amplicon data tutorials as there are many more things to do with metagenomics data. But here are some places to start. If you know of others that are helpful please pass along so we can add them here 🙂

- Mike's recovering genomes from metagenomes [tutorial](https://astrobiomike.github.io/metagenomics/metagen_anvio)
	- [Here is a simplified example of how we recover genomes from metagenomes](/images/gen_from_metagen_slide.png); you can download a <a href="https://ndownloader.figshare.com/files/12367211">keynote slide here</a> and a <a href="https://ndownloader.figshare.com/files/12367226">powerpoint version here</a> if you'd like to use or edit 🙂  


- A nice workflow leading up to and including recovering genomes can be found [here](http://merenlab.org/tutorials/infant-gut/) at the [anvi'o site](http://merenlab.org/software/anvio/) (along with other very informative/helpful tutorials and blogs)
