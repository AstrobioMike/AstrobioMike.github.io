---
layout: main
title: Amplicon and metagenomics overview 
permalink: /misc/amplicon_and_metagen
---  

{% include _amp_and_meta_toc.html %}

>This page presents a broad-level overview of amplicon sequencing and metagenomics as applied to microbial ecology. Both of these methods are most often applied for exploration and hypothesis generation and should be thought of as steps in the process of science rather than end-points – like all tools of science 🙂 

<center><a href="../images/amplicon-vs-metagenomics-overview-w-background.png"><img src="../images/amplicon-vs-metagenomics-split-overview.png"></a></center>

<hr style="height:10px; visibility:hidden;" />
<br>

**Amplicon sequencing**  

Amplicon sequencing of marker-genes (e.g. 16S, 18S, ITS) involves using specific primers that target a specific gene or gene fragment. It is one of the first tools in the microbial ecologist's toolkit. It is most often used as a broad-level survey of community composition used to generate hypotheses based on differences between recovered gene-copy numbers between samples.

**Metagenomics**  

Shotgun metagenomic sequencing aims to amplify all the accessible DNA of a mixed community. It uses random primers and therefore suffers much less from pcr bias (discussed below). Metagenomics provides a window into the taxonomy and functional potential of a sample. Recently, the recovery of representative genomes from metagenomes has become a very powerful approach in microbial ecology, drastically expanding the known Tree of Life by granting us genomic access to as-yet unculturable microbial populations (e.g. [Hug et al. 2016](https://www.nature.com/articles/nmicrobiol201648); [Parks et al. 2017](https://www.nature.com/articles/s41564-017-0012-7)). 

Here we'll discuss some of the things each is useful and not useful for, and then look at some general workflows for each. 


## Amplicon sequencing utility
* **Useful for:**
    * one metric of community composition
        * can say something about relative abundance of gene copies recovered
    * can track changes in community structure (as interpreted by recovered gene copy numbers) in response to a treatment and/or across environmental gradients/time
    * can help generate hypotheses and provide support for further investigation of things 

* **_Not_ useful for:**
    * abundance of organisms (or relative abundance of organisms)
        * recovered gene copies ≠ counts of organisms
            * gene-copy number varies per genome/organism (16S sequence can vary per genome)
            * pcr bias (small scale) -> under/over representation based on primer-binding efficiency
            * pcr bias (large scale) -> "universal" primers, only looking for what we know and they don't even catch all of that 
            * cell-lysis efficiencies
    * function
        * even if we can highly resolve the taxonomy of something from an amplicon sequence, it is still only one fragment of one gene
        * hypothesis generation is okay, e.g. speculating about observed shifts in gene-copy numbers based on what's known about nearest relatives in order to guide further work
        * but, for example, writing a metagenomics paper based on 16S data would likely not be a good idea
        * There is no strain-level resolving capability for a single gene, all that tells you is how similar those genes are.
        
        
As noted above, amplicon data can still be very useful. Most often when people claim it isn't, they are assessing it based on things it's not supposed to do anyway, e.g., this type of question:

> "Why are you doing 16S sequencing? That doesn't tell you anything about function."  

To me is like this type of question:

> "Why are you measuring nitrogen-fixation rates? That doesn't tell you anything about the proteins that are doing it."  

**We shouldn't assess the utility of a tool based on something it's not supposed to do anyway 🙂**
        
## Metagenomics utility
* **Useful for:**
    * functional potential
    * insights into the genomes of as-yet unculturable microbes
    * much better for "relative" abundance due to no confounding copy-number problem and no drastic PCR bias (still not true abundance)
        * still some caveats, like cell-lysis efficiencies  

* **_Not_ useful for:**
    * abundance of organisms
    * "activity"
        * neither is transcriptomics or proteomics for that matter – each gives us insight into cellular regulation at different levels

<challengeBlock>
<center><b>QUICK QUESTION!</b></center>

With all that said, do you think we should expect relative abundance information from amplicon sequencing to match up with relative abundance from metagenomic sequencing? 
<br>

<div class="wrap-collabsible">
  <input id="q1" class="toggle" type="checkbox">
  <label for="q1" class="lbl-toggle">Solution</label>
  <div class="collapsible-content">
    <div class="content-inner">
		
No, and that's not a problem if we understand that neither are meant to tell us a true abundance anyway. They are each providing a different, foggy view into the system we sequenced. And the relative abundance metrics they do provide can still be informative when comparing multiple samples generated the same way 🙂

    </div>
  </div>
</div>
</challengeBlock>


# General workflows

## Amplicon overview

<center><a href="../images/amplicon_overview.png"><img src="../images/amplicon_overview.png"></a></center>

<p align="right"><a href="https://ndownloader.figshare.com/files/15628100">PDF download</a></p>

<br>

#### A Note on OTUs vs ASVs
All sequencing technologies make mistakes, and (to a lesser extent) polymerases make mistakes as well. Our initial methods of clustering similar sequences together based on percent similarity thresholds (generating **O**perational **T**axonomic **U**nits; OTUs) emerged as one way to mitigate these errors and to summarize data – along with a recognized sacrifice in resolution. What wasn't so recognized or understood at first, is that when processing with traditional OTU-clustering methods, these mistakes (along with a hefty contribution from chimeric sequences) artificially increase the number of unique sequences we see in a sample (what we often call "richness", though keep in mind we are counting *sequences*, not organisms) – often by **a lot** (e.g. [Edgar and Flyvbjerg 2015](https://academic.oup.com/bioinformatics/article/31/21/3476/194979){:target="_blank"}, [Callahan et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/){:target="_blank"}, [Edgar 2017](https://peerj.com/articles/3889/){:target="_blank"}, [Prodan et al. 2020](https://doi.org/10.1371/journal.pone.0227434){:target="_blank"}). 

Over the past decade, and particularly with greater frequency the past ~5 years, single-nucleotide resolution methods that directly try to better "denoise" the data (deal with errors) have been developed, e.g. *Oligotyping* [(Eren et al. 2013)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12114){:target="_blank"}, *Minimum Entropy Decomposition* [(Eren et al. 2014)](https://www.nature.com/articles/ismej2014195), *UNOISE* [(Edgar and Flyvbjerg 2015)](https://academic.oup.com/bioinformatics/article/31/21/3476/194979){:target="_blank"}, DADA2 [(Callahan et al. 2016)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/){:target="_blank"}, and deblur [(Amir et al. 2017)](https://msystems.asm.org/content/2/2/e00191-16){:target="_blank"}. These single-nucleotide resolution methods generate what we refer to as **A**mplicon **S**equence **V**ariants (ASVs). The field as a whole is moving towards using solely ASVs, and – in addition to being specifically designed to better deal with errors and successfully drastically reducing the number of spurious recovered sequences – there are pretty good practical reasons for this also. This [Callahan et al. 2017 paper](https://www.nature.com/articles/ismej2017119){:target="_blank"} nicely lays out the practical case for this, summarized in the following points:  

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
    1. attempt to identify the original biological sequences by taking into account error
        * **\+** enables single-nucleotide resolution
        * **\+** can compare across studies
        * **\+** can capture novel diversity

If you happen to work with amplicon data, and are unsure of what's going on with this whole hubbub between OTUs and ASVs, I highly recommend digging into the [Callahan et al. 2017 paper](https://www.nature.com/articles/ismej2017119){:target="_blank"} sometime as a good place to start 🙂

## Metagenomics overview

<center><a href="../images/metagenomics_overview.png"><img src="../images/metagenomics_overview.png"></a></center>

<p align="right"><a href="https://ndownloader.figshare.com/files/15628103">PDF download</a></p>

<br>

---
<br>
