---
layout: main
title: Phylogenomics
categories: [genomics, metagenomics, phylogenomics]
tags: [phylogenomics,phylogenetics,genomics,metagenomics,assembly]
permalink: /genomics/phylogenomics
---

{% include _genomics_phylogenomics_toc.html %}

{% include _side_tab_genomics.html %}

Here we are going to talk about what phylogenomics is (as the term is most often used today) and what goes into generating a phylogenomic tree. Then we're going to create some with [GToTree](https://github.com/AstrobioMike/GToTree/wiki/what-is-gtotree%3F){:target="_blank"} 🙂

<center><a href="../images/hug-tree-of-life.png"><img width="100%" src="../images/hug-tree-of-life.png"></a></center>

<p align="right" style="font-size:0.9em">the instantly iconic "new view of the Tree of Life" from <a href="https://www.nature.com/articles/nmicrobiol201648" target="_blank">Hug et al. 2016</a></p>

<hr style="height:10px; visibility:hidden;" />
---
---
<br>
# What is phylogenomics?
**Put into one over-simplified, slightly misleading, but nonetheless conceptually usefull sentence: *phylogenomics* is attempting to infer evolutionary relationships at the genomic level.** This is over-simplified and slightly misleading because in practice we are never using the "entire" genomes of all of the organisms we wish to focus on (and depending on the breadth of diversity we are considering, it would be impossible and/or meaningless to use their entire genomes because they might be too different). So really, it is more appropriate to say: ***phylogenomics* is attempting to infer evolutionary relationships at something closer to the genome-level than an individual gene phylogeny gets us** (like a 16S rRNA gene tree). 

Most *phylogenetic* trees that biologists are used to seeing and working with are visual representations of the estimated evolutionary relationships of various copies of a single gene-type (like the 16S rRNA gene). When we do this, if we are trying at all to think on the organism level, we are using that gene as a proxy to stand in for the organism itself, and as such we are assuming the evolutionary relationships of those proxies do indeed tell us something meaningful about the evolutionary relationships of their source organisms.

We are doing the same thing with phylogenomics, just instead of a single gene, we are using multiple genes. But in the end, we are still using sequences as proxies. So an even *more* appropriate definition might be: ***phylogenomics* is attempting to infer evolutionary relationships between sequences comprised of multiple concatenated genes, while assuming those inferred evolutionary relationships tell us something meaningful with regard to the evolutionary relationships of their source genomes**. But that can sound pretty convoluted if we're just starting to get into this mental space, hence the over-simplified, slightly misleading first version 🙂

Now, with those caveats laid out there in that last definition, I want to immediately follow up with saying for the most part, it *is* a safe assumption that the inferred evolutionary relationships of the multi-gene sequences we are getting *are* meaningfully representative of their source genomes' evolutionary relationships. And this generally becomes a safer assumption when more *appropriate* genes are included, and a less safer assumption if more *inappropriate* genes are included. Which and how many genes are appropriate entirely depends on the breadth of diversity we are interested in spanning. In general, for a given set of organisms of interest, it'd be a good idea to use as many *single-copy core genes* as that set of organisms has (though this certainly may not always be necessary in order to get a trusted signal that meets our needs).

## What are single-copy core genes? 
One of my responses in a [GToTree issue](https://github.com/AstrobioMike/GToTree/issues/3#issuecomment-573283418):
>One of the built-in assumptions in phylogenetics in general is that the genes being compared are under similar evolutionary pressures, which is tenuous to begin with even with single-copy genes, but gets a little worse when there are paralogs (multiple copies) within the same genome. 

Another [GToTree issue thread](https://github.com/AstrobioMike/GToTree/issues/3#issuecomment-581085307):
>So there are a couple things going on here. When GToTree has clade-specific marker sets, they are specifically for conserved, single-copy genes across that clade (more details on that in numbers 1-4 under the table here). That's different than these HMM sets that are available at PhyloFacts. It seems those are HMMs of all CDSs found in all Bacillus subtilis incorporated reference genomes (as you've pointed out). But we most often wouldn't want to use all CDSs when trying to make a phylogenomic tree, because one of the assumptions built into phylogenetics is that when considering the same gene from multiple organisms, it is presumed be under similar evolutionary pressures in those organisms. This is already a pretty hard thing to ever know for certain, and is probably never actually true, but to be safe, it is best to use genes that are present only in a single copy across the genomes being considered, and not use any that are in multiple copies within individual genomes.

## Overview of the whole process
Genomes -> 1:1 orthologs -> align genes -> trim -> concatenate -> tree


# Let's do some phylogenomics!

## Installing GToTree
A [binder](https://mybinder.org/){:target="_blank"} is available to work in by clicking this badge:  (info on getting into the binder command-line environment can be found [here](/unix/getting-started#accessing-our-command-line-environment){:target="_blank"}, but be sure to activate the binder by clicking on the badge on this page, not that one). If using the binder, we can skip the GToTree installation step as it is already installed there.

If working on your system, we can create a [conda](/unix/conda-intro){:target="_blank"} environment to install [GToTree](https://github.com/AstrobioMike/GToTree/wiki/what-is-gtotree%3F){:target="_blank"} and activate it like so: 

```bash
conda create -y -n gtotree -c conda-forge -c bioconda -c defaults -c astrobiomike gtotree

conda activate gtotree
```

## Alteromonas example

## Tree of Life example

## Mixed-model example

## Tree visualization and manipulation

# The level of diversity of our arbitrary window matters

