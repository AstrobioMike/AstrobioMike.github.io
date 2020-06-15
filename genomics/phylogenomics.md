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

Most *phylogenetic* trees that biologists are used to seeing and working with are visual representations of the estimated evolutionary relationships of various copies of a single gene-type (like the 16S rRNA gene). When we do this, if we are trying at all to think on the organism level, we are using that gene as a proxy to stand in for the organism itself. As such, we are assuming the evolutionary relationships of those proxies do indeed tell us something meaningful about the evolutionary relationships of their source organisms.

We are doing the same thing with phylogenomics, just instead of a single gene, we are using multiple genes. But in the end, we are still using sequences as proxies. So an even *more* appropriate definition might be: ***phylogenomics* is attempting to infer evolutionary relationships between sequences comprised of multiple concatenated genes, while assuming those inferred evolutionary relationships tell us something meaningful with regard to the evolutionary relationships of their source genomes**. But that can sound pretty convoluted if we're just starting to get into this mental space, hence the over-simplified, slightly misleading version first 🙂

Now, with those caveats laid out there in that last definition, I want to immediately follow up with saying for the most part, it *is* a safe assumption that the inferred evolutionary relationships of the multi-gene sequences we are getting *are* meaningfully representative of their source genomes' evolutionary relationships. And this generally becomes a safer assumption when a greater number of *appropriate* genes are included, and a less safer assumption if fewer, or *inappropriate*, genes are included. Which and how many genes are appropriate entirely depends on the breadth of diversity we are interested in spanning (more on that below). In general, for a given set of organisms of interest, it'd be a good idea to use as many *single-copy core genes* as that set of organisms has (though this certainly may not always be necessary in order to get a trusted signal that meets our current needs).

## What are single-copy genes? 
*Single-copy core* genes, or just **s**ingle-**c**opy **g**enes (SCGs), are genes that are present in exactly 1 copy in most of the organisms we happen to be talking about. To compare genes, we need the different organisms we are considering to actually have the gene, that's why we need greater than 0 copies. But we also want genes in single-copy (rather than genes that tend to exist in multiple copies within our target organisms) because *one of the built-in assumptions in phylogenetics in general is that the genes being considered are under similar evolutionary pressures*. This is tenuous to begin with even with single-copy genes (and is probably never actually entirely true), but it becomes much more likely we are violating that assumption if there are paralogs (multiple gene copies that are the result of duplication). So that's why when we are talking about phylogenomics in general, SCGs play such a predominant role. 

## What genes should we use?
The organisms we are considering dictate the genes


## Overview of the whole process
Genomes -> 1:1 orthologs -> align genes -> trim -> concatenate -> tree


# Let's do some phylogenomics!

## Installing GToTree
A [binder](https://mybinder.org/){:target="_blank"} is available to work in by clicking this badge here [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AstrobioMike/binder-happy-belly-phylogenomics/master?urlpath=lab){:target="_blank"} and an example of getting into the binder command-line environment can be found [here](/unix/getting-started#accessing-our-command-line-environment){:target="_blank"} if needed (but be sure to activate the binder by clicking on the badge above, not the one on the other page). 

> **Note on using binder for this tutorial**  
> Binder is amazing, but those environments aren't built for speed, and most laptops will far outperform the binder environment on many tasks. If wanting to actively follow along with the following examples, they will go about 4x faster if working on your own system rather than the Binder environment. None of the examples are large or memory intensive, so if you have access to a [Unix-like environment](https://astrobiomike.github.io/unix/getting_unix_env){:target="_blank"}, I recommend working on your own system and installing `GToTree` as shown below.
> 
> Last note, if working in the binder, be sure to add the `-P` flag to the `GToTree` commands. Binder can't utilize `ftp`, so adding the `-P` flag uses `http` instead.

If using the binder, we can skip the GToTree installation step as it is already installed there.

If working on our system, we can create a [conda](/unix/conda-intro){:target="_blank"} environment, install [GToTree](https://github.com/AstrobioMike/GToTree/wiki/what-is-gtotree%3F){:target="_blank"}, and activate it like so: 

```bash
conda create -y -n gtotree -c conda-forge -c bioconda -c defaults -c astrobiomike gtotree

conda activate gtotree
```

## Alteromonas example

```bash
esearch -query 'Alteromonas[ORGN] AND "latest refseq"[filter] AND "complete genome"[filter] AND (latest[filter] AND all[filter] NOT anomalous[filter])' -db assembly | esummary | xtract -pattern DocumentSummary -def "NA" -element AssemblyAccession > alteromonas_refseq_accessions.txt
```

```bash
# getting "MAG"
curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/271/865/GCA_002271865.1_ASM227186v1/GCA_002271865.1_ASM227186v1_genomic.fna.gz | gunzip - > GCA_002271865.1.fa
# getting an Alpha to allow rooting
curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/365/GCF_000011365.1_ASM1136v1/GCF_000011365.1_ASM1136v1_genomic.gbff.gz | gunzip - > GCF_000011365.1.gbff
```

```bash
printf "GCA_002271865.1.fa\tOur_Alteromonas_MAG\nGCF_000011365.1.gbff\tGCF_000011365.1_Alpha_Outgroup\n" > genome-to-id-map.tsv
```

```bash
ls *.fa > fasta_files.txt
ls *.gbff > genbank_files.txt
```

```bash
GToTree -a alteromonas_refseq_accessions.txt -g genbank_files.txt -f fasta_files.txt -H Gammaproteobacteria -t -L Species,Strain -m genome-to-id-map.tsv -o Alt -P -j 4

### took 18 minutes on binder...
```

```bash
GToTree -a alteromonas_refseq_accessions.txt -g genbank_files.txt -f fasta_files.txt -H Gammaproteobacteria -t -L Species,Strain -m genome-to-id-map.tsv -o Alt -P -j 4

### took 6.5 minutes on my laptop 
```



## Tree of Life example

## Mixed-model example

## Tree visualization and manipulation

# The level of diversity of our arbitrary window matters

