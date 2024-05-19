---
layout: main
title: Accessing data from NCBI with EDirect
categories: [unix, tutorial]
tags: [unix,bash,bioinformatics,tools,ncbi,download,data,tutorial]
permalink: /unix/ncbi_eutils
---

{% include _unix_downloading_from_ncbi_toc.html %}

{% include _side_tab_unix.html %}

[NCBI](https://www.ncbi.nlm.nih.gov/){:target="_blank"} is definitely pretty awesome. But sometimes it can be a little tricky to figure out how to download the data we want – particularly when it's a lot of things and we want and/or need to do it at the command-line rather than at the site. There are some convenient tools available that may help in some situations depending on our needs. For searching by taxonomy and downloading genomes (assemblies), [@kaiblin](https://twitter.com/kaiblin){:target="_blank"} and some others have put together [a great tool](https://github.com/kblin/ncbi-genome-download){:target="_blank"} called `ncbi-genome-download`, and they have [one for downloading individual sequences by accession](https://github.com/kblin/ncbi-acc-download){:target="_blank"} called `ncbi-acc-download`. I also put together a small program specifically for downloading assemblies by accession when making [GToTree](https://github.com/AstrobioMike/GToTree/wiki/what-is-gtotree%3F){:target="_blank"}, and I found it useful in general so it's also now part of my [Bioinf Tools](https://github.com/AstrobioMike/bioinf_tools#bioinformatics-tools-bit){:target="_blank"} as `bit-dl-ncbi-assemblies`. 

But there have been times when I needed more than these could offer, which required me spending a decent amount of time getting used to using NCBI's
[EDirect tools](https://www.ncbi.nlm.nih.gov/books/NBK179288/){:target="_blank"}. **EDirect is a set of tools NCBI provides to enable accessing the vast amount of information stored at NCBI from the command line.** Learning to use it at all was painful at times for me – and still is when I'm trying to figure out new stuff with it. But like a lot of things, once I had a little familiarity with it, I started to appreciate how useful and powerful it can be. 

<hr style="height:15px; visibility:hidden;" />
<center><h3>Why use EDirect?</h3></center>

Being able to access data and info from NCBI at the command line can allow us to: automate and document things well (we can give the exact command used to retrieve information and the date it was executed, rather than "pulled from NCBI"); download directly to a server rather than our local computer; pull more specific information than we can on the site; and more 🙂

There is a lot of info on this at [the main NCBI EDirect page](https://www.ncbi.nlm.nih.gov/books/NBK179288/){:target="_blank"}, and some more [here](https://dataguide.nlm.nih.gov/edirect/overview.html){:target="_blank"}. Those are where I've learned from (along with *a lot* of trial and error). But I still have trouble finding examples when I need them (and I always need them when using this). And I often end up digging through several of my old log files to find things, so...

This page holds some of the ways I've used EDirect, both to serve as a handy archive for myself, and to hopefully help others 🙂 It won't be as comprehensive as most other things on this site, as it's extremely expansive and it's still nowhere near intuitive for me 🤷‍♂️ But these examples may do what is needed, and if not they at least might provide good starting points for building the code that will do what is needed. 

>**NOTE:** This stuff can look messy, mostly because it is; it scares me too. **This page assumes you already have some familiarity with working at the command line.** If you don't yet, then definitely consider running through the [Unix crash course](/unix/unix-intro){:target="_blank"} first 🙂

<hr style="height:15px; visibility:hidden;" />

---
---
<br>

<h1>Using NCBI's EDirect</h1>

For a while, this was also kind of a huge pain to install on some systems. But thanks to the gloriousness of [conda](/unix/conda-intro){:target="_blank"}, that's no longer the case 🙂

<hr style="height:15px; visibility:hidden;" />
## Conda install

```bash
conda create -y -n edirect -c conda-forge -c bioconda -c defaults -c defaults entrez-direct
conda activate edirect
```

<hr style="height:25px; visibility:hidden;" />

---
<br>
## Accessing genome assemblies and info

**Getting all *Alteromonas* assembly accessions** (e.g. to input into [GToTree](https://github.com/AstrobioMike/GToTree/wiki/what-is-gtotree%3F){:target="_blank"} like the [example here](https://github.com/AstrobioMike/GToTree/wiki/example-usage#alteromonas-example){:target="_blank"}!)

```bash
esearch -db assembly -query '"Alteromonas"[Organism] AND latest[filter] AND \
        (all[filter] NOT anomalous[filter] AND all[filter] NOT "derived from \
        surveillance project"[filter])' | esummary | xtract -pattern \
        DocumentSummary -element AssemblyAccession > Alteromonas-assembly-accs.txt
```

>**NOTE:** We can build the search string at the NCBI website and copy and paste it from there (there's a little "Search Details" box at the right side of the search page that adds in things as we modify our search on the site).

<hr style="height:10px; visibility:hidden;" />

**Getting all RefSeq Bacteria assembly accessions, taxids, assembly status, number of contigs, L50, N50, and total assembly length** (took ~5 minutes to get 166,566 records as accessed on 1-Sep-2019)

```bash
esearch -db assembly -query '"Bacteria"[Organism] AND "latest refseq"[filter] AND \
        (all[filter] NOT anomalous[filter] AND all[filter] NOT "derived from \
        surveillance project"[filter])' | esummary | xtract -pattern DocumentSummary \
        -def "NA" -element AssemblyAccession,Taxid,assembly-status -block Stat \
        -if Stat@category -equals contig_count -or Stat@category -equals contig_l50 \
        -or Stat@category -equals contig_n50 -or Stat@category -equals total_length \
        -sep ":" -def "NA" -element Stat@category,Stat \
        > All-bacteria-refseq-complete-assembly-info.tsv
```

<hr style="height:15px; visibility:hidden;" />

**Downloading genomes by accession**  

I typically do this part using `bit-dl-ncbi-assemblies` from my [Bioinf Tools](https://github.com/AstrobioMike/bioinf_tools#bioinformatics-tools-bit){:target="_blank"} after getting the accessions like in the examples above. Here's what that would look like following the [*Alteromonas* example above](/unix/ncbi_eutils#accessing-genome-assemblies-and-info) to download those 160 assemblies (here in fasta format; took ~1 minute):

```bash
bit-dl-ncbi-assemblies -w Alteromonas-assembly-accs.txt -f fasta -j 10
```

But here's an example using EDirect to pull the sequence data for a RefSeq accession:

```bash
esearch -db assembly -query GCF_006538345.1 | elink -target nucleotide -name \
        assembly_nuccore_refseq | efetch -format fasta > GCF_006538345.1.fa
```

And one for a GenBank accession:

```bash
esearch -db assembly -query GCA_006538345.1 | elink -target nucleotide -name \
        assembly_nuccore_insdc | efetch -format fasta > GCA_006538345.1.fa
```

Note the change in the `-name` parameter between those two. "*assembly_nuccore_insdc*" is for GenBank, while "*assembly_nuccore_refseq*" is for RefSeq. Also note that I have no idea what the underlying infrastructure is here, and have to trial-and-error things whenever I'm trying to find something for the first time. Two places to look are with the `einfo -dbs` command and at [this site here](https://www.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html){:target="_blank"}.

<hr style="height:15px; visibility:hidden;" />

---
<br>

### Protein annotations from assembly accession

```bash
esearch -db assembly -query GCA_006538345.1 | elink -target nuccore -name \
        assembly_nuccore_insdc | elink -target protein | efetch -format gb \
        -mode xml | xtract -pattern GBSeq -element GBSeq_accession-version \
        -block GBQualifier -if GBQualifier_name -equals product \
        -element GBQualifier_value > GCA_006538345.1-annotations.tsv
```

<hr style="height:15px; visibility:hidden;" />

---
<br>

## Accessing individual genes/proteins
### Sequences and accessions
**Getting amino acid sequences based on protein-name text search**

```bash
esearch -db protein -query '"nosZ"[Protein name]' | efetch -format fasta > nosZ.faa
```

**Getting only unique sequences from the [Identical Protein Groups](https://www.ncbi.nlm.nih.gov/ipg/docs/about/){:target="_blank"} database**

```bash
esearch -db IPG -query '"nosZ"[Protein name]' | efetch -format fasta > nosZ-IPG.faa
```

**Getting nucleotide coding sequences for proteins based on protein-name text search**

```bash
esearch -db protein -query '"nosZ"[Protein name]' | efetch -format fasta_cds_na > nosZ.fa
```

**Getting protein accessions based on protein-name text search**

```bash
esearch -db protein -query '"nosZ"[Protein name]' | esummary | xtract -pattern \
        DocumentSummary -element AccessionVersion > nosZ-accs.txt
```

**Getting a single protein sequence by accession**

```bash
efetch -db protein -format fasta -id ABA21534.1
```

**Getting a single nucleotide sequence by accession** 

```bash
efetch -db nucleotide -format fasta -id NC_006572.1
```


<hr style="height:15px; visibility:hidden;" />

**Getting many protein sequences by accession**  
We can't provide an input file to the **`efetch`** command, but we can to **`epost`** first. Assuming we have our target accessions in a single-column file, this can be done like so:

```bash
echo -e "ABA21534.1\nWP_013322114.1\nWP_015207051.1" > accs.txt

epost -input accs.txt -db protein | efetch -format fasta
```

>**NOTE:** If doing this with tens of thousands of targets, we need to split things up. There is a section at the end of this page covering this. 

<hr style="height:10px; visibility:hidden;" />

**Getting the nucleotide coding sequence for a protein from the protein accession**

```bash
efetch -db protein -format fasta_cds_na -id ABA21534.1
```

<hr style="height:15px; visibility:hidden;" />

**Getting protein sequences from an assembly (genome) accession**

```bash
esearch -db assembly -query GCA_006538345.1 | elink -target nuccore -name \
        assembly_nuccore_insdc | elink -target protein | efetch -format fasta \
        > GCA_006538345.1.faa
```

<hr style="height:20px; visibility:hidden;" />

---
<br> 
### GIs from accessions
NCBI has unique GI IDs for things in addition to their accessions. I needed GIs specifically once in order to block them from a command-line BLAST search as it can take them to block, but not a list of accessions.

**Getting a single GI from a protein accession**

```bash
esearch -db protein -query AEE52072.1 | esummary | xtract -pattern DocumentSummary \
        -element Gi
```

<hr style="height:10px; visibility:hidden;" />

**Getting many GIs from protein accessions**  
We can use **`epost`** similar to above:

```bash
epost -input accs.txt -db protein | esummary | xtract -pattern DocumentSummary \
      -element Gi
```

<hr style="height:20px; visibility:hidden;" />

<center><b>IMPORTANT</b></center>

**`epost`** will not return things in the same order you input them. If you need to preserve the order, you should consider returning the accession as well, e.g.:

```bash
epost -input accs.txt -db protein | esummary | xtract -pattern DocumentSummary \
      -element AccessionVersion,Gi
```

And then we can sort things based on the accession to make the input and outputs retain order. 

<hr style="height:15px; visibility:hidden;" />

---
<br> 
### Protein annotations from accessions

**Getting a single protein annotation from its accession**  

```bash
efetch -db protein -format gb -mode xml -id ABA21534.1 | xtract -pattern GBSeq \
       -element GBSeq_accession-version -block GBQualifier -if GBQualifier_name \
       -equals product -element GBQualifier_value
```

<hr style="height:10px; visibility:hidden;" />

**Getting many protein annotations from their accessions**  

```bash
epost -input accs.txt -db protein | efetch -format gb -mode xml | xtract -pattern \
      GBSeq -element GBSeq_accession-version -block GBQualifier -if GBQualifier_name \
      -equals product -element GBQualifier_value
```

<hr style="height:25px; visibility:hidden;" />

---
<br>

## Accessing taxonomy info from protein accessions

### Source-organism taxid and lineage from protein accession

**Running for a single protein accession**

```bash
efetch -db protein -id WP_041587175.1 -format xml | xtract -pattern Seq-entry \
       -def "NA" -element Textseq-id_accession,Textseq-id_version \
       -block BioSource_org -def "NA" -element Object-id_id,OrgName_lineage
```

**Running for many protein accessions**

```bash
epost -input accs.txt -db protein | efetch -format xml | xtract -pattern Seq-entry \
      -def "NA" -element Textseq-id_accession,Textseq-id_version \
      -block BioSource_org -def "NA" -element Object-id_id,OrgName_lineage
```

<hr style="height:25px; visibility:hidden;" />

---
<br>

# Breaking up lots of queries

The other thing we have to address is that the [Entrez site notes](https://www.ncbi.nlm.nih.gov/books/NBK179288/#chapter6.Automation){:target="_blank"} that we shouldn't submit more queries than the site can handle for the database we are targeting. This doesn't matter when we're doing a single search that is grabbing lots of things (like the example above to get all RefSeq bacterial accessions), but rather if we are submitting lots of things that need to be queried individually. 

I needed to do this at some point, but I can't remember exactly for what. So this example will be a little contrived, but will show how I did it. Let's first grab all the protein accessions from an *Alteromonas* with assembly accession GCA_006538345.1:

```bash
esearch -db assembly -query GCA_006538345.1 | elink -target nuccore -name \
        assembly_nuccore_insdc | elink -target protein | esummary | xtract -pattern \
        DocumentSummary -element AccessionVersion > GCA_006538345.1-prot-accs.txt
```

This is 4,135 accessions. And here's a nice little bash function to help do this that comes from the [NCBI EDirect page](https://www.ncbi.nlm.nih.gov/books/NBK179288/){:target="_blank"}:

```bash
JoinIntoGroupsOf() {
  xargs -n "$@" echo |
  sed 's/ /,/g'
}
```

This will take chunks of our accessions and provide them separated by commas to our input. It's also a good idea to provide our email address when we are submitting a lot of things. That way if we cause a problem on their end, they can contact us to tell us and we won't keep breaking the system until we get banned 😬

```bash
cat GCA_006538345.1-prot-accs.txt | JoinIntoGroupsOf 200 | xargs -n 1 sh -c \
    'efetch -email "MyEmail@gmail.com" -db protein -id $0 -format fasta' \
    > GCA_006538345.1.faa
``` 

That was submitting chunks of 200 at a time and took about 20 seconds for me. Looking at [this page here of the Entrez database links](https://www.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html){:target="_blank"} we can see that *nuccore_protein* says maximum items processed of 5,000. So I think we can set this value as high as that (which in this case is higher than our total targets, but this is how we'd do it if we had more):

```bash
cat GCA_006538345.1-prot-accs.txt | JoinIntoGroupsOf 5000 | xargs -n 1 sh -c \
    'efetch -email "MyEmail@gmail.com" -db protein -id $0 -format fasta' \
    > GCA_006538345.1.faa
``` 

And that took 6 seconds. 

<hr style="height:25px; visibility:hidden;" />

---
---
<br>

Again, sorry this isn't a comprehensive or a fundamental explanation of any of this. I'm just not at that level with this particular skillset. `xtract` on its own is a pretty expansive language to learn, on top of needing to know the structure of how NCBI links all of its stuff together. It's a lot! I'd like to know it better someday, and then would love to put together a page for it 🤞 

For now, hopefully these examples being stored here will help more than just me 🙂
