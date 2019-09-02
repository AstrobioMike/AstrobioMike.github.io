---
layout: main
title: Accessing data from NCBI with EDirect
categories: [unix, tutorial]
permalink: /unix/ncbi_eutils
---

{% include _unix_downloading_from_ncbi_toc.html %}

{% include _side_tab_unix.html %}

[NCBI](https://www.ncbi.nlm.nih.gov/){:target="_blank"} is definitely pretty awesome. But sometimes it can be a little tricky to figure out how to download the data we want – particularly when it's a lot of things and we want and/or need to do it at the command-line rather than at the site. There are some convenient tools available that may help in some situations depending on our needs. For searching by taxonomy and downloading genomes (assemblies), [@kaiblin](https://twitter.com/kaiblin){:target="_blank"} and some others have put together [a great tool](https://github.com/kblin/ncbi-genome-download){:target="_blank"} called `ncbi-genome-download`, and they have [one for downloading individual sequences by accession](https://github.com/kblin/ncbi-acc-download){:target="_blank"} called `ncbi-acc-download`. I also put together a small program specifically for downloading assemblies by accession when making [GToTree](https://github.com/AstrobioMike/GToTree/wiki/what-is-gtotree%3F){:target="_blank"}, and I found it useful in general so it's also now part of my [BioInf Tools](https://github.com/AstrobioMike/bioinf_tools#bioinformatics-tools-bit){:target="_blank"} as `bit-dl-ncbi-assemblies`. 

But there have been times when I needed more than these could offer, which required me spending a decent amount of time getting used to using NCBI's
[EDirect tools](https://www.ncbi.nlm.nih.gov/books/NBK179288/){:target="_blank"}. **EDirect is a set of tools NCBI provides to enable accessing the vast amount of information stored at NCBI from the command line.** Learning to use it at all was painful at times for me – and still is when I'm trying to figure out new stuff with it. But like a lot of things, once I had a little familiarity with it, I started to appreciate how useful and powerful it can be.

There is a lot of info on this at [the main NCBI page for this here](https://www.ncbi.nlm.nih.gov/books/NBK179288/){:target="_blank"}, and some more [here](https://dataguide.nlm.nih.gov/edirect/overview.html){:target="_blank"}. Those are where I've learned from (along with a lot of trial and error). But I still have trouble finding examples when I need them (and I always need them when using this). And I often end up digging through several of my old log files to find things, so...

This page holds some of the ways I've used EDirect, both to serve as a handy archive for myself, and to hopefully help others 🙂 It won't be as comprehensive as most other things on this site, as it's extremely expansive and it's still nowhere near intuitive for me 🤷‍♂️ But these examples may do what is needed, and if not they at least might provide good starting points for building the code that will do what is needed. 

<hr style="height:15px; visibility:hidden;" />

---
---
<br>

<h1>Using NCBI's EDirect</h1>

For a while, this was also kind of a huge pain to install on some systems. But thanks to the gloriousness of [conda](https://conda.io/docs/){:target="_blank"}, that's no longer the case 🙂

<hr style="height:15px; visibility:hidden;" />
## Conda install

```bash
conda install -y -c conda-forge -c bioconda -c defaults entrez-direct
```

<hr style="height:15px; visibility:hidden;" />
## Accessing genome assemblies and info

* **All *Alteromonas* assembly accessions** (for instance to input into [GToTree](https://github.com/AstrobioMike/GToTree/wiki/what-is-gtotree%3F){:target="_blank"} like the [example here](https://github.com/AstrobioMike/GToTree/wiki/example-usage#alteromonas-example){:target="_blank"}!)

```bash
esearch -db assembly -query '"Alteromonas"[Organism] AND latest[filter] AND \
        (all[filter] NOT anomalous[filter] AND all[filter] NOT "derived from \
        surveillance project"[filter])' | esummary | xtract -pattern \
        DocumentSummary -element AssemblyAccession > Alteromonas-assembly-accs.txt
```

>**NOTE:** We can build the search string at the NCBI website and copy and paste it from there (there's a little "Search Details" box at the right side of the search page that adds in things as we modify our search on the site). Doing the search and download with EDirect still helps a lot because it lets us: automate and document things well; download directly to a server rather than our local computer; pull more specific information than we can on the site; and more 🙂

<hr style="height:10px; visibility:hidden;" />

* **All RefSeq Bacteria assembly accessions, taxids, assembly status, number of contigs, L50, N50, and total assembly length** (took ~5 minutes to get 166,566 records as accessed on 1-Sep-2019)

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

* **Downloading genomes by accession**  

I'd typically do this part with one of the tools listed above in the intro, after getting the accessions like the examples above, but here's an example pulling the sequence data for a RefSeq accession:

```bash
esearch -db assembly -query GCF_006538345.1 | elink -target nucleotide -name \
        assembly_nuccore_refseq | efetch -format fasta > GCF_006538345.1.fa
```

And one for a GenBank accession:

```bash
esearch -db assembly -query GCA_006538345.1 | elink -target nucleotide -name \
        assembly_nuccore_insdc | efetch -format fasta > GCA_006538345.1.fa
```

Note the change in the `-name` parameter between those two. "*assembly_nuccore_insdc*" is for GenBank, while "*assembly_nuccore_refseq*" is for RefSeq. Also note that I have no idea what the underlying infrastructure is here, and have to trial-and-error things whenever I'm trying to find something for the first time. Two places to look are with the `einfo -dbs` command and at [this site here](https://www.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html){:target="_blank"}. Wish I could be more helpful than that, believe me!

Example downloading with [BioInf Tools](https://github.com/AstrobioMike/bioinf_tools#bioinformatics-tools-bit){:target="_blank"}:

```bash
echo -e "GCF_006538345.1\nGCF_006538325.1\nGCF_006538305.1" > assembly-accs.txt

bit-dl-ncbi-assemblies -w assembly-accs.txt -f fasta
```

<hr style="height:15px; visibility:hidden;" />

---
<br>

## Accessing proteins
### Sequences
* **Single protein by accesssion**

```bash
efetch -db protein -format fasta -id ABA21534.1
```
<hr style="height:15px; visibility:hidden;" />

* **Many**  
We can't provide an input file to the **`efetch`** command, but we can to **`epost`** first. Assuming we have our target accessions in a single-column file, this can be done like so:

```bash
echo -e "ABA21534.1\nWP_013322114.1\nWP_015207051.1" > accs.txt

epost -input accs.txt -db protein | efetch -format fasta
```

>**NOTE:** If doing this with tens of thousands of targets, we need to split things up. There is a section at the end of this page covering this. 

<hr style="height:10px; visibility:hidden;" />

* **Nucleotide coding sequence for protein from accession**

```bash
efetch -db protein -format fasta_cds_na -id ABA21534.1
```

<hr style="height:15px; visibility:hidden;" />

* **Protein sequences from assembly accession**

```bash
esearch -db assembly -query GCA_006538345.1 | elink -target nuccore -name \
        assembly_nuccore_insdc | elink -target protein | efetch -format fasta \
        > GCA_006538345.1.faa
```

<hr style="height:20px; visibility:hidden;" />

---
<br> 
### GIs from accessions
I needed to do this to block them from a blast search as it can take a list of GIs to block, but not a list of accessions.

**Single**

```bash
esearch -db protein -query AEE52072.1 | esummary | xtract -pattern DocumentSummary \
        -element Gi
```

<hr style="height:10px; visibility:hidden;" />

**Many**  
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

* **Single**  

```bash
efetch -db protein -format gb -mode xml -id ABA21534.1 | xtract -pattern GBSeq \
       -element GBSeq_accession-version -block GBQualifier -if GBQualifier_name \
       -equals product -element GBQualifier_value
```

<hr style="height:10px; visibility:hidden;" />

* **Many**  

```bash
epost -input accs.txt -db protein | efetch -format gb -mode xml | xtract -pattern \
      GBSeq -element GBSeq_accession-version -block GBQualifier -if GBQualifier_name \
      -equals product -element GBQualifier_value
```

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

<hr style="height:25px; visibility:hidden;" />

---
<br>
# Breaking up lots of queries

The other thing we have to address is that the [Entrez site notes](https://www.ncbi.nlm.nih.gov/books/NBK179288/#chapter6.Automation){:target="_blank"} that we shouldn't submit queries than the site can handle for the database we are targeting. This doesn't matter when we're doing a single search that is grabbing lots of things (like the example above to get all RefSeq bacterial accessions), but rather if we are submitting lots of things that need to be queried individually. 

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

Again, I realize this isn't comprehensive or a fundamental understanding of any of this, and I apologize for that. `xtract` on its own is a pretty expansive language to learn, on top of needing to know the structure of how NCBI links all of its stuff together. It's a lot! But hopefully these examples being stored here will help more than just me 🙂