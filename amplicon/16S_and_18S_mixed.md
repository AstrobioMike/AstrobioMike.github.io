---
layout: main
title: 16S and 18S mixed together
categories: [amplicon, tutorial]
tags: [amplicon,16S,18S,metabarcoding,dada2]
permalink: amplicon/16S_and_18S_mixed
---

{% include _amplicon_16S_and_18S_toc.html %}

{% include _side_tab_amplicon.html %}

So if you're using primers [like these that capture 16S and 18S pretty well](https://onlinelibrary.wiley.com/doi/abs/10.1111/1462-2920.13023){:target="_blank"} and then writing cool papers [like this one](https://www.nature.com/articles/nmicrobiol20165){:target="_blank"}, then you might be wondering if/how you can do this while taking advantage of all the awesomeness that is DADA2. The 18S amplified fragments are typically too long to have any overlap occur with the 2x250bp sequencing that is often performed with these "V4V5" primers, so [one method that has been employed](https://www.nature.com/articles/nmicrobiol20165){:target="_blank"} was to process the reads that didn't successfully merge as the 18S reads. Starting from within DADA2, I considered playing with the `mergePairs()` step and messing with the "justConcatenate" flag. The problem is that we also need to account for those that fail the merge just because of poor quality or whatever other reason they don't merge (rather than it just being because they didn't overlap at all). 

So my good buddy [Josh](https://twitter.com/KlingJoshua){:target="_blank"} and I spent a few hours trying to figure out how we could filter the merged objects to separate out which failed-to-merge sequences failed because they were likely not overlapping (and therefore anticipated to be 18S which we wanted to find) from those that just failed to merge because of sequencing error or other black magic failures (reads we wanted to remove and/or consider likely to be 16S). And, not surprisingly, one can end up in quite the rabbit hole about this looking at nmismatches and nmatches and such in the `mergePairs()` resulting objects. But ultimately it seemed like even if we found optimal parameters for one dataset, they might be different for the next one. 

So I decided to bail on that approach and wanted to see if there was an efficient way (i.e. could run in a reasonable amount of time even on my laptop) to parse the reads into 16S and 18S *before* putting them into DADA2. Here's what I came up with that so far has worked well:

||Step|What we're doing|
|:--:|:--------:|----------|
|1|`Magic-BLAST`|blast all reads to the [PR2 database](https://github.com/pr2database/pr2database#pr2database-){:target="_blank"} (here done with v4.14.0)|
|2|`Filtering Magic-BLAST output`|based on % ID and % of query sequence aligned (of both reads)|
|3|`Splitting 16S/18S reads`|based on the Magic-BLAST filtering|
|4|`Processing both in DADA2`|processing independently and merging at the end|

>**NOTE:** Big thanks to Soluna Sales ([@Yo_mis_mamente](https://twitter.com/Yo_mis_mamente){:taret="_blank"}) for writing in and helping to improve this page 🙂

<hr>
<br>

# 16S/18S example data
For an example of this process, we're going to work with a couple of samples from [the paper](https://www.nature.com/articles/nmicrobiol20165){:target="_blank"} I mentioned above by [David Needham](https://twitter.com/animalkewls){:target="_blank"} and [Jed Fuhrman](https://twitter.com/JedFuhrman){:target="_blank"}. If you'd like to follow along, you can download a small directory containing the data, the [PR2](https://figshare.com/articles/PR2_rRNA_gene_database/3803709){:target="_blank"} database tsv file initially downloaded and used here, and the code that follows (including custom scripts for formatting the PR2 database and splitting the suspected 16S and 18S reads into separate fastq files) with these commands:

```bash
cd ~
curl -L -o dada2_16S_18S_ex.tar.gz https://ndownloader.figshare.com/files/29063763
tar -xzvf dada2_16S_18S_ex.tar.gz
rm dada2_16S_18S_ex.tar.gz
cd dada2_16S_18S_ex
```

The two samples in this working directory are [ERR1018543](https://www.ncbi.nlm.nih.gov/sra/ERR1018543){:target="_blank"} and [ERR10185469](https://www.ncbi.nlm.nih.gov/sra/ERR1018546){:target="_blank"} which were both filtered marine water, size fraction 1µm–80µm, sequenced 2x300 with primers targeting the V4V5 region of the 16S rRNA gene which also capture 18S (515f: 5'-GTGCCAGCMGCCGCGGTAA-3' to 926r: 5'CCGYCAATTYMTTTRAGTTT-3'; [Parada et al. 2015](https://onlinelibrary.wiley.com/doi/abs/10.1111/1462-2920.13023){:target="_blank"}). 

---
<br>


# Conda environment
Here's how we can create a [conda](/unix/conda-intro){:target="_blank"} environment for the work done on this page if wanted. 

```bash
conda create -n hb-16S-18S-example -c conda-forge -c bioconda -c defaults -c astrobiomike \
             python=3 magicblast=1.5.0 cutadapt=3.4 bit=1.8.33 r-base rstudio \
             r-biocmanager r-tidyverse=1.3.1 r-readxl=1.3.1 bioconductor-dada2=1.20.0 \
             bioconductor-decipher=2.20.0

conda activate hb-16S-18S-example
```

<br>

---
<br>

# Getting and formatting the PR2 database
The great folks working on the [PR2 database](https://github.com/pr2database/pr2database#pr2database-){:target="_blank"} have lots of formats available for us. The general table is in xlsx format though ("pr2_version_4.14.0_merged.xlsx"), and I was having trouble manipulating it at the command line due to some cells having line-breaks within them (e.g., entry AB695498.1.1724_U). So we are going to read it into R to parse it down to just 18S seqs and to create the fasta file we need for generating a blast database. 

Downloading at the command-line (it also came with the above working directory download if you pulled that:

```bash
curl -LO https://github.com/pr2database/pr2database/releases/download/v4.14.0/pr2_version_4.14.0_merged.xlsx
```

And now we're going to manipulate it in R. If you created and entered the conda environment above, you can enter RStudio by typing `rstudio` at the command line. 


```R
## in R ##
library(readxl)
tab <- read_xlsx("pr2_version_4.14.0_merged.xlsx")
sub_tab <- tab %>% filter(gene == "18S_rRNA")
# making fasta
headers <- sub_tab %>% pull(pr2_accession)
seqs <- sub_tab %>% pull(sequence)

fasta <- c(rbind(paste0(">", headers), seqs))

# writing out
write(fasta, "pr2_version_4.14.0_euk_18S_only.fasta")
```

And that's our 18S reference fasta file we're going to blast the reads against.

---
<br>

# Trimming our primers 
It's best to cut off the highly conserved primer regions which allowed the amplification of all 3 domains in the first place, as those bases certainly won't help us distinguish things. Here we are using [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html){:target="_blank"}.

```bash
# making samples list
ls *_1.fastq.gz | cut -f1 -d "_" > samples.txt

# trimming primers
for sample in $(cat samples.txt)
do

    printf "\n  Trimming sample: ${sample}\n"
    cutadapt -g GTGCCAGCMGCCGCGGTAA -G CCGYCAATTYMTTTRAGTTT \
             -o ${sample}_1_trimmed.fastq.gz -p ${sample}_2_trimmed.fastq.gz \
             --discard-untrimmed ${sample}_1.fastq.gz ${sample}_2.fastq.gz \
             > ${sample}-cutadapt.log 2>&1

done
```

<br>

---
<br>

# Magic-BLAST
[NCBI's Magic-BLAST](https://ncbi.github.io/magicblast/){:target="_blank"} is a tool based on general [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi){:target="_blank"} principles but built to: 1) deal with high-throughput data (like Illumina reads); 2) consider paired-reads; and 3) work with fastq files. So yeah, I hadn't heard of this before looking for this particular problem, but it's perfect for it 🙂

First we need to make our blast database:

```bash
makeblastdb -in pr2_version_4.14.0_euk_18S_only.fasta \
            -dbtype nucl -parse_seqids -out pr2-magicblast-db
```

And now we can run the blast of each individual sample's paired-reads against our 18S database:

```bash
for sample in $(cat samples.txt)
do

    printf "\n  Doing sample: ${sample}\n"
  
    magicblast -db pr2-magicblast-db -query "${sample}"_1_trimmed.fastq.gz \
             -query_mate "${sample}"_2_trimmed.fastq.gz -infmt fastq \
             -out "${sample}"_mblast_out.txt -outfmt tabular \
             -num_threads 2 -splice F -no_unaligned
             
done
```

<br>

---
<br>

# Filtering Magic-BLAST output
In messing with how to filter this output to be most useful, I ended on requiring greater than 35% of the query to be aligned at greater than 90% ID (for both forward and reverse reads). If only one read passed these criteria, the fragment they originated from was considered not to be a successful hit to the 18S database. This seemed to work well on these samples, but if you're trying this on your own data, I encourage you to mess with the filtering parameters we use below to ensure you are happy with what you're getting.
There's no set values that are known to always be good here, these are just what I ended up using.

First, we are going to trim down the magicblast output tables, and add a new column that has the percent of the query-read that aligned (so we can use that column to filter based on the 35% I mentioned):

```bash
for sample in $(cat samples.txt)
do 

    printf "\n  Doing sample: ${sample}\n"
  
    cut -f 1,2,3,7,8,16 "${sample}"_mblast_out.txt | sed '1d' | sed '1d' | \
      sed 's/# Fields: //' | tr " " "_" | \
      awk -F $'\t' ' BEGIN { OFS=FS } NR==1 { $7="%_query_aln"; print $0 } NR>1 { print $0, ($5-$4)/$6*100 } ' \
      > "${sample}"_mblast_out_mod.txt
      
done
```

Now we are going filter the hits. As mentioned, as coded here, we are requiring that the percent identity was > 90 (column 3 in our tables), and that > 35% of the query aligned (column 7 in our tables). Then we are cutting down to just the first column, which are the read IDs, and removing any that only appear once. This is because the forward and reverse reads have the same ID in these files, so anything that only appears once means only one of the paired reads passed our criteria – and we want both.  

```bash
for sample in $(cat samples.txt)
do 

    printf "\n  Doing sample: ${sample}\n"
  
    awk ' $3 > 90 && $7 > 35 ' "${sample}"_mblast_out_mod.txt | \
        cut -f 1 | uniq -d > "${sample}"_18S_headers.txt
      
done
```

Now for each sample we have a file of headers for the reads that we believe are of 18S origin, e.g.:

```bash
head ERR1018543_18S_headers.txt | sed 's/^/# /'
# ERR1018543.1
# ERR1018543.21
# ERR1018543.28
# ERR1018543.31
# ERR1018543.41
# ERR1018543.43
# ERR1018543.58
# ERR1018543.63
# ERR1018543.64
# ERR1018543.77
```

<br>

---
<br>

# Splitting fastq files into 16S/18S
Next I wrote a little ad hoc python script to split the fastq files for each sample into those that hit our 18S database and those that did not. This script is included the downloaded working directory above, and is also available [here](https://figshare.com/articles/software/split_16S_18S_reads_py/15123234){:target="_blank"}. 

It takes as input the forward and reverse fastq files (gzipped required as written), and it needs the file of 18S-read headers we just made. It gives us back 4 fastq files: 16S and 18S forward reads, and 16S and 18S reverse reads. The help menu can be seen with `python split_16S_18S_reads.py -h`, but here is how I ran it looped for all samples:

```bash
for sample in $(cat samples.txt)
do 

    printf "\n  Doing sample: ${sample}\n"
  
    python split_16S_18S_reads.py -f ${sample}_1_trimmed.fastq.gz \
                                  -r ${sample}_2_trimmed.fastq.gz \
                                  -E ${sample}_18S_headers.txt
  
done
```

Now we can see the four output files for each sample and how they're labeled:

```bash
ls -t | head -n 8 | sed 's/^/# /'
# 16S_ERR1018546_2_trimmed.fastq.gz
# 18S_ERR1018546_2_trimmed.fastq.gz
# 16S_ERR1018546_1_trimmed.fastq.gz
# 18S_ERR1018546_1_trimmed.fastq.gz
# 16S_ERR1018543_2_trimmed.fastq.gz
# 18S_ERR1018543_2_trimmed.fastq.gz
# 16S_ERR1018543_1_trimmed.fastq.gz
# 18S_ERR1018543_1_trimmed.fastq.gz
```

<br>

---
<br>

# Processing both in R with DADA2
This part won't be heavily annotated, as these steps are already annotated and laid out in the [full DADA2 example workflow](/amplicon/dada2_workflow_ex) and this doing the same except for 16S and 18S separately now. But here is the R code (also in the downloaded working directory as "16S_18S_processing.R"). Don't forget, if you installed the conda environment above, you can access this RStudio by running `rstudio` at the command line. Isn't [conda](/unix/conda-intro){:target="_blank"} awesome?? 🙂

---
<br>

## Setting up R environment
```R
library(dada2)
library(tidyverse)
library(DECIPHER)

samples <- scan("samples.txt", what="character")
```

<br>

---
<br>

## 18S processing
```R
forward_18S_reads <- paste0("18S_", samples, "_1_trimmed.fastq.gz")
reverse_18S_reads <- paste0("18S_", samples, "_2_trimmed.fastq.gz")


filtered_forward_18S_reads <- paste0("18S_", samples, "_1_filtered.fastq.gz")
filtered_reverse_18S_reads <- paste0("18S_", samples, "_2_filtered.fastq.gz")

plotQualityProfile(forward_18S_reads) # median (green line) seems to cross Q30 around 230 bases
plotQualityProfile(reverse_18S_reads) # median crosses Q30 around 200


filtered_out_18S <- filterAndTrim(forward_18S_reads, filtered_forward_18S_reads,
                                  reverse_18S_reads, filtered_reverse_18S_reads,
                                  maxEE = c(2,2), rm.phix = TRUE, multithread = TRUE,
                                  truncLen = c(230,200))

plotQualityProfile(filtered_forward_18S_reads)
plotQualityProfile(filtered_reverse_18S_reads)


err_forward_18S_reads <- learnErrors(filtered_forward_18S_reads, multithread = TRUE)
err_reverse_18S_reads <- learnErrors(filtered_reverse_18S_reads, multithread = TRUE)

plotErrors(err_forward_18S_reads, nominalQ = TRUE)
plotErrors(err_reverse_18S_reads, nominalQ = TRUE)

derep_forward_18S <- derepFastq(filtered_forward_18S_reads, verbose = TRUE)
names(derep_forward_18S) <- samples
derep_reverse_18S <- derepFastq(filtered_reverse_18S_reads, verbose = TRUE)
names(derep_reverse_18S) <- samples

dada_forward_18S <- dada(derep_forward_18S, err = err_forward_18S_reads, multithread = TRUE)
dada_reverse_18S <- dada(derep_reverse_18S, err = err_reverse_18S_reads, multithread = TRUE)


  # justConcatenate=TRUE
merged_18S <- mergePairs(dada_forward_18S, derep_forward_18S, dada_reverse_18S,
                         derep_reverse_18S, justConcatenate = TRUE)

seqtab_18S <- makeSequenceTable(merged_18S)
dim(seqtab_18S)[2] # 290
sum(seqtab_18S) # 4580


seqtab.nochim_18S <- removeBimeraDenovo(seqtab_18S, method = "consensus",
                                        multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim_18S)[2] # 236

sum(seqtab.nochim_18S) / sum(seqtab_18S) # 0.97

## looking at counts throughout
getN <- function(x) sum(getUniques(x))

track_18S <- data.frame(row.names = samples, dada2_input = filtered_out_18S[,1],
                        filtered = filtered_out_18S[,2],
                        denoised = sapply(dada_forward_18S, getN),
                        merged = sapply(merged_18S, getN), table=rowSums(seqtab_18S),
                        no_chimeras = rowSums(seqtab.nochim_18S),
                        "perc_reads_survived" = round(rowSums(seqtab.nochim_18S) / filtered_out_18S[,1] * 100, 1))

track_18S
#            dada2_input filtered denoised merged table no_chimeras perc_reads_survived
# ERR1018543        6289     3714     3474   3423  3423        3326                52.9
# ERR1018546        2440     1430     1205   1157  1157        1115                45.7



### Taxonomy
## creating a DNAStringSet object from the ASVs
dna_18S <- DNAStringSet(getSequences(seqtab.nochim_18S))

## downloading silva DECIPHER database
download.file("http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", "SILVA_SSU_r138_2019.RData")
# loading ref taxonomy object
load("SILVA_SSU_r138_2019.RData")


tax_info_18S <- IdTaxa(dna_18S, trainingSet, strand = "both", processors = NULL)


## making and writing out standard output files:
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs_18S <- colnames(seqtab.nochim_18S)

asv_headers_18S <- vector(dim(seqtab.nochim_18S)[2], mode = "character")
for (i in 1:dim(seqtab.nochim_18S)[2]) {
  asv_headers_18S[i] <- paste(">ASV_18S", i, sep = "_")
}

# fasta:
asv_fasta_18S <- c(rbind(asv_headers_18S, asv_seqs_18S))
write(asv_fasta_18S, "18S_ASVs.fa")

# count table:
asv_tab_18S <- t(seqtab.nochim_18S) %>% data.frame
row.names(asv_tab_18S) <- sub(">", "", asv_headers_18S)
asv_tab_18S <- asv_tab_18S %>% rownames_to_column("ASV_ID")
write.table(asv_tab_18S, "18S_ASVs_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# tax table:

    # creating vector of desired ranks
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

    # creating table of taxonomy and setting any that are unclassified as "NA"
tax_tab_18S <- t(sapply(tax_info_18S, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
}))

colnames(tax_tab_18S) <- ranks
row.names(tax_tab_18S) <- NULL
tax_tab_18S <- data.frame("ASV_ID" = sub(">", "", asv_headers_18S), tax_tab_18S, check.names = FALSE)

write.table(tax_tab_18S, "18S_ASVs_taxonomy.tsv", sep = "\t", quote = F, row.names = FALSE)

## saving the seqtab.nochim_18S object
saveRDS(seqtab.nochim_18S, "seqtab.nochim_18S.rds")
```

<br>

---
<br>

## 16S processing
```R
forward_16S_reads <- paste0("16S_", samples, "_1_trimmed.fastq.gz")
reverse_16S_reads <- paste0("16S_", samples, "_2_trimmed.fastq.gz")

filtered_forward_16S_reads <- paste0("16S_", samples, "_1_filtered.fastq.gz")
filtered_reverse_16S_reads <- paste0("16S_", samples, "_2_filtered.fastq.gz")

plotQualityProfile(forward_16S_reads) # median drops below Q30 around 260
  # the primers span 515-926, we cut off about 40 bps when removing the primers, so
  # our target amplicon now is about 370
  # with 260 from forward, the reverse would need to be 110 minimum to reach
plotQualityProfile(reverse_16S_reads) # median drops below Q30 around 200

  # when doing the trimming step, it's important to make sure we aren't trimming them
  # so short that they cannot overlap, which would cause problems when we try to merge later
  # trimming the forward to 260 and reverse to 200 would leave us with around 90 bps overlap

filtered_out_16S <- filterAndTrim(forward_16S_reads, filtered_forward_16S_reads,
                                  reverse_16S_reads, filtered_reverse_16S_reads,
                                  maxEE=c(2,2), rm.phix=TRUE, multithread=TRUE,
                                  truncLen=c(260,200))

plotQualityProfile(filtered_forward_16S_reads)
plotQualityProfile(filtered_reverse_16S_reads)

err_forward_16S_reads <- learnErrors(filtered_forward_16S_reads, multithread = TRUE)
err_reverse_16S_reads <- learnErrors(filtered_reverse_16S_reads, multithread = TRUE)

plotErrors(err_forward_16S_reads, nominalQ = TRUE)
plotErrors(err_reverse_16S_reads, nominalQ = TRUE)

derep_forward_16S <- derepFastq(filtered_forward_16S_reads, verbose = TRUE)
names(derep_forward_16S) <- samples
derep_reverse_16S <- derepFastq(filtered_reverse_16S_reads, verbose = TRUE)
names(derep_reverse_16S) <- samples

dada_forward_16S <- dada(derep_forward_16S, err = err_forward_16S_reads, multithread = TRUE)
dada_reverse_16S <- dada(derep_reverse_16S, err = err_reverse_16S_reads, multithread = TRUE)

  # doing a temp merge without changing the minimum overlap to get a look
  # at the distribution of overlap values
temp_merged_16S <- mergePairs(dada_forward_16S, derep_forward_16S,
                              dada_reverse_16S, derep_reverse_16S)

quantile(temp_merged_16S[[1]]$nmatch, probs=seq(0,1,0.05))
#   0%   5%  10%  15%  20%  25%  30%  35%  40%  45%  50%  55%  60%  65%  70%  75%  80%  85%  90%  95% 100% 
#   36   83   84   85   86   86   86   86   87   87   88   90   90   90   91   91   91   91   92   93  107 
    
    # okay, going to use 80 as min overlap, as that captures >95% of the sequences in there

rm(temp_merged_16S)
merged_16S <- mergePairs(dada_forward_16S, derep_forward_16S, dada_reverse_16S,
                         derep_reverse_16S, minOverlap = 80)


seqtab_16S <- makeSequenceTable(merged_16S)
dim(seqtab_16S)[2] # 624
sum(seqtab_16S) # 39286


seqtab.nochim_16S <- removeBimeraDenovo(seqtab_16S, method = "consensus",
                                        multithread = TRUE, verbose = TRUE)

dim(seqtab.nochim_16S)[2] # 510

sum(seqtab.nochim_16S) / sum(seqtab_16S) # 0.95


## looking at counts throughout
getN <- function(x) sum(getUniques(x))

track_16S <- data.frame(row.names = samples, dada2_input = filtered_out_16S[,1],
                        filtered = filtered_out_16S[,2],
                        denoised = sapply(dada_forward_16S, getN),
                        merged = sapply(merged_16S, getN), table=rowSums(seqtab_16S),
                        no_chimeras = rowSums(seqtab.nochim_16S),
                        "perc_reads_survived" = round(rowSums(seqtab.nochim_16S) / filtered_out_16S[,1] * 100, 1))

track_16S
#            dada2_input filtered denoised merged table no_chimeras perc_reads_survived
# ERR1018543       16033    10818    10401   8101  8101        8042                50.2
# ERR1018546       18351    13407    13019  11725 11725       11582                63.1

### Taxonomy
## creating a DNAStringSet object from the ASVs
dna_16S <- DNAStringSet(getSequences(seqtab.nochim_16S))

# loading ref taxonomy object we downloaded above
load("SILVA_SSU_r138_2019.RData")

tax_info_16S <- IdTaxa(dna_16S, trainingSet, strand = "both", processors = NULL)


### Making and writing out standard output files:
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs_16S <- colnames(seqtab.nochim_16S)

asv_headers_16S <- vector(dim(seqtab.nochim_16S)[2], mode = "character")
for (i in 1:dim(seqtab.nochim_16S)[2]) {
  asv_headers_16S[i] <- paste(">ASV_16S", i, sep = "_")
}

# fasta:
asv_fasta_16S <- c(rbind(asv_headers_16S, asv_seqs_16S))
write(asv_fasta_16S, "16S_ASVs.fa")

# count table:
asv_tab_16S <- t(seqtab.nochim_16S) %>% data.frame
row.names(asv_tab_16S) <- sub(">", "", asv_headers_16S)
asv_tab_16S <- asv_tab_16S %>% rownames_to_column("ASV_ID")
write.table(asv_tab_16S, "16S_ASVs_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# tax table:

    # creating vector of desired ranks
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

    # creating table of taxonomy and setting any that are unclassified as "NA"
tax_tab_16S <- t(sapply(tax_info_16S, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
}))

colnames(tax_tab_16S) <- ranks
row.names(tax_tab_16S) <- NULL
tax_tab_16S <- data.frame("ASV_ID" = sub(">", "", asv_headers_16S), tax_tab_16S, check.names = FALSE)

write.table(tax_tab_16S, "16S_ASVs_taxonomy.tsv", sep = "\t", quote = F, row.names = FALSE)

## saving the seqtab.nochim_16S object
saveRDS(seqtab.nochim_16S, "seqtab.nochim_16S.rds")
```

<br>

---
<br>

## Combining objects in R if wanted
The tables and fasta files we wrote out of R could be combined however we want now, but if we wanted them combined in an object DADA2 like we had for either of them done individually, here's one way we can do that.

```R
class(seqtab.nochim_16S) 
class(seqtab.nochim_18S)
dim(seqtab.nochim_16S)[2] # 510
dim(seqtab.nochim_18S)[2] # 236

head(colnames(seqtab.nochim_18S)) # the sequences are column names in these tables (matrices)
head(rownames(seqtab.nochim_18S)) # and the samples are the rows

  # we can combine them with cbind:
seqtab.nochim <- cbind(seqtab.nochim_16S, seqtab.nochim_18S)
dim(seqtab.nochim) # 2 746
  # this seqtab.nochime object is the same as the individual ones were that we created above following the typical DADA2 processing
  # so we can do whatever to it the same as if we had made it as one the whole time
  # for instance, we can run the taxanomy on the whole thing together the same way we did before merging them
    # (this is the same database we used above, so they will come out the same since we already ran this)

## creating a DNAStringSet object from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

tax_info <- IdTaxa(dna, trainingSet, strand = "both", processors = NULL)

  # and we can write out a combined fasta, count table, and taxonomy table

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)

asv_headers <- vector(dim(seqtab.nochim)[2], mode = "character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}

# fasta:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "All_ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim) %>% data.frame
row.names(asv_tab) <- sub(">", "", asv_headers)
asv_tab <- asv_tab %>% rownames_to_column("ASV_ID")
write.table(asv_tab, "All_ASVs_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# tax table:

    # creating vector of desired ranks
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

    # creating table of taxonomy and setting any that are unclassified as "NA"
tax_tab <- t(sapply(tax_info, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
}))

colnames(tax_tab) <- ranks
row.names(tax_tab) <- NULL
tax_tab <- data.frame("ASV_ID" = sub(">", "", asv_headers), tax_tab, check.names = FALSE)

write.table(tax_tab, "All_ASVs_taxonomy.tsv", sep = "\t", quote = F, row.names = FALSE)
```

<br>

---
<br>

# Evaluating the outcome
From the above processing, we ended up with 237 unique ASVs believed to be derived from 18S fragments, and 511 unique ASVs believed to be derived from 16S fragments. 

That whole process was based on MagicBLAST alignments of the reads to the [v4.14.0 PR2 database](https://github.com/pr2database/pr2database/releases/tag/v4.14.0){:target="_blank"}, so a decent, quick evaluation of if this was all nonsense or not can be done by looking at the taxonomy assigned by [DECIPHER](http://www2.decipher.codes/Classification.html){:target="_blank"} against the [r138 silva database](http://www2.decipher.codes/Classification/TrainingSets/){:target="_blank"} we used in R to do taxonomic classification. So let's look at those taxonomy outputs:

```bash
head -n 5 18S_ASVs_taxonomy.tsv | column -ts $'\t' | sed 's/^/# /'
# ASV_ID     domain     phylum      class               order         family        genus        species
# ASV_18S_1  Eukaryota  Ciliophora  Intramacronucleata  Spirotrichea  Oligotrichia  Strombidium  NA
# ASV_18S_2  Eukaryota  Ciliophora  Intramacronucleata  Spirotrichea  Oligotrichia  Strombidium  NA
# ASV_18S_3  NA         NA          NA                  NA            NA            NA           NA
# ASV_18S_4  NA         NA          NA                  NA            NA            NA           NA

head -n 5 16S_ASVs_taxonomy.tsv | column -ts $'\t' | sed 's/^/# /'
# ASV_ID     domain    phylum         class           order             family             genus    species
# ASV_16S_1  Bacteria  Cyanobacteria  Cyanobacteriia  Synechococcales   Cyanobiaceae       NA       NA
# ASV_16S_2  Bacteria  Bacteroidota   Bacteroidia     Flavobacteriales  Flavobacteriaceae  Formosa  NA
# ASV_16S_3  Bacteria  Bacteroidota   Bacteroidia     Flavobacteriales  NS9 marine group   NA       NA
# ASV_16S_4  Bacteria  Cyanobacteria  Cyanobacteriia  Chloroplast       NA                 NA       NA
```

These are formatted such that the second column holds the domain, so let's see how many each has of Eukaryota, Bacteria, or Archaea in there:

```bash
sed '1d' 18S_ASVs_taxonomy.tsv | cut -f2 | sort | uniq -c | sed 's/^/# /'
#  217 Eukaryota
#   19 NA

sed '1d' 16S_ASVs_taxonomy.tsv | cut -f2 | sort | uniq -c | sed 's/^/# /'
#   12 Archaea
#  345 Bacteria
#  153 NA
```

Not bad! Each has some NAs in there, but there were no Archaea or Bacteria in the 18S split, and there were no Eukarya in the 16S split 👍

---
---
