---
layout: main
title: Useful one-liners I always have handy
categories: [bash, tutorial]
permalink: /bash/one_liners
---

{% include _bash_one_liners_toc.html %}

{% include _side_tab_bash.html %}

Here's a collection of *bash* one-liners I have stored in a text document on my desktop. I've added things to it over time that I either didn't use often enough to have memorized, or that solved some problem for me that was particularly frustrating and so I saved the code there for posterity purposes.  

I have started [this repository]() which holds a lot of these that I use all the time in more user-friendly scripts. I'll get back to dealing with this mess later :)  

<br>

---
---
<br>

<br>
<center><img src="{{ site.url }}/images/under_construction.jpeg"></center>
<center><h3>UNDER CONSTRUCTION</h3></center>
<br>
<br>


```bash
  ############# removing added linewraps from fasta files:
$ awk '!/^>/ { printf "%s", $0; n="\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' in.fa > out.fa 


  ############# using shell variables within awk
$ for i in $(cat gORF_gene_IDs); do awk -v id="$i" '$3 == id { print $2 }' genes_in_splits.txt; done > gORF_gene_splits


  ############# splitting fastq to fasta with sed:
$ sed -n '1~4s/^@/>/p;2~4p' in.fq

  ############# fixing ^M carriage returns that excel likes to put in, so they are newline characters instead:
$ tr "\r" "\n" < paired_dists.csv > paired_dists_fixed.csv

  ############# check out column names to pick which to cut ##############
$ head -n1 pr2_version_4.10.0_merged.tsv | tr "\t" "\n" | cat -n


  ############# print from line "x" to end of file (here from line 2 to leave off a header) 
  		# thanks to Xabier Vázquez-Campos for showing me the "+" usage with `tail` as I was doing this ridiculously convoluted before that :) 
$ tail -n +2 test.txt
  

  ############# delete blank lines with mac (darwin) sed:
$ sed '/^$/d' in.txt

  ############# quick de-interleaving fastq files:
$ paste - - - - - - - - < UW179A_orange_QC_reads.fastq | tee >(cut -f 1-4 | tr '\t' '\n' > reads1.fastq) | cut -f 5-8 | tr '\t' '\n' > reads2.fastq

  ############# quick fastq to fasta:
$ paste - - - - < file.fq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > file.fa

  ############# replacing newline characters (and some dashes in this case) with mac's weird sed:
$ sed -e ':a' -e 'N' -e '$!ba' -e 's/--\n//g' test.faa | head

  # count number of bases in a fasta file:
$ grep -v ">" file.fa | wc | awk '{print $3-$1}'

  # add my typical blast output header to a delimited output file (on osx bash):
$ cat <(printf "qseqid\tqlen\tsseqid\tslen\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tscore\n") ggpS_blastp_out.txt > temp && mv temp ggpS_blastp_out.txt

 ## replace full word only with sed on mac (enclose string with "[[:>:]]...[[:<:]]"):
$ sed "s/[[:<:]]1[[:>:]]/A/" contig_num 

  ## add '>' to every other line:
$ sed 's/^/>/;n' 

  ## interleave two files line by line:
$ paste -d'\n' temp_headers temp_seqs > new_fasta.fa

  ## how to modify parts of the variable you're looping through:
    # this example is grabbing all files with R1, then changing just the R1 to R2, then changing just the R1 to merged (within the brace, specifying the variable, the pattern to be replaced, the pattern to replace it with):
$ for R1 in *_R1.fq; do echo $R1; echo ${R1/_R1/_R2}; echo ${R1/_R1/_merged}; done


  ### sed print the 2nd line of every 4 lines:
$ sed -n '2~4p' 1-F1_S1.fq
```

# Coming soon...

Feel free to [bug me](https://twitter.com/AstrobioMike) to expedite 😊
