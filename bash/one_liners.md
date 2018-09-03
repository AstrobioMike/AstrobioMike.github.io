---
layout: main
title: Useful one-liners I always have handy
categories: [bash, tutorial]
permalink: /bash/one_liners
---

{% include _bash_one_liners_toc.html %}

{% include _side_tab_bash.html %}

Here are some *bash* one-liners I have stored in a text document on my desktop. I've added things to it over time that I either didn't use often enough to have memorized, or that solved some problem for me that was particularly frustrating.  

I have started [this repository at github](https://github.com/AstrobioMike/bioinf_tools){:target="_blank"} which holds a lot of these that I use all the time in more user-friendly scripts, and installation help for that can be found [here](/bash/installing_tools#my-bioinformatics-tools-bit){:target="_blank"}. 

<br>

```bash
  ############# removing added linewraps from fasta files #############
$ awk '!/^>/ { printf "%s", $0; n="\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' in.fa > out.fa 

  ############# using shell variables within awk #############
$ for i in $(cat gORF_gene_IDs); do awk -v id="$i" '$3 == id { print $2 }' genes_in_splits.txt; done > gORF_gene_splits

  ############# fixing ^M carriage returns that excel likes to put in #############
$ tr "\r" "\n" < paired_dists.csv > paired_dists_fixed.csv

  ############# check out column names to pick which to cut ##############
$ head -n1 table.tsv | tr "\t" "\n" | cat -n

  ############# print from line "x" to end of file (here from line 2 to leave off a header) #############
  		# thanks to Xabier Vázquez-Campos for showing me the "+" usage with `tail` as I was doing this ridiculously convoluted before that :) 
$ tail -n +2 test.txt

  ############# delete blank lines #############
$ sed '/^$/d' in.txt

  ############# de-interleave fastq files #############
$ paste - - - - - - - - < interleaved.fq | tee >(cut -f 1-4 | tr '\t' '\n' > R1.fq) | cut -f 5-8 | tr '\t' '\n' > R2.fq

  ############# quick fastq to fasta #############
$ paste - - - - < file.fq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > file.fa

  ############# count number of bases in a fasta file #############
$ grep -v ">" file.fa | wc | awk '{print $3-$1}'

 ############# replace full word only with sed on mac (enclose string with "[[:>:]]...[[:<:]]") #############
$ sed "s/[[:<:]]1[[:>:]]/A/" contig_num 

  ############# add '>' to every other line #############
$ sed 's/^/>/;n' 

  ############# interleave two files line by line #############
$ paste -d '\n' temp_headers temp_seqs > new_fasta.fa

  ############# modify parts of the variable you're looping through #############
    # this example is grabbing all files that end in "_R1.fq",
      # then changing just the "R1" to "R2",
        # then changing just the "R1" to "merged": 
$ for file in *_R1.fq
> do
>   echo $file
>   echo ${file/_R1/_R2}
>   echo ${file/_R1/_merged}
> done

  ############# sed print the 2nd line of every 4 lines #############
$ sed -n '2~4p' in.fq
```
