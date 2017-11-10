---
layout: main
title: Why is this all worth it?
categories: [bash]
permalink: /bash/why
---

{% include _bash_why_toc.html %}

{% include _side_tab_bash.html %}


It can be tough to see why these things are useful before you are actually using them regularly, so here I'm going to provide some real-life examples of how I use *bash* everyday. It is assumed you're already comfortable with *bash* and some of the most common commands, but if not, be sure to check out the [basics](/bash/basics) and [six glorious commands](/bash/six_commands). Without already being comfortable with what's covered in those pages, some of these one-liners might look way more complicated than they actually are.  
<br>

# Getting started  
If you'd like to follow along here, you can get the mock data we'll be working with by coping and pasting these commands into your terminal. If you're not sure what the following commands are doing and are curious, you can find an explanation [here]({{ site.url }}/bash/basics#bottom).

```
cd ~
curl -O https://AstrobioMike.github.io/tutorial_files/bash_why_temp.tar.gz
tar -xvf bash_why_temp.tar.gz
rm bash_basics_why_temp.tar.gz
cd bash_why_temp
```

<br>
# Example 1 - Renaming 1,000 files
So the setup is we have 1,000 files that have names like "Sample-1.fq", "Sample-2.fq", etc., but to use the program we want to use we need the filenames to contain an underscore instead of a dash, like "Sample_1.fq", "Sample_2.fq", etc. Hmm, time to get an undergrad?? Noooo, it's *bash* to the rescue!

Our example_1 directory holds our 1,000 problematic files. 

```
ls | head
ls | wc -l
```

<center><img src="{{ site.url }}/images/1000_files_head.png"></center> 

<br>
Here's the approach. We covered in [*bash* basics](/bash/basics) how to rename files with the `mv` command, so we already know we could easily change one file like this: `mv Sample-1.fq Sample_1.fq`. But of course that doesn't scale, and is just as bad as having the undergrad do this in the Finder window. Instead, we are going to make and do this with a *bash* script. A *bash* script isn't anything special, it's just a bunch of individual *bash* commands one line after another. To demonstrate that, let's make a script out of the two commands we just ran when looking at our files.

`echo` is a command that will spit back out whatever you tell it:

```
echo "Hi Mike"
```

<center><img src="{{ site.url }}/images/echo.png"></center> 

<br>
Seems a little useless at first (unless you're lonely, of course 😢 ). But it's utility makes more sense when you use redirectors to send things somewhere else. For example, instead of using `nano` like we've done before to make a new file, we can make our *bash* script with the `echo` command:

```
echo "ls | head" > test.sh
echo "ls | wc -l" >> test.sh
head test.sh
```

<center><img src="{{ site.url }}/images/make_test_sh.png"></center> 

<br>
The first command redirected "ls | head" into and created the file named "test.sh" (the .sh extension is for "shell"). And the second command appended "ls | wc -l" to that file. Then, as we see with `head`, those are the two lines in the file. We can now run this as a *bash* script. That's done by putting `bash` in front of the file name that is the script:

```
bash test.sh
```

<center><img src="{{ site.url }}/images/test_sh.png"></center> 

<br>
And that ran those commands just the same as if we entered them one at a time actively. This is sometimes referred to as running things in "batch". Also notice our `ls | wc -l` command gave us 1,001 files this time, because there is now also that "test.sh" file in there. If we run `ls *.fq | wc -l`, we'll get our 1,000 again. 

So back to the task at hand, we are going to make a *bash* script that has our 1,000 `mv` commands in it. Here's the whole process:

```
ls *.fq > orig_filenames
tr "-" "_" < orig_filenames > new_filenames
paste -d ' ' orig_filenames new_filenames > both_filenames
sed 's/^/mv /' both_filenames > rename_all.sh
rm *filenames
bash rename_all.sh
``` 

And that's it. If we peek at our files now we see they all have underscores instead of dashes:

<center><img src="{{ site.url }}/images/ex1_done.png"></center> 

<br>

---
<br>
# Example 2 - Making a barcode fasta file to demultiplex samples
In this case we have a mapping file from the sequencing facility that looks like this: 

<center><img src="{{ site.url }}/images/mapping_ex.png"></center> 

<br>
Some sequencing facilities send you back your samples all mixed together, and you need to separate the sequences based on the barcode they have that ties them to which sample they came from. The barcode in this file is the second column, and the sample is specified in the first column. Some programs that do this demultiplexing step require a fasta file to be the "mapping file" to tell it which barcodes belong to which samples. We're going to make that fasta file from this starting file with just a few commands strung together:

```
cut -f1,2 mapping.txt | sed '1d' | sed 's/^/>/' | tr "\t" "\n" > mapping.fa
```

<center><img src="{{ site.url }}/images/making_mapping_fasta.png"></center> 

<br>

---
<br>

# Example 3 - Counting specific amino acids in lots of tables
Recently my good 'ol buddy Josh Kling was yabbering about some paper as he usually does, but for some odd reason I was actually listening this time. It turns out [this paper by Pittera et al.](https://www.nature.com/articles/ismej2016102) was about the thermostability of some *Synechococcus* proteins. They hypothesized that a particular amino acid (position number 43) would more often be alanine in warmer temperature waters, and glycine in colder temperature waters. Having been working on genomics and pangenomics of Syn, I had mapped metagenomic data from about 100 samples from the [TARA Oceans global sampling project](https://www.embl.de/tara-oceans/start/). So we realized I had the data we needed to test this hypothesis just sitting in the computer. All the information was already in the mapping files, we just needed to pull it out. The first part that made that easy is thanks to [Anvi'o](http://merenlab.org/software/anvio/), as that will parse your mapping files for you and give you nice tables like this that we have in our example_3 working directory:

```
column -t ANW_141_05M_24055_AA_freqs.txt | less -S
```

<center><img src="{{ site.url }}/images/aa_counts_ex.png"></center> 

<br>
To get this we just told Anvi'o which gene we were interested in and it spit out this nice table for that specific gene for every environmental sample we had. This is telling us what the Syn population at a specific sample site uses for amino acids at each position within this gene. In the table, every row is an amino acid position in this one gene, and the columns towards the right tell us the total coverage for that position, and the frequency of each amino acid at that position. 

So we had lots of these tables, because there were 32 reference genomes and each had a copy of the gene. Then there were about 32 samples that had greater than 100X coverage that we looked at. In our example here we're only working with one of these samples, so we have just 32 tables for the 32 copies of genes, but the principle is the same.

Back to the task at hand, we want to know how often Syn uses an alanine at amino acid position 43 vs how often it uses glycine at that position. To answer that we needed to pull the appropriate row from all of these tables, and then sum the "coverage", "Ala", and "Gly" columns for comparison. These are columns 7, 8, and 15, which we can double check like so: 

<center><img src="{{ site.url }}/images/aa_counts_2.png"></center> 

<br>
Now that we know we have the right columns, here's the fun part:

```
grep -w "^42" * | cut -f7 | awk '{sum += $1} END {print sum}'
grep -w "^42" * | cut -f8 | awk '{sum += $1} END {print sum}'
grep -w "^42" * | cut -f15 | awk '{sum += $1} END {print sum}'
```

<center><img src="{{ site.url }}/images/aa_counts_results.png"></center> 

<br>
And we see at this particular site the coverage of this specific amino acid position was 1,670, it was an alanine 1,605 of those times, and a glycine only 30 of those times. And sure enough this was a warmer site and we saw the exact opposite when looking at colder sites. Pretty cool! Props to [Piterra et al.](https://www.nature.com/articles/ismej2016102).
<br>

---
<br>
# Example 4 - Converting NCBI taxon IDs to organism lineages  
Recently some collaborators did some sequening with Nanopore's MinION sequencer. This is a "real-time" sequencer that is trying to be as user-friendly as possible so just about anyone can sequence DNA anywhere. They have a whole workflow setup so you can just stream the data to a server as you sequence and it will run things through a taxonomy classification, BUT it just returns you tables with NCBI taxonomy IDs and no actual organism information. 


