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
# Example 3 - Converting NCBI taxon IDs to organism lineages  
Recently some collaborators did some sequening with Nanopore's MinION sequencer. This is a "real-time" sequencer that is trying to be as user-friendly as possible so just about anyone can sequence DNA anywhere. They have a whole workflow setup so you can just stream the data to a server as you sequence and it will run things through a taxonomy classification, BUT it just returns you tables with NCBI taxonomy IDs and no actual organism information. 


<br>

---
<br>
# Example 4 - Counting amino acids in thousands of tables
