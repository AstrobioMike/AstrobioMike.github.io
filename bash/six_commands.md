---
layout: main
title: Six commands worth getting to know
categories: [bash, tutorial]
permalink: /bash/six_commands
---

{% include _bash_6_commands_toc.html %}

{% include _side_tab_bash.html %}

<center><img src="{{ site.url }}/images/under_construction.jpeg"></center>
<center><h3>UNDER CONSTRUCTION</h3></center>
<br>
<br>

Here I would like to introduce you to six glorious commands of *bash* that are absolutely worth having handy in your toolkit. The problem with a lot of these things is that sometimes it's hard to see why exactly something would be useful at first. And of course before you know that some specific tool exists and how it works, you can't exactly realize all the times it would help you. That's how these commands were for me; I didn't know how useful they were, or how often I would use them, until I was well on my way to using them every day.  

I never sat down at any point and specifically tried to learn them. I sort of knew they existed, and occasionally I'd be stuck trying to figure out how to do something and I'd pick up a little more experience with one here-and-there as google led me along. And that's totally fine, but if I had the choice, the first day I started doing anything in the terminal I would have started playing with these commands. 

Here we're going to go over what it is exactly that these commands do, and run some simple examples of each with mock files that you can download if you'd like to follow along. These examples are mostly designed to get you acclimated to just the basic usage of these commands individually, but we will also slowly begin to incorporate [pipes]({{ site.url }}/bash/basics#pipes-and-redirectors) to get an idea of how easy and useful it is to stick multiple commands together. Then at the end, we'll see some real-life examples of how I use these types of commands every day. This is to hopefully help bridge the gap from first exposure to these things, to seeing their actual utility.

I'll note again that basic Unix commands like these are for manipulating [plain text files]({{ site.url }}/bash/basics#whats-a-plain-text-file) only. And also keep in mind that each of these commands is much more expansive than what is presented here.  
<br>

---
<br>
# Getting started  
As with the [bash basics]({{ site.url }}/bash/basics) module, if you'd like to follow along, copy and paste these commands into your terminal to download some small example files to work with.

```
cd ~
curl -O https://AstrobioMike.github.io/tutorial_files/bash_six_commands_temp.tar.gz
tar -xvf bash_six_commands_temp.tar.gz
rm bash_basics_six_commands_temp.tar.gz
cd bash_six_commands_temp
```
<br>

---
<br>
# cut  
`cut` is a command that let's you parse a file by fields (aka columns). To be able to do this, the plain text file needs to be "delimited" by something – meaning there needs to be some character that separates each column. This is most often a tab or a comma, and in our example file in this case, "example_gene_annotations.txt", it is a tab. Let's take a peek at the file with the `head` command:

<center><img src="{{ site.url }}/images/cut_head.png"></center> 

<br>
Yikes, it's kind of hard to tell what's going on due to the line wraps. These particular ones will scale with your terminal window (as they are just "soft wraps", so in some cases it may help you to adjust the size, but sometimes that won't help. Let's get a better look with the `less` command. Running `less example_gene_annotations.txt` would give us a similar view, but if we provide the flag `-S` to the command, it won't wrap lines that reach the edge of the terminal:


```
less -S example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/cut_less.png"></center> 

<br>
Okay, that's a little cleaner. We can see there are some column names in the first row. Don't worry about columns not perfectly lining up, typically the terminal doesn't display things in that fashion (though there are commands that can do that if you'd like). As a reminder, `q` will get you out of `less`. Let's look at just the first row using the `head` command and only pulling out 1 line so we can see the column names:


```
head -n1 example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/cut_head_n1.png"></center> 

<br>
So we can see here this file has 8 columns: "gene_ID", "PC_ID", "genome", "KO_ID", "KO_annotation", "COG_ID", "COG_annotation", and "rRNA". `cut` works based on numerical order. So, column 2 for example, woulc be "PC_ID", (it counts starting at 1, not 0). By default `cut` expects the delimiter to be tabs, so we don't need to specify that in this case. All we need to do is tell cut which fields we want (columns), and which file we want it to cut them from. So let's try pulling just column 2 from the file:


```
cut -f2 example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/cut_f2.png"></center> 

<br>
Pretty simple. You can also request multiple fields by listing them separated by commas:


```
cut -f2,4,6 example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/cut_246.png"></center> 

<br>
Or if the columns you want to pull are contiguous, you can provide a range with a `-`:

```
cut -f1-4 example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/cut_1_4.png"></center> 

<br>
For a quick example of how specifying the delimiter works, let's look at the same file, only delimited by commas instead of tabs. This one is called "example_gene_annotations.csv", the ".csv" extension being for "comma-separated values", but keep in mind it could just as easily have a ".txt" extension and still be delimited by a comma or anything else. 

```
less -S example_gene_annotations.csv
```

<center><img src="{{ site.url }}/images/cut_comma_less.png"></center> 

<br>
Now, if we want to cut the first four columns like just above, we need to tell it not to use the default tab character, but instead to use a comma. We do this with the `-d` argument.

```
cut -d "," -f1-4 example_gene_annotations.csv
```

<center><img src="{{ site.url }}/images/cut_comma_1_4.png"></center> 

<br>
Try running cut on the comma-delimited file without specifying that the delimiter is a comma, `cut -f1-4 example_gene_annotations.csv`, and try to understand why you get back what you do.  
<br>

---
<br>
# grep  
`grep` is a pattern recognition tool. Apparently it stands for "Global Regular Expression Print", as I just learned from google. In its default usage, `grep` will search for whatever string of characters you ask for, in the file you specify, and then return the entire lines that contain that string of characters. For one example of how this works, we're going to search our example gene annotation file for something. This file has [KEGG](http://www.genome.jp/kegg/) and [COG](https://www.ncbi.nlm.nih.gov/COG/) annotations. For the moment, let's pretend we're interested in genes that encode for epoxyqueuosine reductase. We know from above that one of our columns holds KEGG identifiers, and if we search KEGG for this we'd find that it's identifier is [K18979](http://www.genome.jp/dbget-bin/www_bget?ko:K18979). So let's try to pull that out of our file with `grep`:

```
grep "K18979" example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/grep_ko.png"></center> 

<br>
So this may be soft-wrapping, but you can see it printed to the terminal the two lines that had epoxyqueuosine reductase. This seems pretty trivial on such a tiny file, but you can imagine this being useful when working with a file 3 million lines long rather than 10. Say we just wanted to know which genomes had this gene (we saw above in the header to this file that genomes were listed in column 3). We can "pipe" the result of the `grep` into the `cut` command like so:

```
grep "K18979" example_gene_annotations.txt | cut -f3
```

<center><img src="{{ site.url }}/images/grep_cut_ko.png"></center> 

<br>
Let's take a second to understand what's going on here, because this is the foundation of the magic of *bash*. Remember from before what `cut` did and how we entered it. When we entered `cut -f2 example_gene_annotations.txt`, it printed the entire second column to the screen of the file we specified. Here, we specify the file in the original `grep` command, which grabs all of the lines in the file that contain "K18979", then we "pipe" `|` it into the `cut` command and grab the 3rd column. Note that we didn't explicitly specify the file we're acting on in the call to the `cut` command. Baseline *bash* commands like these will automatically act on the output of the previous command when they are preceded by a `|`. The terms you'll see for this type of output and input are "stdout" (for standard out) and "stdin" (for standard in). The `|` character takes the stdout from the previous command and makes it the stdin for the next command.  
Just to close the loop here, try running cut on column 3 of the full file (i.e. without the `grep` command in front): `cut -f3 example_gene_annotations.txt`. 

<center><img src="{{ site.url }}/images/cut_full_genomes.png"></center> 

<br>
Notice that column 3 from *all* of the lines was printed, including the header, which wasn't in our output when we had first used `grep` to find the lines that contained "K18979". It's worth taking a second to think about why that's the case.

Moving on, we are working with genes here, so let's say we want to make a file that has only the gene IDs (column 1), the genome they come from (column 3), the KO identifier (column 4), and the COG identifier (column 6), and only for the genes annotated as encoding for epoxyqueuosine reductase. We can use a combination of `grep` and `cut` strung together with the `|` character to do this.

```
grep "K18979" example_gene_annotations.txt | cut -f1,3,4,6
```

<center><img src="{{ site.url }}/images/grep_cut_ko_parse.png"></center> 

<br>
Right now, the output (stdout) is just printed to the terminal. Instead, let's save this into a new file with `>` redirector. Here we'll send this into a file called "epoxyque_reduct_genes.txt". (Revisit the [pipes and redirectors]({{ site.url }}/bash/basics#pipes-and-redirectors) section of [*bash* basics]({{ site.url }}/bash/basics). 

```
grep "K18979" example_gene_annotations.txt | cut -f1,3,4,6 > epoxyque_reduct_genes.txt
head epoxyque_reduct_genes.txt
```

<center><img src="{{ site.url }}/images/grep_cut_ko_redirect.png"></center> 

<br>
The last thing to note here is that we don't have a header column isn't present in our new file

Now this is a tiny file for example purposes of course, but you can imagine how this could be useful when working with a file 800,000 lines long rather than 8.


<br>
<br>

---
<br>
# paste
`paste` is a command that sticks columns together Now another file we have in our working directory here
<br>
<br>

---
<br>

# awk  
<br>
<br>

---
<br>
# sed  
<br>
<br>

---
<br>

# tr  
<br>
<br>

---
<br>

# Some real-life examples
So yeah, this is all well and good, but I know firsthand that it's really hard to care about these things until you are already using them to make your life easier. But it's also kinda hard to get to that point *before* you care about using them. This is the conundrum.  
The versatility of how these tools can be strung together is what makes them so powerful and useful, but that also makes it difficult to demonstrate; there's no *one* "typical" instance that shows why this is worth it.  
I now use them virtually all day every day. So here are a couple examples of how
