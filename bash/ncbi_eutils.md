---
layout: main
title: Downloading from NCBI with E-Utils
categories: [bash, tutorial]
permalink: /bash/ncbi_eutils
---

{% include _bash_downloading_from_ncbi_toc.html %}

{% include _side_tab_bash.html %}

[NCBI](https://www.ncbi.nlm.nih.gov/){:target="_blank"} is pretty damn awesome. But the first few times I wanted to download a massive amount of reference sequences or genomes I found myself struggling a bit. If that has happened to you, then hopefully this page helps out. NCBI's [Entrez Direct E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK179288/){:target="_blank"} offers one avenue to be able to download data in bulk at the command-line. If you don't have it yet, that link above provides installation instructions and/or you can walkthrough the process [here](/bash/installing_tools#ncbis-e-utilities){:target="_blank"}.  

I don't use this toolset much, and when I do it's only to pull genes or genomes, so this certainly won't be exhaustive at all. Make sure to look over the full functionality [here](https://www.ncbi.nlm.nih.gov/books/NBK25499/){:target="_blank"} sometime. For now, here are some examples of using the `efetch` command to pull sequence data.  

## The efetch command
The `efetch` command let's you pull all kinds of data from NCBI. If you run `efetch -help`, you can look at lots of parameters and types of info you can pull. Here, to get an idea of how the command works, let's just pull one amino acid coding sequence with one accession number for an alkaline phosphatase:

```bash
efetch -db protein -format fasta -id AEE52072.1
```

And after a second the sequence should print out to the screen:

<center><img src="{{ site.url }}/images/eutils_efetch1.png"></center>

<br>
These are some of the typical flags you need to supply to `efetch` or other E-utils commands: `-db` to specify which database; `-format` to tell it how you want the data; and -id to provide the desired accession numbers or unique IDs. Note that the default behavior just prints the output to the terminal, so to save the output you need to [redirect it](/bash/basics#pipes-and-redirectors){:target="_blank"} to a file.  

The efetch command can also take multiple IDs separated by commas. Here's an example pulling two sequences and writing the output to a new file:

```bash
efetch -db protein -format fasta -id AEE52072.1,ADV47642.1 > my_seqs.faa
```

In practice of course we can download one or two from the site though, and we're only using this because we want a lot. Unfortunately you can't provide the `-id` argument a file of accession numbers, and even if you were to actively enter them on the command line, the [Entrez site notes](https://www.ncbi.nlm.nih.gov/books/NBK179288/#chapter6.Automation){:target="_blank"} that you shouldn't do more than blocks of 200 at a time due to server limitations. So next let's look at first how to generate a large list of accession numbers, and then how to solve all our problems with the magic of *bash* :)

## Pulling lots of sequences
For an example, let's imagine we want all the amino acid sequences of the *phoD*-type alkaline phosphatases available in [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/){:target="_blank"} for bacteria (because Euks are too hard). While this is focused on amino acid coding sequences, the same principles apply if you wanted nucleotide sequences of genes or whole genomes. The only things that would change would be how you search your accessions and which options you specify to `efetch`.   

### Generating accessions list
As we just saw, to use `efetch` at the command line we first need to generate a list of accession numbers (or gene IDs). This can be done at the command line too with the `esearch` command, but I personally just do it on their web page. Here are the steps I just took to get the desired accessions for bacterial *phoD*-type amino acid sequences:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; • went to [NCBI](https://www.ncbi.nlm.nih.gov/){:target="_blank"}  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; • changed the search database from "All Databases" to "Protein"  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; • searched for "alkaline phosphatase phoD"  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; • limited the search to only RefSeq by clicking it under "Source databases" on the left side  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; • limited the search to only bacteria by clicking that in the top right  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; • clicked "Send to:", "File", changed format to "Accession List", clicked "Create File"  

At the time of my doing this, this was a total of 10,249 accessions that were written to a file called "sequence.seq":  

<center><img src="{{ site.url }}/images/eutils_efetch2.png"></center>

<br>
For the sake of this example we don't need that many, so I'm going to cut that down to about a 10th and store them in a file called "wanted_accessions.txt":

```bash
head -n 1025 sequence.seq > wanted_accessions.txt
```

### Formatting for bulk download
Remember from the example above that `efetch` can take multiple accessions separated by commas. To see how we can format our accessions list properly, first let's use *bash* to build up an `efetch` command that will run on just the first 10 seqs:

```bash
head wanted_accessions.txt | tr "\n" "," | sed 's/,$//' > ten_formatted.txt
```
Here I used the `head` command to just grab the first 10 accessions, then used the `tr` command to change all newline characters to commas, and then `sed` to remove the very last trailing comma (see the [*bash* basics](/bash/basics){:target="_blank} and [six glorious commands](/bash/six_commands){:target="_blank} pages if you're not yet familiar with these commands).  

<center><img src="{{ site.url }}/images/eutils_efetch3.png"></center>

<br>
That's softwrapped in the image, but we can see that the 10 accessions are now all on one line and separated by commas. Now we simply need to add the rest of the `efetch` command in front of that. The following code replaces the start of every line (here just one) with the `efetch` command we need in front of the comma-delimited accessions, and then writes the output to a new file called "ten_accessions.sh":

```bash
sed 's/^/efetch -db protein -format fasta -id /' ten_formatted.txt > ten_accessions.sh
```

<center><img src="{{ site.url }}/images/eutils_efetch4.png"></center>

<br>
Now all we need to do is call that file as a *bash* script and redirect the output to a new file:

```bash
bash ten_accessions.txt > ten_phoDs.faa
```

<center><img src="{{ site.url }}/images/eutils_efetch5.png"></center>

<br>
Great. Now that we see how we can format one set of accessions for an `efetch` command, that just leaves: splitting the large file of accessions into multiple smaller files; building the formatted `efetch` command for all of them; and then throwing them all into a shell script together. Here I am going to use the `split` command to split up our large accessions file with 1,025 accessions, "wanted_accessions.txt", into as many 200-line files as are needed:

```bash
split -l 200 wanted_accessions.txt temp_block_
``` 

<center><img src="{{ site.url }}/images/eutils_efetch6.png"></center>

<br>
Here the `split` command made 6 files with the prefix we provided as the last positional argument, and all of them have 200 lines except the last which has the remaining 25. Now we can just loop through those to generate the properly formatted shell script like we did above for the individual one:  

```bash
for block in `ls temp_block_*`; 
  do 
  tr "\n" "," < $block | sed 's/,$//' | sed 's/^/efetch -db protein -format fasta -id /'; 
done > pull_more_phoD_seqs.sh
```

And here's what the newly created "pull_more_phoD_seqs.sh" file looks like (in `less` with no softwrapping, so the accession list runs off to the right):

<center><img src="{{ site.url }}/images/eutils_efetch7.png"></center>

<br>
And after running the script, which took about 15 seconds for these 1,025 sequences, we have our references 🙂 

<center><img src="{{ site.url }}/images/eutils_efetch8.png"></center>

<br>

---
---
<br>
I hope this was helpful. As mentioned above, to pull different types of data you would need to run your search appropriately to get the desired accessions, and then you'd need to change the options of `efetch`, like which `-db` and which `-format`. 
