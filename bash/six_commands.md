---
layout: main
title: Six commands worth getting to know
categories: [bash, tutorial]
permalink: /bash/six_commands
---

{% include _bash_6_commands_toc.html %}

{% include _side_tab_bash.html %}

If you aren't already moderately comfortable with working at the command line, I recommended you run through the [*bash* basics](/bash/basics) page first. And if you see something here you aren't familiar with, you will be able to find it there. 

Here I would like to go a touch beyond the basics and introduce you to six glorious commands of *bash* that are absolutely worth having handy in your toolkit. The problem with a lot of these things is that sometimes it's hard to see why exactly something would be useful at first. And of course before you know that some specific tool exists and how it works, you can't exactly realize all the times it would help you. That's how these commands were for me; I didn't know how useful they were, or how often I would use them, until I was well on my way to using them every day.  

To try to help with that, first we're going to go over what it is exactly that these commands do, and run some simple examples of each with mock files that you can download if you'd like to follow along. These are mostly designed to get you acclimated to just the basic usage of these commands individually, but we will also slowly begin to incorporate [pipes and redirectors](/bash/basics#pipes-and-redirectors) to get an idea of how easy and useful it is to stick multiple commands together and write results to a file. Once you're comfortable with these things, we'll go over some more-involved [real-life examples](/bash/why) of how I use these commands every day – to hopefully help bridge the gap between first exposure to these things, to seeing how they're actually useful.

I'll note again that basic Unix commands like these are for manipulating [plain text files](/bash/basics#whats-a-plain-text-file) only. And also keep in mind that each of these commands is much more expansive than what is presented here, so always explore at will!  
<br>

---
---
<br>
# Getting started  
As with the [*bash* basics]({{ site.url }}/bash/basics) module, if you'd like to follow along, copy and paste these commands into your terminal to download some small example files to work with and make sure we're all in the same place. If you're not sure what the following commands are doing and are curious, you can find an explanation [here](/bash/basics#bottom).

```
cd ~
curl -O https://AstrobioMike.github.io/tutorial_files/bash_six_commands_temp.tar.gz
tar -xvf bash_six_commands_temp.tar.gz
rm bash_six_commands_temp.tar.gz
cd bash_six_commands_temp
```

We'll mostly be messing with one file here called "example_gene_annotations.txt", which is a table (a tiny one, for practical purposes here). To give us an idea of what this file looks like in a way you're already used to, here's a screenshot of it in Excel: 
<center><img src="{{ site.url }}/images/excel_tab.png"></center> 

<br>
We can see this file has 8 columns: "gene_ID", "PC_ID", "genome", "KO_ID", "KO_annotation", "COG_ID", "COG_annotation", and "rRNA", and there are 10 rows – the first of which is a header with our column names. At times we will find ourselves You will at times find yourself working with files that have too many rows or columns to be opened in Excel, or may just require more memory to fully open than your computer can handle. So let's try to deduce the same information by looking at this file in the terminal.  

Let's first take a peek at the file with the `head` command:

```
head example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/cut_head.png"></center> 

<br>
Yikes, it's kind of hard to tell what's going on due to the line wraps. These particular ones will scale with your terminal window (as they are just "soft wraps", so in some cases it may help to simply adjust the size of the window, but if the lines are just too long that won't help either. Let's get a better look with the `less` command. Running `less example_gene_annotations.txt` would give us a similar view (with lines soft-wrapped), but if we provide the `-S` flag to the command, it won't wrap lines that reach the edge of the terminal:

```
less -S example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/cut_less.png"></center> 

<br>
Okay, that's a little cleaner. We can now see clearly that the first row is a header with names for each column. Don't worry about columns not perfectly lining up, typically the terminal doesn't display things in that fashion (though there are commands that can do that if you'd like). As a reminder, `q` will get you out of `less`. Let's look at just the first row using the `head` command and only pulling out 1 line so we can isolate the column names:

```
head -n1 example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/cut_head_n1.png"></center> 

<br>
And now we can more clearly see what our columns are named. Last, let's also see how many lines we have:

```
wc -l example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/cut_num_lines.png"></center> 

<br>
Now that we know something about the file we're working with, let's get to manipulating it. As we get into these commands, we'll also begin writing results out to files with the `>` redirector and stringing multiple commands together with the `|` character. If you need a refresher on these, swing by the [pipes and redirectors]({{ site.url }}/bash/basics#pipes-and-redirectors) section of the [*bash* basics]({{ site.url }}/bash/basics) module.  
<br>

---
<br>
# cut  
`cut` is a command that let's you parse a file by fields (aka columns). To be able to do this, the plain text file needs to be "delimited" by something – meaning there needs to be some character that separates each column. This is most often a tab or a comma, and in our example file in this case, "example_gene_annotations.txt", it is a tab. By default `cut` expects the delimiter to be tabs, so we don't need to specify that in this case. All we need to do is tell cut which fields we want (columns), and which file we want it to cut them from. So let's try pulling just column 2 from our mock file:

```
cut -f2 example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/cut_f2.png"></center> 

<br>
Pretty straightforward. We can also request multiple fields by listing them separated by commas. Here let's imagine we want just the "PC_ID", "KO_ID", and "COG_ID" columns:

```
cut -f2,4,6 example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/cut_246.png"></center> 

<br>
Or, if the columns we want to pull are contiguous, we can provide a range with a `-`:

```
cut -f1-4,6 example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/cut_1_4.png"></center> 

<br>
Great. Printing something to the terminal like this is useful when we just want to check for something real quick, but most of the time we want that information for something else. So now we're going to use our `>` redirector to write this out to a file. 

```
cut -f1-4,6 example_gene_annotations.txt > genes_and_identifiers_only.txt
```

This created the file "genes_and_identifiers_only.txt" and put the output from our `cut` command in it. Let's take a look:

<center><img src="{{ site.url }}/images/cut_new_file_head.png"></center> 

<br>
Before moving on from `cut`, let's just look at a quick example of how specifying the delimiter works. There is another file in our current working directory here that is the same as the one we've been working with, only delimited by commas instead of tabs. This one is called "example_gene_annotations.csv", with the ".csv" extension being for "comma-separated values". But keep in mind it could just as easily have had a ".txt" extension and still be delimited by a comma or anything else. And on that note, you will sometimes see tab-delimited files with the extension ".tsv", for "tab-separated values", rather than ".txt". Part of the information we're gathering when we peek at a file before getting to work on it is seeing how it's delimited. Let's look at this new one with `less` first:

```
less -S example_gene_annotations.csv
```

<center><img src="{{ site.url }}/images/cut_comma_less.png"></center> 

<br>
Now if we want to cut columns 1-4 and 6, like we did above, we need to tell it not to use the default tab character, but instead to use a comma. In `cut`, we do this with the `-d` flag:

```
cut -d "," -f1-4,6 example_gene_annotations.csv
```

<center><img src="{{ site.url }}/images/cut_comma_1_4.png"></center> 

<br>
Try running the same `cut` command on this comma-delimited file, only without specifying that the delimiter is a comma, `cut -f1-4,6 example_gene_annotations.csv`. Take a second to think about why you get back what you do.  
<br>

---
<br>
# grep  
`grep` is a pattern recognition tool. Apparently it stands for "Global Regular Expression Print", as I just learned from google. In its default usage, `grep` will search for whatever string of characters you ask for, in whichever file(s) you specify, and then return the entire lines that contain that string of characters. For one example of how this works, we're going to search our example gene annotation file for something. As we've seen, this file has [KEGG](http://www.genome.jp/kegg/) and [COG](https://www.ncbi.nlm.nih.gov/COG/) annotations. For the moment, let's pretend we're interested in genes that encode for epoxyqueuosine reductase. We know from above that one of our columns holds KEGG identifiers, and if we search KEGG for this gene we find that there are two in the KEGG database with identifiers [K09765](http://www.genome.jp/dbget-bin/www_bget?ko:K09765) and [K18979](http://www.genome.jp/dbget-bin/www_bget?ko:K18979). So let's try to pull them out of our file with `grep`:

```
grep "K09765" example_gene_annotations.txt
```

When `grep` returns nothing like this and you're just given your prompt back, that means it didn't find the pattern we searched for. When it finds the pattern in the file, it will print the containing lines to the terminal (unless otherwise directed). Let's try the next one:

```
grep "K18979" example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/grep_ko.png"></center> 

<br>
So these lines may be soft-wrapping, but you can see it printed to the terminal the two lines that had epoxyqueuosine reductase. This seems pretty trivial on such a tiny file, but you can imagine this being useful when working with a file 3 million lines long rather than 10.

This is also a good time to note something about the syntax of the `grep` command. We could have searched for both of these terms at the same time, by separating them with backslash and a pipe ( `\|` ). The backslash character is an "escape" character, that tells *bash* you want to use the following pipe as a special character (rather than actually searching for a pipe in your text. And the pipe as a special character in `grep` can be thought of as an "or" operator. So we would pull any lines that have "K09765" **or** "K198979" by entering the command like this:

```
grep "K09765\|K18979" example_gene_annotations.txt
```

Since only one of these patterns are in our file, we still get the output containing only those lines with "K18979".

<center><img src="{{ site.url }}/images/grep_ko.png"></center> 

<br>
Now say we just wanted to know which genomes had this gene (we saw above in the header to this file that genomes were listed in column 3). We can "pipe" ( `|` ) the result of the `grep` into the `cut` command like so:

```
grep "K18979" example_gene_annotations.txt | cut -f3
```

<center><img src="{{ site.url }}/images/grep_cut_ko.png"></center> 

<br>
This right here is the foundation of the magic of *bash*, so let's take a second to think about it. Remember from before what `cut` did and how we entered it. When we entered `cut -f2 example_gene_annotations.txt`, it printed the entire second column to the screen of the file we specified. Here, we specify the file in the original `grep` command, which grabs all of the lines in the file that contain "K18979", then we "pipe" ( `|` ) the output of that into the `cut` command to grab only the 3rd column. Note that we didn't explicitly specify the file we're acting on in the call to the `cut` command. Baseline *bash* commands like these will automatically act on the output of the previous command when they are preceded by a `|`. The terms you'll see for this type of output and input are "stdout" (for standard out) and "stdin" (for standard in). The `|` character takes the stdout from the previous command and makes it the stdin for the next command.  
Just to close the loop here, try running cut on column 3 of the full file (i.e. without the `grep` command in front): `cut -f3 example_gene_annotations.txt`. 

<center><img src="{{ site.url }}/images/cut_full_genomes.png"></center> 

<br>
Notice that column 3 from *all* of the lines was printed, including the header, which wasn't in our output when we had first used `grep` to find the lines that contained "K18979". It's also worth taking a second to think about why that's the case.

Taking this a step further, we are working with genes here, so let's say we want to make a file of only those genes annotated as encoding for epoxyqueuosine reductase that has the gene IDs (column 1), the genome they come from (column 3), the KO identifier (column 4), and the COG identifier (column 6). We can use a combination of `grep` and `cut` strung together with the `|` character to parse out what we want, and then write it to a new file with the `>` redirector at the end like we did above. This will redirect the "stdout" into the file we specify, rather than printing the result to the terminal. Let's put it all together and then look at the new file with `head`:

```
grep "K18979" example_gene_annotations.txt | cut -f1,3,4,6 > epoxyque_reduct_genes.txt
head epoxyque_reduct_genes.txt
```

<center><img src="{{ site.url }}/images/grep_cut_ko_redirect.png"></center> 

<br>
The last thing to note here is that we don't have a header column in our new file. Again, this is a small file, so it wouldn't be difficult to do this manually, but it might be difficult if this were a huge file. So here's one way we can do this with *bash* very easily using just the commands we've already covered and one additional one called `cat`. `cat` will print the entire contents of whatever file(s) you specify. This is useful for when you want to stick files together (vertically – like one file, and then the next file, but now both are in the same new file). This is easier to understand in practice than in terminology so let's just use it and then think about it. 

Here's what we're gonna do:  
1) use `head` to pull just the first row of column names from our original file, like we did above;  
2) use `cut` on the output of that to get the column names we want, the same we specified to make our subset file (1,3,4,6);  
3) write that output to a file with the `>` redirector;  
4) use `cat` to stick our new header file together with our subset file;  
5) use `rm` to delete our header file;  
6) use `mv` to replace our old subset file (the one with no header) with our new subset file (that has our header).  

Listing it out like that makes it seem like a lot, but here's what it looks like:

```
head -n1 example_gene_annotations.txt | cut -f1,3,4,6 > temp_header
cat temp_header epoxyque_reduct_genes.txt > temp_new_epoxy_genes
rm temp_header
mv temp_new_epoxy_genes epoxyque_reduct_genes.txt
```

And now our subset file has header names:  

<center><img src="{{ site.url }}/images/new_epoxy_head.png"></center> 

<br>
Look over the code again and the steps as listed above to make sure you follow what's going on here.

But we can simplify this even further. Remember from above that we can use `grep` with multiple search terms by separating them with an escape slash and a pipe ( `|\` ), such that it will pull lines that contain any individual search term. Since we know the header has "gene_ID" at the start of the line, we can use this information to `grep` out the header in addition to the epoxyqueuosine reductase genes, then use the `cut` command to cut things down to just the columns we want with the appropriate header included:

```
grep "gene_ID\|K18979" example_gene_annotations.txt | cut -f1,3,4,6 > epoxyque_reduct_genes.txt
```

<center><img src="{{ site.url }}/images/new_epoxy_head2.png"></center> 

<br>
This is another good time to mention a handy flag for `grep`. If we add the flag `-c`, instead of printing out the lines containing the pattern you searched for, `grep` will simply count how many lines contain that pattern and print out that number. For example, we now know there are two lines containing the KEGG identifier "K18979", so if we run `grep` with the `-c` flag it should return the number 2:

```
grep -c "K18979" example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/grep_c.png"></center> 

<br>
I bring this up now, because if any other line in our file had the text string "gene_ID", then our plot above to also pull out the header line would have been flawed. One way we could have checked our assumption was correct (that the string "gene_ID" **only** appeared in the first row and nowhere else) would be running `grep` to count how many were in the file like so:

```
grep -c "gene_ID" example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/grep_c1.png"></center> 

<br>
And by this returning "1", we can be confident we're not accidently grabbing any extra lines we didn't mean to get. 

Many commands in *bash* support what are known as "regular expressions". We've already gotten a little taste of this in the [wildcards]({{ site.url }}/bash/basics#wildcards) section of [bash basics]({{ site.url }}/bash/basics), but as you find yourself wondering if and how you can search for more complicated patters with `grep`, you may want to start googling and learning more about *bash* regular expressions.  
<br>

---
<br>
# paste
We just saw that `cat` sticks things together "vertically". `paste` is our command for sticking things together "horizontally". It treats the things you're sticking together like individual columns and adds some delimiter in between them. 

We have another file we have in our working directory called "aa_lengths.txt". Let's look at it:

```
head aa_lengths.txt
```

<center><img src="{{ site.url }}/images/head_lengths.png"></center> 

<br>
Here we have 2 columns, the first is our gene IDs, the second is how many amino acids long they are. First, just for an example of how `paste` works, let's add these columns to our "genes_and_identifiers_only.txt" file:

```
paste genes_and_identifiers_only.txt aa_lengths.txt
```

<center><img src="{{ site.url }}/images/paste1.png"></center> 

<br>
If we wanted to specify a different delimiter, the `paste` command used the `-d` flag just like `cut`. 

Now that we've seen how `paste` works, let's say we want a file like our "epoxyque_reduct_genes.txt" but that also has this this amino acid length information in it, such that the output looks like this: 

<center><img src="{{ site.url }}/images/head_aa_added.png"></center> 

<br>
There are lots of ways to arrive at this file as an end product. Before looking at how I did it below, see if you can create it on your own by using the things we've covered so far.

...  

...  

...  

...  

...  

...  

Okay, enough ellipses, if you're still reading by this point you're just trying to cheat. 

```
paste example_gene_annotations.txt aa_lengths.txt | cut -f1,3,4,6,10 | grep "gene_ID\|K18979" > epoxyque_reduct_genes.txt
```
<br>

---
<br>

# sed  
`sed` stands for "stream editor". It is a command that does a lot of different things, but probably its most common use is to find a specified string of text in a file, and replace it with something else. You have probably used a search-and-replace like this in a text editor or in excel. `sed` is another *bash* command where "regular expressions" become extremely powerful when you get the hang of them, but for now we'll just look at some very basic usage.

Let's look at our "genes_and_identifiers_only.txt" file again:

```
head genes_and_identifiers_only.txt
```

<center><img src="{{ site.url }}/images/cut_new_file_head.png"></center> 

<br>
And now let's pretend we got updated information about our reference genomes, and what we were previously identifying as genome "GEYO" is now supposed to be "MLee09". The syntax of `sed` may seem a little strange at first, so let's run it and then break it down. 

```
sed 's/GEYO/MLee09/' genes_and_identifiers_only.txt
```

<center><img src="{{ site.url }}/images/sed1.png"></center> 

<br>
Here, the `sed` command is followed by an expression within single quotes. I'll note there is a difference between single quotes and double quotes that you should look into as you begin using them more. But for now it's probably better to use single quotes and when you come across a case where things aren't happening the way you expect, think about the quotes first and look into them more at that point. 

Back to our `sed` search-and-replace expression, within the quotes we actually have 4 items delimited by forward slashes: the first is an `s`, which is for "substitute"; the second is our string we'd like to find, `GEYO`; the third is what we'd like to replace it with, `MLee09`; and the fourth is actually empty right now. To see what the fourth item does, let's look at another example.

Now, let's say we want to replace every "NA" with "<NA>". Let's try the command as we ran it above first:

```
sed 's/NA/<NA>/' genes_and_identifiers_only.txt
```

<center><img src="{{ site.url }}/images/sed2.png"></center> 

<br>
Notice that only the first instance of "NA" in each line was changed. This is the behavior of `sed` when we leave that 4th item empty. In order to replace all items in each line that match our search pattern, we need to provide a `g` for "global" as that 4th item:

```
sed 's/NA/<NA>/g' genes_and_identifiers_only.txt
```

<center><img src="{{ site.url }}/images/sed3.png"></center> 

<br>
`sed` in particular has a ton of functionality and flexibility, we'll see more in-depth applications of it when we get to some real-life examples.

<br>

---
<br>
# awk  
`awk`, like `sed`, is also very expansive. `awk` is the command I go to when I want to do any sort of filtering based on numeric values, or if I want to do any sort of calculations like sum a column.

The syntax of `awk` can also take a little getting used to. This is one of the commands that whenever I need to use it for something new and even remotely complicated, I usually pull up an old file I have or do some googling to remind myself of how to do whatever I'm trying to figure out (this is normal! none of this is about memorizing these things). 

For some examples of `awk`, we're going to work with a typical BLAST output table. Let's take a look at it first with `head`:

<center><img src="{{ site.url }}/images/head_blast.png"></center> 

<br>
And actually now is a good time to show an example of the `column` command. This will keep columns organized together: 

```
column -t example_blast_output.txt
```

<center><img src="{{ site.url }}/images/head_blast_col.png"></center> 

<br>
So we see here we have 6 columns: "qseqid" is the query sequence (what we're blasting); "qlen" is the length of the query sequence; "sseqid" is the subject sequence (what our query hit in our blast database we were blasting against); "slen" is subject length; "pident" is percent identity; and "length" is the length of the alignment.

After blasting you usually want to filter the output by some criteria. For the first example, let's just say we only want to keep hits that were greater than 90% identical. For `awk`, we specify which columns we want to act on with a `$` followed by the column number. So in this case, where the percent identity is in column 5 in our blast table, if we provide an expression like this, `$5 > 95`, `awk` will by default print all of the rows that pass this criterion. Let's give it a shot:

```
awk '$5 > 95' example_blast_output.txt 
```

<center><img src="{{ site.url }}/images/awk1.png"></center> 

<br>
And actually before running that I didn't remember how `awk` was going to treat the text of the first row, but it seems to return that also so we don't have to worry about adding it back in. If it hadn't, we would have done something similar to above to keep our header if we had been writing this output to a new file.

We can also do this sort of filtering based on multiple columns by connecting multiple conditions with "and/or" logical operators. An "and" in `awk` is provided with `&&`, and an "or" is provided with `||`. Just for example's sake, let's say we want our 95% ID criterion, but we also want the query length to be greater than 1000:

```
awk '$5 > 95 && $2 > 1000' example_blast_output.txt
```

<center><img src="{{ site.url }}/images/awk2.png"></center> 

<br>
Percent identity alone doesn't tell us enough information though, as you can have just a tiny portion of a query sequence align to a subject sequence, but still have a high percent identity over that short alignment. So another filtering metric that might typically be incorporated is the alignment length (column 6 here) and the length of the query sequence (column 2 here), to say something like: "only keep the hits that are greater than 95% identical AND have alignments that are greater than 90% of the length of the query sequence. That would look something like this:

```
awk '$5 > 95 && $6 > $2*0.9' example_blast_output.txt
```

<center><img src="{{ site.url }}/images/awk3.png"></center> 

<br>
And although there'd be no reason to do so in this case, let's sum a column just to have the example here. Let's assume we want to sum the length of all the alignments (column 6): 

```
awk '{sum += $6} END {print sum}' example_blast_output.txt
```

<center><img src="{{ site.url }}/images/awk4.png"></center> 

<br>
Here we've introduced two new components: 1) squiggly brackets which are there to separate the two actions of adding up the column ("sum" is just a variable name here, it could be anything), vs printing the variable "sum" at the end; and 2) the new concept of telling `awk` when to "END" whatever it is doing, before going onto the next step. Let's look at what happens when we don't provide that "END" component:

```
awk '{sum += $2} {print sum}' aa_lengths.txt
```

<center><img src="{{ site.url }}/images/awk5.png"></center> 

<br>
That time it printed the value of our "sum" variable after it did each line, rather than waiting until it was done with the whole input. Again the syntax of `awk` takes some getting used to, and pretty much everytime I use it I look at an old log file or google things for it. But there are times I wouldn't know what to do without out, like we'll see in the real-life examples we'll get to soon.  
<br>

---
<br>
# tr
`tr` is for "translate". It is useful for changing single characters to another character, and more often I end up using it to change special characters to other special characters (when `sed` is problematic). 

For example, when you export a table from Excel as a tab-delimited file, it does it with odd characters for the end of a line. Normally, unix-like systems interpret a `\n` as the end of a line. But Excel for some reason exports things with `\r` instead. We can see this if we run `less` on an Excel exported file. We have one of these in our working directory, called "example_gene_annotations_excel_exported.txt". If we run the `less` command on that, we'll get only one "line" returned:

<center><img src="{{ site.url }}/images/excel_exp.png"></center> 

<br>
If you scroll to the right with the right arrow key, you'll see this is all on one "line" as far as `less` is concerned. That's because Excel places `\r` characters at the end of each line (which show up as highlighted "^M" characters), rather than `\n` characters (which are the normal "newline" characters). It may feel like "who cares?" at this point, but unix does. And I would love to hear from a Microsoft developer exactly why they do this! Anyway, we can fix the problem using `tr` like this:

```
tr "\r" "\n" < example_gene_annotations_excel_exported.txt > example_gene_annotations_fixed.txt
```

Notice that this command is different that the previous ones we've seen. `tr` requires that the input file is specified by the `<` character. Then, here we are redirecting the output to a new file with the `>` character. 

Now if we look at this file, we see that it is formatted in a way that the rest of the world can interpret:

```
less -S example_gene_annotations_fixed.txt
```

<center><img src="{{ site.url }}/images/excel_exp2.png"></center> 

<br>
And now we could begin to work on this file as though it weren't an alien.  
<br>

---
---
<br>
# Ok, so where's the advertised gloriousness?
I know, this is all well and good, but again I know firsthand that it's really hard to care about these things until you are already using them to make your life easier. But the problem is that it's also kinda hard to get to that point *before* you care enough about using them to learn them. This is the conundrum we face 😕  

To hopefully help bridge that gap, let's hop over to [Why is this all worth it?]({{ site.url }}/bash/why) for some real-life examples of how I use these things every day.
