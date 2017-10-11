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

Here I would like to introduce you to six glorious commands of *bash* that are awesome to have handy in your toolkit. The problem with a lot of these things is that sometimes it's hard to see why exactly something would be useful at first. And of course before you know that some specific tool exists, you can't exactly realize all the times it would help you. That's how these commands were for me; I didn't know how useful they were, or how often I would use them, until I was well on my way to using them every day.  

I never sat down at any point and specifically tried to learn them. I sort of knew they existed, and occasionally I'd be stuck trying to figure out how to do something and I'd pick up a little more experience with one here-and-there as google led me along. And that's totally fine, but if I had the choice, the first day I started doing anything in terminal I would have started playing with these commands. 

Here we're going to go over what it is exactly that these commands do, and run some simple examples of each with some mock files that you can download if you'd like to follow along. These examples are designed to get you acclimated to just the basic usage of these commands individually. Later when we get to [pipes and redirectors]({{ site.url }}/bash/pipes_redirectors) we'll see how we can easily string multiple commands together to quickly parse files in much more complicated ways, which is really where the utility of these starts to take off. And if you'd like to see some examples of how I end up using these pretty much every day, check out the [Why is this all worth it?]({{ site.url }}/bash/why) post for some 'real-life' examples.  

I'll note again that basic Unix commands like these are for manipulating [plain text files]({{ site.url }}/bash/basics#what-is-a-plain-text-file) only. And also keep in mind that each of these commands are much more expansive than what is presented here.  
<br>

---
<br>
# Getting started  
As with the [bash basics]({{ site.url }}/bash/basics) module, if you'd like to follow along copy and paste these commands into your terminal to download some small example files to work with.

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
`cut` is a command that let's you parse a file by 'fields' (columns). To be able to do this, the plain text file needs to be 'delimited' by something – meaning there needs to be some character that separates each column. This is most often a tab or a comma, and in our example file in this case, "sample_gene_annotations.txt", it is a tab. Let's take a peek at the file with the `head` command:

<center><img src="{{ site.url }}/images/cut_head.png"></center> 

<br>
Yikes, it's kind of hard to tell what's going on due to the line wraps. Let's get a better look with the `less` command. Running `less sample_gene_annotations.txt` would give us a similar view, but if we provide the flag `-S` to the command, it won't wrap lines that reach the edge of the terminal:


```
less -S example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/cut_less.png"></center> 

<br>
Okay, that's a little cleaner. We can see there are some column names in the first row. Don't worry about columns not perfectly lining up, typically the terminal doesn't display things in that fashion (though there are commands that can). As a reminder, `q` will get you out of `less`. Let's look at just the first row using the `head` command and only pulling out 1 line so we can see the column names:


```
head -n1 example_gene_annotations.txt
```

<center><img src="{{ site.url }}/images/cut_head_n1.png"></center> 

<br>
So we can see here this file has 8 columns: "gene_ID", "PC_ID", "genome", "KO_ID", "KO_annotation", "COG_ID", "COG_annotation", and "rRNA". `cut` works based on numerical order. So column 2 would be "PC_ID", for example (it counts starting at 1, not 0). By default `cut` expects the delimiter to be tabs, so we don't need to specify that in this example. All we need to do is tell cut which fields we want (columns), and which file we want it to cut them from. So let's try pulling just column 2 from the file:


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
I realize this is a little anti-climactic at the moment as we haven't even gone over how we can store any of this yet. 


<br>

---
<br>
# grep  
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

# paste  
<br>
<br>

---
<br>
# tr  
<br>
<br>

---
<br>
