---
layout: main
title: R basics
categories: [R, tutorial]
permalink: /R/basics
---

{% include _R_basics_toc.html %}

{% include _side_tab_R.html %}

This module is designed for those that are either completely new to R, or have some experience but maybe don't feel as solid about some of the fundamentals as they'd like. It will run through the very basics such as setting up your working environment, assigning variables, "indexing" (subsetting data), reading in and writing out data, and installing packages. Some relevant terminology is presented [here](/R/index) if you find yourself seeing some words that are unfamilar to you. This page is meant to be a quick-start to get you into and using the R environment. There is, of course, an incredible amount of functionality we won't be touching on here, so be sure to peruse the extensive [R intro documentation available here](https://cran.r-project.org/doc/manuals/r-release/R-intro.html) when you're starting to get your bearings and want to go further 🙂 
<br>
<br>

---
---
<br>

# Accessing our R environment
Before we get started, we need an R environment to work in. You can work on your own computer if you'd like, or you can work in a "Binder" that's been created for this page, see below. 

## On your computer
It is possible your computer already has R, if you are unsure, you can check by opening a terminal window ([a unix-like terminal](/bash/getting_bash_env){:target="_blank"}) and typing `R`. If this launches R rather than giving an error message, you should be good to go (enter `q()` to exit the R environment). If you do not have R, you can download it from here for Mac: [https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/){:target="_blank"}. And if you have a relatively newer Mac, you may also need to install XQuartz which you can get from here: [https://www.xquartz.org/](https://www.xquartz.org/){:target="_blank"}. 

Lastly, I highly, *highly*, **highly** recommend installing the free version of [RStudio](https://www.rstudio.com/){:target="_blank"} if you don't already have it. RStudio is an interface for R that not only makes everything you will do in R easier and more organized, but it's also invaluable for reproducibility of your analyses as it makes it second-nature to generate and save R scripts of everything you're doing, while you're doing it – which is very helpful when you want to look back and see what worked out of the 20 things you just tried 🙂. You can download an RStudio installer from [here](https://www.rstudio.com/products/rstudio/download/#download){:target="_blank"}. See the next section for the typical layout of RStudio.  

## With Binder
[Binder](https://mybinder.org/){:target="_blank"} is an incredible project with incredible people behind hit. I'm still pretty new to it, but the general idea is it makes it easier to setup and share specific working environments in support of open science. What this means for us here is that we can just click this little badge – [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AstrobioMike/binder-R-basics/master?urlpath=rstudio){:target="_blank"} – and it'll open up the proper R environment with all our needed example files in a web-browser ready to rock... how awesome is that?!? So yeah, if you want to work in the binder, click it already! 

When that page finishes loading (it may take a couple minutes), you will see a screen like this:

<center><img src="{{ site.url }}/images/binder-R-launch.png"></center>
<br>

If we then click the icon with the green plus sign at the top-left, and then "R Script", it will open up our text editor area also, and look something like this (minus the labels added here): 

<center><img src="{{ site.url }}/images/binder-R-labeled.png"></center>

<br>

---
<br>
# RStudio layout
RStudio has 4 main panes, as numbered above: 1) console; 2) source; 3) Environment; and 4) Files, Plots, Packages, etc. 

1. The "console" is where you can run commands just as though you were working in an R environment at a command line, and it is also where results will print out. 
2. The "source" pane acts as a sort-of interactive text editor within which you can write out and save all of your commands, and then call them when you'd like. This is one of the reasons R Studio is great, you're constantly building up your notes file as you work, without any added effort needed. If you want to run a command written in the source file, while on the line of the command press `Cmd + Enter` or `Ctrl + Enter`, or you can highlight a section and do the same to run that section. 
3. The "Environment" pane displays all of the variables and data structures you currently have stored.  
4. "Files/Plots/Packages/etc." allows you to navigate your computer in the typical Finder fashion, displays any plots you generate, and serves as your help window. 

Here we're going to be doing our work in the "console" area. To start, let's see how we can get help on a function. To do this in R, we just place a `?` in front of the function name. For example, here is how to see the help info for the function to see what our current working directory is in R:

```R
?getwd
```

<center><img src="{{ site.url }}/images/binder-R-getwd-help.png"></center>
<br>
And notice the pane at the bottom right now shows our help info for this function.  
<br>

---
<br>
# Some practice data
**If you are not using the binder environment**, but want to follow along with this page, copy and paste the following commands into your terminal to get set up with a small, temporary working directory that has the files I'll be working with. If you're unfamiliar with working at the command line, and would like to get to know it better, consider running through [this page](/bash/bash_intro_binder){:target="_blank"} when you can 🙂

<center><b>SKIP THESE COMMANDS IF YOU ARE WORKING IN THE BINDER ENVIRONMENT SHOWN IN THE PICTURE ABOVE</b></center><br>

```
cd ~
curl -O https://AstrobioMike.github.io/tutorial_files/R_basics_temp.tar.gz
tar -xvf R_basics_temp.tar.gz
rm R_basics_temp.tar.gz
cd R_basics_temp
```

<br>

---
<br>

# Setting up our working environment
Just like when working at the command line, or pointing to files in a graphical user interface program, we need to be aware of "where" we are in our computer when working in R. The `getwd()` and `setwd()` help us do this in R. Commands in R typically take this structure, with the command being followed by parentheses (so that's how we'll be listing them for the most part here). Inside those parentheses is where the arguments to the command would go if there are any. In the case of `getwd()` no arguments are needed as it is just going to tell us where we currently are:

```R
getwd()
```

<center><img src="{{ site.url }}/images/binder-R-getwd.png" width="90%"></center>
<br>

Note that this is just a zoomed in image of the "console" pane. 

However, in the case of `setwd()`, we then need to tell it where we want to go. Here, we are saying we want to move into the directory holding our example data, and then checking we moved with `getwd()` again:  

```R
setwd("~/R_basics_temp/")
getwd()
```

<center><img src="{{ site.url }}/images/binder-R-setwd.png" width="90%"></center>
<br>

>**NOTE:** Here we are providing the [absolute path](/bash/bash_intro_binder#absolute-vs-relative-path){:target="_blank"} to where we want to go, and it needs to be within quotations. Putting something within quotations in R is what tells it to take the text just as it is, rather than trying to look for a variable that is named with that text. We'll get into that more later, but just to note now why the quotes are important 🙂  

Now that we are in the correct location that should contain the file for our tutorial here, let's check that it is actually here. At the command line this would be done with the `ls` command, here in R we do it with the `list.files()` function.

```R
list.files()
```

<center><img src="{{ site.url }}/images/binder-R-list-files.png" width="90%"></center>
<br>

---
<br>
# Basic calculations
At its core, R is built for statistics, so with it we can do baseline arithmetic operations like the following examples. (Don't forget we are working in the "console" pane as pictured above, these blocks can be copy-and-pasted there if wanted.)

```R
4 + 4
4 / 2
4 * 4
2 ^ 4 
```

<center><img src="{{ site.url }}/images/R_basic.png" width="90%"></center> 
<br>
>**NOTE:** The `[1]` we're seeing to the left of the output is an "index" number. As we go further we'll see why these are useful, but for now just know it's a counter for how many iterms are returned by a command – for each row of printed output it lists the "index" number of the first item of that row, here we only have 1 row of output, so it's always just showing 1.  

<br>

---
<br>
# Variables in R
Most of the time R acts on things that are stored in variables. In R, things that are stored in variables are referred to as "objects". Objects can be of different types, and the type we are working with, and the type of data within an object, determines what you can and can't do with that object. This can be super-confusing at first, at least it was for me, so don't worry about it too much right now. We'll cover this a little bit moving forward here, but mostly as you start to spend more time bumping into error messages and googling what's up, you'll find that these concepts pretty quickly come into focus.

In R, we assign values to variables and name them using the assignment operator ( `<-` ) as shown in the following code block. Here we're naming the variable "x" (this name could be anything), and giving it a value of 4.  

```R
x <- 4
```

After executing the command, our "environment" pane should now show this variable there.  

Now that we have the value "4" stored in the variable "x", we can use the variable name in functions. Here's some examples doing the same calculations we performed above, but now with our variable.

```R
x
x + 4
x / 2
x * x
2 ^ x
```

<center><img src="{{ site.url }}/images/R_var_arith.png" width="90%"></center> 
<br>
We can also check what type of data is contained within this variable with the `class()` function:

```R
class(x)
```

And we find out it is of class "numeric". Let's try storing a different type, like a word:

```R
w <- "europa" 
class(w)
```

>**NOTE:** Here, the class is "character". Notice again that we need to use quotes when working with characters, like when we set the working directory above. Without the quotes, R will try to find a variable with that name, rather than treating it as plain text. This isn't the case for numbers, like when we set "x" to 4 above (this is also why you can't name a variable starting with a number in R). 

We also can store multiple items into a single variable. A one-dimensional object holding multiple items that are of the same "type" (e.g. numeric, or character) is known as a *vector*. To put multiple items into one object, we can use the `c()` function, stemming from "concatenate". Here we'll make a vector of numbers:

```R
y <- c(5, 6, 7)
y
```

<center><img src="{{ site.url }}/images/R_vector.png" width="90%"></center> 
<br>
Note that this is still of class "numeric" by checking with `class(y)`. It's good practice to get used to actively being aware of what type of objects we are working with. 

Variables can also hold tables. Here we're going to make another vector with 3 numbers, and then combine it with our previous vector in order to make what's known as a *dataframe* (Again, don't worry about remember all this terminology right away! This is just about exposure right now.) We'll do this with the `data.frame()` function, creating a variable called "our_table".

```R
z <- c(8, 9, 10)
our_table <- data.frame(y, z)

our_table

class(our_table)
```

<center><img src="{{ site.url }}/images/R_table.png" width="90%"></center> 
<br>
Dataframes are two-dimensional objects of rows and columns. Here we can see that the default behavior of the `data.frame()` function took our two vectors and put them in a table where each original vector now represents one column. Another similar, but distinct, table structure in R is a "matrix". You will sometimes find you need to convert a dataframe to a matrix or vice versa depending on what you are trying to do with it. Keep this in mind as one of the things to look at first when you run into an error.  

<br>

---
<br>

# The wonderful world of indexing
One of the most powerful things about R is how easy it makes it to subset vectors and tables down to whatever you are interested in via what's known as "indexing". Here we'll look at a couple of the ways we can specify what we would like to subset, and we'll see these in practice on a larger scale below. 

## Subsetting by position
Looking back at our vector stored in variable `y`, it contains 3 values: 5, 6, and 7. These values exist in the object in this order as positions 1, 2, and 3 of the variable `y`. One way we can subset specific values involves using this position information – this position information for each value is known as that value's "index". If we specify the vector name, and then put in brackets `[ ]` the position(s) we are interested in, R will return the value(s) (things following the `#` are "comments" that the program ignores).

```R
y # the whole vector
y[1] # the first item
y[2] # second item
y[3] # third item
```

<center><img src="{{ site.url }}/images/R_vec_index.png" width="90%"></center> 
<br>

>**NOTE:** It's good to think about a way to read this syntax that makes sense to us. The variable we are subsetting from comes first, "y" above, then within brackets we are stating what parts of it we want. Here just by index number, so we're saying something like 'from object "y", give us the first item' (`y[1]`). 


We can also ask for multiple by using the `c()` function we saw above:

```R
y[c(1,3)] # specifying items 1 and 3
```

<center><img src="{{ site.url }}/images/R_vec_index_2.png" width="90%"></center> 
<br>

Ok, so that's how we can subset by saying which positions we want. But in practice we often won't actually know which positions of a vector hold the values we are interested in – meaning we usually won't know the "index" number needed to pull out a specific value. This is where another type of indexing comes into play.

## Subsetting by conditional statements
Another way to subset via indexing in R makes use of conditional statements. For example let's say we wanted all values of our vector `y` that were greater than or equal to 6. We can subset those values by putting `y >= 6` within our subsetting brackets like so:

```R
y # the whole vector
y[y >= 6] # returns just the last two values
```

<center><img src="{{ site.url }}/images/R_vec_index_ge_6.png" width="90%"></center> 
<br>

The way I read the expression `y[y >= 6]` in my head is: "Give me all the values of `y` where `y` is greater than or equal to 6." What R is actually doing under the hood here goes a little further into the weeds than would be helpful right now, so for the moment we're going to move on. But understanding fundamentally how subsetting with conditional statements works is extremely powerful and important, so feel free to dive deeper into this when you can by visiting the [going deeper with indexing page](/R/more_indexing){:target="_blank"} 🙂   

One last thing we're going to introduce here is the `!` character, which inverts the interpretation of conditional expression provided. As we've seen, `y[y >= 6]` will return all values within 'y' that are greater than or equal to 6. But if we add in the `!` point, it will return the opposite for us:

```R
y # whole vector

y[y >= 6] # returns only 6 and 7
y[!y >= 6] # returns only 5
```

<center><img src="{{ site.url }}/images/R_vec_index_ge_6_not1.png" width="90%"></center> 
<br>

The use of the `!` character like this may seem a little unnecessary in the case of strictly numerical conditional expressions like this, but it's very handy for other types of conditional statements. We'll see a somewhat more complicated example below where inverting the `!` logical vector is the only way to actually get at what we want. Again, for now we are glancing over the fundamentals of how R handles indexing with conditional statements like this, but this underlies so many of the ways that we parse down data in R to get what we want, so if you'd like to get into it a little more, be sure to visit [going deeper with indexing](/R/more_indexing){:target="_blank"} when you can 🙂  

So far we've been dealing with subsetting just one-dimensional vectors, but similar rules apply to subsetting two-dimensional tables.

## Subsetting tables
As we've seen, vectors are one-dimensional objects, so when we want to subset from one we only need to specify details for one coordinate for which item(s) we want. But tables are 2-dimensional objects, so we need to provide instructions for handling two coordinates (one for which rows we'd like and one for which columns). In R, this is still done with the same subsetting brackets ( `[ ]` ), but now providing two values within them separated by a comma. **The first value we enter in the brackets specifies which *rows* you'd like, and the the second value (separated by a comma) specifies which *columns*.** Using the table we made above, stored in the variable "our_table", let's run through some examples of what this looks like:

```R
our_table # whole table

our_table[2, 2] # subset value in the second row and second column only
```

<center><img src="{{ site.url }}/images/R_tab_index_1.png" width="90%"></center> 
<br>

If we provide nothing for either the row or the column position, we will get all values for it: 
  
```R
our_table[ , 2] # subset all rows, but only the second column

our_table[3, ] # only row 3, but both columns
```

<center><img src="{{ site.url }}/images/R_tab_index_2.png" width="90%"></center> 
<br>

Notice that when subsetting returns only one column, but multiple rows (as in the first example there, `our_table[ , 2]`), it returns a numeric *vector*. But when subsetting returns one row, but multiple columns (as in the second example there, `our_table[3, ]`), it returns a dataframe:

```R
class(our_table[ , 2])

class(our_table[3, ])
```

<center><img src="{{ site.url }}/images/R_tab_index_3.png" width="90%"></center> 
<br>

This hints at something fundamental about R – that it treats rows and columns differently. This is another detail we don't need to worry about remembering, but just having seen it once may help troubleshoot faster if we happen to run into it sometime 🙂  

If we want, we can tell R to retain the dataframe structure by adding the optional argument `drop=F` like so:

```R
our_table[ , 2]
class(our_table[ , 2])

our_table[ , 2, drop=F]
class(our_table[ , 2, drop=F])
```

<center><img src="{{ site.url }}/images/R_tab_index_4.png" width="90%"></center> 
<br>

Another way we can pull out a specific column from a dataframe as a vector is by the column header/name, in this case we have 2 columns with the names 'y' and 'z'. The function `colnames()` can tell us this:

```R
colnames(our_table)  
```

<center><img src="{{ site.url }}/images/R_tab_index_5.png" width="90%"></center> 
<br>

To specify a column we want to pull from a dataframe, we enter the table variable name (here, "our_table), followed by a `$`, followed by the column name we want:

```R
our_table$z
```

<center><img src="{{ site.url }}/images/R_tab_index_6.png" width="90%"></center> 
<br>

We can also do this with the bracket format of subsetting, and therefore combine it with rows by index. Here we are saying we want rows 2 and 3, and specify the column by name instead of its index number:

```R
our_table[c(2,3), "z"]
```
<center><img src="{{ site.url }}/images/R_tab_index_7.png" width="90%"></center> 
<br>


Indexing in R can definitely seem pretty confusing at first, especially the logical indexing covered in the [going deeper with indexing page](/R/more_indexing){:target="_blank"}, but it's also one of the things that makes R so awesome! 

<br>

---
<br>
# Reading in and writing out data
Most of the time when working with R you're going to want to read in some data, do some stuff to it, and then write out something else to a new file that will then go on to live a wonderous and full life beyond the R environment. Here we're going to cover the basics of reading in and writing out files. 

## Checking out the data in the terminal first
Before we try to read data into R, it's a *really* good idea to know what we're expecting. Let's get some idea of what our example file, "gene_annotations.txt", looks like in the terminal with some of the tools introduced in the [*bash* intro](/bash/bash_intro_binder){:target="_blank"} page.  

We can work at the terminal in RStudio too, if we click the "Terminal" tab at the top of the "source" pane (which is the bottom left one in our binder environment):

<center><img src="{{ site.url }}/images/R_basics_terminal.png" width="90%"></center> 
<br>

>**NOTE:** If there is a conda error message that pops up before the prompt appears like shown in the image above, we can ignore that. 

So in our terminal window, let's change into the directory holding our example file, and then take a peek at it with `less` first:

```bash
cd ~/R_basics_temp/

less -S gene_annotations.txt # the `-S` prevents lines from wrapping
```

<center><img src="{{ site.url }}/images/R_basics_terminal_less.png" width="90%"></center> 
<br>

From this we can see that it's a tab-delimited file, and that it has a header with column names for each column. We can exit `less` by pressing the **`q`** key. 

Let's take a look just at the column names:

```bash
head -n 1 gene_annotations.txt
```

<center><img src="{{ site.url }}/images/R_basics_terminal_head.png" width="90%"></center> 
<br>

We can also quickly check how many rows we should be expecting:

```bash
wc -l gene_annotations.txt
```

<center><img src="{{ site.url }}/images/R_basics_terminal_wc_l.png" width="90%"></center> 
<br>

Ok. So now instead of being blind to what the file holds, we know that it's tab-delimited, it has a header with column names, and it has 8 columns and 84,785 rows (including the header). Awesome. There are some parameters we need to set when we read a file into R, and know these things will help us check to make sure things are working like we want. Now let's get it into R! 

Be sure to switch back to the "Console" tab at the bottom left now, away from the "Terminal" tab, so that our pane looks like this again:

<center><img src="{{ site.url }}/images/R_tab_index_7.png" width="90%"></center> 
<br>

## read.table()
One of the most common ways of reading tables into R is to use the `read.table()` function. To start, let's try reading our "gene_annotations.txt" table into R with no arguments other than specifying the file name:

```R
gene_annotations_tab <- read.table("gene_annotations.txt")
```

<center><img src="{{ site.url }}/images/read_table_err.png" width="90%"></center> 
<br>

Yay our first error! Many error messages may seem a little cryptic at first, but you'll be surprised at how many of them magically start to make sense over time. The important part in this one is at the end where it says "line 1 did not have 22 elements". We know from our exploration in the terminal above that our table should have 8 columns. This is a sign there is something up with how R is trying to split each line into columns. 

If we take a look at the help menu for this function with `?read.table`:

```R
?read.table
```
The help shows up in our bottom right pane. And scanning through there for anything about specifing the delimiter, we can find the argument "sep". And it seems that by default the "sep" argument is set to act on all white space, which includes tabs AND blank spaces:

<center><img src="{{ site.url }}/images/read_table_help_sep.png" width="90%"></center> 
<br>

If we remember looking at our "gene_annotations.txt" file in the terminal with `less`, in addition to it being tab-delimited, there were also spaces within the KO and COG annotation columns. 

<center><img src="{{ site.url }}/images/R_basics_terminal_less.png" width="90%"></center> 
<br>  

So `read.table()` by default is making a new column everywhere there is a space, and then coming back to us and saying "Hey, your first line doesn't have all the columns it should have based on the rest of your file!" Which is nice of it, because it's letting us know something is probably wrong 🙂

Let's try running the command again, but this time stating that the delimiter should only be tabs (tab characters are specified with an backslash followed by a "t" like so: `\t`.

```R
gene_annotations_tab <- read.table("gene_annotations.txt", sep="\t")
```

This works without any errors, let's take a look at it with the `head()` function in R:

```R
head(gene_annotations_tab)
```

<center><img src="{{ site.url }}/images/R_basics_read_table_head.png" width="90%"></center> 
<br>

We can ignore that things are wrapping a little funny because it's wider than the panel can allow, but it put our column names in the first row and added new column names ("V1", "V2", etc.). Looking at the help menu for `read.table()` some more in our bottom right pane, we find there is an argument for "header", which is by default set to `FALSE`. So let's try again but this time we'll specify that there is a header:

```R
gene_annotations_tab <- read.table("gene_annotations.txt", sep="\t", header=TRUE)

head(gene_annotations_tab) 
```

<center><img src="{{ site.url }}/images/R_basics_read_table_head_2.png" width="90%"></center> 
<br>

That looks better. Let's also check our column names and the size of the table:

```R
colnames(gene_annotations_tab)

dim(gene_annotations_tab)
```

<center><img src="{{ site.url }}/images/R_basics_read_table_dim.png" width="90%"></center> 
<br>

>**NOTE:** Now that our vector of column names is longer than the window, our index numbers are printed on the left for each row ("[1]", "[4]", "[7]"). That is the index (positional number) of each row's first item.

So our table is 84,784 rows by 8 columns, which is great as that's what we expect based on what we saw when investigating the file in the terminal. 

Now let's generate a new table so we can practice writing out to a file from R. You may have noticed there are some NAs in our "gene_annotations_tab.txt" table, which are special values to R. These are present in the KEGG and COG annotation and ID columns as "NA" for those genes which weren't annotated. Here, let's pretend we want to subset our full table down to include only those genes that *were* annotated by KEGG. R's `is.na()` function can help us do this. The `is.na()` function will return whether or not each item of an object contains an `NA`, but in this case we are interested in those that are *not* `NA`, meaning we want those that actually contain values (KEGG identifiers in this case). So we need to return the opposite of this using the `!` character like we did above with `y[!y >= 6]`.  

This combines a few concepts, so let's run the code and then we'll break it down 🙂

```R
KEGG_only_tab <- gene_annotations_tab[!is.na(gene_annotations_tab$KO_ID), ]
```

Let's look closely at this first line here. We start off with the name of the new table we're creating "KEGG_only_tab", then we have our assignment character `<-` , and then the part that is actually specifying what we want from the original table. Same as above with subsetting our smaller table, we are giving the variable name that holds the table we want to subset from ("gene_annotations_tab"), followed by the subsetting brackets ( `[ ]` ) that enclose 1) which rows we want to subset *before* the comma, and 2) which columns we want to subset *after* the comma. In this case we have nothing after the comma which means we are taking all columns, but the rows position in these brackets is where the magic is happening. With R it can sometimes be easier to read from the inside and work our way out. Looking at it that way, the inner parentheses enclose `gene_annotations_tab$KO_ID`. If we were to run `gene_annotations_tab$KO_ID` by itself, it would print all 84,784 values within that column as a vector. Wrapping that with the function `is.na()` like we're doing here will return whether or not each item is an `NA` or not. **This would by itself give us a table of all rows where the column KO_ID had an `NA`. But to get all of those that were annotated by KEGG, we flip this with the `!` character in front of the function.**  

And if we peek at our new table with `head()`, we se all top 6 have KEEG annotations, where as before some where NA:

```R
head(KEGG_only_tab)
```

<center><img src="{{ site.url }}/images/R_basics_kegg_head.png" width="90%"></center> 
<br>


And we can also look at how many genes we dropped that didn't have a KEGG annotation assigned to them:

```R
dim(gene_annotations_tab) # 84,784 genes
dim(KEGG_only_tab) # 37,319 had KEGG annotations assigned
```

<center><img src="{{ site.url }}/images/R_basics_kegg_dim.png" width="90%"></center> 
<br>

As mentioned above, there is a little more going on under the hood here that is covered at the [going deeper with indexing page](/R/more_indexing){:target="_blank"}. But it's definitely a little strange at first. Now let's write out our new table so it is free from the constraints of R!

## write.table()

Now, let's write out our new table of only those genes that were annotated by KEGG to a new tab-delimited file called "all_KEGG_annotated_genes.txt". We can do this with the `write.table()` function. If we glance at the help menu for this with `?write.table`, we see that the default delimiter (what separates the columns) seems to be a blank space, so we need to be sure to specify that we instead want it to be tab-delimited by providing the argument `sep="\t"`. We also don't want to keep the row names R added, so we need to set `row.name=F`. And by default it will also write out quotation marks around character strings (like our annotation columns) which we also don't want, so we can set `quote=F`. How would we know all these things at first? We wouldn't, ha. Just like above when we were reading in the file, we can write it out and check if things look how we want, and then look up how to get what we want (this is typical to have to do for me even though I use R pretty regularly). So let's add in these additional arguments to make the file to our liking.

```R
write.table(KEGG_only_tab, "KEGG_annotated.tsv", sep="\t", row.names=F, quote=FALSE)
list.files() # checking it is there now
```

<center><img src="{{ site.url }}/images/R_basics_kegg_write_out.png" width="90%"></center> 
<br>

And as mentioned, it's good practice to peek at the output in the terminal when we are configuring the options to write something out to make sure it's doing what we think it's doing. So we can switch back to the "Terminal" tab in the bottom left pane, and check our new file with `less`:

```bash
less -S KEGG_annotated.tsv
```

<center><img src="{{ site.url }}/images/KEGG_only_tab_less.png" width="90%"></center> 
<br>

---
---
<br>
<h1>Congrats on getting through the basics of R!</h1>
R was not immediately intuitive to me, but it is extremely powerful and *many* statistical tools are designed to work within it. So it is well worth the time getting to know it if working with big data or running complex statistical analyses are part of your work 🙂
