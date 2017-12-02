---
layout: main
title: R basics
categories: [R, tutorial]
permalink: /R/basics
---

{% include _R_basics_toc.html %}

{% include _side_tab_R.html %}

This module is designed for those that are either completely new to R, or have some experience but maybe don't feel as solid about some of the fundamentals as they'd like. It will run through the very basics such as setting up your working environment, assigning variables, "indexing" (subsetting data), reading in and writing out data, and installing packages. Some relevant terminology is presented [here](/R/index) if you find yourself seeing some words that are unfamilar to you. This page is meant to be a quick-start to get you into and using the R environment. There is, of course, an incredible amount of functionality we won't be touching on here, so be sure to peruse the extensive [R intro documentation available here](https://cran.r-project.org/doc/manuals/r-release/R-intro.html) when you're starting to get your bearings and want to go further, and there are also some great courses available at [DataCamp](https://www.datacamp.com/) that are worth checking out after you have a solid foundation. 
<br>
<br>

---
---
<br>

# Installation
It is possible your computer already has R, if you are unsure, you can check by opening a terminal window and typing `R`. If this launches R rather than giving an error message, you should be good to go (enter `q()` to exit the R environment). If you do not have R, you can download it from here for Mac: [https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/). And if you have a relatively newer Mac, you may also need to install XQuartz which you can get from here: [https://www.xquartz.org/](https://www.xquartz.org/). 

Lastly, I highly, *highly*, **highly** recommend installing the free version of [RStudio](https://www.rstudio.com/) if you don't already have it. RStudio is an interface for R that not only makes everything you will do in R easier and more organized, but it's also invaluable for reproducibility of your analyses as it makes it second-nature to generate and save R scripts of everything you're doing, while you're doing it – which is very helpful when you want to look back and see what worked out of the 20 things you just tried 🙂 . You can download an RStudio installer from [here](https://www.rstudio.com/products/rstudio/download/#download).  
<br>

---
<br>
# RStudio layout
RStudio has 4 main panes: 1) console; 2) source; 3) Environment, History, Connections; and 4) Files, Plots, Packages, etc. The "console" is where you can run commands just as though you were working in an R environment at a command line, and it is also where results will print out. The "source" pane acts as a sort-of interactive text document within which you can write out and save all of your commands, and then call them when you'd like. This is one of the reasons R Studio is great, you're constantly building up your notes file as you work, without any added effort needed. If you want to run a command written in the source file, while on the line of the command press `Cmd + Enter`. The "Environment/History/Connections" pane displays all of the variables and data structures you currently have stored. And the last pane, "Files/Plots/Packages/etc.", allows you to peruse your computer in the typical Finder fashion, displays any plots you generate, and serves as your help window. To get help on a function in R, just place a `?` in front of the function name, and leave off the trailing `()`. For example, try entering `?getwd`.  
<br>

---
<br>
# Some practice data
As with the [bash basics](/bash/basics) module, if you'd like to follow along with this tutorial, copy and paste the following commands into your terminal to get set up with a small, temporary working directory that has the files I'll be working with. (If you're unsure of what these mean, you should run through [that page](/bash/basics) and/or just check out [here](/bash/basics#bottom) for some explanation.)

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
# Setting up your working environment
Just like when working at the command line in *bash*, we need to be aware of "where" we are in our computer when working in R. The `getwd()` and `setwd()` help us do this in R. Commands in R typically take this structure, with the command being followed by parentheses. Inside those parentheses is where the arguments to the command would go if there are any. In the case of `getwd()` no arguments are needed as it is just going to tell you where you currently are. But in the case of `setwd()`, we then need to tell it where we want to go. 

```R
  # check where you are in R by default when it boots up
getwd()
  # and let's set our working directory to where we want it to be (where we put our practice data)
setwd("~/R_basics_temp/")
```

Now that we are in the correct location that contains the files for our tutorial here, let's check that that are actually there. In *bash* this would be done with the `ls` command, here in R we do it with the `list.files() function.

```R
list.files()
```

And this should reveal our practice file, "gene_annotations.txt". Just like when working at the command line, in R you can "point" to other places in your computer by providing the relative or absolute path to whichever file/directory you're looking for in the same manner as described [here](/bash/basics#moving-around).  
<br>

---
<br>
# Basic calculations
At its core, R is just a big calculator. So we can do baseline arithmetic operations like the following examples. (Don't forget if you're entering these in the "source" pane of RStudio, you need to press `Cmd + Enter` to execute the command.)

```R
4 + 4
4 / 2
4 * 4
2 ^ 4 # the ^ symbol is for raising to a power
```

<center><img src="{{ site.url }}/images/R_basic.png"></center> 
<br>
And, just to note, the `[1]` you're seeing to the left of the output is an "index" number. As we go further you'll see why these are useful, but for now just know it's a counter for how many iterms are returned by a command, and that for each row of output it lists the "index" number of the first item of that row – and here we only have 1 row of output.  
<br>

---
<br>
# Assigning variables
Most of the time R acts on things that are stored in variables; in R, these specific things that are stored in variables are referred to as "objects". These objects can be of different types, and the type of object you are working with, and the type of data within an object, determines what you can and can't do with that object. This should all seem pretty nebulous at first, at least it was for me, so don't worry about it too much right now. We'll cover this a little bit moving forward here, but mostly as you start to spend more time bumping into error messages and googling what's up, you'll find that these concepts pretty quickly come into focus.

In R, we assign values to variables and name them using the assignment operator ( `<-` ) as shown in the following code block. Here we're naming the variable "x" (arbitrarily), and giving it a value of 4. (Feel free to do this in the "console" or in the "source"
pane.)

```R
x <- 4
```

When you execute the command, it should print in the console pane, and now in your "environment" pane you should see there is a variable there called "x". For a note on syntax, it doesn't matter if you put spaces or not between the variable name, the assignment operator ( `<-` ), and the value you're assigning to it. But generally it makes the code easier to read if you do, and it is good practice to make things as easy to read as possible for the sake others (and yourself) who may be looking at it in the future. 

Now that we have the value "4" stored in the variable "x", we can use the variable name in functions. Here's some examples doing the same calculations we performed above, but now with our variable.

```R
x
x + x
x / 2
x * x
2 ^ x
```

<center><img src="{{ site.url }}/images/R_var_arith.png"></center> 
<br>
We can also check what type of data is contained within this variable with the `class()` function:

```R
class(x)
```

And we find out it is of class "numeric". Let's try storing a different type, like a word:

```R
word <- "europa" 
class(word)
```

Here, the class is "character". It's important to understand why we used quotes here, so let's take a second and look at it. Above, "europa" is a known as a string of characters, and we need to tell that to R by surrounding it in quotes like we did. If we didn't provide the quotes, R would be looking for a variable named `europa`, and it would give us an error saying it wasn't found, give it a shot without the quotes to see. 

We also can store multiple items into a variable. A one-dimensional object holding multiple items that are of the same "type" (e.g. numeric, character) is known as a vector. To put multiple items into one object, we can use the `c()` function, stemming from "concatenate". Here we'll make a vector of numbers:

```R
y <- c(5, 6, 7)
y
```

<center><img src="{{ site.url }}/images/R_vector.png"></center> 
<br>
Note that this is still of class "numeric" by checking with `class(y)`. It's good practice to get used to actively being aware of what type of objects you are working with. 

Variables can also hold tables. Here we're going to make another vector with 3 numbers, and then combine it with our previous vector in order to make what's known as a dataframe. We'll do this with the `data.frame()` function, creating a variable called "our_table".

```R
z <- c(8, 9, 10)
our_table <- data.frame(y, z)
  # let's take a look at it
our_table
  # and let's see what class it is
class(our_table)
```

<center><img src="{{ site.url }}/images/R_table.png"></center> 
<br>
Dataframes are two-dimensional objects of rows and columns. Here we can see that the default behavior of the `data.frame()` function took our two vectors and put them in a table where each original vector now represents one column. Another similar, but distinct, table structure in R is a "matrix". You will sometimes find you need to convert a dataframe to a matrix or vice versa depending on what you are trying to do with it. Keep this in mind as one of the things to look at first when you run into an error.  
<br>

---
<br>
# The wonderful world of indexing
One of the most powerful things about R is how easy it makes it to subset vectors and tables down to whatever you are interested in via what's known as "indexing". Here we'll look at a couple of the ways we can specify what we would like to subset, and we'll see these in practice on a larger scale below. 

## Subsetting by position
Looking back at our vector stored in variable `y`, it contains 3 values: 5, 6, and 7. These values exist in the object in this order as positions 1, 2, and 3 of the variable `y`. One way we can subset specific values involves using this position information – this position information for each value is known as that value's "index". If we specify the vector name, and then put in brackets `[ ]` the position(s) we are interested in, R will return the value(s).

```R
y # the whole vector
y[1] # the first item
y[2] # second item
y[3] # third item
y[c(1,3)] # specifying items 1 and 3 using the c() function
```

<center><img src="{{ site.url }}/images/R_vec_index.png"></center> 
<br>
Ok, so that's how we can subset by saying which positions we want. But in practice we often won't actually know which positions of a vector hold the values we are interested in – meaning we usually won't know the "index" number needed to pull out a specific value. This is where another type of indexing comes into play.

## Subsetting by conditional statements

Another way to subset via indexing in R makes use of conditional statements. For example let's say we wanted all values of our vector `y` that were greater than or equal to 6. We can subset those values by putting `y >= 6` within our subsetting brackets like so:

```R
y # the whole vector
y[y >= 6] # returns just the last two values
```

<center><img src="{{ site.url }}/images/R_vec_index_TF5.png"></center> 
<br>
The way I read the expression `y[y >= 6]` in my head is: "Give me all the values of `y` where `y` is greater than or equal to 6." What R is actually doing under the hood here goes a little further into the weeds than would be helpful right now, so for the moment we're going to move on. But understanding fundamentally how subsetting with conditional statements works is extremely powerful and important, so feel free to dive deeper into this when you can by visiting the [going deeper with indexing page](/R/more_indexing).

One last thing we're going to introduce here is the `!` character, which inverts the interpretation of `TRUE` and `FALSE`. As we've seen, `y[y >= 6]` will return all values within 'y' that are greater than or equal to 6. But if we add in the `!` point, it will return the opposite for us:

```R
y # whole vector

y[y >= 6] # returns only 6 and 7
y[!y >= 6] # returns only 5
```

<center><img src="{{ site.url }}/images/R_vec_index_TF6_not1.png"></center> 
<br>
The use of the `!` character like this may seem a little unnecessary in the case of strictly numerical conditional expressions like this, but it's very handy for other types of conditional statements. We'll see a somewhat more complicated example below where inverting the `!` logical vector is the only way to actually get at what we want. As mentioned, we are glancing over the fundamentals of how R handles indexing with conditional statements like this for now, but this underlies so many of the ways that we parse down data in R to what we want, and it's going to help immensely if you understand the fundamentals of it, so be sure to get to [going deeper with indexing](/R/more_indexing) when you can!

So far we've been dealing with just one-dimensional vectors, but the same rules apply to working with two-dimensional tables.

## Subsetting tables

As we've seen, vectors are one-dimensional objects, so when we want to subset from one we only need to specify details for one coordinate for which item(s) we want. But tables are 2-dimensional objects, so we need to provide instructions for handling two coordinates (one for which rows we'd like and one for which columns). In R, this is still done with the same subsetting brackets ( `[ ]` ), but now providing two values within them separated by a comma. **The first value we enter in the brackets specifies which *rows* you'd like, and the the second value (separated by a comma) specifies which *columns*.** Using the table we made above, "our_table", let's run through some examples of what this looks like:

```R
  # whole table
our_table

  # subset the value in the second row and second column only
our_table[2, 2]

    ### if we provide nothing for either the row position or the column position, we will get all values
  
  # subset all rows, but only the second column
our_table[ , 2] # notice when subsetting returns only one column, it returns a vector

  # only row 3, but both columns
our_table[3, ] # notice when subsetting returns one row, but more than one column, it still returns a dataframe, these intricacies will become important when things are wrapped up in more code

  # we can force the original way to retain the table structure by adding the optional argument drop=F, like so:
our_table[ , 2, drop=F]

   ## another way we can pull out a specific column from a dataframe as a vector is by the column header/name, in this case we have 2 columns with the names 'y' and 'z'
colnames(our_table)  

  # to specify a column we want to pull from a dataframe, we enter the table variable name, followed by a `$`, followed by the column name or names of interest
our_table$z
  
    ## we can also do this with the bracket format of subsetting
  # here, we are still saying we want all rows, but instead of using the number 2 to specify the position of the column we want, we are going to use the name `z`
our_table[ , "z"]
```

<center><img src="{{ site.url }}/images/R_tab_index.png"></center> 
<br>
Indexing in R can seem pretty confusing at first, especially the logical indexing covered in the [going deeper with indexing page](/R/more_indexing), but it's also one of the things that makes R so awesome. Like everything, exposure and practice with the process will make you much more comfortable, and as you begin working with large tables and combining subsetting metrics you'll quickly begin to see how useful this fundamental system is. 
<br>

---
<br>
# Reading in and writing out data
Most of the time when working with R you're going to want to read in some data, do some stuff to it, and then write out something else to a new file that will then go on to live a wonderous and full life beyond the R environment. Here we're going to cover the basics of reading in and writing out files. 

## Checking out the data in the terminal first
Before we try to read data into R, it's a *really* good idea to know what we're expecting. Let's get some idea of what our practice "gene_annotations.txt" file looks like in the terminal with some of the skills we picked up from the [*bash* basics](/bash/basics#working-with-plain-text-files-and-directors) page. 

So go back to your terminal window where you downloaded the practice data for the moment, and let's take a peek with `less -S gene_annotations.txt`: 

<center><img src="{{ site.url }}/images/terminal_table_less.png"></center> 
<br>
From this we can see that it's a tab-delimited file, and that it has a header with column names for each column. Let's take a look just at the column names:

<center><img src="{{ site.url }}/images/R_terminal_colnames.png"></center> 
<br>
And another easy thing to check is how many rows we should be expecting:

<center><img src="{{ site.url }}/images/R_terminal_lines.png"></center> 
<br>
Ok. So now instead of being blind to what the file holds, we know that it's tab-delimited, it has a header with column names, and it has 8 columns and 84,785 rows (including the header). Awesome. Now let's get it into R! 

## read.table()

Probably the most common method of reading in tables is to use the `read.table()` function. To start, let's try reading our gene_annotations.txt table into R with no arguments other than specifying the file name:

```R
gene_annotations_tab <- read.table("gene_annotations.txt")
```

<center><img src="{{ site.url }}/images/read_table_err.png"></center> 
<br>
Yay our first error! Many error messages may seem a little cryptic at first, but you'll be surprised at how many of them magically start to make sense over time. The important part in this one is at the end where it says "line 1 did not have 22 elements", and we know from our exploration in the terminal above that our table should have 8 columns. This is a sign there is something up with how R is trying to split each line into columns. 

If we take a look at the help menu for this function with `?read.table`, and scan for anything about specifing the delimiter, we can find the argument "sep". And it seems that by default the "sep" argument is set to act on all white space, which includes tabs AND blank spaces:

<center><img src="{{ site.url }}/images/read_table_help_sep.png"></center> 
<br>

If we switch back to our terminal and look at our gene_annotations.txt file again with `less -S gene_annotations.txt`, we see that in addition to it being tab-delimited, there are also spaces within the KO and COG annotation columns. 

<center><img src="{{ site.url }}/images/terminal_table_less.png"></center> 
<br>
So `read.table()` by default is making a new column everywhere there is a space, and then coming back to us and saying "Hey, your first line doesn't have all the columns it should have based on the rest of your file. Get lost."

Let's try running the command again, but this time specifying the delimiter should only be tabs (tab characters are specified with an backslash followed by a "t" like so: `\t`.

```R
gene_annotations_tab <- read.table("gene_annotations.txt", sep="\t")
  # and let's take a peek at the table with head():
head(gene_annotations_tab) 
```

That worked, but it put our column names in the first row and added new column names ("V1", "V2", etc.). Looking at the help menu for `read.table()` some more we find there is an argument for "header", which is by default set to `FALSE`. So let's try again but this time we'll specify that there is a header:

```R
gene_annotations_tab <- read.table("gene_annotations.txt", sep="\t", header=TRUE)
  # now let's take a peek at the table with head():
head(gene_annotations_tab) 
  # and check what our column names are:
colnames(gene_annotations_tab)
  # and check the size of the table:
dim(gene_annotations_tab)
```

<center><img src="{{ site.url }}/images/read_table_dim.png"></center> 
<br>
So our table is 84,784 rows by 8 columns, which when not counting the header is exactly what we'd expect from our `wc -l` in the terminal. Notice that since the table doesn't fit all the way across my "console" pane it wraps twice when we viewed it with the `head` command. Also notice that R by default gave it row names that are numbered. This is usually fine, but in this case since our gene_IDs are so similar, we're going to rename the row names to match the gene_ID for each row (to avoid any possible confusion when glancing at the table). We can view and set the row names of our table like with the `row.names()` function. As we saw above, there are a lot of rows here, so we're going to wrap this function with the `head()` function to keep the output manageable.

```R
  # first let's look at what they are currently, going to "wrap
head(row.names(gene_annotations_tab)) # right now they start at 1
  # let's see where they end with the tail() function
tail(row.names(gene_annotations_tab)) # not surprisingly they end at 84,784 - the total number of rows we had

    ## now let's set the row names to equal what the gene_ID is for each row
  # we do this with the assignment character like we saw above, and we'll specify the gene_ID column using the `$` notation we saw earlier
row.names(gene_annotations_tab) <- gene_annotations_tab$gene_ID

  # now let's check the start and end of the row names again:
head(row.names(gene_annotations_tab)) # starts with 0 now
tail(row.names(gene_annotations_tab)) # ends with 84,786 (so some genes have apparently been filtered out at some point, this could have possibly caused confusion had we left the row names as straight numbers from 1-84,784)
```

<center><img src="{{ site.url }}/images/new_table.png"></center> 
<br>
It's also a good idea to run the `str()` function (structure) on a table when you first read it in to make sure the columns are all the appropriate classes – i.e. in this case, `str(gene_annotations_tab)`. When we run that on our table, it tells us the object is a dataframe, its size, and lists each column by name following the `$` character with information about that column. We can see there are many columns that are labeled as "Factors": 

<center><img src="{{ site.url }}/images/str.png"></center> 
<br>
Another default setting of `read.table()` is to convert columns with strings in them (text characters) to factors (categorical text values). Sometimes this is fine, but sometimes it's not what you want. Be sure to pay attention to it. You can change the default behavior of `read.table()` by adding the `stringsAsFactors=FALSE` argument.

<a id="cond_example"></a>

Now let's generate a new table so we can practice writing out to a file from R. You may have noticed there are some NAs in our "gene_annotations_tab" table, which are special values to R. These are present in the KEGG and COG annotation and ID columns as "\<NA\>" for those genes which weren't annotated. Here, let's pretend we want to subset our full table down to include only those genes that *were* annotated by KEGG. R's `is.na()` function can help us do this. The `is.na()` function will return whether or not each item of an object contains an `NA`, but in this case we are interested in those that are *not* `NA`, we want those that actually contain values (KEGG identifiers in this case). So we need to return the opposite of this using the `!` character like we did above with `y[!y >= 6]`.  

This combines a few concepts, so let's run the code and then we'll break it down.

```R
all_KEGG_gene_annotations_tab <- gene_annotations_tab[!is.na(gene_annotations_tab$KO_ID), ]

  # if we peek at the table with head() we see all top 6 have KEGG annotations, whereas before some were NA in the head() output
head(all_KEGG_gene_annotations_tab)

  # we can also look at how many genes we dropped that didn't have a KEGG annotation assigned
dim(gene_annotations_tab) # 84,784 genes
dim(all_KEGG_gene_annotations_tab) # 37,319 had KEGG annotations assigned
```

<center><img src="{{ site.url }}/images/KEGG_only_tab.png"></center> 
<br>
Let's look closely at what's going on in that first line here. We start off with the name of the new table we're creating "all_KEGG_gene_annotations_tab", then we have our assignment character `<-` , and then the part that is actually specifying what we want from the original table. Same as above we are giving the variable name that holds the table we want to subset from ("gene_annotations_tab"), followed by the subsetting brackets ( `[ ]` ) that enclose 1) which rows we want to subset *before* the comma, and 2) which columns we want to subset *after* the comma. In this case we have nothing after the comma which means we are taking all columns, but the rows position in these brackets is where the magic is happening. With R it can sometimes be easier to read from the inside and work your way out. Looking at it that way, the inner parentheses enclose `gene_annotations_tab$KO_ID`. If we were to run `gene_annotations_tab$KO_ID` by itself, it would print all 84,784 values within that column as a vector. Wrapping that with the function `is.na()` like we're doing here will return whether or not each item is an `NA` or not. This would by itself give us a table of all rows where the column KO_ID had an `NA`. But to get all of those that were annotated by KEGG, we flip this with the `!` character in front of the function.  

As mentioned above, there is more going on under the hood here that you can dig into at will at the [going deeper with indexing page](/R/more_indexing) that will help you get a better handle on how R does this. If it seems a little abstract and hard to follow at first, no worries at all! It gets better quickly and only becomes more glorious as you go 🙂

## write.table()

Now, let's write out our new table of only those genes that were annotated by KEGG to a new tab-delimited file called "all_KEGG_annotated_genes.txt". We can do this with the `write.table()` function. If we glance at the help menu for this with `?write.table`, we see that the default delimited seems to be a blank space, so we need to be sure to specify that we instead want it to be tab-delimited by providing the argument `sep="\t"`. But we also don't want to keep our row names anymore ( `row.name=F` ), and by default it will write out quotation marks around character strings (like our annotation columns) which we also don't want ( `quote=F` ). So let's add in these additional arguments to make the file to our liking (check the help menu again for details on the ones we've added here).

```R
write.table(all_KEGG_gene_annotations_tab, "all_KEGG_annotated_genes.txt", sep="\t", row.names=F, quote=FALSE)
```

And it's good practice to peek at the output in the terminal when you are configuring the options to write something out for the first time to make sure it's doing what you think it's doing:

<center><img src="{{ site.url }}/images/KEGG_only_tab_less.png"></center>
<br>

---
---
<br>
<h1>Congrats on getting through the basics of R!</h1>
There is of course so much more to it, but this is a great foundation. And once you have these fundamentals down, everything after is cake 🙂
