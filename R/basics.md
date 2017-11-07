---
layout: main
title: R basics
categories: [R, tutorial]
permalink: /R/basics
---

{% include _R_basics_toc.html %}

{% include _side_tab_R.html %}

This module is designed for those that are either completely new to R, or have some experience but maybe don't feel as solid about some of the fundamentals as they'd like. It will run through the very basics such as setting up your working environment, assigning variables, reading in and writing out data, and installing packages. Some relevant terminology is presented [here]({{ site.url }}/R/index) if you find yourself seeing some words that are unfamilar to you. This page is meant to be solely a quick jump start to get you into and using the R environment. But there is, of course, an incredible amount of functionality we won't be touching on here. Be sure to peruse the extensive [R intro documentation available here](https://cran.r-project.org/doc/manuals/r-release/R-intro.html) when you're starting to get your bearings and want to go further.
<br>
<br>

---
---
<br>

# Installation
It is possible your computer already has R, if you are unsure, you can check by opening a terminal window and typing `R`. If this launches R rather than giving an error message, you should be good to go (enter `q()` to exit the R environment). If you do not have R, you can download it from here for Mac: [https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/). And if you have a relatively newer Mac, you may also need to install XQuartz which you can get from here: [https://www.xquartz.org/](https://www.xquartz.org/). 

Lastly, I highly, *highly*, **highly** recommend installing RStudio if you don't already have it and use it. RStudio is an interface for R that not only makes everything you will do in R easier and more organized, but it's also invaluable for reproducibility of your analyses as it makes it second-nature to generate and save R scripts of everything you're doing, while you're doing it – which is very helpful when you want to look back and see what worked of the 20 things you just tried 🙂 . You can download an RStudio installer from [here](https://www.rstudio.com/products/rstudio/download/#download).

<br>

---
<br>
# RStudio layout
RStudio has 4 main panes: 1) console; 2) source; 3) Environment, History, Connections; and 4) Files, Plots, Packages, etc. You can run commands straight in the console, which would be just like using R in a command line environment, but I personally almost never directly enter something there. The source pane acts as a text document from which you can write out all of your commands and call them when you'd like. This is one of the reasons R Studio is great, you're constantly building up your notes file as you work without any added effort needed. If you want to run a command written in the source file, while on the line of the command press `Cmd + Enter`. The Environment/History pane displays all of the variables and data structures you have currently stored. And the last pane, Files/Plots/etc., allows you to peruse your computer in the typical Finder fashion, displays any plots you generate, and serves as your help window. (To get help on a function in R, just place a `?` in front of the function name. For example, try entering `?getwd`.)  
<br>

---
<br>
# Some practice data
As with the [bash basics]({{ site.url }}/bash/basics) module, if you'd like to follow along with this tutorial, copy and paste the following commands into your terminal to get set up with a small, temporary working directory that has the files I'll be working with. 

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
Just like when working at the command line in *bash*, we need to be aware of "where" we are in our computer when working in R. The `getwd()` and `setwd()` help us do this in R. Commands in R typically take this structure, with the command being followed by parentheses. Inside those parentheses is where the arguments to the command would go if there are any. In the case of `getwd()` no arguments are needed as it is just going to tell you where you currently are. But in the case of `setwd`, we then need to tell it where we want to go. 

```R
  # check where you are in R by default when it boots up
getwd()
  # and let's set our working directory to where we want it to be
setwd("~/R_basics_temp/")
```

Now that we are in the correct location that contains the files for our tutorial here, let's check that that are actually there. In *bash* this would be done with the `ls` command, for list, here in R we do it with the `list.files() function.

```R
list.files()
```

And this should reveal our practice files. Just like when working at the command line, you can "point" to other places in your computer by providing the relative or absolute path to whichever file/directory you're looking for in the same manner as described [here]({{ site.url }}/bash/basics#moving-around).  
<br>

---
<br>
# Basic calculations
As mentioned, R is just a big calculator at its core. So we can do baseline arithmetic operations like the following examples. (Don't forget if you're entering these in the "source" pane of RStudio, you need to press `Cmd + Enter` to execute the command.)

```R
4 + 4
4 / 2
4 * 4
2 ^ 4 # the ^ symbol is for raising to a power
```

<center><img src="{{ site.url }}/images/R_basic.png"></center> 
<br>
Just to note, the "[1]" you're seeing to the left of the output is an "index" number. As we go further you'll see why these are useful, but for now just know it's a counter for how many iterms are returned by a command, and that for each row of output it lists the "index" number of the first item of that row.  
<br>

---
<br>
# Assigning variables
R acts on objects that are stored in variables. These objects can hold different "classes" of data that determine what they can do and what can be done to them, and we'll get into a little bit as we go. We assign values to variables, and name them, as shown in the following code block. Here we're naming the variable "x", and giving it a value of 4.

```R
x <- 4
```

When you execute the command, it should print in the console pane, and now in your environment pane you should see there is a variable there called "x" (named arbitrarily). For a note on syntax, it doesn't matter if you put spaces or not between the variable name, the assignment operator ( `<-` ), and the value you're assigning to it. But generally it makes the code easier to read if you do, and it is good practice to make things as easy to read as possible for others who may be looking at it in the future. 

Now that we have the value 4 stored in the variable "x", we can do stuff with it. Here's some examples doing the same calculations we performed above, but now with our variable.

```R
x
x + x
x / 2
x * x
2 ^ x
```
<center><img src="{{ site.url }}/images/R_var_arith.png"></center> 
<br>
We can also check what class of data is contained within this variable with the `class()` function:

```R
class(x)
```

And we find out it is of class "numeric". Let's try storing a different type, like a word:

```R
word <- "europa" 
class(word)
```

Here, the class is "character". It's important to understand why we used quotes here, so let's take a second and look at it. Above, "europa" is a known as a string of characters, and we need to tell that to R by surrounding it in quotes like we did. If we didn't provide the quotes, R would be looking for an object named europa, and it would give us an error saying it wasn't found, give it a shot without the quotes to see. 

We also can store multiple items into a variable known as a vector. We'll do this here with the `c()` function. This concatenates items together. Here we'll make a vector of numbers:

```R
y <- c(5, 6, 7)
y
```

<center><img src="{{ site.url }}/images/R_vector.png"></center> 
<br>
Note that this is still of class "numeric" by checking with `class(y)`. It's good practice to get used to actively being aware of what class of objects you are working with. 

Variables also hold tables. Here we're going to make another vector with 3 numbers, and then combine them to make what's known as a dataframe. We'll do this with the `data.frame()` function, creating a variable called "our_table".

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
Dataframes are two-dimensional tables of rows and columns. Here we can see that the default behavior of the `data.frame()` function took our two vectors and put them in a table where each represents one column. Another similar, but distinct, table structure in R is a "matrix". You sometimes find you need to convert a dataframe to a matrix or vice versa depending on what you are trying to do with it. Keep this in mind as one of the things to look at first when you run into an error.  
<br>

---
<br>
# The gloriousness of indexing
One of the most powerful things about R is how easy it makes it to parse vectors and tables down to whatever you are interested in via what's known as indexing. Here we'll look at a couple of the ways we can specify what we would like to subset, and we'll see these in practice on a larger scale below. 

Looking back at our vector stored in variable y, it contains 3 values: 5, 6, and 7. One way we can subset specific values from this involves using their position in the vector. If we specify the vector name, and then put in brackets the position we are interested in, R will return the values for just that (or those) position(s) we specify.

```R
y # the whole vector
y[1] # the first item
y[2] # second item
y[3] # third item
y[c(1,3)] # specifying items 1 and 3 using the c() function
```

<center><img src="{{ site.url }}/images/R_vec_index.png"></center> 
<br>
Another way to subset by indexing makes use of `TRUE` and `FALSE` values. Meaning, subsetting with brackets will also return only those values that evaluate to `TRUE`. What makes each one `TRUE` or `FALSE` is dependent upon the conditions we set. For instance, let's say we only want the values from our "y" vector that are greater than or equal to 6. The following code returns 6 and 7 to us as it should:

```R
y # entire vector is returned
y[y >= 6] # only the values that meet our criterion of being greater than or equal to 6 are returned
```

<center><img src="{{ site.url }}/images/R_vec_index_TF.png"></center> 
<br>
The way I read the expression `y[y >= 6]` in my head is: "Give me all the values of vector 'y', where 'y' is greater than or equal to 6." It's imporant to take a second and think about what's going on under the hood here. Let's see what happens when we just run the conditional part that we put in the brackets.

```R
y >= 6
```

<center><img src="{{ site.url }}/images/R_vec_index_TF2.png"></center> 
<br>
This returns a vector of class "logical", which you can check by running `class(y >= 6)`. So what's actually happening when we run `y[y >= 6]` is: 1) is R is looking at each value in the vector – 5, 6, and 7 – and returning a vector that has `TRUE` wherever the condition was met, and `FALSE` wherever it wasn't met; and 2) This T/F vector is what's actually doing the subsetting for us in the command. To R, this is the same as if we were to subset vector 'y' by providing the T/F vector directly, as in this example:

```R
y[c(FALSE, TRUE, TRUE)]
```

<center><img src="{{ site.url }}/images/R_vec_index_TF3.png"></center> 
<br>
One last thing we're going to see here is the `!` character, which inverts the interpretation of `TRUE` and `FALSE`. As we've seen, `y[y >= 6]` will return all values within 'y' that are greater than or equal to six, but if we add in the `!` point, it will return the opposite for us. We can see this just by putting this in front of the portion of the expression that generates the T/F vector:

```R
y >= 6 # returns F, T, T
!y >= 6 # returns T, F, F

y[y >= 6] # returns 6 and 7
y[!y >= 6] # returns 5
```

<center><img src="{{ site.url }}/images/R_vec_index_TF4_not.png"></center> 
<br>
We'll see a somewhat more complicated example of using the `!` character to reverse the "logical" vector we are using to index a bit later.  

We can also subset tables in a similar way. Vectors are one-dimensional objects, so when we index them we only need to provide one value. But tables are 2-dimensional objects, so we need to provide another value (one for row, and one for column). The syntax for how this is done is by providing the variable name first followed by the subsetting brackets (as with vectors), but then the first value you enter in the brackets specifies which **rows** you'd like, and the the second value (separated by a comma) specifies which **columns**. Let's look at our table again and then see some examples of this:

```R
  # whole table
our_table

  # subset the value in the second row and second column only
our_table[2, 2]

    ### if we provide nothing for either the row position or the column position, we will get all values
  
  # subset all rows, but only the second column
our_table[ , 2] # notice when subsetting returns only one column, it returns a vector

  # only row 3, but both columns
our_table[3, ] # notice when subsetting returns one row, but more than one column, it still returns a dataframe, this is important if it is wrapped up in larger code

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
Indexing is *much* more expansive than this, but this is a good start. If any of this seems a little abstract at first, don't worry at all, that's normal (at least it was for me). Like everything, exposure and practice with the process will make you comfortable with it very soon. And as you begin working with large tables and combining subsetting metrics you'll quickly begin to see how useful this fundamental system is. We'll have several larger-scale examples below.  
<br>

---
<br>
# Reading in and writing out data
Now let's read in our practice data table. Rather than give you something straightforward, we're going to intentionally run into a couple of the most frequent snags you'll hit when trying to read in data. A common method of reading in tables is the `read.table()` command. To start, let's try reading our gene_annotations.txt table into R with no arguments other than specifying the file name: 

```R
gene_annotations_tab <- read.table("gene_annotations.txt")
```

<center><img src="{{ site.url }}/images/read_table_err.png"></center> 
<br>
Yay errors! Many of these may seem a little cryptic at first, but you'll be surprised at how many of them start to make sense pretty quickly with just some time. The important part in this one is the end that says "line 1 did not have 22 elements". This is a sign there is something up with how R is trying to split each line into columns.

When there are many things that can vary, (like the settings for how to read in and treat a table), it's important to scan the default parameters of the functions you are using. If we look at the help menu for this command with `?read.table`, we find that the default delimiter ("sep" argument) acts on all white space, which includes tabs AND blank spaces.

<center><img src="{{ site.url }}/images/read_table_help_sep.png"></center> 
<br>
If we switch back to our terminal and look at our gene_annotations.txt file `less -S gene_annotations.txt`, we see that it is tab-delimited, but there are also spaces within the KO and COG annotation columns. 

<center><img src="{{ site.url }}/images/terminal_table_less.png"></center> 
<br>
With the default settings of `read.table()`, R is trying to split columns on blank spaces as well as tabs, this is causing the error we saw as it leads to rows with different numbers of columns. Let's try running the command again, but this time specifying the delimiter should only be tabs (tab characters are specified with an backslash followed by a t: `\t`.

```R
gene_annotations_tab <- read.table("gene_annotations.txt", sep="\t")
  # and let's take a peek at the table with head():
head(gene_annotations_tab) 
```

That worked, but we have a header with column names in the first line, but we didn't tell R that so it assigned names (V1, V2, etc.). Let's read it in again but this time specify that there is a header:

```R
gene_annotations_tab <- read.table("gene_annotations.txt", sep="\t", header=TRUE)
  # now let's take a peek at the table with head():
head(gene_annotations_tab) 
  # and check what our column names are:
colnames(gene_annotations_tab) # looks good
  # and let's check out the size of the table:
dim(gene_annotations_tab)
```

<center><img src="{{ site.url }}/images/read_table_dim.png"></center> 
<br>
So our table is 84,784 rows by 8 columns. Notice that since the table doesn't fit all the way across my Console screen, it wraps twice when we viewed it with the `head` command. Also notice that R by default gave it column names that are numbered. This is usually fine, but in this case since our gene_IDs are so similar, we're going to rename the row names to match the gene_ID for each row (to avoid any possible confusion when glancing at the table). We can view and set the row names of our table like with the `row.names()` function. As we saw above, there are a lot of rows here, so we're going to wrap this function with the `head()` function to keep the output manageable.

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

Now let's generate a new table so we can practice writing out to a file from R. You may have noticed there are some NAs in our gene_annotations_tab table. These are present in the KEGG and COG annotation and ID columns for those genes which weren't annotated. Let's subset our full table to include only those genes that were annotated by KEGG. R's `is.na()` function can help us subset accordingly by generating a "logical" T/F vector just like we did above with our `y[y >= 6]` example above. This will act on each element in a vector and return `TRUE` if it contains an NA value at the respective position, which we can then use to index the full table as R will only pull out the rows where the value was found to be `TRUE`. However, we don't want the NAs, we want all of those that have KEGG annotations. So we need to invert the "logical" vector with the `!` character also like we saw above. This combines a few concepts, so let's run the code and then break down the syntax.

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
Now let's look closely at what's going on in the code we used to generate our subset of only KEGG-annotated genes. We start off with the name of the new table "all_KEGG_gene_annotations_tab", then we have our assignment character `<-` , then the part that is actually specifying what we want from the original table. Same as above we are giving the variable name that holds the table ("gene_annotations_tab"), followed by brackets that enclose which rows we want before the comma, and which columns we want after the comma. We have nothing after the comma so we are taking all columns, but the rows position in these brackets is where the magic is happening. With R it can be sometimes easier to read from the inside and work your way out. So we specific the gene_annotations_tab table, and specifically the column named "KO_ID" by entering that following the `$` . That by itself would print all 84,784 values within that column as a vector. But we have that vector wrapped with the `is.na()` function, which is creating a vector that says `TRUE` wherever the value is equal to NA (which is a special value in R), and `FALSE` wherever the value is not equal to NA. That alone would return us all the rows that have NA in the KO_ID column, but we actually want the opposite of that. So the last part is that we added the `!` character right in front of this expression, which basically inverts the T/F vector and tells the subsetting brackets to take only the rows that do **not** have an NA value for KO_ID, and for those rows take all columns.

Again, if this seems a little abstract and hard to follow at first, no worries at all! It gets better quickly as long as you take the time to try to break it down when you're unsure. 

Now, let's write out our new table of only those genes that were annotated by KEGG to a new tab-delimited file called "all_KEGG_annotated_genes.txt". We can do this with the `write.table()` function. If we glance at the help menu for this with `?write.table`, we see that the default delimited seems to be a blank space, so we need to be sure to specify that. But we also don't want to keep our row names anymore, and be default it will write out quotation marks around character strings (like our annotation columns) which we also don't want. So let's add in some additional arguments to make the file to our liking (check the help menu again for details on the ones we've added here).

```R
write.table(all_KEGG_gene_annotations_tab, "all_KEGG_annotated_genes.txt", sep="\t", quote=FALSE, row.names=F)
```

And it's good practice to peek at the output in the terminal when you are configuring the options to write something out for the first time to make sure it's doing what you think it's doing.

<center><img src="{{ site.url }}/images/KEGG_only_tab_less.png"></center>
<br>

---
<br>
# Installing packages
And the last little bit we'll cover here is how to install and load packages for R. 

## install.packages()
Most often, you will be able to install packages with the `install.packages()` function. For example, if you want to install the [ggplot2 package](https://cran.r-project.org/web/packages/ggplot2/index.html), you would simply enter:

```R
install.packages("ggplot2")
```

And you'll get some info printed to your screen such as this:

<center><img src="{{ site.url }}/images/ggplot2_install.png"></center>
<br>
And all is well with the world and you are ready to load the package with `library("ggplot2")`, at which point you're ready to rock.

But occasionally when using `install.packages()` you will get a message like the following:

<center><img src="{{ site.url }}/images/phyloseq_install_packages.png"></center>
<br>
But do not despair! This is usually just a consequence of the package not having been updated to install this particular way, and you can almost always get around it by installing from bioconductor. When you do run into this, you should head right on over to google and search with terms for the package name and bioconductor, and you'll most likely find a way to install with `biocLite()`. 

## biocLite()
For example, searching for ["phyloseq bioconductor"](https://www.google.com/search?q=bioconductor+phyloseq+R&oq=bioconductor+phyloseq+R&aqs=chrome..69i57j69i60.5566j0j7&sourceid=chrome&ie=UTF-8) returns the [bioconductor homepage of the *phyloseq* package](http://bioconductor.org/packages/release/bioc/html/phyloseq.html) as the top hit. And when you head over there, there are directions to install via bioconductor:

```R
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
```

And sure enough that seems to work just swell:

<center><img src="{{ site.url }}/images/phyloseq_bioconductor_install_run.png"></center>
<br>

And again, after installing, you'd need to load the library in order to access the functions it contains: `library("phyloseq")`.

## install_github()
Occasionally bioconductor may also not workout. At that point I usually turn to searching for the package on github as it may be hosted there, alongside instructions for how to install it via devtools and the `install_github()` function. 

For example, a little googlation for ["tidyr github"](https://www.google.com/search?ei=6BwBWqfgLca6jwOwsIXACg&q=tidyr+github&oq=tidyr+github&gs_l=psy-ab.3..0.1251.3680.0.3814.18.15.3.0.0.0.141.1174.11j3.14.0....0...1.1.64.psy-ab..1.17.1182...0i67k1j0i131k1j0i10k1j0i22i10i30k1.0.Xit6NDyEZS0) returns as the top hit the [tidyr package github page](https://github.com/tidyverse/tidyr), and if you scroll down a little bit there are installation instructions that include how to install.

```R
install.packages("devtools")
devtools::install_github("tidyverse/tidyr") 

  # and after that finishes up, you'd just need to load the library for use
library("tidyr")
```

<br>

---
---
<br>
<h1>Congrats on getting through the basics of R!</h1>
There is of course so much more to it, but this is a great foundation. And once you have these fundamentals down, everything after is cake 🙂
