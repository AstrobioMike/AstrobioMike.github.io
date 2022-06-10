---
layout: main
title: R
categories: [R]
permalink: /R/
---

{% include _side_tab_R.html %}

# R is our friend
[R](https://www.r-project.org/) is generally a big, glorious calculator. More specifically, it's a programming language and working environment for statistical analyses and figure generation, and it's pretty sweet. What really makes R powerful is that it is completely open source and it has a tremendous number of incredible people that contribute "packages" to it for all to use. Packages are bundles of code that perform specific tasks, and they are kinda like Apps in the sense that for most of the things you'll want to do, "there's a package for that". Often you will end up searching for a package that you know exists already either from hearing about it from someone, or seeing it used in a paper, but you can also search the [Comprehensive R Archive Network (CRAN) for packages](https://cran.r-project.org/web/packages/) directly. R serves as the foundation upon which you can utilize this large swath of tools that people all over the world have developed and contributed, and it also is invaluable for parsing tables and creating figures. It's definitely worth the initial time investment it may take to get comfortable with it.  

There is extensive documentation on R at the CRAN [Introduction to R site](https://cran.r-project.org/doc/manuals/r-release/R-intro.html), and it would probably be worthwhile going through it at some point. Consistent with the general approach of the site here, I try to distill things down to just the baseline skills to start in the [R basics]({{ site.url }}/R/basics) page.  

<br>

---
<br>
# R pages
* [Introduction to R](/R/basics)  
* [Installing R packages](/R/installing_packages)  
* [Managing R and RStudio with conda](/R/managing-r-and-rstudio-with-conda)  
* [Other great resources](/R/other_resources)  

<hr style="height:10px; visibility:hidden;" />

---
<br>
<h3>Some terminology</h3>

**Variable**  
A variable is something you define that stores data within it.

**Object**  
An object is a data structure that has specific, known attributes to it that allow you to manipulate it more easily. 

**Vector**  
A one-dimensional structure holding more than one item.

**Data frame**  
A two-dimensional table that can hold different types of data, like categorical, qualitative values and continuous, numerical values in the same table.

**Matrix**  
A two-dimensional table that holds only the same type of data.

**Indexing**  
A way of subsetting R objects down to whatever you want. There is plenty to get yourself comfortable with indexing in R [here](/R/basics#the-wonderful-world-of-indexing) and [here](/R/more_indexing). 

**Function**  
A function is a collection of code that performs a desired task. For example, the function `sum()` adds up all of the components of something and gives the sum back to us. 

**Package (Library)**  
A package (also called library) is a collection of functions. We will often want to use functions already written and shared by wonderful people that are not a part of what comes standard with R, and [downloading and installing packages](/R/installing_packages) makes this easy for us to do.
