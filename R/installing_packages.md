---
layout: main
title: Installing R packages
categories: [R, tutorial]
permalink: /R/installing_packages
---

{% include _R_installing_packages_toc.html %}

{% include _side_tab_R.html %}

Part of what makes R so valuable is that there is an enormous community of people developing software packages for it. People share bundles of code that perform specific tasks through what are known as "packages". Packages are typically maintained at the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/web/packages/){:target="_blank"} and/or at [Bioconductor](https://www.bioconductor.org/){:target="_blank"}. To use a package, we simply need to install it and then load it. Here we'll cover the 3 main ways for installing packages.
<br>
<br>

---
---
<br>

## install.packages()
Often, we will be able to install packages with the `install.packages()` function. For example, if we want to install the [ggplot2 package](https://cran.r-project.org/web/packages/ggplot2/index.html){:target="_blank"}, we would simply enter:

```R
install.packages("ggplot2")
```

And we get some info printed to our screen such as this:

<center><img src="../images/ggplot2_install.png"></center>
<br>
And all is well with the world and we are ready to load the package with `library("ggplot2")`, at which point we're ready to rock.

But occasionally when using `install.packages()` we will get a message like the following:

<center><img src="../images/phyloseq_install_packages.png"></center>
<br>
But do not despair! This is usually just a consequence of the package not being maintained to be installable this particular way, and we can almost always get around it by installing from [Bioconductor](https://bioconductor.org/){target="_blank"} directly. When we do run into this, we should head right on over to google and search with terms for the package name and bioconductor, and we'll most likely find a way to install with `BiocManager::install()`.  
<br>

---
<br>

## BiocManager::install()
BiocManager handles all of the packages hosted on [Bioconductor](https://bioconductor.org/){target="_blank"}. And is my usual next attempt if `install.packages()` does not work. If we don't have it yet, this can be installed with:

```R
install.packages("BiocManager")
```

Then, for example, searching for ["phyloseq bioconductor"](https://www.google.com/search?q=bioconductor+phyloseq+R&oq=bioconductor+phyloseq+R&aqs=chrome..69i57j69i60.5566j0j7&sourceid=chrome&ie=UTF-8) returns the [bioconductor homepage of the *phyloseq* package](http://bioconductor.org/packages/release/bioc/html/phyloseq.html){:target="_blank"} as the top hit. And when we head over there, there are directions to install directly from Bioconductor:

```R
BiocManager::install("phyloseq")
```

And again, after installing, we'd need to load the library in order to access the functions it contains: `library("phyloseq")`.  
<br>

---
<br>

## install_github()
Occasionally Bioconductor may also not workout. At that point I usually turn to searching for the package on github as it may be hosted there, and if so it usually contains instructions for how to install it via devtools and the `install_github()` function. 

For example, a little googlation for ["tidyr github"](https://www.google.com/search?ei=6BwBWqfgLca6jwOwsIXACg&q=tidyr+github&oq=tidyr+github&gs_l=psy-ab.3..0.1251.3680.0.3814.18.15.3.0.0.0.141.1174.11j3.14.0....0...1.1.64.psy-ab..1.17.1182...0i67k1j0i131k1j0i10k1j0i22i10i30k1.0.Xit6NDyEZS0) returns as the top hit the [tidyr package github page](https://github.com/tidyverse/tidyr){:target="_blank"}, and if we scroll down a little bit there are installation instructions that include how to install.

```R
install.packages("devtools")
devtools::install_github("tidyverse/tidyr") 

  # and after that finishes up, we just need to load the library for use
library("tidyr")
```
