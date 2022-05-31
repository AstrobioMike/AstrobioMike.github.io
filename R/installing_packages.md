---
layout: main
title: Installing R packages
categories: [R, tutorial]
permalink: /R/installing_packages
---

{% include _R_installing_packages_toc.html %}

{% include _side_tab_R.html %}

Part of what makes R so valuable is that there is an enormous community of people developing software packages for it. People share bundles of code that perform specific tasks through what are known as "packages". Packages are typically maintained at the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/web/packages/){:target="_blank"} and/or at [Bioconductor](https://www.bioconductor.org/){:target="_blank"}. To use a package, we fist need to install, if we don't have it yet, and then load it. Here we'll cover the 3 main ways for installing packages.

> Keep in mind that the best way to install something is generally going to be the way the package documentation recommends. So that's the first place we should start whenever possible. For instance, the dada2 package has it's recommended method of installation right at the top of [its documentation](https://benjjneb.github.io/dada2/){:target="_blank"}. 

<br>

---
---
<br>

## install.packages()
Often, we will be able to install packages with the `install.packages()` function. For example, if we want to install the glorious [tidyverse package](https://tidyverse.tidyverse.org/){:target="_blank"}, we would run this:

```R
install.packages("tidyverse")
```

This should print a bunch of info to the screen while it's installing, then afterward we would load the library with:

```R
library("tidyverse")
```

But occasionally when using `install.packages()` we will get a message like the following:

```R
install.packages("dada2")
```

```
Warning in install.packages :
  package ‘dada2’ is not available for this version of R
```


But do not despair! This is usually just a consequence of the package not being maintained to be installable with this particular method. Doing a google search for [install dada2](https://www.google.com/search?q=install+dada2&rlz=1C5GCEM_enUS1005US1006&oq=install+dada2&aqs=chrome.0.69i59j0i512l2j0i22i30l4j69i60.1680j0j4&sourceid=chrome&ie=UTF-8){:target="_blank"} will allow us to find the [installation instructions page](https://benjjneb.github.io/dada2/dada-installation.html){:target="_blank"}, which tells us to use `BiocManager::install()` to install this particular package.  
<br>

---
<br>

## BiocManager::install()
BiocManager handles all of the packages hosted on [Bioconductor](https://bioconductor.org/){:target="_blank"}. We do first need to install this before we can use it, using the above method like so: 

```R
install.packages("BiocManager")
```

Then, as we saw above, the [dada2 installation instructions](https://benjjneb.github.io/dada2/dada-installation.html){:target="_blank"} tells us to install it like so: 


```R
BiocManager::install("dada2")
```

And again, after installing, we'd need to load the library in order to access the functions it contains, and we do that with `library("dada2")`.  
<br>

---
<br>

## install_github()
A semi-common way I also find myself using is to install the development version of something that is available on GitHub. We might want to do this, for example, if a bug was fixed in the code but it isn't updated in the packaged release that is available with the above methods. 

For one to do as an example, a little googlation for ["tidyr github"](https://www.google.com/search?q=tidyr+github&rlz=1C5GCEM_enUS1005US1006&oq=tidyr+github&aqs=chrome..69i57j0i433i512j0i512l8.1952j0j7&sourceid=chrome&ie=UTF-8) returns as the top hit the [tidyr package github page](https://github.com/tidyverse/tidyr){:target="_blank"}, and if we scroll down a little bit there are installation instructions that include how to install from github:

```R
install.packages("devtools")
devtools::install_github("tidyverse/tidyr") 
```

And after that finishes up, again, in order to be able to use it, we just need to load it: 

```R
library("tidyr")
```

<br>

---
<br>

> These are the 3 most common ways to install things in R. Finding the package documentation, and doing what they recommend is almost always the best way to go. If that isn't possible for some reason, I'd try the above methods probably in the order they are presented above.
