---
layout: main
title: An introduction to Conda
categories: [unix, tutorial]
tags: [unix,bash,conda,bioinformatics,tools,program,installation,tutorial]
permalink: /unix/conda-intro
---

{% include _unix_conda_toc.html %}

{% include _side_tab_unix.html %}

---
<br>
**This is a page about [Conda](https://conda.io/docs/){:target="_blank"} in the style of this site, but there is also excellent documentation available from the Conda developers and community [here](https://conda.io/projects/conda/en/latest/user-guide/index.html){:target="_blank"}. Thank you, Conda team 🙂**

---
<br>
# What is Conda and why do we love it?
[Conda](https://conda.io/docs/){:target="_blank"} is a package and environment manager that is **by far the easiest way to handle installing most of the tools we want to use** in bioinformatics. Being "conda-installable" requires that someone (could be the developer, could be others) has gone through the trouble of making it that way, so not *everything* is available, but almost everything we're likely to want to use is. Going hand-in-hand with making things easier to install is conda's other value, that it **handles different environments very nicely for us**. Sometimes `Program A` will depend on a specific version of `Program B`. But then, `Program C` will depend on a different version of `Program B`, and this causes problems. **Conda lets us easily create and manage separate environments to avoid these types of version conflicts**, and automatically checks for us when we try to install something new (so we find out now, before we break something somewhere under the hood and have no idea what happened). The benefits go further, like *really* helping with reproducibility too, but let's get into it! 

>**NOTE:** This page assumes already having some familiarity with working at the command line. If that's not the case yet, then consider running through the [Unix crash course](/Unix/unix-intro){:target="_blank"} first 🙂

<hr style="height:10px; visibility:hidden;" />

---
---
<br>
# Binder available
We'll most likely want to be doing this on our own system eventually, but if we just want a temporary system to run through this tutorial, we can open a [Binder](https://mybinder.org/){:target="_blank"} by clicking this badge – [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AstrobioMike/binder-conda-intro/master?urlpath=lab){:target="_blank"}. After the screen loads, we can click the "Terminal" icon under "Other" to launch our command-line environment:

<center><img src="../images/binder-conda-app-launch.png" width="90%"></center>
<br>

<hr style="height:10px; visibility:hidden;" />

---
<br>
# Getting and installing Conda
Conda comes in [two broad forms](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda){:target="_blank"}: Anaconda and Miniconda. Anaconda is large with lots of stuff, Miniconda is more lightweight and we install just what we want. For this reason, we're gonna move forward with Miniconda. The download page for Miniconda is [here](https://conda.io/en/latest/miniconda.html){:target="_blank"}, and we should pick the one appropriate for our operating system. Those labeled Miniconda**2** are using python 2 as the base, while those labeled Miniconda**3** are using python 3. Both can still build environments in the other's python version, but more than likely we should just take python 3. 

Most likely, we will want Miniconda3, and for a 64-bit system. Choosing for our appropriate system here (working in the binder linked above), we are going to want a Linux version, so this is one way to download it to our system (from copying the link under "Miniconda Linux 64-bit":

```bash
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

> **NOTE:** Don't forget to get the link that works for your operating system, and modify it above as needed.

Now we can install that by running it as a `bash` command:

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

**We need to interact with the following during installation:**

1. We will be asked to review the licence agreement. **After pressing return/enter, we can view it and hit the `space` bar to advance. At the bottom, we need to type in `yes` to accept.** 
2. It will then ask if the location is ok, **hitting enter will place it in our current home directory**, which is usually great. 
3. At the end it will ask if we want to initialize Miniconda3. There might be reasons to not want to do this, but until we hit one of those reasons, it probably makes things easier to **say `yes` to that**.  

Great! Now we have conda installed. We just need to reload our current shell session, which we can do like so:

```bash
source ~/.bashrc
```

And now we can see a `(base)` at the start of our prompt, which tells us we are in the "base" conda environment. 

<hr style="height:10px; visibility:hidden;" />

---
<br>
# Creating and navigating environments
Here, "environments" represent the the infrastructure our computer is currently operating in. This includes things like if certain variables exist and how they are set (like our [PATH](https://astrobiomike.github.io/unix/modifying_your_path){:target="_blank"}). 

## Base environment  
The "base" conda environment is, like it sounds, kind of our home base inside conda. We wouldn't want to install lots of complicated programs here, as the more things added, the more likely something is going to end up having a conflict. But the base environment is somewhere we might want to install smaller programs that we tend to use a lot (example below). 

## Making a new environment
The simplest way we can create a new conda environment is like so:

```bash
conda create -n new-env
```

Where the base command is `conda create`, then we are specifying the name of our new environment with `-n` (here "new-env"). It will check some things out and tell us where it is going to put it, **when we hit `y` and enter**, it will be created. 

## Entering an environment

To enter that environment, we need to execute:

```bash
conda activate new-env
```

And now we can see our prompt has changed to have `(new-env)` at the front, telling us we are in that environment. 

If we had forgotten the name, or wanted to see all of our environments, we can do so with:

```bash
conda env list
```

Which will print out all of the available conda environments, and have an asterisk next to the one we are currently in.

## Exiting an environment
We can exit whatever conda environment we are currently in by running:

```bash
conda deactivate
```

## Making an environment with a specific python version
By default, the `conda create` command will use the python version that the base conda environment is running. But we can specify a different one in the command if we'd like:

```bash
conda create -n python-v2.7 python=2.7
```

>**Breakdown**
>* `conda create` – this is our base command
>* `-n python-v2.7` – we are naming the environment "python-v2.7"
>* `python=2.7` – here we are specifying the python version to use within the environment

**After entering `y` to confirm and complete that**, we can see this worked by checking our python version inside and out of these environments.

```bash
conda activate base
python --version

conda activate new-env
python --version # same as base

conda activate python-v2.7
python --version # different than base
```

If we were trying to use a program that was only available in python 2, we could create a conda environment that was prepared to work with it just like this without messing up all the things we use that depend on python 3! (I'll note it's kind of hard to appreciate this sort of thing unless you've fought with the consequences before 🙂)

## Removing an environment
And here is how we can remove an environment, by providing its name to the `-n` flag:

```bash
conda deactivate # we can't be inside the environment we want to remove

conda env remove -n python-v2.7
```

<hr style="height:10px; visibility:hidden;" />

---
<br>
# Finding and installing packages
To install the packages (tools/programs) we want, we use the `conda install` command. But we need to tell conda where to look for it (**conda stores and looks for packages in different "channels"**), and we need to know what it is called in the conda system (usually this is just the program name we are familiar with). Let's imagine we want to install the protein-coding gene-prediction program [prodigal](https://github.com/hyattpd/Prodigal){:target="_blank"}. 

Let's change into our "new-env" and check if `prodigal` is there already:

```bash
conda activate new-env

prodigal -v
```

And we get a "command not found" error if this isn't installed on our system and accessible in our current environment. So let's look at a common way to find it and install it through conda!

## Searching for packages
The first thing I usually do, whether I know if a program is conda-installable or not, is just search in a web-browser for "conda install" plus whatever program I am looking for. For example, searching for "[conda install prodigal](https://www.google.com/search?q=conda+install+prodigal&oq=conda+install+prodigal&aqs=chrome..69i57.2079j0j7&sourceid=chrome&ie=UTF-8){:target="_blank"}" on google brings the top hit back as the [prodigal package page](https://anaconda.org/bioconda/prodigal){:target="_blank"} on [anaconda.org](https://anaconda.org/){:target="_blank"} (which is where these are all hosted). [Anaconda.org](https://anaconda.org/){:target="_blank"} is also a great place to look if we are having trouble finding a program.

## Installing packages
Looking at the [prodigal package page](https://anaconda.org/bioconda/prodigal){:target="_blank"}, we can see there are instructions for installing through conda. These specify where they are specifying the appropriate channel to look in with the `-c` flag. So we can just run this line to install `prodigal`:

```bash
conda install -c bioconda prodigal
```

This prints out such as what and what versions of things are going to be installed, and where they are going to be installed (which will end with the name of the environment we created). **After inputting `y` and pressing return/enter**, the program is installed in this environment, and we can check for `prodigal` again, only this time successfully:

```bash
prodigal -v
```

And now since the program is installed in this environment, we get the version information printed out as we should. If we move out of this environment, it will not be found again (unless it is installed there too):

```bash
conda deactivate

prodigal -v
```

> **NOTE:** Notice that the [prodigal github page](https://github.com/hyattpd/Prodigal){:target="_blank"} doesn't list the conda installation there! This is not all that uncommon, as developers aren't always the folks that make their programs conda-installable. Not seeing conda installation instructions on a particular program's installation page, does not mean it is not available through conda. Be sure to look for it like shown above!

## Uninstalling a package
If we wanted to remove prodigal, we would want to first be inside the environment that has it, so let's switch back in to our new environment:

```bash
conda activate new-env
```

And then we would run:

```bash
conda uninstall prodigal
```

**Once we confirm that with `y` and hit return**, `prodigal` is once again gone:

```bash
prodigal -v
```

## Installing a specific version of a package
If we wanted to specify a different version of `prodigal` to install, we first need to know it exists in conda. We can see all the version types available like so:

```bash
conda search -c bioconda prodigal
```

This lists 2 different versions available (2.6.2 and 2.6.3, with 3 different "builds" for each. A build here is when the program isn't changed, but how it was packaged for conda has changed. If we wanted to install v2.6.2 instead of v2.6.3 like we got above, we could specify it like so:

```bash
conda install -c bioconda prodigal=2.6.2
```

Just like above when we specified the python version we wanted, we can place an equals sign and then the version we want of any program we are installing. 

**Confirming that by entering `y` and hitting return** finishes the installation, and we can see we have this older version now:

```bash
prodigal -v
```

<hr style="height:10px; visibility:hidden;" />

---
<br>
## A note on channels
When we do something in conda, it automatically searches our stored channels (and it does so in [a specific order](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html){:target="_blank"}). We had to specify `-c bioconda` in the command above when we installed `prodigal` because that channel wasn't yet in our stored channels. We can either set these ahead of time, or we can specify them when we install something. Most of the things we biologists will want to use are found in the [bioconda channel](https://bioconda.github.io/index.html){:target="_blank"}. The [prodigal package page instructions on anaconda.org](https://anaconda.org/bioconda/prodigal){:target="_blank"} instructions are actually not 100% ideal. Looking at the [bioconda documentation](https://bioconda.github.io/user/install.html#set-up-channels){:target="_blank"} we can see that we should specify our channel priority such that `conda-forge` is searched before `bioconda` which is searched before `defaults`. 

Even though it worked ok with `prodigal` when we only gave it the `-c bioconda` channel, other programs might run into an issue. So it's best to follow the [bioconda documentation](https://bioconda.github.io/user/install.html#set-up-channels){:target="_blank"} and specify all the channels in the proper order. Here are two ways we can do that. 

**We can specify these channels in the installation command like so:**

```bash
conda activate new-env # making sure we are in our new environment

conda install -c conda-forge -c bioconda -c defaults prodigal
```

(Though in this case it is going to ask us if we want to update the version we have to a newer version, **we can just hit `n` and return**.)

**Or we can set the channels ahead of time** like the [bioconda documentation](https://bioconda.github.io/user/install.html#set-up-channels){:target="_blank"} demonstrates (when done this way, the last one we add has the highest priority):

```bash
conda config --add channels defaults # no need to worry about the warning
conda config --add channels bioconda
conda config --add channels conda-forge
```

And we can see our channel list and priority order with the following:

```bash
conda config --get channels
```

Now we'd be able to install things from bioconda without needing to specify the channels if we wanted, e.g.:

```bash
conda install prodigal
```

(Again this will ask us if we want to update. If we had the same version as the latest, it would just say it is installed alrady. **We can just hit `n` and return again.**)

I tend to list out the channels when installing something just out of habit now, like in the first example above (as it's almost always bioconda in my case). This has the added benefit that I very frequently accidentally type out "conda-forage" 🙂

<hr style="height:10px; visibility:hidden;" />

---
<br>

# An example of something we might want to install in our base environment
As mentioned above, the base environment is somewhere we might want to install smaller programs that we tend to use a lot. For example, I have a collection of small scripts I use very frequently stored in [a package here](https://github.com/AstrobioMike/bioinf_tools#bioinformatics-tools-bit){:target="_blank"} called `bit` – for **b**io**i**nformatics **t**ools ¯\\\_(ツ)_/¯

Let's switch into our base environment and install it (the `-y` flag makes is so we don't need to confirm anything):

```bash
conda activate base

conda install -y -c conda-forge -c bioconda -c defaults -c astrobiomike bit
```

> **Breakdown**
> * `conda install` – our base command
> * `-y` – says not to ask us for any confirmation
> * `-c ...` – specifies the channels
> * `bit` – this positional argument is specifying the program name as it is within conda

This might take about a minute, once it is finished we can try one of them:

```bash
bit-summarize-assembly random-assembly.fa
```

Because I typically use this and other scripts in here all over the place (we can see them all by typing `bit-` then hitting tab twice), it makes sense for me to have them in my base environment. But we wouldn't want to put too many large things into one environment unless it is needed, as it increases the likelihood of version conflicts.

<hr style="height:10px; visibility:hidden;" />

---
<br>
# Creating an environment and installing packages in one command
We can also create an environment and install tools at the same time. Here is creating an environment for [GToTree](https://github.com/AstrobioMike/GToTree/wiki){:target="_blank"} and installing it at the same time:

```bash
conda create -y -n gtotree -c conda-forge -c bioconda -c astrobiomike gtotree python=3.7
```

> **Breakdown**
> * `conda create` – our base command
> * `-y` – says not to ask us for any confirmation
> * `-n gtotree` – creates the environment named "gtotree"
> * `-c ...` – specifies the channels
> * `gtotree` – this positional argument is specifying the program name as it is within conda
> * `python=3.7` – specifies the python version that GToTree expects to operate in

When that is done, we can change into the environment, and see that it is there:

```bash
conda activate gtotree

GToTree -v

gtt-hmms
```

Again, it might be hard to appreciate if we haven't fought with installations, versions, and environments before, but [Conda](https://conda.io/docs/){:target="_blank"} is stellar 🙂

<hr style="height:10px; visibility:hidden;" />

---
---
<br>
# Other Conda Resources
* [Main Conda documentation](https://docs.conda.io/en/latest/){:target="_blank"}
* [Main Conda user guide](https://conda.io/projects/conda/en/latest/user-guide/index.html){:target="_blank"}
* [Getting started tutorial from Conda folks](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html){:target="_blank"}
* [Conda cheat sheet](https://docs.conda.io/projects/conda/en/latest/user-guide/cheatsheet.html){:target="_blank"}
* [Anaconda.org](https://anaconda.org){:target="_blank"} (A good place to search for packages if a regular web-search is failing us)
