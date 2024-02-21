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
<hr style="height:10px; visibility:hidden;" />

**This is a page about [Conda](https://conda.io/docs/){:target="_blank"} in the style of this site, but there is also excellent documentation available from the Conda developers and community [here](https://conda.io/projects/conda/en/latest/user-guide/index.html){:target="_blank"}. Thank you, Conda team 🙂**

---
<br>

# What is Conda and why do we love it?
[Conda](https://conda.io/docs/){:target="_blank"} is a package and environment manager that is **by far the easiest way to handle installing most of the tools we want to use** in bioinformatics. Being "conda-installable" requires that someone (could be the developer, could be others) has gone through the trouble of making it that way, so not *everything* is available, but almost everything we're likely to want to use is. Going hand-in-hand with making things easier to install is conda's other value, that it **handles different environments very nicely for us**. Sometimes `Program A` will depend on a specific version of `Program B`. But then, `Program C` will depend on a different version of `Program B`, and this causes problems. **Conda lets us easily create and manage separate environments to avoid these types of version conflicts**, and automatically checks for us when we try to install something new (so we find out now, before we break something somewhere under the hood and have no idea what happened). The benefits go further, like helping with reproducibility too, but let's get into it! 

>**NOTE:** This page assumes already having some familiarity with working at the command line. If that's not the case yet, then consider running through the [Unix crash course](/unix/unix-intro){:target="_blank"} first 🙂

<hr style="height:10px; visibility:hidden;" />

---
---
<br>

# Binder available
We'll most likely want to be doing this on our own system eventually, but if we just want a temporary system to run through this tutorial, we can open a [Binder](https://mybinder.org/){:target="_blank"} by clicking this badge – [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AstrobioMike/binder-conda-intro/v2?urlpath=lab){:target="_blank"}. After the screen loads, we can click the "Terminal" icon under "Other" to launch our command-line environment:

<center><img src="../images/binder-conda-app-launch.png" width="90%"></center>
<br>

<hr style="height:10px; visibility:hidden;" />

---
<br>

# Getting and installing Conda
Conda comes in [two broad forms](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda){:target="_blank"}: Anaconda and Miniconda. Anaconda is large with lots of programs in it already, Miniconda is more lightweight and then we can install just what we want. For this reason, I always use Miniconda, and that's what we're gonna move forward with here. 

The download page for Miniconda is [here](https://conda.io/en/latest/miniconda.html){:target="_blank"}, and we should pick the one appropriate for our operating system, with one minor exception. If working on a new Mac with the M1/2/etc. chip, I currently think (original date of this suggestion is Jun-2022; still holds Nov-2023) it's best to install the regular Intel version, rather than the available Apple M1 version on the [download page](https://conda.io/en/latest/miniconda.html). This is because not a lot of packages in conda have M1 versions built yet, and the M1 computers come with software (called [rosetta](https://support.apple.com/en-us/HT211861){:target="_blank"}) that will ask to be installed the first time it's needed, but then tries to make it so intel-based things work on the M1. And so far this has worked splendidly for me, while installing an M1 based version of conda led immediately to packages not being available. Over time, as more developers put up versions for the M1 chip, this will of course change.

**If working in the binder linked above or on a Linux system**, which is almost definitely the case if you are logged into a server, we are going to want a Linux version, so this is one way to download it to our system (from copying the link under "Miniconda3 Linux 64-bit"):

```bash
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

**If working on a Mac**, we would want the link under "Miniconda3 macOS Intel x86 64-bit bash" (regardless of if we have the new M1 version or not), and the command we would run would be:  

```bash
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```

**If working on Windows**, I believe you want to be in the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install){:target="_blank"} environment if possible (though I think being in something like mobaXterm will work the same). And we'd want to install the "Miniconda3 Linux 64-bit" (though I don't have a windows machine available to test this more thoroughly on 😞 ). Same as the first example above, that would be this command to download it:

```bash
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Now we can install that by running it as a `bash` command:

```bash
bash Miniconda3-latest-*.sh
```

**We need to interact with the following during installation:**

1. We will be asked to review the licence agreement. **After pressing return/enter, we can view it and hit the `space` bar to advance. At the bottom, we need to type in `yes` to accept.** 
2. It will then ask if the location is ok, **hitting enter will place it in our current home directory**, which is usually great. 
3. At the end it will ask if we want to update our shell profile to automatically initialize conda. There might be reasons to not want to do this, but until we hit one of those reasons, for most users it's easist just to **say `yes` to that**.  

Great! Now we have conda installed. We just need to reload our current shell session, which we can usually do like so:

```bash
source ~/.bash_profile || source ~/.bashrc
  # if on mac, and using zsh, you might need to run
# source ~/.zshrc
  # or you can just open a new terminal session
```

And after that we should see a `(base)` at the start of our prompt, which tells us we are in the "base" conda environment. If the above didn't work, you can also just open a new terminal window, and it should be updated there 👍

<hr style="height:10px; visibility:hidden;" />

---
<br>

# Creating and navigating environments
Here, "environments" represent the the infrastructure our computer is currently operating in. This includes things like if certain variables exist and how they are set (like our [PATH](https://astrobiomike.github.io/unix/modifying_your_path){:target="_blank"}). 

<hr style="height:10px; visibility:hidden;" />

## Base environment  
The "base" conda environment is, like it sounds, kind of our home base inside conda. We wouldn't want to install lots of complicated programs here, as the more things added, the more likely something is going to end up having a conflict. But the base environment is somewhere we might want to install smaller programs that we tend to use a lot (example below). 

<hr style="height:10px; visibility:hidden;" />

## Making a new environment
The simplest way we can create a new conda environment is like so:

```bash
conda create -n new-env
```

Where the base command is `conda create`, then we are specifying the name of our new environment with `-n` (here "new-env"). It will check some things out and tell us where it is going to put it, **when we hit `y` and enter**, it will be created. 

<hr style="height:10px; visibility:hidden;" />

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

<hr style="height:10px; visibility:hidden;" />

## Exiting an environment
We can exit whatever conda environment we are currently in by running:

```bash
conda deactivate
```

<hr style="height:10px; visibility:hidden;" />

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

<hr style="height:10px; visibility:hidden;" />

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

<hr style="height:10px; visibility:hidden;" />

## Searching for packages
The first thing I usually do, whether I know if a program is conda-installable or not, is just search in a web-browser for "conda install" plus whatever program I am looking for. For example, searching for "[conda install prodigal](https://www.google.com/search?q=conda+install+prodigal&oq=conda+install+prodigal&aqs=chrome..69i57.2079j0j7&sourceid=chrome&ie=UTF-8){:target="_blank"}" on google brings the top hit back as the [prodigal package page](https://anaconda.org/bioconda/prodigal){:target="_blank"} on [anaconda.org](https://anaconda.org/){:target="_blank"} (which is where these are all hosted). [Anaconda.org](https://anaconda.org/){:target="_blank"} is also a great place to look if we are having trouble finding a program.

<hr style="height:10px; visibility:hidden;" />

## Installing packages
Looking at the [prodigal package page](https://anaconda.org/bioconda/prodigal){:target="_blank"}, we can see there are instructions for installing through conda. These specify where they are specifying the appropriate channel to look in with the `-c` flag. So we can just run this line to install `prodigal`:

```bash
conda install -c bioconda prodigal
```

> **IMPORTANT NOTE**  
> There is more on this below, but because it's important, we're going to note it here too. Despite what the install examples show on the anaconda pages like the one above, for many packages from bioconda to install properly, more channels need to be set in a specific hierarchy (e.g. see this in the [bioconda documentation here](https://bioconda.github.io/#usage). Some installs work fine regardless (as is the case with `prodigal` here), but others will fail, and even worse, some will succeed but have more insidious issues. This is because even though `prodigal` may be in the bioconda channel, some of its dependencies may be found elsewhere or in multiple places, and having those 3 channels set in the proper priority helps to ensure the correct versions of everything are being used. Since I use conda on different systems all the time, I find it much easier and safer just to run the install command specifiying all the channels needed for it in one line like so:
>
> `conda install -c conda-forge -c bioconda -c defaults prodigal`
>
> That way we get what's needed each time, and we don't have to worry about changing the underlying conda channel hierarchy on each system we use it on like the [bioconda documentation](https://bioconda.github.io) demonstrates. 


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

<hr style="height:10px; visibility:hidden;" />

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

<hr style="height:10px; visibility:hidden;" />

## Installing a specific version of a package
If we wanted to specify a different version of `prodigal` to install, we first need to know it exists in conda. We can see all the version types available like so:

```bash
conda search -c bioconda prodigal
```

This lists 2 different versions available (2.6.2 and 2.6.3, with 3 different "builds" for each. A build here is when the program isn't changed, but how it was packaged for conda has changed. If we wanted to install v2.6.2 instead of v2.6.3 like we got above, we could specify it like so:

```bash
conda install -c conda-forge -c bioconda -c defaults prodigal=2.6.2
```

Just like above when we specified the python version we wanted, we can place an equals sign and then the version we want of any program we are installing. 

**Confirming that by entering `y` and hitting return** finishes the installation, and we can see we have this older version now:

```bash
prodigal -v
```

<hr style="height:10px; visibility:hidden;" />

---
<br>

# A note on channels
When we do something in conda, it automatically searches our stored channels (and it does so in [a specific order](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html){:target="_blank"}). We had to specify `-c bioconda` in the command above when we installed `prodigal` because that channel wasn't yet in our stored channels. We can either set these ahead of time, or we can specify them when we install something. Most of the things we biologists will want to use are found in the [bioconda channel](https://bioconda.github.io/index.html){:target="_blank"}. The [prodigal package page instructions on anaconda.org](https://anaconda.org/bioconda/prodigal){:target="_blank"} instructions are actually not 100% ideal. Looking at the [bioconda documentation](https://bioconda.github.io/#usage){:target="_blank"} we can see that we should specify our channel priority such that `conda-forge` is searched before `bioconda` which is searched before `defaults`. 

Even though it worked ok with `prodigal` when we only gave it the `-c bioconda` channel, other programs might run into an issue. So it's best to follow the [bioconda documentation](https://bioconda.github.io/#usage){:target="_blank"} and specify all the channels in the proper order. Here are two ways we can do that. 

**We can specify these channels in the installation command like so:**

```bash
conda activate new-env # making sure we are in our new environment

conda install -c conda-forge -c bioconda -c defaults prodigal
```

(Though in this case it is going to ask us if we want to update the version we have to a newer version, **we can just hit `n` and return**.)

As mentioned above, I prefer doing things this way so I don't have to ever bother with changing the stored channel configuration on any system I use conda on. Plus, this has the added benefit that I very frequently accidentally type out "conda-forage" 🙂

**Or we can set the channels ahead of time** like the [bioconda documentation](https://bioconda.github.io/#usage){:target="_blank"} demonstrates (when done this way, the last one we add has the highest priority):

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

Let's switch back to our base environment before moving forward:

```bash
conda deactivate
```


<hr style="height:10px; visibility:hidden;" />

---
<br>

# BONUS: [mamba](https://github.com/mamba-org/mamba#mamba){:target="_blank"} (No. 5)

[mamba](https://github.com/mamba-org/mamba#mamba){:target="_blank"} is a drop-in replacement for conda that works to improve upon some aspects of the conda infrastructure. It can very frequently perform installations much faster (helping loads with the "solving environment" steps). 

`conda` needs to be installed first like we did above, then we install `mamba` with `conda`. Trust me, this is likely very worth it for you, I do it on virtually every system I put conda on 🙂

So here's how we can install `mamba`, remember we want to be in the "base" conda environment when we run this (see just above), but we can also explicitly state that we want that environment as done here:

```bash
conda install -n base -c conda-forge mamba
```

Entering 'y' when prompted, and then forevermore we can use `mamba` in place of `conda` when *installing* things. **But it has seemed more stable/rubust to me to still use `conda` when activating/deactivating environments.** Example in this next section...

## Creating an environment and installing packages in one command with mamba

For this example, we'll create an environment with `mamba` while also specifying the channels needed and the program we want to install in one command. I'll shamelessly use my own package of bioinformatic tools, [bit](https://github.com/AstrobioMike/bioinf_tools#bioinformatics-tools-bit){:target="_blank"}  ¯\\\_(ツ)\_/¯

We want to be in our base environment, and then we can create the environment like so: 

```bash
mamba create -n bit -y -c conda-forge -c bioconda -c defaults -c astrobiomike bit
```

> **Breakdown**
> * `mamba create` – our base command (we just swapped `conda` with `mamba`)
> * `-n` – here is where we provide the name we want the environment to have
> * `-y` – says not to ask us for any confirmation
> * `-c ...` – each one of these specifies the channels in the appropriate order
> * `bit` – this positional argument is specifying the program name

This might take about about a minute, then we can activate it like so (note, as mentioned above I tend to use `conda` still for activate/deactivate, and `mamba` for all installations/env. creations):


```bash
conda activate bit
```

Then we can test it with `bit-version`, and run an example command like so:

```bash
bit-version

# downloading tiny example fasta file:
curl -L -o random-assembly.fa https://ndownloader.figshare.com/files/23842415

# and here is using a program in the bit package to summarize it
bit-summarize-assembly random-assembly.fa
```


<hr style="height:10px; visibility:hidden;" />

---
<br>


## Creating an environment from a yaml file
In the examples above, we created a new environment and provided the programs/versions we wanted in the actual command we typed out. Sometimes it's convenient to instead use a [yaml file](https://en.wikipedia.org/wiki/YAML){:target="_blank"} that holds all the information about the programs/versions we want, and we can just give that file to the `conda/mamba create` command. 

Here is an example of a [conda-formatted yaml](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-file-manually){:target="_blank"} file:

```yaml
name: py3.5-with-prodigal
channels:
    - conda-forge
    - bioconda
    - defaults
dependencies:
    - python=3.5
    - prodigal
```

Here, we are specifying the 'name', which becomes the name of the conda environment (meaning what we will give to the `conda activate ...` command), the channels (remember from above, that things from the bioconda channel should have [all 3 of those specified in that order](https://bioconda.github.io/#usage){:target="_blank"}), and then the "dependencies" are the programs and/or versions we want. For this example, just to show it, we are specifying python version 3.5, and no version specified for prodigal. 

We can quickly download that yaml file from above to our working location by running this command:

```bash
curl -L -o example-conda.yaml https://figshare.com/ndownloader/files/35843666
```

And here is how we can create an environment from a yaml file:

```bash
conda env create -f example-conda.yaml
```

> **Breakdown**
> * `conda env create` – our base command, when providing a yaml file, we need to add in the 'env' part here, unlike above where we just used `conda create` or `mamba create` (`mamba env create` would work here too)
> * `-f` – here is where we provide the file that holds the info on the environment we want to create


When that's done, we can activate it the same way:

```bash
conda activate py3.5-with-prodigal
```

And we will have the python version we specified, as well as prodigal in there.

Since we don't really need it, remember we can remove an environment like shown below, we just need to deactivate it first:

```bash
conda deactivate
conda env remove -n py3.5-with-prodigal
```

<hr style="height:10px; visibility:hidden;" />

---
---

>Again, [conda](https://conda.io/docs/){:target="_blank"} might be hard to appreciate if we haven't fought with installations, versions, and environments before, but I promise, it really is stellar 🙂

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
