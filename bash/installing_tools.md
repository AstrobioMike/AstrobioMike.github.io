---
layout: main
title: Installing tools at the command line
categories: [bash, tutorial]
permalink: /bash/installing_tools
---

{% include _bash_install_tools_toc.html %}

{% include _side_tab_bash.html %}

As we've seen in the [basics](/bash/basics), [six glorious commands](/bash/six_commands), and [real-life examples](/bash/why) pages, you can do some amazing things with *bash* as far as manipulating [plain text files](/bash/basics#whats-a-plain-text-file) goes. But that's just the start of things. If you're working towards more bioinformatic-leaning applications, you're certainly also going to need to download and/or install the lots of tools that don't come standard. What's required on our end to get something working properly varies by the tool, and most that are highly used by people have excellent documentation on how to properly download and install them. This page will be an ongoing list of the installation process for any of the tools used on tutorials from this site, and they will hopefully serve as examples to help guide you through installing other things.  
<br>

---
---
<br>
# What's a binary?
You may have heard or read an expression before referring to whether there is a *binary* or *executable* available for download for a particular program. Programs stored in these formats can just be downloaded and then utilized usually with no other work required – they are already "built" when you download them. This is in contrast to programs that you need to download and then "compile from source" – you need to "build" the program first in order for it to work. As usual there is a give and take here, binaries are easier to grab and get running, but you have more freedom and often access to the newest developments when installing something from its source code. There is a good explanation on [stack exchange](https://unix.stackexchange.com/questions/152346/what-is-the-difference-between-building-from-source-and-using-an-install-package) about some of these differences. The examples here so far are all binaries (as I've been selecting programs like that on purpose for ease of use in tutorials), but I will try to get an example up soon of installing from source.  
<br>

---
<br>
# A happy bin
Often a big part of getting things to work properly is having them in a location on the computer that you can access no matter [where you are](/bash/basics#moving-around). As we covered [here](/bash/modifying_your_path), a list of directories that are scanned for programs automatically by your computer is stored in the *bash* special variable called "PATH". For the sake of simplicity, any tools we use in tutorials on this site are going to be put in a directory we made and added to our PATH in the [modifying your PATH walkthrough](/bash/modifying_your_path) called `~/happy_bin`. You can check to see if this is in your PATH already like this:

```bash
echo $PATH | tr ":" "\n" | grep "happy_bin"
```

If this returns something like this:

<center><img src="{{ site.url }}/images/checking_for_happy_bin.png"></center>

<br>
You're good to go. If nothing was returned, then you can create this directory and add it to your PATH like such:

```bash
mkdir ~/happy_bin
export PATH="$PATH:/Users/Mike_Lee/happy_bin" >> ~/.bash_profile
```

And if you're unsure of what's going on here, be sure to visit the [modifying your PATH page](/bash/modifying_your_path).  
<br>

---
<br>
# Installing tools used here
Everything we install here we will put in our `~/happy_bin`, feel free to change that location to where you'd like of course, but then you will need to modify any code accordingly.  

## vsearch
There are instructions to get vsearch up and running [on its github](https://github.com/torognes/vsearch), but these commands should work for you if you're on a Mac **(if you're not, you'll have to download a different version you can find following the above link)**.

```bash
cd ~/happy_bin
curl -LO https://github.com/torognes/vsearch/releases/download/v2.5.1/vsearch-2.5.1-macos-x86_64.tar.gz
tar -xzvf vsearch-2.5.1-macos-x86_64.tar.gz
cp vsearch-2.5.1-macos-x86_64/bin/vsearch .
rm vsearch-2.5.1-macos-x86_64.tar.gz
```

Here we changed into our `~/happy_bin` directory, downloaded the vsearch tool with `curl`, unpacked it and unzipped things with `tar`, copied the main executable file into our `~/happy_bin` directory so that it is in our [PATH](/bash/modifying_your_path) and can be called from anywhere, then finally we deleted the compressed downloaded file.  

## usearch
To get the free version of usearch, you first need to go to [https://www.drive5.com/usearch/download.html](https://www.drive5.com/usearch/download.html), and fill out a (very) short form in order to have a download link sent to you. This usually happens virtually instantly. After getting the link and downloading the file, assuming you're working on a Mac and you downloaded the same version as noted above (v10.0.240) into your default download directory, the following commands will move usearch into our `~/happy_bin` directory and change its properties so we can execute it (if you download a different version, adjust the commands accordingly):

```bash
cd ~/happy_bin
mv ~/Downloads/usearch10.0.240_i86osx32 usearch
chmod +x usearch
```

Here we changed into our `~/happy_bin` directory, moved the downloaded file from our downloads directory to here and renamed it to simply "usearch" with the `mv` command, and then ran the `chmod` command to make the program executable so we can actually run it. Sometimes this `chmod` step is needed and sometimes it isn't depending on how the developer packaged the program. If we had tried to run `usearch` before running the `chmod` command, we would have gotten an error saying permission denied. We can check the permissions and properties ascribed to a file or directory by using the `ls` command with the `-l` optional argument: 

<center><img src="{{ site.url }}/images/ls_l_permission.png"></center>

<br>
The details of what all of the letters and their positions mean here are a bit further into the weeds than we need to worry about right now, but quickly, an "r" is for readable, a "w" is for writable, and an "x" is for executable. Note that the vsearch executable file has an "x" in the position I highlighted for the usearch file. After running `chmod +x usearch`, the output of `ls -l` looks like this:

<center><img src="{{ site.url }}/images/ls_l_permission2.png"></center>

<br>
And now the usearch program can be executed (run) without throwing any errors.  
<br>

---
---
<br>
# More to come
I'll be adding more examples here for each program we use, and feel free to contact me about any tricky ones you run into that might serve as good examples 🙂
