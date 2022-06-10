---
layout: main
title: Getting a Unix-like environment
categories: [unix]
tags: [unix,bash,bioinformatics,tools,program,installation,tutorial]
permalink: /unix/getting_unix_env
---

{% include _side_tab_unix.html %}

> Many of the tutorials on this site have [binder](https://mybinder.org/){:target="_blank"} environments available to be able to jump right in. But here are some ways to access a Unix-like environment on your system. 

# Mac or Linux
If you are working on a Mac or Linux computer, this is already taken care of and you can just do a search for "terminal".  

<hr style="height:10px; visibility:hidden;" />

---
<br>

# Windows
If you are working on a Windows computer, you will need to set up a working command-line environment. There are many ways to do this, below are listed just two if you are looking for suggestions :) 

## Windows 10 or later
* If the operating system is Windows 10 or later, it includes a way to run a Unix-like command-line environment. You can find info on how to install that [here](https://docs.microsoft.com/en-us/windows/wsl/install).

## Earlier than Windows 10
* If the operating system is an earlier version, one way you can set up a command-line working environment is through [Git for Windows](https://gitforwindows.org/). Here are the steps that I followed when trying this with the version noted below.  
  * Download the "Git-2.22.0-64-bit.exe" file from the Git for Windows [download page](https://github.com/git-for-windows/git/releases/tag/v2.22.0.windows.1) (scroll down a little bit to "Assets")  
  * After it is finished downloading, run the installer by opening the file and proceed through the installation:  
    * installing in the default folder location is fine
    * for "Which components should be installed?", make sure the following boxes are checked: "On the Desktop"; "Git Bash Here"; "Git GUI Here"; "Associate .git* configuration files with the default text editor"; and "Associate .sh files to be run with Bash"
    * the shortcuts default location is fine
    * change the default Git editor to "Nano"
    * on the "Adjusting your PATH environment" screen, select "Use Git from Git Bash only"
    * just click "Next" on the remaining configuration windows, and "Install" at the final one
  * When the installation is finished, you should be able to open a terminal window by launching Git Bash from your desktop. 
