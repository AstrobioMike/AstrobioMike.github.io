---
layout: main
title: Getting a Unix environment
categories: [unix]
permalink: /unix/getting_unix_env
---

{% include _side_tab_unix.html %}

# Mac or Linux
If you are working on a Mac or Linux computer, this is already taken care of and you can just do a search for "terminal".  

<hr style="height:10px; visibility:hidden;" />

---
<br>
# Windows
If you are working on a Windows computer, you will need to set up a working command-line environment.

#### Windows 10 or later
* If your operating system is Windows 10 or later, the operating system includes a way to run a command-line environment. You can follow the [helpful steps outlined here](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/){:target="_blank"} to get set up.


#### Earlier than Windows 10
* If you have an earlier Windows operating system, there are a few options. I've seen good results using [Git for Windows](https://gitforwindows.org/){:target="_blank"}, so if you're open to suggestions of what to try, I'd say start with that. To give that a shot you can follow these instructions (updating if there is a newer version available).  
  * Download the "Git-2.18.0-32-bit.exe" file from the [Git for Windows download page](https://github.com/git-for-windows/git/releases/tag/v2.18.0.windows.1){:target="_blank"} (look near the bottom of the page in the "Assets" section).  
  * After it is finished downloading, run the installer by opening the file and proceed through the installation, paying attention to these notes:  
    * installing in the default folder location is fine
    * for "Which components should be installed?", make sure the following boxes are checked: "On the Desktop"; "Git Bash Here"; "Git GUI Here"; "Associate .git* configuration files with the default text editor"; and "Associate .sh files to be run with Bash"
    * the shortcuts default location is fine
    * change the default Git editor to "Nano"
    * on the "Adjusting your PATH environment" screen, select "Use Git from Git Bash only"
    * just click "Next" on the remaining configuration windows, and "Install" at the final one
  * When the installation is finished, you should be able to open a terminal window by launching Git Bash from your desktop. 
