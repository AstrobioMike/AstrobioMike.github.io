---
layout: main
title: bash basics
categories: [bash, tutorial]
permalink: /bash/basics
---

{% include _bash_basics_toc.html %}

{% include _side_tab_bash.html %}

Here we are going to cover the very basics of working in the terminal. For quick definitions of what things like *bash* and 'the terminal' are, see [here]({{ site.url }}/bash/).  

We'll start at the very beginning with the basic formula for how to run commands, then move on to how to navigate around your computer from within the terminal, and end with how to look at, create, and manipulate plain text documents. This is really designed for those with very little to no experience working at the command line, so if you feel comfortable with these tasks already, consider jumping ahead to one of the sections you're less familiar with.  

To start, let's open up a terminal window. If you have a Mac, this is easy as you can just do a spotlight search for "Terminal" and you're ready to rock. If you have a PC, I'm afraid this is where the trip ends for you:  
<img align="right" src="{{ site.url }}/images/oompa.jpg">  
<br>
<br>
<br>
<br>
<br>


Okay, not really. Things do get a little more complicated, but there are programs you can download to get the appropriate *bash* environment on a PC. Unfortunately I have no experience with that though so you'll have visit the almighty [google](https://www.google.com/search?source=hp&q=running+terminal+on+pc&oq=running+terminal+on+pc&gs_l=psy-ab.3..0i22i30k1l2.4571.9293.0.10104.25.22.0.0.0.0.175.1978.12j9.21.0....0...1.1.64.psy-ab..4.20.1827.0..0j0i131k1j0i10k1j33i22i29i30k1.RZTO4OlOhZk).  
<br>

---
<br>
# First things first!
For right now, and only for right now, I would like you to blindly copy and paste the commands here into your terminal window to set up the (very tiny) temporary files you'll need if you want to actively follow along below. Just reading through is of course fine, but if you are totally new to this, I'd recommend doing things with me. I promise this will be the only time you are just copying and pasting without an understanding of what you are doing. Besides, by the time you get to the end of this page, you'll have a much better understanding of what's going on in here.

For now, copy and paste this block of commands into your terminal window, (you can copy and paste all of them at once), and be sure to hit enter after pasting them to trigger the last command.

```
cd ~
curl -O https://AstrobioMike.github.io/tutorial_files/bash_basics_temp.tar.gz
tar -xvf bash_basics_temp.tar.gz
rm bash_basics_temp.tar.gz
cd bash_basics_temp
```

Now, let's get started!  
<br>

---
<br>

# Running a command
The terminal 'prompt' by default is typically a line with some information ending with a dollar sign:

<center><img src="{{ site.url }}/images/blank_prompt.png"></center>  

<br>
This line and the terminal itself are both customizable so things may look a little different on yours, but this is where the magic happens. 


To run a command in the terminal, you simply need to type the command after the prompt. Some commands don't require any further input and after being typed can be executed by hitting the return key. For example, there is a command called `date` that will output the time and date:

<center><img src="{{ site.url }}/images/terminal_date.png"></center>

<br>
And as you can see, after hitting return, a new prompt line appears and the cursor is waiting for the next command. 

The `date` command didn't require anything else, so it worked just fine by itself, but some other commands require what are known as arguments. In these cases the format would be `command argument` (separated by a space). **Some arguments are optional, and some are mandatory**. Probably the most common type of argument is telling *bash* which file you want to do something to. 

For example, the `head` command prints out the first few lines of a file to the terminal, so you need to tell it which file you want it to act on. Take for example this regular text file named "text.txt":

<center><img src="{{ site.url }}/images/text.txt.png"></center> 

<br>
You have this file now too, so if we run the `head` command on it, the terminal prints out the first 10 lines and then returns to the prompt:

<center><img src="{{ site.url }}/images/head_example.png"></center> 

<br>
In this case, providing a file to the `head` command is required; the mandatory argument for this is telling it what file we want it to act on. But this command also has optional arguments. For example, by default it prints out only the first 10 lines, but we can tell it to do however many lines we want by adding what's known as a 'flag'. These are often led by a single dash, followed by a character representing which argument you're specifying, followed by the value (or file) you want to give it. Here is an example where we tell the `head` command we want the first 15 lines of the document, overriding the default 10:


```bash
head -n 15 text.txt
```



<center><img src="{{ site.url }}/images/test.txt_head.png"></center>  

<br>
And that's really it. Those are the fundamentals that govern running virtually any individual command in the terminal. Many commands have a dizzying amount of options that I (and I presume most others) never fully appreciate or explore. You can pull up a 'manual' for any command by using the command `man` followed by the command you're intersted in. For instance, `man head` would bring up the manual for `head` (you can exit the manual by pressing `q` ). Many commands, though not all, also have a help menu of sorts that can usually be accessed by typing the command followed by `-h` or `--help`. I personally most often google for examples of whatever I'm trying to use. But importantly, you don't need to memorize these things. Some things might become second nature depending on how frequently you use them, but the real benefit comes from just knowing what can be done (meaning knowing what tools there are out there), and having a baseline understanding of how to work in the terminal environment. Then when you come across something you need to do, you know what to look for to iron out the details. 

Now that you're familiar with this baseline formula of entering commands, let's move on to looking at some more and learn how to navigate around your computer from inside the terminal.  
<br>

---
<br>
# Moving around

<h4><i>Commands presented in this section:</i></h4>

|Command     |Function          |
|:----------:|------------------|
|`pwd`       |tells you where you are in the computer (print working directory)|
|`ls`        |lists contents of a directory (list)|
|`cd`|changes directories|

Your computer stores files in a hierarchical structure like a tree. You are likely already used to this just by how you would navigate through by clicking on various folders (directories) and finding your way to a file or program. 

When you are working in the terminal, you are always sitting in some directory. Here is the directory we're currently working in, with some generic items inside it, as viewed from the Finder window:

<center><img src="{{ site.url }}/images/directory_example.png"></center>  

<br>
At the top of that image you can see we are in a directory called "bash_basics_temp", and inside that directory there is a subdirectory called "another_directory" and two text files. Additionally at the very bottom there is a line that tells us 'where' we are in the computer. Yours will be different, but this is what mine looks like this: `Macintosh_HD/Users/Mike_Lee/bash_basics_temp`

<br>
This line of directories, delimited by forward slashes, tells us where we are; it's an address. And in the computer world it's called a 'path'.

When we are working in a terminal we need to be aware of where we are in the computer. We can get all of the same information in the terminal by using the commands `pwd` (print working directory – to view the address of the directory we are in) and `ls` (list, to list the contents of the directory we're sitting in):

<center><img src="{{ site.url }}/images/terminal_directory_example.png"></center>  

<br>
It is important to be comfortable thinking about where you are in your computer when working in the terminal. One of the most common errors/easiest mistakes to make is trying to do something to a file that isn't where you think it is. Let's go back to our example above where we used the `head` command on the "text.txt" file, and then let's try it on another file:

```bash
head text.txt
head yet_another_text_file.txt
```

<center><img src="{{ site.url }}/images/file_location_error.png"></center>

<br>
Here the `head` command works fine on "text.txt", but we get an error message when we call it on "yet_another_text_file.txt". Intrepreting error messages in some cases will be tricky. As usual, google is your friend, and most problems you'll run into will be things others have already talked out somewhere – love you, [stackoverflow](https://stackoverflow.com). Fortunately this error message happens to be one of the more straightforward ones. It gives us the command that was used, the file we attempted to call it on, and tells us "No such file or directory". And if we enter the `ls` command just like we did above, we can see the computer is absolutely right (spoiler alert: it usually is). There is no file in the current directory named "yet_another_text_file.txt". And when you enter a file name without any other information, the computer only looks in the exact directory you are sitting in. We'll see what this means in a second. 

In this case the file we are looking for is actually in the directory "another_directory", which is a 'subdirectory' to the one we are sitting in, which we can also see when entering `ls` with no arguments. Further, if we instead enter:  

```bash
ls another_directory/
```

We are now providing an argument to the `ls` command, and asking it to list the contents of the directory named "another_directory".

<center><img src="{{ site.url }}/images/ls_another_dir.png"></center>

<br>
So we see the file we tried to call `head` on isn't in our current working directory, but is actually in a subdirectory just one layer deeper than we are. We can also call `head` on the file by specifying the path (address) of the file like so:

```bash
head another_directory/yet_another_text_file.txt
```

<center><img src="{{ site.url }}/images/head_other_dir_example.png"></center>

<br>
As a side note here, don't be distracted a the line wrap in the terminal window (as seen in mine here). Commands you are entering will sometimes get pretty long and will automatically wrap to the next line without interfering with what you're doing.  

Moving on, there are actually two ways to provide the path to where something is: you can give what is known as the **relative path** or you can give an **absolute path**.  

What we did in the example just above is known as a **relative path** because it takes off from where we entered the command. If we were in a location that didn't have the subdirectory "another_directory" in it, then we would have also gotten an error message because the computer again wouldn't be able to find the file we were pointing to.  

We also could provide the **absolute path** however, which isn't relative to our current location because it takes off from a specific location in the computer, rather than taking off from where we call the command. We actually saw our absolute path earlier when we called the `pwd` command:

<center><img src="{{ site.url }}/images/pwd.png"></center>

<br>
Now, instead of calling `head` on "yet_another_text_file.txt" by providing the relative PATH as we did last time, we will do it using the absolute path. But remember, what `pwd` just showed us is actually where *we* are, and we need to add the subdirectory "another directory" to the end of it in order to specify where the file actual is. To do this, we can copy and paste the output from `pwd`, and then add `/another_directory/yet_another_text_file.txt` to the end of it. Yours will be different, but this is what mine looks like: 

<center><img src="{{ site.url }}/images/head_abs_path_example.png"></center>

<br>
Now that we've covered how to do something to a file in a different location than our current working directory, let's look at how we can also just move ourselves to where that file is. We change directories with the `cd` command. Let's change directories into the subdirectory called "another_directory", check where we are with `pwd`, see if the file we are looking for is actually in the directory with `ls`, and then run `head` on it:

```
cd another_directory/
pwd
ls
head yet_another_text_file.txt
```

<center><img src="{{ site.url }}/images/cd_pwd_ls_head_example.png"></center>

<br>
Great. But now how do we get back 'up' to the directory above us? One way would be to provide the absolute path of where we would like to go, which for me would look like this: 
```
cd /Users/Mike_Lee/bash_basics_temp
``` 

But that would get old fast. A better way to do it involves using special characters:

```
cd ..
```

Here, simply entering two periods as the destination argument to the `cd` command tells the computer we want to go 'up' one level:

<center><img src="{{ site.url }}/images/cd_up_shortcut_example.png"></center>

<br>
Fortunately there are many special characters in *bash*, some of which allow us to navigate around much more easily. Two periods as we just saw refers to the directory just above you. Using just a single period specifies your current directory, which will be useful when we look at how to copy and move files around below. Another special one is simply a lone foreward slash `/`, which we've actually already seen when we provided the absolute path above (e.g. `/Users/Mike_Lee/bash_basics_temp`). What the `/` in the front is actually telling the computer is to start at the "root" directory. This is sort of like the home base of the operating system structure; it is the 'highest' directory level. If you type `ls /`, you will get a list of the directories and files located in the root directory. And if you wanted, you could move yourself into the root directory by entering `cd /`.  

Another special character related to location is the `~` symbol. This shortcut points to your "home" directory. A home directory is like a more personal spot on the computer. If you're on a server with multiple users, it would be your own personal location. And even if you're just working on your own computer, you likely login as a specific user when you turn it on. For instance, my home directory is located at `/Users/Mike_Lee`. We can see that if I change into my home directory and run `pwd`:

<center><img src="{{ site.url }}/images/home_dir.png"></center>

<br>
Another extremely useful one for the `cd` command is simply a dash (`-`). This will change you back to the last directory you were in:

<center><img src="{{ site.url }}/images/change_dir_back.png"></center>

<br>
Having some concept of where you are and how to navigate around the computer via the terminal window alone is an essential skill that you'll develop very quickly. If you'd like, at first you can practice by also having a Finder window open and try 'clicking' around to the same places you are moving through with *bash* commands to help you visualize the structure.  
<br>

---
<br>
# What is a plain text file?

The baseline tools of *bash* are mostly useful for what are known as plain text files, also commonly referred to as 'flat' files. And I'm just realizing now that formally defining what a 'plain text file' is isn't all that simple. There are a few definitions you could check out at the [wiki](https://en.wikipedia.org/wiki/Plain_text) if you are interested, but a simple, working defintion might be something like: a text file that doesn't contain any special formatting characters and can be properly viewed and edited with any standard text editor.

Bioinformaticians (and lots of others who get to play with big data regularly) work with plain text files so much because not being constrained to any special (arbitary) characteristics means methods for interacting with them can be standardized. So when you learn one command that works on a flat file, it works on all of them the same way. Common formats are files with extensions like ".txt", ".csv" for comma-separated values, ".tsv" for tab-seperated values. More specialized extensions like ".docx" or ".xlsx" are not plain text files. Delimited files, like ".csv" and ".tsv", are a simple way to store tables as each row is delimited by a newline and each column by whichever delimiter is specified. And it's this simple formula that makes them very easy to work with. In the next section [(6 commands worth getting to know)](/bash/six_commands) you'll see and work with some examples of tables, and see why *bash* is invaluable for manipulating them. But for now just know that the text files we've been looking at so far are plain text files.  
<br>

---
<br>
# Working with plain text files and directories

<h4><i>Commands presented in this section:</i></h4>

|Command     |Function          |
|:----------:|------------------|
|`head`      |prints the first few lines of a file|
|`tail`      |prints the last few lines of a file|
|`less`      |allows you to browse a file (exit with "q" key)|
|`wc`       |count lines, words, and characters in a file|
|`cp`      |copy a file or directory (use with caution)|
|`mv`      |mv a file or directory (use with caution)|
|`rm`      |delete a file or directory (use with caution)|
|`mkdir`       |create a directory|
|`rmdir`     |delete a directory|
|`nano`     |create and edit plain text files|


<h4>Ways to probe plain text files</h4>
We've already used a very common tool for peeking at files, the `head` command. There is also a corresponding `tail` version that prints the last 10 lines of a file by default:

<center><img src="{{ site.url }}/images/tail_ex.png"></center>

<br>
Commands like `head` and `tail` are particularly useful when you have a very large file that might take a lot of memory to fully open (and therefore could be sluggish). Usually a file is formatted consistently and often just peeking at the first few lines tells you the structure, so you can then pull out just what you are interested in. Along the same lines when you are trying to parse a large file a certain way, you can just do a subset of it with the `head` command to make sure things are working before running it on the entire file. 

Another useful command for just viewing a file is `less`. This opens a searchable reader that allows you to scroll through the document if you just want to read something over. Our example documents are kind of small for `less` to be useful here, but it would be run as such:

```bash
less text.txt
```
<center><img src="{{ site.url }}/images/less_ex.png"></center>

<br>
To exit the `less` program you need to press the "q" key. 

The `wc` command is useful for counting how many lines, words, and characters there are in a file (wanting this information comes up more than you'd expect). For example, running it on the "text.txt" file with no options specified gives us all three:

<center><img src="{{ site.url }}/images/wc_ex.png"></center>
<br>
I personally find myself most often using the `wc` command with the optional argument `-l`, which tells it I only want to know the number of lines in a file:

<center><img src="{{ site.url }}/images/wc_line_ex.png"></center>
<br>

<h4>Ways to manipulate files and directories</h4>

<div class="warning">
<center><h2>WARNING!</h2></center>
Using commands that do things like create, copy, and move and rename files/directories in the terminal <b>will overwrite</b> files/directories that already exist <b>if they have the same name</b>. And using commands that delete things will by default do so without any warning or confirmation. Caution is required until you get used to working with them – and then forever after.
</div>


The commands `cp` and `mv` (copy and move) both function under a similar syntax. The command entered needs to be followed by 2 **positional arguments**. Positional arguments are arguments that understood by the computer to be something specific based on where they come following the command. In the case of the `cp` and `mv` commands the first positional argument is the file you want to act on (the source), and the second positional argument is where you want it to go (the target). 

Let's take a quick look again at what's in our current working directory with `ls`:

<center><img src="{{ site.url }}/images/ls_ex.png"></center>
<br>
We see there are currently 3 items – a directory and two text files. Now, let's make a copy of "text.txt":

```bash
cp text.txt text_copy.txt
```
And now when we list the files in our current working directory they are both there: 

<center><img src="{{ site.url }}/images/cp_ex.png"></center>

<br>
Similarly, we can specify our source or target (first or second arguments) to be somewhere other than our current working directory. Here we are going to make a copy of the file "yet_another_text_file.txt" from the subdirectory "another_directory" and put the copy in our current working directory. To specify the current working directory as the target location, we can simply provide a single period ( `.` ) as the target:

```bash
cp another_directory/yet_another_text_file.txt .
```

<center><img src="{{ site.url }}/images/cp2_ex.png"></center>

<br>
And now we have a copy of that file in our current working directory. 

Notice that in the first `cp` example we provided a new name in the target location, ("text_copy.txt"), while in the second example when we used `.` to specify copying to the current working directory it simply copied the file with the same name as the original. If we had wanted to copy to our current location but also give it a new name, we would run it with the desired new name as the target (second argument):

```bash
cp another_directory/yet_another_text_file.txt yet_another_text_file_copy.txt
```

<center><img src="{{ site.url }}/images/cp3_ex.png"></center>

<br>
Now we have the original copy we made, "yet_another_text_file.txt", and the copy that we renamed, "yet_another_text_file_copy.txt". Note again that if we don't provide a path, relative or absolute, the computer looks in the current working directory.

The `mv` command is used to move files **and** to rename files, and as mentioned it works the same as the `cp` command – requiring a source argument and a target argument that need to be entered in that order. Here we will move a file from our current working directory into our subdirectory, and demonstrate how easy it is to accidentally overwrite something.  

Currently our current working directory and our subdirectory contain these files:  

<center><img src="{{ site.url }}/images/ls_both_dir_ex1.png"></center>

<br>
After we move the "text.txt" file as follows, it is only present in the subdirectory:

```bash
mv text.txt another_directory/
```

<center><img src="{{ site.url }}/images/mv_ls_ex.png"></center>

<br>
Note here, that we didn't provide a file name for the target of the `mv` command, only the location of the subdirectory we wanted to move it into. Just like when we used `cp` to copy a file to our current working directory with a `.` above, this results is keeping the name of the original source file.  

Let's take a quick look again at what's in these two files that are now in the subdirectory:

<center><img src="{{ site.url }}/images/mv_overwrite_ex1.png"></center>

<br>
And now let's see an example of how easy it is to overwrite a file accidentally when using the `mv` command, by specifying the name of the output to be "yet_another_text_file.txt" (first we are getting a copy of "text.txt" back into our working directory, as we just moved it):

```bash
cp another_directory/text.txt .
mv text.txt another_directory/yet_another_text_file.txt
```

<center><img src="{{ site.url }}/images/mv_overwrite_ex2.png"></center>

<br>
Note how the contents of "yet_another_text_file.txt" have changed. It is now the same as "text.txt". Its original contents would be lost forever if this were the only location the file were saved in. 

To delete files (intentionally) there is the `rm` command (remove). This requires at least one argument specifying the file you want to delete. For an example, we'll delete the "yet_another_text_file_copy.txt" file:

```bash
rm yet_another_text_file_copy.txt
```

<center><img src="{{ site.url }}/images/rm_ex.png"></center>

<br>

<h3>How to make and delete directories</h3>

To create new directory we use the command `mkdir` followed by one argument for the directory name:

```bash
mkdir our_new_directory
```
<center><img src="{{ site.url }}/images/mkdir_ex.png"></center>

<br>
And similarly, directories can be deleted with `rmdir`:

```bash
rmdir our_new_directory
```
<center><img src="{{ site.url }}/images/rmdir_ex.png"></center>

<br>
Though if the target directory is not empty, `rmdir` will give you an error as a (rare) safety measure. There are ways to force this action or to use the regular `rm` command on directories that you can find by looking further into the commands and the optional arguments you can provide them.

<h3>Making and editing plain text files</h3>
It is often very useful to be able to generate new plain text files quickly at the command line, or make some changes to an existing one, and there are many ways to do this. One of the many ways involves using a text editor that operates on the command line, and there are several common text editors that fit this bill. I think the easiest to use at first is `nano`. So here we're just going to go over a quick example of that, but you should know there are more and better out there if you're willing to dedicate time to *another* steep learning curve. (I'm still dragging my feet a bit on that one too, so don't feel bad.)  

When we run the command `nano` it will open the program (the "nano" text editor), and our terminal window changes from our regular view, with our prompt that we've seen so far, to a text editor interface. The command itself can be run without any arguments and it will simply open a new file that you can then modify and save before you exit. But I like to start the editor with a file name, which is what the first positional argument is if you give it one:

```bash
nano test.txt
```
<center><img src="{{ site.url }}/images/blank_nano_ex.png"></center>

<br>
Now our view has changed as we are in the program. You can see in mine there is a header bar telling me which version of "nano" I am in and the name of the file we are working on (I gave it "test.txt"). And there are some keyboard shortcuts listed at the bottom. 

Now we can simply type as we normally would to add a bit of text:  

<center><img src="{{ site.url }}/images/nano_text_ex.png"></center>

<br>
And then there are a few ways to save and exit. For me, I press `Ctrl + x`, at which point you will be asked if you want to "Save modified buffer?". Here pressing the `y` key signifies yes we want to save. And then it asks what we'd like to name the file, and has our current file name in place already. So long as we don't want to change the name of the file, we can just hit `return`.  

So altogether one way to get out of "nano" and save the file: `Ctrl + x`, `y`, `return`.  

Now we can see that we've made the new file with `ls` and take a peek at it with `head`:

<center><img src="{{ site.url }}/images/nano_head_ex.png"></center>

<br>
<br>

---
<br>
# Pipes and redirectors
Here we're going to quickly touch upon what makes the Unix command-line environment so powerful. A pipe `|`, is used to connect commands. Basically it takes in the output from the previous command and pipes it into the input of the following command. Here is a simple example of this where we are going to list the contents of the current working directory with `ls`, and then count how many items there are by counting the number of lines with `wc -l`:


```bash
ls | wc -l
```


A simple example using `>` is we're going to list the contents of the current working directory with `ls` and redirect it into a file (rather than printing it to the terminal). 

<br>
<br>

---
<br>
# Tab-completion is your friend!

There is a very important habit you need to develop if you haven't yet, and it is using the `tab` key to complete files and directories at the terminal. You may have noticed there is quite a bit of typing out file names, and spanning multiple directories if what you need to point to is located somewhere other than your current working directory. Fortunately, we don't need to type things out all the way.  

If we are trying to point to a file that's in our current working directory, we can begin typing the filename and then press `tab` to complete it. And this works the same for directories if we are trying to get to a file that is located somewhere else. Whether or not it will complete the file or directory depends on if there is more than one that starts with the same characters as the those you've entered so far.  

Let's look at this example where we want to run `head` on the file "test.txt" like we did above:

<center><img src="{{ site.url }}/images/tab_complete_ex.png"></center>
<br>
Here, I typed `head` then the first two letters of the file I wanted "te", then I hit `tab`. At first it just gave me a notification sound and nothing happened (depending on your setup you may not hear a notification sound). But then I pushed it a second time and it revealed all of the possible files I might be trying to specify. In this case that is two: "test.txt" and "text_copy.txt". I then added an "s" and pressed `tab` again, and that time it finished the file name for me as there was no longer any ambiguity about which file I may mean.  

This may seem like a simple convenience at first, and it certainly is, but tab-completion is also incredibly valuable for making sure you aren't entering anything wrong. We silly humans like to make mistakes like crazy. Any time we can take human-error out of the equation, we should. If you spend any time at the terminal you will quickly realize how invaluable it is to tab-complete everywhere you can, as it assures that things are where you think they are (because if something doesn't tab-complete, it means you are in the wrong location).
<br>
<br>

---
<br>
<h1>Congrats on getting through the basics!</h1>
Some of these things may seem trivial as used in these examples here, and I know that can make the effort required to learn to use them harder to come by, but these basics really are the foundation to do everything at the command line. Understanding the general rules of commands, their arguments, and the required syntax will make everything you try to do easier.  

There are a lot of base commands in *bash*, and a dizzying number of optional arguments for most of them. Again, google is our friend. If you end up working in the terminal window enough, you will remember some things, but also you will often do a quick search to remember what the flag is for a specific argument, or how exactly a specific command works. This really isn't about memorization overall.

Remember the set of commands you ran at the beginning to get the files needed for this walkthrough?

```
cd ~
curl -O https://AstrobioMike.github.io/tutorial_files/bash_basics_temp.tar.gz
tar -xvf bash_basics_temp.tar.gz
rm bash_basics_temp.tar.gz
cd bash_basics_temp
```

I promised you'd have a better idea of what's going on in here by the end, so let's break it down.

We start by changing directories to our 'home' location, and we do this with the special character `~`. Then we use two commands you aren't familiar with yet, `curl` and `tar`. We won't get into these here other than to say `curl` allows you to download things from the internet from the command line, and `tar` is a tool for packing and unpacking files and directories. So with the `curl` command we download the files in a compressed format, and then with `tar` and some options we expand them into a directory called "bash_basics_temp". Then all we do is delete the uncompressed files ("bash_basics_temp.tar.gz"), and change directories into "bash_basics_temp".  

As you can see here, while the commands change, the general structure of how to operate in the command line stays the same. This is part of why the learning curve is steep at the start but it gets better. Once you become comfortable with the general framework, everything is gets easier.  

Now that you're comfortable with the basics, head on over to [6 commands worth getting to know](/bash/six_commands) for a walkthrough of some stellar commands I wouldn't want to live without.
