---
layout: page
categories: [bash, tutorial]
exclude_from_nav: true
---

<h2><i>bash</i> basics</h2>

Here we are going to cover the very basics of working in the terminal. For quick definitions of what things like *bash* and 'the terminal' are, see [here]({{ site.url }}/bash.html).  

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
<center><img src="{{ site.url }}/images/under_construction.jpeg"></center>
<center><h3>UNDER CONSTRUCTION</h3></center>
<br>


### <u>First things first!</u>
For right now, and only for right now, I would like you to blindly copy and paste the commands here to set up the tiny temporary files you'd need if you want to actively follow along below. Just reading through is of course fine, but if you are truly new to this, I'd recommend doing things with me. I promise this will be the only time you are just copying and pasting without an understanding of what you are doing. Also, I like the idea of you running these commands to generate the files, rather than just downloading them from somewhere, because by the time you reach the end of this page you will understand everything you're doing here anyway. So it will be nice for you to look back at this afterwards and have a better understanding of how these steps work.  

For now, copy and paste this block of commands into the terminal window, you can copy and paste all of them at once:

```bash
cd ~
mkdir happy_belly_temp
cd happy_belly_temp
nano


### <u>Running a command at the prompt</u>
The terminal 'prompt' by default is typically a line with some information ending with a dollar sign:

<center><img src="{{ site.url }}/images/blank_prompt.png"></center>

<br>
This line and the terminal itself are both customizable so things may look a little different on yours, but this is where the magic happens. 

Throughout this tutorial I will be presenting blocks of code that will follow a dollar sign:


```bash
$ 
```


That dollar sign just represents your prompt, and shouldn't be typed when replicating the commands being used.

To run a command in the terminal, you simply need to type the command after the prompt. Some commands don't require any further input and can be executed by hitting the return key. For example, there is a command called `date` that will output the time and date:
<<<<<<< HEAD

<center><img src="{{ site.url }}/images/terminal_date.png"></center>
=======
<center><img src="{{ site.url }}/images/terminal_date.png"></center>

![Alt text]({{ site.url }}/images/terminal_date.png "date")
>>>>>>> 048bdf81a01ab4b4b5f5bd7e9574424456d6f9c4

<br>
And as you can see, after hitting return, a new prompt line appears and the cursor is waiting for the next command. 

The `date` command didn't require anything else, so it worked just fine by itself, but some other commands require what are known as arguments. In these cases the format would be `command` `argument` (separated by a space). Some arguments are optional, and some are mandatory. Probably the most common type of argument is telling *bash* which file you want to do something to. 

For example, the `head` command prints out the first few lines of a file to the terminal, so you need to tell it which file you want it to act on. Take for example this regular text file named 'test.txt':

<center><img src="{{ site.url }}/images/text.txt.png"></center> 

<br>
If I run the 'head' command on this file, the terminal prints out the first 10 lines and then returns to the prompt:
<center><img src="{{ site.url }}/images/head_example.png"></center> 
<br>
In this case, providing a file to the `head` command is required; a mandatory 'argument' for this is telling it what file we want it to act on. But this command also has optional arguments. For example, by default it prints out only the first 10 lines, but we can tell it to do however many lines we want by adding what's known as a 'flag'. These are often led by a single dash, followed by a character representing which argument you're specifying, followed by the value (or file) you want to give it. Here is an example where we tell the `head` command we want the first 15 lines of the document, overriding the default 10:


```bash
$ head -n 15 text.txt
```



<center><img src="{{ site.url }}/images/test.txt_head.png"></center>  

<br>
And that's really it. Those are the fundamentals that govern running virtually any individual command in the terminal. Now that you're familiar with this baseline formula, let's move on to looking at some more commands and learn how to navigate around your computer from inside the terminal.  
<br>

---
<br>
### <u>System structure and moving around your computer within the terminal</u>

#### <i>Commands presented in this section:</i>  

|Command     |Function          |
|:----------:|------------------|
|`pwd`       |tells you where you are in the computer (print working directory)|
|`ls`        |lists contents of a directory (list)|
|`cd`|changes directories|

Your computer stores files in a hierarchical structure like a tree. You are likely already used to this just by how you would navigate through by clicking on various folders (directories) and finding your way to a file. 

When you are working in the terminal, you are always sitting in some directory. Here is the directory I happen to be working in, with some generic items inside it, as viewed from the Finder window:

<center><img src="{{ site.url }}/images/directory_example.png"></center>  

<br>
At the top of that image you can see we are in a directory called "temp", and inside that directory there is a subdirectory called "another_directory" and two text files. Additionally at the very bottom there is a line that tells us 'where' we are in the computer:

<center><code>Macintosh_HD/Users/Mike_Lee/Documents/web_content/temp</code></center>

<br>
This line of directories, delimited by forward slashes, tells us where we are; it's an address. And in the computer world it's called a PATH.

When we are working in a terminal we need to be aware of where we are in the computer. We can get all of the same information in the terminal by using the commands `pwd` (print working directory – to view the address of the directory we are in) and `ls` (list, to list the contents of the directory we're sitting in):

<center><img src="{{ site.url }}/images/terminal_directory_example.png"></center>  

<br>
It is important to be comfortable thinking about where you are in your computer when working in the terminal. One of the most common errors/easiest mistakes to make is to be trying to do something to a file that isn't where you think it is. Let's go back to our example above where we used the `head` command on the "text.txt" file, and then let's try it on another file:

```bash
$ head text.txt
$ head yet_another_text_file.txt
```

<center><img src="{{ site.url }}/images/file_location_error.png"></center>

<br>
Here the `head` command works fine on "text.txt", but we get an error message when we call it on "yet_another_text_file.txt". Intrepreting error messages in some cases will be tricky. As usual, google is your friend, and most problems you'll run into will be things others have already talked out somewhere – love you, [stackoverflow](https://stackoverflow.com). Fortunately this error message happens to be one of the more straightforward ones. It gives us the command that was used, the file we attempted to call it on, and tells us "No such file or directory". And if we enter the `ls` command just like we did above, we can see the computer is absolutely right (spoiler alert: it usually is). There is no file in the current directory named "yet_another_text_file.txt". And when you enter a file name without any other information, the computer only looks in the exact directory you are sitting in.  

In this case the file we are looking for is actually in the directory "another_directory", which is a 'subdirectory' to the one we are sitting in, which we can also see when entering `ls` with no arguments. Further, if we instead enter:  

```bash
$ ls another_directory/
```

We are now providing an argument to the `ls` command, and asking it to list the contents of the directory named "another_directory".

<center><img src="{{ site.url }}/images/ls_another_dir.png"></center>

<br>
So we see the file we tried to call `head` on isn't in our current working directory, but is actually in a subdirectory just one layer deeper than we are. We can also call `head` on the file by specifying the PATH (address) of the file like so:

```bash
$ head another_directory/yet_another_text_file.txt
```

<center><img src="{{ site.url }}/images/head_other_dir_example.png"></center>

<br>
There are two ways to provide the PATH to where something is: you can give what is known as the "relative PATH"; or you can give an "absolute PATH".  

What we did in the example just above is known as a relative PATH because it takes off from where we entered the command. If we were in a location that didn't have the subdirectory "another_directory" in it, then we would have also gotten an error message because the computer again wouldn't be able to find the file we were pointing to.  

We also could provide the absolutely PATH however, which isn't relative to our current location because it takes off from a specific location in the computer, rather than taking off from where we call the command. We actually saw our absolute PATH earlier when we called the `pwd` command:

<center><img src="{{ site.url }}/images/pwd.png"></center>

<br>
Now, instead of calling `head` on "yet_another_text_file.txt" by providing the relative PATH as we did last time, we will do it using the absolute PATH. But remember, what `pwd` just showed us is where *we* are, and we need to add the subdirectory "another directory" to the end of it in order to be specifying where the file actual is:

```bash
$ head /Users/Mike_Lee/Documents/web_content/temp/another_directory/yet_another_text_file.txt
```

<center><img src="{{ site.url }}/images/head_abs_path_example.png"></center>

<br>
Don't be distracted by the line wrap in the PATH. Commands you are entering will sometimes get pretty long and will automatically wrap to the next line without interfering with what you're doing.  


Now that we've covered how to do something to a file in a different location than our current working directory, let's look at how we can also just move ourselves to where that file is. We change directories with the `cd` command. Let's change directories into "another_directory", check where we are with `pwd`, check that the filing we are looking for is actually in the directory with `ls`, and then run `head` on it:

```
$ cd another_directory/
$ pwd
$ ls
$ head yet_another_text_file.txt
```

<center><img src="{{ site.url }}/images/cd_pwd_ls_head_example.png"></center>

<br>
Great. But now how do we get back 'up' to the directory above us? One way would be to provide the absolute PATH of where we would like to go:

```
cd /Users/Mike_Lee/Documents/web_content/temp
```

But that would get old fast. A better way to do it involves special characters:

```
cd ..
```

Here, simply entering two periods as the destination argument to the `cd` command tells the computer we want to go 'up' one level:

<center><img src="{{ site.url }}/images/cd_up_shortcut_example.png"></center>

<br>
Fortunately there are many special characters in *bash*, some of which allow us to navigate around much more easily. Two periods as we just saw refers to the directory just above you. Another special one is simply a lone foreward slash `/`, which we've actually already seen when we provided the absolute PATH above. What the `/` is actually telling the computer is to start at the "root" directory. This is sort of like the home base of the operating system structure; it is the 'highest' directory level. If you type `ls /`, you will get a list of the directories and files located in the root directory. And if you wanted, you could move yourself into the root directory by entering `cd /`.  

Another special character related to location is the `~` symbol. This shortcut points to your "home" directory. A home directory is like a more personal spot on the computer. If you're on a server with multiple users, it would be your own personal location. And even if you're just working on your own computer, you likely login as a specific user when you turn it on. For instance, my home directory is located at `/Users/Mike_Lee`. We can see that if I change into my home directory and run `pwd`:

<center><img src="{{ site.url }}/images/home_dir.png"></center>

<br>
Another really useful one for the `cd` command is simply a dash (`-`). This will change you back to the last directory you were in:

<center><img src="{{ site.url }}/images/change_dir_back.png"></center>

<br>
Having some concept of where you are and how to navigate around the computer via the terminal window alone is an essential skill that you'll develop very quickly. If you'd like, at first you can practice by also having a Finder window open and try 'clicking' around to the same places you are moving through with *bash* commands to help you visualize the structure.  
<br>

---
<br>
### <u>Probing files</u>

#### <i>Commands presented in this section:</i>  

|Command     |Function          |
|:----------:|------------------|
|`head`      |prints the first few lines of a file|
|`tail`      |prints the last few lines of a file|
|`less`      |allows you to browse a file (exit with "q" key)|
| `wc`       |count lines, words, and characters in a file|



less, head, tail, wc
<br>

### Manipulating files
cp, mv, rm, mkdir, rmdir 
<br>

### Making plain text files
nano, vim, emacs

{% highlight rouge %}
$ for i in `cat samples`; do echo $i; wc -l "$i".txt; done
{% endhighlight %}


```bash
$ for i in `cat samples`; do echo $i; wc -l "$i".txt; done
```
