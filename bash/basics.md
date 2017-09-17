---
layout: page
categories: [bash, tutorial]
exclude_from_nav: true
---

<h2><i>bash</i> basics</h2>

Here we are going to cover the very basics of working in the terminal. For quick definitions of what things like *bash* and 'the terminal' are, see [here]({{ site.url }}/bash.html).

We'll start at the very beginning with the basic formula for how to run commands, then move on to how to navigate around your computer from within the terminal, and end with how to manipulate and create documents. This is probably suitable for those with little to no experience, so if you feel comfortable with these tasks already, consider jumping ahead to one of the sections you're less familiar with. To start, let's open up a terminal window. If you have a Mac, this is easy as you can just do a spotlight search for Terminal and you're ready to rock. If you have a PC, I'm afraid this is where the trip ends for you:  
<img align="right" src="{{ site.url }}/images/oompa.jpg">  
<br>
<br>
<br>
<br>
<br>


Okay, not really. Things do get a little more complicated, but there are programs you can download to get the appropriate bash environment on a PC. Unfortunately I have no experience with that though so you'll have visit the almighty [google](https://www.google.com/search?source=hp&q=running+terminal+on+pc&oq=running+terminal+on+pc&gs_l=psy-ab.3..0i22i30k1l2.4571.9293.0.10104.25.22.0.0.0.0.175.1978.12j9.21.0....0...1.1.64.psy-ab..4.20.1827.0..0j0i131k1j0i10k1j33i22i29i30k1.RZTO4OlOhZk).


<br>
<center><img src="{{ site.url }}/images/under_construction.jpeg"></center>
<center><h3>UNDER CONSTRUCTION</h3></center>
<br>


### <u>Running a command at the prompt</u>
The terminal prompt by default is typically a line with some information ending with a dollar sign:
<center><img src="{{ site.url }}/images/blank_prompt.png"></center> 
<br>
This line and the terminal itself are both customizable so things may look a little different on yours, but this is where the magic happens. 

Throughout this tutorial I will be presenting blocks of code that will follow a dollar sign:


```bash
$ 
```


That dollar sign just represents your prompt, and shouldn't be typed when replicating the commands being used. 

To run a command in the terminal, you simply need to type the command after the prompt. Some commands don't require any further input and can be executed by hitting the return key. For example, there is a command called `date` that will output the time and date:
<center><img src="{{ site.url }}/images/terminal_date.png"></center> 
<br>
And as you can see, after hitting return, a new prompt line appears and the cursor is waiting for the next command. 

The `date` command didn't require anything else, so it worked just fine by itself, but some other commands require what are known as arguments. In these cases the format would be `command` `argument` (separated by a space). Some arguments are optional, and some are mandatory. Probably the most common type of argument is telling *bash* which file you want to do something to. 

For example, the `head` command prints out the first few lines of a file to the terminal, so you need to tell it which file you want it to act on. Take for example this regular text file named 'test.txt':
<center><img src="{{ site.url }}/images/test.txt.png"></center> 
<br>
If I run the 'head' command on this file, the terminal prints out the first 10 lines and then returns to the prompt:
<center><img src="{{ site.url }}/images/head_example.png"></center> 
<br>
In this case, providing a file to the 'head' command is required. But this command also has optional arguments. For example, by default it prints out only the first 10 lines, but we tell it to do however many lines we want by adding what's known as a 'flag'. These are often led by a single dash, followed by a character representing which argument you're specifying, followed by the value (or file) you want to give it. Here is an example where we tell the 'head' command we want the first 4 lines of the document, overriding the default 10:


```bash
$ head -n 4 text.txt
```



<center><img src="{{ site.url }}/images/test.txt_head4.png"></center>  

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
In that image you can see we are in a directory called "temp", and inside that directory there is a sub-directory called "another_directory" and two text files. Additionally at the very bottom there is a line that tells us 'where' we are in the computer:
<center><code>Macintosh_HD/Users/Mike_Lee/Documents/web_content/temp</code></center>
<br>
This line of directories delimited by forward slashes tells us where we are; it is just like an address. 

When we are working in a terminal we need to know where we are as well. We can get all of the same information in the terminal by using the commands `pwd` (print working directory – to view the address of the directory we are in) and `ls` (list, to list the contents of the directory we're sitting in):

<center><img src="{{ site.url }}/images/terminal_directory_example.png"></center>  
<br>

It is important to be comfortable thinking about where you are in your computer when working in the terminal. One of the most common errors/easiest mistakes to make is to be trying to do something to a file that isn't where you think it is. Let's go back to our example above where we used the `head` command:

<center><img src="{{ site.url }}/images/file_location_error.png"></center>
<br>

Need to point at other directory, then merge last two sections







special characters for changing directories '~'  '/'  '.' '..'
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

