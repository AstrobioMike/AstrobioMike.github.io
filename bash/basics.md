---
layout: page
categories: [bash, tutorial]
exclude_from_nav: true
---

<h2><i>bash</i> basics</h2>
Here we are going to cover the very basics of working in the terminal. For quick definitions of what things like *bash* and 'the terminal' are, see [here]({{ site.url }}/bash.html). We'll start at the very beginning with the basic formula for how to run commands, then  move onto how to navigate around your computer with these commands, and end with how to manipulate and create documents. This is probably suitable for those with little to no experience, so if you feel comfortable with these tasks already, consider jumping ahead to one of the sections you're less familiar with. 


<br>
<center><img src="{{ site.url }}/images/under_construction.jpeg"></center>
<center><h3>UNDER CONSTRUCTION</h3></center>
<br>


### <u>Running a command at the prompt</u>
The terminal prompt by default is typically a dollar sign:
<center><img src="{{ site.url }}/images/blank_prompt.png"></center> 
<br>
This line and the terminal itself are both customizable so things may look a little different on yours, but this is where the magic happens. 

Throughout this tutorial I will be presenting blocks of code that will follow a dollar sign:

{% highlight rouge %}
$
{% endhighlight %}

That dollar sign just represents your prompt, and shouldn't be typed when replicating the commands being used. 

To run a command in the terminal, you simply need to type the command after the prompt. Some commands don't require any further input and can be executed by hitting the return key. For example, there is a command called `date` that will output the time and date:
<center><img src="{{ site.url }}/images/terminal_date.png"></center> 
<br>
And as you can see, after hitting return, a new prompt line appears and the cursor is waiting for the next command. 

The `date` command didn't require anything else, so it worked just fine by itself, but some other commands require what are known as arguments. In these cases the format would be `command` `arguments` (separated by a space). Some arguments are optional, and some are mandatory. Probably the most common type of argument is which file you want to do something to. 

For example, the `head` command prints out the first 10 lines of a file to the terminal, so you need to tell it which file you want it to act on. Take for example this regular text file named 'test.txt':
<center><img src="{{ site.url }}/images/test.txt.png"></center> 
<br>
If I run the 'head' command on this file, the terminal prints out the first 10 lines and then returns to the prompt:
<center><img src="{{ site.url }}/images/head_example.png"></center> 
<br>
In this case, providing a file to the 'head' command is required. But this command also has optional arguments. For example, by default it prints out only the first 10 lines, but we tell it to do however many lines we want by adding what's known as a 'flag'. These are often led by a single dash, followed by a character representing which argument you're specifying, followed by the value (or file) you want to give it. Here is an example where we tell the 'head' command we want the first 4 lines of the document, overriding the default 10:

{% highlight rouge %}
$ head -n 4 text.txt
{% endhighlight %}

<center><img src="{{ site.url }}/images/test.txt_head4.png"></center>  

<br>
And that's really it. Those are the fundamentals that govern running virtually any individual command in the terminal. Now that you're familiar with this baseline formula, let's move on to looking at more commands and learn how to navigate around your computer from inside the terminal.  
<br>

---
<br>
### <u>System structure</u>
file tree structure, pwd, ls
<br>  


{% highlight rouge %}
$ for i in `cat samples`; do echo $i; wc -l "$i".txt; done
{% endhighlight %}

### Moving around  
cd, special characters for changing directories '~'  '/'  '.'
<br>

### Probing files
less, head, tail, wc
<br>

### Manipulating files
cp, mv, rm, mkdir, rmdir 
<br>

### Making plain text files
nano, vim, emacs


