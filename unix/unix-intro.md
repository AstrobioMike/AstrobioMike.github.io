---
layout: main
title: Unix crash course
categories: [unix]
tags: [unix,bash,bioinformatics,tutorial]
permalink: /unix/unix-intro
---

{% include _side_tab_unix.html %}

**Unix is very likely the most foundational skillset we can develop for bioinformatics** (and much more than bioinformatics). Many of the most common and powerful bioinformatics approaches happen in this text-based environment, and having a solid foundation here can make everything we’re trying to learn and do much easier. This is a set of 5 introductory tutorials to help us get from being completely new to Unix up to being great friends with it 🙂 

<hr style="height:10px; visibility:hidden;" />

---
<br>

# What is Unix?
The term "Unix" encompasses a family of operating systems that share a common ancestor that was initially developed in the late 1960s. It was built based off a modular design philosophy, where relatively simple tools can be strung together seemlessly in order to accomplish more complicated tasks. In short, it has spread like crazy throughout the computational world over the past 5 decades. 

Many successful offshoots of Unix have been developed, which is why we might see the term "Unix-like", in order to encompass all the systems that aren't technically "Unix". 

Here are some terms that are often used interchangeably – not because it's important to remember them or any differences (it's not for must of us), but just to have them laid out somewhere in case you here people saying them.

| Term     | What it is          |
|:-------------:|------------------|
| **`shell`** | what we use to talk to the computer; anything where we are pointing and clicking with a mouse is a **G**raphical **U**ser **I**nterface (**GUI**) shell; something with text only is a **C**ommand **L**ine **I**nterface (**CLI**) shell |  
| **`command line/terminal`** | a text-based environment capable of taking input and providing output |  
| **`bash`** | the most common programming language used at a Unix command-line |  
| **`Unix`** | a family of operating systems (we also use the term "Unix-like" because one of the most popular operating systems derived from Unix is specifically named as [*not* being Unix](https://en.wikipedia.org/wiki/GNU){:target="_blank"}) |  


<hr style="height:10px; visibility:hidden;" />

---
<br>

# Why learn Unix?
Getting familiar with working at a "Unix-like command-line” is one of the most fundamental skillsets we can develop for bioinformatics, but also much, much more. As Brian Kerrigan (a team member of the original Unix team) puts it in his 2019 book [*Unix: A history and a memoir*](https://www.cs.princeton.edu/~bwk/memoir.html){:target="_blank"}:

> "*Unix and its derivatives aren't widely known outside a particular technical community, but they are at the heart of any number of systems that are part of everyone's world. Google, Facebook, Amazon, and plenty of other services are powered by Unix-like operating systems. If you have a cell phone or a Mac, it runs on some version of Unix. If you have gadgets like Alexa at home or navigation software in your car, they're powered by Unix-like systems too.*"

Being the framework for so much of our world, learning to speak its language also gives us greater access to things like remote servers and cloud-computing. It can allow us to access and manipulate large datasets we otherwise couldn't, and enable us to use programs we otherwise couldn't.

Which brings us back to it being foundational to bioinformatics. Many of the most common and powerful bioinformatics approaches happen in this text-based environment, and having a solid foundation here can make everything we’re trying to learn and do *much* easier. 

**Sooo, here are just a few reasons why it's worth it to learn Unix:**  

* it’s the foundation for most of bioinformatics (and much more)  
* enables the use of non-GUI (Graphical User Interface) tools  
* improves reproducibility (GUI's are super-convenient for lots of things, but they are not ideal when it comes to reproducibility)  
* enables things like quickly performing operations on large files (without needing to read them into memory)  
* can allow us to programmatically access data  
* helps automate repetitive tasks (need to rename 1,000 files?)  
* enables use of higher-powered computers elsewhere (servers/cloud-computing)  


<hr style="height:10px; visibility:hidden;" />

---
<br>

# So let's get to it!

>**NOTE**  
> Keep in mind while going through these pages is that this is all about <i>exposure</i>, not memorization or mastering anything. Don't worry about the details!  

1. [Getting started](/unix/getting-started)
2. [Working with files and directories](/unix/working-with-files-and-dirs)
3. [Redirectors and wildcards](/unix/wild-redirectors)
4. [Six glorious commands](/unix/six-glorious-commands)
5. [Variables and For loops](/unix/for-loops)  

<hr style="height:10px; visibility:hidden;" />

---
---

<h5><a href="/unix/getting-started" style="float: right"><b>Next:</b> 1. Getting started</a></h5>
