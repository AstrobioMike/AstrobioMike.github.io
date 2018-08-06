---
layout: main
title: Unix
categories: [bash]
permalink: /bash/
---

{% include _side_tab_bash.html %}

**Shells and Unix and terminals, oh my!**  
Along with my disclaimers about not being an expert as laid out [here](/about/#Disclaimers), I'm not going to try to give you technical explanations of all the underlying infrastructure of how working in a terminal environment works. I simply don't have a strong enough understanding myself. I've read about all these things in the many O'Reilly books I have; I get curious from time to time and google things and remember stuff (at least what I understand) for a short while afterwards; and I occasionally get to corner someone who really knows their stuff (like [Evan Boylen](https://twitter.com/ebolyen)) and ask them a bunch of questions. But I don't seem to retain many of the details, and while we all want to know everything, we also all need to get the current analysis done so we can get to the rest of our never-ending list of things we're behind on.  
<br>
So if you'd like to get started right away, head on over to [Unix basics](/bash/basics) and get started now. You can always revisit these concepts when your curiousity outweighs your motivation to do work. I've also found that (at least for me personally) it's a big help to be able to throttle back the instinct to want to understand everything *right now*. And that it's totally okay (and often better) to just let a deeper understanding grow as we develop our skills in this sort of expanding sphere into all of these different elements. But if you want to read over a few terms you'll hear and see a lot, here I'll give just the layman understanding that I have – with the unspoken agreement between us that, if you were so inclined, just starting at the wiki for any of these things can take you much further.  
<br>

---  
<br>

<h3>Some terminology</h3>

**Unix**  
Unix is a type of operating system, or more precisely, a family of operating systems as the [wikipedia](https://en.wikipedia.org/wiki/Unix) words it. This is why you will hear people say things like "Unix-like", and refer to Unix being on both Mac and Linux computers (they run Unix-like operating systems), but not on PCs (they do not 😞 ). 

**Terminal**  
A terminal is just a text-based, command-line environment where you give input and get output. It's a place where you can talk to your operating system through what is known as a *shell*. You indirectly interact with your operating system all the time, like when you open the internet, or are editing figures, etc., but through a terminal you have much more freedom. For instance, in a typical graphical-user-interface program like we use all the time, you can only ask the computer to do explicitly pre-defined tasks. Whereas in the terminal, you can ask much more personal questions and have much more interesting conversations with your operating system. Seriously though, the real value of this is that it results in your having more capabilities to do what *you* need to do, because you're no longer limited to doing only things that have been precisely laid out before. 

**Shell**  
This is what runs in a terminal. A shell is your ambassador to your Unix-like operating system. It allows you to interface with the system you are working on by taking what you want, translating it into something the operating system understands, and then translating the output back into something you can understand. There are multiple types of shells (*bash*, *sh*, *ksh*, *csh*, bla bla bla, and lots more), but by far the most common these days is *bash*. 


***Bash***  
*Bash* is currently the default shell in most environments you will work in, which is great. There are, however, multiple version of *bash*, which means if you regularly work on different systems (say your own computer as well as on various servers), even though you will likely be working in *bash* on all of them, you will still occasionally run into minor differences in syntax that you'll have to navigate from time to time. Fortunately these are usually easily solved with minimal googlation effort as everyone else has hit those exact same snags before. 
<br>  

---

<br>
So there you have it. In 4 paragraphs you now know as much about these terms as I know after 4 years of being in the weeds with them. Hmm, saying it that way I'm not sure if you should feel better or if I should feel worse... 🤔  

Either way, head on over to [*bash* basics](/bash/basics) to get started!
