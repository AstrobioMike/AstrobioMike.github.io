---
layout: main
title: The wonderful world of loops
categories: [bash, tutorial]
permalink: /bash/loops_to_help
---

{% include _bash_loops_toc.html %}

{% include _side_tab_bash.html %}

Loops are extremely powerful in all programming languages. They are what let us write out a command or operation once, and have it run on all of our samples or files or whatever we want to act on. Not only is this powerful, but it also helps with keeping our code more concise and readable, and it helps elmininate some more of our mortal enemy (human error). Again, this is because we only need to write out what we want to do once and it will then be done the same way to all. There are multiple types of loops, but here we are going to cover what is probably the most common type: the for loop. If you are new to the command line, this will be much easier to follow if you've run through the [intro to bash](/bash/basics){:target="_blank"} first. 

To get into loops, we quickly first need to more explicitly introduce bash variables.  
<br>

# Variables in bash
We already saw the special bash variable "PATH" in the [modifying your PATH page](/bash/modifying_your_path#the-path-demystified){:target="_blank"}. This is a special variable in all Unix-like systems that holds all of the locations (directories) that your computer looks for programs. When we wanted to access our PATH, we had to put a dollar sign in front of it. Variables in bash are preceded with a `$` so it knows to interpret the following string of attached text as a variable, rather than just as plain text. For example, in your terminal window, try running `echo PATH`. Now try running `echo $PATH`. With the first bash just prints out "PATH" to the screen, but with the second it interprets it as a variable and prints out your PATH (a colon-delimited list of all the directories your computer scans for programs). If you'd like you can [pipe](/bash/basics#pipes-and-redirectors){:target="_blank"} that into the `tr` command to make it more human-readable: `echo $PATH | tr ":" "\n"`.

To set a variable in bash we provide the variable name we want, an equals sign, and then what we want the variable to hold (with no spaces in between any of that). Let's try it:

```
my_var="Europa"
```

Nothing prints out when a variable is set, but the value "Europa" has been stored in the variable "my_var". We can see this if we try to do something with it:

```
echo $my_var
```

If we wanted to set a variable that held spaces, we could surround it in quotations to tell bash it should be considered as one thing:

```
my_new_var="Europa is awesome."

echo $my_new_var
```

It's important to keep in mind that these variables are just holding text. This is why we've been checking them with `echo`. If we just enter `$my_var` at the command line, when it is holding the text "Europa", we would get an error telling us "Europa: command not found". This is because to bash, us entering `$my_var` like this is equivalent to us entering the text "Europa". And when that was the only thing we provided on the command line, it was expected to be a command. We can also see this if we try to run the `cat` command on the variable:

```
cat $my_var
```

This gives us the error "Europa: No such file or directory", because there it is serving as the argument being passed to the `cat` command. If a file existed with the name "Europa" in our current working directory, then that command wouldn't give an error. Let's make a quick file and see how that would work:

```
echo $my_new_var > $my_var

cat $my_var
```

Here we `echo` the text stored in `$my_new_var` and [redirect](/bash/basics#pipes-and-redirectors){:target="_blank"} that output into a new file named with the text stored in `$my_var`. This made a file called "Europa" that holds the text within it "Europa is awesome."  

This may seem a little confusing at first, but having a good grasp of variables will help a lot when starting to work with loops. So it's worth taking a second to make sure this makes sense. Make some of your own variables and try using them with some commands.  

<br>
# For loops
We're going to be making some mock files in the following examples, so if you'd like, consider making the following directory for us to work in:

```
mkdir ~/for_loop_temp
cd ~/for_loop_temp
```

<br>
## The 4 magic words
There are 4 special words in the syntax of a bash for loop: `for`, `in`, `do`, and `done`. 

|magic word|purpose|
|:-----:|:------|
|`for`|set the loop variable name|
|`in`|specify whatever it is we are looping over|
|`do`|tell bash what we want to do with each item|
|`done`|tell bash we are done telling it what to do with each item|

Let's see what this looks like in practice. Here we are going to: name the variable "word"; loop over 3 words (cat, dog, and ukulele); and we're going to just `echo` each word, which will just print each item to the terminal. 

```
for word in cat dog ukulele
do
  echo $word
done
```

Just to note, we don't need to put these on separate lines, and we don't need to indent over the "body" of the loop like we did here, but both can help with readability. We could also enter it like so on one line:

```
for word in cat dog ukulele; do echo $word; done
```

We can also do multiple things within the body of the loop. Here we'll start with the same loop, but add another line that also sends the word into a file we'll call "words.txt":
```
for word in cat dog ukulele
do
  echo $word
  echo $word >> words.txt
done
```

Now we have a file that holds these words. Note that we used `>>` as the [redirector](/bash/basics#pipes-and-redirectors){:target="_blank"}, and not just `>`. Remember `>>` will append to the file whereas `>` will overwrite it. If we used `>`, each time through the loop we would be overwriting it and at the end we would have a file with one line of "ukulele" only.

```
cat words.txt
```

This is most often the way I specify the things I want to loop over, when they are stored in a file. As we'll see down below, I will typically make a file that holds all of the sample or genome names I'm working with for a given project. That way anytime I want to do something to all of them (which is most of the time), it's easy to do. Let's look at looping through the lines of a file.

<br>
## Looping through lines of a file
Instead of typing out the elements we want to loop over, we can execute a command in such a way that the output of that command becomes the list of things we are looping over. The syntax of how to do this may seem a little odd at first, but here is an example with our "words.txt" we just made:

```
for word in $(cat words.txt)
do
  echo $word
done
```

Here, when we say `$(cat words.txt)`, bash will perform that operation first, and then put the output in its place. In the context of items of the loop, this will make an item out of each string of characters separated by a space or a newline character. So when you have a file with just one column like we have here, each line becomes an item in the loop having the same effect as when we typed them all out like above.

#### BONUS: Brace expansion 101

We're going to use another bash trick called "brace expansion" in order to make a bunch of example files. One way brace expansion is useful is for expanding a range of numbers or letters separated by two periods (`..`). This works like this:

```
echo {1..5}
 # outputs 1 2 3 4 5

echo {a..e}
 # outputs a b c d e
```

And we can add text before or after the braces, like so:

```
echo sample_{1..5}
 # outputs sample_1 sample_2 sample_3 sample_4 sample_5

echo sample_{A..E}
 # outputs sample_A sample_B sample_C sample_D sample_E 
```

So we're going to use this to make 100 (blank) sample fasta files using the `touch` command. The `touch` command either creates an empty file if the file doesn't exist yet, or if it does exist the command updates the last modified timestamp):

```
  # currently only "words.txt" is in our working directory
ls 

  # creating 100 blank files with the common fasta extension ".fa" 
touch sample_{1..100}.fa

ls
ls *.fa | wc -l
```

As mentioned above, the first thing I typically do with a new project is make a file that holds all of the sample names, here is one way to do that:

```
ls *.fa | cut -f1 -d "." > samples
head samples
```

Right now that's sorted in computer speak, but here's one way we can change it into a more human-friendly sort:

```
sort -t "_" -nk 2 samples > s_temp
mv s_temp samples
```

(You may want to review the [bash intro](/bash/basics){:target="_blank"} and [six glorious commands](/bash/six_commands){:target="_blank"} pages if what we're doing here is unclear.)

Now the file "samples" holds all of our sample names on individual lines, in an order more intuitive for us:

```
head samples
```

We can now use this file to loop through all of our samples and do something to each of them. Notice that we just have the base sample names in there (with no extention), this means we need to provide that in the loop if we want to act on those fasta files, but it is also easier to do more. For example, say we wanted to make individual directories for each sample, and then move that sample's corresponding fasta file into it:

```
for sample in $(cat samples)
do
  mkdir $sample
  mv "$sample".fa $sample
done
```

Now we have 100 subdirectories, one for each sample, and within each is that sample's fasta file:

```
ls

ls sample_1/
```

<br>
## Nested loops
A nested loop is a loop within a loop (cue Inception theme music... yea, I don't remember it either actually). Nested loops are helpful when you want to iterate over two groups of things at the same time. This example is going to be a little contrived because we are working with blank files here, but imagine for every sample we have, we wanted to add a sequence to that sample's fasta file for each of the items in our "words.txt" file: 

```
for sample in $(cat samples)
do 
  echo "On $sample" # tell us which sample we're on
  
  for word in $(cat words.txt)
  do
    echo $word # to tell us which word we're on 
    echo ">$word" >> "$sample"/"$sample".fa
    echo "ATGCATGC" >> "$sample"/"$sample".fa
  done
    
  echo "" # adding a blank space between each sample in the output
    
done
```

Now if we look at any of our sample fasta files, they each have a sequence added for each of our words (the DNA from the ukulele might be contamination).

```
cat sample_1/sample_1.fa
```

<br>
## Setting variables within a loop
To do some more complicated operations, it is sometimes useful to set and use additional variables inside the loop. For this example, let's say we wanted to count how many total sequences were in each sample and then write this information out to a tab-delimited table that had two columns: our sample name; and how many total sequences. Each sample has the same number of sequences right now, so let's add a few sequences to a couple samples just so we that see the output actually reflects a difference:

for sample in sample_2 sample_5 sample_7 sample_99
do 
  echo "On $sample" # tell us which sample we're on
  
  for word in $(cat words.txt)
  do
    echo "  $word" # to tell us which word we're on 
    echo ">$word" >> "$sample"/"$sample".fa
    echo "ATGCATGC" >> "$sample"/"$sample".fa
  done
    
  echo "" # adding a blank space between each sample in the output
    
done

Now if we look at one of those samples, we see they have 6 sequences:

```
cat sample_2/sample_2.fa
```

So now let's count how many sequences are in each, save that as a variable that is changing on each iteration of the loop, and write out the sample and sequence count as a table in a new file:

```
for sample in $(cat samples)
do 
  echo "On $sample" # tell us which sample we're on
  
    # setting new variable that is the count of sequences in each fasta file
  num_seqs=$(grep -c ">" "$sample"/"$sample".fa)
  echo "  $sample" has "$num_seqs" sequences # print out the count of sequences for each as we go for the heck of it
  
  printf "$sample\t$num_seqs\n" >> seqs_per_sample.tsv
      
  echo "" # adding a blank space between each sample in the output
  
done
```

Now we have a new tab-delimited table that tells us how many sequences are in each sample:

```
head seqs_per_sample.tsv
```
<br>

___
<br>
Like the other pages here, this is just barely scratching the surface. BUT even though loops can get much more complicated as needed, getting these foundational skills down is all that's needed to start harnessing their awesome power 🙂
