---
layout: main
title: Going deeper with indexing
categories: [R, tutorial]
permalink: /R/more_indexing
---

{% include _R_indexing_toc.html %}

{% include _side_tab_R.html %}

Okay, so some of the awesome utility of R is how easy it makes it to subset things down to whatever you want. The fundamentals of this were covered [here](/R/basics#the-wonderful-world-of-indexing) (and if you're completely new to R be sure to run through the [basics](/R/basics) before this), but now we're going to go a bit deeper into the weeds of how exactly subsetting by [conditional statements](/R/basics#subsetting-by-conditional-statements) works. Having a better grasp of what R's doing under the hood will greatly improve your abilities to write and make use of more complicated subsetting expressions.
<br>
<br>

---
---
<br>

# Subsetting by position

Subsetting by index number (the actual sequential position of a value) is pretty straightforward as we covered [here](/R/basics#subsetting-by-position). Let's make a simple vector to work with and then look at that again: 

```R
y <- c(5, 6, 7)
y
```

<center><img src="{{ site.url }}/images/R_vector.png"></center> 
<br>
And now let's subset it by index values:

```R
y # returns the whole vector
y[1] # returns the first item
y[2] # returns the second item
y[3] # returns the third item
y[c(1,3)] # returns items 1 and 3
```

<center><img src="{{ site.url }}/images/R_vec_index.png"></center> 
<br>

Ok. So each value in a vector has an index number which is just the order of things, and we can use that index number to pull out the value that is held at that position. And we saw how this worked for 2-dimensional structures like tables [here](/R/basics#subsetting-tables).  
<br>

---
<br>
# Subsetting by conditional statements

We also saw how to subset using a conditional statement like so:

```R
y # the whole vector
y[y >= 6] # returns just the last two values
```

<center><img src="{{ site.url }}/images/R_vec_index_TF5.png"></center> 
<br>

Now we're going to break down exactly what is going on here, and it involves the use of `TRUE` and `FALSE`. `TRUE` and `FALSE` are special terms in R that are of type "logical". You can see this if you run `class(TRUE)`.  

### Conditional statements resolve to TRUE/FALSE

When subsetting by index number above we provided R a *numeric* vector within the subsetting brackets: `y[c(1,3)]`. When given a numeric vector like this, R interprets those numbers as which index positions we want to pull out.  

If we instead place a *logical* vector of `TRUE` and `FALSE` values within the subsetting brackets ( `[ ]` ), R will return only the values corresponding to index positions where our `TRUE/FALSE` vector holds `TRUE`.  

This is definitely a little abstract at first, but it's integral to how more-advanced subsetting is done in R, so it's absolutely worth fighting through the initial head-scratching phase! Let's look at this in practice:

```R
y # the whole vector
y[c(FALSE, TRUE, TRUE)] # just the last two values
```

<center><img src="{{ site.url }}/images/R_vec_index_TF1.png"></center> 
<br>
We can see here R only returned the last two values of the vector stored in variable `y`, where the `TRUE/FALSE` vector we provided within the subsetting brackets contained the values `TRUE`.  

Of course in this case we still specified exactly what we wanted by making those positions `TRUE` in our subsetting logical vector, and so at this point this isn't any more useful than giving the exact index positions like above. But we actually don't need to explicity provide the `TRUE/FALSE` vector. **Instead, we can provide a *conditional expression* that will resolve to a `TRUE/FALSE` vector.** Conditional expressions involve things like greater than `>`, less than `<`, equal to `==`, greater than or equal to `>=`, and less than or equal to `<=`. To see what this means in practice, let's look at how R handles some conditional expressions:

```R
2 > 1 # FALSE 
2 > 2 # FALSE
2 >= 2 # TRUE
```

<center><img src="{{ site.url }}/images/R_vec_index_TF2.png"></center> 
<br>
R returns "logical" values (i.e. `TRUE/FALSE` values) when presented with conditional expressions like this. As such `2 > 1` resolves to `FALSE`, `2 > 2` resolves to `FALSE`, and `2 >= 2` resolves to `TRUE`. Ok, great, let's walk this just a few more steps forward.  

If we provide a vector instead of a single digit, R will check the conditional expression on each item of the vector. Let's see this in action (remember the `c()` function from above to create a vector):

```R
c(5, 6, 7) >= 6 # FALSE, TRUE, TRUE 
```

<center><img src="{{ site.url }}/images/R_vec_index_TF3.png"></center> 
<br>
Here, R checked if each of the values within the vector (5, 6, and 7) resolved to `TRUE` or `FALSE` based on the conditional `>= 6`, and returned a *logical* vector of `FALSE`, `TRUE`, and `TRUE`. We can also provide a vector stored in a variable and R acts upon it the same way:

```R
y # our vector from above
y >= 6 # FALSE, TRUE, TRUE
```

<center><img src="{{ site.url }}/images/R_vec_index_TF4.png"></center> 
<br>
Now, also remember from above that we directly provided a `TRUE/FALSE` vector for our first example of how to subset using a "logical" vector:  

```R
y # the whole vector
y[c(FALSE, TRUE, TRUE)] # just the last two values
```

<center><img src="{{ site.url }}/images/R_vec_index_TF1.png"></center> 
<br>
So you might be able to see where this is going now, when we provided that conditional statement `y >= 6` within our subsetting brackets, R was doing this and returning only the positions where the condition resolved to `TRUE`:

```R
y # the whole vector
y[y >= 6] # returns just the last two values
```

<center><img src="{{ site.url }}/images/R_vec_index_TF5.png"></center> 
<br>
The way I read the expression `y[y >= 6]` in my head is: "Give me all the values of vector 'y', where 'y' is greater than or equal to 6." This fundamental concept is key to what makes indexing in R so powerful!  

One last thing we're going to introduce here is the `!` character, which inverts the interpretation of `TRUE` and `FALSE`. As we've seen, `y[y >= 6]` will return all values within 'y' that are greater than or equal to 6. But if we add in the `!` point, it will return the opposite for us. We can see this just by putting this in front of the portion of the expression that generates the T/F vector:

```R
y

y >= 6 # returns FALSE, TRUE, TRUE
!y >= 6 # returns TRUE, FALSE, FALSE

y[y >= 6] # returns only 6 and 7
y[!y >= 6] # returns only 5
```

<center><img src="{{ site.url }}/images/R_vec_index_TF6_not.png"></center> 
<br>
The use of the `!` character in this case may seem a little unnecessary when we could just switch around our equality expression, but it's very handy for other types of conditional statements like the one we ran through [here](/R/basics#cond_example), where inverting the `TRUE/FALSE` vector was the only way we could pull out what we wanted. 
<br>

---
---
<br>
If you're still reading this then kudos to you for getting through that! This concept of indexing with logical vectors underlies so many of the ways that we parse down data in R, and it's going to help immensely if you understand the fundamentals of it.
