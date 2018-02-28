---
layout: main
title: Amplicon analysis
categories: [amplicon]
permalink: /amplicon/
---

{% include _side_tab_amplicon.html %}

<div class="warning">
<h2>WARNING!</h2>
Don't get too lost in the weeds when working with marker-gene data. More often than not the kind of detailed questions that will land you in that situation can't be answered with this type of data anyway!</div>

**Amplicons and marker-genes and tags, oh my!**  
Most often a marker-gene analysis is the microbial ecologist's first tool in a vast toolkit. It is primarily used as a broad survey of community structure. As the warning notes above, when first beginning to work with this type of data it can be easy to get caught spinning your wheels about a subtle component in your processing pipeline that ultimately has a negligible impact compared to the noise we are working through. What I mean by this is, generally speaking, tag sequencing is not the appropriate tool to answer really meticulous questions. It is a tool for comparing baseline *proxies* of metrics about microbial communities. It is a tool of exploration and hypothesis generation, not hypothesis confirmation.  

There are a lot of things to keep in mind regarding what tag sequencing means, what it doesn't mean, what it can tell you, and what it can't. And people more immersed in the intricacies than I have written and published lots about this already. But at some point I'll give my two cents on all of these things in the [Caveat Central]({{ site.url}}/amplicon/caveats) post for those interested. For now, if you're new to this type of analysis, below you can find some common terminology defined as I understand things, and you can walkthrough a complete [example workflow]({{ site.url }}/amplicon/workflow_ex) for processing and analyzing a tag dataset.  
<br>  

---  
<br>

<h3>Some terminology</h3>
As with most things, trying to pin down precise, one-answer definitions can be difficult sometimes. This is because language is fluid and we are constantly trying to find better ways to communicate more efficiently and with more clarity. That's a good thing, but it can add to the confusion when first trying to find your footing in a new area. Here are some terms you might hear and what they mean most often in my experience (though there can certainly be other interpretations).

**Amplicon**  
The resulting sequence of a targeted amplification of genetic material. Targeted meaning primers were used 

**Marker gene**  
A gene that can be useful for delineating organisms, like a fingerprint.

**OTU**  
Operational Taxonomic Unit. OTUs are an artificial, arbitrary construct useful for grouping sequences together into units to help us summarize and analyze things, and they also intended to help deal with  sequencing error intrinsic to the technology.

**ASV**  
Amplicon Sequence Variant. Resulting sequences from newer processing methodologies that attempt to take into account sequencing error rates and are believed to represent true biological sequences. The case for moving towards using ASVs over OTUs is made very well in [this paper](https://www.nature.com/articles/ismej2017119) by [@bejcal](https://twitter.com/bejcal), [@joey711](https://twitter.com/joey711), and [@SherlockpHolmes](https://twitter.com/SherlockpHolmes).

**Tag sequencing**  
In my opinion, the most useful and most often used meaning of "tag sequencing" is when you are using primers that target something specific, like a marker gene. So the opposite of shotgun sequencing which uses random hexamer or nonamer primers. But on rare occassions I have heard this used to mean barcoding samples.

**Barcodes**  
Barcodes refer to the sequences ligated to your individual samples' genetic material before they get all mixed together to be sequenced together. These barcodes are then unique to each sample, so you can afterwards identify which sequences came from which samples.

**Demultiplex**  
Demultiplexing refers to the step in processing amplicon sequence data where you'd use the barcode information in order to know which sequences came from which samples after they had all been sequenced together.
