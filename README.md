# ACGTrie
**[UPDATE: For the latest and fastest code, check out the dev branch by Jim Bailey!]**

Make exploring DNA composition easy!

Sequenced DNA has two inherent properties - it's position in the genome/transcriptome from which it came, and the composition of the DNA itself.
Whilst there are *many* tools for mapping reads to positions on the genome and deriving information there, far fewer tools exist that ignore positional information altogether and focus solely on the DNA itself.

Of the few tools that do focus soley on the DNA, almost all of them build fixed length "k-mer" lookup tables - tables which can be thought of as storing a fixed amount of DNA per row, as well as the frequency of the fragment in the input data. For example:

    ACGTA 39
    CATGA 4
    TTACG 13
    CGAGC 21
    ...

These tables make it quick and easy to look up how much "ACGTA" there is in some input data, and are extremely useful in a number of specific circumstances like de novo genome assembly. However, if you wish to look at another length of DNA, say "ACGTAG" or perhaps "ACGT", you would have to re-run the whole analysis again, often on specialist hardware with a lot of RAM. You would also end up with two output files which, although they share no overlapping information per-se, could be stored in another much more efficient single format.

But none of these issues are really deal-breakers on their own; you can always buy a bigger computer or more disk space afterall... no, what inspired the ACGTrie project and what is argubly the most restrictive part of k-mer analysis is that if Researchers have to choose a fixed-length k-mer before they can even start analysing their DNA composition, their entire downstream analysis is dominated by an aribitary choice of k-mer length - a choice which they cannot possibly know in advance how significant it's effect will be on the end result. In short, making big choices early on and not incentivising curiosity restricts the space that DNA composition analysis has to grow. Curiosities or hunches go uninvestigated, because exploring them is simply too inconvenient..

So the ACGTrie project has two primary goals: 
- To make getting the DNA composition statistics as stress-free as possible.
- To make the output format as convenient as possible for others to pick up and hack on - to find interesting and new ways to analyse, compare and contrast DNA composition information from one or more samples.

To read more on how we aim to go about doing that, check out the OUTPUT, PREPROCESSORS, and UP2BIT readmes in the docs directory :)
