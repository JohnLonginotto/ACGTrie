# Preprocessors

ACGTrie is designed to be as general-purpose and generic as possible, always taking its input data from the standard unix input 
pipe (**stdin**). DNA sequences however are generally kept locked up in a variety of different specalized  Bioinformatic formats - 
**BAM** files, **SAM** files, **FASTA** files, **SQL** databases, etc - and so to access that data and send it to ACGTrie one has to find a 
preprocessor specific to their sequencing data format. 

Different preprocessors can be found in the *preprocessor* directory, so if you just want to build an ACGTrie as quickly as 
possible, go there, grab a pre-processor (and a copy of ACGtrie, presumably the *FAST* version), and start making DNA 
composition tables :)

If however you cannot find a pre-processor for your DNA input format, or you want to write a custom preprocessor for a specific 
analysis (say, just reads from Chromosome 1 between positions x and y, or perhaps a pre-processor that reverse-compliments DNA 
on the reverse strand before sending it to ACGTrie), then you may have to get your hands dirty and write one yourself. 
Don't worry, a simple pre-processor can be written in just 5 lines of code!
There can be, however, a lot more to pre-processors than just getting DNA and piping it to ACGTrie. 
A smart pre-processor can significantly speed up the time taken to make the output. 

The rest of this document focuses on the known optimizations that preprocessors can do to speed things up, which is something 
you will want to read about if you intend on making your own pre-processors. They broadly fall into the 7 catagories:

- run length encoding,
- fragmentation,
- sorting
- parallelization & memory reduction,
- clipping, and
- filtering/deleting

# Run Length Encoding

Often in genomic datasets, DNA is sorted by position in the genome. This means that very often we will be feeding into ACGTrie the
exact same fragments multiple times over. For example, if the *genome* was:

                                               ACGTTTAACGTTGCA

We might expect our sorted sequences to be:

                  ACGT ACGT ACGTT CGT CGTTT CGTTT TTTAA TAACG AACGT ACGTT ACGTT TTGC TGCA

Note the repitition of some adjacent DNA sequences such as the first two reads. Now we COULD pipe these reads directly into ACGTrie,
and it would to do **13** insertions in total. However, we could also keep a little running buffer in the preprocessor which stores
the fragments and their frequency, and then pass the fragments and their counts to ACGTrie in CSV format. For example, if our
buffer stored 2 fragments at a time:

                  ACGT,2  ACGTT,1  CGT,1  CGTTT,2  TTTAA,1  TTTGC,1  AACGT,1  ACGTT,2  TTGC,1  TGCA,1

Now ACGTrie only has to insert **10** items, because as it updates the COUNT value for each fragment, instead of adding 1 it simply adds 2.

In fact, if we processed ALL the data first and then added it to the trie, we would see that *ACGTT* actually appears *3* times, and
we could get it down to just **9** insertions - but of course if we could do that for a whole sequencing file, we wouldn't need ACGTrie
in the first place because we have essentially created a hash table :) But since hash tables will not fit into memory, there is a 
trade off for how big the buffer should be, and how often we send buffered data to ACGTrie. 
This practice is commonly known as Run Length Encoding or RLE for short, and it is the first optimization your pre-processor can do.
You do not have to do anything special to get RLE for ACGTrie other than send it comma-delimited, newline-sepurated text rather than
just standard newline sepurated text.

# Fragmentation

While it might be nice to make an ACGTrie of just input sequences/reads, most of the time we actually are interested in the DNA composition. To go from DNA fragments/reads to DNA composition, all we need to do is fragment our read into all the possible sub-fragments (note: I know the terminiology is confusing because a read is already a fragment of the genome, so we are "fragmenting fragments to get sub-fragments"!). By way of example, if we had the DNA 'ACGT', we would fragment it like so:

                                   'A', 'C', 'G', 'T', 'AC', 'CG', 'GT', 'ACG', 'CGT', 'ACGT'

This can either be done by either the pre-processor, or by ACGTrie it itself if you use the --fragment parameter. While the latter is certainly simpler, and allows for situations such as piping data straight out of a file or database, it is actually often beneficial to have the pre-processor do it for a number of reasons. The first is that we saw previously that keeping a buffer of fragments can reduce the amount of repetition ACGTrie has to do. Well keeping a buffer of sub-fragments is *even more* beneficial, because subfragments are far less likely to be unqiue compared to fragments, and the more we can pile up before sending to ACGTrie, the less work ACGTrie has to do!

But there is another, even better reason. In the above example, ACGTrie would have to find the A row, then +1. Then later find the AC row and +1. Then the ACG row and +1. And finally the ACGT row a +1. That requires 4 sepurate trips to the trie to complete. A **much** more efficient method would be to fragment the data in a special way:

                                                   'ACGT', 'CGT', 'GT', 'T'

and then tell ACGTrie to use the special --walk parameter. This means that for every row in the table that ACGTrie visits to get to its destination row, it +1s along the way. This means instead of 10 trips to the trie, we only need 4. The downside of this is that the pre-processor has to fragment in this very specific way. Fortunately, it is not that complicated - all you do is take a bite off the left of the string 1 character at a time until there is nothing left.
The best thing about this method is that it also works seemlessly with Run Length Encoding, as ACGTrie can +X to every row as it walks the trie.

# Sorting

Unfortunately, ACGTrie has 1 property that concerns its developers and that is that the order in which sequences are added to the trie effect both the time taken to add data, and the final structure of the trie. Put simply, if we add the longest subfragments first, those subfragments will have to be broken up into two or more rows when smaller, partially matching subfragments are added later. Alternatively, adding the shortest subfragments first means no splitting of subfragments is ever needed. You may intuitively think that the former is slower than the latter, but actually you would be mistaken - it is often faster to add the long subfragments first, then split them later if you have to, than it is to always have to hop through many rows to add the long fragment at the end.
So why not always add longer fragments first? Well unfortunately there is a downside to that approche too, and that is that spliting rows up a lot leads to a messy and sometimes less efficient trie than if rows are only ever added when needed. The result is that longest-first tries tend to be larger in their final size than shortest-first, but are made quicker.
In the future will will probably have ACGTrie do some sort of normalization on the trie so that the input order does not effect the result... but that will have for a rainy day :) For now, just know that if you change the order of the data you feed to ACGTrie, you likely change its resultant file size, and will almost definitely change its MD5 checksum!

# Parallelization & Memory Reduction

It is said that you can't build a trie in a parallelized way, because two processes might interfear with each other when they try to read/write rows to the table - a situation known as a Race Condition. Fortunately, this is where our pre-processor can step in and split the workload up into *branches* for multiple processors.

Again, it is really quite simple. After creating a buffer of subfragments, if we only select subfragments which start with, for example, an *A*, we will only create 1/4th of the total trie. Row 0 will only have 1 warp pipe - to an A. If our preprocessor was to only select subfragments starting with A **and** before feeding those subfragments to ACGTrie cut off that initial A, the resultant trei would have a root row which would be essentially identical to the first A row of the full trie.

So you can probably see where this is going... if we have the pre-processor create subfragments in its buffer as always, but pipe certain subfragments to different ACGTrie processes (depending on what they start with), making sure to cut off whatever DNA it is that defines which ACGTrie process the subfragment should go to, we can create multiple tries in parallel and simply concatinate them all back up together at the end by creating a root node that points to the first line of each file segment.
In fact, there is no reason to only split the data 4-ways. We could split it based on the first 2 bases (16-ways) or 3 bases (81-ways), etc. At which level to split the branches depends on how many processors you want to use simultaniously, or how much total memory you have to commit to ACGTrie. In theory, extremly large tries could be built with very limited total RAM if the input data was split in this way and read many times over.

But again, this only works if *the preprocessor is doing the fragmenting*, or if *no fragmentation at all* is desired.

## Clipping

Depending on the questions and input data, 99% of the information that can be gained from looking at DNA composition can be seen in the first 20 or 30 hops away from the root node. Beyond that you will see extreamly long, but not very high frequency, sequences. If output filesize is important to you, one option is to have the pre-processor create subfragments in the usual way, but only keep the first X bases of the subfragment in the buffer. This will reduce the total output size considerably, but
should still give you more than enough information to get an idea of the overall DNA composition. It is very to remember to clip your sequences *after fragmentation*, and not clip before fragmentation, as this will change everything!

##  Filtering / Deleting

Again, this is just another technique to reduce the amount of data added to the trie, or prune the trie after creation.
In our initial tests, more than half of a final trie contains rows with less than 5 COUNTs in total. Unlike the clipping technique however, while these rows may predominantly be at the terminating (no warp pipe) rows of the trie, it is not guarenteed, so you end up with a tree in which you dont really know what has been deleted. For some statistics you may with to compare a rows count to the total count of that row's depth, and by deleting indiscriminatly rows with a COUNT of less than X, it may effect this statistics. There is also the option of deleting rows in the preprocessor before they even make it to ACGTrie based on probabalstic models much like the Bloom filters of typical k-mer analysis tools, but for now that level of complexity and uncertainty does not seem to be worth the cost-benefit of simpler, although perhaps slower, methods of pruning the trie.
