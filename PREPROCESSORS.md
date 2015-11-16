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
- parrallelization,
- maximum-memory reduction,
- clipping,
- filtering, and
- sorting

# Run Length Encoding

Often in genomic datasets, DNA is sorted by position in the genome. This means that very often we will be feeding into ACGTrie the
exact same fragments multiple times over. For example, if the *genome* was:

                                               ACGTTTAACGTTGCA

We might expect our sorted sequences to be:

                  ACGT ACGT ACGTT CGT CGTTT CGTTT TTTAA TTTGC AACGT ACGTT ACGTT TTGC TGCA

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

So, say we read some DNA such as 'ACGT' from a sequencing file. At some point, we are going to have to break it down into 
all the possible sub-fragments that it could be. In full, this would look like:

                                   'A', 'C', 'G', 'T', 'AC', 'CG', 'GT', 'ACG', 'CGT', 'ACGT' 

We call this fragmenting, and it can either be done by either the pre-processor, or by ACGTrie it itself if you use the --fragment 
parameter. While the latter is certainly simpler, and allows for situations such as piping data straight out of a file or 
database, it is actually often beneficial to have the pre-processor do it for a number of reasons. The first is that keeping a 
buffer of sub-fragments in the pre-processor is far superior to keeping a buffer of unique fragments, because even if two consecutive
reads are not 100% identical, their subfragments might be - and since ACGTrie adds subfragments and not actually whole fragments,
buffering subfragments is ideal.
