#!/bin/sh
true = True
## The polyglot below just runs the code in pypy if it is avalible, else python.
if true :'''';then if [ $(which pypy) ];then pypy $0 $@;else python $0 $@;fi;fi;exit;'''

import csv
import sys
import json
import cffi
import time
import ctypes
import argparse
import tempfile
import datetime
import itertools
import collections

#from numba import jit

# We only use NumPy if this code is executed with pypy, otherwise NumPy is so slow 
# that genetic drift of the sequenced organism makes most results irrelevent...
usingNumpy = False
if "__pypy__" in sys.builtin_module_names:
    try:
        import numpy; usingNumpy = True; print "   [ Using NumPyPy! ]"
    except:                              print "   [ Using PyPy :) ]"
else:                                    print "   [ Using Python ]"

## FORCE USING NUMPY FOR THIS CYTHON TEST:
# usingNumpy = True
# import numpy
#cimport numpy

## Parse user-supplied command line options:
parser = argparse.ArgumentParser(     description="Put in a BAM file or fragments of DNA, get out an ACGTrie table.")
parser.add_argument("-o", "--output", metavar='/path/to/output.trie', help='Required. Output filename.')
parser.add_argument("--rows",         default=10000000, type=int,     help="Required. How big to make the ACGTrie.")
parser.add_argument("--walk",         action='store_true',            help="MANUAL: Tells ACGTrie to incrementally add to trie.")
parser.add_argument("--fragment",     action='store_true',            help="MANUAL: Tells ACGTrie to make subfragments itself.")
args = parser.parse_args()

if args.output == None: print '''
    You need to provide a path/filename for your output, as well as
    an estimate for the number of rows the final trie will contain!
    Actually, your pre-processor should have done all this for you,
    are you using one? Pre-processors are designed specifically for
    popular sequencing data formats so that DNA can pass quickly in
    to one or more ACGTrie processes in an optimized way. Visit the
    ac.gt website to find a pre-processor for your data format (BAM
    or FASTA or whatever), or, if we haven't made one for your data
    format yet, write one and we will put it up on the site for all :)
'''; exit()


#############
## ACGTrie ##
##########################################################################################################
##                                                                                                      ##
## Whats up bio-bro? Oh, you're looking to count the DNA sequences in your read data? Well you came to  ##
## the right place! Probably. I don't know. I'm not a very good salesmen...                             ##
##                                                                                                      ##
## So it's pretty simple really; if we had a sequence like 'ACGACCC', and we want to ask some questions ##
## about the "composition" of the DNA, we might be interested in a table like the one below:            ##
'''                                                                                                     ##
##                                   row    Seq      Count                                              ##
##                                    0     C        4                                                  ##
##                                    1     A        2                                                  ##
##                                    2     G        1                                                  ##
##                                    3     AC       2                                                  ##
##                                    4     CC       2                                                  ##
##                                    5     CG       1                                                  ##
##                                    6     GA       1                                                  ##
##                                    7     ACC      1                                                  ##
##                                    8     ACG      1                                                  ##
##                                    9     CCC      1                                                  ##
##                                    10    CGA      1                                                  ##
##                                    11    GAC      1                                                  ##
##                                    12    CGAC     1                                                  ##
##                                    13    ACCC     1                                                  ##
##                                    14    ACGA     1                                                  ##
##                                    63    GACC     1                                                  ##
##                                    32    CGACC    1                                                  ##
##                                    17    ACGAC    1                                                  ##
##                                    18    GACCC    1                                                  ##
##                                    19    CGACCC   1                                                  ##
##                                    20    ACGACC   1                                                  ##
##                                    21    ACGACCC  1                                                  ##
'''                                                                                                     ##
## Then we can ask questions like 'how prevelent is motif X?' or 'how much seq data covers repeats?',   ##
## or, if we have two of these tables, 'what are the key differences between sample A and sample B?'.   ##
## But at the end of the day it's really nothing more than a table of counts for all the different      ##
## subfragments of DNA within our sample. Note that 'CC' has a count of 2 due to the 'CCC', but you get ##
## the basic idea. Problem is, if we tried to make the table above for some 5Gb+ seq data, there would  ##
## be SO MANY combinations of DNA that we would run out of memory/RAM before we processed even 1% of    ##
## chromosome 1. If we tried to build the table directly onto the hard drive, it would take a VERY long ##
## time because compared to RAM, disk access is really slow. Bioinformaticians have tried to solve this ##
## problem by *approximating* the DNA composition (instead of getting exact results) which allows for a ##
## lot less RAM to be used during counting. However, these very advanced algorithums used to do this    ##
## approximation (Google "k-mer analysis") come with some really frustrating consequences:              ##
## 1) k-mer counting tools are some of the most complex tools in the whole field of Bioinformatics,     ##
## well above anything a Biologist might be expected to know or understand intuitively. It's not just   ##
## the core algorithums that are unfamiliar (e.g. Bloom filters), but also all the little hacks their   ##
## authors use to squeeze the most performance out of their code. If a Biologist cannot understand the  ##
## methods, how can they be expected to understand the results? Furthermore, all this complexity has a  ##
## knock-on effect on the user experience, as less development time was spent on the interface. Even as ##
## a Computer Scientist I struggled to figure out which parameters I needed to run most k-mer analysis  ##
## tools, often because some options (like how much RAM to use) isn't known until AFTER the analysis is ##
## completed! While some authors address these issues and try to help you out (like the authors of the  ##
## excellent Jellyfish tool), other tools who's names I won't mention will just exit with no error if   ##
## the right parameters aren't specified. During tests, one very popular tool even exited without error ##
## halfway through the analysis, but had left an unfinished (but still valid) output file!?             ##
## 2) Due to randomness inherent in these fast algorithums, k-mer tools often give different results    ##
## every time the same data is analysed (with the same tool). This would be acceptible for a quick look ##
## at the data, but what we really want for this project is a constant, reproducable "result" file for  ##
## a given seq data input file.                                                                         ##
## 3) Their outputs are in complex binary formats which do not lend themselves to manual analysis,      ##
## or analysis by other tools in high-level languages like Python. A good format should make it easy    ##
## for others to start hacking on it quickly and easily, allowing them to 'discover' new analyses.      ##
## 4) They take many hours to run, and for me, anything that takes more than 2 hours to run may as well ##
## take 20. I'm not going to wait around by the computer for the result either way, so the time savings ##
## of a few hours are not worth all of this added complexity, randomness, and downstream inflexability. ##
##                                                                                                      ##
## So in an attempt to address these issues, this program focuses on two things: simplicity, and a      ##
## really small and easy-to-work-with output format. This is what we came up with:                      ##                                                                        ##
'''                                                                                                     ##
##                               row  A   C   T   G   #   Seq                                           ##
##                              -----------------------------                                           ##
##                                0   1   2   0   3   7                                                 ##
##                                1   0   5   0   4   2   C                                             ##
##                                2   0   7   0   6   4                                                 ##
##                                3   0   0   0   0   1   ACCC                                          ##
##                                4   0   0   0   0   1   ACCC                                          ##
##                                5   0   0   0   0   1   C                                             ##
##                                6   0   0   0   0   1   ACCC                                          ##
##                                7   0   8   0   0   2                                                 ##
##                                8   0   0   0   0   1                                                 ##
'''                                                                                                     ##
## So it's just a table, with 6 columns: A, C, T, G, #, and Seq ("row" is just the row number, it's not ##
## actually stored anywhere). You start on row 0, and use one of the A C T G "warp pipe" columns to go  ##
## to a new row. If you have ever read one of those "create your own adventure" books, you'll recognize ##
## this format immediately :)                                                                           ##
## The "#" column tells us how many times the fragment that this row represents was seen in the data.   ##
## The Seq column allows us to group multiple fragments together into 1 row. For example, if there were ##
## many fragments that started with AAAAA, but only one of which ended in 'TCGA', that different ending ##
## could be stored in a single row, rather than 4 rows (one for each base). It might be confusing at    ##
## first, particularly if you forget that the row itself represents a base due to the warp pipe used to ##
## get to it, but once you "get it" you'll see it's really really simple.                               ##
##                                                                                                      ##
## This data structure is refered to in Computer Science as a 'trie', hence "ACGTrie" :)                ##
##                                                                                                      ##
## So say you were interested in counting how many "A"s there were in the data. You'd start on row 0,   ##
## and from there goto the A warp pipe, which has a "1" in it. This means to go to row 1. Then from row ##           
## 1, which represents "A", we can see in the # column there is a "2". This means there are 2 As in the ##
## input data. There's also some Seq in this row, a "C", which means that if we wanted to look up "AC"  ##
## we would also end up here, and that would also be a 2. You can confirm this on the original table up ##
## top. In fact it's worth taking a minute or two to try out some of the subfragments in the original   ##
## table, and seeing if you can get the same result in the ACGTrie format by-hand :)                    ##
##                                                                                                      ##
## But where this format really shines is when you have a lot of data. My estimates put it at 170x      ##
## smaller than the big table for real-world data, and much higher when your input data has very long   ##
## reads. And 170x smaller data means 170x less RAM (or 170x more data can be processed before you run  ##
## out of RAM) so you can take all that RAM-money and spend it on cool stuff like food and clothes.     ##
## I love food and clothes - good luck reading the rest of the code :)                                  ##
##########################################################################################################


###################
## up2bit format ##
##########################################################################################################
##                                                                                                      ##
## If you are a Biologist, brace yourself; winter is coming and you're about to go over the wall :T     ##
##                                                                                                      ##
## So, in the above table you can see DNA sequences like 'ACCC'. In the actual table that is saved on   ##
## your hard drive, we can't store little A/C/T/Gs, because computers only understand 1s and 0s.        ##
## We would need like, 2s and 3s as well. And hard drives dont support 2s and 3s. Pretty sure they dont ##
## anyways. To be honest, I have never really looked. So whatever, we have to encode our DNA in 1s and  ##
## 0s somehow. One mechanism is to use "ASCII". You might have heard of ASCII. It's a pretty big deal.  ##
## The idea is that you always use 8 bits of binary for every letter you want to store, a bit like how  ##
## DNA coding for amino acids always use 3 'bits' of DNA. 8 bits of binary gives you 256 possible       ##
## combinations. 00000000, 000000001, 00000010 ... 11111111. That's enough to store all the letters on  ##
## a Latin keyboard, and a few others letters and symbols too. The combination of 1s and 0s that define ##
## each letter/symbol is set out by the ASCII "encoding standard". This very file is encoded in ASCII!  ##
## See this --> A <-- ? That totally took up 8 bits of binary. Specifically "01000001". But if we only  ##
## needed to save 4 letters, A,C,T, and G, we don't really need to use 8 bits per letter, just like if  ##
## an organism only had 15 amino acids and 1 stop codon, its codons could be just 2 bases of DNA long.  ##
## Now i'm not one to size-shame large encodings - data is beautiful at any size - but it is kind of    ##
## wasteful to use 8 bits when we only need 2 - 00, 01, 10, and 11. This format actually already exists ##
## and is called "2bit". Maybe you heard of it when you downloaded a genome fasta file. It's basically  ##
## the same as ASCII but 4x smaller by storing DNA sequences in this compact encoding.                  ##
## But 2bit has a problem - we can't use it when the space to store our data is fixed, like our 64bit   ##
## ACGTrie table columns. No data in a 32bit column might be "00000000000000000000000000000000", but in ##
## 2bit that would mean "AAAAAAAAAAAAAAAA". It's like our row is literally screaming at us because it   ##
## has no way to encode silence. Poor row. So I made a format which supports a variable amount of DNA   ##
## info in a fixed space called 'up2bit'. Get it? While 2bit can encode 32 letters of DNA in 64bits, it ##
## cannot encode 31 or 30 or even 0 letters of DNA. But up2bit supports 'up too' 31 letters in 64bits.  ##
## Basically you sacrifice 1 letter for the ability to store a variable amount of letters.              ##
## It's very simple really - reading RIGHT TO LEFT (which is also how binary is counted) just encode up ##
## to 63 letters in regular 2bit format, then stick a '01' on the left-end. This is also a C, but that  ##
## isn't important. Any unused bits to the left of the 01 cap are 0s.                                   ##
## So "ACTG" in up2bit format would look like:                                                          ##
'''                                          ' -  G  T  C  A'                                           ##
##                                           '01 11 10 01 00'                                           ##
'''                                                                                                     ##
## Which when stored in a 32bit space, looks like:                                                      ##
'''                                                                                                     ##
##                           '00 00 00 00 00 00 00 00 00 00 00 01 11 10 01 00'                          ##
'''                                                                                                     ##
## Another example, here's 15 Ts:                                                                       ##
'''                                                                                                     ##
##                           '01 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10'                          ##
'''                                                                                                     ##
## Whilst no data at all would look like:                                                               ##
'''                                                                                                     ##
##                           '00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 01'                          ##
'''                                                                                                     ##
## Which is just the number 1 in binary. If you see in the code if statements asking if Seq == 1, this  ##
## is why - we're actually asking if it's empty or not.                                                 ##
##                                                                                                      ##
## So i'm sure you have questions. The first - "Why are we reading it right to left?"                   ##
## It's because in binary both "0111100100" and "00000000000000000000000111100100" are the same thing.  ##
## They both equal 484. You can add as many 0s to the left, and you'll still get the same number. This  ##
## means that an up2bit DNA sequence always has the same value whether you store it in 32bits, 64bits,  ##
## or even 3234983488bits. And 'no data' will also always equal '1'.                                    ##
##                                                                                                      ##
## Next question - "Why use A-00, C-01, T-10, G-11? It's not in alphabetical order, nor is it the order ##
## of the letters of your website!!?"                                                                   ##
## Well firstly, it is the order of the letters for the website (try and see!). Secondly, it's because  ##
## there's a neat trick you can do to convert those letters in ASCII directly into binary values. Here  ##
## is the binary encodings for A, C, G, and T in ASCII:                                                 ##
'''        -vv-                                                                                         ##
    A - 1000001                                                                                         ##
    C - 1000011                                                                                         ##
    G - 1000111                                                                                         ##
    T - 1010100                                                                                         ##
'''                                                                                                     ##
## See the 2nd and 3rd from the right values? 00, 01, 11, 10. In some programming languages like C, we  ##
## can cut out just those numbers very quickly. Unfortunately, in python, the code to do that, which is ##
'''                                                                                                     ##
    ord('T')>>1 &3                                                                                      ##
'''                                                                                                     ##
## is actually slower than just looking up the value in a dictionary. Also, i'd file this under "little ##
## hacks authors use to speed up their code", which make the code faster, but confuse the hell out of   ##
## anyone who doesn't know what it means. Having said that, I stuck with it to help another brother out ##
## who might use up2bit files and want's to make use of this trick in their own code.                   ##
##                                                                                                      ##
## Penultimate Question: Why use '01' as the capping bits? Why not '11' or '10'?                        ##
## This one is kind of long-winded, so unless you really care i'd skip it...                            ##
## Well sometimes, unfortunately, computers need to store negative numbers. How do you store a negative ##
## number when you only have 1s and 0s? Well, you say 'these numbers can be negative, and if they are I ##
## will make the very first bit a 1.' This means you can only count HALF as high, as the max number in  ##
## say, 8bit, would be 01111111, but you can also count all those numbers in negative too, so the range ##
## is still the same. Also you can make a number negative really quickly, by just flipping the first 0  ## 
## to a 1, or visa-versa. So when we define an encoding, say 8bit, 32bit, 64bit, etc, we first have to  ##
## say if it will be "signed" (can have negative numbers, first bit determines +/-), or "unsigned",     ##
## meaning they can only count 0 or higher (but twice as high as signed). Sometimes unsigned numbers    ##
## will also tell you about this cool indie band they're in, but it's best to just ignore them.         ##
## So why do we care? Well, in python (and other high-level programming languages), maths is generally  ##
## done using signed 64bit numbers. Therefore, if we want to store a python number in a NumPy array     ##
## that is anything other than signed 64bit, we're going to have to do some checking to make sure the   ##
## number will fit (not too high, or negative, etc), and this check takes time. However, if we are just ##
## stuffing signed 64bit numbers into a signed 64bit space, we don't need to do it. Thus, NumPy is alot ##
## faster working with signed 64bit spaces than it is unsigned 64bit spaces (or any other sized spaces) ##
## If the cap is therefore '01', we can never have a negative number, even if all the space is used up  ##
## by Gs (11). No negative numbers means we do not have weird issues like adding 1 to a negative number ##
## and the number getting smaller. Also it means you can put up2bit in both signed and unsigned spaces  ##
## without needing to check if it will fit (it always will).                                            ##
##                                                                                                      ##
## Final question: "If ACGTrie uses 64bits to store DNA, which means 1bp encoded in the warp pipe used  ##
## to get to that row, and 31bp of DNA in the Seq column, what happens if the DNA to be added is more   ##
## than 32 letters long?"                                                                               ##
## You just simply cut your sequence at the 32nd base, the 33rd base is encoded in the next warp pipe,  ##
## and the extra N bases are added to the next rows seq, and so on. So a 99bp fragment would need       ##
## atleast 4 rows to store all its data. However that is only true for 64bit columns. There's nothing   ##
## stopping you from editing the code and using 128bits or even higher for the seq column - then you    ##
## could store a lot more DNA before we have to break and go to a new row. While 64bit columns are good ##
## for fragments around 50bp in length. I highly recommend increasing this to something higher if your  ##
## fragments are longer than 50bp. Really, it depends on your data. If you do change the NumPy size,    ##
## just remember to change all occurences in this file of '32' to half of whatever you choose your new  ##
## size will be, and then all occurences of '64' in this file to your new size. Everything else, as     ##
## Steve Jobs would say, should "just work".                                                            ##
##                                                                                                      ##
## This concludes all the information I can provide about the up2bit encoding. If you can think of an   ##
## improvement to the format let me know, and we can break backwards compatibility together! \o/        ##
##                                                                                                      ##
##########################################################################################################



####################################################
## Optimization 1 : Pre-processing & adding rows ##
##########################################################################################################
##                                                                                                      ##
## Now we get to the good stuff! The first thing you need to know about how ACGTrie works is that it    ##
## *ALWAYS* takes it's data on the standard in (stdin) as strings of DNA, which it will put into a trie ##
## in memory, before writing it to disk. The program that feeds ACGTrie this data is called a           ##
## "pre-processor". When you give ACGTrie a SAM/BAM file to anaylse with "--input ./some.bam", you are  ##
## actually asking ACGTrie to be a pre-processor, which will start another "--cpu N" ACGTrie processes  ##
## to accept its data. The reason we do this is so you can write your own pre-processors very easily,   ##
## for example a FASTA file pre-processor, or perhaps something that reads an SQL database and prints   ##
## out DNA fragments. Whatever you have, a pre-processor can be made in just a couple of lines of code. ##
## However, pre-processors can also do a whole lot more. For example they could feed only DNA starting  ##
## with an 'A' to an ACGTrie process, making the output file only 1/4er of the size of the full trie,   ##
## and thus only 1/4er of the RAM to make. Repeat the process for 'C', 'G', and 'T', before combining   ##
## all these sub-tries in a smart way, and huge tries can be built on modest hardware. This is how we   ##
## get parallelism with --cpu N.                                                                        ##
##                                                                                                      ##
## By giving you the ability to write your own pre-processors, you have maximum flexability in what is  ##
## added to your ACGTries. For example, if you wanted to make a trie where the reads on the reverse     ##
## strand of the BAM file get reverse-complimented before being added to the trie (to look for motifs), ##
## this is really easy! If you want to clip reads to be no longer than 20bp long to make a more compact ##
## ACGTrie, no problem! Just feed in the DNA you want to count, and sepurate read counting complexity   ##
## from read gathering complexity.                                                                      ##
##                                                                                                      ##
## But there are other, more significant optimizations that pre-processing can do. For example, if our  ##
## read was 'ACGT', we wouldn't want to just find/create the ACGT row and +1 it's count - we would      ##
## also want to add counts to all the sub-fragments of the read too, which for 'ACGT' would look like:  ##
'''                      'A', 'C', 'G', 'T', 'AC', 'CG', 'GT', 'ACG', 'CGT', 'ACGT'                    '''
## So while one option is for your pre-processor to pass the above subfragments to ACGTrie, there is a  ##
## much better way, which is get your pre-processor to fragment the string like this:                   ##
'''                                      'ACGT', 'CGT', 'GT', 'T'                                      '''
## (which is to basically trim bases off the left-hand side on-by-one) then pass just those 4 fragments ##
## to an ACGTrie process with the "--walk" parameter set. This will make ACGTrie use a different adding ##
## function which instead of just incrementing the last row in the trie that represents the DNA, will   ##
## actually incriment all the rows it visits as it walks through the trie to find the last row.         ##
## So for 'CGT', the first C row also gets +1 and then the CG row gets a +1, and then finally the end   ##
## CGT row too. Fragmenting and adding in this way has the same result as adding all 10 sub-fragments   ##
## individually, but is obviously much faster since we only walk the trie 4 times and not 10. Also, the ##
## longer your reads, the more efficient the --walk optimization becomes, although of course this all   ##
## only applies if you want to count subfragments. Maybe you actually just want to count fragments with ##
## no subfragmenting at all, in which case just use ACGTrie without --walk and not do any fragmenting   ##
## in your pre-processor. Easy :)                                                                       ##
##                                                                                                      ##
## Or maybe you do want to count all subfragments, but you don't have control over your pre-processor.  ##
## It's just a 'dumb' SQL database churning our DNA. In this case you can get ACGTrie to do all the     ##
## subfragmenting for you, with the --fragment parameter, which will also automatically use the --walk  ##
## optimization. With --fragment, you would just pass your whole read, in this case 'ACGT', to the      ##
## ACGTrie process. Just remember that if you do use --fragment and a dumb pre-processor, you can not   ##
## also have your pre-processor only pass the fragments that start with an A (as in the above RAM       ##
## optimization example) because the sub-fragments are not guarenteed to also start with an A.          ##
##                                                                                                      ##
## The code to do --fragment's subfragmenting, as well as the two variants of adding logic (with and    ##
## without the --walk optimization) can be found below.                                                 ##
##                                                                                                      ##
##########################################################################################################

## Framents the input sequence for addRowWalk. Used when ACGTrie called with --fragment 
def subfragment(DNA,count):
    for l in range(len(DNA)-1,-1,-1):
        seqChunks[DNA[l:]] += count   # You will read about seqChunks next

## Adds DNA sequence/counts to the trie, and all subfragments it visits along the way.
## This is by far the most complicated part of the entire program. If you have any questions 
## please don't hesitate to ask me in IRC at http://ac.gt/chat

#@jit
def addRowWalk(DNA,count,nextRowToAdd,subfragments):
    dna = []
    for char in DNA: dna.append(('A','C','T','G').index(char))
    #dna = [('A','C','T','G').index(char) for char in DNA]                           ## Turn DNA into a list of 2bit numbers.
    subfragments += len(DNA)                                                        ## Just a counter that keeps a tally of how many fragments went into the trie. Goes in header.
    row = 0                                                                         ## We always start on row 0.
    o = 0
    while True:
        #if seqLen == 0: break                                                      ## This is here to remind myself never to impliment it (sometimes you need to split even if seq = [])
        if SEQ[row] == 1:                                                           ## There is no sequence data in this row, and thus 3 possible options going forward.
                                                                                    ## 1: We also have no data to add (seq == []). Just do nothing and we'll break out of the loop.
                                                                                    ## 2: We have data to add (seq != []) and will need to TAKE a pipe to the next row.
                                                                                    ## 3: We have data to add (seq != []) but will need to MAKE a pipe and a next row and all rows thereafter.                                                                              
            COUNT[row] += count                                                     ## ( But all of them start by adding +1 to this rows's # count )
            if o == len(dna): return nextRowToAdd, subfragments
            temp = int((A,C,T,G)[dna[o]][row])                                  ## Getting the value for thisRow[seq[0]] takes time. We only want to get it once.
            if temp != 0:                                                       ## Type 2. Happens 97% of the time, which is why its the first thing we try.
                row = temp
                o  += 1                                                         ## Note the seq becomes a bit shorter, because we need to lose a base following the pipe.
                continue
            else:                                                               ## Type 3. Happens just 2.3% of the time. 
                for y in xrange(0, len(dna)-o, 32):                             ##         For every 32-letter chunk of sequence,
                    (A,C,T,G)[dna[o+y]][row] = nextRowToAdd                     ##         Make a new pipe in the old row
                    COUNT[nextRowToAdd]       = count                           ##         Add count to the new row
                    tmp = 0
                    for t,z in enumerate(xrange(y+1,y+32 if y+32 < len(dna)-o else len(dna)-o)): tmp += dna[o+z]<<(2*t)
                    SEQ[nextRowToAdd]         = tmp + (1<< (31 if 31 < len(dna)-o-1 else len(dna)-o-1)*2)
                    row                       = nextRowToAdd                                          ##         Set row to the newly created row, and repeat.
                    nextRowToAdd             += 1
                return nextRowToAdd, subfragments
                                                                                    ## as checking would waste time 99.3% of the time. If the above loop fails because seq[0] cant be found, we
                                                                                    ## must have a Type 1 situation, so the best thing to do is just break whenever we hit that error.

        else:                                                                       ## There IS sequence data in this row, now we have three possible options going forward:
                                                                                    ## 1: The DNA in the row is longer than the DNA in the fragment. (45%)
                                                                                    ## 2: The DNA in the row is shorter than the DNA in the fragment. (56%)
                                                                                    ## 3: The DNA in the row is the same length as the DNA in the fragment.
                                                                                    ## Within these three branches, the shorter DNA can match the longer, or, it can not.
            temp = int(SEQ[row])
            up2bit = []
            for x in range(0,temp.bit_length()-2,2): up2bit.append((temp >> x) & 3)
            #up2bit = [((temp >> x) & 3) for x in range(0,temp.bit_length()-2,2)]
            if len(dna)-o < len(up2bit):                                            ## Type 1. DNA in the row is longer than DNA in the fragment. This gives us two possibilities:
                for x in xrange(0,len(dna)-o):
                    if up2bit[x] != dna[o+x]:
                        A[nextRowToAdd] = A[row]
                        C[nextRowToAdd] = C[row]
                        T[nextRowToAdd] = T[row]
                        G[nextRowToAdd] = G[row]
                        COUNT[nextRowToAdd] = COUNT[row]
                        tmp = 0
                        for t,z in enumerate(xrange(x+1,len(up2bit))): tmp += up2bit[z]<<(2*t)
                        SEQ[nextRowToAdd] = tmp + (1 << (len(up2bit)-x-1)*2)
                        A[row] = 0
                        C[row] = 0
                        T[row] = 0
                        G[row] = 0
                        (A,C,T,G)[up2bit[x]][row] = nextRowToAdd
                        (A,C,T,G)[dna[o+x]][row] = nextRowToAdd+1
                        COUNT[row] += count
                        tmp = 0
                        for t,z in enumerate(xrange(0,x)): tmp += up2bit[z]<<(2*t)
                        SEQ[row]  = tmp + (1 << x*2)
                        COUNT[nextRowToAdd+1] = count                                                    ## This second new row already has 0,0,0,0 for its pipes, so we just have to set the
                        tmp = 0
                        for t,z in enumerate(xrange(o+x+1,len(dna))): tmp += dna[z]<<(2*t)
                        SEQ[nextRowToAdd+1]  = tmp + (1 << (len(dna)-o-x-1)*2)
                        nextRowToAdd           += 2
                        return nextRowToAdd, subfragments
                A[nextRowToAdd]                     = A[row]
                C[nextRowToAdd]                     = C[row]
                T[nextRowToAdd]                     = T[row]
                G[nextRowToAdd]                     = G[row]
                COUNT[nextRowToAdd]                 = COUNT[row]
                tmp = 0
                for t,z in enumerate(xrange(len(dna)-o+1,len(up2bit))): tmp += up2bit[z]<<(2*t)
                SEQ[nextRowToAdd]                   = tmp + (1 << (len(up2bit)-1+o-len(dna))*2)
                A[row]                              = 0
                C[row]                              = 0
                T[row]                              = 0
                G[row]                              = 0
                COUNT[row]                         += count
                tmp = 0
                for t,z in enumerate(xrange(0,len(dna)-o)): up2bit[z]<<(2*t)
                SEQ[row]                            = tmp + (1 << (len(dna)-o)*2)
                (A,C,T,G)[up2bit[len(dna)-o]][row]  = nextRowToAdd
                nextRowToAdd                       += 1
                return nextRowToAdd, subfragments
            elif len(dna)-o > len(up2bit):                                                ## Type 2. The DNA in the fragment is longer than DNA in the row. This gives us again two possibilities:
                #if dna[o:len(up2bit)+o] == up2bit:                                                     ## 1) The row's DNA, although shorter than the fragments's DNA, matches the fragment.
                for x in xrange(0,len(up2bit)):
                    if dna[o+x] != up2bit[x]:
                        A[nextRowToAdd] = A[row]
                        C[nextRowToAdd] = C[row]
                        T[nextRowToAdd] = T[row]
                        G[nextRowToAdd] = G[row]
                        COUNT[nextRowToAdd] = COUNT[row]
                        tmp = 0
                        for t,z in enumerate(xrange(x+1,len(up2bit))): tmp += up2bit[z]<<(2*t)
                        SEQ[nextRowToAdd] = tmp + (1 << (len(up2bit)-x-1)*2)
                        A[row] = 0
                        C[row] = 0
                        T[row] = 0
                        G[row] = 0
                        (A,C,T,G)[up2bit[x]][row] = nextRowToAdd
                        (A,C,T,G)[dna[o+x]][row]    = nextRowToAdd+1
                        COUNT[row] += count
                        tmp = 0
                        for t,z in enumerate(xrange(0,x)): tmp += up2bit[z]<<(2*t)
                        SEQ[row]  = tmp + (1 << x*2)
                        nextRowToAdd  += 1
                        o += x                                                                      ## We know we need to make a new row for the DNA the fragment had that the original row's
                        for y in xrange(0, len(dna)-o, 32):                                         ## Seq did not, but we don't know how long the DNA we need to add is. We might need to 
                            (A,C,T,G)[dna[o+y]][row] = nextRowToAdd                              ## make several new rows to fit it all in, and so that is what this loop is doing.
                            COUNT[nextRowToAdd]      = count
                            tmp = 0
                            for t,z in enumerate(xrange(y+1,y+32 if y+32 < len(dna)-o else len(dna)-o)): tmp += dna[o+z]<<(2*t)
                            SEQ[nextRowToAdd]        = tmp + (1<< (31 if 31 < len(dna)-o-1 else len(dna)-o-1)*2)
                            row                      = nextRowToAdd
                            nextRowToAdd            += 1
                        return nextRowToAdd, subfragments
                COUNT[row] += count
                o += len(up2bit)                                                                    ## cut our fragment's DNA to be just the stuff the fragment has extra,
                temp = (A,C,T,G)[dna[o]][row]                                                       ## and check to see if this row has a warp pipe to where we want to go next.
                if temp == 0:                                                                       ## If there is no warp pipe...
                    for y in xrange(0, len(dna)-o, 32):                                             ## We create one (or as many as we need in a chain)
                        (A,C,T,G)[dna[o+y]][row] = nextRowToAdd
                        COUNT[nextRowToAdd]      = count
                        tmp = 0
                        for t,z in enumerate(xrange(y+1,y+32 if y+32 < len(dna)-o else len(dna)-o)): tmp += dna[o+z]<<(2*t)
                        SEQ[nextRowToAdd]        = tmp + (1<< (31 if 31 < len(dna)-o-1 else len(dna)-o-1)*2)
                        row                      = nextRowToAdd
                        nextRowToAdd            += 1
                        return nextRowToAdd, subfragments
                else:                                                                               ## But if there is a warp pipe,
                    row = temp                                                                      ## we just take it :)
                    o += 1
                    continue
            else:                                                                   ## Type 3. The DNA in the fragment is the same length as the DNA in the row. This gives us again two possibilities:
                for x in xrange(0,len(up2bit)):
                    if dna[o+x] != up2bit[x]:
                        A[nextRowToAdd]           = A[row]
                        C[nextRowToAdd]           = C[row]
                        T[nextRowToAdd]           = T[row]
                        G[nextRowToAdd]           = G[row]
                        COUNT[nextRowToAdd]       = COUNT[row]
                        tmp = 0
                        for t,z in enumerate(xrange(x+1,len(up2bit))): tmp += up2bit[z]<<(2*t)
                        SEQ[nextRowToAdd]         = tmp + (1 << (len(up2bit)-x-1)*2)
                        A[row]                    = 0
                        C[row]                    = 0
                        T[row]                    = 0
                        G[row]                    = 0
                        (A,C,T,G)[up2bit[x]][row] = nextRowToAdd
                        (A,C,T,G)[dna[o+x]][row]  = nextRowToAdd+1
                        COUNT[row]               += count
                        tmp = 0
                        for t,z in enumerate(xrange(0,x)): tmp += up2bit[z]<<(2*t)
                        SEQ[row]                  = tmp + (1 << x*2)
                        COUNT[nextRowToAdd+1]     = count
                        tmp = 0
                        for t,z in enumerate(xrange(o+x+1,len(dna))): tmp += dna[z]<<(2*t)
                        SEQ[nextRowToAdd+1]       = tmp + (1 << (len(dna)-o-x-1)*2)
                        nextRowToAdd             += 2                                                                   ## row, so we don't have to do the loop we did at the end of Type 2.
                        return nextRowToAdd, subfragments
                COUNT[row] += count                                                   ## 1) Very simply, we just increment the count and we're done. 
                return nextRowToAdd, subfragments

############################################
## Optimization 2 : Multiplying & Sorting ##
##########################################################################################################
##                                                                                                      ##
## When reading BAM files and other chromosomal-position-sorted data, it is quite common to see the     ##
## exact same DNA fragment being read several times over as multiple reads pile up in the same location ##
## Even if the reads themselves are not identical duplicates, their sub-fragments almost certainly will ##
## if they span the same genomic region. So another optimization method is to have your pre-processor   ##
## keep a buffer of fragments and counts, such that when the same fragment is seen more than once, it's ##
## count is incremented before it is even handed over to ACGTrie. When it is handed over, instead of    ##
## using the usual standard input format of:                                                            ##
'''                                     'ACATCGA\nACATCGA\nCATCGAG'                                    '''
## your pre-processor uses the CSV format of:                                                           ##
'''                                       'ACATCGA,2\nCATCGAG,1'                                       '''
## Where \n is the new-line character or, in laymans terms, the invisible character you put in when you ##
## want to make a new line by pressing the enter key.                                                   ##
## This is automatically detected by ACGTrie, and doing so means we only have to do two trie-walks      ##
## instead of three, because in the first walk we just +2 to every row rather than +1. This works       ##
## regardless of any other pre-optimizations you do or do not use such as --fragment or --walk, and for ##
## data like BAM files can improve performances by basically as much as your read depth. It might even  ##
## be a good idea to read you whole BAM file first to make this buffer, before passing anything at all  ##
## to ACGTrie. The size of your buffer is up to your pre-processor! With --input on SAM/BAM files, the  ##
## default is to make a buffer of 10,000 reads, which are stored in the "seqChunks" dictionary. By the  ##
## way, if you were a Computer Scientist, you might call this optimization Run Length Encoding or "RLE" ##
##                                                                                                      ##
## Finally, there is one more trick we can use to increase our trie-creation speed, and that is to sort ##
## the fragments so the largest fragments are added first. In all honesty, I don't really know why this ##
## works - I just tried it one day to see if it made a difference and it did. Computers are weird.      ##
## Fortunately, because we are already creating a buffer for the RLE, a quick sort isn't hard :)        ##
## Again though, this sorting optimization does not work if you use --fragment since it has to be done  ##
## by the pre-processor!                                                                                ##
##                                                                                                      ##
##########################################################################################################

## The seqChunks dictionary stores sub-fragments and counts and is used whenever ACGTrie is called with
## --fragment, or when an --input SAM/BAM file is provided and ACGTrie needs to be a pre-processor. 
seqChunks = collections.defaultdict(int)

## emptyCache is called whenever there is enough buffered DNA. Does longest-fragment-first sorting too.
def emptyCache(add,nextRowToAdd,subfragments):
    global seqChunks
    for fragment,count in sorted(seqChunks.items(), reverse=True, key=lambda t: len(t[0])):
        nextRowToAdd,subfragments = add(fragment,count,nextRowToAdd,subfragments)
    seqChunks = collections.defaultdict(int)

## Just a little function to figure out how much RAM we're using and display it in a human-readable way.
def getRAM():
    if usingNumpy: num = A.nbytes*7
    else:
        # num = ffi.sizeof(A)
        num = ctypes.sizeof(A)*7
    for unit in [' ',' K',' M',' G',' T']:
        if abs(num) < 1024.0:
            return "%3.1f%sb" % (num, unit)
        num /= 1024.0

#############################
## Putting it all together ##
##########################################################################################################
##                                                                                                      ##
## Here the actual logic of the program begins. We first figure out if we are just accepting data for   ##
## adding to a trie, or if we are going to act as a pre-processor for a SAM/BAM file...                 ##
##                                                                                                      ##
##########################################################################################################

# Our table starts off at 280Mb and grows by 280Mb every time it runs out of space :)
if usingNumpy:
    A     = numpy.zeros(args.rows, dtype='uint32')
    C     = numpy.zeros(args.rows, dtype='uint32')
    G     = numpy.zeros(args.rows, dtype='uint32')
    T     = numpy.zeros(args.rows, dtype='uint32')
    COUNT = numpy.zeros(args.rows, dtype='uint32')
    SEQ   = numpy.zeros(args.rows, dtype='int64' )
else:
    Array32 = ctypes.c_uint32 * args.rows
    Array64 = ctypes.c_int64  * args.rows
    A,C,G,T,COUNT,SEQ = Array32(),Array32(),Array32(),Array32(),Array32(),Array64()
# else:
#     ffi = cffi.FFI()
#     A = ffi.new("uint32_t [args.rows]")
#     C = ffi.new("uint32_t [args.rows]")
#     G = ffi.new("uint32_t [args.rows]")
#     T = ffi.new("uint32_t [args.rows]")
#     COUNT = ffi.new("uint32_t [args.rows]")
#     SEQ = ffi.new("int64_t [args.rows]")

A[0],C[0],G[0],T[0],COUNT[0],SEQ[0] = 0,0,0,0,0,1 ## Row 0, the root row/node.
nextRowToAdd = 1
print '   [ RAM used @ ' + getRAM() + ' ]'

# Read first row to see if DNA comes with counts. Create appropriate stdin generator.
firstFragment = sys.stdin.readline().strip().split(',')
if   len(firstFragment) == 1: stdin = zip(sys.stdin, itertools.repeat(1)) ; firstFragment.append(1)
elif len(firstFragment) == 2: stdin = csv.reader(sys.stdin, delimiter=',') ; firstFragment[1] = int(firstFragment[1])
else: print 'ERROR: I do not understand this kind of stdin format :U'; exit()

# What kind of adding function to use?
# if args.walk or args.fragment: add = addRowWalk
# else:                          add = addRow
add = addRowWalk

startTime = time.time()
linesRead = 1
fragmentAvg = len(firstFragment[0])
subfragments = 0
if args.fragment:
    subfragment(firstFragment[0],firstFragment[1])
    for DNA,count in stdin:
        linesRead += 1
        fragmentAvg += len(DNA)
        DNA,count = DNA.rstrip(),int(count)
        subfragment(DNA,count)
        if linesRead % 10000 == 0:
            emptyCache(add,nextRowToAdd,subfragments)
    emptyCache(add,nextRowToAdd,subfragments)
else:
    nextRowToAdd,subfragments = add(firstFragment[0],firstFragment[1],nextRowToAdd,subfragments)
    for DNA,count in stdin:
        if nextRowToAdd + 100 > len(A):
            if usingNumpy:
                # Originally I wrote to disk then pulled it back because most methods to
                # extend C structured array requires having both the old and new array in
                # memory at once for some period of time, but since our columns are individual
                # arrays, we only require 28% of the RAM that doing it as a single array takes. 
                A = numpy.concatenate((A,numpy.zeros(10000000, dtype='uint32')))
                C = numpy.concatenate((C,numpy.zeros(10000000, dtype='uint32')))
                G = numpy.concatenate((G,numpy.zeros(10000000, dtype='uint32')))
                T = numpy.concatenate((T,numpy.zeros(10000000, dtype='uint32')))
                COUNT = numpy.concatenate((COUNT,numpy.zeros(10000000, dtype='uint32')))
                SEQ = numpy.concatenate((SEQ,numpy.zeros(10000000, dtype='int64')))
            else:
                newSize = A._length_+10000000
                ctypes.resize(A, ctypes.sizeof(A._type_)*newSize)
                A = (A._type_*newSize).from_address(ctypes.addressof(A))
                ctypes.resize(C, ctypes.sizeof(C._type_)*newSize)
                C = (C._type_*newSize).from_address(ctypes.addressof(C))
                ctypes.resize(G, ctypes.sizeof(G._type_)*newSize)
                G = (A._type_*newSize).from_address(ctypes.addressof(G))
                ctypes.resize(T, ctypes.sizeof(T._type_)*newSize)
                T = (A._type_*newSize).from_address(ctypes.addressof(T))
                ctypes.resize(COUNT, ctypes.sizeof(COUNT._type_)*newSize)
                COUNT = (A._type_*newSize).from_address(ctypes.addressof(COUNT))
                ctypes.resize(SEQ, ctypes.sizeof(SEQ._type_)*newSize)
                SEQ = (SEQ._type_*newSize).from_address(ctypes.addressof(SEQ))
            print '   [ RAM used @ ' + getRAM() + ' ]'
        linesRead += 1
        fragmentAvg += len(DNA)
        DNA,count = DNA.rstrip(),int(count)
        #######################
        nextRowToAdd,subfragments = add(DNA,count,nextRowToAdd,subfragments)

print 'Done in: ', time.time()-startTime
## Trie is now built. Time to add a header and call it a day:
print sum(A),sum(C),sum(T),sum(G),sum(COUNT),sum(SEQ)
exit()


header = 'HEADER_START\n'
headerInfo = {
    #'structs': [(x,str(y[0])) for x,y in sorted(trie.dtype.fields.items(),key=lambda k: k[1])],     ## List of name/structs used in NumPy array. Readers of ACGTrie files should check this as might not always be standard.
    'fragments': linesRead,                                                                         ## Total number of lines in the SAM/BAM file read
    'subfragments': subfragments,                                                                   ## The number of subfragments added by ACGTrie
    'fragmentAvgLen': int(fragmentAvg/linesRead),                                                   ## The average length of each fragment (not subfragment)
    'rows': int(nextRowToAdd),                                                                      ## The number of rows in the ACGTrie table
    'analysisTime': str(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')),                     ## When the table was made.
    'countOverflow': {},                                                                            ## ?
    'warpOverflow': {}                                                                              ## ?
}
theJson = json.dumps(headerInfo,sort_keys=True, indent=4)       ## The header format here is exactly the
if theJson.count('\n') < 100: header += theJson                 ## same as all the header formats in all
else: header += json.dumps(headerInfo,sort_keys=True)           ## AC.GT projects. JSON data, typically
while header.count('\n') < 99: header += '\n'                   ## laid out nicely to be human-readable,
header += 'HEADER_END\n'                                        ## with HEADER_END on the 100th line so
                                                                ## "head -100 ./file" or "tail +101 ./file"
fileA     = open(args.output + '.A', 'wb');     fileA.write(header)
fileC     = open(args.output + '.C', 'wb');     fileC.write(header)
fileG     = open(args.output + '.G', 'wb');     fileG.write(header)
fileT     = open(args.output + '.T', 'wb');     fileT.write(header)
fileCOUNT = open(args.output + '.COUNT', 'wb'); fileCOUNT.write(header)
fileSEQ   = open(args.output + '.SEQ', 'wb');   fileSEQ.write(header)

if usingNumpy:
    if sys.byteorder == 'big':
        print '   [ Flipping eggs. ]'; ## Little Endians 4 lyfe yo.
        A.byteswap(True);C.byteswap(True);G.byteswap(True);T.byteswap(True);COUNT.byteswap(True);SEQ.byteswap(True)
    fileA.write(A[:nextRowToAdd])
    fileC.write(C[:nextRowToAdd])
    fileG.write(G[:nextRowToAdd])
    fileT.write(T[:nextRowToAdd])
    fileCOUNT.write(COUNT[:nextRowToAdd])
    fileSEQ.write(SEQ[:nextRowToAdd])
else:
    fin32 = ctypes.c_uint32 * int(nextRowToAdd)
    fin64 = ctypes.c_int64 * int(nextRowToAdd)
    fileA.write(fin32.from_address(ctypes.addressof(A)))
    fileC.write(fin32.from_address(ctypes.addressof(C)))
    fileG.write(fin32.from_address(ctypes.addressof(G)))
    fileT.write(fin32.from_address(ctypes.addressof(T)))
    fileCOUNT.write(fin32.from_address(ctypes.addressof(COUNT)))
    fileSEQ.write(fin64.from_address(ctypes.addressof(SEQ)))

fileA.close()
fileC.close()
fileG.close()
fileT.close()
fileCOUNT.close()
fileSEQ.close()

print linesRead, subfragments, nextRowToAdd
'''

rstrip vs [:-1] ?

pip update numpy?

Depth column? See how much time it adds.
Why RAM used so much more than what gets written to disk?
numpy.array instead of uint32?
(pytables/array.array?)
--interpreter /path/to/pypy/or/cython   (must have hts-python not pysam. Maybe would work with samtools only.)

OVERFLOWS

Speed test against other Trie-making programs. Perhaps have a speed and RAM competition.
Consider using 5x array.array - one for each column - and pypy for fast access.
Also pytables instead of numpy?
Get addRow to work.
'''


'''
Future ideas:
DEPTH array? Has the depth of the row/node in the trie. uint8 as depth unlikely to be bigger than 256.
BACK array ? Has the row number for the parent row. uint32.
Any way to compress this data but still randomly-access it's contents? Numpy's memmap only works on uncompressed. Maybe use HD5F? Overhead? (No pytables support in pypy)
Function just like addNodeWalk but only +1s to the last node, not the intermediate nodes.
Sort table. I think compression will work a lot better if the table is sorted by SEQ or something. Branches of trei dont need sorting, just array whilst correcting A/C/G/T/BACK pipes.
Compress table. The table is actually not at optimal size after being made. Rows 'deeper' in the trie can be merged to nodes higher up which now have space because they were once.
Delete rows. Delete rows from the table with a count lower than X. Delete noise and save a lot of space? Pipes need to be maintained.
'''