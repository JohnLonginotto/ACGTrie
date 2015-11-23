#!/usr/bin/env python

####################
## For Biologists ##
##########################################################################################################
##                                                                                                      ##
## Hello :)                                                                                             ##
## This file is designed to help you understand how ACGTrie works on a practical level. Nothing about   ##
## ACGTrie is hard to understand, but for Scientists unfamiliar with coding/python this should provide  ##
## a solid understanding of exactly what is happening to your data, and hopefully from there you can    ##
## experiment with both the algorithm and the outputs (since you will know exactly how they are made).  ##
##                                                                                                      ##
## The downside to this simplicity and verbosity however is that, because we don't specify every little ##
## detail in how to do things, the Python intepreter (that turns written Python code into the low-level ##
## machine code required by the CPU) has to make some educated guesses - and it doesn't always get it   ##
## right. The result is code which, although clear for humans, is alot slower than it could be. Athough ##
## ACGTrie_LEARN.py will give the exact same output as ACGTrie_FAST.py, it will take longer to do so.   ##
## Thus, for your actual analysis, use ACGTrie_FAST! Use ACGTrie_LEARN.py to learn the algorithum :)    ##
##                                                                                                      ##
##########################################################################################################

#####################
## For Programmers ##
##########################################################################################################
##                                                                                                      ##
## This script is an executable lab book.                                                               ##
## It is not Object-Oriented, Modularized, or PEP 8 compliant. And that is OK.                          ##
## Push requests should focus on improving readibility for Biologists and not optimising the code for   ##
## speed, memory usage, collaboration, etc.                                                             ##
##                                                                                                      ##
##########################################################################################################


## First we import some extra modules that come with python that we will be using:
import csv
import sys
import json
import time
import argparse
import datetime
import collections

## We use high-performance C arrays to store the data, but accessing C arrays from Python can
## be done using three different modules: cffi, ctypes, and numpy. ACGTrie can use any, however
## cffi or numpy is usually the fastest if you have it. ACGTrie_FAST will determine the best array
## for your system, but in LEARN you can set it manually to see the difference for yourself.

arrayKind = 'cffi' #    <-- change this if you like.  valid options are 'cffi', 'numpy' and 'ctypes'
if arrayKind == 'cffi':   import cffi    ; print '   [ Using cffi ]'
if arrayKind == 'numpy':  import numpy   ; print '   [ Using numpy ]'
if arrayKind == 'ctypes': import ctypes  ; print '   [ Using ctypes ]'

## This code can be run in both the standard Python intepreter (called CPython) or the PyPy interpreter.
## PyPy is a lot better at figuring out how to turn Python into machine code, and so it is a lot faster than
## the standard CPython intepreter. However, PyPy does not yet support all the external modules that CPython does. 
## In ACGTrie_FAST this is all handled automatically, but since _LEARN lets you choose the C array module yourself, 
## we need to do some checks to make sure you haven't changed the above to something incompatible with your system.
if "__pypy__" in sys.builtin_module_names:
    if arrayKind == 'numpy':
        try:    import numpy
        except: print 'ERROR: You are using pypy, but you havent installed numpy for pypy (numpypy) yet! Go to the site and grab the latest version to continue :)'; exit()
    if arrayKind == 'ctypes': 
        print 'WARN: The pypy ctypes module does not support ctypes.resize(). Please make sure you use enough RAM (--rows) when you start, because we cant add more later :('
else:
    if arrayKind == 'numpy':
        try:    import numpy
        except: print 'ERROR: You do not have numpy installed! Grab it via pip install numpy :)'; exit()
    if arrayKind == 'cffi':
        try:    cffi.FFI().memmove
        except: print 'WARN: The version of cffi you have does not support memmove(). Please make sure you use enough RAM (--rows) when you start, because we cant add more later :('

## Parse user-supplied command line options:
parser = argparse.ArgumentParser(     description="Put in a BAM file or fragments of DNA, get out an ACGTrie table.")
parser.add_argument("-o", "--output", metavar='/path/to/output.trie', help='Required. Output filename.')
parser.add_argument("-q", "--quiet",  action='store_true',            help="Optional. Do not print anything to stdout (except for errors)")
parser.add_argument("--rows",         default=10000000, type=int,     help="Optional. How big to make the ACGTrie table initially.")
parser.add_argument("--walk",         action='store_true',            help="Optional. Tells ACGTrie to incrementally add to trie.")
parser.add_argument("--fragment",     action='store_true',            help="Optional. Tells ACGTrie to make subfragments itself.")
parser.add_argument("--debug",        action='store_true',            help="Optional. Debug mode. Builds hashtable as well as trie.")
args = parser.parse_args()

if args.output == None: print '''
    You havn't provided a path or filename for your output!
    Actually, your pre-processor should have done this for you.
    Are you using a preprocessor? Preprocessors are designed to get
    sequencing data out of the popular data formats and pass it quickly 
    and efficiently into one or more ACGTrie processes in an optimized way. 
    Please visit the ac.gt website to find a pre-processor for your data format
    (BAM, FASTA or whatever), or, if we have not made one for your data format yet,
    write one and we will put it up on the website for all to use and benefit from! :-)
'''; exit()


##############################################
## Two algorithums, with and without --walk ##
##########################################################################################################
##                                                                                                      ##
## So you have read the docs on the up2bit format, the output format, and preprocessors - (if not, you  ##
## are going to want to do that first!) - now you want to know how DNA is actually added to the trie.   ##
## You already know that there are 2 algorithums, one that's highly optimized for DNA composition (all  ##
## subfragments in the input fragment) activated with --walk, and one that just adds the fragment and   ##
## not the subfragments. What you might not know is that the algorithum for --walk is simpler than the  ##
## algorithum that just adds the fragment. You will see why as we continue, but for now we will start   ##
## with the --walk algorithum.                                                                          ##
##                                                                                                      ##
##########################################################################################################


## This function turns a string of DNA into a list of 2bit numbers.
def dna_list(DNA):
    return [ ('A','C','T','G').index(char) for char in DNA ]

## This function turns the up2bit number in the SEQ column of the trie into a list of 2bit values like the above.
def up2bit_list(value):
    value = int(value) # Numbers in cffi arrays don't have the .bit_length property, so we convert to int.
    return [((value >> x) & 3) for x in range(0,value.bit_length()-1,2)] # This uses bit-shifting to get the pairs of 2 bits at a time.

## Once we have both the incoming DNA fragment and the DNA in the row as 2bit lists, we can compare them base-by-base.
## This function sees if the fragment DNA and the SEQ DNA matches or not, and if not, where in the list they mismatch.
def firstNonMatching(longerDNA,shorterDNA):
    for x,y in enumerate(longerDNA):
        if y != shorterDNA[x]:
            return  x           ## If they dont match, x is the base where the mismatch happens
    else:   return -1           ## If there is no DNA (because we used it up while traveling the warp to get to this row), x is -1
    return x                    ## If they match, x is the length of the shorter DNA.

## This function turns out list of 2bit values into a single up2bit number (so it can be added to the trie)
def list_up2bit(twobits):
    up2bit = 1                             ## The '01' cap.
    for twobit in reversed(twobits):       ## Kind of similar to string_up2bit except we already start with a list of
        up2bit = (up2bit << 2) + twobit    ## 2bit numbers and we just need to cap it and stack them all together.
    return up2bit

## OK thats all the low-level stuff covered. With the above we can take our first DNA fragment, turn it to a list, then turn it to an up2bit number,
## and add it to the trie. But what about the second fragment? As there is already DNA in the trie, we have to add our new DNA in such a way that the old
## rows still make sense. If the new fragment is the same but longer than the old, we have to extend the old data with a new warp pipe and new rows.
## If the new fragment is the same but shorter than the old fragment, we have to split a row in the old data so that new COUNT values can be added without
## effecting the COUNT values of the extra stuff in the old row.
## If the new DNA is a partial match, but also has some unique DNA of its own, we have to split an existing row into three rows: shared DNA, path to old DNA,
## and a path to the new DNA.
## The following functions deal with splitting rows:  

## If we have to split, we start by copying the data that is unique to the existing row into a new row:
def copyRow(row,nextRowToAdd,up2bit,x):
    A[nextRowToAdd] = A[row]
    C[nextRowToAdd] = C[row]
    T[nextRowToAdd] = T[row]
    G[nextRowToAdd] = G[row]
    COUNT[nextRowToAdd] = COUNT[row]
    SEQ[nextRowToAdd] = list_up2bit(up2bit[x+1:len(up2bit)]) ## Only DNA after the mismatch (x+1:) is copied over (the mismatch itself will
                                                             ## be encoded in the warp pipe).

## We then edit the original row to have either 1 or 2 warps to pipes pointing to the 1 or 2 new rows.
## The first new row is the one we copied above, the old data. The second hasn't been made yet, but if
## we have to make it, we already know where it will be - nextRowToAdd!
## the row above.
def freshRow(row,count,DNA1,DNA2,x):
    global nextRowToAdd
    A[row] = 0
    C[row] = 0
    T[row] = 0
    G[row] = 0
    COUNT[row] += count
    (A,C,T,G)[DNA1[x]][row] = nextRowToAdd         ## first pipe
    SEQ[row] =  list_up2bit(DNA1[:x])              ## original row DNA before the mismatch ([:x])
    nextRowToAdd += 1
    if DNA2 != None:
        (A,C,T,G)[DNA2[x]][row] = nextRowToAdd  ## second pipe (if needed)

## And finally, if we did make two warp pipes before, we need to make the second row (DNA unique to the incoming fragment):
def secondRow(count,seq,x):
    global nextRowToAdd
    COUNT[nextRowToAdd] = count
    SEQ[nextRowToAdd]   = list_up2bit(seq[x+1:len(seq)])
    nextRowToAdd += 1

## Sometimes we know we can add DNA to the trie without interfearing with any other paths.
## However, we might have DNA that is longer than what the SEQ column can hold (usually 32bp of DNA if the SEQ array is 64bits).
## When this happens we use this function to just keep adding as many new rows as it takes until all the DNA is added to the trie.
def restOfSeq(seq,row,count):
    global nextRowToAdd
    for y in xrange(0, len(seq), 32):
        (A,C,T,G)[seq[y]][row]    = nextRowToAdd
        COUNT[nextRowToAdd]       = count
        SEQ[nextRowToAdd]         = list_up2bit(seq[y+1:y+32])
        row                       = nextRowToAdd
        nextRowToAdd             += 1

## Right, that's all the mid-level stuff covered. Argubly that was the hardest bit. Now comes the part
## where we wire it all together to build the trie using the --walk algorithum. Here we go!
def addRowWalk(DNA,count):
    global nextRowToAdd
    seq = dna_list(DNA)                                                         ## Turn fragment DNA into a list of 2bit numbers.
    row = 0                                                                     ## We always start on row 0.
    while True:                                                                 ## Here we start our main loop, where we walk through rows in the trie reading and or adding data until we hit a 'return'.
        if SEQ[row] == 1:                                                       ## There is no sequence data in this row, and thus 3 possible options going forward:
            COUNT[row] += count                                                 ## ( But all of them start by adding +1 to this rows's COUNT )
            if not seq: return                                                  ## Type 1: We have no more DNA in our fragment, so we're totally done!
            warp_pipe = (A,C,T,G)[seq[0]][row]                                  ## Otherwize we get the value from the next warp pipe we need to take.
            if warp_pipe:                                                       ## Type 2: If there is a warp pipe (value != 0):
                row = warp_pipe                                                 ##         the row becomes the row the warp pipe pointed to.
                seq = seq[1:]                                                   ##         and our seq becomes a bit shorter, as we need to lose a base when following the pipe.
                continue                                                        ##         and finally 'continue', meaning we repeat the loop again starting on the new row.
            else:                                                               ## Type 3: If there is no warp pipe to take:
                restOfSeq(seq,row,count)                                        ##         We need to just make as many rows as needed to get all our fragment DNA into the trie.
                                                                                ##
        else:                                                                   ## If there IS sequence data in this row, things become a bit more complicated. We have another 3 possible options:
            up2bit = up2bit_list(SEQ[row])                                      ## First turn the SEQ DNA in this row into a list like the fragment DNA already is ('seq')
            if len(seq) < len(up2bit):                                          ## Type 1: The DNA in the fragment is shorter than the DNA in the row. This gives us 2 possibilities:
                if up2bit[:len(seq)] == seq:                                    ##         1) The Fragment's DNA, although shorter than the row's DNA, matches the row.
                    copyRow(row,nextRowToAdd,up2bit,len(seq))                   ##            So we just copy the extra DNA in the row out to a new row, with the original COUNT value
                    freshRow(row,count,up2bit,None,len(seq))                    ##            and fix up the original row to have just the shared DNA in the SEQ, an incremented COUNT value, and a warp to the new row.
                else:                                                           ##         2) The Fragment's DNA is shorter AND doesn't match the existing SEQ data.
                    x = firstNonMatching(seq,up2bit)                            ##            So we will have to find the mismatch,
                    copyRow(row,nextRowToAdd,up2bit,x)                          ##            Copy out to a new the matching (if any) DNA.
                    freshRow(row,count,up2bit,seq,x)                            ##            Update the original row with either 1 or 2 warp pipes.
                    secondRow(count,seq,x)                                      ##            And add the second new row, containing the DNA that was in the fragment but not the original SEQ.
            elif len(seq) > len(up2bit):                                        ## Type 2: The DNA in the fragment is longer than DNA in the row. This gives us again 2 possibilities:
                if seq[:len(up2bit)] == up2bit:                                 ##         1) The row's DNA, although shorter than the fragments's DNA, matches the fragment.
                    COUNT[row] += count                                         ##            In which case we can just increment the count of this row,
                    seq         = seq[len(up2bit):]                             ##            cut our fragment's DNA to be just the stuff the fragment has extra,
                    warp_pipe   = (A,C,T,G)[seq[0]][row]                        ##            and check to see if this row has a warp pipe to where we want to go next.
                    if not warp_pipe:                                           ##            If it doesnt...
                        restOfSeq(seq,row,count)                                ##            we add as many rows as required to get all our fragment's DNA into the trie
                    else:                                                       ##            But if there is a warp pipe,
                        row = warp_pipe                                         ##            we just take it as before :)
                        seq = seq[1:]                                           ##
                        continue                                                ##
                else:                                                           ##         2) Otherwise, the rows's DNA and the fragments's DNA do NOT match, and we have to make 2 new rows...
                    x = firstNonMatching(up2bit,seq)                            ##            First we find the mismatch,
                    copyRow(row,nextRowToAdd,up2bit,x)                          ##            make a new row for the DNA which was not in the fragment's DNA 
                    freshRow(row,count,up2bit,seq,x)                            ##            update the original row to reflect this, and create two warp pipes. One to the row we made above,
                    restOfSeq(seq[x:],row,count)                                ##            and one to the first of potentially many rows we need to make to get all our fragment's DNA into the trie.
            else:                                                               ## Type 3. So here, the DNA in the fragment is the same length as the DNA in the row. Again 2 possibilities:
                if up2bit == seq:                                               ##         1) We have an identical match,
                    COUNT[row] += count                                         ##            so we just increment the count and we're done!
                else:                                                           ##         2) Otherwise, we have to make 2 new rows, just like in the above.
                    x = firstNonMatching(seq,up2bit)                            ##            This could have also been x = firstNonMatching(up2bit,seq)
                    copyRow(row,nextRowToAdd,up2bit,x)                          ##            yadda
                    freshRow(row,count,up2bit,seq,x)                            ##            yadda
                    secondRow(count,seq,x)                                      ##            yadda
        return                                                                  ## And we're done! Exit the While loop to get a new fragment from the stdin :) 




## OK that was all the code required to make the trie using --walk.
## Somewhat counter-intuitively, the code for adding just the fragments and not the sub-fragments is a bit
## more complicated. We basically have to treat the very last base in the fragment as special, and add it at the end because
## it always has to have its own row with no other DNA in it's SEQ. All the DNA before the last base has to have a count of 0,
## or at least we +0 to each row we walk over to get to that last, special, base.
## I will only comment the important differences between this algorithum and the above.

## First, we use this function to add the last base. It has to know if it just adds, or splits, etc.
## It's basically a combination of copyRow and freshRow.
def addLastBase(row,lastBase,count):
    global nextRowToAdd
    warp_pipe = (A,C,T,G)[lastBase][row]
    if warp_pipe:
        row = warp_pipe
        if SEQ[row] == 1:
            COUNT[row] += count
        else:
            up2bit = up2bit_list(SEQ[row])
            copyRow(row,nextRowToAdd,up2bit,0)
            A[row] = 0
            C[row] = 0
            T[row] = 0
            G[row] = 0
            COUNT[row] += count
            (A,C,T,G)[up2bit[0]][row] = nextRowToAdd
            SEQ[row] =  1
            nextRowToAdd += 1
    else:
        (A,C,T,G)[lastBase][row] = nextRowToAdd
        SEQ[nextRowToAdd] = 1
        COUNT[nextRowToAdd] = count
        nextRowToAdd += 1


## And this is the main program logic:
def addRow(DNA,count):
    global nextRowToAdd
    seq = dna_list(DNA)
    lastBase = seq[-1]                                                  ## So we start but grabbing that last base,
    seq = seq[:-1]                                                      ## and chopping it off the fragment.
    row = 0                                                             ## Now we can add the fragment to the trie as usual (with a count of 0)
    while True:
        if SEQ[row] == 1:
            if len(seq) == 0:
                addLastBase(row,lastBase,count)                         ## but when we would normally return, we add that last base.
                return
            warp_pipe = (A,C,T,G)[seq[0]][row]
            if warp_pipe:
                row = warp_pipe
                seq = seq[1:]
                continue
            else:
                restOfSeq(seq,row,0)
                addLastBase(nextRowToAdd-1,lastBase,count)              ## here.
                return
        else:
            up2bit = up2bit_list(SEQ[row])
            if len(seq) < len(up2bit):
                if up2bit[:len(seq)] == seq: 
                    copyRow(row,nextRowToAdd,up2bit,len(seq))
                    freshRow(row,0,up2bit,None,len(seq))
                    addLastBase(row,lastBase,count)                     ## here.
                else:
                    x = firstNonMatching(seq,up2bit)
                    copyRow(row,nextRowToAdd,up2bit,x)
                    freshRow(row,0,up2bit,seq,x)
                    secondRow(0,seq,x)
                    addLastBase(nextRowToAdd-1,lastBase,count)          ## here.
            elif len(seq) > len(up2bit):
                if seq[:len(up2bit)] == up2bit:
                    seq         = seq[len(up2bit):]
                    warp_pipe   = (A,C,T,G)[seq[0]][row]
                    if not warp_pipe:
                        restOfSeq(seq,row,0)
                        addLastBase(nextRowToAdd-1,lastBase,count)      ## here.
                        return
                    else:
                        row = warp_pipe
                        seq = seq[1:]
                        continue
                else:
                    x = firstNonMatching(up2bit,seq)
                    copyRow(row,nextRowToAdd,up2bit,x)
                    freshRow(row,0,up2bit,seq,x)
                    restOfSeq(seq[x:],row,0)
                    addLastBase(nextRowToAdd-1,lastBase,count)          ## here.
                    return
            else:
                if up2bit == seq:
                    addLastBase(row,lastBase,count)                     ## here.
                else:
                    x = firstNonMatching(seq,up2bit)
                    copyRow(row,nextRowToAdd,up2bit,x)
                    freshRow(row,0,up2bit,seq,x)
                    secondRow(0,seq,x)
                    addLastBase(nextRowToAdd-1,lastBase,count)          ## and finally, here :)
        break


####################
## Bonus Features ##
##########################################################################################################
##                                                                                                      ##
## So we've done all the hard stuff - and to be honest it is pretty confusing. It took me a good solid  ##
## month to figure out both algorithums and how to impliment them effectively (in FAST). So if you feel ##
## like you don't really get straight away, well thats partially unavoidable, and partially my fault!   ##
##                                                                                                      ##
## Fortunately the rest of the code is dead-easy. Its just helper functions, reading the stdin to pass  ##
## fragments to the functions above, and a little bit of debugging stuff for developers. Good luck! :)  ##
##                                                                                                      ##
##########################################################################################################

## Figures out how much RAM we are currently using and returns it in a human-readable way.
def getRAM():
    if args.quiet: return # we dont print anything when --quiet is used...
    else:
        if   arrayKind == 'cffi':   num = ffi.sizeof(A)*7
        elif arrayKind == 'numpy':  num = A.nbytes*7
        elif arrayKind == 'ctypes': num = ctypes.sizeof(A)*7
        for unit in [' ',' K',' M',' G',' T']:
            if abs(num) >= 1024.0: num /= 1024.0
            else: print "   [ RAM used @ %3.1f%sb ]" % (num, unit); return

## These two functions are used when ACGTrie is called with --fragment. Since we are fragmenting 
## ourselves here, it also makes sense to keep a very basic buffer called seqChunks going too.
def subfragment(DNA,count):
    for l in range(len(DNA)-1,-1,-1):
        seqChunks[DNA[l:]] += count

## Which we empty with emptyCache.
def emptyCache(add,nextRowToAdd):
    global seqChunks
    for fragment,count in sorted(seqChunks.items(), reverse=True, key=lambda t: len(t[0])):
        add(fragment,count)
    seqChunks = collections.defaultdict(int)

## If the user wants to --debug the code, we do this by also creating a hashtable whilst
## building the tree. Of course, hash tables can only store a few thousand fragments before
## they fill up all the RAM, but that is enough to test the algorithum works.
## The way we add data to the hashtables for the two algorithums is described below:
def addRowDebug(DNA,count):
    hashTable[DNA] += count
    addRow(DNA,count)
def addRowWalkDebug(DNA,count):
    for l in range(1,len(DNA)+1):
        hashTable[DNA[:l]] += count
    addRowWalk(DNA,count)

## On start up, we allocate a block of RAM memory for the trie to use. If we add so many
## rows that we start running out of space, we need to allocate more memory.
## Unfortunately, the only way to "extend" a C array is to liturally allocate a new, bigger block,
## copy the old data over to it, and delete the old block.
## Fortunately, this can be done very fast (no need to re-compute anything), and because we 
## break each column of our trie-table into its own little array, it doesn't require 2x the total RAM used
## to do it, just 2x the largest column, which is the SEQ column (which is why we copy it first).
def growTrie():
    global A
    global C
    global G
    global T
    global COUNT
    global SEQ
    growBy = 10000000     ## You can change this if you want to experiment with growing your arrays in smaller/larger chunks.
    currentSize = len(A)
    newSize = currentSize + growBy
    if currentSize == 4294967295:
        print '''
    ERROR: MAXIMUM ROW LIMIT REACHED!
           First of all, let me congratulate you on having more than 120Gb of RAM. Thats crazy.
           Unfortunately, a single ACGTrie output cannot have more than this number of rows unless
           it was prefix-split and rejoined by a preprocessor. If you are already prefix-splitting,
           rerun the analysis at a higher prefix-split depth. If this is not possible, please
           talk to one of the developers about your situation and how you want to proceed, and we'll
           work with you to figure out the best course of action :)
              '''; exit()
    if newSize > 4294967295: newSize = 4294967295; growBy = 4294967295 - currentSize
    if arrayKind == 'numpy':
        SEQ = numpy.concatenate((SEQ,numpy.zeros(growBy, dtype='int64')))
        A = numpy.concatenate((A,numpy.zeros(growBy, dtype='uint32')))
        C = numpy.concatenate((C,numpy.zeros(growBy, dtype='uint32')))
        G = numpy.concatenate((G,numpy.zeros(growBy, dtype='uint32')))
        T = numpy.concatenate((T,numpy.zeros(growBy, dtype='uint32')))
        COUNT = numpy.concatenate((COUNT,numpy.zeros(growBy, dtype='uint32')))
    elif arrayKind == 'ctypes':
        ctypes.resize(SEQ, ctypes.sizeof(SEQ._type_)*newSize)
        SEQ = (SEQ._type_*newSize).from_address(ctypes.addressof(SEQ))
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
    elif arrayKind == 'cffi':
        try:
            tempSEQ   = SEQ;   SEQ   = ffi.new("int64_t[]",  newSize); ffi.memmove( SEQ,   tempSEQ,   currentSize*(ffi.sizeof(SEQ)/newSize)   )
            tempA     = A;     A     = ffi.new("uint32_t[]", newSize); ffi.memmove( A,     tempA,     currentSize*(ffi.sizeof(A)/newSize)     )
            tempC     = C;     C     = ffi.new("uint32_t[]", newSize); ffi.memmove( C,     tempC,     currentSize*(ffi.sizeof(C)/newSize)     )
            tempG     = G;     G     = ffi.new("uint32_t[]", newSize); ffi.memmove( G,     tempG,     currentSize*(ffi.sizeof(G)/newSize)     )
            tempT     = T;     T     = ffi.new("uint32_t[]", newSize); ffi.memmove( T,     tempT,     currentSize*(ffi.sizeof(T)/newSize)     )
            tempCOUNT = COUNT; COUNT = ffi.new("uint32_t[]", newSize); ffi.memmove( COUNT, tempCOUNT, currentSize*(ffi.sizeof(COUNT)/newSize) )
        except AttributeError:
            print 'ERROR! Your version of python does not fully support cffi! Please upgrade, or set the intial amount of RAM to something larger.'
            exit()
    getRAM()

## Stats are nice. Stats in your face is better.
## We record some stats so that anything unusual will be displayed on program exit.
class fileStats:
    def __init__(self,DNA,count):
        self.start = time.time()
        self.linesRead = 1
        self.linesProcessed = count
        self.datetime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    def add(self,DNA,count):
        self.linesRead += 1
        self.linesProcessed += count
    def result(self):
        return (self.linesRead,self.linesProcessed,time.time()-self.start,self.datetime)


######################
## Create the table ##
##########################################################################################################
##                                                                                                      ##
## There are many ways to create an in-memory table. Some are faster than others. Some are part of the  ##
## standard library and don't require any other modules to be downloaded. Some only work on the latest  ##
## versions of python. But after we have created the table, we read/write data in exactly the same way  ##
## The three ways we currently support are numpy, ctypes and cffi.                                      ##
##                                                                                                      ##
##########################################################################################################

# Our table starts off at 280Mb (or the number of rows provided by --rows) :)
if arrayKind == 'numpy':
    A     = numpy.zeros(args.rows, dtype='uint32')
    C     = numpy.zeros(args.rows, dtype='uint32')
    G     = numpy.zeros(args.rows, dtype='uint32')
    T     = numpy.zeros(args.rows, dtype='uint32')
    COUNT = numpy.zeros(args.rows, dtype='uint32')
    SEQ   = numpy.zeros(args.rows, dtype='int64' )

elif  arrayKind == 'ctypes':
    Array32 = ctypes.c_uint32 * args.rows
    Array64 = ctypes.c_int64  * args.rows
    A = Array32()
    C = Array32()
    G = Array32()
    T = Array32()
    COUNT = Array32()
    SEQ = Array64()

elif  arrayKind == 'cffi':
    ffi = cffi.FFI()
    A = ffi.new("uint32_t[]", args.rows)
    C = ffi.new("uint32_t[]", args.rows)
    G = ffi.new("uint32_t[]", args.rows)
    T = ffi.new("uint32_t[]", args.rows)
    COUNT = ffi.new("uint32_t[]", args.rows)
    SEQ = ffi.new("int64_t[]", args.rows)

## Create row 0, the root row/node:
A[0],C[0],G[0],T[0],COUNT[0],SEQ[0] = 0,0,0,0,0,1
nextRowToAdd = 1
countOverflow = {} ## The maximum value in a uint32 column is 4294967295, which means we cant increase the COUNT or visit a row via a warp above this value  
warpOverflow = {}  ## without changing all columns to uint64, practically doubling the amount of memory needed. So instead of doing that, we can just keep
                   ## track of the few rows that do go over (for the COUNT) and break up out trie into smaller prefix-tries to totally avoid the warp-pipe issue.
                   ## But for now both are only done by the preprocessors, not ACGTrie itself. This may change soon.
getRAM()           ## Print the initial amount of ram being used.


#########################
## Start reading stdin ##
##########################################################################################################
##                                                                                                      ##
## First we need to figure out what format the input data is coming in as - just DNA with newlines, or  ##
## in CSV format with counts too. Once we've done that, we need to check the execution parameters, and  ##
## enter the appropriate stdin read loop, with the right adding function.                               ##
##                                                                                                      ##
##########################################################################################################

# Read first row to see if DNA comes with counts.
firstFragment = sys.stdin.readline().strip().split(',')
if len(firstFragment) == 1:                                                             ## No counts!
    firstFragment.append(1)                                                             ## Add a 1 to the firstFragment
    stdin = ((line.rstrip(),1) for line in sys.stdin)                                   ## and add a 1 for every time we read a fragment from the stdin.
elif len(firstFragment) == 2:                                                           ## There are counts!
    firstFragment[1] = int(firstFragment[1])                                            ## Convert the count from a string to a number
    stdin = ( (DNA,int(count)) for DNA,count in csv.reader(sys.stdin, delimiter=',') )  ## and do the same every other time we read from the stdin.
else:
    print 'ERROR: I do not understand this kind of stdin format :U'
    exit()

# What kind of adding function to use?
if args.walk or args.fragment: add = addRowWalk                                         ## addRowWalk is used to efficiently create DNA composition tries.
else:                          add = addRow                                             ## addRow is used to make treis of fragments and not DNA composition.

stats = fileStats(firstFragment[0],firstFragment[1])                                    ## Initiliaze the stats counter class with data from the first fragment. Kicks off the timer too.

# If we are using --debug, we need a hash table too
if args.debug:
    hashTable = collections.defaultdict(int)
    if add == addRowWalk: add = addRowWalkDebug
    else:                 add = addRowDebug

## If ACGTrie has to fragment the reads to get DNA composition itself:
if args.fragment:
    seqChunks = collections.defaultdict(int)
    subfragment(firstFragment[0],firstFragment[1])
    for DNA,count in stdin:
        stats.add(DNA,count)
        subfragment(DNA,count)
        if len(seqChunks) > 100000:
            if nextRowToAdd + 100000 > len(A): growTrie()
            emptyCache(add,nextRowToAdd)
    emptyCache(add,nextRowToAdd)

## Else, we add data as-is without fragmenting:
else:
    add(*firstFragment)
    for DNA,count in stdin:
        stats.add(DNA,count)
        add(DNA,count)
        if nextRowToAdd + 100 > len(A): growTrie()

## OK, all done - lets get the stats:
linesRead, linesProcessed, duration, startDatetime = stats.result()

## Print some output:
if not args.quiet:
    print 'Done in: ', duration
    print 'Lines read: ', linesRead
    print 'Fragments processed: ', linesProcessed
    print 'Sum of each array (A/C/G/T/COUNT/SEQ): ',sum(A),sum(C),sum(T),sum(G),sum(COUNT),sum(SEQ)
if linesProcessed > 4294967295 and (args.walk or args.fragment):
    print '''
    ERROR: The largest COUNT value a row can have is 4294967295, but you have added more fragments than that to the trie!
           This means that some of the COUNT values - certainly the root node - have rolled over and started on 0 again. 
           There are three fixes: re-run the data using prefix-branch splitting (see ACGTrie_BAM) to split the data into 
           smaller chunks which dont overflow, (and let the preprocessor figure out the overflows); manually find the rows
           whos COUNT has rolled over and fix the header manually; ask the authors of ACGTrie to impliment a new addRow or 
           addRowWalk function that checks for overflows on every update.
           I'm sorry for the inconvenience, I just really didn't expect this to ever happen, and to fix it for the 0.0001% 
           of time happens will slow everything down for the other 99.9999% of the time where its not needed :(
    '''

########################
## Write out the data ##
##########################################################################################################
##                                                                                                      ##
## Like all ACGT software, it's very important that output files contain a JSON header with metadata    ##
## about the analysis. Having said that, im not a huge fan of 6 files all with (almost) the same header ##
## ... but it is much safer for ACGTrie to always add a header, rather than trust the preprocessors to  ##
## do the right thing and make one. If preprocessors HAVE to parse it out to concat, they're likely to  ##
## put it back in again properly too.                                                                   ##
##                                                                                                      ##
##########################################################################################################

header32 = 'HEADER_START\n';
header64 = 'HEADER_START\n';
uint32_head = {
    'structs': 'uint32',
    'linesRead': linesRead,
    'fragments': linesProcessed,
    'rows': nextRowToAdd,
    'analysisTime': startDatetime,
    'analysisDuration': duration,
    'countOverflow': countOverflow,
    'warpOverflow': warpOverflow
}
int64_head = dict(uint32_head); int64_head['structs'] = 'int64'
uint32_json = json.dumps(uint32_head,sort_keys=True, indent=4)       ## The header format here is exactly the
int64_json = json.dumps(int64_head,sort_keys=True, indent=4)       ## The header format here is exactly the
if uint32_json.count('\n') < 100:
    header32 += uint32_json                 ## same as all the header formats in all
    header64 +=  int64_json                  ## same as all the header formats in all
else: 
    header32 += json.dumps(uint32_json,sort_keys=True)                 ## same as all the header formats in all
    header64 += json.dumps( int64_json,sort_keys=True)                 ## same as all the header formats in all
while header32.count('\n') < 99:
    header32 += '\n'
    header64 += '\n'
header32 += 'HEADER_END\n'                                        ## with HEADER_END on the 100th line so
header64 += 'HEADER_END\n'                                        ## "head -100 ./file" or "tail +101 ./file"
fileA     = open(args.output + '.A', 'wb');     fileA.write(header32)
fileC     = open(args.output + '.C', 'wb');     fileC.write(header32)
fileG     = open(args.output + '.G', 'wb');     fileG.write(header32)
fileT     = open(args.output + '.T', 'wb');     fileT.write(header32)
fileCOUNT = open(args.output + '.COUNT', 'wb'); fileCOUNT.write(header32)
fileSEQ   = open(args.output + '.SEQ', 'wb');   fileSEQ.write(header64)

## We write the data in the arrays to disk - however we only write up to nextRowToAdd, otherwise
## we would be writing a whole bunch of 0s to disk ;-)
if arrayKind == 'numpy':
    fileA.write(A[:nextRowToAdd])
    fileC.write(C[:nextRowToAdd])
    fileG.write(G[:nextRowToAdd])
    fileT.write(T[:nextRowToAdd])
    fileCOUNT.write(COUNT[:nextRowToAdd])
    fileSEQ.write(SEQ[:nextRowToAdd])
elif arrayKind == 'ctypes':
    fin32 = ctypes.c_uint32 * int(nextRowToAdd)
    fin64 = ctypes.c_int64 * int(nextRowToAdd)
    fileA.write(fin32.from_address(ctypes.addressof(A)))
    fileC.write(fin32.from_address(ctypes.addressof(C)))
    fileG.write(fin32.from_address(ctypes.addressof(G)))
    fileT.write(fin32.from_address(ctypes.addressof(T)))
    fileCOUNT.write(fin32.from_address(ctypes.addressof(COUNT)))
    fileSEQ.write(fin64.from_address(ctypes.addressof(SEQ)))
elif arrayKind == 'cffi':
    fin32 = nextRowToAdd * (ffi.sizeof(A)/len(A))
    fin64 = nextRowToAdd * (ffi.sizeof(SEQ)/len(A))
    fileA.write(ffi.buffer(A,fin32))
    fileC.write(ffi.buffer(C,fin32))
    fileG.write(ffi.buffer(G,fin32))
    fileT.write(ffi.buffer(T,fin32))
    fileCOUNT.write(ffi.buffer(COUNT,fin32))
    fileSEQ.write(ffi.buffer(SEQ,fin64))

fileA.close()
fileC.close()
fileG.close()
fileT.close()
fileCOUNT.close()
fileSEQ.close()


###############
## Debugging ##
##########################################################################################################
##                                                                                                      ##
## Bugs are a natural part of life. The following functions will help you search out and squash them    ##
## if you call ACGTrie with --debug :)                                                                  ##
##                                                                                                      ##
##########################################################################################################

## Turn an up2bit number into a string of DNA letters
def up2bit_string(value):
    up2bits = [((value >> x) & 3) for x in range(0,value.bit_length()-1,2)]
    return ''.join([('A','C','T','G')[up2bit] for up2bit in up2bits])

## If the function is run, it will print out the trie in the same format as you saw in the docs.
def printTrie():
    print '\nrow     A    C    T    G  COUNT  SEQ\n-------------------------------------------'
    for row in xrange(0,nextRowToAdd):
        print str(  row ).rjust(4),
        print str(A[row]).rjust(4),
        print str(C[row]).rjust(4),
        print str(T[row]).rjust(4),
        print str(G[row]).rjust(4),
        print str(COUNT[row]).rjust(4) + '   ',
        print up2bit_string(SEQ[row])
    print '-------------------------------------------'
    ## We also print out 2 rows after the last row that was added.
    ## Why? Because bugs dont play by the rules, and anything is possible.
    for row in xrange(nextRowToAdd,nextRowToAdd+2):
        print str(  row ).rjust(4),
        print str(A[row]).rjust(4),
        print str(C[row]).rjust(4),
        print str(T[row]).rjust(4),
        print str(G[row]).rjust(4),
        print str(COUNT[row]).rjust(4) + '   ',
        print up2bit_string(SEQ[row])
    print '-------------------------------------------'

## A little function to get the counts of some DNA (and all of its prefixes)
## out of the trie. If you just needed the final value, this could be a lot simplier...
def getScore(string):
    row = 0                                                             ## We always start on row 0
    counts = []                                                         ## We will return a list of counts for each row we visit
    seq = [('A','C','T','G').index(char) for char in string]            ## Turn the string into a list of 2bit numbers
    while True:
        rowCount = COUNT[row]
        counts.append(rowCount)                                         ## Straight away we will add this row's # value to our rowCount list
        seqLen = len(seq)
        if seqLen == 0: break                                           ## If we have no more DNA in our input string, we're done :)
        up2bit = up2bit_list(SEQ[row])                                  ## Else we get some more data..
        up2bitLen = len(up2bit)
        nextPipe = (A,C,T,G)[seq[0]][row]
        if nextPipe and up2bitLen == 0:                                 ## If we have no more DNA in this row's Seq,
            row = nextPipe                                              ## take the warp pipe to the next row.
            seq = seq[1:]
            continue
        elif seqLen <= up2bitLen:                                       ## If we have more DNA in our row than our string...
            if up2bit[:seqLen] == seq: counts += [rowCount] * seqLen    ## 1) Check they match up. If so add this row's count N times.
            else:
                for x,y in enumerate(seq):                              ## 2) If they dont match up, find where they diverge
                    if y != up2bit[x]: break                            ##    and just add this row's count for that many times
                counts += [rowCount] *  x
        else:                                                           ## Finally, to be here we must have more DNA in our string
            if seq[:up2bitLen] == up2bit:                               ## than in our row. Thus as before we see if they match.
                counts += [rowCount] * up2bitLen                        ## If they do, add N counts to the counts list, and take
                row = (A,C,T,G)[seq[up2bitLen]][row]                         ## the next warp pipe out (if available)
                if row != 0:
                    seq = seq[up2bitLen+1:]
                    continue
            else:                                                       ## Otherwise find out where the two sequences diverge and
                for x,y in enumerate(up2bit):                           ## just add the appropriate number of counts to the counts
                    if y != seq[x]: break                               ## list before bottoming out of the while loop and breaking.
                counts += [rowCount] *  x
        break
    return counts

## If ACGTrie run with --debug, whilst the trie is being built we also build
## a hash table. The hashtable and the tree should match up. Note that this
## doesn't make sure the trie has things the hash table doesn't (although this
## is very unlikely), and of course making hashtables like this uses a LOT of RAM,
## so only run debug on a small (< 1 million fragments) of data :)
## Works whether you --fragment, --walk, or not.
def debug():
    for DNA,count in sorted(hashTable.items()):
        trieScore = getScore(DNA)
        if count != trieScore[-1] or len(DNA) != len(trieScore)-1:
            if count != trieScore[-1]: print 'DEBUG: COUNT mismatch found between trie and hash-table'
            if len(DNA) != len(trieScore)-1: print 'DEBUG: LENGTH mismatch found between trie and hash-table'
            print DNA,count,trieScore
            printTrie()
            return
    print 'Both trie and hashtable match up! :)'

if args.debug:
    sys.stdin = open('/dev/tty') ## We take data from the stdin, which gets closed when the preprocessor/etc closes the pipe. We need to re-open it for the user to type stuff.
    while True:
        userInput = raw_input('Type debug, print, or enter DNA to get its count starting from the root node: ')
        if userInput.strip('ACGT') != '' and userInput != 'debug' and userInput != 'print': 
            print 'Sorry, I didnt understand that - please try again.'
            continue
        if userInput == 'print': printTrie()
        elif userInput == 'debug': debug()
        else: print getScore(userInput)


'''
Future ideas:
DEPTH array? Has the depth of the row/node in the trie. uint8 as depth unlikely to be bigger than 256.
BACK array ? Has the row number for the parent row. uint32.
Any way to compress this data but still randomly-access it's contents? Numpy's memmap only works on uncompressed. Maybe use HD5F? Overhead? (No pytables support in pypy)
Function just like addNodeWalk but only +1s to the last node, not the intermediate nodes.

Post-processing ideas:
Sort table. I think compression will work a lot better if the table is sorted by SEQ or something. Branches of trei dont need sorting, just array whilst correcting A/C/G/T/BACK pipes.
Compact table. The table is actually not at optimal size after being made. Rows 'deeper' in the trie can be merged to nodes higher up which now have space because they were split.
Delete rows. Delete rows from the table with a count lower than X or a depth higher than X. Pipes need to be maintained.
'''
