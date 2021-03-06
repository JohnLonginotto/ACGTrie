#!/usr/bin/env python

import csv
import sys
import json
import time
import argparse
import tempfile
import datetime
import itertools
import collections

## Valid options are 'cffi', 'ctypes' and 'numpy'. 
arrayKind = 'cffi'

## We test to see if we can use the above C array module.
if arrayKind == 'cffi': import cffi      ; print '   [ Using cffi ]'
if arrayKind == 'numpy': import numpy    ; print '   [ Using numpy ]'
if arrayKind == 'ctypes': import ctypes  ; print '   [ Using ctypes ]'
if "__pypy__" in sys.builtin_module_names:
    if arrayKind == 'ctypes': print 'WARN: The pypy ctypes module does not support ctypes.resize(). Please make sure you use enough RAM (--rows) when you start, because we cant add more later :('
    if arrayKind == 'numpy':
        try: import numpy
        except: print 'ERROR: You are using pypy, but you havent installed numpy for pypy (numpypy) yet! Go to the site and grab the latest version to continue :)'; exit()
else:
    if arrayKind == 'numpy':
        try: import numpy
        except: print 'ERROR: You do not have numpy installed! Grab it via pip install numpy :)'; exit()
    if arrayKind == 'cffi':
        try: cffi.FFI().memmove
        except: print 'WARN: The version of cffi you have does not support memmove(). Please make sure you use enough RAM (--rows) when you start, because we cant add more later :('

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

###################
## Trie creation ##
##########################################################################################################
##                                                                                                      ##
## Below is the code we use to add data to the trie. Unlike ACGTrie_LEARN, addRow and addRowWalk do not ##
## call other functions when running. This helps the Python compiler to inline the code to effient C.   ##
## We also try to keep operations which copy memory (like slicing or summing) to a minimum, since high  ##
## performance python intepreters like pypy tend to have a more efficient (less aggressive) garbage     ##
## collection of unused things in memory.                                                               ##
##                                                                                                      ##
##########################################################################################################

## Finally, the main logic of the whole program - how to use the above to add data to the trie:

def addRowWalk(DNA,count):
    global nextRowToAdd
    dna = [ ord(char)>>1 &3 for char in DNA ]
    row = 0                                                                         ## We always start on row 0.
    o   = 0
    while True:
        #if seqLen == 0: break                                                      ## This is here to remind myself never to impliment it (sometimes you need to split even if seq = [])
        if SEQ[row] == 1:                                                           ## There is no sequence data in this row, and thus 3 possible options going forward.
                                                                                    ## 1: We also have no data to add (seq == []). Just do nothing and we'll break out of the loop.
                                                                                    ## 2: We have data to add (seq != []) and will need to TAKE a pipe to the next row.
                                                                                    ## 3: We have data to add (seq != []) but will need to MAKE a pipe and a next row and all rows thereafter.                                                                              
            COUNT[row] += count                                                     ## ( But all of them start by adding +1 to this rows's # count )
            if o == len(dna): return
            warp_pipe = int((A,C,T,G)[dna[o]][row])                                  ## Getting the value for thisRow[seq[0]] takes time. We only want to get it once.
            if warp_pipe:                                                       ## Type 2. Happens 97% of the time, which is why its the first thing we try.
                row = warp_pipe
                o  += 1                                                         ## Note the seq becomes a bit shorter, because we need to lose a base following the pipe.
                continue
            else:                                                               ## Type 3. Happens just 2.3% of the time. 
                for y in xrange(0, len(dna)-o, 32):                             ##         For every 32-letter chunk of sequence,
                    (A,C,T,G)[dna[o+y]][row]  = nextRowToAdd                    ##         Make a new pipe in the old row
                    COUNT[nextRowToAdd]       = count                           ##         Add count to the new row
                    SEQ[nextRowToAdd]         = sum([dna[o+z]<<(2*t) for t,z in enumerate(xrange(y+1,y+32 if y+32 < len(dna)-o else len(dna)-o))],(1<< (31 if 31 < len(dna)-o-1 else len(dna)-o-1)*2))
                    row                       = nextRowToAdd                                          ##         Set row to the newly created row, and repeat.
                    nextRowToAdd             += 1
                return
                                                                                    ## as checking would waste time 99.3% of the time. If the above loop fails because seq[0] cant be found, we
                                                                                    ## must have a Type 1 situation, so the best thing to do is just break whenever we hit that error.

        else:                                                                       ## There IS sequence data in this row, now we have three possible options going forward:
                                                                                    ## 1: The DNA in the row is longer than the DNA in the fragment. (45%)
                                                                                    ## 2: The DNA in the row is shorter than the DNA in the fragment. (56%)
                                                                                    ## 3: The DNA in the row is the same length as the DNA in the fragment.
                                                                                    ## Within these three branches, the shorter DNA can match the longer, or, it can not.
            temp = int(SEQ[row])
            up2bit = [((temp >> x) & 3) for x in range(0,temp.bit_length()-2,2)]
            if len(dna)-o < len(up2bit):                                            ## Type 1. DNA in the row is longer than DNA in the fragment. This gives us two possibilities:
                for x in xrange(0,len(dna)-o):
                    if up2bit[x] != dna[o+x]:
                        A[nextRowToAdd]           = A[row]
                        C[nextRowToAdd]           = C[row]
                        T[nextRowToAdd]           = T[row]
                        G[nextRowToAdd]           = G[row]
                        COUNT[nextRowToAdd]       = COUNT[row]
                        SEQ[nextRowToAdd]         = sum([up2bit[z]<<(2*t) for t,z in enumerate(xrange(x+1,len(up2bit)))],(1 << (len(up2bit)-x-1)*2))
                        A[row]                    = 0
                        C[row]                    = 0
                        T[row]                    = 0
                        G[row]                    = 0
                        (A,C,T,G)[up2bit[x]][row] = nextRowToAdd
                        (A,C,T,G)[dna[o+x]][row]  = nextRowToAdd+1
                        COUNT[row]               += count
                        SEQ[row]                  = sum([up2bit[z]<<(2*t) for t,z in enumerate(xrange(0,x))],(1 << x*2))
                        COUNT[nextRowToAdd+1]     = count                                                    ## This second new row already has 0,0,0,0 for its pipes, so we just have to set the
                        SEQ[nextRowToAdd+1]       = sum([dna[z]<<(2*t) for t,z in enumerate(xrange(o+x+1,len(dna)))],(1 << (len(dna)-o-x-1)*2))
                        nextRowToAdd             += 2
                        return
                A[nextRowToAdd]                     = A[row]
                C[nextRowToAdd]                     = C[row]
                T[nextRowToAdd]                     = T[row]
                G[nextRowToAdd]                     = G[row]
                COUNT[nextRowToAdd]                 = COUNT[row]
                SEQ[nextRowToAdd]                   = sum([up2bit[z]<<(2*t) for t,z in enumerate(xrange(len(dna)-o+1,len(up2bit)))],(1 << (len(up2bit)-1+o-len(dna))*2))
                A[row]                              = 0
                C[row]                              = 0
                T[row]                              = 0
                G[row]                              = 0
                COUNT[row]                         += count
                SEQ[row]                            = sum([up2bit[z]<<(2*t) for t,z in enumerate(xrange(0,len(dna)-o))],(1 << (len(dna)-o)*2))
                (A,C,T,G)[up2bit[len(dna)-o]][row]  = nextRowToAdd
                nextRowToAdd                       += 1
                return
            elif len(dna)-o > len(up2bit):                                                ## Type 2. The DNA in the fragment is longer than DNA in the row. This gives us again two possibilities:
                for x in xrange(0,len(up2bit)):
                    if dna[o+x] != up2bit[x]:
                        A[nextRowToAdd]           = A[row]
                        C[nextRowToAdd]           = C[row]
                        T[nextRowToAdd]           = T[row]
                        G[nextRowToAdd]           = G[row]
                        COUNT[nextRowToAdd]       = COUNT[row]
                        SEQ[nextRowToAdd]         = sum([up2bit[z]<<(2*t) for t,z in enumerate(xrange(x+1,len(up2bit)))],(1 << (len(up2bit)-x-1)*2))
                        A[row]                    = 0
                        C[row]                    = 0
                        T[row]                    = 0
                        G[row]                    = 0
                        (A,C,T,G)[up2bit[x]][row] = nextRowToAdd
                        (A,C,T,G)[dna[o+x]][row]  = nextRowToAdd+1
                        COUNT[row]               += count
                        SEQ[row]                  = sum([up2bit[z]<<(2*t) for t,z in enumerate(xrange(0,x))],(1 << x*2))
                        nextRowToAdd             += 1
                        o                        += x                                               ## We know we need to make a new row for the DNA the fragment had that the original row's
                        for y in xrange(0, len(dna)-o, 32):                                         ## Seq did not, but we don't know how long the DNA we need to add is. We might need to 
                            (A,C,T,G)[dna[o+y]][row] = nextRowToAdd                                 ## make several new rows to fit it all in, and so that is what this loop is doing.
                            COUNT[nextRowToAdd]      = count
                            SEQ[nextRowToAdd]        = sum([dna[o+z]<<(2*t) for t,z in enumerate(xrange(y+1,y+32 if y+32 < len(dna)-o else len(dna)-o))], (1<< (31 if 31 < len(dna)-o-1 else len(dna)-o-1)*2))
                            row                      = nextRowToAdd
                            nextRowToAdd            += 1
                        return
                COUNT[row] += count
                o          += len(up2bit)                                                           ## cut our fragment's DNA to be just the stuff the fragment has extra,
                temp        = (A,C,T,G)[dna[o]][row]                                                ## and check to see if this row has a warp pipe to where we want to go next.
                if temp == 0:                                                                       ## If there is no warp pipe...
                    for y in xrange(0, len(dna)-o, 32):                                             ## We create one (or as many as we need in a chain)
                        (A,C,T,G)[dna[o+y]][row] = nextRowToAdd
                        COUNT[nextRowToAdd]      = count
                        SEQ[nextRowToAdd]        = sum([dna[o+z]<<(2*t) for t,z in enumerate(xrange(y+1,y+32 if y+32 < len(dna)-o else len(dna)-o))], (1<< (31 if 31 < len(dna)-o-1 else len(dna)-o-1)*2))
                        row                      = nextRowToAdd
                        nextRowToAdd            += 1
                        return
                else:                                                                               ## But if there is a warp pipe,
                    row = temp                                                                      ## we just take it :)
                    o  += 1
                    continue
            else:                                                                   ## Type 3. The DNA in the fragment is the same length as the DNA in the row. This gives us again two possibilities:
                for x in xrange(0,len(up2bit)):
                    if dna[o+x] != up2bit[x]:
                        A[nextRowToAdd]           = A[row]
                        C[nextRowToAdd]           = C[row]
                        T[nextRowToAdd]           = T[row]
                        G[nextRowToAdd]           = G[row]
                        COUNT[nextRowToAdd]       = COUNT[row]
                        SEQ[nextRowToAdd]         = sum([up2bit[z]<<(2*t) for t,z in enumerate(xrange(x+1,len(up2bit)))],(1 << (len(up2bit)-x-1)*2))
                        A[row]                    = 0
                        C[row]                    = 0
                        T[row]                    = 0
                        G[row]                    = 0
                        (A,C,T,G)[up2bit[x]][row] = nextRowToAdd
                        (A,C,T,G)[dna[o+x]][row]  = nextRowToAdd+1
                        COUNT[row]               += count
                        SEQ[row]                  = sum([up2bit[z]<<(2*t) for t,z in enumerate(xrange(0,x))],(1 << x*2))
                        COUNT[nextRowToAdd+1]     = count
                        SEQ[nextRowToAdd+1]       = sum([dna[z]<<(2*t) for t,z in enumerate(xrange(o+x+1,len(dna)))],(1 << (len(dna)-o-x-1)*2))
                        nextRowToAdd             += 2                                                                   ## row, so we don't have to do the loop we did at the end of Type 2.
                        return
                COUNT[row] += count                                                   ## 1) Very simply, we just increment the count and we're done. 
                return


####################
## Bonus Features ##
##########################################################################################################
##                                                                                                      ##
## Below are some functions that are useful to have around, but not necessarily important or even used. ##
##                                                                                                      ##
##########################################################################################################

## Figures out how much RAM we are using and returns it in a human-readable way.
def getRAM():
    if   arrayKind == 'cffi':   num = ffi.sizeof(A)*7
    elif arrayKind == 'numpy':  num = A.nbytes*7
    elif arrayKind == 'ctypes': num = ctypes.sizeof(A)*7
    for unit in [' ',' K',' M',' G',' T']:
        if abs(num) >= 1024.0: num /= 1024.0
        else: print "   [ RAM used @ %3.1f%sb ]" % (num, unit); return

## These two functions are used when ACGTrie is called with --fragment:
def subfragment(DNA,count):
    for l in range(len(DNA)-1,-1,-1):
        seqChunks[DNA[l:]] += count

def emptyCache(add,nextRowToAdd):
    global seqChunks
    for fragment,count in sorted(seqChunks.items(), reverse=True, key=lambda t: len(t[0])):
        add(fragment,count)
    seqChunks = collections.defaultdict(int)

def growTrie():
    # Originally I wrote to disk then pulled it back because most methods to
    # extend C structured array requires having both the old and new array in
    # memory at once for some period of time, but since our columns are individual
    # arrays, we only require 28% of the RAM that doing it as a single array takes.
    # This is also why we extend SEQ first.
    global A
    global C
    global G
    global T
    global COUNT
    global SEQ

    if arrayKind == 'numpy':
        SEQ = numpy.concatenate((SEQ,numpy.zeros(10000000, dtype='int64')))
        A = numpy.concatenate((A,numpy.zeros(10000000, dtype='uint32')))
        C = numpy.concatenate((C,numpy.zeros(10000000, dtype='uint32')))
        G = numpy.concatenate((G,numpy.zeros(10000000, dtype='uint32')))
        T = numpy.concatenate((T,numpy.zeros(10000000, dtype='uint32')))
        COUNT = numpy.concatenate((COUNT,numpy.zeros(10000000, dtype='uint32')))
    elif arrayKind == 'ctypes':
        newSize = A._length_+10000000
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
            tempA     = A;     A     = ffi.new("uint32_t[]", len(A)+10000000);     ffi.memmove(A,tempA,len(tempA)*(ffi.sizeof(A)/len(A)))
            tempC     = C;     C     = ffi.new("uint32_t[]", len(C)+10000000);     ffi.memmove(C,tempC,len(tempC)*(ffi.sizeof(C)/len(C)))
            tempG     = G;     G     = ffi.new("uint32_t[]", len(G)+10000000);     ffi.memmove(G,tempG,len(tempG)*(ffi.sizeof(G)/len(G)))
            tempT     = T;     T     = ffi.new("uint32_t[]", len(T)+10000000);     ffi.memmove(T,tempT,len(tempT)*(ffi.sizeof(T)/len(T)))
            tempCOUNT = COUNT; COUNT = ffi.new("uint32_t[]", len(COUNT)+10000000); ffi.memmove(COUNT,tempCOUNT,len(tempCOUNT)*(ffi.sizeof(COUNT)/len(COUNT)))
            tempSEQ   = SEQ;   SEQ   = ffi.new("int64_t[]", len(SEQ)+10000000);    ffi.memmove(SEQ,tempSEQ,len(tempSEQ)*(ffi.sizeof(SEQ)/len(SEQ)))
        except AttributeError:
            print 'ERROR! Your version of python does not fully support cffi! Please upgrade, or set the intial amount of RAM to something larger.'
            exit()
    getRAM()

class fileStats:
    def __init__(self,DNA,count):
        self.start = time.time()
        self.linesRead = count
        self.fragmentAvg = len(DNA)*count
    def add(self,DNA,count):
        self.linesRead += count
        self.fragmentAvg += len(DNA)*count
    def result(self):
        return (self.linesRead,int(self.fragmentAvg/self.linesRead),time.time()-self.start)

######################
## Create the table ##
##########################################################################################################
##                                                                                                      ##
## There are many ways to create an in-memory table. Some are faster than others. Some are part of the  ##
## standard library and don't require any other modules to be downloaded. Some only work on the latest  ##
## versions of python. But after we have created the table, we read/write data in exactly the same way  ##
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
countOverflow = {}
warpOverflow = {}
getRAM()

#########################
## Start reading stdin ##
##########################################################################################################
##                                                                                                      ##
## First we need to figure out what format the input data is coming in as - just DNA with newlines, or  ##
## in CSV format with counts too. Once we've done that, we need to check the execution parameters, and  ##
## enter the appropriate stdin read loop, with the right adding function.                               ##
##                                                                                                      ##
##########################################################################################################

startTime = str(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

# Read first row to see if DNA comes with counts. Create appropriate stdin generator.
firstFragment = sys.stdin.readline().strip().split(',')
if   len(firstFragment) == 1:
    stdin = ((line.rstrip(),1) for line in sys.stdin)
    firstFragment.append(1)
elif len(firstFragment) == 2:
    stdin = ( (DNA,int(count)) for DNA,count in csv.reader(sys.stdin, delimiter=',') )
    firstFragment[1] = int(firstFragment[1])
else:
    print 'ERROR: I do not understand this kind of stdin format :U'; exit()

#What kind of adding function to use?
if args.walk or args.fragment: add = addRowWalk
else:                          add = None #addRow function not yet written.

stats = fileStats(firstFragment[0],firstFragment[1])

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

## Else, we just go straight into it.
else:
    add(*firstFragment)
    for DNA,count in stdin:
        stats.add(DNA,count)
        add(DNA,count)
        if nextRowToAdd + 100 > len(A): growTrie()

linesRead, fragmentAvg, duration = stats.result()

print 'Done in: ', duration
print 'Lines read: ', linesRead
print 'Average fragment length: ', fragmentAvg
print sum(A),sum(C),sum(T),sum(G),sum(COUNT),sum(SEQ)

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
    'fragments': linesRead,
    'fragmentAvgLen': fragmentAvg,
    'rows': nextRowToAdd,
    'analysisTime': startTime,
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

if sys.byteorder == 'big':
    print '   [ Flipping eggs. ]'; ## Little Endians 4 lyfe yo.
    if arrayKind == 'numpy': A.byteswap(True);C.byteswap(True);G.byteswap(True);T.byteswap(True);COUNT.byteswap(True);SEQ.byteswap(True)
    else:
        print 'Support for byteswap in cffi is comming soon! Until then, youll have to swap manually in numpy on loading the data on another machine.'

print len(A)
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

'''
Determine 
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
