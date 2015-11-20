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
## As we add new DNA to the trie, we have to check if the new DNA fits into an existing row, or needs   ##
## to become a new row of it's own. The following functions are used in the creation of the trie. Note  ##
## that this script only serves to explain the algorithum we use, and is not gaurenteed to be quick!!   ##
##                                                                                                      ##
##########################################################################################################

## We use this function to turn the up2bit number in the SEQ column into a list of 2bit values (to compare to the list of 2bit numbers of the new DNA)
def up2bit_list(value):
    value = int(value)
    up2bits = [((value >> x) & 3) for x in range(0,64,2)]
    return up2bits[:-up2bits[::-1].index(1)-1]                 ## A fancy way of saying 'reverse the list, and cut off the '01' cap and everything before it'.

## We start by using this function to see if the existing/new DNA lists matches or not:
def firstNonMatching(longerDNA,shorterDNA):
    for x,y in enumerate(longerDNA):
        if y != shorterDNA[x]:
            return  x           ## If they dont match, x is the base where the mismatch happens
    else:   return -1           ## If there is no DNA (because we used it up while traveling the warp to get to this row), x is -1
    return x                    ## If they match, x is the length of the shorter DNA.


## If it isn't a perfect match, we will have to split the existing row into 2 or 3 rows (1 or 2 new rows).
## We make 1 row if the new fragment matches but isn't the full length of the existing row's SEQ, 
## and we make 2 rows when the new fragment matches somewhat the existing row, but also contains some unique DNA of its own.

## If we have to split, we start by copying the data that is unique to the existing row to a new row:
def copyRow(row,nextRowToAdd,up2bit,x):
    A[nextRowToAdd] = A[row]
    C[nextRowToAdd] = C[row]
    T[nextRowToAdd] = T[row]
    G[nextRowToAdd] = G[row]
    COUNT[nextRowToAdd] = COUNT[row]
    SEQ[nextRowToAdd] = list_up2bit(up2bit[x+1:len(up2bit)])

## We then edit the original row to have either 1 or 2 warps to pipes pointing to the 1 or 2 new rows.
## The COUNT is increased of course, and the SEQ becomes whatever came before the mismatch.
def freshRow(firstPipe,secondPipe,row,count,DNAJ,x):
    global nextRowToAdd
    A[row] = 0
    C[row] = 0
    T[row] = 0
    G[row] = 0
    COUNT[row] += count
    SEQ[row] =  list_up2bit(DNAJ[:x])
    (A,C,T,G)[firstPipe][row] = nextRowToAdd       ## first pipe
    nextRowToAdd += 1
    if secondPipe != None:
        (A,C,T,G)[secondPipe][row] = nextRowToAdd  ## second pipe

## And finally, if needed, make a second row (DNA unique to the incoming fragment):
def secondRow(count,seq,x):
    global nextRowToAdd
    COUNT[nextRowToAdd] = count
    SEQ[nextRowToAdd]   = list_up2bit(seq[x+1:len(seq)])
    nextRowToAdd += 1

## Sometimes we have DNA to add that is longer than what the SEQ column can hold. This happens irrespective of
## splitting or just adding from nothing. When this happens, we keep adding as many new rows as it takes.
def restOfSeq(seq,row,count):
    global nextRowToAdd
    for y in xrange(0, len(seq), 32):
        (A,C,T,G)[seq[y]][row]    = nextRowToAdd
        COUNT[nextRowToAdd]       = count
        SEQ[nextRowToAdd]         = list_up2bit(seq[y+1:y+32])
        row                       = nextRowToAdd
        nextRowToAdd             += 1

## Oh, and to turn our lists back to up2bit numbers, we use this function:
def list_up2bit(twobits):
    up2bit = 1                             ## The '01' cap.
    for twobit in reversed(twobits):       ## Kind of similar to string_up2bit except we already start with a list of
        up2bit = (up2bit << 2) + twobit    ## 2bit numbers and we just need to cap it and stack them all together.
    return up2bit

## Finally, the main logic of the whole program - how to use the above to add data to the trie:
def addRowWalk(DNA,count):
    global nextRowToAdd
    seq = [('A','C','T','G').index(char) for char in DNA]                           ## Turn DNA into a list of 2bit numbers.
    row = 0                                                                         ## We always start on row 0.
    while True:
        if SEQ[row] == 1:                                                           ## There is no sequence data in this row, and thus 3 possible options going forward.
            COUNT[row] += count                                                     ## ( But all of them start by adding +1 to this rows's # count )
            if not seq: return
            warp_pipe = (A,C,T,G)[seq[0]][row]                                  ## Getting the value for thisRow[seq[0]] takes time. We only want to get it once.
            if warp_pipe:                                                       ## Type 2. Happens 97% of the time, which is why its the first thing we try.
                row = warp_pipe
                seq = seq[1:]                                                   ## Note the seq becomes a bit shorter, because we need to lose a base following the pipe.
                continue
            else: restOfSeq(seq,row,count)
        else:                                                                       ## There IS sequence data in this row, now we have three possible options going forward:
            up2bit = up2bit_list(SEQ[row])
            if len(seq) < len(up2bit):                                                  ## Type 1. DNA in the row is longer than DNA in the fragment. This gives us two possibilities:
                if up2bit[:len(seq)] == seq:                                                              ## 1) The Fragment's DNA, although shorter than the row's DNA, matches the row.
                    copyRow(row,nextRowToAdd,up2bit,len(seq))
                    freshRow(up2bit[len(seq)],None,row,count,up2bit,len(seq))
                else:                                                                                   ## 2) The Fragment's DNA and the row's DNA do NOT match. So now we have to make 2 new rows...
                    x = firstNonMatching(seq,up2bit)
                    copyRow(row,nextRowToAdd,up2bit,x)
                    freshRow(up2bit[x],seq[x],row,count,up2bit,x)
                    secondRow(count,seq,x)
            elif len(seq) > len(up2bit):                                                ## Type 2. The DNA in the fragment is longer than DNA in the row. This gives us again two possibilities:
                if seq[:len(up2bit)] == up2bit:                                                         ## 1) The row's DNA, although shorter than the fragments's DNA, matches the fragment.
                    COUNT[row] += count
                    seq         = seq[len(up2bit):]                                                             ## cut our fragment's DNA to be just the stuff the fragment has extra,
                    warp_pipe   = (A,C,T,G)[seq[0]][row]                                                  ## and check to see if this row has a warp pipe to where we want to go next.
                    if not warp_pipe: restOfSeq(seq,row,count)
                    else:                                                                               ## But if there is a warp pipe,
                        row = warp_pipe                                                                 ## we just take it :)
                        seq = seq[1:]
                        continue
                else:                                                                                   ## 2) The rows's DNA and the fragments's DNA do NOT match. So now we have to make 2 new rows...
                    x = firstNonMatching(up2bit,seq)
                    copyRow(row,nextRowToAdd,up2bit,x)
                    freshRow(up2bit[x],seq[x],row,count,up2bit,x)
                    restOfSeq(seq[x:],row,count)
            else:                                                                   ## Type 3. The DNA in the fragment is the same length as the DNA in the row. This gives us again two possibilities:
                if up2bit == seq:
                    COUNT[row] += count                                                   ## 1) Very simply, we just increment the count and we're done. 
                else:
                    x = firstNonMatching(seq,up2bit)
                    copyRow(row,nextRowToAdd,up2bit,x)
                    freshRow(up2bit[x],seq[x],row,count,up2bit,x)
                    secondRow(count,seq,x)
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
