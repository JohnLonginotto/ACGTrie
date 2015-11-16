#!/usr/bin/env python
import hts
import time
import argparse
import subprocess
import collections

## Parse user-supplied command line options:
parser = argparse.ArgumentParser(                              description="Put in a BAM file or fragments of DNA, get out an ACGTrie table.")
parser.add_argument("-o", "--output", metavar='/path/to/output.trie', help='Required. Output filename.')
parser.add_argument("-i", "--input",  metavar='/path/to/file.bam',    help='Required for SAM/BAM analysis. If no input file is provided, newline-sepurated DNA can be taken via stdin (called MANUAL mode, see code for more info...)')
parser.add_argument("--cpu",          metavar='1',default=1,          help="Optional. Number of processes/cores you want to use.")
parser.add_argument("--acgtrie",      default='ACGTrie',	          help="Optional. Path to ACGTrie if it cant be found automatically.")
parser.add_argument("-q", '--quiet',  action='store_true',            help='Optional. No status bars in output. Good for logs.')
args = parser.parse_args()

if args.input == None or args.output == None: print '''
	Oops.
    You need to provide a path/filename for your inputs and outputs!
    E.g. ./ACGTrie_BAM.py --input myInput.bam --output myOutput
'''; exit()

###############################
## Pre-processor for SAM/BAM ##
##########################################################################################################
##                                                                                                      ##
## So really that is the end of ACGTrie. The rest is just an example pre-processor for BAM files that   ##
## just so happens to be in the ACGTrie code aswell. In the future I will probably move it out and have ##
## a list of pre-processors for different file formats                                                  ##
## So as I mentioned before, you can do a lot of things with a smart pre-processor. Here is a demo of a ##
## pre-processor for SAM/BAM files. I've tried to get it as optimzed as possible, so its certainly not  ##
## the simplest example, but it sure is efficient :)                                                    ##
##                                                                                                      ##
##########################################################################################################

inputData = hts.Bam(args.input)
seqChunks = collections.OrderedDict()
# seqChunks = collections.defaultdict(int)
lostChildren = collections.defaultdict(int)
command = "'" + args.acgtrie + "' --rows 27844500 --walk --output " + args.output + '.64.AZ'
#command = "'" + args.acgtrie + "' --rows 27844500 "
print command
firstSubprocess = subprocess.Popen(command, stdin=subprocess.PIPE, shell=True, executable='/bin/bash')
for totalReads,line in enumerate(inputData):
    if totalReads == 1000000: break
    seq = line.seq
    if 'A' not in seq or 'N' in seq: continue # may want to split on N and treat as more than 1 read, or just throw away subfrags with N...?
    idx = -1
    while True:
        idx = seq.find('A',idx+1)
        if idx == -1: break
        try: seqChunks[seq[idx+5:idx+25]] += 1
        except KeyError: seqChunks[seq[idx+5:idx+25]] = 1
        # seqChunks[seq[idx+5:]] += 1
    if len(seqChunks) > 100000:
        #print totalReads, len(seqChunks)
        if '' in seqChunks:	
            lostChildren['A'] += seqChunks['']
            del seqChunks['']

        # Type 1 - Takes the oldest 50000 lines, sorts by len. 
        temp = {}
        while len(seqChunks) > 50000:
            k,v = seqChunks.popitem(last=False)
            temp[k] = v
        for fragment,count in sorted(temp.items(), reverse=True):
            firstSubprocess.stdin.write(fragment + ',' + str(count) + '\n') 

        ## Type 2 - Swap to this if it turns out Type 1 is too slow, because there isn't much in it.
        # x = 0
        # temp2 = []
        # for fragment,count in sorted(seqChunks.items(), reverse=True):
        #     x+=1
        #     temp2.append(fragment)
        #     firstSubprocess.stdin.write(fragment + ',' + str(count) + '\n')
        #     if x == 30000: # leaving 20000 in the dict
        #         for t in temp2: del seqChunks[t]
        #         break

## Finish up:
if '' in seqChunks: lostChildren['A'] += seqChunks['']; del seqChunks['']
'''
for fragment,count in sorted(seqChunks.items(), reverse=False, key=lambda t: len(t[0])):  ## smallest.
for fragment,count in sorted(seqChunks.items(), reverse=True): ## fastest
'''
for fragment,count in sorted(seqChunks.items(), reverse=True):
    firstSubprocess.stdin.write(fragment + ',' + str(count) + '\n') 

firstSubprocess.stdin.close()
firstSubprocess.wait()

#print 'lost children: ', lostChildren

'''
Count the occurence of each branch and re-arrange to suit!
--interpreter /path/to/pypy/or/cython   (must have hts-python not pysam. Maybe would work with samtools only.)
fragments up to Xbp in length - then uint32?
- How many processors?
- So each processor can have a table up to X bytes of RAM. This determines number of rows in numpy.zeros() (unless we can find out if we have MTU or empty())
- So need to process data no larger than X rows per cpu.
- MD5 file and get estimate for abundance of 5mers, 4mers, 3mers & average/max fragment length:
                    def MD5andCount(inputFile,samtools):
                        readCounting = subprocess.Popen('"' + samtools + '" view -c -', stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='/bin/bash')
                        md5 = hashlib.md5()
                        chunkSize = 128*md5.block_size
                        ping.change('#',os.path.getsize(inputFile))
                        with open(inputFile,'rb') as f: 
                            for x,chunk in enumerate(iter(lambda: f.read(chunkSize), b'')): 
                                md5.update(chunk)     # Read the file once, but 
                                readCounting.stdin.write(chunk)
                                ping.pong(x*chunkSize)
                        readCounting.stdin.close()
                        readCounting.wait()
                        return [int(readCounting.stdout.read()), md5.hexdigest()]
Use info on xmers to decide how to group 5/4/3mers chunks that will fit into RAM (or if RAM will be an issue?)
Estimate 4mers from 5mers?
Estimate 3mers from 4 mers?
Use some xmer plus avg/max fragment length to rows converter

# Preprocessor has to add: 
    # Responsible & Contact
    # New Structs?
    # sampleHash
    # sampleName
    # samplePath
    # sampleSize
    # sampleReads?
    # sampleCreated
    # analysisFilters
    # analysisCreated
    ## AND MORE - check rawSeQL's makeHeader, and SeQCs data for BAMs
'''

'''
Get --input to work, spawning processes to do 256 branches (first four levels)
Root row is the sequence of the branch.
For values that exceed # dtype, # == 0, signalling to check header for actual count.
Code to combine multiple tries made with --input.
Also try reading the struct data in JS land. Maybe this is where HDF5 will come in handy? Maybe numpy arrays for javascript exist.
Speed test against other Trie-making programs. Perhaps have a speed and RAM competition.
Consider using 5x array.array - one for each column - and pypy for fast access.
Also pytables instead of numpy?
Consider having trie trim all DNA to 20 characters after fragmentation of full Xbp read, perhaps top of addRow/Walk.
Consider 6th column for "trie depth". Need only be 8 bits. Also consider # column being 32 bits, and Seq being 64.
Get addRow to work.
'''

# --memory bytes ?

'''
orderedDict - try using the last 50%, and try delete/inserting after every +=1
Try cutting all substrings to be no longer than 31bp, or even 25. (I cant think of a use of >30bp)
Speedtest against existing trie-software.
Read 1024th of the input file, in this process, using all possible memory for trie.
From result, estimate total time for completion, (time*1024)/(cpu).
Work out amount of memory needed for a 1024th (based on number of rows before seq is 0 not 3)
Work out how to get as few readthroughs as possible, with the --cpus and ram, 
processing 1024ths, or 256ths, or 64ths, or 32ths, or 4ers, or just 1 cpu.
--level (don't do this check, just jump straight in to splitting at level X)

Better to use more CPUs and than less memory
Better to do less read-throughs than less memory

mores cpus > less readthroughs >= less memory
'''
