## ACGTrie
 Hi! In this doc we explain what ACGTrie is doing in a bit more depth. To truly understand how the    
 code works, you need to look at the code itself. We provide two versions, ACGTrie_FAST, and          
 ACGTrie_LEARN. Both work in the same way and both get the same result, it is simply that LEARN is    
 easy to read python, while FAST is highly optimized for performance :)

So it's pretty simple really; if we had a single sequence like 'ACGACCC', and we want to ask some questions about it's "DNA composition", we might be interested in a table like the one below:

                                               row    Seq      Count                                              
                                                0     C        4                                                  
                                                1     A        2                                                  
                                                2     G        1                                                  
                                                3     AC       2                                                  
                                                4     CC       2                                                  
                                                5     CG       1                                                  
                                                6     GA       1                                                  
                                                7     ACC      1                                                  
                                                8     ACG      1                                                  
                                                9     CCC      1                                                  
                                                10    CGA      1                                                  
                                                11    GAC      1                                                  
                                                12    CGAC     1                                                  
                                                13    ACCC     1                                                  
                                                14    ACGA     1                                                  
                                                63    GACC     1                                                  
                                                32    CGACC    1                                                  
                                                17    ACGAC    1                                                  
                                                18    GACCC    1                                                  
                                                19    CGACCC   1                                                  
                                                20    ACGACC   1                                                  
                                                21    ACGACCC  1                                                  

Now we can ask questions like 'how prevelent is motif X?' or 'how much seq data covers repeats?', or, if we have two of these  tables, 'what are the key differences between sample A and sample B?'. But at the end of the day it's really nothing more than a table of counts for all the different subfragments of DNA within our sample. Note that 'CC' has a count of 2 due to the 'CCC, but you get the basic idea...

Problem is, if we tried to make the table above for some 5Gb+ seq data, there would be **so many** combinations of DNA that we would run out of memory/RAM before we processed even 1% of chromosome 1. If we tried to build the table directly onto the hard drive, it would take a **very long time** because compared to RAM, disk access is really slow. 

Bioinformaticians have tried to solve this problem with two methods: only storing subfragments of a specific size, and approximating the DNA composition (instead of getting exact results) using a variety of ever more ingenious (and complex) computer science algorithums. At ACGT we took a different approache, to keep our algorithums as simple as possible so that even Biologists with a passing interest in how their analysis tool works can "get it" it in just a few minutes, and to make an output file that is really easy to work with. No need to visit Wikipedia for entries on Bloom filters or complex graph theory - everything is right here the docs and code :)

*(also, our tools are pretty darn fast)*

So in an attempt to make a small and simple, yet comprehensive output format, this is what we have come up with:
                                                                                                                
                                          row  A   C   T   G COUNT SEQ                                           
                                         ----------------------------------                                         
                                           0   1   2   0   3   7                                                   
                                           1   0   5   0   4   2   C                                             
                                           2   0   7   0   6   4                                                   
                                           3   0   0   0   0   1   ACCC                                          
                                           4   0   0   0   0   1   ACCC                                          
                                           5   0   0   0   0   1   C                                             
                                           6   0   0   0   0   1   ACCC                                          
                                           7   0   8   0   0   2                                                 
                                           8   0   0   0   0   1                                                 
                                                                                                     
 So it's just a table, with 6 columns: **A**, **C**, **T**, **G**, **COUNT**, and **SEQ** ("row" is just the row number, it's not actually stored anywhere). To find the count of some fragment of DNA, you start on row 0, and use one of the **A** **C** **T** **G** *"warp pipe"* columns to go to a new row. If you have ever read one of those *create your own adventure* books, you'll recognize this format immediately :)
 
 The **COUNT** column tells us how many times the fragment that this row represents was seen in the data. The **SEQ** column allows us to group multiple fragments together into 1 row. For example, if there were many fragments that started with AAAAA, but only one of which ended in TCGA, that different ending could be stored in a single row, rather than 4 rows (one for each base). 
 
It might be confusing at first, particularly if you forget that the row itself represents a base due to the warp pipe used to  get to it, but once you "get it" you'll see it's really really simple. This data structure is refered to in Computer Science as a 'trie', hence "ACGTrie" :)
 
So say you were interested in counting how many "A"s there were in the data. You'd start on row 0, and from there goto the **A** warp pipe, which has a *1* in it. This means to go to row 1. Then from row 1, which represents "A", we can see in the **COUNT** column there is a *2*. This means there are 2 As in the input data. There's also some **SEQ** data in this row, a *C*, which means that if we wanted to look up "AC" we would also end up here, and that would also be a 2. You can confirm this on the original table up top. In fact it's worth taking a minute or two to try out some of the subfragments in the original table, and seeing if you can get the same result in the ACGTrie format by-hand :)

So as you can already see, our ACGTrie table takes less data to store the same information as in the larger table at the top of this doc - but where this format really shines is when you have a **lot** of data. My estimates put it at 170x smaller than the big table for real-world data, and much higher when your input data has very long strings of DNA.
170x smaller data means 170x less RAM needed when constructing the table, which is how we can create a **full** DNA composition table in **less** space than a typical k-mer tool can.
Furthermore, while a fixed-length k-mer tool can build a DNA composition table for a given k-mer size much faster than ACGTrie can for all lengths of DNA, if it had to be run 20 times to get just fragments up to length 20, ACGTrie would have already finished and would include fragments of theoretically infinite length :)

To see examples of how to read ACGTrie tables from within python (and hopefully other languages soon as people contribute to the project!) check the samples directory for example code. To see how specifically the ACGTrie tables are built, check the ACGTrie_LEARN.py code - it's well commented I promise ;) 
