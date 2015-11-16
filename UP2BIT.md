## up2bit format 

Here we define the **up2bit** format - a binary encoding used to store *variable-length* DNA in a *fixed* number of bits. 
Basically, we use it to turn a string of DNA (of any length) into a number, which can then be stored very efficiently.

Although ACGTrie wouldn't be possible without the up2bit format, in truth the implimentation specifics are not something most 
Biologists will care about, so it is probably safe to skip over it entirely ;)

If you wish to keep reading, *great*! But brace yourself... winter is coming and you are about to go over the wall :T

As you may have seen in the *OUTPUT_FORMAT* docs, we have a column called **SEQ** which contains DNA like "ACGT".
In the actual table that is saved on your hard drive, we can't store little As, Cs, Ts, and Gs, because computers only 
understand 1s and 0s. We would need like, 2s and 3s as well, and hard drives dont support 2s and 3s. Pretty sure they dont anyways.
To be honest, I have never really looked... 

So whatever, we have to encode our DNA in 1s and 0s somehow.
One mechanism is to use **ASCII**. You might have heard of ASCII. It's a pretty big deal.
The idea is that you always use 8 bits of binary for every letter you want to store, a bit like how DNA coding for amino acids 
always use 3 'bits' of DNA. 8 bits of binary gives you 256 possible combinations. *00000000*, *000000001*, *00000010* ... *11111111*.
Thats enough combinations to store all the letters on a Latin keyboard, and a few extra letters and symbols too. 
Which combination of 1s and 0s maps to which letter/symbol is set out by the ASCII encoding standard. 
So *A* in ASCII is "01000001", *T* is '01010100', *a* is '01100001', *@* is '01000000', etc.

Now i'm not one to size-shame large encodings - data is beautiful at any size - but if we only need to store four letters: 
*A*, *C*, *T* and *G*, it is kind of wasteful to use 8bits when we really only need 2bits - *00*, *01*, *10*, and *11*.

This format actually already exists and is called **2bit**. Maybe you heard of it when you downloaded a genome fasta file.
It's basically the same as ASCII but 4x smaller by storing DNA sequences in this compact encoding.
But 2bit has a problem - we can't use it when the space to store our data is fixed, like our 64bit ACGTrie table columns.
No data in a 32bit column might be *00000000000000000000000000000000*, but in 2bit that would mean *AAAAAAAAAAAAAAAA*.
It's like our row is literally screaming at us because it has no way to encode silence. Poor row. So we made a format which 
supports a variable amount of DNA info in a fixed space called up2bit. Get it? While 2bit can encode *32* letters of 
DNA in *64bits*, it cannot encode *31* or *30* or even *0* letters of DNA. But up2bit supports **up too** 31 letters in 64bits.
You sacrifice 1 letter for the ability to store a variable amount of letters. 

It's very simple really - reading **right to left** (which is also how binary is counted) just encode all your letters in 2 bits,
then stick a *01* on the left end (this is also a C, but that is a coincidence). Any unused bits in your  fixed space can then just be
left as 0s.

So *ACTG* in up2bit format would look like:

                                               -  G  T  C  A
                                              01 11 10 01 00

Which when stored in a 32bit space, looks like:

             00 00 00 00 00 00 00 00 00 00 00 01 11 10 01 00

Here's another example, 15 Ts:                                                                       

             01 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10

Whilst no data at all would look like:                                                               

             00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 01

Which is just the number **1** in binary. If you see in the code if statements asking *if SEQ == 1*, 
this is why; we are actually asking if it's empty or not.

And that's one of the great things about the up2bit format - because we read right to left and all unused spaces are added to 
the left, whatever fixed space you store up2bit data in, it wont change the value that the binary represents. *ACTG* in up2bit is **always** 
 *228*, whether you store it in 8bits, 16bits, 32bits, 64bits, or even 3234983488bits; and no data is always *1*.
You now know pretty much everything needed to impliment and decode the up2bit format, but here are some Frequently Asked Questions for the curious:

*"Why use A-00, C-01, T-10, G-11? It's not in alphabetical order, nor is it the order of the letters of your website!!?"*

Well firstly, it is the order of the letters for the website (check http://ac.tg and see!). Secondly, it's because  
 there's a neat trick you can do to convert those four letters of ASCII directly into up2bit values. 
 Here are the binary encodings for A, C, G, and T in ASCII:                                                 

           -vv-                                                                                         
    A - 1000001                                                                                         
    C - 1000011                                                                                         
    G - 1000111                                                                                         
    T - 1010100                                                                                         

See the 2nd and 3rd from the right values? They just happen to be *00*, *01*, *11*, and *10*. In some programming languages like
C, we can cut out just those numbers very quickly. Unfortunately, in python, the code to do that, which is:

    ord('T')>>1 &3

Is actually slower than just looking up the value in a dictionary/list. Also, unless you know what you're looking at, it would be
confusing as hell to anyone who didnt know the trick. For this reason, we dont use it in our code, but we stuck with the format 
to maybe help another brother out who might use up2bit files and want's to make use of this trick in their own code.

*"Why use '01' as the capping bits? Why not '11' or '10'?"*


Well sometimes, unfortunately, computers need to store negative numbers. How do you store a negative number when you only have 
1s and 0s? Well, you say in the code *'these numbers can be negative, and if they are I will make the very first bit a 1'*. 
This means you can only count *half as high*, as the max number in say, 8bit, would be 01111111, but you can also count all 
those numbers in negative too, so the range is still 256, but the max number is 128.

So when we define an encoding, say 8bit, 32bit, 64bit, etc, we first have to say if it will be "signed" (can have negative 
numbers, first bit determines +/-), or "unsigned", meaning they can only count 0 or higher (but twice as high as signed). 
Sometimes unsigned numbers will also tell you about this cool indie band they're in, but it's best to just ignore them.

So why do we care? Well, in python (and other high-level programming languages), maths is generally done using signed 64bit 
numbers. Therefore, if we want to store a python number in a NumPy array that is anything other than signed 64bit, we're going 
to have to do some checking to make sure the number will fit (not too high, or negative, etc), and this check takes time. 
However, if we are just stuffing native signed 64bit numbers into a signed 64bit space, we don't need to do the checks. 
Thus, NumPy is alot faster working with signed 64bit arrays than it is with unsigned 64bit arrays (or anything else).
If the cap is **01**, we can *never have a negative number*.
No negative numbers means we do not have weird issues like adding 1 to a negative number and the number getting *smaller*.
It also means that we do not need to do any checks when copying up2bit numbers from unsigned spaces to signed spaces, because 
they will always fit, and always represent the exact same number. win-win :)

*"If ACGTrie uses 64bits to store DNA, which means 1bp encoded in the warp pipe used to get to that row, and 31bp of DNA in the Seq column, what happens if the DNA to be added is more than 32 letters long?"*

This one is easy - you just simply cut your sequence at the 32nd base, the 33rd base is encoded in the next warp pipe, and the 
extra N bases are added to the next rows seq, and so on. So a 99bp fragment would need at least 4 rows to store all its data.
However that is only true for 64bit columns. There's nothing stopping you from editing the code and using 128bits or even 
higher for the seq column - then you could store a lot more DNA before we have to break and go to a new row. While 64bit columns
are good for fragments around 50bp in length, I highly recommend increasing this to something higher if your fragments are longer
than 50bp. Really, it depends on your data. If you do change the array size, everything else, as Steve Jobs would say, should "just work", since the file headers contain the array sizes, and all programs that we write to parse ACGTrie data is compatible with arrays of variable size :)
                                                                                                      
 This concludes all the information I can provide about the up2bit encoding. If you can think of an   
 improvement to the format let me know, and we can break backwards compatibility together! \o/
