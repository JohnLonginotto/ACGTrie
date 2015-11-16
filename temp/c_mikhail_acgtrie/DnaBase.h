#pragma once

#ifdef WINDOWS
#  include "debug_utils.h"
#else
#  define FUNC_GUARD
#endif

// TODO? to use masks arrays for 3 << (charInd * 2) and binary 111111s ((c_one << (charCount * 2)) - 1)

class CSequenceUp2Bit
{
public:
	typedef unsigned long long TSequenceUp2BitStorage;
	static const int c_maxLen = sizeof(TSequenceUp2BitStorage) * 4 - 1;

	TSequenceUp2BitStorage up2BitCode;

	CSequenceUp2Bit()
		: up2BitCode(1)
	{	}

	int GetLength() const;     // This method is optimized, but not very fast
	int GetChar2Bits(int charInd) const
	{
		assert(charInd >= 0 && charInd < GetLength());
		return (up2BitCode >> (charInd * 2)) & 3;
	}
	CSequenceUp2Bit GetSubsequence(int startCharInd, int charCount) const;
	//TSequenceUp2BitStorage GetBits // Without first 01

private:
	static const TSequenceUp2BitStorage c_one = 1;
};


// DNA characters are packed into 2-bit arrays. 
// Let's try to pack them into 8-byte integers, as 01-prefixed up2bit.
// In big bits cut operations this will give simpler code and maybe a very little slowdown
// on 32-bit processors. Simple operations like getting 2 bits from the array
// we'll optimize separately
typedef CSequenceUp2Bit::TSequenceUp2BitStorage TDna2BitsPortion;        
//// DNA characters are packed into 2-bit arrays. One group was made 4 byte (16 characters)
//// despite the fact that CSequenceUp2Bit was made 8-byte.
//// This was done in order not to do superfluous work when processing 64-bit integers. 
//// I suspect 64-bit integers processing can be slower than 32-bit's even on 64-bit processors
//// and is apparently slower on 32-bit ones
//typedef unsigned int TDna2BitsPortion;               
static const int c_dnaCharsInPortion = sizeof(TDna2BitsPortion) * 4;

class CDna2Bits
{
public:
	CDna2Bits()
		: m_len(0)
	{	}

	void AssignFromString(const std::string &dnaStr);
	int GetLength() const      {  return m_len;  }
	int GetChar2Bits(int charInd) const
	{
		assert(charInd >= 0 && charInd < GetLength());
		return (m_bits[charInd / c_dnaCharsInPortion] >> 
						 ((charInd % c_dnaCharsInPortion) * 2)) & 3;
	}
	const TDna2BitsPortion *GetBits() const;  //-?
	void GetSubsequenceUp2Bit(CSequenceUp2Bit &seqUp2Bit, int startCharInd, int charCount) const;
	int GetEqualCharCount(CSequenceUp2Bit &seqUp2Bit, int startCharInd) const;
	            // Returns number of equal DNA characters between DNA's subsequence and 
					    // the given sequence
	//bool AreSubsequencesEqual(const CSequenceUp2Bit &seqUp2Bit, 
	//	                    int thisStartCharInd, // int seqStartCharInd, 
	//											int charCount) const;

	static void RunUnitTests();

private:
	static const TDna2BitsPortion c_one = 1;
	
	typedef std::vector<TDna2BitsPortion> TDna2Bits;

	TDna2Bits m_bits;
	int m_len;   // Number of 2-bit DNA characters

	static int GetFirstNonZeroCharacter(TDna2BitsPortion codesDiff);
};


void DnaBase_RunUnitTests();
