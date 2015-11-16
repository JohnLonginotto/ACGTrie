#include "StdAfx.h"
#include "DnaTrieBuilder.h"

int CSequenceUp2Bit::GetLength() const
{
	int zeroBitCount = 0;
	unsigned int curBlock;

	assert(sizeof(TSequenceUp2BitStorage) == 8);             // This method implemented for 8-byte codes
	assert((up2BitCode >> 32) == ((int *)&up2BitCode)[1]);   // Sometimes >> 32 works incorrectly
	                                                         // We imply that machine architecture is Little Endian
	if ((up2BitCode >> 32) == 0)
	{
		zeroBitCount = 32;
		curBlock = unsigned int(up2BitCode);
	}
	else 
		curBlock = up2BitCode >> 32;

	for (int mid = 16; mid != 1; mid /= 2)
	{
		if ((curBlock >> mid) == 0)
		{
			zeroBitCount += mid;
		//	curBlock &= (1U << mid) - 1;
		}
		else 
			curBlock >>= mid;
	}
	assert(curBlock == 1);
	assert((up2BitCode >> (sizeof(TSequenceUp2BitStorage) * 8 - zeroBitCount)) == 0);
	assert((up2BitCode >> (sizeof(TSequenceUp2BitStorage) * 8 - zeroBitCount - 2)) == 1);
	return (sizeof(TSequenceUp2BitStorage) * 8 - zeroBitCount) / 2 - 1;
}

CSequenceUp2Bit CSequenceUp2Bit::GetSubsequence(int startCharInd, int charCount) const
{
	CSequenceUp2Bit newSeq;

	assert(startCharInd >= 0 && charCount >= 0);
	assert(startCharInd + charCount <= GetLength());
	newSeq.up2BitCode = ((up2BitCode >> (startCharInd * 2)) &
					((c_one << (charCount * 2)) - 1)) |
				(c_one << (charCount * 2));
	return newSeq;
}


void CDna2Bits::AssignFromString(const std::string &dnaStr)
{
	m_len = int(dnaStr.length());

	int portionCount = (m_len - 1) / c_dnaCharsInPortion + 1;
	
	if (portionCount > int(m_bits.size()))
		m_bits.resize(portionCount);

	for (int portionInd = 0; portionInd < portionCount; portionInd++)
	{
		TDna2BitsPortion curPortion = 0;
		int curCharCount = c_dnaCharsInPortion;

		if (portionInd == portionCount - 1)
			curCharCount = (m_len - 1) % c_dnaCharsInPortion + 1;
		
		for (int i = 0; i < curCharCount; i++)
		{
			char ch = dnaStr[portionInd * c_dnaCharsInPortion + i];

			switch (ch)
			{
				case 'A':
				case 'C':
				case 'T':
				case 'G':
					break;
				default:
					throw std::exception("Invalid input character");
			}
			curPortion |= TDna2BitsPortion((ch >> 1) & 3) << (i * 2);
		}

		m_bits[portionInd] = curPortion;
	}
}

void CDna2Bits::GetSubsequenceUp2Bit(CSequenceUp2Bit &seqUp2Bit, int startCharInd, int charCount) const
{
	// TODO? to process empty sequences separately and to return assert(charCount > 0
	assert(charCount >= 0 && charCount <= CSequenceUp2Bit::c_maxLen);
	assert(sizeof(CSequenceUp2Bit::TSequenceUp2BitStorage) <= 8);   // This method implemented for up to 8-byte codes
	assert(sizeof(TDna2BitsPortion) == 8);                          // and 4-byte processing portions

	// Subsequence will touch at maximum 2 portions: ... <8-byte portion N> <8-byte portion N + 1>...
	//                                                            ^<subsequence>^
	int startPortionInd = startCharInd / c_dnaCharsInPortion;
	TDna2BitsPortion finalPortion = m_bits[startPortionInd] >> (startCharInd % c_dnaCharsInPortion * 2);
	int lastPortionInd = (startCharInd + charCount - 1) / c_dnaCharsInPortion;

	if (lastPortionInd > startPortionInd)     // TODO? to remove this conditional check by adding one more m_bits cell with value 0
	{
		int startPortionCharCount = c_dnaCharsInPortion - startCharInd % c_dnaCharsInPortion;
	
		finalPortion |= (m_bits[lastPortionInd] & 
						((c_one << ((charCount - startPortionCharCount) * 2)) - 1)) <<
				  (startPortionCharCount * 2);
	}

	// TODO: to test by means of assert(GetSubsequenceUp2Bit_Slow( ), 
  //                                  seqUp2Bit.up2BitCode == finalPortion | (c_one << (charCount * 2)));
	seqUp2Bit.up2BitCode = finalPortion | (c_one << (charCount * 2));
}

int CDna2Bits::GetEqualCharCount(CSequenceUp2Bit &seqUp2Bit, int startCharInd) const
{
	assert(sizeof(CSequenceUp2Bit::TSequenceUp2BitStorage) == sizeof(TDna2BitsPortion)); 
	//assert(startCharInd >= 0 && startCharInd < GetLength());
	assert(startCharInd >= 0 && startCharInd <= GetLength());

	int seqLen = seqUp2Bit.GetLength();
	int restCharCount = GetLength() - startCharInd;
	CSequenceUp2Bit subseqUp2Bit;

	if (restCharCount > CSequenceUp2Bit::c_maxLen)
		restCharCount = CSequenceUp2Bit::c_maxLen;
	//assert(restCharCount * 4 <= sizeof(TDna2BitsPortion));
	GetSubsequenceUp2Bit(subseqUp2Bit, startCharInd, restCharCount);

	// Leading 01s will not disturb
	CSequenceUp2Bit::TSequenceUp2BitStorage code1 = subseqUp2Bit.up2BitCode;
	CSequenceUp2Bit::TSequenceUp2BitStorage code2 = seqUp2Bit.up2BitCode;
	//// Removing leading 01s and getting pure 2-bit arrays of DNA characters
	//CSequenceUp2Bit::TSequenceUp2BitStorage code1 = subseqUp2Bit.up2BitCode ^ 
	//			(c_one << (restCharCount * 2));
	//CSequenceUp2Bit::TSequenceUp2BitStorage code2 = seqUp2Bit.up2BitCode ^ 
	//			(c_one << (seqLen * 2));

	// XOR between 2-bit arrays will reveal differences
	CSequenceUp2Bit::TSequenceUp2BitStorage codesDiff = code1 ^ code2;
	int diffCharInd = GetFirstNonZeroCharacter(codesDiff);

	return std::min(diffCharInd, std::min(restCharCount, seqLen));
}

int CDna2Bits::GetFirstNonZeroCharacter(TDna2BitsPortion codesDiff) 
{
	int zeroBitCount = 0;
	unsigned int curBlock;

	assert(sizeof(TDna2BitsPortion) == 8);             // This method implemented for 8-byte codes
	if ((codesDiff & 0xFFFFFFFF) == 0)
	{
		zeroBitCount = 32;
		curBlock = unsigned int(codesDiff >> 32);
	}
	else 
		curBlock = unsigned int(codesDiff);
	
	if (curBlock == 0)
		return (zeroBitCount + 32) / 2;    // I suppose this is very typical case

	for (int mid = 16; mid != 1; mid /= 2)
	{
		if ((curBlock & ((1 << mid) - 1)) == 0)
		{
			zeroBitCount += mid;
			curBlock >>= mid;
		}
	}

	return zeroBitCount / 2;
}


//bool CDna2Bits::AreSubsequencesEqual(const CSequenceUp2Bit &seqUp2Bit, 
//		                    int thisStartCharInd, int charCount) const
//{
//	assert(sizeof(CSequenceUp2Bit::TSequenceUp2BitStorage) == sizeof(TDna2BitsPortion)); 
//
//
//
//	return false;
//}
//

void CDna2Bits::RunUnitTests()
{
	if (GetFirstNonZeroCharacter(0) != sizeof(TDna2BitsPortion) * 4)
		throw std::exception("CDna2Bits unit test failed");
	assert(sizeof(TDna2BitsPortion) == 8);  
	for (int i = 0; i < 64; i++)
	{
		if (GetFirstNonZeroCharacter(c_one << i) != i / 2)
			throw std::exception("CDna2Bits unit test failed"); 
		if (GetFirstNonZeroCharacter((c_one << i) | (c_one << 62)) != i / 2)
			throw std::exception("CDna2Bits unit test failed"); 
	}	
}


void DnaBase_RunUnitTests()
{
	CDna2Bits::RunUnitTests();
	

}