#include "StdAfx.h"
#include "DnaTrieBuilder.h"

CDnaTrieBuilder::CDnaTrieBuilder()
	: m_rowCount(1)
{
	m_table.resize(4);

	//TTrieRow &row = m_table[0];
	//
	//row.
}

CDnaTrieBuilder::~CDnaTrieBuilder()
{
}

void CDnaTrieBuilder::ResizeTable(int rowCount)
{
	assert(rowCount > m_rowCount);
	m_table.resize(rowCount);
	printf("%d MB of RAM allocated\n", int((rowCount * sizeof(TTrieRow)) >> 20));
}

std::string g_curDna;  //d_

void CDnaTrieBuilder::AddDna(const std::string &dnaStr, int count)
{
	g_curDna = dnaStr;

	FUNC_GUARD

	if (m_rowCount + 10 + dnaStr.length() / 2 >= m_table.size())
		ResizeTable(int(m_table.size() * 2));

	int rowInd = 0;
	int o = 0;

	m_curDna.AssignFromString(dnaStr);

	int dnaLen = m_curDna.GetLength();
	//const TDna2BitsPortion *pDnaBits = m_curDna.GetBits();

	assert(dnaLen == dnaStr.length());

	while (true)
	{
		TTrieRow &row = m_table[rowInd];
		int rowSeqLen = row.seq.GetLength();

		//if (rowSeqLen == 0)                                    // if SEQ[row] == 1:
		//{
		//	row.count += count;
		//	if (o == dnaLen)
		//		break;

		//	int curChar2Bits = m_curDna.GetChar2Bits(o);         // curChar2Bits = oseq[o]
		//	int nextRowInd = row.GetPipeRowIndex(curChar2Bits);  // temp = int((A,C,T,G)[oseq[o]][row])

		//	if (nextRowInd != TTrieRow::c_emptyRowInd)           // if temp != 0:
		//	{
		//		assert(nextRowInd > rowInd && nextRowInd < m_rowCount);
		//		rowInd = nextRowInd;                               // row = temp
		//		o++;
		//		//printf("%s: row %d, o %d\n", dnaStr.c_str(), rowInd, o);  //d_
		//		continue;
		//	}
		//	else
		//		AttachNewSequence(rowInd, o, count);
		//}
		//else
		{
			int equalCharCount = m_curDna.GetEqualCharCount(row.seq, o);
			int restCharCount = dnaLen - o;

			assert(equalCharCount >= 0 && equalCharCount <= dnaLen - o && equalCharCount <= rowSeqLen);
			if (equalCharCount == rowSeqLen &&
				  restCharCount == rowSeqLen)                        // if up2bit == oseq[o:]: 
			{
				row.count += count;                                  // COUNT[row] += count
        break;					                           
			}

			if (equalCharCount < rowSeqLen)
			{		
				// Splitting existing sequence

				int nextRowToAdd = m_rowCount;
				TTrieRow &nextRow = m_table[nextRowToAdd];

				nextRow = row;                                     // A[nextRowToAdd]     = A[row]   C[...
			//	m_curDna.GetSubsequenceUp2Bit(nextRow.seq, o, rowSeqLen);
				nextRow.seq = row.seq.GetSubsequence(equalCharCount + 1, rowSeqLen - equalCharCount - 1);
				row.SetEmptyRowInds();                             // A[row]              = 0   C[row] =...
				//row.count += count;
				row.SetPipeRowIndex(row.seq.GetChar2Bits(equalCharCount), nextRowToAdd);     // (A,C,T,G)[up2bit[len(oseq)-o]][row] = nextRowToAdd     
				row.seq = row.seq.GetSubsequence(0, equalCharCount);    // SEQ[row] = sum([up2bit[z]<<(2*t) 
				m_rowCount++;
			}

			row.count += count;
			if (equalCharCount < restCharCount)
			{
				//assert(restCharCount > rowSeqLen);
				o += equalCharCount;

				// Very similar to the first block (if (rowSeqLen == 0))
				//if (o == dnaLen)
				//	break;

				int curChar2Bits = m_curDna.GetChar2Bits(o);         // curChar2Bits = oseq[o]
				int nextRowInd = row.GetPipeRowIndex(curChar2Bits);  // temp = int((A,C,T,G)[oseq[o]][row])

				if (nextRowInd != TTrieRow::c_emptyRowInd)           // if temp != 0:
				{
					assert(nextRowInd < m_rowCount);
					//assert(nextRowInd > rowInd && nextRowInd < m_rowCount);
					rowInd = nextRowInd;                               // row = temp
					o++;
					//printf("%s: row %d, o %d\n", dnaStr.c_str(), rowInd, o);  //d_
					continue;
				}
				else
					AttachNewSequence(rowInd, o, count);
			}  
		}
		break;
	}

	//if (m_rowCount % (1 << 18) == 0)
	//{
	//	printf("DNA %s table:\n", dnaStr.c_str());
	//	PrintTable();
	//}
}

void CDnaTrieBuilder::AttachNewSequence(int rowInd, int startCharInd, int count)
{
	int curRowInd = rowInd;
	int curCharInd = startCharInd;
	int dnaLen = m_curDna.GetLength();

	do                                                       // for y in xrange(0, len(oseq)-o, 32):
	{
		assert(curCharInd < dnaLen);
		assert(m_rowCount < int(m_table.size()));
		
		int newSeqLen = dnaLen - curCharInd - 1;
		
		if (newSeqLen > CSequenceUp2Bit::c_maxLen)
			newSeqLen = CSequenceUp2Bit::c_maxLen;

		TTrieRow &row = m_table[rowInd];
		int curChar2Bits = m_curDna.GetChar2Bits(curCharInd);  // oseq[o+y]
		int nextRowToAdd = m_rowCount;
		TTrieRow &nextRow = m_table[nextRowToAdd];

		row.SetPipeRowIndex(curChar2Bits, nextRowToAdd);       // (A,C,T,G)[oseq[o+y]][row] = nextRowToAdd
		nextRow.count = count;                                 // COUNT[nextRowToAdd]       = count
		m_curDna.GetSubsequenceUp2Bit(nextRow.seq, curCharInd + 1, newSeqLen);      // SEQ[nextRowToAdd]         = ...
		m_rowCount++;

		curCharInd += 1 + newSeqLen;
	}
	while (curCharInd != dnaLen);
}                        

#define WRITE_ARRAY_TO_FILE(arr, arrFileName) \
	f = fopen((fileNameBegin + "." #arrFileName).c_str(), "wb"); \
	for (i = 0; i < m_rowCount; i++) \
		fwrite(&(m_table[i].arr), sizeof(m_table[0].arr), 1, f); \
	fclose(f);

void CDnaTrieBuilder::WriteToFiles(const std::string &fileNameBegin) const
{
	FUNC_GUARD

	FILE *f;
	int i;

	// TODO: conversion from Big Endian
	WRITE_ARRAY_TO_FILE(a, A);
	WRITE_ARRAY_TO_FILE(c, C);
	WRITE_ARRAY_TO_FILE(t, T);
	WRITE_ARRAY_TO_FILE(g, G);
	WRITE_ARRAY_TO_FILE(count, COUNT);
	WRITE_ARRAY_TO_FILE(seq.up2BitCode, SEQ);
}

#define PRINT_ARRAY_SUM(arr) sum = 0; \
	for (i = 0; i < m_rowCount; i++) \
		sum += m_table[i].arr; \
	printf("%lld", sum);

void CDnaTrieBuilder::PrintCheckSum()
{
	CSequenceUp2Bit::TSequenceUp2BitStorage sum;
	int i;

	printf("Sums: ");
	PRINT_ARRAY_SUM(a);
	printf(" "); 
	PRINT_ARRAY_SUM(c);
	printf(" "); 
	PRINT_ARRAY_SUM(t);
	printf(" "); 
	PRINT_ARRAY_SUM(g);
	printf(" "); 
	PRINT_ARRAY_SUM(count)
	printf(" "); 
	PRINT_ARRAY_SUM(seq.up2BitCode);
	printf("\n"); 
}

#define PRINT_ARRAY(arr, arrOutName) printf("%s: ", #arrOutName); \
	for (i = 0; i < m_rowCount - 1; i++) \
	  printf("%lld, ", (unsigned long long)(m_table[i].arr)); \
	printf("%lld\n", (unsigned long long)(m_table[i].arr));

void CDnaTrieBuilder::PrintTable()
{
	int i;
	
	PRINT_ARRAY(a, A);
	PRINT_ARRAY(c, C);
	PRINT_ARRAY(t, T);
	PRINT_ARRAY(g, G);
	PRINT_ARRAY(count, Cnt);
	PRINT_ARRAY(seq.up2BitCode, Seq);
}
