#pragma once

#include "DnaBase.h"

struct TTrieRow
{
	static const int c_emptyRowInd = 0;

	int a, c, t, g;       // Child nodes' row indices
	int count;            // TODO? to store seq.GetLength()
	CSequenceUp2Bit seq;

	TTrieRow()
		: a(c_emptyRowInd), c(c_emptyRowInd), t(c_emptyRowInd), g(c_emptyRowInd),
	    count(0)
	{	}

	void SetEmptyRowInds()
	{
		a = c_emptyRowInd;
		c = c_emptyRowInd;
		t = c_emptyRowInd;
		g = c_emptyRowInd;
	}

	void CheckAllRowIndsEmpty()
	{
		assert(a == c_emptyRowInd);
		assert(c == c_emptyRowInd);
		assert(t == c_emptyRowInd);
		assert(g == c_emptyRowInd);
	}

	void SetPipeRowIndex(int pipeInd, int rowInd)
	{
		assert(pipeInd >= 0 && pipeInd < 4);
		assert(rowInd >= 0);
		(&a)[pipeInd] = rowInd;
	}

	int GetPipeRowIndex(int pipeInd) const
	{
		assert(pipeInd >= 0 && pipeInd < 4);
		return (&a)[pipeInd];
	}

};


class CDnaTrieBuilder
{
public:
	CDnaTrieBuilder();
	~CDnaTrieBuilder();
	void ResizeTable(int rowCount);
	void AddDna(const std::string &dnaStr, int count);
	void WriteToFiles(const std::string &fileNameBegin) const;

	void PrintCheckSum();
	void PrintTable();

private:
	typedef std::vector<TTrieRow> TTrieTable;

	int m_rowCount;
	TTrieTable m_table;

	CDna2Bits m_curDna;

	void AttachNewSequence(int rowInd, int startCharInd, int count);
};

