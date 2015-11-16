#include "stdafx.h"
#include "DnaTrieBuilder.h"

void TrimRight(std::string &st)
{
	size_t i;

	i = st.size () - 1;
	while (i>0 && (st[i-1]<=' '))
		i--;
	st.resize (i);
}

void LoadDataToTrie(CDnaTrieBuilder &trieBuilder, const std::string &inputFileName)
{
	static const int c_maxStrLen = 16384;
	FILE *f = stdin;

	if (!inputFileName.empty())
		f = fopen(inputFileName.c_str(), "r");

	char buf[c_maxStrLen];
	std::string str;
	int readStrCount = 0;

	if (!f)
		throw std::exception("Couldn't open input file");

	buf[0] = 0;
	while (true)
	{
		FUNC_GUARD

		char *res = fgets(buf, c_maxStrLen, f);

		str = buf;
		TrimRight(str);
		if (!res || str.empty())
		{
			printf("%d lines read\n", readStrCount);
			trieBuilder.PrintCheckSum();
			if (!inputFileName.empty())
				fclose(f);
			return;
		}

		std::string::size_type commaPos = str.rfind(',');
		int count = -1;

		if (commaPos == std::string::npos)
			throw std::exception("No comma in input string");
		if (sscanf(str.c_str() + commaPos + 1, "%d", &count) != 1 || count <= 0)
			throw std::exception("Invalid DNA count");
	
		str.resize(commaPos);

		trieBuilder.AddDna(str, count);

		readStrCount++;
		if (readStrCount % (1 << 17) == 0)
			printf("%d lines read\n", readStrCount);
	}
}

int main(int argc, char* argv[])
{
	try
	{
		printf("DNA Trie Builder v. 2.0\n\n");
		if (argc <= 2)
		{
			printf("Usage: DnaTrieBuilder.exe [-i <input text file>] -o <output file> <options>\n\n"
				     "Parameters and options can be in any order.\n"
				     "Output file sets beginning of produced output files, "
						 "actual files will be named *.A, *.C, *.G, *.T, *.COUNT and *.SEQ.\n"
						 "Options: --rows - number of table rows which are allocated at start");
			return -1;
		}

		DnaBase_RunUnitTests();

		std::string inputFileName, outputFileNameBegin;
		int allocRowCount = 10000000;
		CDnaTrieBuilder trieBuilder;

		for (int i = 1; i < argc; i++)
		{
			char *str = argv[i];
			bool nextArgExists = (i + 1 < argc);
			char *nextStr = argv[i + 1];

			if (stricmp(str, "-i") == 0)
			{
				if (!nextArgExists)
					throw std::exception("Input file is not specified");
				inputFileName = nextStr;
				i++;
			}
			else  if (stricmp(str, "-o") == 0)
			{
				if (!nextArgExists)
					throw std::exception("Output file is not specified");
				outputFileNameBegin = nextStr;
				i++;
			}
			else  if (stricmp(str, "--rows") == 0)
			{
				if (!nextArgExists)
					throw std::exception("Rows value is not specified");
				if (sscanf(nextStr, "%d", &allocRowCount) != 1 || allocRowCount < 1)
					throw std::exception("Invalid rows value");
				i++;
			}
		}

		trieBuilder.ResizeTable(allocRowCount);

		if (outputFileNameBegin.empty())
			throw std::exception("No output file specified");
		
		LoadDataToTrie(trieBuilder, inputFileName);
		trieBuilder.WriteToFiles(outputFileNameBegin);
		printf("\nDone          \n");

#ifdef FUNC_TIME_MEASUREMENT
		time_func_guard::outputFunctionsStatistics();
#endif
	}
  catch (const std::exception &e)
	{
    printf("\nError: %s", e.what());
	//	_getch();
		return -50;
	}

//	_getch();
	return 0;
}
