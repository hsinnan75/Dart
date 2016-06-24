#include "structure.h"

int iChromsomeNum;
map<int64_t, int> ChrLocMap;
int64_t GenomeSize, TwoGenomeSize;
vector<Chromosome_t> ChromosomeVec;

static const char CapBaseMap[255] =
{
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /*   0 -   9 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /*  10 -  19 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /*  20 -  29 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /*  30 -  39 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /*  40 -  49 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /*  50 -  59 */
	'N', 'N', 'N', 'N', 'N', 'A', 'N', 'C', 'N', 'N', /*  60 -  69 */
	'N', 'G', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /*  70 -  79 */
	'N', 'N', 'N', 'N', 'T', 'N', 'N', 'N', 'N', 'N', /*  80 -  89 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N', 'C', /*  90 -  99 */
	'N', 'N', 'N', 'G', 'N', 'N', 'N', 'N', 'N', 'N', /* 100 - 109 */
	'N', 'N', 'N', 'N', 'N', 'N', 'T', 'N', 'N', 'N', /* 110 - 119 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /* 120 - 129 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /* 130 - 139 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /* 140 - 149 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /* 150 - 159 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /* 160 - 169 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /* 170 - 179 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /* 180 - 189 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /* 190 - 199 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /* 200 - 209 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /* 210 - 219 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /* 220 - 229 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /* 230 - 239 */
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', /* 240 - 249 */
	'N', 'N', 'N', 'N', 'N'                           /* 250 - 254 */
};

void Convert2Cap(string& seq)
{
	for (string::iterator iter = seq.begin(); iter != seq.end(); iter++)
		*iter = CapBaseMap[(int)(*iter)];
}

int IdentifyHeaderBoundary(string& str)
{
	int i = 0;
	for (string::iterator iter = str.begin(); iter != str.end(); iter++)
	{
		if (*iter == ' ' || *iter == '\t')
		{
			i = iter - str.begin();
			break;
		}
	}
	if (i == 0) i = (int)str.length();

	return i - 1;
}

bool GetNextEntry(fstream& file, string& header, string& seq)
{
	char ch;
	string str;

	getline(file, header);
	if (header == "" || file.eof()) return false;
	else
	{
		if (FastQFormat)
		{
			getline(file, seq); getline(file, str), getline(file, str);
		}
		else
		{
			seq.clear();
			while (!file.eof())
			{
				getline(file, str); seq += str; file.get(ch);
				if (file.eof()) break;
				else file.unget();
				if (ch == '>') break;
			}
		}
		return true;
	}
}

int GetNextChunk(bool bSepLibrary, fstream& file, fstream& file2, ReadItem_t* ReadArr)
{
	int i, p, len, iCount = 0;
	string header, seq, quality;

	while (!file.eof())
	{
		if (!GetNextEntry(file, header, seq)) break;
		//header
		p = IdentifyHeaderBoundary(header);
		ReadArr[iCount].header = new char[p + 1]; ReadArr[iCount].header[p] = '\0';
		strncpy(ReadArr[iCount].header, header.c_str() + 1, p);

		//sequence
		Convert2Cap(seq);
		ReadArr[iCount].rlen = len = (int)seq.length();
		ReadArr[iCount].seq = new char[len + 1];
		strcpy(ReadArr[iCount].seq, seq.c_str());

		ReadArr[iCount].EncodeSeq = new uint8_t[len];
		for (i = 0; i != len; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];

		iCount++;

		if (bSepLibrary) GetNextEntry(file2, header, seq);
		else if (!GetNextEntry(file, header, seq)) break;

		//header
		p = IdentifyHeaderBoundary(header);
		ReadArr[iCount].header = new char[p + 1]; ReadArr[iCount].header[p] = '\0';
		strncpy(ReadArr[iCount].header, header.c_str() + 1, p);
		//sequence
		Convert2Cap(seq);
		ReadArr[iCount].rlen = len = (int)seq.length();
		ReadArr[iCount].seq = new char[len + 1];
		if (bPairEnd) GetComplementarySeq(len, (char*)seq.c_str(), ReadArr[iCount].seq);
		else strcpy(ReadArr[iCount].seq, seq.c_str());

		ReadArr[iCount].EncodeSeq = new uint8_t[len];
		for (i = 0; i != len; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];

		iCount++;
		if (iCount == ReadChunkSize) break;
	}

	return iCount;
}

bool CheckBWAIndexFiles(string IndexPrefix)
{
	fstream file;
	string filename;
	bool bChecked=true;

	filename = IndexPrefix + ".ann"; file.open(filename.c_str(), ios_base::in);
	if(!file.is_open()) return false; else file.close();
	
	filename = IndexPrefix + ".amb"; file.open(filename.c_str(), ios_base::in);
	if(!file.is_open()) return false; else file.close();

	filename = IndexPrefix + ".pac"; file.open(filename.c_str(), ios_base::in);
	if(!file.is_open()) return false; else file.close();

	return bChecked;
}
