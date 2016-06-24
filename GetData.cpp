#include "structure.h"

int iChromsomeNum;
map<int64_t, int> ChrLocMap;
int64_t GenomeSize, TwoGenomeSize;
vector<Chromosome_t> ChromosomeVec;

void Convert2Cap(string& seq)
{
	for (string::iterator iter = seq.begin(); iter != seq.end(); iter++)
	{
		if (islower(*iter)) *iter = toupper(*iter);
	}
}

int IdentifyHeaderBoundary(string& str)
{
	int i, len = (int)str.length();
	for (i = 1; i < len; i++) if (str[i] == ' ' || str[i] == '\t') break;

	return i;
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
				getline(file, str); seq += str;
				file.get(ch);
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
	string header, seq, quality;
	int i, p, len, iCount = 0, BaseCount = 0;

	while (!file.eof())
	{
		if (!GetNextEntry(file, header, seq)) break;
		//header
		p = IdentifyHeaderBoundary(header);
		ReadArr[iCount].header = new char[p + 1]; ReadArr[iCount].header[p] = '\0';
		strncpy(ReadArr[iCount].header, header.c_str(), p);

		//sequence
		Convert2Cap(seq);
		ReadArr[iCount].rlen = len = (int)seq.length();
		ReadArr[iCount].seq = new char[len + 1]; ReadArr[iCount].seq[len] = '\0';
		strcpy(ReadArr[iCount].seq, seq.c_str());

		ReadArr[iCount].EncodeSeq = new uint8_t[len];
		for (i = 0; i != len; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];

		iCount++; BaseCount += len;

		if (bSepLibrary) GetNextEntry(file2, header, seq);
		else if (!GetNextEntry(file, header, seq)) break;

		//header
		p = IdentifyHeaderBoundary(header);
		ReadArr[iCount].header = new char[p + 1]; ReadArr[iCount].header[p] = '\0';
		strncpy(ReadArr[iCount].header, header.c_str(), p);
		//sequence
		Convert2Cap(seq);
		ReadArr[iCount].rlen = len = (int)seq.length();
		ReadArr[iCount].seq = new char[len + 1]; ReadArr[iCount].seq[len] = '\0';
		if (bPairEnd) GetComplementarySeq(len, (char*)seq.c_str(), ReadArr[iCount].seq);
		else strcpy(ReadArr[iCount].seq, seq.c_str());

		ReadArr[iCount].EncodeSeq = new uint8_t[len];
		for (i = 0; i != len; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];

		iCount++; BaseCount += len;
		if (iCount == ReadChunkSize || BaseCount >= 1000000) break;
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
