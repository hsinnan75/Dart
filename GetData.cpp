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

int IdentifyHeaderBoundary(char* str, int len)
{
	int i;

	for (i = 1; i < len; i++)
	{
		if (str[i] == ' ' || str[i] == '\t') return i;
	}
	return len - 1;
}

ReadItem_t GetNextEntry(FILE *file)
{
	ssize_t len;
	size_t size = 0;
	ReadItem_t read;
	char *buffer = NULL;

	read.header = read.seq = NULL; read.rlen = 0;

	if ((len = getline(&buffer, &size, file)) != -1)
	{
		len = IdentifyHeaderBoundary(buffer, len) - 1; read.header = new char[len + 1];
		strncpy(read.header, (buffer + 1), len); read.header[len] = '\0';
		if (FastQFormat)
		{
			if ((read.rlen = getline(&buffer, &size, file)) != -1)
			{
				read.seq = new char[read.rlen];
				strncpy(read.seq, buffer, read.rlen);
				read.rlen -= 1; read.seq[read.rlen] = '\0';
				getline(&buffer, &size, file); getline(&buffer, &size, file);
			}
			else read.rlen = 0;
		}
		else
		{
			string seq;
			while (true)
			{
				if ((len = getline(&buffer, &size, file)) == -1) break;
				if (buffer[0] == '>')
				{
					fseek(file, 0 - len, SEEK_CUR);
					break;
				}
				else
				{
					buffer[len - 1] = '\0'; seq += buffer;
				}
			}
			if ((read.rlen = (int)seq.length()) > 0)
			{
				read.seq = new char[read.rlen + 1];
				strcpy(read.seq, (char*)seq.c_str());
				read.seq[read.rlen] = '\0';
			}
		}
	}
	free(buffer);

	return read;
}

int GetNextChunk(bool bSepLibrary, FILE *file, FILE *file2, ReadItem_t* ReadArr)
{
	char* rseq;
	int i, iBase = 0, iCount = 0;

	while (true)
	{
		if ((ReadArr[iCount] = GetNextEntry(file)).rlen == 0) break;
		ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		//printf("%s\n%s\n", ReadArr[iCount].header, ReadArr[iCount].seq);

		for (i = 0; i != ReadArr[iCount].rlen; i++)
		{
			//if (islower(ReadArr[iCount].seq[i])) ReadArr[iCount].seq[i] = toupper(ReadArr[iCount].seq[i]);
			ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		}
		iBase += ReadArr[iCount].rlen; iCount++;

		if (bSepLibrary) ReadArr[iCount] = GetNextEntry(file2);
		else if ((ReadArr[iCount] = GetNextEntry(file)).rlen == 0) break;
		//printf("%s\n%s\n", ReadArr[iCount].header, ReadArr[iCount].seq);

		if (bPairEnd)
		{
			rseq = new char[ReadArr[iCount].rlen];
			GetComplementarySeq(ReadArr[iCount].rlen, ReadArr[iCount].seq, rseq);
			copy(rseq, rseq + ReadArr[iCount].rlen, ReadArr[iCount].seq);
			delete[] rseq;
		}
		ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		for (i = 0; i != ReadArr[iCount].rlen; i++)
		{
			//if (islower(ReadArr[iCount].seq[i])) ReadArr[iCount].seq[i] = toupper(ReadArr[iCount].seq[i]);
			ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		}
		iBase += ReadArr[iCount].rlen; iCount++;
		if (iCount == ReadChunkSize || iBase > 1000000) break;
	}
	return iCount;
}

ReadItem_t gzGetNextEntry(gzFile file)
{
	int len;
	ReadItem_t read;
	char buffer[1024];

	read.header = read.seq = NULL; read.rlen = 0;

	if (gzgets(file, buffer, 1024) != NULL)
	{
		len = IdentifyHeaderBoundary(buffer, strlen(buffer)) - 1;
		if (len > 0 && (buffer[0] == '@' || buffer[0] == '>'))
		{
			read.header = new char[len + 1];
			strncpy(read.header, (buffer + 1), len); read.header[len] = '\0';
			gzgets(file, buffer, 1024); read.rlen = strlen(buffer) - 1; read.seq = new char[read.rlen + 1]; read.seq[read.rlen] = '\0';
			strncpy(read.seq, buffer, read.rlen);

			if (FastQFormat)
			{
				gzgets(file, buffer, 1024); gzgets(file, buffer, 1024);
			}
		}
	}
	return read;
}

int gzGetNextChunk(bool bSepLibrary, gzFile file, gzFile file2, ReadItem_t* ReadArr)
{
	char* rseq;
	int i, iBase = 0, iCount = 0;

	while (true)
	{
		if ((ReadArr[iCount] = gzGetNextEntry(file)).rlen == 0) break;
		ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		for (i = 0; i != ReadArr[iCount].rlen; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		iBase += ReadArr[iCount].rlen; iCount++;

		if (bSepLibrary) ReadArr[iCount] = gzGetNextEntry(file2);
		else if ((ReadArr[iCount] = gzGetNextEntry(file)).rlen == 0) break;

		if (bPairEnd)
		{
			rseq = new char[ReadArr[iCount].rlen];
			GetComplementarySeq(ReadArr[iCount].rlen, ReadArr[iCount].seq, rseq);
			copy(rseq, rseq + ReadArr[iCount].rlen, ReadArr[iCount].seq);
			delete[] rseq;
		}
		ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		for (i = 0; i != ReadArr[iCount].rlen; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		iBase += ReadArr[iCount].rlen; iCount++;

		if (iCount == ReadChunkSize || iBase > 1000000) break;
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
