#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <string>
#include <map>

using namespace std;

int GetNextStop(string& CIGAR, int pos)
{
	int iRetPos;

	for (iRetPos = pos + 1; iRetPos<(int)CIGAR.length(); iRetPos++)
	{
		if (isalpha(CIGAR[iRetPos])) break;
	}
	return iRetPos;
}

pair<int, int> CalSeqIdy(string& aln1, string& aln2)
{
	int i, idy, len;
	for (len = (int)aln1.length(), idy = 0, i = 0; i < len; i++) if (aln1[i] == aln2[i]) idy++;

	return make_pair(idy, len);
}

pair<int, int> CalSeqIdentity(int rlen, int chrlen, int gPos, string& CIGAR, string& qseq, string& rseq)
{
	string aln1, aln2, tmp;
	int len, iPos1, iPos2, rPos, num, idy;

	len = (int)CIGAR.length();
	for (rPos = iPos1 = iPos2 = 0;;)
	{
		if ((iPos2 = GetNextStop(CIGAR, iPos1)) == len) break;

		tmp = CIGAR.substr(iPos1, iPos2 - iPos1);
		num = atoi(tmp.c_str());

		if ((CIGAR[iPos2] == 'M' || CIGAR[iPos2] == 'I' || CIGAR[iPos2] == 'S') && rPos + num > rlen) break;
		if ((CIGAR[iPos2] == 'M' || CIGAR[iPos2] == 'D') && gPos + num > chrlen) break;

		if (CIGAR[iPos2] == 'I')
		{
			aln1.append(qseq.substr(rPos, num)); rPos += num;
			aln2.insert(aln2.end(), num, '-');
		}
		else if (CIGAR[iPos2] == 'D')
		{
			aln1.insert(aln1.end(), num, '-');
			aln2.append(rseq.substr(gPos, num)); gPos += num;
		}
		else if (CIGAR[iPos2] == 'S')
		{
			rPos += num;
		}
		else if (CIGAR[iPos2] == 'N')
		{
			gPos += num;
		}
		else if (CIGAR[iPos2] != 'H')
		{
			aln1.append(qseq.substr(rPos, num)); rPos += num;
			aln2.append(rseq.substr(gPos, num)); gPos += num;
		}
		fflush(stdout);

		iPos1 = iPos2 + 1;
		if (iPos2 == (int)CIGAR.length() - 1) break;
	}
	return CalSeqIdy(aln1, aln2);
}

int main(int argc, char* argv[])
{
	fstream file;
	stringstream ss;
	pair<int, int> IntPair;
	map<string, string> RefSeqMap;
	int64_t iTotalIdy, iIdenticalBase, iTotalLength;
	string str, tmp, pre_qname, qname, qseq, chrname, chrseq, CIGAR;
	int iFlag, rPos, gPos, MAPQ, num, rlen, chrlen, iTotal, iAln, iHit;

	iTotal = iAln = 0; iTotalIdy = iIdenticalBase = iTotalLength = 0;
	file.open("hg38s.fa", ios_base::in); // read reference sequences
	while (!file.eof())
	{
		getline(file, str);
		if ((int)str.length()>0 && str[0] == '>')
		{
			if (chrseq != "") RefSeqMap.insert(make_pair(chrname, chrseq));
			chrname = str.substr(1); chrseq = "";
		}
		else chrseq.append(str);
	}
	file.close();

	RefSeqMap.insert(make_pair(chrname, chrseq));

	file.clear(); file.open(argv[1], ios_base::in);
	while (!file.eof())
	{
		getline(file, str); if (str == "") break;
		if (str[0] == '@') continue;

		ss.clear(); ss.str(str);
		ss >> qname >> iFlag >> chrname >> gPos >> MAPQ >> CIGAR >> tmp >> tmp >> tmp >> qseq;

		if (pre_qname != qname)
		{
			iHit = 1;
			pre_qname = qname;
		}
		else if (++iHit > 2) continue;

		iTotal++;
		if (CIGAR == "*" || (gPos -= 1) < 0 || RefSeqMap.find(chrname) == RefSeqMap.end()) continue;

		iAln++; rlen = (int)qseq.length(); chrlen = (int)RefSeqMap[chrname].length();
		for (string::iterator iter = qseq.begin(); iter != qseq.end(); iter++) if (islower(*iter)) *iter = toupper(*iter);

		IntPair = CalSeqIdentity(rlen, chrlen, gPos, CIGAR, qseq, RefSeqMap[chrname]);
		if (IntPair.second > 0)
		{
			iIdenticalBase += IntPair.first; iTotalLength += IntPair.second;
			iTotalIdy += (1000 * IntPair.first / IntPair.second);
		}
		//else printf("%s\n", str.c_str());
		if (iAln % 10000 == 0) fprintf(stderr, "\rsensitivity = %d / %d = %.3f, AvgSeqIdy = %.3f, AvgIdenticalBase = %.3f, AvgAlnLen = %.3f", iAln, iTotal, (1.0*iAln / iTotal + 0.0005), (1.0*iTotalIdy / iAln / 1000.0 + 0.0005), 1.*iIdenticalBase / iAln, 1.*iTotalLength / iAln);
	}
	file.close();
	
	if(iAln > 0) fprintf(stderr, "\rsensitivity = %d / %d = %.3f, AvgSeqIdy = %.3f, AvgIdenticalBase = %.3f, AvgAlnLen = %.3f\n\n", iAln, iTotal, (1.0*iAln / iTotal + 0.0005), (1.0*iTotalIdy / iAln / 1000.0 + 0.0005), 1.*iIdenticalBase / iAln, 1.*iTotalLength / iAln);
	else fprintf(stderr, "\rsensitivity = 0, AvgSeqIdy = 0\n");
	//if (iTotalLength > 0) fprintf(stderr, "Identical_BasePairs / TotalAlnLength = %lld / %lld = %.3f\n", iIdenticalBase, iTotalLength, 1.0*iIdenticalBase / iTotalLength);

	return 0;
}
