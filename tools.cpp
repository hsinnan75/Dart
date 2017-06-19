#include "structure.h" 

static const char ReverseMap[255] =
{
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*   0 -   9 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  10 -  19 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  20 -  29 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  30 -  39 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  40 -  49 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  50 -  59 */
		'\0', '\0', '\0', '\0', '\0', 'T',  '\0', 'G',  '\0', '\0', /*  60 -  69 */
		'\0', 'C',  '\0', '\0', '\0', '\0', '\0', '\0', 'N',  '\0', /*  70 -  79 */
		'\0', '\0', '\0', '\0', 'A',  'A',  '\0', '\0', '\0', '\0', /*  80 -  89 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', 'T',  '\0', 'G',  /*  90 -  99 */
		'\0', '\0', '\0', 'C',  '\0', '\0', '\0', '\0', '\0', '\0', /* 100 - 109 */
		'N',  '\0', '\0', '\0', '\0', '\0', 'A',  'A',  '\0', '\0', /* 110 - 119 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 120 - 129 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 130 - 139 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 140 - 149 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 150 - 159 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 160 - 169 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 170 - 179 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 180 - 189 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 190 - 199 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 200 - 209 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 210 - 219 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 220 - 229 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 230 - 239 */
		'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 240 - 249 */
		'\0', '\0', '\0', '\0', '\0'                                /* 250 - 254 */
};

void GetComplementarySeq(int len, char* seq, char* rseq)
{
	int i, j;

	for (j = len - 1, i = 0; i < j; i++, j--)
	{
		rseq[i] = ReverseMap[(int)seq[j]];
		rseq[j] = ReverseMap[(int)seq[i]];
	}
	if (i == j) rseq[i] = ReverseMap[(int)seq[i]];
}

int CalFragPairNonIdenticalBases(int len, char* frag1, char* frag2)
{
	int i, c;

	for (c = 0, i = 0; i < len; i++) if (frag1[i] != frag2[i]) c++;

	return c;
}

int CalFragPairIdenticalBases(int len, char* frag1, char* frag2)
{
	int i, c;

	for (c = 0, i = 0; i < len; i++) if (frag1[i] != frag2[i]) c++;

	return (len - c);
}

int AddNewCigarElements(string& str1, string& str2, vector<pair<int,char> >& cigar_vec)
{
	char state = '*';
	int i, c, len, score = 0;

	for (len = (int)str1.length(), c = 0, i = 0; i < len; i++)
	{
		if (str1[i] == '-')
		{
			if (state == 'D') c++;
			else
			{
				if (c > 0)
				{
					cigar_vec.push_back(make_pair(c, state));
					//if (bDebugMode) printf("%d%c", c, state);
				}
				c = 1; state = 'D';
			}
		}
		else if (str2[i] == '-')
		{
			if (state == 'I') c++;
			else
			{
				if (c > 0)
				{
					cigar_vec.push_back(make_pair(c, state));
					//if (bDebugMode) printf("%d%c", c, state);
				}
				c = 1; state = 'I';
			}
		}
		else
		{
			if (str1[i] == str2[i]) score++;

			if (state == 'M') c++;
			else
			{
				if (c > 0)
				{
					cigar_vec.push_back(make_pair(c, state));
					//if (bDebugMode) printf("%d%c", c, state);
				}
				c = 1; state = 'M';
			}
		}
	}
	if (c > 0)
	{
		cigar_vec.push_back(make_pair(c, state));
		//if (bDebugMode) printf("%d%c", c, state);
	}
	return score;
}

void ShowSeedLocationInfo(int64_t MyPos)
{
	int64_t gPos;

	map<int64_t, int>::const_iterator iter = ChrLocMap.lower_bound(MyPos);
	if (MyPos < GenomeSize) gPos = MyPos - ChromosomeVec[iter->second].FowardLocation;
	else gPos = iter->first - MyPos;
	printf("\t\t\t\t\tChr [%s, %lld]\n", ChromosomeVec[iter->second].name, gPos);
}

void ShowSeedInfo(vector<SeedPair_t>& SeedPairVec)
{
	for (vector<SeedPair_t>::const_iterator iter = SeedPairVec.begin(); iter != SeedPairVec.end(); iter++)
	{
		if (iter->rLen > 0 || iter->gLen > 0)
		{
			printf("\t\tseed#%d: R[%d-%d]=%d G[%lld-%lld]=%d Diff=%lld %s\n", (int)(iter - SeedPairVec.begin() + 1), iter->rPos, iter->rPos + iter->rLen - 1, iter->rLen, iter->gPos, iter->gPos + iter->gLen - 1, iter->gLen, iter->PosDiff, (iter->bSimple ? "Simple" : "Normal"));
			if(iter->gPos < GenomeSize) ShowSeedLocationInfo(iter->gPos);
			else ShowSeedLocationInfo(iter->gPos + iter->gLen - 1);
		}
	}
	printf("\n\n"); fflush(stdout);
}

int ProcessNormalSequencePair(char* seq, SeedPair_t& sp, vector<pair<int, char> >& cigar_vec)
{
	int score = 0;

	if (sp.PosDiff == -1)
	{
		cigar_vec.push_back(make_pair(sp.rLen, 'S'));
	}
	else if (sp.rLen == 0 || sp.gLen == 0)
	{
		if (sp.rLen > 0) cigar_vec.push_back(make_pair(sp.rLen, 'I')); // insertion
		else if (sp.gLen > 0) cigar_vec.push_back(make_pair(sp.gLen, 'D')); // deletion
	}
	else
	{
		string frag1, frag2;
		frag1.resize(sp.rLen); strncpy((char*)frag1.c_str(), seq + sp.rPos, sp.rLen);
		frag2.resize(sp.gLen); strncpy((char*)frag2.c_str(), RefSequence + sp.gPos, sp.gLen);

		if (sp.rLen == sp.gLen && (sp.rLen - (score = CalFragPairIdenticalBases(sp.rLen, (char*)frag1.c_str(), (char*)frag2.c_str())) <= (sp.rLen < 25 ? 5 : (int)(sp.rLen*.2))))
		{
			cigar_vec.push_back(make_pair(sp.rLen, 'M'));
		}
		else
		{
			PairwiseSequenceAlignment(sp.rLen, frag1, sp.gLen, frag2);
			score = AddNewCigarElements(frag1, frag2, cigar_vec);
		}
		//if (bDebugMode) printf("NormalPair:\n%s #read[%d-%d]=%d\n%s #chr[%lld-%lld]=%d\nScore=%d\n\n", frag1.c_str(), sp.rPos, sp.rPos + sp.rLen - 1, sp.rLen, frag2.c_str(), sp.gPos, sp.gPos + sp.gLen - 1, sp.gLen, score);
	}
	return score;
}

int ProcessHeadSequencePair(char* seq, SeedPair_t& sp, vector<pair<int, char> >& cigar_vec)
{
	int score;
	string frag1, frag2;

	frag1.resize(sp.rLen); strncpy((char*)frag1.c_str(), seq + sp.rPos, sp.rLen);
	frag2.resize(sp.gLen); strncpy((char*)frag2.c_str(), RefSequence + sp.gPos, sp.gLen);

	//if (bDebugMode) printf("Head:\n%s #read[%d-%d]=%d\n%s #chr[%lld-%lld]=%d\n\n\n", frag1.c_str(), sp.rPos, sp.rPos + sp.rLen - 1, sp.rLen, frag2.c_str(), sp.gPos, sp.gPos + sp.gLen - 1, sp.gLen), fflush(stdout);

	if (sp.rLen == sp.gLen && (sp.rLen - (score = CalFragPairIdenticalBases(sp.rLen, (char*)frag1.c_str(), (char*)frag2.c_str()))) <= (sp.rLen < 25 ? 5 : (int)(sp.rLen*0.2))) cigar_vec.push_back(make_pair(sp.rLen, 'M'));
	else
	{
		PairwiseSequenceAlignment(sp.rLen, frag1, sp.gLen, frag2);

		//Case1: -X..X vs XX..X (leading gaps in the read block)
		int p = 0; while (frag1[p] == '-') p++;
		if (p > 0) // shrink the genome block
		{
			frag1.erase(0, p); frag2.erase(0, p);
			sp.gPos += p; sp.gLen -= p;
		}
		//Case2: XX..X vs -X..X (leading gaps in the genome block)
		p = 0; while (frag2[p] == '-') p++;
		if (p > 0) // shrink the read block
		{
			frag1.erase(0, p); frag2.erase(0, p);
			sp.rPos += p; sp.rLen -= p; cigar_vec.push_back(make_pair(p, 'S'));
		}
		score = AddNewCigarElements(frag1, frag2, cigar_vec);
		//if (bDebugMode) printf("H2:\n%s #read[%d-%d]=%d\n%s #chr[%lld-%lld]=%d\nScore=%d\n\n", frag1.c_str(), sp.rPos, sp.rPos + sp.rLen - 1, sp.rLen, frag2.c_str(), sp.gPos, sp.gPos + sp.gLen - 1, sp.gLen, score);
	}
	return score;
}

int ProcessTailSequencePair(char* seq, SeedPair_t& sp, vector<pair<int, char> >& cigar_vec)
{
	int score;
	string frag1, frag2;

	frag1.resize(sp.rLen); strncpy((char*)frag1.c_str(), seq + sp.rPos, sp.rLen);
	frag2.resize(sp.gLen); strncpy((char*)frag2.c_str(), RefSequence + sp.gPos, sp.gLen);

	//if (bDebugMode) printf("Tail:\n%s #read[%d-%d]=%d\n%s #chr[%lld-%lld]=%d\n\n\n", frag1.c_str(), sp.rPos, sp.rPos + sp.rLen - 1, sp.rLen, frag2.c_str(), sp.gPos, sp.gPos + sp.gLen - 1, sp.gLen), fflush(stdout);

	if (sp.rLen == sp.gLen && (sp.rLen - (score = CalFragPairIdenticalBases(sp.rLen, (char*)frag1.c_str(), (char*)frag2.c_str()))) <= (sp.rLen < 25 ? 5 : (int)(sp.rLen*0.2))) cigar_vec.push_back(make_pair(sp.rLen, 'M'));
	else
	{
		PairwiseSequenceAlignment(sp.rLen, frag1, sp.gLen, frag2);

		//Case1: X..X- vs X..XX (tailing gaps in the read block)
		int p = (int)frag1.length() - 1, c = 0; while (frag1[p--] == '-') c++;
		if (c > 0) // shrink the genome block
		{
			frag1.resize((int)frag1.length() - c); frag2.resize((int)frag2.length() - c);
			sp.gLen -= c;
			//if (bDebugMode) printf("find %d tailing gaps in the read block\n", c);
		}
		//Case2: X..XX vs X..X- (tailing gaps in the genome block)
		p = (int)frag2.length() - 1, c = 0; while (frag2[p--] == '-') c++;
		if (c > 0) // shrink the read block
		{
			frag1.resize((int)frag1.length() - c); frag2.resize((int)frag2.length() - c);
			sp.rLen -= c;
			if (bDebugMode) printf("find %d tailing gaps in the genome block\n", c);
		}
		score = AddNewCigarElements(frag1, frag2, cigar_vec); if (c > 0) cigar_vec.push_back(make_pair(c, 'S'));
		//if (bDebugMode) printf("T2:\n%s #read[%d-%d]=%d\n%s #chr[%lld-%lld]=%d\nScore=%d\n\n", frag1.c_str(), sp.rPos, sp.rPos + sp.rLen - 1, sp.rLen, frag2.c_str(), sp.gPos, sp.gPos + sp.gLen - 1, sp.gLen, score);
	}
	return score;
}
