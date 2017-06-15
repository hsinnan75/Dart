#include "structure.h"

#define MinIntronSize 5

map<int64_t, int> ExonMap;
static pthread_mutex_t ExonLoc;
int ShiftArr[19] = { 0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8, 9, -9 };

typedef struct
{
	int p;
	int left_ext;
	int right_ext;
} GappedExtension_t;

bool CompByPosDiff(const SeedPair_t& p1, const SeedPair_t& p2)
{
	if (p1.PosDiff == p2.PosDiff) return p1.rPos < p2.rPos;
	else return p1.PosDiff < p2.PosDiff;
}

bool CompByGenomePos(const SeedPair_t& p1, const SeedPair_t& p2)
{
	if (p1.gPos == p2.gPos) return p1.rPos < p2.rPos;
	else return p1.gPos < p2.gPos;
}

bool CompByReadPos(const SeedPair_t& p1, const SeedPair_t& p2)
{
	return p1.rPos < p2.rPos;
}

bool CompByFirstInt(const pair<int, int>& p1, const pair<int, int>& p2)
{
	return p1.first < p2.first;
}

string GenerateCIGAR(vector<pair<int, char> >& cigar_vec)
{
	int i, num, c;
	char state, buf[8];
	string cigar_str;

	for (state = '\0', num = (int)cigar_vec.size(), c = 0, i = 0; i != num; i++)
	{
		if (cigar_vec[i].second != state)
		{
			if (c > 0) sprintf(buf, "%d%c", c, state), cigar_str += buf;

			c = cigar_vec[i].first; state = cigar_vec[i].second;
		}
		else c += cigar_vec[i].first;
	}
	if (c > 0) sprintf(buf, "%d%c", c, state), cigar_str += buf;

	//if (bDebugMode) printf("CIGAR=%s\n\n\n", cigar_str.c_str());

	return cigar_str;
}

void ShowSpliceJunctions(char* header, Coordinate_t& coor)
{
	int p, len;
	int64_t gPos = coor.gPos;
	string tmp, CIGAR = coor.CIGAR;

	printf("\nSplice junctions [%s]: %ld\n", header, coor.gPos);
	while (CIGAR != "")
	{
		for (p = 0; p < (int)CIGAR.length(); p++) if (isalpha(CIGAR[p])) break;

		tmp = CIGAR.substr(0, p); len = atoi(tmp.c_str());
		if (CIGAR[p] == 'M' || CIGAR[p] == 'N' || CIGAR[p] == 'D') gPos += len;
		//printf("%d%c --> %ld\n", len, CIGAR[p], gPos);
		if (CIGAR[p] == 'N') ShowSeedLocationInfo(gPos);
		CIGAR = CIGAR.substr(p + 1);
	}
	printf("\n\n");
}

Coordinate_t GenCoordinateInfo(bool bFirstRead, int64_t gPos, int64_t end_gPos, vector<pair<int, char> >& cigar_vec)
{
	Coordinate_t coor;
	map<int64_t, int>::iterator iter;

	if (gPos < GenomeSize) // forward strand
	{
		if (bFirstRead) coor.bDir = true;
		else coor.bDir = false;

		iter = ChrLocMap.lower_bound(gPos);
		coor.ChromosomeIdx = iter->second;
		coor.gPos = gPos + 1 - ChromosomeVec[coor.ChromosomeIdx].FowardLocation;
	}
	else
	{
		if (bFirstRead) coor.bDir = false;
		else coor.bDir = true;

		reverse(cigar_vec.begin(), cigar_vec.end());

		iter = ChrLocMap.lower_bound(gPos);
		coor.gPos = iter->first - end_gPos + 1; coor.ChromosomeIdx = iter->second;
		//if(bDebugMode) printf("matched chr=%s, loc=%ld, gPos: %ld -> %ld\n", ChromosomeVec[coor.ChromosomeIdx].name, ChromosomeVec[coor.ChromosomeIdx].ReverseLocation, gPos, coor.gPos);
	}
	//if (bDebugMode) printf("gPos: %ld --> %ld %s\n", gPos, coor.gPos, (coor.bDir ? "Forward" : "Reverse"));

	coor.CIGAR = GenerateCIGAR(cigar_vec);

	return coor;
}

string ReverseCIGAR(string& CIGAR)
{
	string RevCIGAR;
	int len, pos1, pos2;

	len = (int)CIGAR.length();

	for (pos1 = 0, pos2 = 1; pos1<len; pos2++)
	{
		if (isalpha(CIGAR[pos2]))
		{
			RevCIGAR.insert(0, CIGAR.substr(pos1, pos2 - pos1 + 1));
			pos1 = pos2 + 1;
		}
	}
	return RevCIGAR;
}

bool CheckCoordinateValidity(vector<SeedPair_t>& SeedVec)
{
	bool bValid = true;
	int64_t gPos1 = 0, gPos2 = TwoGenomeSize;

	for (vector<SeedPair_t>::iterator iter = SeedVec.begin(); iter != SeedVec.end(); iter++)
	{
		if (iter->gLen > 0)
		{
			gPos1 = iter->gPos;
			break;
		}
	}
	for (vector<SeedPair_t>::reverse_iterator iter = SeedVec.rbegin(); iter != SeedVec.rend(); iter++)
	{
		if (iter->gLen > 0)
		{
			gPos2 = iter->gPos + iter->gLen - 1;
			break;
		}
	}
	if ((gPos1 < GenomeSize && gPos2 >= GenomeSize) || (gPos1 >= GenomeSize && gPos2 < GenomeSize) || ChrLocMap.lower_bound(gPos1)->second != ChrLocMap.lower_bound(gPos2)->second)
	{
		bValid = false;
		//if (bDebugMode) printf("%ld and %ld are not in the same chromosome!\n", gPos1, gPos2);
	}
	return bValid;
}

bool CheckCandidateValidity(vector<SeedPair_t>& SeedPairVec)
{
	int i, j, num;
	bool bValidity = true;

	for (num = (int)SeedPairVec.size(), i = 0, j = 1; j < num; i++, j++)
	{
		if ((SeedPairVec[i].rLen > 0 && SeedPairVec[j].rLen > 0 && SeedPairVec[i].rPos + SeedPairVec[i].rLen - SeedPairVec[j].rPos > 0) || (SeedPairVec[i].gLen > 0 && SeedPairVec[j].gLen > 0 && SeedPairVec[i].gPos + SeedPairVec[i].gLen - SeedPairVec[j].gPos > 0))
		{
			bValidity = false;
			break;
		}
	}
	return bValidity;
}

vector<SeedPair_t> IdentifySeedPairs(int rlen, uint8_t* EncodeSeq)
{
	SeedPair_t SeedPair;
	int i, pos, end_pos;
	vector<SeedPair_t> SeedPairVec;
	bwtSearchResult_t bwtSearchResult;

	SeedPairVec.clear(); pos = 0, end_pos = rlen - 13;
	SeedPair.bAcceptorSite = false; SeedPair.bSimple = true;

	while (pos < end_pos)
	{
		if (EncodeSeq[pos] > 3) pos++;
		else
		{
			bwtSearchResult = BWT_Search(EncodeSeq, pos, rlen);
			//if (bDebugMode) printf("Pos=%d, Freq=%d, Len=%d\n", pos, bwtSearchResult.freq, bwtSearchResult.len);
			if (bwtSearchResult.freq > 0)
			{
				SeedPair.rPos = pos; SeedPair.rLen = SeedPair.gLen = bwtSearchResult.len;
				for (i = 0; i != bwtSearchResult.freq; i++)
				{
					SeedPair.PosDiff = (SeedPair.gPos = bwtSearchResult.LocArr[i]) - SeedPair.rPos;
					SeedPairVec.push_back(SeedPair);
				}
				delete[] bwtSearchResult.LocArr;
				pos += (bwtSearchResult.len);
			}
			else pos++;
		}
	}
	sort(SeedPairVec.begin(), SeedPairVec.end(), CompByGenomePos);

	return SeedPairVec;
}

void MergeAdjacentSimplePairs(vector<SeedPair_t>& SeedVec)
{
	int i, j, k, shift, num;

	for (num = (int)SeedVec.size(), i = 0; i < num;)
	{
		if (SeedVec[i].rLen == 0) continue;
		//if (bDebugMode) printf("\tr[%d-%d] g[%ld-%ld], len=%d, PosDiff=%ld\n", SeedVec[i].rPos, SeedVec[i].rPos + SeedVec[i].rLen - 1, SeedVec[i].gPos, SeedVec[i].gPos + SeedVec[i].gLen - 1, SeedVec[i].rLen, SeedVec[i].PosDiff);
		for (j = i + 1; j < num; j++)
		{
			if (SeedVec[j].PosDiff != SeedVec[i].PosDiff || SeedVec[j].rPos - SeedVec[j - 1].rPos > 1) break;
			//if (bDebugMode) printf("\tr[%d-%d] g[%ld-%ld], len=%d, PosDiff=%ld\n", SeedVec[j].rPos, SeedVec[j].rPos + SeedVec[j].rLen - 1, SeedVec[j].gPos, SeedVec[j].gPos + SeedVec[j].gLen - 1, SeedVec[j].rLen, SeedVec[j].PosDiff);
		}
		if (j - i > 1)
		{
			shift = SeedVec[j - 1].rPos - SeedVec[i].rPos;
			SeedVec[i].rLen += shift; SeedVec[i].gLen += shift;
			//if (bDebugMode) printf("Merged!!r[%d-%d] g[%ld-%ld], len=%d, PosDiff=%ld\n\n", SeedVec[i].rPos, SeedVec[i].rPos + SeedVec[i].rLen - 1, SeedVec[i].gPos, SeedVec[i].gPos + SeedVec[i].gLen - 1, SeedVec[i].rLen, SeedVec[i].PosDiff);
			for (k = i + 1; k < j; k++) SeedVec[k].rLen = SeedVec[k].gLen = 0;
		}
		i = j;
	}
}

vector<AlignmentCandidate_t> GenerateAlignmentCandidate(int rlen, vector<SeedPair_t>& SeedPairVec)
{
	int64_t PosDiff;
	int i, j, k, thr, num;
	AlignmentCandidate_t AlignmentCandidate;
	vector<AlignmentCandidate_t> AlignmentVec;

	if ((num = SeedPairVec.size() == 0)) return AlignmentVec;
	else AlignmentCandidate.PairedAlnCanIdx = -1;

	thr = (int)(rlen*0.3);

	//if (bDebugMode) printf("\n\nRaw seeds:\n"), ShowSeedInfo(SeedPairVec);
	 num = (int)SeedPairVec.size();

	i = 0; while (i < num && SeedPairVec[i].PosDiff < 0) i++;
	for (; i < num;)
	{
		AlignmentCandidate.Score = SeedPairVec[i].rLen;
		AlignmentCandidate.SeedVec.resize(1); AlignmentCandidate.SeedVec[0] = SeedPairVec[i];
		//if (bDebugMode) printf("Leading seed:  r[%d-%d] g[%ld-%ld], len=%d, PosDiff=%ld\n", SeedPairVec[i].rPos, SeedPairVec[i].rPos + SeedPairVec[i].rLen - 1, SeedPairVec[i].gPos, SeedPairVec[i].gPos + SeedPairVec[i].gLen - 1, SeedPairVec[i].rLen, SeedPairVec[i].PosDiff);
		for (j = i, k = i + 1; k < num; k++)
		{
			PosDiff = abs(SeedPairVec[k].PosDiff - SeedPairVec[j].PosDiff);
			if (PosDiff < MaxGaps || (PosDiff < MaxIntronSize && SeedPairVec[k].rPos > SeedPairVec[j].rPos))
			{
				//if (bDebugMode) printf("add seed (PosDiff=%ld): r[%d-%d] g[%ld-%ld], len=%d, PosDiff=%ld\n", PosDiff, SeedPairVec[k].rPos, SeedPairVec[k].rPos + SeedPairVec[k].rLen - 1, SeedPairVec[k].gPos, SeedPairVec[k].gPos + SeedPairVec[k].gLen - 1, SeedPairVec[k].rLen, SeedPairVec[k].PosDiff);
				AlignmentCandidate.Score += SeedPairVec[k].rLen;
				AlignmentCandidate.SeedVec.push_back(SeedPairVec[k]);
				j = k;
			}
			else break;
		}
		if (AlignmentCandidate.Score > thr)
		{
			if ((AlignmentCandidate.PosDiff = AlignmentCandidate.SeedVec[0].PosDiff) < 0) AlignmentCandidate.PosDiff = 0;
			//sort(AlignmentCandidate.SeedVec.begin(), AlignmentCandidate.SeedVec.end(), CompByGenomePos);
			AlignmentVec.push_back(AlignmentCandidate);
			//if (bDebugMode)
			//{
			//	printf("Candidate score = %d\n", AlignmentCandidate.Score);
			//	ShowSeedInfo(AlignmentCandidate.SeedVec);
			//}
		}
		i = k;
	}
	return AlignmentVec;
}

void RemoveShortSeeds(vector<SeedPair_t>& SeedVec, int thr)
{
	for (vector<SeedPair_t>::iterator iter = SeedVec.begin(); iter != SeedVec.end();)
	{
		if (iter->rLen <= thr) iter = SeedVec.erase(iter);
		else iter++;
	}
}

void RemoveNullSeeds(vector<SeedPair_t>& SeedVec)
{
	for (vector<SeedPair_t>::iterator iter = SeedVec.begin(); iter != SeedVec.end();)
	{
		if (iter->rLen == 0) iter = SeedVec.erase(iter);
		else iter++;
	}
}

void ShowFragmentPair(char* seq, SeedPair_t& SeedPair)
{
	int i, glen, shift;
	string frag1, frag2;
	int64_t LeftGPos, RightGPos;

	printf("r[%d-%d], g[%ld-%ld], Len=%d\n", SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.rLen);

	LeftGPos = SeedPair.gPos - 50;
	RightGPos = SeedPair.gPos + SeedPair.gLen + 50;

	if (LeftGPos < 0) LeftGPos = 0; glen = RightGPos - LeftGPos + 1; shift = SeedPair.gPos - LeftGPos;
	frag1.resize(SeedPair.rLen); frag2.resize(glen);
	strncpy((char*)frag1.c_str(), seq + SeedPair.rPos, SeedPair.rLen);
	strncpy((char*)frag2.c_str(), RefSequence + LeftGPos, glen);
	printf("%s%s\n%s\n\n", string().assign(shift, ' ').c_str(), frag1.c_str(), frag2.c_str());
}

void ShowAlnStatus(string& frag1, string& frag2, string& aln)
{
	int i, len = (int)frag1.length(); aln.resize(len);

	for (i = 0; i < len; i++)
	{
		if (frag1[i] == frag2[i]) aln[i] = ' ';
		else aln[i] = '*';
	}
}

pair<int, int> IdentifyBestUnGappedPartition(char* seq, int rGaps, SeedPair_t& LeftSeed, SeedPair_t& RightSeed)
{
	int64_t gPos;
	vector<int> Rvec, Lvec;
	int i, p, s, rPos, max_score = 0;

	Rvec.resize(rGaps + 1); s = 0; rPos = LeftSeed.rPos + LeftSeed.rLen; gPos = LeftSeed.gPos + LeftSeed.gLen;
	for (i = 1; i <= rGaps; rPos++, gPos++, i++)
	{
		if (nst_nt4_table[(int)seq[rPos]] == nst_nt4_table[(int)RefSequence[gPos]]) s++;
		Rvec[i] = s; // ext_len=i, score=s
	}
	//string frag1, frag2, aln;
	//frag1.resize(rGaps); strncpy((char*)frag1.c_str(), seq + LeftSeed.rPos + LeftSeed.rLen, rGaps);
	//frag2.resize(rGaps); strncpy((char*)frag2.c_str(), RefSequence + LeftSeed.gPos + LeftSeed.gLen, rGaps);
	//ShowAlnStatus(frag1, frag2, aln);
	//printf("Extension -->\n %s\n %s\n %s\n", frag1.c_str(), frag2.c_str(), aln.c_str());
	//for (i = 0; i <= rGaps; i++) printf("%d", Rvec[i] % 10);
	//printf("\n\n");

	Lvec.resize(rGaps + 1); s = 0; rPos = RightSeed.rPos - 1; gPos = RightSeed.gPos - 1;
	for (i = 1; i <= rGaps; rPos--, gPos--, i++)
	{
		if (nst_nt4_table[(int)seq[rPos]] == nst_nt4_table[(int)RefSequence[gPos]]) s++;
		Lvec[rGaps - i] = s; // ext_len=i, score=s
	}
	//strncpy((char*)frag1.c_str(), seq + RightSeed.rPos - rGaps, rGaps);
	//strncpy((char*)frag2.c_str(), RefSequence + RightSeed.gPos - rGaps, rGaps);
	//ShowAlnStatus(frag1, frag2, aln);
	//printf("Extension <--\n%s\n%s\n%s\n", frag1.c_str(), frag2.c_str(), aln.c_str());
	//for (i = 0; i <= rGaps; i++) printf("%d", Lvec[i] % 10);
	//printf("\n\n");

	for (i = p = 0; i <= rGaps; i++)
	{
		s = Rvec[i] + Lvec[i];
		//if (s > max_score || (s == max_score && abs(rGaps - p*2) < abs(rGaps - i*2)))
		if (s > max_score)
		{
			max_score = s;
			p = i;
		}
	}	
	//printf("BestPartition s= %d / %d, score=%d\n\n", p, rGaps, max_score);

	return make_pair(p, max_score);
}

GappedExtension_t IdentifyBestGappedPartition(char* seq, int rGaps, SeedPair_t& LeftSeed, SeedPair_t& RightSeed)
{
	int64_t gPos;
	vector<int> Rvec, Lvec;
	GappedExtension_t GappedExtension;
	int i, len, p, q, s, max_score = 0;
	string frag1, frag2, frag3, frag4, aln;

	frag1.resize(rGaps); strncpy((char*)frag1.c_str(), seq + LeftSeed.rPos + LeftSeed.rLen, rGaps);
	frag2.resize(rGaps); strncpy((char*)frag2.c_str(), RefSequence + LeftSeed.gPos + LeftSeed.gLen, rGaps);
	PairwiseSequenceAlignment(rGaps, frag1, rGaps, frag2);
	// identify tailing gaps
	len = (int)frag1.length(); i = len - 1; while (frag2[i] == '-') i--;
	for (i += 1, gPos = LeftSeed.gPos + LeftSeed.gLen + rGaps; i < len; i++,gPos++) frag2[i] = RefSequence[gPos];

	Rvec.resize(rGaps + 1);
	for (p = s = 0, i = 0; i < len; i++)
	{
		if (frag1[i] == frag2[i]) s++;
		//else if (s > 0 && (frag1[i] == '-' || frag2[i] == '-')) s--;
		if (frag1[i] != '-') p++;
		Rvec[p] = s; // ext_len = p, score=s
	}
	//ShowAlnStatus(frag1, frag2, aln); printf("Gapped alignment -->\n%s\n%s\n%s\n", frag1.c_str(), frag2.c_str(), aln.c_str());
	//for (p = i = 0; i < len; i++)
	//{
	//	if (frag1[i] != '-') p++;
	//	printf("%d", Rvec[p] % 10);
	//}
	//printf("\n\n");

	frag3.resize(rGaps); strncpy((char*)frag3.c_str(), seq + LeftSeed.rPos + LeftSeed.rLen, rGaps);
	frag4.resize(rGaps); strncpy((char*)frag4.c_str(), RefSequence + RightSeed.gPos - rGaps, rGaps);
	PairwiseSequenceAlignment(rGaps, frag3, rGaps, frag4);
	// identify heading gaps
	i = 0; while (frag4[i] == '-') i++;
	for (i -= 1, gPos = RightSeed.gPos - rGaps; i >= 0;i--,gPos--) frag4[i] = RefSequence[gPos];

	len = (int)frag3.length(); Lvec.resize(rGaps + 1);
	for (p = 0, s = 0, i = len - 1; i >= 0; i--)
	{
		if (frag3[i] == frag4[i]) s++;
		//else if (s > 0 && (frag3[i] == '-' || frag4[i] == '-')) s--;
		if (frag3[i] != '-') p++;
		Lvec[rGaps - p] = s; // ext_len = p, score = s
	}
	//ShowAlnStatus(frag3, frag4, aln); printf("Gapped alignment <--\n%s\n%s\n%s\n", frag3.c_str(), frag4.c_str(), aln.c_str());
	//for (p = i = 0; i < len; i++)
	//{
	//	printf("%d", Lvec[p] % 10);
	//	if (frag3[i] != '-') p++;
	//}
	//printf("\n\n");

	for (GappedExtension.p = i = 0; i <= rGaps; i++)
	{
		s = Rvec[i] + Lvec[i];
		if (s > max_score)
		{
			max_score = s;
			GappedExtension.p = i;
		}
	}
	for (GappedExtension.right_ext = 0, p = GappedExtension.p, i = 0; p > 0; i++)
	{
		if (frag1[i] != '-') p--;
		if (frag2[i] != '-') GappedExtension.right_ext++;
	}
	for (GappedExtension.left_ext = 0, p = rGaps - GappedExtension.p, i = (int)frag3.length() - 1; p > 0; i--)
	{
		if (frag3[i] != '-') p--;
		if (frag4[i] != '-') GappedExtension.left_ext++;
	}
	return GappedExtension;
}

void FillGapsBetweenAdjacentSeeds(char* seq, SeedPair_t& left_seed, SeedPair_t& right_seed, vector<SeedPair_t>& Vec)
{
	int64_t gPos;
	SeedPair_t seed;
	int rPos, rGaps;
	pair<int, int> Partition;

	seed.bAcceptorSite = false; seed.rLen = seed.gLen = 0;

	rGaps = right_seed.rPos - (left_seed.rPos + left_seed.rLen);
	//if (bDebugMode)
	//{
	//	printf("\n\nrGaps=%d\n", rGaps);
	//	ShowFragmentPair(seq, left_seed); ShowFragmentPair(seq, right_seed);
	//}
	Partition = IdentifyBestUnGappedPartition(seq, rGaps, left_seed, right_seed);
	if (rGaps < 5 || Partition.second >= (int)ceil(rGaps*0.8))
	{
		// add one or two normal pairs in between
		seed.bSimple = false;
		if (Partition.first > 0)
		{
			left_seed.rLen += Partition.first;
			left_seed.gLen += Partition.first;
			//if (bDebugMode) printf("Case1.1:\nr[%d-%d]=%d, g[%ld-%ld]=%d, PosDiff=%d\n", left_seed.rPos, left_seed.rPos + left_seed.rLen - 1, left_seed.rLen, left_seed.gPos, left_seed.gPos + left_seed.gLen - 1, left_seed.gLen, left_seed.PosDiff);
		}
		if ((rGaps -= Partition.first) > 0)
		{
			right_seed.rPos -= rGaps; right_seed.rLen += rGaps;
			right_seed.gPos -= rGaps; right_seed.gLen += rGaps;
			//if (bDebugMode) printf("Case1.2:\nr[%d-%d]=%d, g[%ld-%ld]=%d, PosDiff=%d\n", right_seed.rPos, right_seed.rPos + right_seed.rLen - 1, right_seed.rLen, right_seed.gPos, right_seed.gPos + right_seed.gLen - 1, right_seed.gLen, right_seed.PosDiff);
		}
	}
	else
	{
		GappedExtension_t GappedExtension = IdentifyBestGappedPartition(seq, rGaps, left_seed, right_seed);
		// add one or two normal pairs in between
		seed.bSimple = false;
		//if (bDebugMode) printf("Case result: GappedExtension.p=%d, left_ext=%d, right_ext=%d\n", GappedExtension.p, GappedExtension.left_ext, GappedExtension.right_ext);
		if (GappedExtension.p > 0)
		{
			if (GappedExtension.p == GappedExtension.right_ext)
			{
				left_seed.rLen += GappedExtension.p;
				left_seed.gLen += GappedExtension.p;
				//if (bDebugMode) printf("Case2.1:\nr[%d-%d]=%d, g[%ld-%ld]=%d, PosDiff=%d\n", left_seed.rPos, left_seed.rPos + left_seed.rLen - 1, left_seed.rLen, left_seed.gPos, left_seed.gPos + left_seed.gLen - 1, left_seed.gLen, left_seed.PosDiff);
			}
			else
			{
				seed.rPos = left_seed.rPos + left_seed.rLen;
				seed.gPos = left_seed.gPos + left_seed.gLen;
				seed.PosDiff = seed.gPos - seed.rPos;
				seed.rLen = GappedExtension.p;
				seed.gLen = GappedExtension.right_ext;
				Vec.push_back(seed);
				//if (bDebugMode) printf("Case2.1 normal pairs between exons:\nr[%d-%d]=%d, g[%ld-%ld]=%d\n", seed.rPos, seed.rPos + seed.rLen - 1, seed.rLen, seed.gPos, seed.gPos + seed.gLen - 1, seed.gLen);
			}
		}
		if ((rGaps -= GappedExtension.p) > 0)
		{
			if (rGaps == GappedExtension.left_ext)
			{
				right_seed.rPos -= rGaps; right_seed.rLen += rGaps;
				right_seed.gPos -= rGaps; right_seed.gLen += rGaps;
				//if (bDebugMode) printf("Case2_2:\nr[%d-%d]=%d, g[%ld-%ld]=%d, PosDiff=%d\n", right_seed.rPos, right_seed.rPos + right_seed.rLen - 1, right_seed.rLen, right_seed.gPos, right_seed.gPos + right_seed.gLen - 1, right_seed.gLen, right_seed.PosDiff);
			}
			else
			{
				seed.rLen = rGaps;
				seed.gLen = GappedExtension.left_ext;
				seed.rPos = right_seed.rPos - seed.rLen;
				seed.gPos = right_seed.gPos - seed.gLen;
				seed.PosDiff = seed.gPos - seed.rPos;
				Vec.push_back(seed);
				//if (bDebugMode) printf("Case2_2 normal pairs between exons:\nr[%d-%d]=%d, g[%ld-%ld]=%d\n", seed.rPos, seed.rPos + seed.rLen - 1, seed.rLen, seed.gPos, seed.gPos + seed.gLen - 1, seed.gLen);
			}
		}
	}
}

void SeedExtension(char* seq, vector<SeedPair_t>& SeedVec)
{
	int64_t gPos;
	SeedPair_t seed;
	vector<SeedPair_t> Vec;
	pair<int, int> Partition;
	int i, rPos, rGaps, PosDiff, num = (int)SeedVec.size();

	for (i = 1; i < num; i++)
	{
		if ((PosDiff = (int)(SeedVec[i].PosDiff - SeedVec[i - 1].PosDiff)) > MinIntronSize && SeedVec[i].rPos > (SeedVec[i - 1].rPos + SeedVec[i - 1].rLen))
			FillGapsBetweenAdjacentSeeds(seq, SeedVec[i - 1], SeedVec[i], Vec);
	}
	if (Vec.size() > 0)
	{
		//printf("\nadd more seed by SeedExtension\n\n");
		copy(Vec.begin(), Vec.end(), back_inserter(SeedVec));
		sort(SeedVec.begin(), SeedVec.end(), CompByGenomePos);
	}
}

SeedPair_t ReseedingWithSpecificRegion(char* seq, int rBegin, int rEnd, int64_t L_Boundary, int64_t R_Boundary)
{
	SeedPair_t seed;
	char *frag1, *frag2;
	vector<SeedPair_t> vec;
	int i, j, rlen, glen, thr;
	int64_t gPos, PosDiff, best_gPos;

	//if (bDebugMode) printf("Perform re-seeding between r[%d-%d] g[%ld - %ld]\n", rBegin, rEnd - 1, L_Boundary, R_Boundary);

	rlen = rEnd - rBegin; glen = R_Boundary - L_Boundary;
	frag1 = new char[rlen + 1]; strncpy(frag1, seq + rBegin, rlen);
	frag2 = new char[glen + 1]; strncpy(frag2, RefSequence + L_Boundary, glen);

	if ((thr = (int)(rlen*0.85)) < KmerSize) thr = KmerSize;
	seed = GenerateLongestSimplePairsFromFragmentPair(rlen, frag1, glen, frag2);
	if(seed.rLen >= thr)
	{
		seed.rPos += rBegin;
		seed.gPos += L_Boundary;
		seed.PosDiff = seed.gPos - seed.rPos;
		//printf("gap seed r[%d-%d] g[%ld-%ld], len=%d\n", seed.rPos, seed.rPos + seed.rLen - 1, seed.gPos, seed.gPos + seed.gLen - 1, seed.rLen);
	}
	else seed.rLen = 0;
	//else printf("no any seed found!\n\n");

	delete[] frag1; delete[] frag2;

	return seed;
}

SeedPair_t ReseedingFromEnd(bool bDir, uint8_t* seq, int rBegin, int rEnd, int64_t OriginalGenomicPos)
{
	SeedPair_t seed;
	int i, j, rlen, glen, thr;
	map<int64_t, int>::iterator ChrIter;
	int64_t gPos, gStop, L_Boundary, R_Boundary, PosDiff, best_gPos;

	//if (bDebugMode) printf("Perform re-seeding with r[%d-%d]\n", rBegin, rEnd);

	seed.rLen = seed.gLen = 0;

	ChrIter = ChrLocMap.lower_bound(OriginalGenomicPos);
	if (bDir) gStop = ChromosomeVec[ChrIter->second].FowardLocation + ChromosomeVec[ChrIter->second].len;
	else gStop = ChromosomeVec[ChrIter->second].FowardLocation;

	if (bDir) // -->
	{
		L_Boundary = OriginalGenomicPos;
		if ((R_Boundary = L_Boundary + MaxIntronSize) >= gStop) R_Boundary = gStop;
	}
	else // <--
	{
		R_Boundary = OriginalGenomicPos;
		if ((L_Boundary = R_Boundary - MaxIntronSize) < gStop) L_Boundary = gStop;
	}
	//if (bDebugMode) printf("Search area: %ld-%ld\n", L_Boundary, R_Boundary), fflush(stdout);

	bwtSearchResult_t bwtSearchResult = BWT_LocalSearch(seq, rBegin, rEnd, L_Boundary, R_Boundary);
	if (bwtSearchResult.freq > 0 && (seed.PosDiff = seed.gPos - seed.rPos) < MaxIntronSize)
	{
		seed.rPos = 0;
		seed.gPos = bwtSearchResult.LocArr[0];
		seed.rLen = seed.gLen = bwtSearchResult.len;
		//if (bDebugMode) printf("Re-seeding result: r[%d-%d] g[%ld-%ld], len=%d\n", seed.rPos, seed.rPos + seed.rLen - 1, seed.gPos, seed.gPos + seed.gLen - 1, seed.rLen), fflush(stdout);

		delete[] bwtSearchResult.LocArr;
	}
	return seed;
}

SeedPair_t IdentifyHeadingSeed(char* seq, uint8_t* EncodeSeq, int nGaps, int64_t gPos)
{
	SeedPair_t seed;
	seed.rLen = seed.gLen = 0;

	if(nGaps >= 10) seed = ReseedingFromEnd(false, EncodeSeq, 0, nGaps, gPos);

	return seed;
}

SeedPair_t IdentifyTailingSeed(char* seq, uint8_t* EncodeSeq, int rPos_begin, int rPos_end, int64_t gPos)
{
	int nGaps;
	SeedPair_t seed;

	seed.rLen = seed.gLen = 0;
	if(rPos_end - rPos_begin >= 10) seed = ReseedingFromEnd(true, EncodeSeq, rPos_begin, rPos_end, gPos);

	return seed;
}

void IdentifyMissingSeeds(int rlen, char* seq, vector<SeedPair_t>& SeedVec)
{
	int64_t gPos;
	SeedPair_t seed;
	pair<int, int> Partition;
	int i, rPos, rGaps, shift, PosDiff, num = (int)SeedVec.size();

	for (i = 1; i < num; i++)
	{
		if ((PosDiff = (int)(SeedVec[i].PosDiff - SeedVec[i - 1].PosDiff)) > MaxGaps && (rGaps = SeedVec[i].rPos - SeedVec[i - 1].rPos - SeedVec[i - 1].rLen) > 20)
		{
			seed = ReseedingWithSpecificRegion(seq, SeedVec[i - 1].rPos + SeedVec[i - 1].rLen, SeedVec[i].rPos, SeedVec[i - 1].gPos + SeedVec[i - 1].gLen, SeedVec[i].gPos);
			if (seed.rLen > 0) SeedVec.push_back(seed);
		}
	}
	if ((int)SeedVec.size() > num) sort(SeedVec.begin(), SeedVec.end(), CompByGenomePos);
}

bool CheckSeqFragment(int64_t LeftGPos, int64_t RightGPos, int shift)
{
	int i;
	bool bIdentical = true;

	if (shift > 0)
	{
		for (i = 0; i < shift; i++, LeftGPos++, RightGPos++)
		{
			if (RefSequence[LeftGPos] != RefSequence[RightGPos])
			{
				bIdentical = false;
				break;
			}
		}
	}
	else
	{
		for (shift = 0 - shift, LeftGPos -= shift, RightGPos-=shift, i = 0; i < shift; i++, LeftGPos++, RightGPos++)
		{
			if (RefSequence[LeftGPos] != RefSequence[RightGPos])
			{
				bIdentical = false;
				break;
			}
		}
	}
	return bIdentical;
}

int IdentifySpliceJunction(int SJtype, SeedPair_t& left_seed, SeedPair_t& right_seed)
{
	int i, j, shift;
	bool bChecked = false;
	int64_t LeftGPos, RightGPos, gPos1, gPos2;

	i = (left_seed.rLen < right_seed.rLen ? left_seed.rLen : right_seed.rLen);
	j = (left_seed.gLen < right_seed.gLen ? left_seed.gLen : right_seed.gLen);

	if (i < j) j = i;
	if (j > 9) j = 9; j <<= 1;

	LeftGPos = left_seed.gPos + left_seed.gLen;
	RightGPos = right_seed.gPos;

	for (i = 0; i <= j; i++)
	{
		shift = ShiftArr[i];
		if (shift != 0 && CheckSeqFragment(LeftGPos, RightGPos, shift) == false) continue;

		gPos1 = LeftGPos + shift; gPos2 = RightGPos - 2 + shift;
		if (RefSequence[gPos1] == SpliceJunctionArr[SJtype][0] && RefSequence[gPos1 + 1] == SpliceJunctionArr[SJtype][1] && RefSequence[gPos2] == SpliceJunctionArr[SJtype][3] && RefSequence[gPos2 + 1] == SpliceJunctionArr[SJtype][4]) break;
	}
	if (i > j) return 10;
	else return shift;
}

int CheckSpliceJunction(int rlen, char* seq, uint8_t* EncodeSeq, vector<SeedPair_t>& SeedVec)
{
	int64_t gPos;
	SeedPair_t seed;
	vector<SeedPair_t> Vec;
	vector<pair<int, int> > vec, best_vec;
	int i, j, rPos, SJtype, shift, c, mis, min_cost, best_type, num = (int)SeedVec.size();

	min_cost = 1000; best_type = -1;
	for (SJtype = 0; SJtype < 4; SJtype++)
	{
		vec.clear(); mis = 0;
		for (c = 0, i = 1; i < num; i++)
		{
			if ((SeedVec[i].PosDiff - SeedVec[i - 1].PosDiff) > MinIntronSize && SeedVec[i-1].bSimple && SeedVec[i].bSimple)
			{
				shift = IdentifySpliceJunction(SJtype, SeedVec[i - 1], SeedVec[i]);
				if (shift != 10) vec.push_back(make_pair(i, shift));
				else mis++;
				c += abs(shift);
			}
		}
		if (vec.size() > 0 && c < min_cost)
		{
			min_cost = c;
			best_type = SJtype;
			best_vec = vec;
		}
		if (mis == 0) break;
	}
	if (best_type != -1)
	{
		for (i = 0; i < (int)best_vec.size(); i++)
		{
			j = best_vec[i].first;
			shift = best_vec[i].second;

			if (shift != 10)
			{
				SeedVec[j].bAcceptorSite = true;
				if (shift != 0)
				{
					SeedVec[j - 1].rLen += shift; SeedVec[j - 1].gLen += shift;
					SeedVec[j].rLen -= shift; SeedVec[j].gLen -= shift;
					SeedVec[j].rPos += shift; SeedVec[j].gPos += shift;
				}
			}
		}
	}
	//if (bDebugMode)
	//{
	//	if (best_type != -1)
	//	{
	//		printf("Splice junction type = %s\n", SpliceJunctionArr[best_type]);
	//		for (num = (int)SeedVec.size(), i = 0; i < num; i++) ShowFragmentPair(seq, SeedVec[i]);
	//	}
	//	else printf("Cannot find SJ!!\n");
	//}
	return best_type;
}

void RemoveTandemRepeatSeeds(vector<SeedPair_t>& SeedVec)
{
	int i, j, k, num;
	bool bTandemRepeat = false;
	vector<pair<int, int> > vec;

	num = (int)SeedVec.size(); if (num < 2) return;
	vec.resize(num);
	for (i = 0; i < num; i++) vec[i] = make_pair(SeedVec[i].rPos, i);
	sort(vec.begin(), vec.end(), CompByFirstInt); //first:rPos, second: gPos rank

	for (i = 0; i < num;) // identify all repetitive rPos
	{
		j = i + 1; while (j < num && vec[j].first == vec[i].first) j++;
		if (j - i > 1)
		{
			bTandemRepeat = true;
			//printf("Tandem repeat found: rPos = %d\n", vec[i].first);
			for (k = i; k < j; k++) SeedVec[vec[k].second].rLen = SeedVec[vec[k].second].gLen = 0;
		}
		i = j;
	}
	//vec.clear(); vector<pair<int, int> >().swap(vec);
	if (bTandemRepeat) RemoveNullSeeds(SeedVec);
	//if (bDebugMode) printf("After tandem repeat checking\n"), ShowSeedInfo(SeedVec);
}

int IdentifyTranslocationRange(int i, int num, vector<pair<int, int> >& vec, vector<SeedPair_t>& SeedVec)
{
	int j, max_idx = vec[i].second;

	for (j = i + 1; j <= max_idx; j++)
	{
		if (vec[j].second > max_idx) max_idx = vec[j].second;
	}
	return max_idx;
}

void RemoveTranslocatedSeeds(vector<SeedPair_t>& SeedVec)
{
	int64_t gPos;
	int i, j, k, s1, s2, num;
	bool bTranslocation = false;
	vector<pair<int, int> > vec;

	num = (int)SeedVec.size(); if (num < 2) return;
	vec.resize(num);
	for (i = 0; i < num; i++) vec[i] = make_pair(SeedVec[i].rPos, i);
	sort(vec.begin(), vec.end(), CompByFirstInt); //first:rPos, second: gPos rank

	// checking translocation
	for (i = 0; i < num; i++)
	{
		if (vec[i].first != SeedVec[i].rPos)
		{
			bTranslocation = true;
			j = IdentifyTranslocationRange(i, num, vec, SeedVec);
			s1 = 0; s2 = 0;
			for (k = i; k <= j; k++)
			{
				if (k < vec[k].second) s1 += SeedVec[vec[k].second].rLen;
				else s2 += SeedVec[vec[k].second].rLen;
			}
			if (s1 > s2)
			{
				for (k = i; k <= j; k++)
					if (k > vec[k].second)
					{
						//printf("set #%d seed to zero length\n", vec[k].second + 1);
						SeedVec[vec[k].second].rLen = SeedVec[vec[k].second].gLen = 0;
					}
			}
			else
			{
				for (k = i; k <= j; k++)
					if (k < vec[k].second)
					{
						//printf("set #%d seed to zero length\n", vec[k].second + 1);
						SeedVec[vec[k].second].rLen = SeedVec[vec[k].second].gLen = 0;
					}
			}
			i = j;
		}
	}
	if (bTranslocation) RemoveNullSeeds(SeedVec);
	//if (bDebugMode) printf("After translocation checking\n"), ShowSeedInfo(SeedVec);
}

bool CheckSeedOverlapping(SeedPair_t& p1, SeedPair_t& p2)
{
	int iOverlap;
	bool bMaster = true;

	//printf("[1]: p1: r[%d-%d]=%d g[%ld-%ld]=%d vs p2: r[%d-%d]=%d g[%ld-%ld]=%d\n", p1.rPos, p1.rPos + p1.rLen - 1, p1.rLen, p1.gPos, p1.gPos + p1.gLen - 1, p1.gLen, p2.rPos, p2.rPos + p2.rLen - 1, p2.rLen, p2.gPos, p2.gPos + p2.gLen - 1, p2.gLen);
	if ((iOverlap = p1.rPos + p1.rLen - p2.rPos) > 0)
	{
		if (p1.rLen < p2.rLen)
		{
			bMaster = false;
			if (p1.rLen > iOverlap)
			{
				p1.gLen = (p1.rLen -= iOverlap);
			}
			else p1.rLen = p1.gLen = 0;
		}
		else
		{
			if (p2.rLen > iOverlap)
			{
				p2.rPos += iOverlap; p2.gPos += iOverlap;
				p2.gLen = (p2.rLen -= iOverlap);
			}
			else p2.rLen = p2.gLen = 0;
		}
	}
	if ((p1.rLen > 0 && p2.rLen > 0) && (iOverlap = p1.gPos + p1.gLen - p2.gPos) > 0)
	{
		if (p1.gLen < p2.gLen)
		{
			bMaster = false;
			if (p1.rLen > iOverlap)
			{
				p1.gLen = (p1.rLen -= iOverlap);
			}
			else p1.rLen = p1.gLen = 0;
		}
		else
		{
			if (p2.rLen > iOverlap)
			{
				p2.rPos += iOverlap; p2.gPos += iOverlap;
				p2.gLen = (p2.rLen -= iOverlap);
			}
			else p2.rLen = p2.gLen = 0;
		}
	}
	//printf("[2]: p1: r[%d-%d]=%d g[%ld-%ld]=%d vs p2: r[%d-%d]=%d g[%ld-%ld]=%d\n\n", p1.rPos, p1.rPos + p1.rLen - 1, p1.rLen, p1.gPos, p1.gPos + p1.gLen - 1, p1.gLen, p2.rPos, p2.rPos + p2.rLen - 1, p2.rLen, p2.gPos, p2.gPos + p2.gLen - 1, p2.gLen);
	return bMaster;
}

int LocateThePreviousSeedIdx(int i, vector<SeedPair_t>& SeedVec)
{
	while (i > 0 && SeedVec[i].rLen == 0) i--;
	if (i < 0) return 0;
	else return i;
}

void CheckOverlappingSeeds(vector<SeedPair_t>& SeedVec)
{
	int64_t gEnd;
	int i, j, num, rEnd;
	bool bNullSeed = false;

	num = (int)SeedVec.size(); if (num < 2) return;
	for (i = 0; i < num;)
	{
		if (SeedVec[i].rLen > 0)
		{
			rEnd = SeedVec[i].rPos + SeedVec[i].rLen - 1;
			gEnd = SeedVec[i].gPos + SeedVec[i].gLen - 1;

			//printf("overlap check for seed#%d, rEnd=%d, gEned=%ld\n", i + 1, rEnd, gEnd);
			for (j = i + 1; j < num; j++)
			{
				if (SeedVec[j].rLen == 0) continue;
				if (rEnd < SeedVec[j].rPos && gEnd < SeedVec[j].gPos) break;
				//printf("\ttest seed#%d\n", j + 1);
				if (CheckSeedOverlapping(SeedVec[i], SeedVec[j]) == false) break;
			}
			if (SeedVec[i].rLen == 0)
			{
				bNullSeed = true;
				i = LocateThePreviousSeedIdx(i - 1, SeedVec);
			}
			else i++;
		}
		else
		{
			bNullSeed = true; i++;
		}
	}
	if (bNullSeed) RemoveNullSeeds(SeedVec);
	//if (bDebugMode) printf("after overlap checking\n"), ShowSeedInfo(SeedVec);
}

void IdentifyNormalPairs(int rlen, char* seq, vector<SeedPair_t>& SeedVec)
{
	SeedPair_t SeedPair;
	int i, j, rGaps, gGaps, num;

	SeedPair.bAcceptorSite = SeedPair.bSimple = false;

	if (SeedVec.size() > 1)
	{
		CheckOverlappingSeeds(SeedVec);

		num = (int)SeedVec.size();
		for (i = 0, j = 1; j < num; i++, j++)
		{
			//if (abs(SeedVec[j].PosDiff - SeedVec[i].PosDiff) > MaxGaps) continue;
			//if ((SeedVec[j].PosDiff - SeedVec[i].PosDiff) > MaxGaps) continue;
			if (SeedVec[j].rPos - SeedVec[i].rPos - SeedVec[i].rLen == 0) continue;

			rGaps = SeedVec[j].rPos - (SeedVec[i].rPos + SeedVec[i].rLen); if (rGaps < 0) rGaps = 0;
			gGaps = SeedVec[j].gPos - (SeedVec[i].gPos + SeedVec[i].gLen); if (gGaps < 0) gGaps = 0; else if (gGaps > 30 && gGaps > (rGaps<<1)) gGaps = 0;

			if (rGaps > 0 || gGaps > 0)
			{
				SeedPair.rPos = SeedVec[i].rPos + SeedVec[i].rLen;
				SeedPair.gPos = SeedVec[i].gPos + SeedVec[i].gLen;
				SeedPair.PosDiff = SeedPair.gPos - SeedPair.rPos;
				SeedPair.rLen = rGaps; SeedPair.gLen = gGaps;
				
				SeedVec.push_back(SeedPair);
				//if (bDebugMode) printf("normal pair:\nR1[%d-%d] G1[%ld-%ld]\nNR[%d-%d]=%d NG[%ld-%ld]=%d\nR2[%d-%d] G2[%ld-%ld]\n", SeedVec[i].rPos, SeedVec[i].rPos + SeedVec[i].rLen - 1, SeedVec[i].gPos, SeedVec[i].gPos + SeedVec[i].gLen - 1, SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen, SeedVec[j].rPos, SeedVec[j].rPos + SeedVec[j].rLen - 1, SeedVec[j].gPos, SeedVec[j].gPos + SeedVec[j].gLen - 1); fflush(stdout);
			}
		}
		if ((int)SeedVec.size() > num) inplace_merge(SeedVec.begin(), SeedVec.begin()+num, SeedVec.end(), CompByGenomePos);
	}
}

void UpdateMyExonMap(map<int64_t, int>& MyExonMap, vector<SeedPair_t>& SeedPairVec)
{
	int i, tail_idx, num = (int)SeedPairVec.size();

	for (tail_idx = num - 1, i = 1; i < num; i++)
	{
		if (SeedPairVec[i].bAcceptorSite)
		{
			//if (bDebugMode) printf("Acceptor site found: %ld\n", SeedPairVec[i].gPos);
			if(i != tail_idx) MyExonMap[SeedPairVec[i].gPos] = SeedPairVec[i].gLen;
			else MyExonMap[SeedPairVec[i].gPos] = 0;
		}
	}
}

void GenMappingReport(bool bFirstRead, ReadItem_t& read, vector<AlignmentCandidate_t>& AlignmentVec)
{
	Coordinate_t coor;
	int i, j, k, g, num;
	map<int64_t, int> MyExonMap;
	vector<pair<int, char> > cigar_vec;
	map<int64_t, int>::iterator iter, ExonMapIter;

	if (bDebugMode) printf("\n\nGenerate alignment for read %s (%d cans)\n", read.header, (int)AlignmentVec.size()), fflush(stdout);
	
	read.score = read.iBestAlnCanIdx = 0;
	if ((read.CanNum = (int)AlignmentVec.size()) > 0)
	{
		read.AlnReportArr = new AlignmentReport_t[read.CanNum];
		for (i = 0; i != (int)AlignmentVec.size(); i++)
		{
			read.AlnReportArr[i].AlnScore = 0;
			read.AlnReportArr[i].PairedAlnCanIdx = AlignmentVec[i].PairedAlnCanIdx;

			if (AlignmentVec[i].Score == 0) continue;
			//if (bDebugMode)
			//{
			//	printf("Original seeds\n");
			//	ShowSeedInfo(AlignmentVec[i].SeedVec);
			//}
			RemoveTandemRepeatSeeds(AlignmentVec[i].SeedVec);
			RemoveTranslocatedSeeds(AlignmentVec[i].SeedVec);
			IdentifyMissingSeeds(read.rlen, read.seq, AlignmentVec[i].SeedVec);
			//printf("After IdentifyMissingSeeds\n"), ShowSeedInfo(AlignmentVec[i].SeedVec);
			SeedExtension(read.seq, AlignmentVec[i].SeedVec);
			//printf("After SeedExtension\n"), ShowSeedInfo(AlignmentVec[i].SeedVec);
			AlignmentVec[i].SJtype = CheckSpliceJunction(read.rlen, read.seq, read.EncodeSeq, AlignmentVec[i].SeedVec);
			//printf("After CheckSpliceJunction\n"), ShowSeedInfo(AlignmentVec[i].SeedVec);
			IdentifyNormalPairs(read.rlen, read.seq, AlignmentVec[i].SeedVec); // fill missing framgment pairs between simple pairs
			//printf("After IdentifyNormalPairs\n"), ShowSeedInfo(AlignmentVec[i].SeedVec);
			//if (bDebugMode)
			//{
			//	printf("Process candidate#%d (Score = %d, SegmentPair#=%d): \n", i + 1, AlignmentVec[i].Score, (int)AlignmentVec[i].SeedVec.size());
			//	ShowSeedInfo(AlignmentVec[i].SeedVec);
			//}
			num = (int)AlignmentVec[i].SeedVec.size();
			if (num > 1 && CheckCoordinateValidity(AlignmentVec[i].SeedVec) == false)
			{
				//if (bDebugMode) printf("CheckCoordinateValidity fails!\n");
				continue;
			}
			cigar_vec.clear();
			for (j = 0; j != num; j++)
			{
				if (AlignmentVec[i].SeedVec[j].rLen == 0 && AlignmentVec[i].SeedVec[j].gLen == 0) continue;
				else
				{
					if (j > 0 && (g = AlignmentVec[i].SeedVec[j].gPos - (AlignmentVec[i].SeedVec[j - 1].gPos + AlignmentVec[i].SeedVec[j - 1].gLen)) > 0) cigar_vec.push_back(make_pair(g, 'N'));

					if (AlignmentVec[i].SeedVec[j].bSimple)
					{
						//if (bDebugMode) ShowFragmentPair(AlignmentVec[i].SeedVec[j]);
						cigar_vec.push_back(make_pair(AlignmentVec[i].SeedVec[j].rLen, 'M'));
						read.AlnReportArr[i].AlnScore += AlignmentVec[i].SeedVec[j].rLen;
					}
					else
					{
						//if (bDebugMode) printf("Check normal pair#%d: R[%d-%d]=%d G[%ld-%ld]=%d\n", j + 1, AlignmentVec[i].SeedVec[j].rPos, AlignmentVec[i].SeedVec[j].rPos + AlignmentVec[i].SeedVec[j].rLen - 1, AlignmentVec[i].SeedVec[j].rLen, AlignmentVec[i].SeedVec[j].gPos, AlignmentVec[i].SeedVec[j].gPos + AlignmentVec[i].SeedVec[j].gLen - 1, AlignmentVec[i].SeedVec[j].gLen), fflush(stdout);
						if (j == 0)
						{
							read.AlnReportArr[i].AlnScore += ProcessHeadSequencePair(read.seq, AlignmentVec[i].SeedVec[0], cigar_vec);
						}
						else if (j == num - 1)
						{
							read.AlnReportArr[i].AlnScore += ProcessTailSequencePair(read.seq, AlignmentVec[i].SeedVec[j], cigar_vec);
						}
						else
						{
							read.AlnReportArr[i].AlnScore += ProcessNormalSequencePair(read.seq, AlignmentVec[i].SeedVec[j], cigar_vec);
						}
					}
				}
			}
			if (num > 0)
			{
				if ((j = AlignmentVec[i].SeedVec[0].rPos) > 0) cigar_vec.insert(cigar_vec.begin(), make_pair(j, 'S'));
				if ((j = read.rlen - (AlignmentVec[i].SeedVec[num - 1].rPos + AlignmentVec[i].SeedVec[num - 1].rLen)) > 0) cigar_vec.push_back(make_pair(j, 'S'));
			}
			//if (bDebugMode) printf("Alignment score = %d (rlen=%d) \n", read.AlnReportArr[i].AlnScore, read.rlen);

			read.AlnReportArr[i].coor = GenCoordinateInfo(bFirstRead, AlignmentVec[i].SeedVec[0].gPos, (AlignmentVec[i].SeedVec[num - 1].gPos + AlignmentVec[i].SeedVec[num - 1].gLen - 1), cigar_vec);
			if (read.AlnReportArr[i].AlnScore > read.score)
			{
				read.iBestAlnCanIdx = i;
				read.sub_score = read.score;
				read.score = read.AlnReportArr[i].AlnScore;
			}
			else if (read.AlnReportArr[i].AlnScore == read.score)
			{
				read.sub_score = read.score;
				//if (!bMultiHit && ChromosomeVec[read.AlnReportArr[i].coor.ChromosomeIdx].len > ChromosomeVec[read.AlnReportArr[read.iBestAlnCanIdx].coor.ChromosomeIdx].len) read.iBestAlnCanIdx = i;
			}
		}
	}
	else
	{
		read.CanNum = 1; read.iBestAlnCanIdx = 0;
		read.AlnReportArr = new AlignmentReport_t[1];
		read.AlnReportArr[0].AlnScore = 0;
		read.AlnReportArr[0].PairedAlnCanIdx = -1;
	}
}
