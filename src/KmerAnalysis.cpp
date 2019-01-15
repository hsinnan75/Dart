#include "structure.h"

bool CompByKmerID(const KmerItem_t& p1, const KmerItem_t& p2)
{
	return p1.wid < p2.wid;
}

bool CompByKmerRPos(const KmerPair_t& p1, const KmerPair_t& p2)
{
	return p1.rPos < p2.rPos;
}

bool CompByKmerGPos(const KmerPair_t& p1, const KmerPair_t& p2)
{
	if (p1.gPos == p2.gPos) return p1.rPos < p2.rPos;
	else return p1.gPos < p2.gPos;
}

bool CompByKmerPosDiff(const KmerPair_t& p1, const KmerPair_t& p2)
{
	if (p1.PosDiff == p2.PosDiff) return p1.rPos < p2.rPos;
	else return p1.PosDiff < p2.PosDiff;
}

uint32_t CreateKmerID(const char* seq, short pos)
{
	uint32_t i, id, end_pos = pos + KmerSize;

	for (id = 0, i = pos; i < end_pos; i++) id = (id << 2) + nst_nt4_table[(int)seq[i]];

	return id;
}

vector<KmerItem_t> CreateKmerVecFromReadSeq(int len, char* seq)
{
	int count, head, tail;
	vector<KmerItem_t> vec;
	KmerItem_t WordPosition;

	tail = 0; count = 0; vec.clear();

	while (count < KmerSize && tail < len)
	{
		if (seq[tail++] != 'N') count++;
		else count = 0;
	}
	if (count == KmerSize) // found the first kmer
	{
		WordPosition.pos = (head = tail - KmerSize); WordPosition.wid = CreateKmerID(seq, head);
		vec.push_back(WordPosition);

		for (head += 1; tail < len; head++, tail++)
		{
			if (seq[tail] != 'N')
			{
				WordPosition.pos = head;
				WordPosition.wid = ((WordPosition.wid & KmerPower) << 2) + nst_nt4_table[(int)seq[tail]];
				vec.push_back(WordPosition);
			}
			else
			{
				// find next kmer without 'N'
				count = 0; tail++;
				while (count < KmerSize && tail < len)
				{
					if (seq[tail++] != 'N') count++;
					else count = 0;
				}
				if (count == KmerSize)
				{
					WordPosition.pos = (head = tail - KmerSize); WordPosition.wid = CreateKmerID(seq, head);
					vec.push_back(WordPosition);
				}
				else break;
			}
		}
		sort(vec.begin(), vec.end(), CompByKmerID);
	}
	return vec;
}

vector<KmerPair_t> IdentifyCommonKmers(vector<KmerItem_t>& vec1, vector<KmerItem_t>& vec2)
{
	int i, num;
	uint32_t wid;
	KmerPair_t KmerPair;
	vector<KmerPair_t> KmerPairVec;
	vector<KmerItem_t>::iterator KmerPtr;

	num = (int)vec1.size();
	for (i = 0; i < num; i++)
	{
		wid = vec1[i].wid; KmerPtr = lower_bound(vec2.begin(), vec2.end(), vec1[i], CompByKmerID);
		while (KmerPtr != vec2.end() && KmerPtr->wid == wid)
		{
			KmerPair.rPos = vec1[i].pos;
			KmerPair.PosDiff = (KmerPair.gPos = KmerPtr->pos) - KmerPair.rPos;
			KmerPairVec.push_back(KmerPair);

			KmerPtr++;
		}
	}
	sort(KmerPairVec.begin(), KmerPairVec.end(), CompByKmerPosDiff);

	return KmerPairVec;
}

SeedPair_t GenerateSimplePairsFromCommonKmers(vector<KmerPair_t>& KmerPairVec)
{
	SeedPair_t SeedPair;
	int i, j, l, max_len, PosDiff, num;

	SeedPair.bSimple = true; SeedPair.rLen = SeedPair.gLen = 0;

	num = (int)KmerPairVec.size(); max_len = 0;
	for (i = 0; i < num;)
	{
		//printf("KmerPair %d: %d PosDiff = %d\n", i+1, KmerPairVec[i].rPos, KmerPairVec[i].PosDiff);
		for (PosDiff = KmerPairVec[i].PosDiff, j = i + 1; j < num; j++)
		{
			if (KmerPairVec[j].PosDiff != PosDiff) break;
		}
		if ((l = KmerSize + (KmerPairVec[j-1].rPos - KmerPairVec[i].rPos)) > max_len)
		{
			SeedPair.rPos = KmerPairVec[i].rPos;
			SeedPair.gPos = KmerPairVec[i].gPos;
			max_len = SeedPair.rLen = SeedPair.gLen = l;
		}
		i = j;
	}
	return SeedPair;
}

SeedPair_t GenerateLongestSimplePairsFromFragmentPair(int len1, char* frag1, int len2, char* frag2)
{
	SeedPair_t SeedPair;
	vector<KmerPair_t> KmerPairVec;
	vector<KmerItem_t> KmerVec1, KmerVec2;
	int i, j, l, s, max_len, PosDiff, num;

	KmerVec1 = CreateKmerVecFromReadSeq(len1, frag1); 
	KmerVec2 = CreateKmerVecFromReadSeq(len2, frag2);
	KmerPairVec = IdentifyCommonKmers(KmerVec1, KmerVec2);
	SeedPair.bSimple = true; SeedPair.bAcceptorSite = false; SeedPair.rLen = SeedPair.gLen = 0;

	num = (int)KmerPairVec.size(); max_len = 0;
	for (s = 1, i = 0; i < num;)
	{
		//printf("KmerPair %d: %d PosDiff = %d\n", i+1, KmerPairVec[i].rPos, KmerPairVec[i].PosDiff);
		for (PosDiff = KmerPairVec[i].PosDiff, j = i + 1; j < num; j++)
		{
			if (KmerPairVec[j].PosDiff != PosDiff) break;
			else s++;
		}
		if ((l = KmerSize + (KmerPairVec[j - 1].rPos - KmerPairVec[i].rPos)) > max_len && s > (l - KmerSize)/2)
		{
			SeedPair.rPos = KmerPairVec[i].rPos;
			SeedPair.gPos = KmerPairVec[i].gPos;
			max_len = SeedPair.rLen = SeedPair.gLen = l;
			s = 1;
		}
		i = j;
	}

	return SeedPair;
}

