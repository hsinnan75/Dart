#include <cmath>
#include <sys/time.h>
#include "structure.h"

#define MAPQ_COEF 30

FILE *output = stdout;
time_t StartProcessTime;
bool bSepLibrary = false;
fstream ReadFileHandler1, ReadFileHandler2;
static pthread_mutex_t LibraryLock, OutputLock;
map<pair<int64_t, int64_t>, SpliceJunction_t> SpliceJunctionMap;
int iTotalReadNum = 0, iUniqueMapped = 0, iUnMapped = 0, iPaired = 0;

void ShowAlignmentCandidateInfo(bool bFirstRead, char* header, vector<AlignmentCandidate_t>& AlignmentVec)
{
	vector<AlignmentCandidate_t>::iterator iter;

	printf("\n%s\n", string().assign(100, '-').c_str());

	printf("Alignment Candidate for read %s /%d\n", header, bFirstRead?1:2);
	for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++)
	{
		if (iter->Score == 0) continue;
		printf("\tcandidate#%d: Score=%d\n", (int)(iter - AlignmentVec.begin() + 1), iter->Score);
		ShowSeedInfo(iter->SeedVec);
	}
	printf("%s\n\n", string().assign(100, '-').c_str());

	fflush(stdout);
}

bool CompByCandidateScore(const AlignmentCandidate_t& p1, const AlignmentCandidate_t& p2)
{
	if (p1.Score == p2.Score) return p1.PosDiff < p2.PosDiff;
	else return p1.Score > p2.Score;
}

void SetSingleAlignmentFlag(ReadItem_t& read)
{
	int i;

	if ((read.score > read.sub_score && ++iUniqueMapped) || bMultiHit == false) // unique mapping or bMultiHit=false
	{
		i = read.iBestAlnCanIdx;
		if (read.AlnReportArr[i].coor.bDir == false) read.AlnReportArr[i].iFrag = 0x10;
		else read.AlnReportArr[i].iFrag = 0;
	}
	else if (read.score > 0)
	{
		for (i = 0; i < read.CanNum; i++)
		{
			if (read.AlnReportArr[i].AlnScore > 0)
			{
				if (read.AlnReportArr[i].coor.bDir == false) read.AlnReportArr[i].iFrag = 0x10;
				else read.AlnReportArr[i].iFrag = 0;
			}
		}
	}
	else
	{
		read.AlnReportArr[0].iFrag = 0x4;
	}
}

void SetPairedAlignmentFlag(ReadItem_t& read1, ReadItem_t& read2)
{
	int i, j;

	//printf("read1:[%d, %d]:#%d, read2:[%d, %d] #%d\n", read1.score, read1.sub_score, read1.iBestAlnCanIdx, read2.score, read2.sub_score, read2.iBestAlnCanIdx); fflush(stdout);
	if (read1.score > read1.sub_score && read2.score > read2.sub_score) // unique mapping
	{
		iUniqueMapped += 2;
		i = read1.iBestAlnCanIdx;
		read1.AlnReportArr[i].iFrag = 0x41; // read1 is the first read in a pair

		j = read2.iBestAlnCanIdx;
		read2.AlnReportArr[j].iFrag = 0x81; // read2 is the second read in a pair

		if (j == read1.AlnReportArr[i].PairedAlnCanIdx) // reads are mapped in a proper pair
		{
			read1.AlnReportArr[i].iFrag |= 0x2;
			read2.AlnReportArr[j].iFrag |= 0x2;
		}
		read1.AlnReportArr[i].iFrag |= (read1.AlnReportArr[i].coor.bDir ? 0x20 : 0x10);
		read2.AlnReportArr[j].iFrag |= (read2.AlnReportArr[j].coor.bDir ? 0x20 : 0x10);
	}
	else
	{
		if ((read1.score > read1.sub_score && ++iUniqueMapped) || bMultiHit == false) // unique mapping or bMultiHit=false
		{
			i = read1.iBestAlnCanIdx;
			read1.AlnReportArr[i].iFrag = 0x41; // read1 is the first read in a pair
			read1.AlnReportArr[i].iFrag |= (read1.AlnReportArr[i].coor.bDir ? 0x20 : 0x10);
			if ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0) read1.AlnReportArr[i].iFrag |= 0x2;// reads are mapped in a proper pair
			else read1.AlnReportArr[i].iFrag |= 0x8; // next segment unmapped
		}
		else if (read1.score > 0)
		{
			for (i = 0; i < read1.CanNum; i++)
			{
				if (read1.AlnReportArr[i].AlnScore > 0)
				{
					read1.AlnReportArr[i].iFrag = 0x41; // read1 is the first read in a pair
					read1.AlnReportArr[i].iFrag |= (read1.AlnReportArr[i].coor.bDir ? 0x20 : 0x10);

					if ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0) read1.AlnReportArr[i].iFrag |= 0x2;// reads are mapped in a proper pair
					else read1.AlnReportArr[i].iFrag |= 0x8; // next segment unmapped
				}
			}
		}
		else
		{
			read1.AlnReportArr[0].iFrag = 0x41; // read1 is the first read in a pair
			read1.AlnReportArr[0].iFrag |= 0x4;

			if (read2.score == 0) read1.AlnReportArr[0].iFrag |= 0x8; // next segment unmapped
			else read1.AlnReportArr[0].iFrag |= (read2.AlnReportArr[read2.iBestAlnCanIdx].coor.bDir ? 0x10 : 0x20);
		}

		if ((read2.score > read2.sub_score && ++iUniqueMapped) || bMultiHit == false) // unique mapping or bMultiHit=false
		{
			j = read2.iBestAlnCanIdx;
			read2.AlnReportArr[j].iFrag = 0x81; // read2 is the second read in a pair
			read2.AlnReportArr[j].iFrag |= (read2.AlnReportArr[j].coor.bDir ? 0x20 : 0x10);

			if ((i = read2.AlnReportArr[j].PairedAlnCanIdx) != -1 && read1.AlnReportArr[i].AlnScore > 0) read2.AlnReportArr[j].iFrag |= 0x2;// reads are mapped in a proper pair
			else read2.AlnReportArr[j].iFrag |= 0x8; // next segment unmapped
		}
		else if (read2.score > 0)
		{
			for (j = 0; j < read2.CanNum; j++)
			{
				if (read2.AlnReportArr[j].AlnScore > 0)
				{
					read2.AlnReportArr[j].iFrag = 0x81; // read2 is the second read in a pair
					read2.AlnReportArr[j].iFrag |= (read2.AlnReportArr[j].coor.bDir ? 0x20 : 0x10);

					if ((i = read2.AlnReportArr[j].PairedAlnCanIdx) != -1 && read1.AlnReportArr[i].AlnScore > 0) read2.AlnReportArr[j].iFrag |= 0x2;// reads are mapped in a proper pair
					else read2.AlnReportArr[j].iFrag |= 0x8; // next segment unmapped
				}
			}
		}
		else
		{
			read2.AlnReportArr[0].iFrag = 0x81; // read2 is the second read in a pair
			read2.AlnReportArr[0].iFrag |= 0x4; // segment unmapped
			if (read1.score == 0) read2.AlnReportArr[0].iFrag |= 0x8; // next segment unmapped
			else read2.AlnReportArr[0].iFrag |= (read1.AlnReportArr[read1.iBestAlnCanIdx].coor.bDir ? 0x10 : 0x20);
		}
	}
}

void EvaluateMAPQ(ReadItem_t& read)
{
	int i, iMap;

	if (read.score == 0 || read.score == read.sub_score) read.mapq = 0;
	else
	{
		if (read.sub_score == 0 || read.score - read.sub_score > 10) read.mapq = 255;
		else
		{
			for (iMap = 0, i = 0; i < read.CanNum; i++)
				if (read.AlnReportArr[i].AlnScore == read.score) iMap++;
			if(iMap > 1) read.mapq = (int)(-10 * log10(1 - (1.0 / iMap)));
			else read.mapq = 0;
		}
	}
}

void OutputPairedSamFile(ReadItem_t& read1, ReadItem_t& read2)
{
	char *seq, *rseq;
	int i, j, dist = 0;

	if (read1.score == 0 || (bUnique && read1.mapq == 0))
	{
		iUnMapped++;
		fprintf(output, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t*\tAS:i:0\tXS:i:0\n", read1.header, read1.AlnReportArr[0].iFrag, read1.seq);
	}
	else
	{
		seq = read1.seq; rseq = NULL;
		for (i = read1.iBestAlnCanIdx; i < read1.CanNum; i++)
		{
			if (read1.AlnReportArr[i].AlnScore > 0)
			{
				if (read1.AlnReportArr[i].coor.bDir == false && rseq == NULL)
				{
					rseq = new char[read1.rlen + 1]; rseq[read1.rlen] = '\0';
					GetComplementarySeq(read1.rlen, seq, rseq);
				}
				if ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0)
				{
					dist = (int)(read2.AlnReportArr[j].coor.gPos - read1.AlnReportArr[i].coor.gPos + (read1.AlnReportArr[i].coor.bDir ? read2.rlen : 0 - read1.rlen));
					if (i == read1.iBestAlnCanIdx) iPaired += 2;
					fprintf(output, "%s\t%d\t%s\t%ld\t%d\t%s\t=\t%ld\t%d\t%s\t*\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read1.header, read1.AlnReportArr[i].iFrag, ChromosomeVec[read1.AlnReportArr[i].coor.ChromosomeIdx].name, read1.AlnReportArr[i].coor.gPos, read1.mapq, read1.AlnReportArr[i].coor.CIGAR.c_str(), read2.AlnReportArr[j].coor.gPos, dist, (read1.AlnReportArr[i].coor.bDir ? seq : rseq), read1.rlen - read1.score, read1.score, read1.sub_score);
				}
				else fprintf(output, "%s\t%d\t%s\t%ld\t%d\t%s\t*\t0\t0\t%s\t*\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read1.header, read1.AlnReportArr[i].iFrag, ChromosomeVec[read1.AlnReportArr[i].coor.ChromosomeIdx].name, read1.AlnReportArr[i].coor.gPos, read1.mapq, read1.AlnReportArr[i].coor.CIGAR.c_str(), (read1.AlnReportArr[i].coor.bDir ? seq : rseq), read1.rlen - read1.score, read1.score, read1.sub_score);
			}
			if (!bMultiHit) break;
		}
		if (rseq != NULL)
		{
			delete[] rseq;
			rseq = NULL;
		}
	}

	if (read2.score == 0 || (bUnique && read2.mapq == 0))
	{
		iUnMapped++;
		fprintf(output, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t*\tAS:i:0\tXS:i:0\n", read2.header, read2.AlnReportArr[0].iFrag, read2.seq);
	}
	else
	{
		rseq = read2.seq; seq = NULL;
		for (j = read2.iBestAlnCanIdx; j < read2.CanNum; j++)
		{
			if (read2.AlnReportArr[j].AlnScore > 0)
			{
				if (read2.AlnReportArr[j].coor.bDir == true && seq == NULL)
				{
					seq = new char[read2.rlen + 1]; seq[read2.rlen] = '\0';
					GetComplementarySeq(read2.rlen, rseq, seq);
				}
				if ((i = read2.AlnReportArr[j].PairedAlnCanIdx) != -1 && read1.AlnReportArr[i].AlnScore > 0)
				{
					dist = 0 - ((int)(read2.AlnReportArr[j].coor.gPos - read1.AlnReportArr[i].coor.gPos + (read1.AlnReportArr[i].coor.bDir ? read2.rlen : 0 - read1.rlen)));
					fprintf(output, "%s\t%d\t%s\t%ld\t%d\t%s\t=\t%ld\t%d\t%s\t*\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read2.header, read2.AlnReportArr[j].iFrag, ChromosomeVec[read2.AlnReportArr[j].coor.ChromosomeIdx].name, read2.AlnReportArr[j].coor.gPos, read2.mapq, read2.AlnReportArr[j].coor.CIGAR.c_str(), read1.AlnReportArr[i].coor.gPos, dist, (read2.AlnReportArr[j].coor.bDir ? seq : rseq), read2.rlen - read2.score, read2.score, read2.sub_score);
				}
				else fprintf(output, "%s\t%d\t%s\t%ld\t%d\t%s\t*\t0\t0\t%s\t*\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read2.header, read2.AlnReportArr[j].iFrag, ChromosomeVec[read2.AlnReportArr[j].coor.ChromosomeIdx].name, read2.AlnReportArr[j].coor.gPos, read2.mapq, read2.AlnReportArr[j].coor.CIGAR.c_str(), (read2.AlnReportArr[j].coor.bDir ? seq : rseq), read2.rlen - read2.score, read2.score, read2.sub_score);
			}
			if (!bMultiHit) break;
		}
		if (seq != NULL)
		{
			delete[] seq;
			seq = NULL;
		}
	}
}

void OutputSingledSamFile(ReadItem_t& read)
{
	if (read.score == 0 || (bUnique && read.mapq == 0))
	{
		iUnMapped++;
		fprintf(output, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t*\tAS:i:0\tXS:i:0\n", read.header, read.AlnReportArr[0].iFrag, read.seq);
	}
	else
	{
		int i;
		char *seq, *rseq;

		seq = read.seq; rseq = NULL;
		for (i = read.iBestAlnCanIdx; i < read.CanNum; i++)
		{
			if (read.AlnReportArr[i].AlnScore == read.score)
			{
				if (read.AlnReportArr[i].coor.bDir == false && rseq == NULL)
				{
					rseq = new char[read.rlen + 1]; rseq[read.rlen] = '\0';
					GetComplementarySeq(read.rlen, seq, rseq);
				}
				fprintf(output, "%s\t%d\t%s\t%ld\t%d\t%s\t*\t0\t0\t%s\t*\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read.header, read.AlnReportArr[i].iFrag, ChromosomeVec[read.AlnReportArr[i].coor.ChromosomeIdx].name, read.AlnReportArr[i].coor.gPos, read.mapq, read.AlnReportArr[i].coor.CIGAR.c_str(), (read.AlnReportArr[i].coor.bDir ? seq : rseq), read.rlen - read.score, read.score, read.sub_score);
				if (!bMultiHit) break;
			}
		}
		if (rseq != NULL)
		{
			delete[] rseq;
			rseq = NULL;
		}
	}
}

void RemoveRedundantCandidates(vector<AlignmentCandidate_t>& AlignmentVec)
{
	bool bAmbiguous = false;
	int thr, score1, score2;
	vector<AlignmentCandidate_t>::iterator iter;

	if (AlignmentVec.size() <= 1) return;
	else
	{
		score1 = score2 = 0;
		for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++)
		{
			if (iter->Score > score2)
			{
				if (iter->Score >= score1)
				{
					score2 = score1;
					score1 = iter->Score;
				}
				else score2 = iter->Score;
			}
			else if (iter->Score == score2) score2 = score1;
		}
		//if (score1 == score2) bAmbiguous = true;

		if (score1 == score2 || score1 - score2 > 20) thr = score1;
		else thr = score2;

		if (bDebugMode) printf("Candidate score threshold = %d\n", thr);
		for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++) if (iter->Score < thr) iter->Score = 0;
	}
	//return bAmbiguous;
}

bool CheckPairedAlignmentCandidates(vector<AlignmentCandidate_t>& AlignmentVec1, vector<AlignmentCandidate_t>& AlignmentVec2)
{
	bool bPairing = false;
	int64_t dist, min_dist;
	int i, j, best_mate, num1, num2;

	num1 = (int)AlignmentVec1.size(); num2 = (int)AlignmentVec2.size();

	if (num1*num2 > 1000)
	{
		RemoveRedundantCandidates(AlignmentVec1);
		RemoveRedundantCandidates(AlignmentVec2);
	}
	for (i = 0; i != num1; i++)
	{
		if (AlignmentVec1[i].Score == 0) continue;

		for (best_mate = -1, min_dist = MaxIntronSize, j = 0; j != num2; j++)
		{
			if (AlignmentVec2[j].Score == 0) continue;

			dist = abs(AlignmentVec2[j].PosDiff - AlignmentVec1[i].PosDiff);
			//printf("#%d (s=%d) and #%d (s=%d) (dist=%ld)\n", i+1, AlignmentVec1[i].Score, j+1, AlignmentVec2[j].Score, dist), fflush(stdout);
			if (dist < min_dist)
			{
				best_mate = j;
				min_dist = dist;
			}
		}
		if (best_mate != -1)
		{
			j = best_mate;
			if (AlignmentVec2[j].PairedAlnCanIdx == -1)
			{
				bPairing = true;
				AlignmentVec1[i].PairedAlnCanIdx = j;
				AlignmentVec2[j].PairedAlnCanIdx = i;
			}
			else if (AlignmentVec1[i].Score > AlignmentVec1[AlignmentVec2[j].PairedAlnCanIdx].Score)
			{
				AlignmentVec1[AlignmentVec2[j].PairedAlnCanIdx].PairedAlnCanIdx = -1;
				AlignmentVec1[i].PairedAlnCanIdx = j;
				AlignmentVec2[j].PairedAlnCanIdx = i;
			}
		}
	}
	return bPairing;
}

void RemoveUnMatedAlignmentCandidates(vector<AlignmentCandidate_t>& AlignmentVec1, vector<AlignmentCandidate_t>& AlignmentVec2)
{
	int i, j, num1, num2;

	num1 = (int)AlignmentVec1.size(); num2 = (int)AlignmentVec2.size();

	for (i = 0; i != num1; i++)
	{
		if (AlignmentVec1[i].PairedAlnCanIdx == -1) AlignmentVec1[i].Score = 0;
		else
		{
			j = AlignmentVec1[i].PairedAlnCanIdx;
			AlignmentVec1[i].Score = AlignmentVec2[j].Score = AlignmentVec1[i].Score + AlignmentVec2[j].Score;
		}
	}
	for (j = 0; j != num2; j++) if (AlignmentVec2[j].PairedAlnCanIdx == -1) AlignmentVec2[j].Score = 0;

	if (bDebugMode)
	{
		for (i = 0; i != num1; i++)
		{
			if ((j = AlignmentVec1[i].PairedAlnCanIdx) != -1)
				printf("#%d(s=%d) and #%d(s=%d) are pairing\n", i + 1, AlignmentVec1[i].Score, j + 1, AlignmentVec2[j].Score);
		}
	}
}

void CheckPairedFinalAlignments(ReadItem_t& read1, ReadItem_t& read2)
{
	bool bMated;
	int64_t dist;
	int i, j, best_mate, s;

	//printf("BestIdx1=%d, BestIdx2=%d\n", read1.iBestAlnCanIdx + 1, read2.iBestAlnCanIdx + 1);
	bMated = read1.AlnReportArr[read1.iBestAlnCanIdx].PairedAlnCanIdx == read2.iBestAlnCanIdx ? true : false;
	if (!bMultiHit && bMated) return;
	if (!bMated && read1.score > 0 && read2.score > 0) // identify mated pairs
	{
		for (s = 0, i = 0; i != read1.CanNum; i++)
		{
			if (read1.AlnReportArr[i].AlnScore > 0 && (j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0)
			{
				bMated = true;
				if (s < read1.AlnReportArr[i].AlnScore + read2.AlnReportArr[j].AlnScore)
				{
					s = read1.AlnReportArr[i].AlnScore + read2.AlnReportArr[j].AlnScore;
					read1.iBestAlnCanIdx = i; read1.score = read1.AlnReportArr[i].AlnScore;
					read2.iBestAlnCanIdx = j; read2.score = read2.AlnReportArr[j].AlnScore;
				}
			}
		}
	}
	if (bMated)
	{
		for (i = 0; i != read1.CanNum; i++)
		{
			if (read1.AlnReportArr[i].AlnScore != read1.score || ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore != read2.score))
			{
				read1.AlnReportArr[i].AlnScore = 0;
				read1.AlnReportArr[i].PairedAlnCanIdx = -1;
				continue;
			}
		}
	}
	else // remove all mated info
	{
		for (i = 0; i != read1.CanNum; i++)
		{
			if (read1.AlnReportArr[i].PairedAlnCanIdx != -1) read1.AlnReportArr[i].PairedAlnCanIdx = -1;
			if (read1.AlnReportArr[i].AlnScore > 0 && read1.AlnReportArr[i].AlnScore != read1.score) read1.AlnReportArr[i].AlnScore = 0;
		}
		for (j = 0; j != read2.CanNum; j++)
		{
			if (read2.AlnReportArr[j].PairedAlnCanIdx != -1) read2.AlnReportArr[j].PairedAlnCanIdx = -1;
			if (read2.AlnReportArr[j].AlnScore > 0 && read2.AlnReportArr[j].AlnScore != read2.score) read2.AlnReportArr[j].AlnScore = 0;
		}
	}
}

void UpdateLocalSJMap(AlignmentCandidate_t& Aln, map<pair<int64_t, int64_t>, SpliceJunction_t>& LocalSJMap)
{
	if (Aln.SJtype == -1) return;

	int64_t g1, g2;
	SpliceJunction_t SpliceJunction;
	int i, num = (int)Aln.SeedVec.size();
	map<pair<int64_t, int64_t>, SpliceJunction_t>::iterator iter;

	for (i = 1; i < num; i++)
	{
		if (Aln.SeedVec[i].bAcceptorSite)
		{
			if (Aln.PosDiff < GenomeSize)
			{
				g1 = Aln.SeedVec[i - 1].gPos + Aln.SeedVec[i - 1].gLen;
				g2 = Aln.SeedVec[i].gPos - 1;
			}
			else
			{
				g1 = TwoGenomeSize - Aln.SeedVec[i].gPos;
				g2 = TwoGenomeSize - 1 - (Aln.SeedVec[i - 1].gPos + Aln.SeedVec[i - 1].gLen);
			}
			if ((iter = LocalSJMap.find(make_pair(g1, g2))) != LocalSJMap.end()) iter->second.iCount++;
			else
			{
				SpliceJunction.iCount = 1;
				SpliceJunction.type = Aln.SJtype;
				LocalSJMap.insert(make_pair(make_pair(g1, g2), SpliceJunction));
			}
		}
	}
}

void UpdateGlobalSJMap(map<pair<int64_t, int64_t>, SpliceJunction_t>& LocalSJMap)
{
	map<pair<int64_t, int64_t>, SpliceJunction_t>::iterator iter, ft;

	for (iter = LocalSJMap.begin(); iter != LocalSJMap.end(); iter++)
	{
		ft = SpliceJunctionMap.find(iter->first);
		if (ft != SpliceJunctionMap.end()) ft->second.iCount += iter->second.iCount;
		else SpliceJunctionMap.insert(make_pair(iter->first, iter->second));
	}
}

void *ReadMapping(void *arg)
{
	bool bAmb1, bAmb2;
	ReadItem_t* ReadArr = NULL;
	int i, j, max_dist, ReadNum, EstDistance;
	vector<SeedPair_t> SeedPairVec1, SeedPairVec2;
	map<pair<int64_t, int64_t>, SpliceJunction_t> LocalSJMap;
	vector<AlignmentCandidate_t> AlignmentVec1, AlignmentVec2;

	ReadArr = new ReadItem_t[ReadChunkSize];
	while (true)
	{
		pthread_mutex_lock(&LibraryLock);
		ReadNum = GetNextChunk(bSepLibrary, ReadFileHandler1, ReadFileHandler2, ReadArr);
		fprintf(stderr, "\r%d %s tags have been processed in %ld seconds...", iTotalReadNum, (bPairEnd? "paired-end":"singled-end"), (time(NULL) - StartProcessTime));
		pthread_mutex_unlock(&LibraryLock);
		
		if (ReadNum == 0) break;
		if (bPairEnd && ReadNum % 2 == 0)
		{
			for (i = 0, j = 1; i != ReadNum; i += 2, j += 2)
			{
				//if (bDebugMode) printf("Mapping paired reads#%d %s (len=%d) and %s (len=%d):\n", i + 1, ReadArr[i].header + 1, ReadArr[i].rlen, ReadArr[j].header + 1, ReadArr[j].rlen);

				SeedPairVec1 = IdentifySeedPairs(ReadArr[i].rlen, ReadArr[i].EncodeSeq); //if (bDebugMode) ShowSeedInfo(SeedPairVec1);
				AlignmentVec1 = GenerateAlignmentCandidate(ReadArr[i].rlen, SeedPairVec1);

				SeedPairVec2 = IdentifySeedPairs(ReadArr[j].rlen, ReadArr[j].EncodeSeq); //if (bDebugMode) ShowSeedInfo(SeedPairVec2);
				AlignmentVec2 = GenerateAlignmentCandidate(ReadArr[j].rlen, SeedPairVec2);

				//if (bDebugMode) ShowAlignmentCandidateInfo(true, ReadArr[i].header+1, AlignmentVec1), ShowAlignmentCandidateInfo(false, ReadArr[j].header+1, AlignmentVec2);
				if (CheckPairedAlignmentCandidates(AlignmentVec1, AlignmentVec2))
				{
					RemoveUnMatedAlignmentCandidates(AlignmentVec1, AlignmentVec2);
				}
				//else 
				//{
				//	AlignmentVec1.clear();
				//	AlignmentVec2.clear();
				//}
				RemoveRedundantCandidates(AlignmentVec1); RemoveRedundantCandidates(AlignmentVec2);
				
				//if (bAmb1 || bAmb2)
				//{
				//	AlignmentVec1.clear();
				//	AlignmentVec2.clear();
				//}
				//if (bDebugMode)
				//{
				//	ShowAlignmentCandidateInfo(1, ReadArr[i].header + 1, AlignmentVec1);
				//	ShowAlignmentCandidateInfo(0, ReadArr[j].header + 1, AlignmentVec2);
				//}
				GenMappingReport(true,  ReadArr[i], AlignmentVec1);
				GenMappingReport(false, ReadArr[j], AlignmentVec2);

				CheckPairedFinalAlignments(ReadArr[i], ReadArr[j]);

				SetPairedAlignmentFlag(ReadArr[i], ReadArr[j]);
				EvaluateMAPQ(ReadArr[i]); EvaluateMAPQ(ReadArr[j]);

				if (SJFileName != NULL && ReadArr[i].mapq == 255) UpdateLocalSJMap(AlignmentVec1[ReadArr[i].iBestAlnCanIdx], LocalSJMap);
				if (SJFileName != NULL && ReadArr[j].mapq == 255) UpdateLocalSJMap(AlignmentVec2[ReadArr[j].iBestAlnCanIdx], LocalSJMap);

				//if (bDebugMode) printf("\nEnd of mapping for read#%s\n%s\n", ReadArr[i].header, string().assign(100, '=').c_str());
			}
		}
		else
		{
			for (i = 0; i != ReadNum; i++)
			{
				//if (bDebugMode) printf("Mapping single read#%d %s (len=%d):\n", i + 1, ReadArr[i].header + 1, ReadArr[i].rlen);
				//fprintf(stdout, "%s\n", ReadArr[i].header + 1); fflush(output);

				SeedPairVec1 = IdentifySeedPairs(ReadArr[i].rlen, ReadArr[i].EncodeSeq); //if (bDebugMode) ShowSeedInfo(SeedPairVec1);
				AlignmentVec1 = GenerateAlignmentCandidate(ReadArr[i].rlen, SeedPairVec1);
				RemoveRedundantCandidates(AlignmentVec1); if (bDebugMode) ShowAlignmentCandidateInfo(1, ReadArr[i].header+1, AlignmentVec1);
				GenMappingReport(true, ReadArr[i], AlignmentVec1);
				SetSingleAlignmentFlag(ReadArr[i]); EvaluateMAPQ(ReadArr[i]);

				if (SJFileName != NULL && ReadArr[i].mapq == 255) UpdateLocalSJMap(AlignmentVec1[ReadArr[i].iBestAlnCanIdx], LocalSJMap);

				//if (bDebugMode) printf("\nEnd of mapping for read#%s\n%s\n", ReadArr[i].header, string().assign(100, '=').c_str());
			}
		}
		pthread_mutex_lock(&OutputLock);
		iTotalReadNum += ReadNum;
		if (bPairEnd && ReadNum %2 == 0) for (i = 0, j = 1; i != ReadNum; i += 2, j += 2) OutputPairedSamFile(ReadArr[i], ReadArr[j]);
		else for (i = 0; i != ReadNum; i++) OutputSingledSamFile(ReadArr[i]);
		pthread_mutex_unlock(&OutputLock);

		for (i = 0; i != ReadNum; i++)
		{
			delete[] ReadArr[i].header;
			delete[] ReadArr[i].seq;
			delete[] ReadArr[i].EncodeSeq;
			delete[] ReadArr[i].AlnReportArr;
		}
	}
	delete[] ReadArr;

	pthread_mutex_lock(&OutputLock);
	UpdateGlobalSJMap(LocalSJMap);
	pthread_mutex_unlock(&OutputLock);

	return (void*)(1);
}

int AbsLoc2ChrLoc(int64_t& g1, int64_t& g2)
{
	int ChrIdx;
	map<int64_t, int>::iterator iter;

	iter = ChrLocMap.lower_bound(g1); ChrIdx = iter->second;
	g1 = g1 + 1 - ChromosomeVec[ChrIdx].FowardLocation;
	g2 = g2 + 1 - ChromosomeVec[ChrIdx].FowardLocation;
	return ChrIdx;
}

void OutputSpliceJunctions()
{
	bool bDir;
	int64_t g1, g2;
	int i, n, ChrIdx;
	int SJtypeNum[4] = { 0, 0, 0, 0 };
	FILE *SJFile = fopen(SJFileName, "w");

	for (map<pair<int64_t, int64_t>, SpliceJunction_t>::iterator iter = SpliceJunctionMap.begin(); iter != SpliceJunctionMap.end(); iter++)
	{
		g1 = iter->first.first; g2 = iter->first.second;
		ChrIdx = AbsLoc2ChrLoc(g1, g2); //SJtypeNum[iter->second.type]++;
		fprintf(SJFile, "%s\t%d\t%d\t%d\n", ChromosomeVec[ChrIdx].name, g1, g2, iter->second.iCount);
	}
	if ((n = (int)SpliceJunctionMap.size()) > 0)
	{
		fprintf(stderr, "\nFound %d splice junctions (file: %s)\n", n, SJFileName);
		//for (i = 0; i < 4; i++) fprintf(stderr, "\t%s (%.2f%%)\n", SpliceJunctionArr[i], (int)(10000 * (1.0*SJtypeNum[i] / n)) / 100.0);
	}
	SpliceJunctionMap.clear(); fclose(SJFile);
}

void Mapping()
{
	int i;

	ReadFileHandler1.open(ReadFileName, ios_base::in);
	if (ReadFileName2 != NULL)
	{
		bSepLibrary = bPairEnd = true;
		ReadFileHandler2.open(ReadFileName2, ios_base::in);
	}
	if (SamFileName != NULL) output = fopen(SamFileName, "w");

	if (bDebugMode) iThreadNum = 1;
	else
	{
		fprintf(output, "@PG\tPN:Dart\tVN:%s\n", VersionStr);
		for (i = 0; i < iChromsomeNum; i++) fprintf(output, "@SQ\tSN:%s\tLN:%ld\n", ChromosomeVec[i].name, ChromosomeVec[i].len);
	}

	//iThreadNum = 1;
	StartProcessTime = time(NULL);
	pthread_t *ThreadArr = new pthread_t[iThreadNum];
	for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, ReadMapping, NULL);
	for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);

	ReadFileHandler1.close(); if (ReadFileName2 != NULL) ReadFileHandler2.close();
	fprintf(stderr, "\rAll the %d %s reads have been processed in %lld seconds.\n", iTotalReadNum, (bPairEnd? "paired-end":"single-end"), (long long)(time(NULL) - StartProcessTime));

	if (SamFileName != NULL) fclose(output);

	delete[] ThreadArr;

	if(iTotalReadNum > 0)
	{
		if (bPairEnd) fprintf(stderr, "\t# of total mapped reads = %d (sensitivity = %.2f%%)\n\t# of paired sequences = %d (%.2f%%)\n", iTotalReadNum - iUnMapped, (int)(10000 * (1.0*(iTotalReadNum - iUnMapped) / iTotalReadNum) + 0.5) / 100.0, iPaired, (int)(10000 * (1.0*iPaired / iTotalReadNum) + 0.5) / 100.0);
		else fprintf(stderr, "\t# of total mapped reads = %d (sensitivity = %.2f%%)\n", iTotalReadNum - iUnMapped, (int)(10000 * (1.0*(iTotalReadNum - iUnMapped) / iTotalReadNum) + 0.5) / 100.0);
		fprintf(stderr, "\t# of unique mapped reads = %d (%.2f%%)\n", iUniqueMapped, (int)(10000 * (1.0*iUniqueMapped / iTotalReadNum) + 0.5) / 100.0);
		if (!bUnique) fprintf(stderr, "\t# of multiple mapped reads = %d (%.2f%%)\n", (iTotalReadNum - iUnMapped - iUniqueMapped), (int)(10000 * (1.0*(iTotalReadNum - iUnMapped - iUniqueMapped) / iTotalReadNum) + 0.5) / 100.0);
		fprintf(stderr, "\t# of unmapped reads = %d (%.2f%%)\n", iUnMapped, (int)(10000 * (1.0*iUnMapped / iTotalReadNum) + 0.5) / 100.0);
	}
	if (SJFileName != NULL) OutputSpliceJunctions();
}
