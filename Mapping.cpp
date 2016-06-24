#include <cmath>
#include "structure.h"

#define MAPQ_COEF 30

FILE *output = stdout;
time_t StartProcessTime;
bool bSepLibrary = false;
int iTotalReadNum = 0, iUnMapped = 0;
fstream ReadFileHandler1, ReadFileHandler2;
static pthread_mutex_t LibraryLock, OutputLock;

void ShowAlignmentCandidateInfo(char* header, vector<AlignmentCandidate_t>& AlignmentVec)
{
	vector<AlignmentCandidate_t>::iterator iter;

	printf("\n%s\n", string().assign(100, '-').c_str());

	printf("Alignment Candidate for read %s\n", header);
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

void SetPairedAlignmentFlag(MappingReport_t& report1, MappingReport_t& report2)
{
	report1.iFrag = 0x41; // read1 is the first read in a pair
	report2.iFrag = 0x81; // read2 is the second read in a pair

	if (report1.bPaired) // the reads are mapped in a proper pair
	{
		report1.iFrag |= 0x2;
		report2.iFrag |= 0x2;

		report1.iFrag |= (report1.Coor.bDir ? 0x20 : 0x10);
		report2.iFrag |= (report2.Coor.bDir ? 0x20 : 0x10);
	}
	else
	{
		if (report1.score == 0) // read1 is unmapped
		{
			report1.iFrag |= 0x4; // segment unmapped
			report2.iFrag |= 0x8; // next segment unmapped

			if (report2.score > 0)
			{
				report1.iFrag |= (report2.Coor.bDir ? 0x20 : 0x10);
				report2.iFrag |= (report2.Coor.bDir ? 0x10 : 0x20);
			}
		}
		else
		{
			report1.iFrag |= (report1.Coor.bDir ? 0x20 : 0x10);
			report2.iFrag |= (report1.Coor.bDir ? 0x10 : 0x20);
		}

		if (report2.score == 0) // read2 is unmapped
		{
			report1.iFrag |= 0x8; // next segment unmapped
			report2.iFrag |= 0x4; // segment unmapped

			if (report1.score > 0)
			{
				report1.iFrag |= (report1.Coor.bDir ? 0x10 : 0x20);
				report2.iFrag |= (report1.Coor.bDir ? 0x20 : 0x10);
			}
		}
		else
		{
			report1.iFrag |= (report2.Coor.bDir ? 0x10 : 0x20);
			report2.iFrag |= (report2.Coor.bDir ? 0x20 : 0x10);
		}
	}
}

void EvaluateMAPQ(int rlen, MappingReport_t& report)
{
	float f;

	if (report.score == 0 || report.score == report.sub_score) report.mapq = 0;
	else
	{
		if (report.sub_score == 0) report.mapq = 255;
		else
		{
			report.mapq = (int)(MAPQ_COEF * (1 - (float)(report.score - report.sub_score)/10.0)*log(report.score) + 0.4999);
			if (report.mapq > 255) report.mapq = 255;
			else if(report.mapq < 0) report.mapq = 0;

			f = 1.0*report.score / rlen;
			report.mapq = (f < 0.95 ? (int)(report.mapq * f * f) : report.mapq);
		}
	}
}

void OutputPairedSamFile(ReadItem_t& read1, MappingReport_t& report1, ReadItem_t& read2, MappingReport_t& report2)
{
	int i, j, dist, paired_idx;

	if (report1.score == 0)
	{
		iUnMapped++;
		fprintf(output, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t*\tAS:i:0\tXS:i:0\n", read1.header + 1, report1.iFrag, read1.seq);
		//if(report2.score == 0) fprintf(output, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t*\tAS:i:0\tXS:i:0\n", read1.header + 1, report1.iFrag, read1.seq);
		//else fprintf(output, "%s\t%d\t%s\t%ld\t0\t*\t=\t%ld\t0\t%s\t*\tAS:i:0\tXS:i:0\n", read1.header + 1, report1.iFrag, ChromosomeVec[report2.Coor.ChromosomeIdx].name, report2.Coor.gPos, report2.Coor.gPos, read1.seq);
	}
	else
	{
		if (report1.bPaired)
		{
			dist = (int)(report2.Coor.gPos - report1.Coor.gPos + (report1.Coor.bDir ? read2.rlen : 0 - read1.rlen));
			fprintf(output, "%s\t%d\t%s\t%ld\t%d\t%s\t=\t%ld\t%d\t%s\t*\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read1.header + 1, report1.iFrag, ChromosomeVec[report1.Coor.ChromosomeIdx].name, report1.Coor.gPos, report1.mapq, report1.Coor.CIGAR.c_str(), report2.Coor.gPos, dist, read1.seq, read1.rlen - report1.score, report1.score, report1.sub_score);
		}
		else
		{
			fprintf(output, "%s\t%d\t%s\t%ld\t%d\t%s\t*\t0\t0\t%s\t*\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read1.header + 1, report1.iFrag, ChromosomeVec[report1.Coor.ChromosomeIdx].name, report1.Coor.gPos, report1.mapq, report1.Coor.CIGAR.c_str(), read1.seq, read1.rlen - report1.score, report1.score, report1.sub_score);
		}
	}
	if (report2.score == 0)
	{
		iUnMapped++;
		if (report1.score == 0) fprintf(output, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t*\tAS:i:0\tXS:i:0\n", read2.header + 1, report2.iFrag, read2.seq);
		//if (report1.score == 0) fprintf(output, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t*\tAS:i:0\tXS:i:0\n", read2.header + 1, report2.iFrag, read2.seq);
		//else fprintf(output, "%s\t%d\t%s\t%ld\t0\t*\t=\t%ld\t0\t%s\t*\tAS:i:0\tXS:i:0\n", read2.header + 1, report2.iFrag, ChromosomeVec[report1.Coor.ChromosomeIdx].name, report1.Coor.gPos, report1.Coor.gPos, read2.seq);
	}
	else
	{
		if (report2.bPaired)
		{
			dist = 0 - dist;
			fprintf(output, "%s\t%d\t%s\t%ld\t%d\t%s\t=\t%ld\t%d\t%s\t*\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read2.header + 1, report2.iFrag, ChromosomeVec[report2.Coor.ChromosomeIdx].name, report2.Coor.gPos, report2.mapq, report2.Coor.CIGAR.c_str(), report1.Coor.gPos, dist, read2.seq, read2.rlen - report2.score, report2.score, report2.sub_score);
		}
		else
		{
			fprintf(output, "%s\t%d\t%s\t%ld\t%d\t%s\t*\t0\t0\t%s\t*\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read2.header + 1, report2.iFrag, ChromosomeVec[report2.Coor.ChromosomeIdx].name, report2.Coor.gPos, report2.mapq, report2.Coor.CIGAR.c_str(), read2.seq, read2.rlen - report2.score, report2.score, report2.sub_score);
		}
	}
}

void OutputSingledSamFile(ReadItem_t& read, MappingReport_t& report)
{
	if (report.score == 0)
	{
		iUnMapped++;
		fprintf(output, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t*\tAS:i:0\tXS:i:0\n", read.header + 1, report.iFrag, read.seq);
	}
	else
	{
		fprintf(output, "%s\t%d\t%s\t%ld\t%d\t%s\t*\t0\t0\t%s\t*\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read.header + 1, report.iFrag, ChromosomeVec[report.Coor.ChromosomeIdx].name, report.Coor.gPos, report.mapq, report.Coor.CIGAR.c_str(), read.seq, read.rlen - report.score, report.score, report.sub_score);
	}
}

void RemoveRedundantCandidates(vector<AlignmentCandidate_t>& AlignmentVec)
{
	int thr;
	vector<AlignmentCandidate_t>::iterator iter;

	if (AlignmentVec.size() <= 1) return;
	else
	{
		int score1, score2;

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
		if (score1 == score2 || score1 - score2 > 20) thr = score1;
		else thr = score2;

		if (bDebugMode) printf("Candidate score threshold = %d\n", thr);
		for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++) if (iter->Score < thr) iter->Score = 0;
	}
}

bool CheckPairedAlignmentCandidates(vector<AlignmentCandidate_t>& AlignmentVec1, vector<AlignmentCandidate_t>& AlignmentVec2)
{
	bool bPaired = false;
	int64_t dist, min_dist;
	int i, j, p, s, max_score, num1, num2;

	num1 = (int)AlignmentVec1.size(); num2 = (int)AlignmentVec2.size(); max_score = 0;
	for (i = 0; i != num1; i++)
	{
		if (AlignmentVec1[i].Score == 0) continue;

		min_dist = 100000;
		for (p = -1, j = 0; j != num2; j++)
		{
			if (AlignmentVec2[j].Score == 0) continue;

			dist = abs(AlignmentVec2[j].PosDiff - AlignmentVec1[i].PosDiff);
			//printf("#%d (s=%d) and #%d (s=%d), dist=%ld\n", i + 1, AlignmentVec1[i].Score, j + 1, AlignmentVec2[j].Score, dist);
			if (dist < min_dist)
			{
				min_dist = dist;
				p = j;
			}
		}
		if (p != -1)
		{
			j = p;
			s = AlignmentVec1[i].Score + AlignmentVec2[j].Score;
			if(bDebugMode) printf("Can#%d (s=%d) and Can#%d (s=%d) are pairing\n", i + 1, AlignmentVec1[i].Score, j + 1, AlignmentVec2[j].Score);

			if (AlignmentVec2[j].PairedAlnCanIdx == -1 || (AlignmentVec1[i].Score > AlignmentVec1[AlignmentVec2[j].PairedAlnCanIdx].Score))
			{
				AlignmentVec1[i].PairedAlnCanIdx = j;
				AlignmentVec2[j].PairedAlnCanIdx = i;

				if (s > max_score) max_score = s;
			}
		}
	}
	if (max_score > 0)
	{
		max_score = (int)(max_score*0.75);
		for (i = 0; i != num1; i++)
		{
			if (AlignmentVec1[i].PairedAlnCanIdx == -1 || (AlignmentVec1[i].Score + AlignmentVec2[AlignmentVec1[i].PairedAlnCanIdx].Score < max_score)) AlignmentVec1[i].Score = 0;
		}
		for (j = 0; j != num2; j++)
		{
			if (AlignmentVec2[j].PairedAlnCanIdx == -1 || (AlignmentVec1[AlignmentVec2[j].PairedAlnCanIdx].Score + AlignmentVec2[j].Score < max_score)) AlignmentVec2[j].Score = 0;
		}
		bPaired = true;
	}
}

void CheckPairingStatus(MappingReport_t& report1, vector<AlignmentCandidate_t>& AlignmentVec1, MappingReport_t& report2, vector<AlignmentCandidate_t>& AlignmentVec2)
{
	int i, j, dist, idx1, idx2;

	if (report1.score > 0 && report2.score > 0 && report1.Coor.ChromosomeIdx == report2.Coor.ChromosomeIdx && (report2.Coor.gPos - report1.Coor.gPos) < 1000000)
	{
		//printf("reads are paired aligned\n");
		report1.bPaired = report2.bPaired = true;
	}
	else
	{
		//printf("unpaired aligned --> %ld vs %ld\n", report1.Coor.gPos, report2.Coor.gPos);
		report1.bPaired = report2.bPaired = true;
	}
}

void RemovedUnPairedCandidates(vector<AlignmentCandidate_t>& AlignmentVec1, vector<AlignmentCandidate_t>& AlignmentVec2)
{
	int i, num;

	for (num = (int)AlignmentVec1.size(), i = 0; i < num; i++)
	{
		if (AlignmentVec1[i].PairedAlnCanIdx != -1 && AlignmentVec2[AlignmentVec1[i].PairedAlnCanIdx].Score == 0)
			AlignmentVec1[i].Score = 0;
	}

	for (num = (int)AlignmentVec2.size(), i = 0; i < num; i++)
	{
		if (AlignmentVec2[i].PairedAlnCanIdx != -1 && AlignmentVec1[AlignmentVec2[i].PairedAlnCanIdx].Score == 0)
			AlignmentVec2[i].Score = 0;
	}
}

void *ReadMapping(void *arg)
{
	ReadItem_t* ReadArr = NULL;
	MappingReport_t* ReportArr = NULL;
	int i, j, max_dist, ReadNum, EstDistance;
	vector<SeedPair_t> SeedPairVec1, SeedPairVec2;
	vector<AlignmentCandidate_t> AlignmentVec1, AlignmentVec2;

	ReadArr = new ReadItem_t[ReadChunkSize];  ReportArr = new MappingReport_t[ReadChunkSize];
	while (true)
	{
		pthread_mutex_lock(&LibraryLock);
		ReadNum = GetNextChunk(bSepLibrary, ReadFileHandler1, ReadFileHandler2, ReadArr);
		fprintf(stderr, "\r%d %s reads have been processed in %ld seconds...", iTotalReadNum, (bPairEnd? "paired-end":"singled-end"), (time(NULL) - StartProcessTime));
		pthread_mutex_unlock(&LibraryLock);
		
		if (ReadNum == 0) break;
		if (bPairEnd && ReadNum % 2 == 0)
		{
			for (i = 0, j = 1; i != ReadNum; i += 2, j += 2)
			{
				if (bDebugMode) printf("Mapping paired reads#%d %s (len=%d) and %s (len=%d):\n", i + 1, ReadArr[i].header + 1, ReadArr[i].rlen, ReadArr[j].header + 1, ReadArr[j].rlen);

				IdentifySeedPairs(ReadArr[i].rlen, ReadArr[i].EncodeSeq, SeedPairVec1); //if (bDebugMode) ShowSeedInfo(SeedPairVec1);
				GenerateAlignmentCandidate(ReadArr[i].rlen, SeedPairVec1, AlignmentVec1);

				IdentifySeedPairs(ReadArr[j].rlen, ReadArr[j].EncodeSeq, SeedPairVec2); //if (bDebugMode) ShowSeedInfo(SeedPairVec2);
				GenerateAlignmentCandidate(ReadArr[j].rlen, SeedPairVec2, AlignmentVec2);

				if (bDebugMode) ShowAlignmentCandidateInfo(ReadArr[i].header+1, AlignmentVec1), ShowAlignmentCandidateInfo(ReadArr[j].header+1, AlignmentVec2);
				if (!CheckPairedAlignmentCandidates(AlignmentVec1, AlignmentVec2)) RemoveRedundantCandidates(AlignmentVec1); RemoveRedundantCandidates(AlignmentVec2);

				if (bDebugMode)
				{
					ShowAlignmentCandidateInfo(ReadArr[i].header + 1, AlignmentVec1);
					ShowAlignmentCandidateInfo(ReadArr[j].header + 1, AlignmentVec2);
				}
				ReportArr[i] = GenMappingReport(true,  ReadArr[i], AlignmentVec1);
				ReportArr[j] = GenMappingReport(false, ReadArr[j], AlignmentVec2);

				CheckPairingStatus(ReportArr[i], AlignmentVec1, ReportArr[j], AlignmentVec2);
				// reverse read sequence if it is reversely mapped
				if (ReportArr[i].score > 0 && ReportArr[i].Coor.bDir == false)
				{
					string myseq = ReadArr[i].seq;
					GetComplementarySeq(ReadArr[i].rlen, (char*)myseq.c_str(), ReadArr[i].seq);
				}
				if (ReportArr[j].score > 0 && ReportArr[j].Coor.bDir == true)
				{
					string myseq = ReadArr[j].seq;
					GetComplementarySeq(ReadArr[j].rlen, (char*)myseq.c_str(), ReadArr[j].seq);
				}
				SetPairedAlignmentFlag(ReportArr[i], ReportArr[j]);
				EvaluateMAPQ(ReadArr[i].rlen, ReportArr[i]); EvaluateMAPQ(ReadArr[j].rlen, ReportArr[j]);

				if (bDebugMode) printf("\nEnd of mapping\n\n\n");
			}
		}
		else
		{
			for (i = 0; i != ReadNum; i++)
			{
				if (bDebugMode) printf("Mapping single read#%d %s (len=%d):\n", i + 1, ReadArr[i].header + 1, ReadArr[i].rlen);
				//fprintf(stdout, "%s\n", ReadArr[i].header + 1); fflush(output);

				IdentifySeedPairs(ReadArr[i].rlen, ReadArr[i].EncodeSeq, SeedPairVec1); //if (bDebugMode) ShowSeedInfo(SeedPairVec1);
				GenerateAlignmentCandidate(ReadArr[i].rlen, SeedPairVec1, AlignmentVec1);
				RemoveRedundantCandidates(AlignmentVec1); if (bDebugMode) ShowAlignmentCandidateInfo(ReadArr[i].header+1, AlignmentVec1);
				ReportArr[i] = GenMappingReport(true, ReadArr[i], AlignmentVec1);

				if (ReportArr[i].score > 0 && ReportArr[i].Coor.bDir == false)
				{
					string myseq = ReadArr[i].seq;
					GetComplementarySeq(ReadArr[i].rlen, (char*)myseq.c_str(), ReadArr[i].seq);
				}
				ReportArr[i].bPaired = false; ReportArr[i].iFrag = (ReportArr[i].score == 0 ? 0x4 : (ReportArr[i].Coor.bDir ? 0x0 : 0x10));
				EvaluateMAPQ(ReadArr[i].rlen, ReportArr[i]);

				if (bDebugMode) printf("\nEnd of mapping --> %ld\n%s\n", (ReportArr[i].score == 0 ? 0 : ReportArr[i].Coor.gPos), string().assign(100, '=').c_str());
			}
		}
		pthread_mutex_lock(&OutputLock);
		iTotalReadNum += ReadNum;
		if (bPairEnd && ReadNum %2 == 0) for (i = 0, j = 1; i != ReadNum; i += 2, j += 2) OutputPairedSamFile(ReadArr[i], ReportArr[i], ReadArr[j], ReportArr[j]);
		else for (i = 0; i != ReadNum; i++) OutputSingledSamFile(ReadArr[i], ReportArr[i]);
		pthread_mutex_unlock(&OutputLock);

		for (i = 0; i != ReadNum; i++)
		{
			delete[] ReadArr[i].header;
			delete[] ReadArr[i].seq;
			delete[] ReadArr[i].EncodeSeq;
		}
	}
	delete[] ReadArr; delete[] ReportArr;
	//fprintf(stderr, "\niSimpleLen=%.2f, iNormalLen=%.2f\n", 1.0*iSimpleLength / (iSimpleLength + iNormalLength), 1.0*iNormalLength / (iSimpleLength + iNormalLength));

	return (void*)(1);
}

void LoadExonInfo()
{
	int exon_len;
	int64_t gPos;
	fstream file;
	string chr, str;
	stringstream ss;

	file.open("300.error_free.exon.txt", ios_base::in);
	while (!file.eof())
	{
		getline(file, str); if (str == "") break;
		ss >> gPos >> exon_len;
		ExonMap.insert(make_pair(gPos, exon_len));
	}
	file.close();
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
		fprintf(output, "@PG\tPN:Kart\tVN:%s\n", VersionStr);
		for (i = 0; i < iChromsomeNum; i++) fprintf(output, "@SQ\tSN:%s\tLN:%ld\n", ChromosomeVec[i].name, ChromosomeVec[i].len);
	}

	//LoadExonInfo();
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
		if (bPairEnd) fprintf(stderr, "\t# of total mapped sequences = %d (sensitivity = %.2f%%)\n", iTotalReadNum - iUnMapped, (int)(10000 * (1.0*(iTotalReadNum - iUnMapped) / iTotalReadNum) + 0.5) / 100.0);
		else fprintf(stderr, "\t# of total mapped sequences = %d (sensitivity = %.2f%%)\n", iTotalReadNum - iUnMapped, (int)(10000 * (1.0*(iTotalReadNum - iUnMapped) / iTotalReadNum) + 0.5) / 100.0);
	}
	//map<int64_t, int>::const_iterator iter, ChrIter;
	//for (iter = ExonMap.begin(); iter != ExonMap.end(); iter++)
	//{
	//	printf("%ld %d\n", iter->first, iter->second);
	//}
}
