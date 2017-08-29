#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

typedef struct
{
	string chr;
	int start;
	int end;
} ExonItem_t;

typedef struct
{
	string chr;
	int start;
	int end;
} SJItem_t;

int iAnnotated = 0;
vector<SJItem_t> AnnSJVec, ReportSJVec;
map<string, pair<int, int> > ChrIdxMap;

bool CompBySJItem(const SJItem_t& p1, const SJItem_t& p2)
{
	if (p1.chr == p2.chr) return p1.start < p2.start;
	else return p1.chr < p2.chr;
}

void LoadExonInfo()
{
	fstream file;
	stringstream ss;
	string str, tmp;
	SJItem_t SJItem;

	file.open("junctions.txt", ios_base::in);
	if (!file.is_open()) fprintf(stderr, "error open file: junctions.txt\n"), exit(0);

	while(!file.eof())
	{
		getline(file, str); if (str == "") break;
		ss.clear(); ss.str(str); ss >> SJItem.chr >> SJItem.start >> SJItem.end;
		AnnSJVec.push_back(SJItem);
	}
	file.close();

	sort(AnnSJVec.begin(), AnnSJVec.end(), CompBySJItem);
	//printf("Get %d splice junctions\n", (int)AnnExonVec.size());
	//for (int i = 0; i < AnnExonVec.size(); i++) printf("%s %d %d\n", AnnExonVec[i].chr.c_str(), AnnExonVec[i].start, AnnExonVec[i].end);
	SJItem.chr = "NULL";
	for (vector<SJItem_t>::iterator iter = AnnSJVec.begin(); iter != AnnSJVec.end(); iter++)
	{
		if (iter->chr != SJItem.chr)
		{
			ChrIdxMap[SJItem.chr].second = iter - AnnSJVec.begin();
			ChrIdxMap[(SJItem.chr = iter->chr)].first = iter - AnnSJVec.begin();
		}
	}
	ChrIdxMap[SJItem.chr].second = AnnSJVec.size();
	//for (map<string, pair<int, int> >::iterator iter = ChrIdxMap.begin(); iter != ChrIdxMap.end(); iter++)
	//	printf("%s: %d - %d\n", iter->first.c_str(), iter->second.first, iter->second.second);
}

void GetReportedExonList(string filename)
{
	fstream file;
	stringstream ss;
	string str, tmp;
	SJItem_t SJItem;

	file.open(filename.c_str(), ios_base::in);
	while (!file.eof())
	{
		getline(file, str); if (str == "") break;
		ss.clear(); ss.str(str); ss >> SJItem.chr >> SJItem.start >> SJItem.end;
		ReportSJVec.push_back(SJItem);
	}
	file.close();
	//printf("Get %d splice junctions\n", (int)ReportExonVec.size());
	//for (int i = 0; i < 10; i++) printf("%s %d %d\n", ReportExonVec[i].chr.c_str(), ReportExonVec[i].start, ReportExonVec[i].end);
}

void CheckSpliceJunctions()
{
	SJItem_t SJItem;
	int i, j, idx1, idx2, num;

	for (num = (int)ReportSJVec.size(), i = 0; i < num; i++)
	{
		SJItem = ReportSJVec[i];
		idx1 = ChrIdxMap[SJItem.chr].first;
		idx2 = ChrIdxMap[SJItem.chr].second;
		//printf("check %s %d %d\n", ExonItem.chr.c_str(), ExonItem.start, ExonItem.end);
		for (j = idx1; j < idx2; j++)
		{
			if (abs(SJItem.start - AnnSJVec[j].start) < 5 && abs(SJItem.end - AnnSJVec[j].end) < 5)
			{
				//printf("SJ: [%d-%d] vs [%d-%d]\n", ExonItem.start, ExonItem.end, AnnExonVec[j].start, AnnExonVec[j].end);
				iAnnotated++;
				break;
			}
		}
	}
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s SJ_file\n", argv[0]);
		exit(0);
	}
	LoadExonInfo();
	GetReportedExonList(argv[1]);
	CheckSpliceJunctions();

	printf("# of SJ = %d\n# of Reported SJ = %d\nAcc = %d (%.2f%%)\n", (int)AnnSJVec.size(), (int)ReportSJVec.size(), iAnnotated, (int)(10000*(1.0*iAnnotated/(int)ReportSJVec.size()))/100.0);

	return 0;
}
