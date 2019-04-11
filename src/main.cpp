#include "structure.h"
#include <sys/stat.h>

bwt_t *Refbwt;
bwaidx_t *RefIdx;
char SJFileName[256];
const char* VersionStr = "1.3.5";
vector<string> ReadFileNameVec1, ReadFileNameVec2;
char *RefSequence, *IndexFileName, *OutputFileName;
int iThreadNum, MaxInsertSize, MaxGaps, MaxIntronSize, OutputFileFormat;
bool bDebugMode, bSilent, bPairEnd, FastQFormat, bMultiHit, bUnique, gzCompressed;
const char* SpliceJunctionArr[4] = { "GT/AG", "CT/AC", "GC/AG", "CT/GC" };

void ShowProgramUsage(const char* program)
{
	fprintf(stdout, "\nDART v%s (Hsin-Nan Lin & Wen-Lian Hsu)\n\n", VersionStr);
	fprintf(stdout, "Usage: %s -i Index_Prefix -f <ReadFile_A1 ReadFile_B1 ...> [-f2 <ReadFile_A2 ReadFile_B2 ...>] -o|-bo Alignment_Output\n\n", program);
	fprintf(stdout, "Options: -t INT        number of threads [4]\n");
	fprintf(stdout, "         -f            files with #1 mates reads\n");
	fprintf(stdout, "         -f2           files with #2 mates reads\n");
	fprintf(stdout, "         -o            alignment filename in SAM format\n");
	fprintf(stdout, "         -bo           alignment filename in BAM format\n");
	fprintf(stdout, "         -j            splice junction output filename [junctions.tab]\n");
	fprintf(stdout, "         -m            output multiple alignments\n");
	fprintf(stdout, "         -p            paired-end reads are interlaced in the same file\n");
	fprintf(stdout, "         -unique       output unique alignments\n");
	fprintf(stdout, "         -intron       the maximal intron size [500000]\n");
	fprintf(stdout, "         -v            version\n");
	fprintf(stdout, "\n");
}

bool CheckOutputFileName()
{
	bool bRet = true;

	if (strcmp(OutputFileName, "output.sam") != 0)
	{
		struct stat s;
		string filename;

		filename = OutputFileName;
		if (stat(OutputFileName, &s) == 0)
		{
			if (s.st_mode & S_IFDIR)
			{
				bRet = false;
				fprintf(stdout, "Warning: %s is a directory!\n", OutputFileName);
			}
			else if (s.st_mode & S_IFREG)
			{
			}
			else
			{
				bRet = false;
				fprintf(stdout, "Warning: %s is not a regular file!\n", OutputFileName);
			}
		}
	}
	return bRet;
}

bool CheckInputFiles()
{
	struct stat s;
	bool bRet = true;

	for (vector<string>::iterator iter = ReadFileNameVec1.begin(); iter != ReadFileNameVec1.end(); iter++)
	{
		if (stat(iter->c_str(), &s) == -1)
		{
			bRet = false;
			fprintf(stderr, "Cannot access file:[%s]\n", (char*)iter->c_str());
		}
	}
	for (vector<string>::iterator iter = ReadFileNameVec2.begin(); iter != ReadFileNameVec2.end(); iter++)
	{
		if (stat(iter->c_str(), &s) == -1)
		{
			bRet = false;
			fprintf(stderr, "Cannot access file:[%s]\n", (char*)iter->c_str());
		}
	}
	return bRet;
}

int main(int argc, char* argv[])
{
	int i;
	string parameter, str;

	MaxGaps = 5;
	iThreadNum = 4;
	bPairEnd = false;
	bDebugMode = false;
	bMultiHit = false;
	bUnique = false;
	bSilent = false;
	MaxIntronSize = 500000;
	OutputFileName = (char*)"output.sam";
	OutputFileFormat = 0; // 0:sam 1:bam
	FastQFormat = true; // fastq:true, fasta:false

	strcpy(SJFileName, "junctions.tab");
	RefSequence = IndexFileName = NULL;

	if (argc == 1 || strcmp(argv[1], "-h") == 0) ShowProgramUsage(argv[0]);
	else if (strcmp(argv[1], "update") == 0)
	{
		system("git fetch; git merge origin/master master;make");
		exit(0);
	}
	else
	{
		for (i = 1; i < argc; i++)
		{
			parameter = argv[i];

			if (parameter == "-i") IndexFileName = argv[++i];
			else if (parameter == "-f")
			{
				while (++i < argc && argv[i][0] != '-') ReadFileNameVec1.push_back(argv[i]);
				i--;
			}
			else if (parameter == "-f2")
			{
				while (++i < argc && argv[i][0] != '-') ReadFileNameVec2.push_back(argv[i]);
				i--;
			}
			else if (parameter == "-t")
			{
				if ((iThreadNum = atoi(argv[++i])) > 16)
				{
					fprintf(stdout, "Warning! Thread number is limited to 16!\n");
					iThreadNum = 16;
				}
			}
			else if (parameter == "-o")
			{
				OutputFileFormat = 0;
				OutputFileName = argv[++i];
			}
			else if (parameter == "-bo")
			{
				OutputFileFormat = 1;
				OutputFileName = argv[++i];
			}
			else if (parameter == "-silent") bSilent = true;
			else if (parameter == "-j") strcpy(SJFileName, argv[++i]);
			else if (parameter == "-p") bPairEnd = true;
			else if (parameter == "-m") bMultiHit = true;
			else if (parameter == "-unique") bUnique = true;
			else if (parameter == "-intron")
			{
				if ((MaxIntronSize = atoi(argv[++i])) < 100000) MaxIntronSize = 100000;
			}
			else if (parameter == "-d" || parameter == "-debug") bDebugMode = true;
			else if (parameter == "-v" || parameter == "--version")
			{
				fprintf(stdout, "DART v%s\n\n", VersionStr);
				exit(0);
			}
			else
			{
				fprintf(stderr, "Error! Unknow parameter: %s\n", argv[i]);
				ShowProgramUsage(argv[0]);
				exit(1);
			}
		}

		if (ReadFileNameVec1.size() == 0)
		{
			fprintf(stderr, "Error! Please specify a valid read input!\n");
			ShowProgramUsage(argv[0]);
			exit(1);
		}
		if (ReadFileNameVec2.size() > 0 && ReadFileNameVec1.size() != ReadFileNameVec2.size())
		{
			fprintf(stderr, "Error! Paired-end reads input numbers do not match!\n");
			fprintf(stderr, "Read1:\n"); for (vector<string>::iterator iter = ReadFileNameVec1.begin(); iter != ReadFileNameVec1.end(); iter++) fprintf(stderr, "\t%s\n", (char*)iter->c_str());
			fprintf(stderr, "Read2:\n"); for (vector<string>::iterator iter = ReadFileNameVec2.begin(); iter != ReadFileNameVec2.end(); iter++) fprintf(stderr, "\t%s\n", (char*)iter->c_str());
			exit(1);
		}
		if (CheckInputFiles() == false || CheckOutputFileName() == false) exit(0);
		if (IndexFileName != NULL && CheckBWAIndexFiles(IndexFileName)) RefIdx = bwa_idx_load(IndexFileName);
		else
		{
			fprintf(stderr, "Error! Please specify a valid reference index!\n");
			ShowProgramUsage(argv[0]);
			exit(1);
		}
		if (RefIdx == 0) fprintf(stderr, "\n\nError! Index files are corrupt!\n");
		else
		{
			Refbwt = RefIdx->bwt;
			RestoreReferenceInfo();
			Mapping();
			bwa_idx_destroy(RefIdx);
			if (RefSequence != NULL) delete[] RefSequence;
		}
	}

	return 0;
}
