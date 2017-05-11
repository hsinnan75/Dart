#include "structure.h"
#include <sys/stat.h>

bwt_t *Refbwt;
bwaidx_t *RefIdx;
char SJFileName[256];
const char* VersionStr = "1.0.0";
bool bDebugMode, bPairEnd, FastQFormat, bMultiHit, bUnique;
int iThreadNum, MaxInsertSize, MaxGaps, MaxIntronSize, iInsertSize;
const char* SpliceJunctionArr[4] = { "GT/AG", "CT/AC", "GC/AG", "CT/GC" };
char *RefSequence, *IndexFileName, *ReadFileName, *ReadFileName2, *SamFileName;

void ShowProgramUsage(const char* program)
{
	fprintf(stderr, "\nDART v%s [Developers: Hsin-Nan Lin and Wen-Lian Hsu]\n\n", VersionStr);
	fprintf(stderr, "Usage: %s aln|index\n\n", program);
	fprintf(stderr, "Command: index		index the reference sequences with FASTA format\n");
	fprintf(stderr, "         aln		read alignment\n");
	fprintf(stderr, "\n");
}

bool CheckReadFile(char* filename, bool bFirst)
{
	char ch;
	fstream file;
	bool bCheck = true;
	string str, header, seq;
	int iCount = 0, iBaseCount = 0;

	file.open(filename, ios_base::in);
	if (!file.is_open()) return false;
	else
	{
		getline(file, header);
		if (header == "") return false;
		else
		{
			if (bFirst)
			{
				if (header[0] == '>') FastQFormat = false;
				else FastQFormat = true;
			}
			else
			{
				if ((FastQFormat && header[0] != '@') || (!FastQFormat && header[0] != '>')) return false;
			}
		}
	}
	file.close();

	if (!bCheck) fprintf(stderr, "Read file: %s cannot be accessed or with incompatible format!\n", filename);

	return bCheck;
}

bool CheckSamFileName()
{
	struct stat s;
	bool bRet = true;

	if (stat(SamFileName, &s) == 0)
	{
		if (s.st_mode & S_IFDIR)
		{
			bRet = false;
			fprintf(stderr, "Warning: %s is a directory!\n", SamFileName);
		}
		else if (s.st_mode & S_IFREG)
		{
			//it's a file
		}
		else
		{
			bRet = false;
			fprintf(stderr, "Warning: %s is not a regular file!\n", SamFileName);
		}
	}
	return bRet;
}

int main(int argc, char* argv[])
{
	int i;
	string parameter, str;

	MaxGaps = 5;
	iInsertSize = 700;
	iThreadNum = 16;
	bPairEnd = false;
	bDebugMode = false;
	bMultiHit = false;
	bUnique = false;
	MaxIntronSize = 500000;
	FastQFormat = true; // fastq:true, fasta:false

	strcpy(SJFileName, "junctions.tab");
	RefSequence = IndexFileName = ReadFileName = ReadFileName2 = SamFileName = NULL;

	if (argc == 1) ShowProgramUsage(argv[0]);
	else if (strcmp(argv[1], "index") == 0)
	{
		string cmd(argv[0]);
		cmd = cmd.substr(0, cmd.find_last_of('/')+1) + "bwt_index";

		for (i = 2; i < argc; i++) cmd += " " + (string)argv[i];
		system((char*)cmd.c_str());
	}
	else if (strcmp(argv[1], "aln") == 0)
	{
		for (i = 2; i < argc; i++)
		{
			parameter = argv[i];

			if (parameter == "-i") IndexFileName = argv[++i];
			else if (parameter == "-f" || parameter == "-q") ReadFileName = argv[++i];
			else if (parameter == "-f2" || parameter == "-q2") ReadFileName2 = argv[++i];
			else if (parameter == "-t")
			{
				if ((iThreadNum = atoi(argv[++i])) > 16)
				{
					fprintf(stderr, "Warning! Thread number is limited to 16!\n");
					iThreadNum = 16;
				}
			}
			else if (parameter == "-o") SamFileName = argv[++i];
			else if (parameter == "-j") strcpy(SJFileName, argv[++i]);
			else if (parameter == "-pair" || parameter == "-p") bPairEnd = true;
			else if (parameter == "-m") bMultiHit = true;
			else if (parameter == "-unique") bUnique = true;
			else if (parameter == "-intron")
			{
				if ((MaxIntronSize = atoi(argv[++i])) < 100000) MaxIntronSize = 100000;
			}
			else if (parameter == "-ins")
			{
				if ((MaxInsertSize = atoi(argv[++i])) < 500) MaxInsertSize = 500;
			}
			else if (parameter == "-d" || parameter == "-debug") bDebugMode = true;
			else fprintf(stderr, "Warning! Unknow parameter: %s\n", argv[i]);
		}

		if (IndexFileName == NULL || ReadFileName == NULL)
		{
			fprintf(stderr, "\n");
			fprintf(stderr, "DART v%s\n", VersionStr);
			fprintf(stderr, "Usage: %s aln [-i IndexFile Prefix] -f|-q ReadFile [-f2|-q2 ReadFile2] > out.sam\n\n", argv[0]);
			fprintf(stderr, "Options: -t INT       number of threads [16]\n");
			fprintf(stderr, "         -f           files with #1 mates reads\n");
			fprintf(stderr, "         -f2          files with #2 mates reads\n");
			fprintf(stderr, "         -m           output multiple alignments\n");
			fprintf(stderr, "         -unique      output uniquely mapped alignments\n");
			fprintf(stderr, "         -o           sam output filename\n");
			fprintf(stderr, "         -j           splice junction output filename\n");
			fprintf(stderr, "         -p           paired-end reads are interlaced in the same file\n");
			fprintf(stderr, "\n");
			exit(0);
		}
		if (CheckReadFile(ReadFileName, true) == false) exit(0);
		if (ReadFileName2 != NULL && CheckReadFile(ReadFileName2, false) == false) exit(0);
		if (SamFileName != NULL && CheckSamFileName() == false) exit(0);

		if (CheckBWAIndexFiles(IndexFileName)) RefIdx = bwa_idx_load(IndexFileName);
		else RefIdx = 0;

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
	else ShowProgramUsage(argv[0]);

	return 0;
}
