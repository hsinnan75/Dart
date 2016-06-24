#include "structure.h"

bwt_t *Refbwt;
bwaidx_t *RefIdx;
const char* VersionStr = "1.0.0";
bool bDebugMode, bPairEnd, FastQFormat;
int iThreadNum, MaxInsertSize, MaxGaps, MaxIntronSize;
char *RefSequence, *IndexFileName, *ReadFileName, *ReadFileName2, *SamFileName;

void ShowProgramUsage(const char* program)
{
	fprintf(stderr, "\nKart-RNA v%s [Developers: Hsin-Nan Lin and Wen-Lian Hsu]\n\n", VersionStr);
	fprintf(stderr, "Usage: %s aln|index\n\n", program);
	fprintf(stderr, "Command: index		index the reference sequences with FASTA format\n");
	fprintf(stderr, "         aln		read alignment\n");
	fprintf(stderr, "\n");
}

bool CheckReadFile(char* filename)
{
	FILE* f = fopen(filename, "r");

	if (f == NULL) return false;
	else return true;
}

bool CheckDataFormat(const char* fn)
{
	fstream file;
	string header;
	bool bCheck = true;

	file.open(fn, ios_base::in); getline(file, header);
	if (header == "" || (FastQFormat && header[0] != '@') || (!FastQFormat && header[0] != '>')) bCheck = false;
	file.close();

	return bCheck;
}

int main(int argc, char* argv[])
{
	int i;
	string parameter, str;

	MaxGaps = 5;
	iThreadNum = 16;
	bPairEnd = false;
	bDebugMode = false;
	MaxIntronSize = 500000;
	FastQFormat = true; // fastq:true, fasta:false
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
			else if (parameter == "-f")
			{
				ReadFileName = argv[++i];
				FastQFormat = false;
			}
			else if (parameter == "-q") ReadFileName = argv[++i];
			else if (parameter == "-f2" || parameter =="-q2") ReadFileName2 = argv[++i];
			else if (parameter == "-t")
			{
				if ((iThreadNum = atoi(argv[++i])) > 16)
				{
					fprintf(stderr, "Warning! Thread number is limited to 16!\n");
					iThreadNum = 16;
				}
			}
			else if (parameter == "-o") SamFileName = argv[++i];
			else if (parameter == "-pair" || parameter == "-p") bPairEnd = true;
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
			fprintf(stderr, "Kart-RNA v%s\n", VersionStr);
			fprintf(stderr, "Usage: %s aln [-i IndexFile Prefix] -f|-q ReadFile [-f2|-q2 ReadFile2] > out.sam\n\n", argv[0]);
			fprintf(stderr, "Options: -t INT\t	number of threads [16]\n");
			fprintf(stderr, "         -f|-f2\t	reads are in the FASTA format\n");
			fprintf(stderr, "         -q|-q2\t	reads are in the FASTQ format\n");
			fprintf(stderr, "         -p	\t	paired-end reads are interlaced in the same file\n");
			fprintf(stderr, "\n");
			exit(0);
		}
		if (CheckReadFile(ReadFileName) == false) fprintf(stderr, "Cannot open the read file: %s\n", ReadFileName), exit(0);
		if (ReadFileName2 != NULL && CheckReadFile(ReadFileName2) == false) fprintf(stderr, "Cannot open the read file: %s\n", ReadFileName2), exit(0);

		if (CheckDataFormat(ReadFileName) == false)
		{
			fprintf(stderr, "%s is not a %s file, please check again.\n", ReadFileName, (FastQFormat ? "fastq" : "fasta"));
			exit(0);
		}
		if (ReadFileName2 != NULL && CheckDataFormat(ReadFileName2) == false)
		{
			fprintf(stderr, "%s is not a %s file, please check again.\n", ReadFileName2, (FastQFormat ? "fastq" : "fasta"));
			exit(0);
		}


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
