#include "structure.h"
#include <sys/stat.h>

bwt_t *Refbwt;
bwaidx_t *RefIdx;
char SJFileName[256];
const char* VersionStr = "1.1.0";
int iThreadNum, MaxInsertSize, MaxGaps, MaxIntronSize, OutputFileFormat;
bool bDebugMode, bPairEnd, FastQFormat, bMultiHit, bUnique, gzCompressed;
const char* SpliceJunctionArr[4] = { "GT/AG", "CT/AC", "GC/AG", "CT/GC" };
char *RefSequence, *IndexFileName, *ReadFileName, *ReadFileName2, *OutputFileName;

void ShowProgramUsage(const char* program)
{
	fprintf(stderr, "\nDART v%s\n\n", VersionStr);
	fprintf(stderr, "Usage: %s -i Index_Prefix -f ReadFile [-f2 ReadFile2]\n\n", program);
	fprintf(stderr, "Options: -t INT        number of threads [16]\n");
	fprintf(stderr, "         -f            files with #1 mates reads\n");
	fprintf(stderr, "         -f2           files with #2 mates reads\n");
	fprintf(stderr, "         -o            alignment filename for output [stdout]\n");
	fprintf(stderr, "         -j            splice junction output filename [junctions.tab]\n");
	fprintf(stderr, "         -m            output multiple alignments\n");
	fprintf(stderr, "         -p            paired-end reads are interlaced in the same file\n");
	fprintf(stderr, "         -unique       output unique alignments\n");
	fprintf(stderr, "         -intron       the maximal intron size [500000]\n");
	fprintf(stderr, "\n");
}

bool Check_gzInputFormat()
{
	string filename, FileExt;
	bool bCompressed = false;

	filename = ReadFileName; FileExt = filename.substr(filename.find_last_of('.') + 1);
	if (FileExt == "gz") bCompressed = true;

	return bCompressed;
}

bool CheckOutputFileName()
{
	struct stat s;
	bool bRet = true;
	string filename, FileExt;

	filename = OutputFileName; FileExt = filename.substr(filename.find_last_of('.') + 1);

	if (FileExt == "gz") OutputFileFormat = 1;

	//if (OutputFileFormat == 0) fprintf(stderr, "OutputFile = %s [format=sam]\n", OutputFileName);
	//else if (OutputFileFormat == 1) fprintf(stderr, "OutputFile = %s [format=sam.gz]\n", OutputFileName);
	//else fprintf(stderr, "OutputFile = %s [format=bam]\n", OutputFileName);

	if (stat(OutputFileName, &s) == 0)
	{
		if (s.st_mode & S_IFDIR)
		{
			bRet = false;
			fprintf(stderr, "Warning: %s is a directory!\n", OutputFileName);
		}
		else if (s.st_mode & S_IFREG)
		{
		}
		else
		{
			bRet = false;
			fprintf(stderr, "Warning: %s is not a regular file!\n", OutputFileName);
		}
	}
	return bRet;
}

bool CheckReadFile(char* filename, bool& bReadFormat)
{
	fstream file;
	string header;
	bool bCheck = true;

	file.open(filename, ios_base::in);
	if (!file.is_open()) return false;
	else
	{
		getline(file, header);
		if (header == "") return false;
		else
		{
			if (header[0] == '@') bReadFormat = true;
			else bReadFormat = false;
		}
	}
	file.close();

	return bCheck;
}

bool CheckCompressedFile(char* filename)
{
	gzFile file;
	bool bCheck = true;

	if ((file = gzopen(filename, "rb")) == Z_NULL) bCheck = false;
	gzclose(file);

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
	bMultiHit = false;
	bUnique = false;
	MaxIntronSize = 500000;
	OutputFileFormat = 0; // 0:sam 1:sam.gz
	FastQFormat = true; // fastq:true, fasta:false

	strcpy(SJFileName, "junctions.tab");
	RefSequence = IndexFileName = ReadFileName = ReadFileName2 = OutputFileName = NULL;

	if (argc == 1 || strcmp(argv[1], "-h") == 0) ShowProgramUsage(argv[0]);
	else
	{
		for (i = 1; i < argc; i++)
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
			else if (parameter == "-o") OutputFileName = argv[++i];
			else if (parameter == "-j") strcpy(SJFileName, argv[++i]);
			else if (parameter == "-p") bPairEnd = true;
			else if (parameter == "-m") bMultiHit = true;
			else if (parameter == "-unique") bUnique = true;
			else if (parameter == "-intron")
			{
				if ((MaxIntronSize = atoi(argv[++i])) < 100000) MaxIntronSize = 100000;
			}
			else if (parameter == "-d" || parameter == "-debug") bDebugMode = true;
			else
			{
				fprintf(stderr, "Warning! Unknow parameter: %s\n", argv[i]);
				ShowProgramUsage(argv[0]);
				exit(0);
			}
		}

		if (IndexFileName == NULL || ReadFileName == NULL)
		{
			fprintf(stderr, "Warning! Please specify a valid index prefix and read files!\n");
			ShowProgramUsage(argv[0]);
			exit(0);
		}
		gzCompressed = Check_gzInputFormat();
		if ((gzCompressed && CheckCompressedFile(ReadFileName) == false) ||
			(!gzCompressed && CheckReadFile(ReadFileName, FastQFormat) == false))
			fprintf(stderr, "Cannot open the read file: %s\n", ReadFileName), exit(0);

		if (ReadFileName2 != NULL)
		{
			bool FastQFormat2 = true;
			if ((gzCompressed && CheckCompressedFile(ReadFileName2) == false) ||
				(!gzCompressed && CheckReadFile(ReadFileName2, FastQFormat2) == false))
				fprintf(stderr, "Cannot open the read file: %s\n", ReadFileName2), exit(0);

			if (FastQFormat2 != FastQFormat) fprintf(stderr, "The input files are not with the same format!\n"), exit(0);
		}
		if (OutputFileName != NULL && CheckOutputFileName() == false) exit(0); 

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

	return 0;
}
