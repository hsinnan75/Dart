#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <string>

using namespace std;

void IdentifyGenomicRegion(string& header, string& r_chr, int64_t& left_gPos, int64_t& right_gPos)
{
	string tmp;
	int p1, p2, p3;

	p1 = header.find_first_of(':');
	p2 = header.find_first_of('-');
	p3 = header.find_first_of('W');

	r_chr = header.substr(0, p1);
	tmp = header.substr(p1 + 1, p2 - p1 - 1); left_gPos = atoi(tmp.c_str());
	tmp = header.substr(p2 + 1, p3 - p2 + 1); right_gPos = atoi(tmp.c_str());

	//printf("header:%s --> chr=%s, l=%ld, r=%ld\n", header.c_str(), r_chr.c_str(), left_gPos, right_gPos);
}

int main(int argc, char* argv[])
{
	fstream file;
	float acc, sen;
	stringstream ss;
	int iTotal, iCor, mapq, iEmpty;
	int64_t gPos, left_gPos, right_gPos;
	string str, header, tmp, p_chr, r_chr, CIGAR;

	file.open(argv[1], ios_base::in); iTotal = iCor = iEmpty = 0;
	while (!file.eof())
	{
		getline(file, str); if (str == "") break;
		if (str[0] == '@') continue;

		ss.clear(); ss.str(str); ss >> header >> tmp >> p_chr >> gPos >> mapq >> CIGAR;
		IdentifyGenomicRegion(header, r_chr, left_gPos, right_gPos);

		if (mapq == 0) continue;
		iTotal++;
		if (p_chr == "*") iEmpty++;
		else if (CIGAR == "*") printf("%s\n", str.c_str());

		if (p_chr == r_chr && gPos >= left_gPos && gPos <= right_gPos) iCor++;
	}
	file.close();

	if (iTotal > 0) sen = (int)(1000 * (1.0*(iTotal - iEmpty) / iTotal + 0.0005)) / 10.0;
	else sen = 0;

	if (iTotal > 0) acc = (int)(1000 * (1.0*iCor / (iTotal - iEmpty) + 0.0005)) / 10.0;
	else acc = 0;

	printf("Sen = %d / %d = %.2f\nAcc = %d / %d = %.2f\n", iTotal - iEmpty, iTotal, sen, iCor, (iTotal - iEmpty), acc);

	return 0;
}
