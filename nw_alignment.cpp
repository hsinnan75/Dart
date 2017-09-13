#include "structure.h"

const short MaxPenalty = -32768;
const short OPEN_GAP = -1;
const short EXTEND_GAP = -1;
const short NEW_GAP = -2; //OPEN_GAP + EXTEND_GAP;

double max(short x, short y)
{
	return x > y ? x : y;
}

double max(short x, short y, short z)
{
	return x > y ? max(x, z) : max(y, z);
}

void nw_alignment(int m, string& s1, int n, string& s2)
{
	int i, j;

	m = m + 1, n = n + 1;

	short** r = new short*[m];
	short** t = new short*[m];
	short** s = new short*[m];

	for (i = 0; i < m; i++)
	{
		r[i] = new short[n];
		t[i] = new short[n];
		s[i] = new short[n];
	}

	// initialization
	r[0][0] = t[0][0] = s[0][0] = 0;
	for (i = 1; i < m; i++)
	{
		r[i][0] = MaxPenalty;
		s[i][0] = t[i][0] = 0;
	}

	for (j = 1; j < n; j++)
	{
		t[0][j] = MaxPenalty;
		s[0][j] = r[0][j] = 0;
	}

	for (i = 1; i < m; i++)
	{
		for (j = 1; j < n; j++)
		{
			r[i][j] = max(r[i][j - 1] + EXTEND_GAP, s[i][j - 1] + NEW_GAP);
			t[i][j] = max(t[i - 1][j] + EXTEND_GAP, s[i - 1][j] + NEW_GAP);
			//if (s1[i - 1] == 'N' || s2[j - 1] == 'N') s[i][j] = max(s[i - 1][j - 1] + 1, r[i][j], t[i][j]);
			//else s[i][j] = max(s[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 1 : -1), r[i][j], t[i][j]);
			s[i][j] = max(s[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 1 : -1), r[i][j], t[i][j]);
		}
	}

	i = m - 1, j = n - 1;
	while (i > 0 || j > 0) {
		if (s[i][j] == r[i][j]) {
			s1.insert(i, 1, '-');
			j--;
		}
		else if (s[i][j] == t[i][j]) {
			s2.insert(j, 1, '-');
			i--;
		}
		else {
			i--, j--;
		}
	}

	for (i = 0; i < m; i++)
	{
		delete[] r[i]; delete[] t[i]; delete[] s[i];
	}

	delete[] r; delete[] t; delete[] s;
}
