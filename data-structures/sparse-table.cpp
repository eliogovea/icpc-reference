#include <bits/stdc++.h>

using namespace std;

const int N = 5005;
const int LN = 15;
const int INF = 1e9;

string s;
int n;
int sum[N];
int maxST[LN][N];
int minST[LN][N];

int lg[N];

void build() {
	for (int i = 1; i <= n; i++) {
		maxST[0][i] = sum[i];
		minST[0][i] = sum[i];
	}
	for (int i = 1; (1 << i) <= n; i++) {
		for (int l = 1; l + (1 << i) - 1 <= n; l++) {
			maxST[i][l] = max(maxST[i - 1][l], maxST[i - 1][l + (1 << (i - 1))]);
			minST[i][l] = min(minST[i - 1][l], minST[i - 1][l + (1 << (i - 1))]);
		}
	}
}

inline int queryMax(int l, int r) {
	int len = r - l + 1;
	return max(maxST[lg[len]][l], maxST[lg[len]][r - (1 << lg[len]) + 1]);
}

inline int queryMin(int l, int r) {
	int len = r - l + 1;
	return min(minST[lg[len]][l], minST[lg[len]][r - (1 << lg[len]) + 1]);
}

int main() {
    
}
