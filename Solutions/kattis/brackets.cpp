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
	ios::sync_with_stdio(false);
	cin.tie(0);

	//freopen("dat.txt", "r", stdin);

	cin >> s;
	n = s.size();

	lg[1] = 0;
	for (int i = 2; i <= n; i++) {
		lg[i] = lg[i / 2] + 1;
	}

	bool good = true;
	for (int i = 1; i <= n; i++) {
		sum[i] = sum[i - 1];
		if (s[i - 1] == '(') {
			sum[i]++;
		} else {
			sum[i]--;
		}
		if (sum[i] < 0) {
			good = false;
		}
	}

	if (sum[n] != 0) {
		good = false;
	}

	if (good) {
		cout << "possible\n";
		return 0;
	}
	build();
	bool possible = false;
	for (int l = 1; l <= n && !possible; l++) {
		for (int r = l; r <= n && !possible; r++) {
			//cerr << "interaval: " << l << " " << r << "\n";
			/// [1, l)
			if (l > 1) {
				int minVal = queryMin(1, l - 1);
				//cerr << "[1, l) " << minVal << "\n";
				if (minVal < 0) {
					continue;
				}
			}
			/// [l, r]
			{
				int minVal = 2 * sum[l - 1] - queryMax(l, r);
				//cerr << "[l, r] " << minVal << "\n";
				if (minVal < 0) {
					continue;
				}
			}
			/// (r, n]
			if (r < n) {
				int minVal = 2 * (sum[l - 1] - sum[r]) + queryMin(r + 1, n);
				//cerr << "(r, n] " << minVal << "\n";
				if (minVal < 0) {
					continue;
				}
			}
			/// check sum
			int sumT;
			if (r < n) {
				sumT = 2 * (sum[l - 1] - sum[r]) + sum[n];
			} else {
				sumT = 2 * sum[l - 1] - sum[n];
			}

			//cerr << "sumT " << sumT << "\n";
			if (sumT != 0) {
				continue;
			}
			possible = true;
		}
	}
	if (possible) {
		cout << "possible\n";
	} else {
		cout << "impossible\n";
	}
}
