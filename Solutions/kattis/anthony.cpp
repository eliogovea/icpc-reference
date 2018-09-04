/// https://open.kattis.com/contests/z7x6z7/problems/anthony

#include <bits/stdc++.h>

using namespace std;

int n, m;
double p[2500];

double dp[1500][1500];

double calc(int a, int b) {
	if (a == n && b == m) {
		return 1.0;
	}
	if (dp[a][b] > -0.5) {
		return dp[a][b];
	}

}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);

	cin >> n >> m;
	for (int i = 0; i < n + m - 1; i++) {
		cin >> p[i];
	}

	dp[0][0] = 1;
	for (int x = 0; x < n; x++) {
		for (int y = 0; y < m; y++) {
			dp[x + 1][y] += dp[x][y] * (1.0 - p[x + y]);
			dp[x][y + 1] += dp[x][y] * p[x + y];
		}
	}
	double answer = 0.0;
	for (int x = 0; x < n; x++) {
		answer += dp[x][m];
	}
	cout.precision(14);
	cout << fixed << answer << "\n";
}
