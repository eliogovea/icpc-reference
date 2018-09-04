// Divide & Conquer Optimization
// dp[i][j] = min(dp[i - 1][k] + C[k][j]) : (k < j)
// sufficient condition: pos[i][j] <= pos[i][j + 1]

#include <bits/stdc++.h>

using namespace std;

const int N = 1005, K = 205, INF = 1e9;
int n, k, u[N][N], sum[N][N];
int dp[K][N], pos[K][N];

inline int cost(int l, int r) {
	return (sum[r][r] - sum[r][l - 1] -
					sum[l - 1][r] + sum[l - 1][l - 1]) >> 1;
}

void solve(int g, int l, int r, int pl, int pr) {
	if (l > r) return; int mid = (l + r) >> 1;
	dp[g][mid] = INF; pos[g][mid] = -1;
	for (int i = pl; i < mid && i <= pr; i++)
		if (dp[g][mid] > dp[g - 1][i] + cost(i + 1, mid)) {
			dp[g][mid] = dp[g - 1][i] + cost(i + 1, mid);
			pos[g][mid] = i;
		}
	solve(g, l, mid - 1, pl, pos[g][mid]);
	solve(g, mid + 1, r, pos[g][mid], pr);
}


int main() {
	cin >> n >> k;
	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++)
			cin >> u[i][j];
	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++)
			sum[i][j] = u[i][j] + sum[i - 1][j]
			+ sum[i][j - 1] - sum[i - 1][j - 1];
	for (int i = 1; i <= n; i++)
		dp[1][i] = cost(1, i), pos[1][i] = 0;
	int ans = dp[1][n];
	for (int i = 2; i <= k; i++)
		solve(i, i, n, i - 1, n),
		ans = min(ans, dp[i][n]);
	cout << ans << "\n";
}
