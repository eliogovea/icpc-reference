// Knuth Optimization
// f[i][j] = min(f[i][k] + f[k][j] + c[i][j]):
// [i][j] -> interval [i, j)
// sufficient condition:
// p[i][j - 1] <= p[i][j] <= p[i + 1][j]
// e.g. Optimal Binary Search Tree

#include <bits/stdc++.h>

using namespace std;

void solve(vector<int> v) {
    int n = v.size();
    vector<int> sv(n + 1);
    for (int i = 1; i <= n; i++) {
        sv[i] = sv[i - 1] + v[i - 1];
    }
    auto c = [&](int l, int r, int m) {
        assert(l <= m && m < r);
        return sv[r] - sv[l] - v[m];
    };  // sv[i] = v[0] + v[1] + ... + v[i - 1]
    for (int i = 0; i < n; i++) {
        f[i][i + 1] = 0;
        p[i][i + 1] = i;
    }
    for (int s = 2; s <= n; s++) {
        for (int l = 0; l + s <= n; l++) {
            int r = l + s;
            f[l][r] = -1;
            for (int m = p[l][r - 1]; m <= p[l + 1][r]; m++) {
                int x = f[l][m] + f[m + 1][r] + c(l, r, m);
                if (f[l][r] == -1 || x <= f[l][r]) {
                    f[l][r] = x;
                    p[l][r] = m;
                }
            }
        }
    }
    return f[0][n];
}

int main() {}
