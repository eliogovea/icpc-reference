// problem link http://coj.uci.cu/24h/problem.xhtml?pid=2643

#include <bits/stdc++.h>

using namespace std;

typedef complex<double> base;

const double PI = acos(-1.0);

void fft(vector<base>& P, bool invert) {
    int n = P.size();
    int ln = 0;
    while ((1 << ln) < n) {
        ln++;
    }
    assert((1 << ln) == n);
    for (int i = 0; i < n; i++) {
        int x = i;
        int y = 0;
        for (int j = 0; j < ln; j++) {
            y = (y << 1) | (x & 1);
            x >>= 1;
        }
        if (y < i) {
            swap(P[y], P[i]);
        }
    }
    for (int len = 2; len <= n; len <<= 1) {
        int half = len >> 1;
        double angle = 2 * PI / len;
        if (invert) {
            angle = -angle;
        }
        base step(cos(angle), sin(angle));
        for (int i = 0; i < n; i += len) {
            base w(1);
            for (int j = 0; j < half; j++) {
                base u = P[i + j];
                base v = w * P[i + j + half];
                P[i + j] = u + v;
                P[i + j + half] = u - v;
                w = w * step;
            }
        }
    }
    if (invert) {
        for (int i = 0; i < n; i++) {
            P[i] = P[i] / (double)n;
        }
    }
}

void fft(vector<vector<base> >& P, bool invert) {
    int n = P.size();
    int m = P[0].size();
    for (int r = 0; r < n; r++) {
        fft(P[r], invert);
    }
    for (int c = 0; c < m; c++) {
        vector<base> Q(n);
        for (int r = 0; r < n; r++) {
            Q[r] = P[r][c];
        }
        fft(Q, invert);
        for (int r = 0; r < n; r++) {
            P[r][c] = Q[r];
        }
    }
}

inline int fix(int n) {
    int x = 1;
    while (x <= n) {
        x <<= 1;
    }
    return x;
}

void solve() {
    int na, ma;
    cin >> na >> ma;

    vector<vector<int> > ha(na, vector<int>(ma));

    for (int r = 0; r < na; r++) {
        for (int c = 0; c < ma; c++) {
            cin >> ha[r][c];
        }
    }

    int nb, mb;
    cin >> nb >> mb;

    vector<vector<int> > hb(nb, vector<int>(mb));

    for (int r = 0; r < nb; r++) {
        for (int c = 0; c < mb; c++) {
            cin >> hb[r][c];
        }
    }

    int n = fix(na);
    int m = fix(ma);

    vector<vector<base> > Pa(n, vector<base>(m, 0));
    for (int r = 0; r < na; r++) {
        for (int c = 0; c < ma; c++) {
            Pa[r][c] = ha[r][c];
        }
    }

    vector<vector<base> > Pb(n, vector<base>(m, 0));
    for (int r = 0; r < nb; r++) {
        for (int c = 0; c < mb; c++) {
            Pb[nb - 1 - r][mb - 1 - c] = hb[r][c];
        }
    }

    fft(Pa, false);
    fft(Pb, false);
    for (int r = 0; r < n; r++) {
        for (int c = 0; c < m; c++) {
            Pa[r][c] = Pa[r][c] * Pb[r][c];
        }
    }
    fft(Pa, true);
    vector<vector<long long> > sha(na + 1, vector<long long>(ma + 1));
    for (int r = 1; r <= na; r++) {
        for (int c = 1; c <= ma; c++) {
            sha[r][c] = (long long)ha[r - 1][c - 1] * ha[r - 1][c - 1] +
                        sha[r - 1][c] + sha[r][c - 1] - sha[r - 1][c - 1];
        }
    }

    long long allb = 0;
    for (int r = 0; r < nb; r++) {
        for (int c = 0; c < mb; c++) {
            allb += (long long)hb[r][c] * hb[r][c];
        }
    }

    auto get = [&sha](int r1, int c1, int r2, int c2) {
        return sha[r2][c2] - sha[r2][c1] - sha[r1][c2] + sha[r1][c1];
    };

    const long long INF = 1e17;

    long long ans1 = INF;
    long long ans2 = INF;

    for (int r = nb - 1; r < na; r++) {
        for (int c = mb - 1; c < ma; c++) {
            long long v = get(r + 1 - nb, c + 1 - mb, r + 1, c + 1) + allb -
                          2LL * round((Pa[r][c].real()));
            if (v < ans1) {
                ans2 = ans1;
                ans1 = v;
            } else if (v < ans2) {
                ans2 = v;
            }
        }
    }

    cout << ans1 << " " << ans2 << "\n";
}
int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);

    solve();
}
