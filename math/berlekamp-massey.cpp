#include <bits/stdc++.h>

using namespace std;

// Given the first elements of a sequence
// returns its recursive equation
// M prime
template <const int M>
struct BerlekampMassey {
    inline int power(int x, int n) {
        int y = 1 % M;
        while (n) {
            if (n & 1) { y = (long long)y * x % M; }
            x = (long long)x * x % M;
            n >>= 1;
        }
        return y;
    }
    inline bool inverse(int x) { return power(x, M - 2); }
    inline vector<int> shift(vector<int>& P, int d) {
        vector<int> Q(d + (int)P.size());
        for (int i = 0; i < (int)P.size(); i++) Q[i + d] = P[i];
        return Q;
    }
    int calc(vector<int>& P, vector<int>& d, int pos) {
        int res = 0;
        for (int i = 0; i < (int)P.size(); i++) {
            res += (long long)d[pos - i] * P[i] % M;
            if (res >= M) res -= M;
        }
        return res;
    }
    vector<int> sub(vector<int> P, const vector<int>& Q) {
        if (Q.size() > P.size()) P.resize(Q.size());
        for (int i = 0; i < (int)Q.size(); i++) {
            P[i] -= Q[i];
            if (P[i] < 0) P[i] += M;
        }
        return P;
    }
    vector<int> scale(const int c, vector<int> P) {
        for (auto& x : P) x = (long long)x * c % M;
        return P;
    }
    vector<int> solve(vector<int> f) {
        int n = f.size();
        vector<int> s(1, 1), b(1, 1);
        for (int i = 1, j = 0, ld = f[0]; i < n; i++) {
            int d = calc(s, f, i);
            if (d == 0) continue;
            int c = (long long)d * inverse(ld) % M;
            if (((int)s.size() - 1) * 2 <= i) {
                auto ob = b;
                b = s;
                s = sub(s, scale(c, shift(ob, i - j)));
                j = i;
                ld = d;
            } else
                s = sub(s, scale(c, shift(b, i - j)));
        }
        while (s.size() > 0 && s.back() == 0) {
            s.pop_back();  // ???
        }
        return s;
    }
};

int main() {
    // TODO test
}
