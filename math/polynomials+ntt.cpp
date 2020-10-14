#include <bits/stdc++.h>

using namespace std;

const int M = 7340033;
vector<int> root, invRoot;
bool ready = false;
inline void add(int& a, int b) {
    a += b;
    if (a >= M) a -= M;
}
inline int mul(int a, int b) { return (long long)a * b % M; }
inline int power(int x, int n) {
    int y = 1;
    while (n) {
        if (n & 1) y = mul(y, x);
        x = mul(x, x);
        n >>= 1;
    }
    return y;
}
void calcRoots() {
    if (ready) return;
    ready = true;
    int a = 2;
    while (power(a, (M - 1) / 2) != M - 1) a++;
    for (int l = 1; (M - 1) % l == 0; l <<= 1) {
        root.push_back(power(a, (M - 1) / l));
        invRoot.push_back(power(root.back(), M - 2));
    }
}
void transform(vector<int>& P, int n, bool invert) {
    if (!ready) calcRoots();
    int ln = 0;
    while ((1 << ln) < n) ln++;
    for (int i = 0; i < n; i++) {
        int x = i, y = 0;
        for (int j = 0; j < ln; j++) y = (y << 1) | (x & 1), x >>= 1;
        if (y < i) swap(P[y], P[i]);
    }
    for (int e = 1; (1 << e) <= n; e++) {
        int len = (1 << e), half = len >> 1;
        int step = invert ? invRoot[e] : root[e];
        for (int i = 0; i < n; i += len) {
            int w = 1;
            for (int j = 0; j < half; j++) {
                int u = P[i + j];
                int v = mul(P[i + j + half], w);
                P[i + j] = u;
                add(P[i + j], v);
                P[i + j + half] = u;
                add(P[i + j + half], M - v);
                w = mul(w, step);
            }
        }
    }
    if (invert) {
        int in = power(n, M - 2);
        for (int i = 0; i < n; i++) P[i] = mul(P[i], in);
    }
}

vector<int> mul(vector<int> P, vector<int> Q) {
    assert(P.size() > 0 && Q.size() > 0);
    int s = P.size() + Q.size() - 1, n = 1;
    while (n < s) n <<= 1;
    P.resize(n);
    Q.resize(n);
    if (P == Q)
        transform(P, n, false), Q = P;
    else
        transform(P, n, false), transform(Q, n, false);
    for (int i = 0; i < n; i++) P[i] = mul(P[i], Q[i]);
    transform(P, n, true);
    P.resize(s);
    return P;
}

vector<int> inverse(vector<int> P) {
    int s = P.size();
    assert(s > 0);
    if (P[0] == 0) return {};
    int n = 1;
    while (n < P.size()) n <<= 1;
    P.resize(n);
    vector<int> Q(2 * n), R(2 * n), S(2 * n);
    R[0] = power(P[0], M - 2);
    for (int k = 2; k <= n; k *= 2) {
        for (int i = 0; i < k; i++) S[i] = R[i];
        for (int i = 0; i < min(k, n); i++) Q[i] = P[i];
        for (int i = min(k, n); i < 2 * k; i++) Q[i] = 0;
        transform(S, 2 * k, false);
        transform(Q, 2 * k, false);
        for (int i = 0; i < 2 * k; i++) S[i] = mul(S[i], mul(S[i], Q[i]));
        transform(S, 2 * k, true);
        for (int i = 0; i < k; i++) add(R[i], R[i]), add(R[i], M - S[i]);
    }
    R.resize(s);
    return R;
}

vector<int> integral(vector<int> P) {
    assert(P.size() > 0);
    P.push_back(0);
    for (int i = P.size() - 1; i > 0; i--)
        P[i] = mul(P[i - 1], power(i, M - 2));
    P[0] = 0;
    return P;
}

vector<int> derivative(vector<int> P) {
    assert(P.size() > 0);
    if (P.size() == 1) return {0};
    for (int i = 0; i < P.size() - 1; i++) P[i] = mul(P[i + 1], i + 1);
    P.pop_back();
    return P;
}

vector<int> log(vector<int> P) {
    int s = P.size();
    assert(s > 0);
    assert(P[0] == 1);
    P = integral(mul(derivative(P), inverse(P)));
    P.resize(s);
    return P;
}

vector<int> exp(vector<int> P) {
    int s = P.size(), n = 1;
    assert(s > 0 && P[0] == 0);
    while (n < s) n <<= 1;
    vector<int> R({1});
    for (int k = 2; k <= n; k <<= 1) {
        vector<int> Q(k);
        R.resize(k);
        for (int i = 0; i < min(k, n); i++) Q[i] = P[i];
        auto logR = log(R);
        for (int i = 0; i < k; i++) add(Q[i], M - logR[i]);
        add(Q[0], 1);
        R = mul(R, Q);
    }
    R.resize(s);
    return R;
}

//  P ^ a, a is a real number
//  for log P[0] == 1, if P[0] != 1 transform
//  P into c * (x ^ d) * Q and Q[0] = 1
//  P ^ a = (c ^ a) * (x ^ (d * a)) * exp(a * log(Q))
// a = 1 / 2, P[0] == 1
vector<int> sqrt(vector<int> P) {
    int n = P.size(), inv2 = power(2, M - 2);
    P = log(P);
    for (int i = 0; i < n; i++) P[i] = mul(P[i], inv2);
    P = exp(P);
    P.resize(n);
    return P;
}
