#include <bits/stdc++.h>

using namespace std;

typedef complex<double> base;

const double PI = acos(-1);
const int N = 1 << 19;

void fft(vector<base> &a, bool invert) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        while (j >= bit) {
            j -= bit;
            bit >>= 1;
        }
        j += bit;
        if (i < j) { swap(a[i], a[j]); }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double angle = 2.0 * PI / len;
        if (invert) { angle = -angle; }
        base wlen(cos(angle), sin(angle));
        for (int i = 0; i < n; i += len) {
            base w(1);
            for (int j = 0; j < len / 2; j++) {
                base u = a[i + j];
                base v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (invert) {
        for (int i = 0; i < n; i++) { a[i] /= n; }
    }
}

void square(vector<long long> &p) {
    int len = 1;
    while (len < p.size()) { len *= 2; }
    len *= 2;
    p.resize(len);
    vector<base> fp(p.begin(), p.end());
    fft(fp, false);
    for (int i = 0; i < len; i++) { fp[i] *= fp[i]; }
    fft(fp, true);
    for (int i = 0; i < len; i++) { p[i] = (long long)(fp[i].real() + 0.5); }
}

vector<long long> cdiv(N, 0);

int st[4 * N];

void build(int x, int l, int r) {
    if (l == r) {
        st[x] = l;
    } else {
        int mid = (l + r) >> 1;
        build(2 * x, l, mid);
        build(2 * x + 1, mid + 1, r);
        if (cdiv[st[2 * x]] > cdiv[st[2 * x + 1]]) {
            st[x] = st[2 * x];
        } else if (cdiv[st[2 * x + 1]] > cdiv[st[2 * x]]) {
            st[x] = st[2 * x + 1];
        } else {
            st[x] = min(st[2 * x], st[2 * x + 1]);
        }
    }
}

int query(int x, int l, int r, int ql, int qr) {
    if (l > qr || r < ql) { return -1; }
    if (l >= ql && r <= qr) { return st[x]; }
    int mid = (l + r) >> 1;
    int q1 = query(2 * x, l, mid, ql, qr);
    int q2 = query(2 * x + 1, mid + 1, r, ql, qr);
    if (q1 == -1) { return q2; }
    if (q2 == -1) { return q1; }
    if (cdiv[q1] > cdiv[q2]) { return q1; }
    if (cdiv[q2] > cdiv[q1]) { return q2; }
    return min(q1, q2);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    for (int i = 1; i <= N; i++) {
        for (int j = i; j <= N; j += i) { cdiv[j]++; }
    }
    square(cdiv);
    build(1, 1, N - 1);
    int t;
    cin >> t;
    while (t--) {
        int a, b;
        cin >> a >> b;
        int pos = query(1, 1, N - 1, a, b);
        cout << pos << " " << cdiv[pos] << "\n";
    }
}
