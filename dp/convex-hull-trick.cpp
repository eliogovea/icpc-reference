// Convex Hull Trick
// DP[i] = min(DP[j] + b[j] * a[i]) : (j < i)
// sufficient condition: b[j] >= b[j + 1]
// (b[j] <= b[j + 1] para minimizar)

#include <bits/stdc++.h>

using namespace std;

const int N = 100005;

struct line {
    long long m, n;
    long long y(long long x) { return m * x + n; }
};

inline long double inter(line a, line b) {
    return (long double)(b.n - a.n) / (long double)(a.m - b.m);
}

struct ConvexHullTrick {
    line ch[N];
    int size;
    void clear() { size = 0; }
    void add(line l) {
        while (size > 1 && inter(l, ch[size - 1]) < inter(l, ch[size - 2]))
            size--;
        ch[size++] = l;
    }
    long long get_min(long long x) {
        int id = 0, lo = 1, hi = size - 1;
        while (lo <= hi) {
            int mid = (lo + hi) >> 1;
            if (ch[mid].y(x) < ch[mid - 1].y(x)) {
                id = mid;
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        return ch[id].y(x);
    }
} CH;

int n;
long long a[N], b[N];
long long dp[N];

int main() {
    cin >> n;
    for (int i = 0; i < n; i++) cin >> a[i];
    for (int i = 0; i < n; i++) cin >> b[i];
    CH.clear();
    CH.add((line){b[0], 0});
    for (int i = 1; i < n; i++) {
        dp[i] = CH.get_min(a[i]);
        CH.add((line){b[i], dp[i]});
    }
    cout << dp[n - 1] << "\n";
}
