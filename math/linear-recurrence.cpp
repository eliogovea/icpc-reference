#include <bits/stdc++.h>

using namespace std;

// linear recurrence solver
// f -> first elements of the secuence
// t -> transitions f[n] = -t[1] * f[n - 1] + -t[2] * f[n - 2] + ... + -t[k] *
// f[n - k] t[0] = 1 O(k * k * log(n)) ensure 0 <= f[i] < M, 0 <= t[i] < M !!!!!
// TODO precalculate the transitions O(k * k * log(MAXN))
int linearRecurrence(vector<int> f, vector<int> t, long long n, const int M) {
    int k = f.size();
    f.resize(2 * k);

    assert(t.size() == k + 1);
    vector<int> nt(2 * k + 1);

    while (n >= k) {
        // calc the next k elements of f using t
        for (int i = k; i < 2 * k; i++) {
            f[i] = 0;
            for (int j = 1; j <= k; j++) {
                f[i] += (long long)(M - t[j]) * f[i - j] % M;
                if (f[i] >= M) {
                    f[i] -= M;
                }
            }
        }

        // get the new transitions T(x) <- T(x) * T(-x);
        for (int i = 0; i <= 2 * k; i++) {
            nt[i] = 0;
        }

        for (int i = 0; i <= k; i++) {
            for (int j = 0; j <= k; j++) {
                int v = (long long)t[i] * t[j] % M;
                if (j & 1) {
                    v = M - v;
                }
                nt[i + j] += v;
                if (nt[i + j] >= M) {
                    nt[i + j] -= M;
                }
            }
        }

        for (int i = 0; i <= k; i++) {
            t[i] = nt[i << 1];
        }

        // update f depending on the parity of n
        for (int i = 0, o = (n & 1); i < k; i++) {
            f[i] = f[2 * i + o];
        }

        n >>= 1;
    }

    return f[n];
}

// accepted
void spoj_seq() {
    const int M = 1000 * 1000 * 1000;

    int t;
    cin >> t;

    while (t--) {
        int k;
        cin >> k;

        vector<int> f(k), t(k + 1);

        for (int i = 0; i < k; i++) {
            cin >> f[i];
        }

        t[0] = 1;
        for (int i = 1; i <= k; i++) {
            cin >> t[i];
            if (t[i] != 0) {
                t[i] = M - t[i];
            }
        }

        int n;
        cin >> n;
        n--;

        cout << linearRecurrence(f, t, n, M) << "\n";
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);

    spoj_seq();
}
