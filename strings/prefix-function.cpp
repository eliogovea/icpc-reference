#include <bits/stdc++.h>

using namespace std;

// tested problem C https://codeforces.com/group/3qadGzUdR4/contest/101710
template <class T>
vector<int> prefixFunction(const T& s, int n) {
    vector<int> pi(n);
    for (int i = 1, l = 0; i < n; i++) {
        while (l > 0 && s[i] != s[l]) { l = pi[l - 1]; }
        if (s[i] == s[l]) { l++; }
        pi[i] = l;
    }
    return pi;
}

// tested problem J https://codeforces.com/group/3qadGzUdR4/contest/101710
// minimal string p, ppp..., contains s as prefix
// e.g. for abcabcab the answer if abc
template <class T>
int minimalStringPeriod(const T& s, int n) {
    auto pi = prefixFunction(s, n);
    int r = n;
    int l = pi[n - 1];
    while (true) {
        int lp = n - l;
        int x = lp - pi[lp - 1];
        if (lp % x != 0) { x = lp; }
        if (lp >= l) { r = min(r, x); }
        if (l == 0) { break; }
        l = pi[l - 1];
    }
    return r;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);

    string s;
    cin >> s;

    auto pi = prefixFunction(s, s.size());

    for (int i = 0; i < pi.size(); i++) {
        cout << pi[i];
        if (i + 1 < pi.size()) { cout << " "; }
    }
    cout << "\n";
}
