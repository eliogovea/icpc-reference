#include <bits/stdc++.h>

using namespace std;

template <class T>
vector<int> normalize(vector<T> v) {
    auto cv = v;
    vector<int> nv(v.size());
    sort(cv.begin(), cv.end());
    cv.erase(unique(cv.begin(), cv.end()), cv.end());
    for (int i = 0; i < v.size(); i++)
        nv[i] = upper_bound(cv.begin(), cv.end(), v[i]) - cv.begin();
    return nv;
}

// 0 < *min_element(v.begin(), v.end()) !!!
vector<int> get_suffix_array(vector<int> v, int alpha = -1) {
    if (alpha == -1) alpha = v.size() + 1;
    v.push_back(0);
    int n = v.size(), classes = alpha;
    vector<int> cnt(max(n, alpha));
    vector<int> p(n), np(n), c(n), nc(n);
    for (int i = 0; i < n; i++) p[i] = i, c[i] = v[i];
    for (int len = 1; len < 2 * n; len <<= 1) {
        int hlen = len >> 1;
        for (int i = 0; i < n; i++) np[i] = (p[i] - hlen + n) % n;
        for (int i = 0; i < classes; i++) cnt[i] = 0;
        for (int i = 0; i < n; i++) cnt[c[i]]++;
        for (int i = 1; i < classes; i++) cnt[i] += cnt[i - 1];
        for (int i = n - 1; i >= 0; i--) p[--cnt[c[np[i]]]] = np[i];
        classes = 0;
        for (int i = 0; i < n; i++) {
            if (i == 0 || c[p[i]] != c[p[i - 1]] ||
                c[(p[i] + hlen) % n] != c[(p[i - 1] + hlen) % n])
                classes++;
            nc[p[i]] = classes - 1;
        }
        for (int i = 0; i < n; i++) c[i] = nc[i];
    }
    p.erase(p.begin());
    return p;
}

vector<int> get_lcp(vector<int>& v, vector<int>& sa) {
    int n = v.size();
    vector<int> rank(n), lcp(n);
    for (int i = 0; i < n; i++) rank[sa[i]] = i;
    for (int i = 0, l = 0; i < n; i++) {
        if (rank[i] == n - 1) {
            l = 0;
            continue;
        }
        l = max(0, l - 1);
        int j = sa[rank[i] + 1];
        while (i + l < n && j + l < n && v[i + l] == v[j + l]) l++;
        lcp[rank[i]] = l;
    }
    return lcp;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);

    vector<int> v = {3, 1, 2, 3, 4, 3, 1};
    auto sa = get_suffix_array(v, -1);
    auto lcp = get_lcp(v, sa);
}
