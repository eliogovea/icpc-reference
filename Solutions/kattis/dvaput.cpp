#include <bits/stdc++.h>

using namespace std;

vector <int> get_suffix_array(const string & s, const int alpha = 256) {
	int n = s.size() + 1, classes = alpha;
	vector <int> cnt(max(n, alpha)), p(n), np(n), c(n), nc(n);
	for (int i = 0; i < n; i++) p[i] = i, c[i] = s[i];
	for (int len = 1; len < 2 * n; len <<= 1) {
		int hlen = len >> 1;
		for (int i = 0; i < n; i++) np[i] = (p[i] - hlen + n) % n;
		for (int i = 0; i < classes; i++) cnt[i] = 0;
		for (int i = 0; i < n; i++ ) cnt[c[i]]++;
		for (int i = 1; i < classes; i++) cnt[i] += cnt[i - 1];
		for (int i = n - 1; i >= 0; i--) p[--cnt[c[np[i]]]] = np[i];
		classes = 0;
		for (int i = 0; i < n; i++) {
			if (i == 0 || c[p[i]] != c[p[i - 1]] || c[(p[i] + hlen) % n] != c[(p[i - 1] + hlen) % n])
				classes++;
			nc[p[i]] = classes - 1;
		}
		for (int i = 0; i < n; i++) c[i] = nc[i];
	}
	for (int i = 0; i < n - 1; i++) p[i] = p[i + 1];
	p.pop_back();
	return p;
}

vector <int> get_lcp(const string & s, const vector <int> & sa) {
	int n = s.size(); vector <int> rank(n), lcp(n);
	for (int i = 0; i < n; i++) rank[sa[i]] = i;
	for (int i = 0, l = 0; i < n; i++) {
		if (rank[i] == n - 1) {l = 0; continue;}
		l = max(0, l - 1); int j = sa[rank[i] + 1];
		while (i + l < n && j + l < n && s[i + l] == s[j + l]) l++;
		lcp[rank[i]] = l;
	}
	return lcp;
}

void check(const string & s, const vector <int> & sa, const vector <int> & lcp) {
	int n = s.size();
	for (int i = 0; i < n; i++) {
		cerr << s.substr(sa[i], n - sa[i]) << " " << lcp[i] << "\n";
	}
}

void solve() {
	int n;
	string s;

	cin >> n >> s;
	vector <int> sa = get_suffix_array(s);
	vector <int> lcp = get_lcp(s, sa);

	// check(s, sa, lcp);

	cout << *max_element(lcp.begin(), lcp.end()) << "\n";
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);

	solve();
}
