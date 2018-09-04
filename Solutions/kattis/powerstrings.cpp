#include <bits/stdc++.h>

using namespace std;

vector <int> prefix_function(const string & s) {
	int n = s.size();
	vector <int> pi(n);
	for (int i = 1, j = 0; i < n; i++) {
		while (j > 0 && s[i] != s[j]) {
			j = pi[j - 1];
		}
		if (s[i] == s[j]) {
			j++;
		}
		pi[i] = j;
	}
	return pi;
}

int calc_power(const string & s) {
	int n = s.size();
	int x = prefix_function(s).back();
	if (n % (n - x)) {
		return 1;
	}
	return n / (n - x);
}

void solve() {
	string s;
	while (cin >> s && s != ".") {
		cout << calc_power(s) << "\n";
	}
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);

	solve();
}
