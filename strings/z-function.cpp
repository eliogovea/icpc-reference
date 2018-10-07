#include <bits/stdc++.h>

using namespace std;

// tested problem D https://codeforces.com/group/3qadGzUdR4/contest/101710
template <class T>
vector <int> zFunction(const T & s, int n) {
	int l = 0, r = 0;
	vector <int> z(n);
	for (int i = 1; i < n; i++) {
		if (i <= r) {
			z[i] = min(z[i - l], r - i + 1);
		}
		while (i + z[i] < n && s[i + z[i]] == s[z[i]]) {
			z[i]++;
		}
		if (i + z[i] - 1 > r) {
			l = i;
			r = i + z[i] - 1;
		}
	}
	return z;
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);

	string s;
	cin >> s;

	auto z = zFunction(s, s.size());

	for (int i = 0; i < z.size(); i++) {
		cout << z[i];
		if (i + 1 < z.size()) {
			cout << " ";
		}
	}
	cout << "\n";
}
