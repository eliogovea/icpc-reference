#include <bits/stdc++.h>

using namespace std;

int val1(string s) {
	char res = 0;
	for (int i = 0; i < s.size(); i++) {
		res ^= s[i];
	}
	return res;
}

const long long B1 = 357;
const long long M1 = 1000 * 1000 * 1000 + 7;

long long v2(string s) {
	long long res = 0;
	for (int i = 0; i < s.size(); i++) {
		res = (B1 * res + (long long)s[i]) % M1;
	}
	return res;
}

const long long B2 = 489;
const long long M2 = (1 << 27) + (1 << 25) + 1;

long long v3(string s) {
	long long res = 0;
	for (int i = 0; i < s.size(); i++) {
		res = (B2 * res + (long long)s[i]) % M2;
	}
	return res;
}

int n;
long long a[555], b[555], c[555];

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);

	string s;

	while (cin >> n && n) {
		getline(cin, s);
		int un = 0;
		int col = 0;
		for (int i = 0; i < n; i++) {
			getline(cin, s);
			a[i] = val1(s);
			b[i] = v2(s);
			c[i] = v3(s);

			bool ok = true;

			for (int j = 0; j < i; j++) {
				if (a[i] == a[j] && b[i] == b[j] && c[i] == c[j]) {
					ok = false;
				}
			}
			if (ok) {
				un++;
			}
			for (int j = 0; j < i; j++) {
				if (a[i] == a[j]) {
					if (b[i] != b[j] || c[i] != c[j]) {
						col++;
					}
				}
			}
		}
		cout << un << " " << col << "\n";
	}
}
