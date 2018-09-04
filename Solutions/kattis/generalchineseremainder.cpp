#include <bits/stdc++.h>

using namespace std;

using LL = long long;

// solve a * x + b * y = gcd(a, b)
LL egcd(LL a, LL b, LL & x, LL & y) {
	if (a == 0) {x = 0; y = 1; return b;}
	LL g = egcd(b % a, a, y, x);
	x -= (b / a) * y; return g;
}

pair <LL, LL> crt(LL r1, LL m1, LL r2, LL m2) {
	// TODO test
	LL d = __gcd(m1, m2);
	if (r1 % d != r2 % d) return {-1, -1};
	LL rd = r1 % d;
	r1 /= d; m1 /= d; r2 /= d; m2 /= d;
	if (m1 < m2) {swap(r1, r2); swap(m1, m2);}
	LL k = (r2 - r1) % m2; if (k < 0) k += m2;
	LL x, y; egcd(m1, m2, x, y);
	x %= m2; if (x < 0) x += m2;
	k *= x; k %= m2;
	return {m1 * m2 * d, (k * m1 + r1) * d + rd};
}

// https://open.kattis.com/problems/generalchineseremainder
void generalchineseremainder() {
	int n;
	cin >> n;
	while (n--) {
		int r1, r2, m1, m2;
		cin >> r1 >> m1 >> r2 >> m2;
		auto ans = crt(r1, m1, r2, m2);
		if (ans.first == -1) {
			cout << "no solution\n";
		} else {
			cout << ans.second << " " << ans.first << "\n";
		}
	}
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);
	generalchineseremainder();
}