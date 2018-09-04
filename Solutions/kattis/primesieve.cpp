#include <bits/stdc++.h>

using namespace std;

typedef long long LL;

inline void add(LL & a, LL b, LL M) {
	a += b; if (a >= M) a -= M;
}

inline LL mul(LL a, LL b, LL m) {
	// return (__int128)a * b % m;
	// a %= m; b %= m;
	if (m <= 2e9) return a * b % m;
	LL q = (long double)a * b / m;
	LL r = a * b - q * m; r %= m;
	if (r < 0) r += m;
	return r;
} // to avoid overflow, m < 1e18

inline LL power(LL x, LL n, LL m) {
	if (x == 0 || m == 1) return 0;
	LL y = 1 % m; x %= m;
	while (n > 0) {
		if (n & 1) y = mul(y, x, m);
		x = mul(x, x, m); n >>= 1;
	}
	return y;
}

// solve a * x + b * y = gcd(a, b)
LL egcd(LL a, LL b, LL & x, LL & y) {
	if (a == 0) {x = 0; y = 1; return b;}
	LL g = egcd(b % a, a, y, x);
	x -= (b / a) * y; return g;
}

LL inverse(LL n, LL m) {
	LL x, y; LL g = egcd(n, m, x, y);
	if (g != 1) return -1;
	x %= m; if (x < 0) x += m;
	assert((n * x % m) == 1); return x;
}

bool ready = false;
const int P = 1000 * 1000; 
bool _isPrime[P];
int minPrime[P];
vector <int> primes;

inline bool fastSieve() { // O(n)
	for (int i = 0; i < P; i++) 
		_isPrime[i] = true;
	_isPrime[0] = _isPrime[1] = false;
	primes.reserve(P); //
	for (int i = 2; i < P; i++) {
		if (_isPrime[i]) primes.push_back(i), minPrime[i] = i;
		for (int j = 0; j < primes.size() && primes[j] * i < P; j++) {
			_isPrime[primes[j] * i] = false;
			minPrime[primes[j] * i] = primes[j];
			if (i % primes[j] == 0) break;
		}
	}
}

// Miller Rabin primality test !!!
inline bool witness(LL x, LL n, LL s, LL d) {
	LL cur = power(x, d, n);
	if (cur == 1) return false;
	for (int r = 0; r < s; r++) {
		if (cur == n - 1) return false;
		cur = mul(cur, cur, n);
	}
	return true;
}

bool isPrime(long long n) { 
	if (!ready) fastSieve(), ready = true; // !!!
	if (n < P) return _isPrime[n];
	if (n < 2) return false;
	for (int x : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29}) {
		if (n == x) return true;
		if (n % x == 0) return false;
	}
	if (n < 31 * 31) return true;
	int s = 0; LL d = n - 1;
	while (!(d & 1)) s++, d >>= 1;
	// for n int: test = {2, 7, 61}
	// for n < 3e18: test = {2, 3, 5, 7, 11, 13, 17, 19, 23}
	static vector <LL> test = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
	for (long long x : test) {
		if (x % n == 0) return true;
		if (witness(x, n, s, d)) return false;
	}
//		for (auto p : primes) { // slow for testing
//			if ((long long)p * p > n) break;
//			if (n % p == 0) return false;
//		}
	return true;
} // ends Miller Rabin primality test !!!

// O(n ^ 0.75) same idea for other functions over primes
LL countPrimes(LL n) {
	int r = sqrt((double)n); 
	while ((LL)(r + 1) * (r + 1) <= n) r++; 
	vector <LL> values; // all [n / x]
	for (LL x = 1; x <= n; x = n / (n / x) + 1) 
		values.push_back(n / x);
	function <int(LL)> getPos = [&](LL x) {
		int pos = (x > r) ? n / x - 1 : values.size() - x;
		/* assert(values[pos] == x); */ return pos; 
	};
	// auto primes = getPrimes(r);
	if (!ready) fastSieve(), ready = true;
	vector <LL> cnt(values.size()); 
	for (auto v : values) cnt[getPos(v)] = v - 1;
	for (auto p : primes) for (auto v : values) {
		if ((LL)p * p > v) break;
		cnt[getPos(v)] -= cnt[getPos(v / p)] - cnt[getPos(p - 1)];
	}
	return cnt[getPos(n)];
} // ends prime counting

void primesieve() {
	int n, q;
	cin >> n >> q;
	
	cout << countPrimes(n) << "\n";
	while (q--) {
		int x;
		cin >> x;
		cout << isPrime(x) << "\n";
	}
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);
	primesieve();
}