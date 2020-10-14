#include <bits/stdc++.h>

using namespace std;

typedef long long LL;

inline void add(LL& a, LL b, LL M) {
    a += b;
    if (a >= M) a -= M;
}

inline LL mul(LL a, LL b, LL m) {
    // return (__int128)a * b % m;
    // a %= m; b %= m;
    if (m <= 2e9) return a * b % m;
    LL q = (long double)a * b / m;
    LL r = a * b - q * m;
    r %= m;
    if (r < 0) r += m;
    return r;
}  // to avoid overflow, m < 1e18

inline LL power(LL x, LL n, LL m) {
    if (x == 0 || m == 1) return 0;
    LL y = 1 % m;
    x %= m;
    while (n > 0) {
        if (n & 1) y = mul(y, x, m);
        x = mul(x, x, m);
        n >>= 1;
    }
    return y;
}

// solve a * x + b * y = gcd(a, b)
LL egcd(LL a, LL b, LL& x, LL& y) {
    if (a == 0) {
        x = 0;
        y = 1;
        return b;
    }
    LL g = egcd(b % a, a, y, x);
    x -= (b / a) * y;
    return g;
}

LL inverse(LL n, LL m) {
    LL x, y;
    LL g = egcd(n, m, x, y);
    if (g != 1) return -1;
    x %= m;
    if (x < 0) x += m;
    assert((n * x % m) == 1);
    return x;
}

bool ready = false;
const int P = 1000 * 1000;
bool _isPrime[P];
int minPrime[P];
vector<int> primes;

inline bool fastSieve() {  // O(n)
    for (int i = 0; i < P; i++) _isPrime[i] = true;
    _isPrime[0] = _isPrime[1] = false;
    primes.reserve(P);  //
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
    if (!ready) fastSieve(), ready = true;  // !!!
    if (n < P) return _isPrime[n];
    if (n < 2) return false;
    for (int x : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29}) {
        if (n == x) return true;
        if (n % x == 0) return false;
    }
    if (n < 31 * 31) return true;
    int s = 0;
    LL d = n - 1;
    while (!(d & 1)) s++, d >>= 1;
    // for n int: test = {2, 7, 61}
    // for n < 3e18: test = {2, 3, 5, 7, 11, 13, 17, 19, 23}
    static vector<LL> test = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    for (long long x : test) {
        if (x % n == 0) return true;
        if (witness(x, n, s, d)) return false;
    }
    //		for (auto p : primes) { // slow for testing
    //			if ((long long)p * p > n) break;
    //			if (n % p == 0) return false;
    //		}
    return true;
}  // ends Miller Rabin primality test !!!

LL nextPrime(LL x) {
    for (x = max(2LL, x + 1); x % 6; x++)
        if (isPrime(x)) return x;
    for (; true; x += 6) {
        if (isPrime(x + 1)) return x + 1;
        if (isPrime(x + 5)) return x + 5;
    }
}

// O(n ^ 0.75) same idea for other functions over primes
LL countPrimes(LL n) {
    int r = sqrt((double)n);
    while ((LL)(r + 1) * (r + 1) <= n) r++;
    vector<LL> values;  // all [n / x]
    for (LL x = 1; x <= n; x = n / (n / x) + 1) values.push_back(n / x);
    function<int(LL)> getPos = [&](LL x) {
        int pos = (x > r) ? n / x - 1 : values.size() - x;
        /* assert(values[pos] == x); */ return pos;
    };
    // auto primes = getPrimes(r);
    if (!ready) fastSieve(), ready = true;
    vector<LL> cnt(values.size());
    for (auto v : values) cnt[getPos(v)] = v - 1;
    for (auto p : primes)
        for (auto v : values) {
            if ((LL)p * p > v) break;
            cnt[getPos(v)] -= cnt[getPos(v / p)] - cnt[getPos(p - 1)];
        }
    return cnt[getPos(n)];
}  // ends prime counting

// Pollard Rho factorization
void rho(LL n, LL c, vector<LL>& factors) {
    if (n == 1) return;
    if (n < P) {                                // use sieve
        if (!ready) fastSieve(), ready = true;  // !!!
        while (n > 1) {
            int p = minPrime[n];
            while (n > 1 && minPrime[n] == p) factors.push_back(p), n /= p;
        }
        return;
    }
    if (isPrime(n)) {
        factors.push_back(n);
        return;
    }
    if (!(n & 1)) {
        factors.push_back(2);
        rho(n / 2, c, factors);
        return;
    }
    LL x = 2, s = 2, p = 1, l = 1;
    function<LL(LL)> f = [&c, &n](LL x) {
        return (LL)(((__int128)x * x + c) % n);
    };
    while (true) {
        x = f(x);
        LL g = __gcd(abs(x - s), n);
        if (g != 1) {
            rho(g, c + 1, factors);
            rho(n / g, c + 1, factors);
            return;
        }
        if (p == l) s = x, p <<= 1, l = 0;
        l++;
    }
}
vector<pair<LL, int>> factorize(LL n) {
    vector<LL> p;
    rho(n, 1, p);
    sort(p.begin(), p.end());
    vector<pair<LL, int>> f;
    for (int i = 0, j = 0; i < p.size(); i = j) {
        while (j < p.size() && p[i] == p[j]) j++;
        f.emplace_back(p[i], j - i);
    }
    return f;
}  // ends pollar rho factorization

LL phi(LL n) {  // euler phi
    auto f = factorize(n);
    LL r = n;
    for (auto p : f) r -= r / p.first;
    return r;
}

LL primitiveRoot(LL n) {
    if (n <= 0) return -1;
    if (n == 1 || n == 2 || n == 4) return n - 1;
    auto f = factorize(n);
    if (f[0].first == 2 && (f[0].second != 1 || f.size() != 2)) return -1;
    if (f[0].first != 2 && f.size() != 1) return -1;
    int phin = phi(n);
    f = factorize(phin);
    for (int g = 2; g < n; g++) {
        if (power(g, phin, n) != 1) continue;
        bool ok = true;
        for (auto& p : f)
            if (power(g, phin / p.first, n) == 1) {
                ok = false;
                break;
            }
        if (ok) return g;
    }
    assert(false);
    return -1;
}

// a ^ x = b (mod c)
LL discreteLogarithm(LL a, LL b, LL m) {
    if (b == 1) return 0;
    unordered_map<LL, LL> M;
    int c = (int)sqrt((double)m) + 1;
    LL v = power(a, c, m), pv = 1, w = b;
    for (int i = 1; i <= c; i++) pv = mul(pv, v, m), M[pv] = i;
    LL ans = -1;
    for (int i = 0; i < c; i++) {
        if (M.find(w) != M.end()) {
            int lg = M[w] * c - i;
            if (ans == -1 || lg < ans) ans = lg;
        }
        w = mul(w, a, m);
    }
    return ans;
}  // TODO swap steps for finding min

// x ^ n = a (mod m)
vector<LL> discreteRoot(LL a, LL n, LL m) {
    LL g = primitiveRoot(m);
    LL sq = (LL)sqrt(m) + 1;
    // TODO terminar
}

LL legendre(LL n, LL p) { return power(n, (p - 1LL) / 2LL, p); }

// Tonelli Shanks TODO test
LL discreteSqrt(LL a, LL p) {  // returns -1 if no solution
    if (a == 0) return 0;
    if (p == 2) return a;
    if (legendre(a, p) != 1) return -1;
    LL d = p - 1;
    int s = 0;
    while (!(d & 1)) s++, d >>= 1LL;
    LL q = 2;
    while (legendre(q, p) != p - 1) q++;
    LL t = power(a, (d + 1) / 2, p);
    LL r = power(a, d, p);
    while (r != 1) {
        int i = 0;
        LL v = 1;
        while (power(r, v, p) != 1) i++, v *= 2LL;
        LL e = power(2, n - i - 1, p);
        LL u = power(q, k * e, p);
        t = mul(t, u, p);
        r = mul(r, mul(u, u, p), p);
    }
    return t;
}  // if y is a solution; then -y is too

pair<LL, LL> crt(LL r1, LL m1, LL r2, LL m2) {
    // TODO test
    LL d = __gcd(m1, m2);
    if (r1 % d != r2 % d) return {-1, -1};
    LL rd = r1 % d;
    r1 /= d;
    m1 /= d;
    r2 /= d;
    m2 /= d;
    if (m1 < m2) {
        swap(r1, r2);
        swap(m1, m2);
    }
    LL k = (r2 - r1) % m2;
    if (k < 0) k += m2;
    LL x, y;
    egcd(m1, m2, x, y);
    x %= m2;
    if (x < 0) x += m2;
    k *= x;
    k %= m2;
    return {m1 * m2 * d, (k * m1 + r1) * d + rd};
}

// todas las fracciones reducidas tal que
// el denominador es menor o igual que n
// a / b y c / d son consecutivas si y solo si:
// b * c - a * d = 1; TODO test
vector<pair<LL, LL>> farey(int n) {
    LL a = 0, b = 1, c = 1, d = n;
    vector<pair<LL, LL>> s;
    s.push_back({0, 1});
    while (c < n) {
        LL k = (n + b) / d;
        LL e = k * c - a;
        LL f = k * d - b;
        a = c;
        b = d;
        c = e;
        d = f;
        s.push_back({a, b});
    }
    return s;
}

// fibonacci sequence O(log(n))
void _fib(LL n, LL m, LL& x, LL& y) {
    if (n == 1) {
        x = 1;
        y = 1;
    } else if (n & 1) {
        _fib(n - 1, m, y, x);
        y += x;
        if (y >= m) y -= m;
    } else {
        LL a, b;
        _fib(n >> 1, m, a, b);
        y = (mul(a, a, m) + mul(b, b, m)) % m;
        x = (mul(a, b, m) + mul(a, b - a + m, m)) % m;
    }
}

LL fib(LL n, LL m) {  // O(log(n))
    assert(n > 0);
    LL x, y;
    _fib(n, m, x, y);
    return x;
}  // ends fibonacci

vector<LL> partitions(int n, LL m) {  // O(n ^ (3 / 2))
    vector<LL> p(n + 1);
    p[0] = 1;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; i >= (3 * j * j - j) / 2; j++) {
            LL u = p[i - (3 * j * j - j) / 2];
            if (!(j & 1)) u = m - u;
            add(p[i], u, m);
            if (i >= (3 * j * j + j) / 2) {
                LL v = p[i - (3 * j * j + j) / 2];
                if (!(j & 1)) v = m - v;
                add(p[i], v, m);
            }
        }
    }
    return p;
}  // end partitions implementation

// O(n * n * log(n))
template <const int M, const int B = 63>  // B bits
struct linearRec {
    int n;
    vector<int> f;           // f={f(1), f(2), ..., f(n)}
    vector<int> t;           // f(k)=t[0]*f(k-1)+t[1]*f(k-2)+...+t[n-1]*f(k-n)
    vector<vector<int>> fn;  // for bin power

    vector<int> add(vector<int>& a, vector<int>& b) {
        vector<int> r(2 * n + 1);
        for (int i = 0; i <= n; i++)
            for (int j = 0; j <= n; j++) {
                r[i + j] += (long long)a[i] * b[j] % M;
                if (r[i + j] >= M) r[i + j] -= M;
            }
        for (int i = 2 * n; i > n; i--) {
            for (int j = 0; j < n; j++) {
                r[i - 1 - j] += (long long)r[i] * t[j] % M;
                if (r[i - 1 - j] >= M) r[i - 1 - j] -= M;
            }
            r[i] = 0;
        }
        r.erase(r.begin() + n + 1, r.end());
        return r;
    }

    linearRec(vector<int>& _f, vector<int>& _t) : f(_f), t(_t) {
        n = f.size();
        vector<int> a(n + 1);
        a[1] = 1;
        fn.push_back(a);
        for (int i = 1; i < B; i++) fn.push_back(add(fn[i - 1], fn[i - 1]));
    }

    int calc(long long k) {
        vector<int> a(n + 1);
        a[0] = 1;
        for (int i = 0; i < B; i++)
            if (k & (1LL << i)) a = add(a, fn[i]);
        int r = 0;
        for (int i = 0; i < n; i++) {
            r += (long long)a[i + 1] * f[i] % M;
            if (r >= M) r -= M;
        }
        return r;
    }
};  // ends linear recurrence solver

// TODO add latice count

// TODO add sum(phi(n)) for(n <= N), N <= 1e12, project euler

// TODO add Berlekamp-Massey algorithm

// TODO Lucas theorem

// TODO Gauss Matrix Solver Modulo

void init() {
    fastSieve();
    ready = true;
}

void testMillerRabin() {
    assert(isPrime(2));
    assert(!isPrime(4));
    assert(isPrime(3));
    assert(!isPrime(141));
    assert(isPrime(1000 * 1000 * 1000 + 7));
    assert(!isPrime(1000 * 1000 * 1000 + 11));
    cerr << nextPrime(1000 * 1000 * 1000) << "\n";
    long long np = nextPrime(100LL * 1000LL * 1000LL * 1000LL * 1000LL);
    init();
    assert(isPrime(np));
    assert(!isPrime(np + 2));
    cerr << "test ok\n";
}

void testPollardRho() {
    LL n = 45;
    // LL n = 12123423422424LL;
    auto f = factorize(n);
    for (auto x : f) {
        cerr << "(" << x.first << ", " << x.second << ")";
    }
    cerr << "\n";
}

void testIsPrime() {
    int t;
    cin >> t;
    while (t--) {
        long long n;
        cin >> n;
        cout << (isPrime(n) ? "YES" : "NO") << "\n";
    }
}

void testFibonacci() {
    for (int i = 1; i < 20; i++) {
        cerr << fib(i, 1000000007) << " ";
    }
    cerr << "\n";
}

void testPartitions() {
    auto p = partitions(10, 1000000007);
    for (auto x : p) {
        cerr << x << " ";
    }
    cerr << "\n";
}

void testCountPrimes() {
    vector<LL> test = {100, 1000000, 100000000000LL};
    for (auto x : test) {
        cerr << x << " " << countPrimes(x) << "\n";
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);

    // testMillerRabin();
    // testPollardRho();
    // testIsPrime();
    // testFibonacci();
    // testPartitions();
    testCountPrimes();
}
