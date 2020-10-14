/*
1 -  Utiles (add, mul, power)
2 -  EGCD
                 2.1 - EGCD
                 2.2 - CRT
                 2.3 - Modular Inverse
3 -  Fast Sieve O(n) [necesario para otras funciones]
4 -  Primes
                 4.1 - Miller Rabin
                 4.2 - Next Prime
                 4.3 - Count Primes O(n ^ 0.75)
5 -  Pollard Rho [usa Fast Sieve + Miller Rabin]
6 -  Euler phi
7 -  Primitive Root
8 -  Discrete Logarithm
9 -  Discrete Root
10 - Discrete sqrt
11 - Farey
12 - Fibonacci O(log(n))
13 - Partitions O(n ^ 1.5)
14 - Linear recurrence [calc k-th in O(n*n*log(k))]
15 - Latice points under al line
*/

/// 1 - Utiles
typedef long long LL;
inline void add(LL& a, LL b, LL M) {
    a += b;
    if (a >= M) a -= M;
}

inline LL mul(LL a, LL b, LL m) {
    // return (__int128)a * b % m; // !!!!
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
// 2 - Extended GCD
// 2.1 - Extended GCD[solve a * x + b * y = gcd(a, b)]
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

// 2.1 - Chinese Remainder Theorem
// solve a system of congruence equations
pair<LL, LL> crt(LL r1, LL m1, LL r2, LL m2) {
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

// 2.2 - Modular Inverse
// inverse(n, m) = x iff m * x % m == 1
LL inverse(LL n, LL m) {
    LL x, y;
    LL g = egcd(n, m, x, y);
    if (g != 1) return -1;
    x %= m;
    if (x < 0) x += m;
    assert((n * x % m) == 1);
    return x;
}

// 3 - Fast Sieve
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
        for (auto p : primes) {
            if (p * i >= P) break;
            _isPrime[p * i] = false;
            minPrime[p * i] = p;
            if (i % p == 0) break;
        }
    }
}
// 4 - Primality
// 4.1 - Miller Rabin primality test (OK)
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
    // if n int:  test = {2, 7, 61}
    // if n<3e18: test = {2, 3, 5, 7, 11, 13, 17, 19, 23}
    static vector<LL> test = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    for (long long x : test) {
        if (x % n == 0) return true;
        if (witness(x, n, s, d)) return false;
    }
    return true;
}  // ends Miller Rabin primality test !!!

// 4.2 - Next Prime
LL nextPrime(LL x) {
    for (x = max(2LL, x + 1); x % 6; x++)
        if (isPrime(x)) return x;
    for (; true; x += 6) {
        if (isPrime(x + 1)) return x + 1;
        if (isPrime(x + 5)) return x + 5;
    }
}
// 4.3 - Count Primes O(n ^ 0.75)
// O(n^0.75) same idea for other functions over primes
LL countPrimes(LL n) {
    int r = sqrt((double)n);
    while ((LL)(r + 1) * (r + 1) <= n) r++;
    vector<LL> values;  // all [n / x]
    for (LL x = 1; x <= n; x = n / (n / x) + 1) values.push_back(n / x);
    function<int(LL)> pos = [&](LL x) {
        if (x > r) return (int)(n / x) - 1;
        return (int)(values.size() - x);
    };
    // auto primes = getPrimes(r);
    if (!ready) fastSieve(), ready = true;
    vector<LL> cnt(values.size());
    for (auto v : values) cnt[pos(v)] = v - 1;
    for (auto p : primes)
        for (auto v : values) {
            if ((LL)p * p > v) break;
            cnt[pos(v)] -= cnt[pos(v / p)];
            cnt[pos(v)] += cnt[pos(p - 1)];
        }
    return cnt[pos(n)];
}  // ends prime counting

// 5 - Factorization
void rho(LL n, LL c, vector<LL>& fp) {
    if (n == 1) return;
    if (n < P) {                                // use sieve
        if (!ready) fastSieve(), ready = true;  // !!!
        while (n > 1) {
            int p = minPrime[n];
            while (n > 1 && minPrime[n] == p) fp.push_back(p), n /= p;
        }
        return;
    }
    if (isPrime(n)) {
        fp.push_back(n);
        return;
    }
    if (!(n & 1)) {
        fp.push_back(2);
        rho(n / 2, c, fp);
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
            rho(g, c + 1, fp);
            rho(n / g, c + 1, fp);
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

// 6 - Euler Phi (counts the nomber of relative primes to n in the interval [1,
// n])
LL phi(LL n) {  // euler phi
    auto f = factorize(n);
    LL r = n;
    for (auto p : f) r -= r / p.first;
    return r;
}

// 7 - Primitive root, ()
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

// 8 - Discrete logarithm, baby step, giant step algorithm
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
}  // TODO swap steps, calc min answer faster

// 9 - Discrete Root, x ^ n = a (mod m)
// Same idea used in discrete logarithm using
// primitive root (g ^ n) ^ x') = g ^ (a')
// + congruence equation n * x' = a' (modulo phi(m))
// gcd(n, phi(m)) solutions

// 10 - Discrete sqrt
LL legendre(LL n, LL p) { return power(n, (p - 1LL) / 2LL, p); }

// Tonelli Shanks
LL discreteSqrt(LL a, LL p) {  // -1 if no solution
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
        LL e = power(2, s - i - 1, p);
        LL u = power(q, d * e, p);
        t = mul(t, u, p);
        r = mul(r, mul(u, u, p), p);
    }
    return t;
}  // if y is a solution; then -y is too

// 11 - Farey
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

// 12 - Fibonacci sequence O(log(n))
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

// 13 - Partitions
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

// 14 - Linear recurrence
// O(n * n * log(n))
template <const int M, const int B = 63>  // B bits
struct linearRec {
    int n;
    vector<int> f;  // f={f(1), f(2), ..., f(n)}
    // f(k)=t[0]*f(k-1)+t[1]*f(k-2)+...+t[n-1]*f(k-n)
    vector<int> t;
    vector<vector<int>> fn;  // for bin power
    vector<int> add(vector<int>& a, vector<int>& b) {
        vector<int> r(2 * n + 1);
        for (int i = 0; i <= n; i++)
            for (int j = 0; j <= n; j++) add(r[i + j], mul(a[i], b[j], M), M);
        for (int i = 2 * n; i > n; i--) {
            for (int j = 0; j < n; j++)
                add(r[i - 1 - j], mul(r[i], t[j], M), M);
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

// 15 - Lattice points below segment :
// solves for sigma [(a + b * i) / m] where 0 <= i < n.
LL solve(LL n, LL a, LL b, LL m) {
    if (b == 0) return n * (a / m);
    if (a >= m) return n * (a / m) + solve(n, a % m, b, m);
    if (b >= m) return (n - 1) * n / 2 * (b / m) + solve(n, a, b % m, m);
    return solve((a + b * n) / m, (a + b * n) % m, m, b);
}
