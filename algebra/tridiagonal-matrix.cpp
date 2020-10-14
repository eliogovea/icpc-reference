#include <bits/stdc++.h>

using namespace std;

// a[i] * x[i - 1] + b[i] * x[i] + c[i] * x[i + 1] = d[i]
// a[0] = 0, c[n - 1] = 0
template <class T>
vector<T> tridiagonal(vector<T> a, vector<T> b, vector<T> c, vector<T> d) {
    int n = d.size();
    c[0] /= b[0];
    for (int i = 1; i < n; i++) c[i] /= (b[i] - a[i] * c[i - 1]);
    d[0] /= b[0];
    for (int i = 1; i < n; i++)
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
    vector<T> x(n);
    x[n - 1] = d[n - 1];
    for (int i = n - 2; i >= 0; i--) x[i] = d[i] - c[i] * x[i + 1];
    return x;
}

int main() {}
