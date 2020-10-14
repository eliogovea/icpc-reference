template <class F, int N = 1 << 20>
double simpson(F f, double l, double r) {
    double h = (r - l) / (double)N, s = 0.0;
    for (int i = 0; i <= N; i++) {
        double x = l + h * i;
        s += f(x) * ((i == 0 || i == N) ? 1 : ((i & 1) == 0) ? 2 : 4);
    }
    s *= h / 3.0;
    return s;
}
