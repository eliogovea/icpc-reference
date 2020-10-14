template <class T>
void xor_fft(vector<T>& P, bool invert) {
    int n = P.size();
    int ln = 0;
    while ((1 << ln) < n) ln++;
    for (int bit = 0; bit < ln; bit++)
        for (int mask = 0; mask < n; mask++)
            if (!(mask & (1 << bit))) {
                int u = P[mask], v = P[mask | (1 << bit)];
                P[mask] = u + v;
                P[mask | (1 << bit)] = u - v;
            }
    if (invert)
        for (auto& x : P) x /= T(n);
}
