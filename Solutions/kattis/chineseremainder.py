# https://open.kattis.com/problems/chineseremainder

def e_gcd(a, b):
    if b == 0:
        return [a, 1, 0]
    g, x, y = e_gcd(b, a % b)
    return [g, y, x - y * (a // b)]

def crt(r1, m1, r2, m2):
    g, x, y = e_gcd(m1, m2)
    m = m1 * m2
    r = (((r1 * m2 * y + r2 * m1 * x) % m) + m) % m
    return [r, m]

t = int(input())

while t > 0:
    t = t - 1
    r1, m1, r2, m2 = map(int, input().split())
    r, m = crt(r1, m1, r2, m2)
    print("%d %d" % (r, m))
