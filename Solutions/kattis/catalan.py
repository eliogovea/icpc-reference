# https://open.kattis.com/problems/catalan

N = 10 * 1000

fact = [0] * (N + 10)
fact[0] = 1
for i in range(1, N + 2):
    fact[i] = fact[i - 1] * i

def catalan(n):
    return fact[2 * n] // fact[n] // fact[n + 1]

t = int(input())

while t > 0:
    t = t - 1
    n = int(input())
    print(catalan(n))
