#include <bits/stdc++.h>

using namespace std;

template <class T>
vector<int> manacher(const T& s, int n) {
    vector<int> len(2 * n - 1);
    for (int i = 0, l = 0, r = -1; i < 2 * n - 1; i++) {
        int x = (i + 1) >> 1, y = i >> 1, z = 0;
        if (x < r) {
            z = min(len[2 * (l + r) - i], r - x);
        }
        while (0 <= x - z - 1 && y + z + 1 < n &&
               s[x - z - 1] == s[y + z + 1]) {
            z++;
        }
        len[i] = z;
        if (y + z >= r) {
            l = x - z;
            r = y + z;
        }
    }
    return len;
}

bool isPalindrome(int l, int r, const vector<int>& len) {
    int x = (l + r + 1) >> 1, y = (l + r) >> 1, z = len[l + r];
    return x - z <= l && r <= y + z;
}

// accepted
void timus_2060() {
    string s;
    cin >> s;

    int n = s.size();

    auto r = manacher(s, n);

    vector<int> start(n);
    vector<int> end(n);

    for (int i = 0; i < 2 * n - 1; i++) {
        int x = (i + 1) >> 1;
        int y = i >> 1;
        int z = r[i];

        if (x - z <= y + z) {
            start[x - z]++;
            if (x != y) {
                start[x]--;
            } else {
                if (x + 1 < n) {
                    start[x + 1]--;
                }
            }

            end[x]++;
            if (y + z + 1 < n) {
                end[y + z + 1]--;
            }
        }
    }

    for (int i = 1; i < n; i++) {
        start[i] += start[i - 1];
        end[i] += end[i - 1];
    }

    long long answer = 0;

    for (int i = 0; i + 1 < n; i++) {
        answer += (long long)end[i] * start[i + 1];
    }

    cout << answer << "\n";
}

int main() { timus_2060(); }
