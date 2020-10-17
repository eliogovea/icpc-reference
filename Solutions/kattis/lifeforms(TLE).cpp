#include <bits/stdc++.h>

using namespace std;

const int MAXS = 2 * 100 * 1006 + 1;

int length[MAXS];
map<char, int> to[MAXS];
int suffixLink[MAXS];

int size;
int last;

int getNew(int _length) {
    length[size] = _length;
    to[size] = map<char, int>();
    suffixLink[size] = -1;
    return size++;
}

int getClone(int from, int _length) {
    length[size] = _length;
    to[size] = to[from];
    suffixLink[size] = suffixLink[from];
    return size++;
}

void init() {
    size = 0;
    last = getNew(0);
}

void add(char c) {
    int p = last;
    if (to[p].find(c) != to[p].end()) {
        int q = to[p][c];
        if (length[p] + 1 == length[q]) {
            last = q;
        } else {
            int clone = getClone(q, length[p] + 1);
            suffixLink[q] = clone;
            while (p != -1 && to[p][c] == q) {
                to[p][c] = clone;
                p = suffixLink[p];
            }
            last = clone;
        }
    } else {
        int cur = getNew(length[p] + 1);
        while (p != -1 && to[p].find(c) == to[p].end()) {
            to[p][c] = cur;
            p = suffixLink[p];
        }
        if (p == -1) {
            suffixLink[cur] = 0;
        } else {
            int q = to[p][c];
            if (length[p] + 1 == length[q]) {
                suffixLink[cur] = q;
            } else {
                int clone = getClone(q, length[p] + 1);
                suffixLink[q] = clone;
                suffixLink[cur] = clone;
                while (p != -1 && to[p][c] == q) { to[p][c] = clone; }
            }
        }
        last = cur;
    }
}

int n;
string s;
bool used[105][MAXS];
bool ok[MAXS];

int freq[MAXS];
int order[MAXS];

int from[MAXS];
char let[MAXS];
vector<string> ans;

void solve() {
    bool firstCase = true;
    while (cin >> n && n) {
        if (!firstCase) { cout << "\n"; }
        firstCase = false;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < size; j++) { used[i][j] = false; }
            ok[i] = false;
        }
        init();
        for (int i = 0; i < n; i++) {
            cin >> s;
            for (int j = 0; j < s.size(); j++) {
                add(s[j]);
                used[i][last] = true;
            }
            last = 0;
        }
        int maxLength = 0;
        for (int i = 0; i < size; i++) {
            maxLength = max(maxLength, length[i]);
        }
        for (int i = 0; i <= maxLength; i++) { freq[i] = 0; }
        for (int i = 0; i < size; i++) { freq[length[i]]++; }
        for (int i = 1; i <= maxLength; i++) { freq[i] += freq[i - 1]; }
        for (int i = 0; i < size; i++) { order[--freq[length[i]]] = i; }
        for (int i = size - 1; i > 0; i--) {
            int x = order[i];
            for (int j = 0; j < n; j++) {
                used[j][suffixLink[x]] |= used[j][x];
            }
        }
        int ansLength = 0;
        for (int i = size - 1; i > 0; i--) {
            int x = order[i];
            if (length[x] < ansLength) { break; }
            int cnt = 0;
            for (int j = 0; j < n; j++) {
                if (used[j][x]) { cnt++; }
            }
            if (cnt > n / 2) {
                ok[x] = true;
                ansLength = length[x];
            }
        }
        if (ansLength == 0) {
            cout << "?\n";
            continue;
        }
        for (int i = 0; i < size; i++) {
            int x = order[i];
            for (map<char, int>::iterator it = to[x].begin(); it != to[x].end();
                 it++) {
                if (length[x] + 1 == length[it->second]) {
                    from[it->second] = x;
                    let[it->second] = it->first;
                }
            }
        }
        ans.clear();
        for (int s = 1; s < size; s++) {
            if (length[s] == ansLength && ok[s]) {
                string tmp;
                int now = s;
                while (now != 0) {
                    tmp += let[now];
                    now = from[now];
                }
                reverse(tmp.begin(), tmp.end());
                ans.push_back(tmp);
            }
        }
        sort(ans.begin(), ans.end());
        for (int i = 0; i < ans.size(); i++) { cout << ans[i] << "\n"; }
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    solve();
}
