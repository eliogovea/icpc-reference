// https://www.hackerrank.com/challenges/largest-rectangle/problem

#include <bits/stdc++.h>

using namespace std;

long long largestRectangle(vector<int> h) {
    int n = h.size();
    vector<int> to_left(n, -1);
    {
        vector<int> s(n);
        int top = 0;
        for (int i = n - 1; i >= 0; i--) {
            while (top > 0 && h[s[top - 1]] > h[i]) {
                to_left[s[top - 1]] = i;
                top--;
            }
            s[top++] = i;
        }
    }
    vector<int> to_right(n, n);
    {
        vector<int> s(n);
        int top = 0;
        for (int i = 0; i < n; i++) {
            while (top > 0 && h[s[top - 1]] > h[i]) {
                to_right[s[top - 1]] = i;
                top--;
            }
            s[top++] = i;
        }
    }
    long long answer = 0;
    for (int i = 0; i < n; i++) {
        answer = max(answer, (long long)h[i] *
                                 (long long)(to_right[i] - to_left[i] - 1LL));
    }
    return answer;
}

int main() {
    int n;
    cin >> n;
    vector<int> h(n);
    for (int h_i = 0; h_i < n; h_i++) {
        cin >> h[h_i];
    }
    long result = largestRectangle(h);
    cout << result << endl;
    return 0;
}
