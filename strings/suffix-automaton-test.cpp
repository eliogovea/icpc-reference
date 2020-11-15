#include <iostream>
#include <string>

#include "suffix-automaton.hpp"

// https://judge.yosupo.jp/problem/number_of_substrings
void NumberOfSubstrings() {
    using SuffixAutomaton::SuffixAutomaton;
    constexpr std::size_t MaxWordLength = 500 * 1000 + 10;
    constexpr std::size_t MaxAutomatonSize = 2 * MaxWordLength;
    static SuffixAutomaton<char, MaxAutomatonSize> suffixAutomaton;
    std::string word;
    std::cin >> word;
    for (const auto& c : word) { suffixAutomaton.Extend(c); }
    std::int64_t answer = 0;
    for (const auto& state : suffixAutomaton) {
        if (state.suffixLink != nullptr) {
            answer += state.maxLength - state.suffixLink->maxLength;
        }
    }
    std::cout << answer << "\n";
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    NumberOfSubstrings();
}