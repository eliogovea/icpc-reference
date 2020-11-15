#pragma once

#include <array>
#include <map>

namespace SuffixAutomaton {

template <typename CharT, std::size_t Maxsize>
class SuffixAutomaton {
   public:
    SuffixAutomaton() : size_ {0}, root_ {nullptr}, last_ {nullptr} {
        root_ = CreateNewState(0);
        last_ = root_;
    }

    template <typename Iterator>
    explicit SuffixAutomaton(Iterator begin, Iterator end) : SuffixAutomaton() {
        static_assert(
          std::is_same<
            CharT, typename std::iterator_traits<Iterator>::value_type>::value);
        for (auto i = begin; i != end; ++i) { Extend(static_cast<CharT>(*i)); }
    }

    template <typename Container>
    explicit SuffixAutomaton(const Container& c)
        : SuffixAutomaton(std::begin(c), std::end(c)) {}

    struct State {
        State(std::size_t maxLength = 0)
            : maxLength {maxLength}, go {}, suffixLink {nullptr} {}
        std::size_t maxLength;
        std::map<CharT, State*> go;
        State* suffixLink;
    };

    void Extend(const CharT& c) {
        if (last_->go.count(c) != 0) {
            auto to = last_->go[c];
            if (last_->maxLength + 1 == to->maxLength) {
                last_ = to;
            } else {
                auto cloned = CreateClonedState(to, last_->maxLength + 1);
                to->suffixLink = cloned;
                while (last_ != nullptr && last_->go[c] == to) {
                    last_->go[c] = cloned;
                    last_ = last_->suffixLink;
                }
                last_ = cloned;
            }
        } else {
            auto newState = CreateNewState(last_->maxLength + 1);
            while (last_ != nullptr && last_->go.count(c) == 0) {
                last_->go[c] = newState;
                last_ = last_->suffixLink;
            }
            if (last_ == nullptr) {
                newState->suffixLink = root_;
            } else {
                auto to = last_->go[c];
                if (last_->maxLength + 1 == to->maxLength) {
                    newState->suffixLink = to;
                } else {
                    auto cloned = CreateClonedState(to, last_->maxLength + 1);
                    to->suffixLink = cloned;
                    newState->suffixLink = cloned;
                    while (last_ != nullptr && last_->go[c] == to) {
                        last_->go[c] = cloned;
                        last_ = last_->suffixLink;
                    }
                }
            }
            last_ = newState;
        }
    }

    template<typename Iterator>
    bool Check(Iterator begin, Iterator end) {
        State* now = root_;
        for (auto i = begin; i != end; ++i) {
            const auto& it = now->go.find(*i);
            if (it == now->go.end()) {
                return false;
            }
            now = it->second;
        }
        return true;
    }

    std::size_t Size() const { return size_; }
    State* GetRoot() const { return root_; }
    State* GetLast() const { return last_; }

    State const* StateAt(const std::size_t& index) const {
        return states_[index];
    }

    std::size_t IndexOf(const State* s) const {
        return std::distance(std::begin(states_), s);
    }

    void ResetLast() { last_ = root_; }

    State* begin() { return root_; }
    State const* begin() const { return root_; }
    State* end() { return &states_[size_]; }
    State const* end() const { return &states_[size_]; }

   private:
    State* CreateNewState(int maxLength) {
        states_[size_] = State(maxLength);
        return &states_[size_++];
    }
    State* CreateClonedState(State* original, std::size_t maxLength) {
        states_[size_] = *original;
        states_[size_].maxLength = maxLength;
        return &states_[size_++];
    }
    State* root_;
    State* last_;
    std::size_t size_;
    std::array<State, Maxsize> states_;
};

}  // namespace SuffixAutomaton