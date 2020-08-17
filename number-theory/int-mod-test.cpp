#include "int-mod.hpp"

int main() {
    constexpr int modulo = 1000 * 1000 * 1000 + 7;
    using int_mod = int32_mod_t<modulo>;

    int_mod x{2}, y{1};
    std::cout << x.value << " " << y.value << "\n";

    auto z = y / x;

    std::cout << z << "\n";
    std::cout << int_mod::power(x, modulo - 2) << "\n";
}
