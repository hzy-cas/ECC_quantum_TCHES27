#include "../cpp/GateManager.h"
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    GateManager gm(argc > 1 ? std::stoull(argv[1]) : 10000000);
    unsigned width, count;
    while (std::cin >> width >> count) {
        gm.clear();
        for (unsigned i = 0; i < count; ++i) {
            char op; unsigned a, b, c;
            std::cin >> op >> a >> b >> c;
            if (op == 'X') gm.add_X(a);
            else if (op == 'C') gm.add_CNOT(a, b);
            else if (op == 'T') gm.add_Toffoli(a, b, c);
            else return 2;
        }
        auto [fd, cd, td, tc, cc] = gm.optimize_and_count(width);
        std::cout << fd << ' ' << cd << ' ' << td << ' ' << tc << ' ' << cc << '\n';
    }
}
