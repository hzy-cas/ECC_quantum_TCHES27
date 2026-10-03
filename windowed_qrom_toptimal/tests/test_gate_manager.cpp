#include "../backend/src/Infrastructure/GateManager.h"
#include <stdexcept>
#include <vector>
int main() {
    for (int input=0;input<8;++input) {
        GateManager gm;
        gm.add_CNOT(0,3); gm.add_CNOT(0,3); // cancellation
        gm.add_Toffoli(0,1,2); gm.add_X(1); gm.add_CNOT(2,3);
        auto counts=gm.optimize_and_count(4);
        if (std::get<3>(counts)!=1 || std::get<4>(counts)!=1) throw std::runtime_error("counts");
        std::vector<uint8_t> state={static_cast<uint8_t>(input&1),static_cast<uint8_t>((input>>1)&1),
            static_cast<uint8_t>((input>>2)&1),0};
        auto initial=state;
        gm.apply_to_basis(state);
        if (state[2]!=(initial[2]^(initial[0]&initial[1])) || state[3]!=state[2])
            throw std::runtime_error("simulation");
        gm.apply_to_basis(state,true);
        if (state!=initial) throw std::runtime_error("inverse");
    }
}
