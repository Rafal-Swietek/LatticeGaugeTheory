
#define CONFIG 7
#include <set>
// #include <cstdint>
// #include <stdexcept>
// #include <armadillo>
// #include <cstdint>
// #include <stdexcept>
// #include <string>
// #include <iostream>
// #include <algorithm>

#include "../../../QHamSolver/include/user_interface.hpp"
#include "../includes/PackedLattice2D/LatticeRepresentation.hpp"
// #include "includes/symmetry_generators2D.hpp"
#include "../includes/hilbert_space.hpp"

int main(const int argc, char* argv[]) 
{
    int Lx = 10;
    int Ly = 10;
    int Nu = 98;
    int Nd = 2;

    PackedLattice2D state(Lx, Ly);
    state.fill_all();
    // state.print_state();
    std::cout << std::endl;

    // QHS::Symmetric2DLattice<QOps::particle::fermion> _hilbert(Lx, Ly, Nu, Nd, v_1d<QOps::generic_operator<cpx, PackedLattice2D>>(), 0 );
    // auto mapy = _hilbert.get_mapping();
    // for(auto& x : mapy){
    //     x.print_state();
    //     std::cout << std::endl;
    // }
    // std::cout << mapy.size() << std::endl;
    QHS::Symmetric2DLattice<QOps::particle::fermion> _hilbert;
    _hilbert.set_dimensions(Lx, Ly, Nu, Nd, 0);
    std::cout << _hilbert.count_dimension() << std::endl;
    // PackedLattice2D state(Lx, Ly);
    // state.set_block(0, 0, 1);
    // state.set_block(0, 1, 1);
    // state.set_block(0, 5, 1);
    // state.set_block(1, 1, 1);
    // state.set_block(2, 2, 1);
    // state.set_block(4, 0, 1);
    // state.set_block(4, 3, 1);
    // state.set_block(4, 5, 1);
    // // state.set_block(5, 0, 1);
    // // state.set_block(5, 1, 1);
    // // state.set_block(5, 2, 1);
    // // state.set_block(1, 1, 2);
    // state.print_state();
    // std::cout << std::endl << std::endl;


    // auto Tx = QOps::_translation_X_symmetry(Lx, Ly, 0, 1);
    // auto Ty = QOps::_translation_Y_symmetry(Lx, Ly, 0, 1);

    // auto [new_state, tmp1] = Tx(state);
    // new_state.print_state();
    // std::cout << std::endl << std::endl;

    // auto [new_state2, tmp2] = Ty(state);
    // new_state2.print_state();
    // std::cout << std::endl << std::endl;
    
    return 0;
}