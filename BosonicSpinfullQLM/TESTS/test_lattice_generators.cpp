
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
    int Lx = 16;
    int Ly = 16;
    int Nu = Lx*Ly-2;
    int Nd = 2;
    int V = Lx * Ly;
    // PackedLattice2D state(Lx, Ly);
    // state.fill_all();
    // // state.print_state();
    // std::cout << std::endl;
    omp_set_num_threads(10);
    // // std::cout << dim_tot << std::endl;
    // std::vector<int> Z2sym = {-1, 1};
    // std::vector<int> noZ2sym = {1};
    // u64 dim_sum = 0;
    // printSeparated(std::cout, "\t", 16, true, "kx", "ky", "px", "py", "pd", "z", "D");
    // for(int kx = 0; kx < Lx; kx++){
    // for(int ky = 0; ky < Ly; ky++){
    //     bool is_piX_sector = Lx % 2? false : kx  == int(Lx / 2);
    //     bool is_piY_sector = Ly % 2? false : ky  == int(Ly / 2);
    //     bool is_kx_real_sector = kx==0 || is_piX_sector;
    //     bool is_ky_real_sector = ky==0 || is_piY_sector;

    //     std::vector<int> px_sec = is_kx_real_sector? Z2sym : noZ2sym;
    //     std::vector<int> py_sec = is_ky_real_sector? Z2sym : noZ2sym;
    //     std::vector<int> z_sec = Nu == Nd? Z2sym : noZ2sym;
    //     // if(is_kx_real_sector && is_ky_real_sector){
    //     //     rot_sec = {-1, 1};
    //     // }
    //     for(int z : z_sec){
    //         for(int px : px_sec){ for(int py : py_sec){     
    //             std::vector<int> pd_sec = (Lx == Ly && kx == ky && px == py)? Z2sym : noZ2sym;
    //             for(int pd : pd_sec){
    //                     v_1d<QOps::generic_operator<cpx, PackedLattice2D>> symmetry_generators;
    //                     if(is_kx_real_sector)
    //                         symmetry_generators.emplace_back(QOps::_mirror_X_symmetry(Lx, Ly, px));
    //                     if(is_ky_real_sector)
    //                         symmetry_generators.emplace_back(QOps::_mirror_Y_symmetry(Lx, Ly, py));

    //                     //<! Symmetries <of square lattice (diagonal mirror symmetries)
    //                     //<! Add mirror >symmetry along diagonal valid for all kx = ky
    //                     if(Lx == Ly && px == py && kx == ky)
    //                         symmetry_generators.emplace_back(QOps::_mirror_diagonal_symmetry(Lx, Ly, pd));
                        
    //                     if(Nu == Nd)
    //                         symmetry_generators.emplace_back(QOps::_spin_inversion_symmetry(Lx, Ly, z));
    //                     QHS::Symmetric2DLattice<QOps::particle::fermion> _hilbert(Lx, Ly, Nu, Nd, symmetry_generators, 0, kx, ky, 1, 1, true, true);
    //                     auto mapy = _hilbert.get_mapping();
    //                     dim_sum += mapy.size();
    //                     printSeparated(std::cout, "\t", 16, true, kx, ky, px, py, pd, z, mapy.size());
    //                 }
    //             }
    //         }
    //     }
    // }      
    // }
    // std::cout << std::endl << std::endl;
    // std::cout << "Summed Dimension: " << dim_sum << std::endl;

    // QHS::Symmetric2DLattice<QOps::particle::fermion> _hilbert(Lx, Ly, Nu, Nd, v_1d<QOps::generic_operator<cpx, PackedLattice2D>>(), 0, 0, 0, 1, 1, true, true);
    // auto mapy = _hilbert.get_mapping();
    // for(auto& x : mapy){
    //     x.print_state();
    //     std::cout << std::endl;
    // }
    // std::cout << mapy.size() << std::endl;

    QHS::Symmetric2DLattice<QOps::particle::fermion> _hilbert;
    _hilbert.set_dimensions(Lx, Ly, Nu, Nd, 0);
    u64 dim_tot = _hilbert.count_dimension();
    std::cout << "Total Dimension: " << dim_tot << std::endl;
    // for(int Nd = 1; Nd < 4; Nd++){
    //     for(int Nu = 1; Nu < 4; Nu++)
    //     {
    //         QHS::Symmetric2DLattice<QOps::particle::fermion> _hilbert(Lx, Ly, Nu, Nd, v_1d<QOps::generic_operator<cpx, PackedLattice2D>>(), 0, 0, 0, 1, 1, true, true);
    //         auto mapy = _hilbert.get_mapping();
    //         std::cout << "N\u2193 = " << Nd << "\t\t" << "N\u2191 = " << Nu << "\t\t D = " <<  mapy.size() << std::endl;
        
    //     }
    // }
    // PackedLattice2D state(Lx, Ly);
    // state.set_block(0, 0, DOWN_MONOMER);
    // state.set_block(1, 0, DOWN_MONOMER);
    // state.set_block(2, 0, DIMER_RIGHT);
    // state.set_block(3, 0, DIMER_LEFT);
    // state.set_block(0, 1, DOWN_MONOMER);
    // state.set_block(1, 1, UP_MONOMER);
    // state.set_block(2, 1, DIMER_UP);
    // state.set_block(3, 1, UP_MONOMER);
    // state.set_block(0, 2, DIMER_UP);
    // state.set_block(1, 2, UP_MONOMER);
    // state.set_block(2, 2, DIMER_DOWN);
    // state.set_block(3, 2, UP_MONOMER);
    // state.set_block(0, 3, DIMER_DOWN);
    // state.set_block(1, 3, UP_MONOMER);
    // state.set_block(2, 3, DOWN_MONOMER);
    // state.set_block(3, 3, UP_MONOMER);
    // state.print_state();
    // std::cout << std::endl << std::endl;


    // auto Tx = QOps::_translation_X_symmetry(Lx, Ly, 0, 1);
    // auto Ty = QOps::_translation_Y_symmetry(Lx, Ly, 0, 1);
    // auto Px = QOps::_mirror_X_symmetry(Lx, Ly, 1);
    // auto Py = QOps::_mirror_Y_symmetry(Lx, Ly, 1);
    // auto D = QOps::_mirror_diagonal_symmetry(Lx, Ly, 1);
    // auto AD = QOps::_mirror_anti_diagonal_symmetry(Lx, Ly, 1);
    // auto C4 = QOps::_rotation_90_symmetry(Lx, Ly, 0);

    // PackedLattice2D new_state(Lx, Ly);
    // cpx tmp;

    // std::cout << "------------------------ Tx * Px" << std::endl;
    // std::tie(new_state, tmp) = Px(state);
    // std::tie(new_state, tmp) = Tx(new_state);
    // new_state.print_state();
    // std::cout << std::endl << std::endl;

    // std::cout << "------------------------ Px * Tx" << std::endl;
    // std::tie(new_state, tmp) = Tx(state);
    // std::tie(new_state, tmp) = Px(new_state);
    // new_state.print_state();
    // std::cout << std::endl << std::endl;
    
    
    return 0;
}