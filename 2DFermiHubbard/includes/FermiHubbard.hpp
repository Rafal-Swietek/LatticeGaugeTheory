#pragma once

//#define _XYZ_SYM

#ifndef _FermiHubbard
#define _FermiHubbard

#include "../../../QHamSolver/include/lattices/_base.hpp"
#include "../../../QHamSolver/include/hilbert_space/symmetries.hpp"
#include "../../../QHamSolver/include/hilbert_space/double_u1.hpp"
// #include "helpers.hpp"

using U1U1Hilbert = QHS::U1U1_hilbert_space<QHS::U1::charge, true>;

#ifdef USE_REAL_SECTORS
    using elem_ty = double;
#else
    using elem_ty = cpx;
#endif

const arma::vec _fermion_number = {0, 1, 1, 2};

// using elem_ty = cpx;
/// @brief Fully anisotropic spin chain (XYZ) with point symmetries
class FermiHubbard : 
#ifdef NO_POINT_SYM
    public QHS::hamiltonian_base<elem_ty, U1U1Hilbert>
#else
    public QHS::hamiltonian_base<elem_ty, QHS::point_symmetric<QOps::particle::fermion>>
#endif
{
    //<! ----------------------------------------------------- INHERIT TYPEDEFs FROM BASE
    #ifdef NO_POINT_SYM
    typedef typename QHS::hamiltonian_base<elem_ty, U1U1Hilbert>::matrix        matrix;
    typedef typename QHS::hamiltonian_base<elem_ty, U1U1Hilbert>::sparse_matrix sparse_matrix;
    #else
    typedef typename QHS::hamiltonian_base<elem_ty, QHS::point_symmetric<QOps::particle::fermion>>::matrix        matrix;
    typedef typename QHS::hamiltonian_base<elem_ty, QHS::point_symmetric<QOps::particle::fermion>>::sparse_matrix sparse_matrix;
    #endif

    //<! ----------------------------------------------------- MODEL PARAMETERS
private:
    lattice::lattice2D _lattice;             // lattice of system (cubic by default)
    v_1d<QOps::genOp> symmetry_generators;    // list of symmetry generators

    //<! Symmetry contained in struct
    struct {
        int kx_sym;
        int px_sym;
        int ky_sym;
        int py_sym;
        int Nup;
        int Ndown;
    } syms;

    double _J = 1;                          // hopping strength
    double _delta = 0.0;                    // staggered field
    double _U = 0.0;                        // interaction
    double _F = 0.0;                       // tilted field
    // double _Fy = 0.0;                       // tilted field
    int Lx = 4;
    int Ly = 4;
    bool _use_symmetries = true;            // [temporary] boolean choose if use symmetries

    //<! ----------------------------------------------------- INITIALIZE MODEL
    virtual void init() override;
    void set_symmetry_generators();

public:
    //<! ----------------------------------------------------- CONSTRUCTORS
    FermiHubbard() = default;
    FermiHubbard(std::istream& os);
    FermiHubbard(int _BC, unsigned int Lx, unsigned int Ly, double J, double U, double delta, double F, int Nup, int Ndown, bool use_syms = true);

    //<! ----------------------------------------------------- HAMILTONIAN BUILDERS
    virtual void create_hamiltonian() override;
    virtual sparse_matrix create_local_hamiltonian(int site) override;
    virtual void set_hamiltonian_elements(u64 k, elem_ty value, u64 new_idx) override;

    //<! ----------------------------------------------------- OVERRIDEN OPERATORS
    virtual std::ostream& write(std::ostream&) const override;
    virtual std::istream& read(std::istream&) override;

    //<! ----------------------------------------------------- ADDITIONAL METHODS
    void print2D_state(u64 state) {
        auto check_spin = QOps::__builtins::get_digit(this->system_size);

        arma::mat state_in_lattice(this->Ly, this->Lx, arma::fill::zeros);
        for(long int j = 0; j < this->system_size; j++)
        {
            int n_particle = check_spin(state, j);
            // if(n_particle > 1) n_particle--; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
            int xcoord = j % this->Lx;
            int ycoord = this->Ly - 1 - j / this->Lx;
            // printSeparated(std::cout, "\t", 20, true, state, j, xcoord, ycoord, n_particle);
            state_in_lattice(ycoord, xcoord) = n_particle;
        }
        std::cout << state << "\t" << boost::dynamic_bitset<>(2*this->system_size, state) << "\n" << state_in_lattice << std::endl << std::endl;
    }
    
};


#endif