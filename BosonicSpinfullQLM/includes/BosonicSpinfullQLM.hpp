#pragma once

//#define _XYZ_SYM

#ifndef _BosonicSpinfullQLM
#define _BosonicSpinfullQLM

// #include "../../../QHamSolver/include/lattices/_base.hpp"
// #include "../../../QHamSolver/include/hilbert_space/symmetries.hpp"
// #include "../../../QHamSolver/include/hilbert_space/double_u1.hpp"
#include "PackedLattice2D/LatticeRepresentation.hpp"
#include "hilbert_space.hpp"
// #include "helpers.hpp"


#ifdef USE_REAL_SECTORS
    using elem_ty = double;
#else
    using elem_ty = cpx;
#endif

const arma::vec _fermion_number = {1, 1, 0, 0, 0, 0, 0};
const arma::vec _fermion_up_number = {0, 1, 0, 0, 0, 0, 0};
const arma::vec _fermion_down_number = {1, 0, 0, 0, 0, 0, 0};
const arma::vec _dimer_number = {0, 0, 1, 1, 1, 1, 0};

//<! Symmetry contained in struct
struct SymmetrySectors{
    int kx_sym = 0;             // translation in the X direction
    int ky_sym = 0;             // translation in the Y direction
    int px_sym = 1;             // mirror symmetry in X direction (around Y axis)
    int py_sym = 1;             // mirror symmetry in Y direction (around X axis)
    int pd_sym = 1;             // reflection symmetry across diagonal (for Lx = Ly only)
    // int pa_sym;              // reflection symmetry across anti-diagonal (for Lx = Ly only) // NOT COMMUTING WITH THE OTHER SYMMETRIES!
    // int rot_sym;             // rotation symmetry with \pi / 2 angle (for Lx = Ly only)     // NOT COMMUTING WITH THE OTHER SYMMETRIES!
    int z_sym = 1;              // spin-flip parity
    int Nup = 1;                // translation in the X direction
    int Ndown = 1;              // translation in the X direction
    
    void print(){
        std::cout << "Sectors:\n";
        printSeparated(std::cout, "\n", 20, false, "kx = " + std::to_string(this->kx_sym), "ky = " + std::to_string(this->ky_sym), "px = " + std::to_string(this->px_sym));
        printSeparated(std::cout, "\n", 20, false, "py = " + std::to_string(this->py_sym), "pd = " + std::to_string(this->pd_sym), "z = " + std::to_string(this->z_sym));
        printSeparated(std::cout, "\n", 20, false, "N\u2191 = " + std::to_string(this->Nup), "N\u2193 = " + std::to_string(this->Ndown));
    }
};

// using elem_ty = cpx;
/// @brief Fully anisotropic spin chain (XYZ) with point symmetries
class BosonicSpinfullQLM : 
    public QHS::hamiltonian_base<elem_ty, QHS::Symmetric2DLattice<QOps::particle::boson>>
{
    //<! ----------------------------------------------------- INHERIT TYPEDEFs FROM BASE
    typedef typename QHS::hamiltonian_base<elem_ty, QHS::Symmetric2DLattice<QOps::particle::boson>>::matrix        matrix;
    typedef typename QHS::hamiltonian_base<elem_ty, QHS::Symmetric2DLattice<QOps::particle::boson>>::sparse_matrix sparse_matrix;
    
    //<! ----------------------------------------------------- MODEL PARAMETERS
private:
    v_1d<QOps::generic_operator<cpx, PackedLattice2D>> symmetry_generators;    // list of symmetry generators

    SymmetrySectors syms;                   // structure holding symmetry sectors

    double _tau = 1;                        // hopping strength
    double _mass = 0.0;                     // mass term
    double _kappa = 0.0;                    // spin inversion term ( breaks U(1) spin symmetry )
    
    int Lx = 4;
    int Ly = 4;
    bool _use_symmetries = true;            // [temporary] boolean choose if use symmetries

    //<! ----------------------------------------------------- INITIALIZE MODEL
    virtual void init() override;
    void set_symmetry_generators();

public:
    //<! ----------------------------------------------------- CONSTRUCTORS
    BosonicSpinfullQLM() = default;
    BosonicSpinfullQLM(std::istream& os);
    BosonicSpinfullQLM(int _BC, unsigned int Lx, unsigned int Ly, double tau, double m, SymmetrySectors syms_in, 
        #ifdef USE_SYMMETRIES
        bool use_syms = true
        #else
        bool use_syms = false
        #endif
        );

    //<! ----------------------------------------------------- HAMILTONIAN BUILDERS
    virtual void create_hamiltonian() override;
    virtual sparse_matrix create_local_hamiltonian(int site) override;
    virtual void set_hamiltonian_elements(u64 k, elem_ty value, u64 new_idx) override;

    //<! ----------------------------------------------------- OVERRIDEN OPERATORS
    virtual std::ostream& write(std::ostream&) const override;
    virtual std::istream& read(std::istream&) override;
    
};


#endif