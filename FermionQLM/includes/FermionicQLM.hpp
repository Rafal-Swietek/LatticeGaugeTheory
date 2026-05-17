#pragma once

//#define _XYZ_SYM

#ifndef _FermionicQLM
#define _FermionicQLM

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
const arma::vec _doublon_number = {0, 0, 1, 1, 1, 1, 0};

// using elem_ty = cpx;
/// @brief Fully anisotropic spin chain (XYZ) with point symmetries
class FermionicQLM : 
    public QHS::hamiltonian_base<elem_ty, QHS::Symmetric2DLattice<QOps::particle::fermion>>
{
    //<! ----------------------------------------------------- INHERIT TYPEDEFs FROM BASE
    typedef typename QHS::hamiltonian_base<elem_ty, QHS::Symmetric2DLattice<QOps::particle::fermion>>::matrix        matrix;
    typedef typename QHS::hamiltonian_base<elem_ty, QHS::Symmetric2DLattice<QOps::particle::fermion>>::sparse_matrix sparse_matrix;
    
    //<! ----------------------------------------------------- MODEL PARAMETERS
private:
    v_1d<QOps::genOp> symmetry_generators;    // list of symmetry generators

    //<! Symmetry contained in struct
    struct {
        int kx_sym;             // translation in the X direction
        int px_sym;             // mirror symmetry in X direction (around Y axis)
        int ky_sym;             // translation in the Y direction
        int py_sym;             // mirror symmetry in Y direction (around X axis)
        int z_sym;              // spin-flip parity
        int r_sym;              // reflection symmetry across diagonal (for Lx = Ly only)
        int rot_sym;            // rotation symmetry with \pi / 2 angle (for Lx = Ly only)
        int Nup;                // translation in the X direction
        int Ndown;              // translation in the X direction
    } syms;

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
    FermionicQLM() = default;
    FermionicQLM(std::istream& os);
    FermionicQLM(int _BC, unsigned int Lx, unsigned int Ly, int Nup, int Ndown, double tau, double m, bool use_syms = true);

    //<! ----------------------------------------------------- HAMILTONIAN BUILDERS
    virtual void create_hamiltonian() override;
    virtual sparse_matrix create_local_hamiltonian(int site) override;
    virtual void set_hamiltonian_elements(u64 k, elem_ty value, u64 new_idx) override;

    //<! ----------------------------------------------------- OVERRIDEN OPERATORS
    virtual std::ostream& write(std::ostream&) const override;
    virtual std::istream& read(std::istream&) override;
    
};


#endif