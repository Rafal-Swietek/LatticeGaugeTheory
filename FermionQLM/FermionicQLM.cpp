#include "includes/config.hpp"
#include "../../QHamSolver/include/QHamSolver.h"
#include "includes/FermionicQLM.hpp"


//<! ---------------------------------------------------------------------------------------------------------------------------------------
//<! ------------------------------------------------------------------------------------------------------------------------ IMPLEMENTATION

//<! ------------------------------------------------------------------------------ CONSTRUCTORS

/// @brief Initialize model dependencies (symetries, hilbert space, ...)
void FermionicQLM::init()
{   
    _debug_start( clk::time_point start = std::chrono::system_clock::now(); )
    // set symmetry generators
    this->set_symmetry_generators();

    this->_hilbert_space = QHS::Symmetric2DLattice<QOps::particle::fermion>(this->Lx, this->Ly, this->syms.Nup, this->syms.Ndown, v_1d<QOps::generic_operator<cpx, PackedLattice2D>>(), this->_boundary_condition );
    this->dim = this->_hilbert_space.get_hilbert_space_size();
    
    // auto mapy = this->_hilbert_space.get_mapping();
    // for(auto& x : mapy){
    //     x.print_state();
    //     std::cout << std::endl;
    // }
    // std::cout << mapy.size() << std::endl;
    _debug_end( std::cout << "\t\tFinished setting generating reduced basis (U(1) x point symmetries) with size:\t dim=" << this->dim << "\tin " << tim_s(start) << " seconds" << std::endl; )

    // create hamiltonian
    _debug_start( start = std::chrono::system_clock::now(); )
    this->create_hamiltonian();
    _debug_end( std::cout << "\t\tFinished generating Hamiltonian in " << tim_s(start) << " seconds" << std::endl; )
}

/// @brief Constructor for GoldenChain model class
/// @param _BC boundary condition
/// @param L system size
/// @param J model parameter
/// @param c coupling to Q5
/// @param ksym quasimomentum symmetry sector
/// @param rsym other Z_2 symmetry sector
/// @param zzsym spin flip in Z symmetry sector
/// @param use_syms use symmetric code?
FermionicQLM::FermionicQLM(int _BC, unsigned int Lx, unsigned int Ly, int Nup, int Ndown, double tau, double m, bool use_syms)
{ 
    CONSTRUCTOR_CALL;

    this->_boundary_condition = _BC;
    this->Lx = Lx;
    this->Ly = Ly; 
    this->system_size = Lx * Ly;
    this->_tau = tau;
    this->_mass = m;
    
    this->_use_symmetries = use_syms;
    //<! symmetries
    // this->syms.k_sym = ksym;
    // this->syms.p_sym = psym;
    // this->syms._Sz = Sz;
    this->syms.Nup = Nup;
    this->syms.Ndown = Ndown;

    // #ifdef USE_REAL_SECTORS
    //     if(this->_boundary_condition == 0){ // only for PBC
    //         bool is_k_sector_real = (std::abs(two_pi * ksym / this->system_size) < 1e-4) || (std::abs(two_pi * ksym / this->system_size - pi) < 1e-4);
    //         _assert_(is_k_sector_real, NOT_ALLOWED_SYM_SECTOR "\n\t\tMatrix type is real due to USE_REAL_SECTORS macro, but quasimomentum sector is complex, i.e. k != 0, pi");
    //     }
    // #endif

    this->init(); 
}

/// @brief Constructor from input stream
/// @param os input stream
FermionicQLM::FermionicQLM(std::istream& os)
    { os >> *this; }

/// @brief Set symmetry generators (among spin flips if fields perpendicular to spin axis are 0)
void FermionicQLM::set_symmetry_generators()
{   
    // parity symmetry
    // this->symmetry_generators.emplace_back(QOps::_parity_symmetry(this->system_size, this->syms.p_sym));
}

//<! ------------------------------------------------------------------------------ HAMILTONIAN BUILDERS
/// @brief Set hamiltonian matrix element given with value and new index
/// @param k current basis state
/// @param value value of matrix element
/// @param new_idx new index to be found in hilbert space
void FermionicQLM::set_hamiltonian_elements(u64 k, elem_ty value, u64 new_idx)
{   
    u64 state;
    elem_ty sym_eig;

    // try {
    //     std::tie(state, sym_eig) = this->_hilbert_space.find_matrix_element(new_idx, this->_hilbert_space.get_norm(k));
        
    //     H(state, k) += value * std::conj(sym_eig);
    // } 
    // catch (const std::exception& err) {
    //     std::cout << "Exception:\t" << err.what() << "\n";
    //     std::cout << "SHit ehhh..." << std::endl;
    //     printSeparated(std::cout, "\t", 14, true, new_idx, this->_hilbert_space(k), value, sym_eig);
    // }
}

/// @brief Method to create hamiltonian within the class
void FermionicQLM::create_hamiltonian()
{
    this->H = sparse_matrix(this->dim, this->dim);
	for (int j = 0; j < this->system_size; j++) 
        this->H += this->create_local_hamiltonian(j);
    
    // std::cout << this->H - this->H.t() << std::endl;
}


/// @brief Method to create hamiltonian within the class
/// @param site site index where the local hamiltonian acts
/// @return the local hamiltonian at site site
typename FermionicQLM::sparse_matrix FermionicQLM::create_local_hamiltonian(int site)
{
    sparse_matrix H_local(dim, dim);
    auto set_loc_hamiltonian_elements = [this, &H_local](u64 k, elem_ty value, const PackedLattice2D& new_state){
        try {
            // #ifdef NO_POINT_SYM
                u64 idx = this->_hilbert_space.find(new_state);
                H_local(idx, k) += value;
                H_local(k, idx) += value;
            // #else
            //     u64 state;
            //     elem_ty sym_eig;
            //     std::tie(state, sym_eig) = this->_hilbert_space.find_matrix_element(new_idx, this->_hilbert_space.get_norm(k));
            //     #ifdef USE_REAL_SECTORS
            //         H_local(state, k) += std::real(value * sym_eig);
            //         H_local(k, state) += std::real(value * sym_eig);
            //     #else
            //         H_local(state, k) += value * sym_eig;
            //         H_local(k, state) += std::conj(value * sym_eig);
            //     #endif
            // #endif
        } 
        catch (const std::exception& err) {
            std::cout << "Exception:\t" << err.what() << "\n";
            std::cout << "SHit ehhh..." << std::endl;
            // printSeparated(std::cout, "\t", 14, true, new_idx, new_idx, this->_hilbert_space(k), value, sym_eig);
        }
    };
    int xcoord = site % this->Lx;
    int ycoord = site / this->Lx;

    #ifdef EXTRA_DEBUG
        std::cout << "\n ------------ Test Neigbhours in lattice -----------------" << std::endl;
        for (QHS::Direction dir : {QHS::Direction::Right, QHS::Direction::Left, QHS::Direction::Up, QHS::Direction::Down}) {
            int nei = this->_hilbert_space.move(site, dir);
            printSeparated(std::cout, "\t", 20, true, site, xcoord, ycoord, nei);
        }   
        std::cout << "\n ---------------------------------------------------------" << std::endl;
    #endif

    std::vector<std::uint8_t> _dimers     = {DIMER_RIGHT, DIMER_LEFT, DIMER_UP, DIMER_DOWN};
    std::vector<std::uint8_t> _dimers_nei = {DIMER_LEFT, DIMER_RIGHT, DIMER_DOWN, DIMER_UP};
    if(this->_boundary_condition && (xcoord >= this->Lx - 1 && ycoord >= this->Ly - 1))
    {
        return sparse_matrix(dim, dim);
    } else {
        for (size_t k = 0; k < this->dim; k++) {
            PackedLattice2D base_state = this->_hilbert_space(k);
            std::uint8_t _bitj = base_state.get_block(xcoord, ycoord);

            //<! Mass terms
            if(_bitj == UP_MONOMER || _bitj == DOWN_MONOMER)
                H_local(k,k) += this->_mass;
            
            //<! kinetic term
            int idx = 0;
            for (QHS::Direction dir : {QHS::Direction::Right, QHS::Direction::Left, QHS::Direction::Up, QHS::Direction::Down}) {
				int nei = this->_hilbert_space.move(site, dir);
				if (nei >= 0) {
                    std::uint8_t _bitnei = base_state.get_block(nei);
                    if( (_bitj == UP_MONOMER && _bitnei == DOWN_MONOMER)
                        ||
                        (_bitj == DOWN_MONOMER && _bitnei == UP_MONOMER)
                    ){
                        PackedLattice2D state = base_state;
                        state.set_block(site, _dimers[idx]);
                        state.set_block(nei, _dimers_nei[idx]);
                        set_loc_hamiltonian_elements(k, -this->_tau, state);
                    }
				}
                idx++;
			}
        }
        return H_local;
    }
}


//<! ------------------------------------------------------------------------------ OVVERRIDEN OPERATORS AND OPERATOR KERNELS
/// @brief Read model parameters from input stream
/// @tparam U1_sector U(1) symmetry sector as teamplate input 
/// @param os input stream to read parameters
std::istream& FermionicQLM::read(std::istream& os)
{
    
    return os;
}

/// @brief Write hamiltonian to stream as human readable
/// @param os input stream to read parameters
std::ostream& FermionicQLM::write(std::ostream& os) const
{
    printSeparated(os, "\t", 16, true, "Model:", "GoldenChain spin chain");
    os << std::endl;
    printSeparated(os, "\t", 16, true, "Hamiltonian:", "H = \u03A3_j P_j (X_j X_j+1 + Y_j Y_j+1 + \u0394 Z_j Z_j+1) P_j+1 ]");
    printSeparated(os, "\t", 16, true, "----------------------------------------------------------------------------------------------------");
    printSeparated(os, "\t", 16, true, "Parameters:");
    printSeparated(os, "\t", 16, true, "L", this->system_size);
    // printSeparated(os, "\t", 16, true, "J", this->_J);
    // printSeparated(os, "\t", 16, true, "\u0394", this->_delta);
    
    // printSeparated(os, "\t", 16, true, "k", this->syms.k_sym);
    // printSeparated(os, "\t", 16, true, "p", this->syms.p_sym);
    // printSeparated(os, "\t", 16, true, "Sz", this->syms._Sz);

    printSeparated(os, "\t", 16, true, "----------------------------------------------------------------------------------------------------");
    printSeparated(os, "\t", 16, true, "----------------------------------------------------------------------------------------------------");

    return os;
}



//<! ------------------------------------------------------------------------------ ADDITIONAL METHODS FOR SYMMETRIC HAMILTONIAN



