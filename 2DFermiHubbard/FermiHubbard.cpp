#include "includes/config.hpp"
#include "../../QHamSolver/include/QHamSolver.h"
#include "includes/FermiHubbard.hpp"


//<! ---------------------------------------------------------------------------------------------------------------------------------------
//<! ------------------------------------------------------------------------------------------------------------------------ IMPLEMENTATION

//<! ------------------------------------------------------------------------------ CONSTRUCTORS

/// @brief Initialize model dependencies (symetries, hilbert space, ...)
void FermiHubbard::init()
{   
    _debug_start( clk::time_point start = std::chrono::system_clock::now(); )
    // set symmetry generators
    this->set_symmetry_generators();

    // // initialize hilbert space with PXP + U(1)
    // auto Translate = QOps::__builtins::translation(this->system_size, 1);
    // auto flip = QOps::__builtins::spin_flip_x(this->system_size);
    // const int Lx = this->system_size;
    // const float Sz = this->syms._Sz;
    // auto some_kernel = [&Translate, &flip, Sz, Lx](u64 n){
    //     // n = std::get<0>(flip(n));
    //     int num_particles = Sz / _Spin + Lx / 2;
    //     // printSeparated(std::cout, "\t", 20, true, n, boost::dynamic_bitset<>(Lx, n), num_particles);
    //     return ( !( (n) & std::get<0>( Translate(n) ) ) ) && (__builtin_popcountll(n) == num_particles);
    //     // return !( (n) & std::get<0>( Translate(n) ) );
    // };

    // auto _hilbert_FH     = U1U1Hilbert(this->system_size, this->syms.Nup, this->syms.Ndown);
    // auto _second_hilbert = this->_use_symmetries? QHS::point_symmetric( this->system_size, this->symmetry_generators, this->_boundary_condition, this->syms.k_sym, 0) :\
    //                                                 QHS::point_symmetric( this->system_size, v_1d<QOps::genOp>(), 1, 0, -1);
    
    // include input argument into symmetric subspace with global sector and then iterate over global sector to find reduced subspace.
    #ifdef NO_POINT_SYM
    this->_hilbert_space = U1U1Hilbert(this->system_size, this->syms.Ndown, this->syms.Nup);
    // this->_hilbert_space = _second_hilbert;
    #else
    this->_hilbert_space = tensor(
                                QHS::point_symmetric<QOps::particle::fermion>( this->system_size, v_1d<QOps::genOp>(), 1, 0, 0 ),
                                U1U1Hilbert(this->system_size, this->syms.Ndown, this->syms.Nup)
                                    );
    #endif
    this->dim = this->_hilbert_space.get_hilbert_space_size();
    // for(int k = 0; k < this->dim; k++){
    //     u64 state = this->_hilbert_space(k);
    //     print2D_state(state);
    // }
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
FermiHubbard::FermiHubbard(int _BC, unsigned int Lx, unsigned int Ly, double J, double U, double delta, double F, int Nup, int Ndown, bool use_syms)
{ 
    CONSTRUCTOR_CALL;

    this->_boundary_condition = _BC;
    this->Lx = Lx;
    this->Ly = Ly; 
    this->system_size = Lx * Ly;
    this->_J = J;
    this->_delta = delta;
    this->_U = U;
    this->_F = F;

    this->_lattice = lattice::lattice2D(this->Lx, this->Ly, 1-this->_boundary_condition);
    std::cout << "\t\tFinished generating lattice of dimension:\t Lx=" << this->Lx << ", Ly=" << this->Ly << std::endl;
    
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
FermiHubbard::FermiHubbard(std::istream& os)
    { os >> *this; }

/// @brief Set symmetry generators (among spin flips if fields perpendicular to spin axis are 0)
void FermiHubbard::set_symmetry_generators()
{   
    // parity symmetry
    // this->symmetry_generators.emplace_back(QOps::_parity_symmetry(this->system_size, this->syms.p_sym));
}

//<! ------------------------------------------------------------------------------ HAMILTONIAN BUILDERS
/// @brief Set hamiltonian matrix element given with value and new index
/// @param k current basis state
/// @param value value of matrix element
/// @param new_idx new index to be found in hilbert space
void FermiHubbard::set_hamiltonian_elements(u64 k, elem_ty value, u64 new_idx)
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
void FermiHubbard::create_hamiltonian()
{
    this->H = sparse_matrix(this->dim, this->dim);
	for (int j = 0; j < this->system_size; j++) 
        this->H += this->create_local_hamiltonian(j);
    
    // std::cout << this->H - this->H.t() << std::endl;
}


/// @brief Method to create hamiltonian within the class
/// @param site site index where the local hamiltonian acts
/// @return the local hamiltonian at site site
typename FermiHubbard::sparse_matrix FermiHubbard::create_local_hamiltonian(int site)
{
    sparse_matrix H_local(dim, dim);
    auto set_loc_hamiltonian_elements = [this, &H_local](u64 k, elem_ty value, u64 new_idx){
        u64 state;
        elem_ty sym_eig;
        try {
            #ifdef NO_POINT_SYM
                u64 idx = this->_hilbert_space.find(new_idx);
                H_local(idx, k) += value;
                H_local(k, idx) += value;
            #else
                std::tie(state, sym_eig) = this->_hilbert_space.find_matrix_element(new_idx, this->_hilbert_space.get_norm(k));
                #ifdef USE_REAL_SECTORS
                    H_local(state, k) += std::real(value * sym_eig);
                    H_local(k, state) += std::real(value * sym_eig);
                #else
                    H_local(state, k) += value * sym_eig;
                    H_local(k, state) += std::conj(value * sym_eig);
                #endif
            #endif
        } 
        catch (const std::exception& err) {
            std::cout << "Exception:\t" << err.what() << "\n";
            std::cout << "SHit ehhh..." << std::endl;
            printSeparated(std::cout, "\t", 14, true, new_idx, new_idx, this->_hilbert_space(k), value, sym_eig);
        }
    };
    int xcoord = site % this->Lx;
    int ycoord = site / this->Lx;
    auto check_spin = QOps::__builtins::get_digit(this->system_size);
    
    auto neis = this->_lattice.get_neighbours(site);
    #ifdef EXTRA_DEBUG
        std::cout << "\n ------------ Test Neigbhours in lattice -----------------" << std::endl;
        for(auto& nei : neis)
            printSeparated(std::cout, "\t", 20, true, site, xcoord, ycoord, nei);
        std::cout << "\n ---------------------------------------------------------" << std::endl;
    #endif
    if(this->_boundary_condition && (xcoord >= this->Lx - 1 && ycoord >= this->Ly - 1))
    {
        return sparse_matrix(dim, dim);
    } else {
        for (size_t k = 0; k < this->dim; k++) {
            u64 base_state = this->_hilbert_space(k);
            int n_i = check_spin(base_state, site);
            int density = _fermion_number(n_i);

            // On-site coulomb repulsion
            if(n_i == 3) H_local(k, k) += this->_U;
            
            // Staggered potential
            int sign1 = xcoord%2 == 0? 1 : -1;
            int sign2 = ycoord%2 == 0? 1 : -1;
            H_local(k, k) += this->_delta / 2.0 * (sign1 + sign2) * density;

            // Linear tilt in both x and y direction
            H_local(k, k) += this->_F * (xcoord + 1.5 * ycoord) * density;

            // Hop hop my boy
            for(auto& nei : neis)
            {
                if(nei >= 0){
                    int n_nei = check_spin(base_state, nei);
                    
                    // Hopping: |00>|10> -> |10>|00>, |01>|10> -> |11>|00>, |01>|11> -> |11>|01> and |00>|11> -> |10>|01>
                    //           |0>|2>  ->  |2>|0>,   |1>|2>  ->  |3>|0>,   |1>|3>  ->  |3>|1>  and  |0>|3>  ->  |2>|1>
                    if( (n_i <= 1 && n_nei >= 2))
                    {
                        auto [val, state_tmp]   = operators::fermions::spin_half::anihilate_up<elem_ty>(base_state, this->system_size, nei);
                        auto [val2, state]      = operators::fermions::spin_half::create_up<elem_ty>(state_tmp, this->system_size, site);   
                        if(std::abs(val * val2) > 0){
                            set_loc_hamiltonian_elements(k, -this->_J * val * val2, state);
                            // printSeparated(std::cout, "\t", 20, true, "-----ONE-----", k, site, nei, n_i, n_nei, val, val2);
                            // print2D_state(base_state);
                            // print2D_state(state_tmp);
                            // print2D_state(state);
                        }
                    }
                    // Hopping: |00>|01> -> |01>|00>, |10>|01> -> |11>|00>, |10>|11> -> |11>|10> and |00>|11> -> |01>|10>
                    //           |0>|1>  ->  |1>|0>,   |2>|1>  ->  |3>|0>,   |2>|3>  ->  |3>|2>  and  |0>|3>  ->  |1>|2>
                    if( (n_i%2 == 0 && n_nei == 1) || (n_i%2 == 0 && n_nei == 3))
                    {
                        auto [val, state_tmp]   = operators::fermions::spin_half::anihilate_down<elem_ty>(base_state, this->system_size, nei);
                        auto [val2, state]      = operators::fermions::spin_half::create_down<elem_ty>(state_tmp, this->system_size, site);   
                        if(std::abs(val * val2) > 0){
                            set_loc_hamiltonian_elements(k, -this->_J * val * val2, state);
                            // printSeparated(std::cout, "\t", 20, true, "-----TWO-----", k, site, nei, n_i, n_nei, val, val2);
                            // print2D_state(base_state);
                            // print2D_state(state_tmp);
                            // print2D_state(state);
                        }
                    }

                    // // Hopping: |00>|10> -> |10>|00>, |01>|10> -> |11>|00>, |01>|11> -> |11>|01> and |00>|11> -> |10>|01>
                    // //           |0>|2>  ->  |2>|0>,   |1>|2>  ->  |3>|0>,   |1>|3>  ->  |3>|1>  and  |0>|3>  ->  |2>|1>
                    // if( (n_nei <= 1 && n_i >= 2))
                    // {
                    //     auto [val, state_tmp]   = operators::fermions::spin_half::anihilate_up<elem_ty>(base_state, this->system_size, nei);
                    //     auto [val2, state]      = operators::fermions::spin_half::create_up<elem_ty>(state_tmp, this->system_size, site);   
                    //     if(std::abs(val * val2) > 0){
                    //         set_loc_hamiltonian_elements(k, -this->_J * val * val2, state);
                    //         printSeparated(std::cout, "\t", 20, true, "-----THREE-----", k, site, nei, n_i, n_nei, val * val2);
                    //         print2D_state(base_state);
                    //         print2D_state(state);
                    //     }
                    // }
                    // // Hopping: |00>|01> -> |01>|00>, |10>|01> -> |11>|00>, |10>|11> -> |11>|10> and |00>|11> -> |01>|10>
                    // //           |0>|1>  ->  |1>|0>,   |2>|1>  ->  |3>|0>,   |2>|3>  ->  |3>|2>  and  |0>|3>  ->  |1>|2>
                    // if( (n_nei%2 == 0 && n_i == 1) && (n_nei%2 == 0 && n_i == 3))
                    // {
                    //     auto [val, state_tmp]   = operators::fermions::spin_half::anihilate_down<elem_ty>(base_state, this->system_size, site);
                    //     auto [val2, state]      = operators::fermions::spin_half::create_down<elem_ty>(state_tmp, this->system_size, nei);
                    //     if(std::abs(val * val2) > 0){
                    //         set_loc_hamiltonian_elements(k, -this->_J * val * val2, state);
                    //         printSeparated(std::cout, "\t", 20, true, "-----FOUR-----", k, site, nei, n_i, n_nei, val * val2);
                    //         print2D_state(base_state);
                    //         print2D_state(state);
                    //     }
                    // }
                }
            }
        }
        return H_local;
    }
}


//<! ------------------------------------------------------------------------------ OVVERRIDEN OPERATORS AND OPERATOR KERNELS
/// @brief Read model parameters from input stream
/// @tparam U1_sector U(1) symmetry sector as teamplate input 
/// @param os input stream to read parameters
std::istream& FermiHubbard::read(std::istream& os)
{
    
    return os;
}

/// @brief Write hamiltonian to stream as human readable
/// @param os input stream to read parameters
std::ostream& FermiHubbard::write(std::ostream& os) const
{
    printSeparated(os, "\t", 16, true, "Model:", "GoldenChain spin chain");
    os << std::endl;
    printSeparated(os, "\t", 16, true, "Hamiltonian:", "H = \u03A3_j P_j (X_j X_j+1 + Y_j Y_j+1 + \u0394 Z_j Z_j+1) P_j+1 ]");
    printSeparated(os, "\t", 16, true, "----------------------------------------------------------------------------------------------------");
    printSeparated(os, "\t", 16, true, "Parameters:");
    printSeparated(os, "\t", 16, true, "L", this->system_size);
    printSeparated(os, "\t", 16, true, "J", this->_J);
    printSeparated(os, "\t", 16, true, "\u0394", this->_delta);
    
    // printSeparated(os, "\t", 16, true, "k", this->syms.k_sym);
    // printSeparated(os, "\t", 16, true, "p", this->syms.p_sym);
    // printSeparated(os, "\t", 16, true, "Sz", this->syms._Sz);

    printSeparated(os, "\t", 16, true, "----------------------------------------------------------------------------------------------------");
    printSeparated(os, "\t", 16, true, "----------------------------------------------------------------------------------------------------");

    return os;
}



//<! ------------------------------------------------------------------------------ ADDITIONAL METHODS FOR SYMMETRIC HAMILTONIAN



