#pragma once

#ifndef UI
#define UI

#include "config.hpp"



// #ifdef NO_POINT_SYM
//     #include "../../../QHamSolver/include/user_interface_dis.hpp"
//     #define FermiUI user_interface_dis<FermionicQLM>
// #else
    #include "../../../QHamSolver/include/user_interface_sym.hpp"
    #define FermiUI user_interface_sym<FermionicQLM>
// #endif  
#include "FermionicQLM.hpp"
// ----------------------------------------------------------------------------- UI QUANTUM SUN -----------------------------------------------------------------------------

namespace FermionicQLM_UI{
    class ui : public FermiUI{
        
    protected:
        double tau, mass, taus, masss;
        int taun, massn;
        
        int Lx, Ly, Lxs, Lys, Lxn, Lyn;

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
        
        /// @brief 
        /// @param ks 
        /// @return 
        bool k_real_sec(int ks) const { return (ks == 0 || ks == this->L / 2.) || (this->boundary_conditions == 1); };
        
        typedef typename FermiUI::model_pointer model_pointer;
        typedef typename FermiUI::element_type element_type;
        
    public:
		// ----------------------------------- CONSTRUCTORS
		ui() = default;
		ui(int argc, char** argv);														// standard constructor
		
        // ----------------------------------- PARSER FUNCTION FOR HELP
		void print_help() const override;
		
        // ----------------------------------- HELPING FUNCIONS
		void set_default() override;													// set default parameters
		void printAllOptions() const override;

		// ----------------------------------- REAL PARSER
		void parse_cmd_options(int argc, std::vector<std::string> argv) override;		// the function to parse the command line

		// ----------------------------------- SIMULATION
		void make_sim() override;														// make default simulation
        virtual model_pointer create_new_model_pointer() override;
	    virtual void reset_model_pointer() override;

        // ----------------------------------- MODEL DEPENDENT FUNCTIONS
        virtual std::string set_info(std::vector<std::string> skip = {}, 
										std::string sep = "_") const override;
        void compare_energies();
        // void compare_hamiltonian();
        // void check_symmetry_generators();

        template <
			typename callable, 
			typename... _types
			> 
			void loopSymmetrySectors(
			callable& lambda, //!< callable function
			_types... args										   //!< arguments passed to callable interface lambda
		) {
            
			const int k_end = (this->boundary_conditions) ? 1 : this->L;
		// #pragma omp parallel for num_threads(outer_threads)// schedule(dynamic)
        // #ifdef USE_REAL_SECTORS
        //     for(int ks : (this->L%2 || this->boundary_conditions? v_1d<int>({0}) : v_1d<int>({0, (int)this->L/2})) ){
        // #else
		// 	for (int ks = 1; ks < k_end; ks++) {
        // #endif
			for (int ks = 0; ks < k_end; ks++) {
				v_1d<int> psec = (k_real_sec(ks))? v_1d<int>({-1, 1}) : v_1d<int>({1});
                std::cout << ks << "\t\t" << psec << std::endl;
                for(auto& ps : psec){
                    //<! create local lambda for multithreading enivorontment
                    auto dummy_lambda = [&lambda](int k, int p, auto... args){
                        lambda(k, p, args...);
                    };
                    dummy_lambda(ks, ps, args...);
                }
			}
		}
        

        void spin_impurity();
        // ----------------------------------- HELPER FUNCTIONS
	    
        // virtual element_type
        // jE_mat_elem_kernel(
        //     const arma::Col<element_type>& state1, 
        //     const arma::Col<element_type>& state2,
        //     int i, u64 k, const QOps::_ifun& check_spin
        //     ) override;
    };
}

#endif