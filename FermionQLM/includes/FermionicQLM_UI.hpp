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

        SymmetrySectors syms; // structure holding symmetry sectors
        
        /// @brief 
        /// @param ks 
        /// @return 
        bool kx_real_sec(int ks) const { return (ks == 0 || (this->Lx % 2? false : ks  == int(this->Lx / 2))) || (this->boundary_conditions == 1); };
        /// @brief 
        /// @param ks 
        /// @return 
        bool ky_real_sec(int ks) const { return (ks == 0 || (this->Ly % 2? false : ks  == int(this->Ly / 2))) || (this->boundary_conditions == 1); };
        

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
            
			const int kx_end = (this->boundary_conditions) ? 1 : this->Lx;
            const int ky_end = (this->boundary_conditions) ? 1 : this->Ly;
		// #pragma omp parallel for num_threads(outer_threads)// schedule(dynamic)
        // #ifdef USE_REAL_SECTORS
        //     for(int ks : (this->L%2 || this->boundary_conditions? v_1d<int>({0}) : v_1d<int>({0, (int)this->L/2})) ){
        // #else
		// 	for (int ks = 1; ks < k_end; ks++) {
        // #endif
            v_1d<int> z_sec = this->syms.Ndown == this->syms.Nup? v_1d<int>({1, -1}) : v_1d<int>({1});
			for (int kxs = 0; kxs < kx_end; kxs++) {
                for (int kys = 0; kys < ky_end; kys++) {
                    v_1d<int> pxsec = (kx_real_sec(kxs))? v_1d<int>({1}) : v_1d<int>({1});
                    v_1d<int> pysec = (ky_real_sec(kys))? v_1d<int>({1}) : v_1d<int>({1});
                    for(auto& zs : z_sec){
                        for(auto& pxs : pxsec){
                            for(auto& pys : pysec)
                            {   
                                std::vector<int> pd_sec = (this->Lx == this->Ly && kxs == kys && pxs == pys)? v_1d<int>({1, -1}) : v_1d<int>({1});
                                for(auto& pds : pd_sec){
                                    //<! create local lambda for multithreading enivorontment
                                    auto dummy_lambda = [&lambda](int kx, int ky, int z, int px, int py, int pd, auto... args){
                                        lambda(kx, ky, z, px, py, pd, args...);
                                    };
                                    dummy_lambda(kxs, kys, zs, pxs, pys, pds, args...);
                                }
                            }
                        }
                    }
                }
            }
		}
        

        void spin_impurity();
        void dimer_dynamics();
        void dynamical_structure_factors();
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