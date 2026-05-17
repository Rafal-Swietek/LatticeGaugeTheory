#include "includes/FermionicQLM_UI.hpp"

int outer_threads = 1;
int num_of_threads = 1;

 /* Overloading * operator */
std::string operator*(std::string a, int b) 
{
    string output = "";
    while (b--) {
        output += a;
    }
    return output;
};

namespace FermionicQLM_UI{



void ui::make_sim(){
    printAllOptions();

    this->ptr_to_model = create_new_model_pointer();
    
    // return;
    // compare_energies(); return;
    
	clk::time_point start = std::chrono::system_clock::now();
    switch (this->fun)
	{
	case 0: 
		diagonalize(); 
		break;
	case 1:
		spin_impurity();
        break;
	case 2:
		spin_impurity();
		break;
	default:
		#define generate_scaling_array(name) arma::linspace(this->name, this->name + this->name##s * (this->name##n - 1), this->name##n)
        #define for_loop(param, var) for (auto& param : generate_scaling_array(var))

        for_loop(system_size, L){ 
            for_loop(taux, tau){    
            for_loop(massx, mass){    
        {
            this->L = system_size;
            this->tau = taux;
            this->mass = massx;
            this->site = this->L / 2.;
            
            this->reset_model_pointer();
            auto start_loop = std::chrono::system_clock::now();
            std::cout << " - - START NEW ITERATION:\t\t par = "; // simulation end
            printSeparated(std::cout, "\t", 16, true, this->L, this->tau, this->mass);
            spin_impurity();

            std::cout << "\t\t - - - - - - FINISHED ITERATION IN : " << tim_s(start_loop) << " seconds\n\t\t\t Total time : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
        }}}}
        std::cout << "Add default function" << std::endl;
	}
	std::cout << " - - - - - - FINISHED CALCULATIONS IN : " << tim_s(start) << " seconds - - - - - - " << std::endl; // simulation end
}



// -------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------- MODEL DEPENDENT FUNCTIONS
/// @brief Compare energie spactra for full model and all symmetry sectors combined
void ui::spin_impurity()
{          
    std::string dir = this->saving_dir + "SpinImpurity_Test" + kPSep;
	createDirs(dir);
	
	size_t dim = this->ptr_to_model->get_hilbert_size();
	std::string info = this->set_info();

    size_t dim_cut = 1e5;
    arma::sp_cx_mat H = cast_cx_sparse(this->ptr_to_model->get_hamiltonian());
	const size_t size = dim > dim_cut? this->l_steps : dim;

	clk::time_point start = std::chrono::system_clock::now();
    clk::time_point start_tot = start;
    if(dim > dim_cut){
        this->ptr_to_model->diag_sparse(this->l_steps, this->l_bundle, this->tol, this->seed);	
    }
    else{
        this->ptr_to_model->diagonalization();
    }
    const arma::vec E = this->ptr_to_model->get_eigenvalues();
    const auto& V = this->ptr_to_model->get_eigenvectors();
    std::cout << " - - - - - - finished diagonalization in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();

    PackedLattice2D initial_state(this->Lx, this->Ly);
    initial_state.fill_all(UP_MONOMER);
    // int N = std::min(this->)
    
    int x0 = this->Lx/2;
    int y0 = this->Ly/2;
    int nd = 1;
    initial_state.set_block(x0, y0, DOWN_MONOMER);
    for(int yp : {-1, 1}){
        for(int xp : {-1, 1}){
            if(nd >= this->syms.Ndown) continue;
            initial_state.set_block(x0 + xp, y0 + yp, DOWN_MONOMER);
            nd++;
        }
    }
    while(nd < this->syms.Ndown){
        initial_state.set_block(site + nd, DOWN_MONOMER);
    }
    

    auto _hilbert_space = this->ptr_to_model->get_model_ref().get_hilbert_space();
    u64 index = _hilbert_space.find(initial_state);
    initial_state.print_state();
    std::cout << "\n" << index << std::endl << std::endl;
    arma::Col<element_type> coeff = V.row(index).t();

    E.save(arma::hdf5_name(dir + info + ".hdf5", "energies"));
	// auto H = this->ptr_to_model->get_dense_hamiltonian();
    // H = arma::abs(H);
    // H.save(arma::hdf5_name(dir + info + ".hdf5", "Hamiltonian", arma::hdf5_opts::append));


    // Lochscmidt echo, total Nd on monomeres, Nd Nd correlation function, nd(k, t) ? energies as function of m

    arma::cx_vec init_state(dim, arma::fill::zeros);
    init_state(index) = 1.0;
    std::cout << " - - - - - - finished setting up initial state in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();
    
    arma::vec _ndown_density(size, arma::fill::zeros);
    arma::vec _nmonomer_density(size, arma::fill::zeros);
    for(int k = 0; k < size; k++)
    {
        arma::Col<ui::element_type> state = this->ptr_to_model->get_eigenState(k);
        arma::mat _ndown_density_state(this->Ly, this->Lx, arma::fill::zeros);
        arma::mat _nmonomer_density_state(this->Ly, this->Lx, arma::fill::zeros);
        for(long alfa = 0; alfa < dim; alfa++)
        {
            auto base_state = _hilbert_space(alfa);
        // #pragma omp parallel for
            for(std::size_t y = 0; y < this->Ly; y++){
                for(std::size_t x = 0; x < this->Lx; x++){
                    auto _conf = base_state.get_block(x, y);
                    _ndown_density_state(y, x) += _fermion_down_number(_conf) * std::abs(state(alfa)) * std::abs(state(alfa));
                    _nmonomer_density_state(y, x) += _fermion_number(_conf) * std::abs(state(alfa)) * std::abs(state(alfa));
                }
            }
        }
	    {
            _ndown_density_state.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDown_State/n = " + std::to_string(k), arma::hdf5_opts::append));
            _nmonomer_density_state.save(arma::hdf5_name(dir + info + ".hdf5", "Monomer_State/n = " + std::to_string(k), arma::hdf5_opts::append));
        }
        _ndown_density(k) = arma::accu( _ndown_density_state ) / (this->Lx * this->Ly);
        _nmonomer_density(k) = arma::accu( _nmonomer_density_state ) / (this->Lx * this->Ly);
    }
    _ndown_density.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDownDensity_states", arma::hdf5_opts::append));
    _nmonomer_density.save(arma::hdf5_name(dir + info + ".hdf5", "MonomerDensity_states", arma::hdf5_opts::append));
    std::cout << " - - - - - - finished eigenstate measures in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();

    arma::vec times = arma::linspace(0, 100, 10001);
    this->dt = times(1) - times(0);
    _ndown_density = arma::vec(times.size(), arma::fill::zeros);
    _nmonomer_density = arma::vec(times.size(), arma::fill::zeros);
    
    arma::vec _ndown_spreading(times.size(), arma::fill::zeros);

    arma::vec return_probability(times.size(), arma::fill::zeros);

    times.save(arma::hdf5_name(dir + info + ".hdf5", "times", arma::hdf5_opts::append));
    for(long t_idx = 0; t_idx < times.size(); t_idx++)
    {
        double time = times(t_idx);
        if(dim <= dim_cut){
            init_state = arma::cx_vec(dim, arma::fill::zeros);
            for(long alfa = 0; alfa < dim; alfa++)
            {
                // auto state = V.col(alfa);
                arma::Col<ui::element_type> state = this->ptr_to_model->get_eigenState(alfa);
                init_state += std::exp(-1i * time * E(alfa)) * state * coeff(alfa);
            }
            init_state = arma::normalise(init_state);
        }
        return_probability(t_idx) = std::abs(init_state(index) * init_state(index));

        arma::mat _ndown_density_t(this->Ly, this->Lx, arma::fill::zeros);
        arma::mat _nmonomer_density_t(this->Ly, this->Lx, arma::fill::zeros);
        for(long alfa = 0; alfa < dim; alfa++)
        {
            auto base_state = _hilbert_space(alfa);
        // #pragma omp parallel for
            for(std::size_t y = 0; y < this->Ly; y++){
                for(std::size_t x = 0; x < this->Lx; x++){
                    auto _conf = base_state.get_block(x, y);
                    _ndown_density_t(y, x) += _fermion_down_number(_conf) * std::abs(init_state(alfa)) * std::abs(init_state(alfa));
                    _nmonomer_density_t(y, x) += _fermion_number(_conf) * std::abs(init_state(alfa)) * std::abs(init_state(alfa));
                    
                    _ndown_spreading(t_idx) += ((x - x0)*(x - x0) + (y - y0)*(y - y0)) * _fermion_down_number(_conf) * std::abs(init_state(alfa)) * std::abs(init_state(alfa));
                }
            }
        }
        _ndown_density(t_idx) = arma::accu( _ndown_density_t ) / (this->Lx * this->Ly);
        _nmonomer_density(t_idx) = arma::accu( _nmonomer_density_t ) / (this->Lx * this->Ly);

        if(dim > dim_cut){
            lanczos::Lanczos<cpx> lancz(H, this->l_steps, -1, this->tol, this->seed, true, false, init_state);
            lancz.time_evolution_step(init_state, this->dt);
        }
        
	    if( t_idx%20 == 0){
            _ndown_density_t.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDownDensity_snapshot/t = " + to_string_prec(time), arma::hdf5_opts::append));
            _nmonomer_density_t.save(arma::hdf5_name(dir + info + ".hdf5", "MonomerDensity_snapshot/t = " + to_string_prec(time), arma::hdf5_opts::append));
        }
        // arma::vec state_abs = arma::abs(init_state);
        // state_abs.save(arma::hdf5_name("GaugeViolation.hdf5", "state: t = " + to_string_prec(time), arma::hdf5_opts::append));
    }
    _ndown_density.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDownDensity", arma::hdf5_opts::append));
    _nmonomer_density.save(arma::hdf5_name(dir + info + ".hdf5", "MonomerDensity", arma::hdf5_opts::append));
    return_probability.save(arma::hdf5_name(dir + info + ".hdf5", "return_probability", arma::hdf5_opts::append));
    _ndown_spreading.save(arma::hdf5_name(dir + info + ".hdf5", "spreading", arma::hdf5_opts::append));

    std::cout << " - - - - - - finished time evolution in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    std::cout << " - - - - - - finished evolution of spin down particle test in : " << tim_s(start_tot) << " s - - - - - - " << std::endl; // simulation end
}


/// @brief 
/// @param skip 
/// @param sep 
/// @return 
std::string ui::set_info(std::vector<std::string> skip, std::string sep) const
{
        std::string name = "Lx=" + std::to_string(this->Lx) + \
            ",Ly=" + std::to_string(this->Ly) + \
            ",Nu=" + std::to_string(this->syms.Nup) + \
            ",Nd=" + std::to_string(this->syms.Ndown) + \
            ",tau=" + to_string_prec(this->tau) + \
            ",m=" + to_string_prec(this->mass);
        // #ifdef USE_SYMMETRIES
        // name += ",Sz=" + std::to_string(this->syms.Sz);
        // if(this->boundary_conditions == 0)      name += ",k=" + std::to_string(this->syms.k_sym);
        // if(this->k_real_sec(this->syms.k_sym))  name += ",p=" + std::to_string(this->syms.p_sym);
        // #else
        //     name += ",w=" + to_string_prec(this->w) + \
        //             ",edge=" + std::to_string((int)this->add_edge_fields) + \
        //             ",pb=" + std::to_string((int)this->add_parity_breaking);
        // #endif

		auto tmp = split_str(name, ",");
		std::string tmp_str = sep;
		for (int i = 0; i < tmp.size(); i++) {
			bool save = true;
			for (auto& skip_param : skip)
			{
				// skip the element if we don't want it to be included in the info
				if (split_str(tmp[i], "=")[0] == skip_param)
					save = false;
			}
			if (save) tmp_str += tmp[i] + ",";
		}
		tmp_str.pop_back();
		return tmp_str;
}

/// @brief Compare energie spactra for full model and all symmetry sectors combined
void ui::compare_energies()
{                

    // auto full_model = std::make_unique<QHS::QHamSolver<FermiHubbard>>(this->boundary_conditions, this->Lx, this->Ly, this->J, this->U, this->delta, this->F, this->syms.Nup, this->syms.Ndown);
    // full_model->diagonalization(false);
    // arma::vec E_dis = full_model->get_eigenvalues();// + this->J1 * (this->L - int(this->boundary_conditions)) * (3 + this->eta1 * this->eta1) / 8.;
    
    // v_1d<double> Esym;
	// v_1d<std::string> symms;
    // auto kernel = [&Esym, &symms, this](int k, int p)
    // {
    //     auto symmetric_model = std::make_unique<QHS::QHamSolver<FermiHubbard>>(this->boundary_conditions, this->Lx, this->Ly, this->J, this->U, this->delta, this->F, this->syms.Nup, this->syms.Ndown);
    //     symmetric_model->diagonalization(false);
    //     arma::vec E = symmetric_model->get_eigenvalues();
        
    //     Esym.insert(Esym.end(), std::make_move_iterator(E.begin()), std::make_move_iterator(E.end()));
    //     v_1d<std::string> temp_str = v_1d<std::string>(E.size(), "k=" + std::to_string(k) + ",p=" + to_string(p));
	// 	symms.insert(symms.end(), std::make_move_iterator(temp_str.begin()), std::make_move_iterator(temp_str.end()));
    // };
    // loopSymmetrySectors(kernel);

    
    // auto permut = sort_permutation(Esym, [](const double a, const double b)
	// 							   { return a < b; });
	// apply_permutation(Esym, permut);
	// apply_permutation(symms, permut);
	// std::cout << std::endl << Esym.size() << std::endl << E_dis.size() << std::endl;
	// printSeparated(std::cout, "\t", 20, true, "symmetry sector", "Energy sym", "Energy total", "difference");
	// for (int k = 0; k < min((int)E_dis.size(), (int)Esym.size()); k++){
    //     // if(std::abs(Esym[k] - E_dis(k)) > 1e-13)
	// 	    printSeparated(std::cout, "\t", 20, true, symms[k], Esym[k], E_dis(k), Esym[k] - E_dis(k));
    // }
}

// -------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------- IMPLEMENTATION OF UI

/// @brief Create unique pointer to model with current parameters in class
typename ui::model_pointer ui::create_new_model_pointer(){
    return std::make_unique<QHS::QHamSolver<FermionicQLM>>(this->boundary_conditions, this->Lx, this->Ly, this->syms.Nup, this->syms.Ndown, this->tau, this->mass);
}

/// @brief Reset member unique pointer to model with current parameters in class
void ui::reset_model_pointer(){
    this->ptr_to_model.reset(new QHS::QHamSolver<FermionicQLM>(this->boundary_conditions, this->Lx, this->Ly, this->syms.Nup, this->syms.Ndown, this->tau, this->mass) );
}

/// @brief 
/// @param argc 
/// @param argv 
ui::ui(int argc, char **argv)
{
    auto input = change_input_to_vec_of_str(argc, argv);			// change standard input to vec of strings
	input = std::vector<std::string>(input.begin()++, input.end()); // skip the first element which is the name of file
	// plog::init(plog::info, "log.txt");						    // initialize logger
	
    if (std::string option = this->getCmdOption(input, "-f"); option != "")
	    input = this->parse_input_file(option); // parse input from file
	
	this->parse_cmd_options((int)input.size(), input); // parse input from CMD directly
}


/// @brief 
/// @param argc 
/// @param argv 
void ui::parse_cmd_options(int argc, std::vector<std::string> argv)
{
    //<! set all general UI parameters
    FermiUI::parse_cmd_options(argc, argv);

    //<! set the remaining UI parameters
	std::string choosen_option = "";	

	#define set_param(name) choosen_option = "-" #name;                                 \
	                        this->set_option(this->name, argv, choosen_option);         \
                                                                                        \
	                        choosen_option = "-" #name "s";                             \
	                        this->set_option(this->name##s, argv, choosen_option);      \
                                                                                        \
	                        choosen_option = "-" #name "n";                             \
	                        this->set_option(this->name##n, argv, choosen_option, true);
    set_param(tau);
    set_param(mass);
    set_param(Lx);
    set_param(Ly);

    this->L = this->Lx * this->Ly;

    //<! SYMMETRIES
    choosen_option = "-Nup";
    this->set_option(this->syms.Nup, argv, choosen_option);
    choosen_option = "-Ndown";
    this->set_option(this->syms.Ndown, argv, choosen_option);
    
    if(this->syms.Ndown < 0)
        this->syms.Ndown = this->L - this->syms.Nup;
    if(this->syms.Nup < 0)
        this->syms.Nup = this->L - this->syms.Ndown;

    choosen_option = "-k";
    // this->set_option(this->syms.k_sym, argv, choosen_option);

    choosen_option = "-p";
    // this->set_option(this->syms.p_sym, argv, choosen_option);

    choosen_option = "-n";
    // this->set_option(this->syms.n, argv, choosen_option);

    // if(this->syms.n == 0){
    //     choosen_option = "-Sz";
    //     this->set_option(this->syms.Sz, argv, choosen_option);
    //     if(this->L % 2 == 1 && this->syms.Sz == 0.0)
    //         this->syms.Sz = _Spin;
    // } else{
    //     this->syms.Sz = _Spin * (this->syms.n - 0.5) * this->L;
    // }
    
    

    //<! FOLDER
    std::string folder = "." + kPSep + "results" + kPSep;
    // #ifdef USE_SYMMETRIES
    //     folder += "symmetries" + kPSep;
    // #else
    //     folder += "disorder" + kPSep;
    // #endif
    switch(this->boundary_conditions){
        case 0: folder += "PBC" + kPSep; break;
        case 1: folder += "OBC" + kPSep; break;
        case 2: folder += "ABC" + kPSep; break;
        default:
            folder += "PBC" + kPSep; 
            break;
        
    }
    folder = folder + this->dir_prefix;
    if (fs::create_directories(folder) || fs::is_directory(folder)) // creating the directory for saving the files with results
    	this->saving_dir = folder;									// if can create dir this is is
}


/// @brief 
void ui::set_default(){
    FermiUI::set_default();
    this->tau = 1.0;
	this->taus = 0.0;
	this->taun = 1;
    
    this->mass = 0.0;
	this->masss = 0.0;
	this->massn = 1;

    // this->syms.k_sym = 0;
    // this->syms.p_sym = 1;
    // this->syms.Sz = (this->L%2) * _Spin;
    // this->syms.n = 0;
}

/// @brief 
void ui::print_help() const {
    FermiUI::print_help();

    printf(" Flags for 2D FermiHubbard model:\n");
    printSeparated(std::cout, "\t", 20, true, "-J", "(double)", "some coupling strength");
    printSeparated(std::cout, "\t", 20, true, "-Js", "(double)", "step in coupling strength sweep");
    printSeparated(std::cout, "\t", 20, true, "-Jn", "(int)", "number of couplings in the sweep");
    printSeparated(std::cout, "\t", 20, true, "-delta", "(double)", "coupling to Q5 charge");
    printSeparated(std::cout, "\t", 20, true, "-deltas", "(double)", "step in coupling to Q5 charge sweep");
    printSeparated(std::cout, "\t", 20, true, "-deltan", "(int)", "number of couplings to Q5 in the sweep");
    // #ifdef USE_SYMMETRIES
        printSeparated(std::cout, "\t", 20, true, "-k", "(int)", "quasimomentum symmetry sector");
        printSeparated(std::cout, "\t", 20, true, "-p", "(int)", "parity symmetry sector");
        printSeparated(std::cout, "\t", 20, true, "-Sz", "(float)", "magnetization sector");

    // #else
    //     printSeparated(std::cout, "\t", 20, true, "-w", "(double)", "disorder strength from uniform distribution");
    //     printSeparated(std::cout, "\t", 20, true, "-ws", "(double)", "step in disorder strength sweep");
    //     printSeparated(std::cout, "\t", 20, true, "-wn", "(int)", "number of disorder in the sweep");

    //     printSeparated(std::cout, "\t", 20, true, "-seed", "(u64)", "randomness in position for coupling to grain");
    //     printSeparated(std::cout, "\t", 20, true, "-edge", "(boolean)", "add edge fields for SUSY (when no disorder, i.e. w=0)");
    //     printSeparated(std::cout, "\t", 20, true, "-pb", "(boolean)", "add parity breaking term on edge (when no disorder, i.e. w=0)");
    // #endif
	std::cout << std::endl;
}

/// @brief 
void ui::printAllOptions() const{
    FermiUI::printAllOptions();
    // std::cout << "H = \u03A3_j H_j" << std::endl << std::endl;
	// std::cout << "H_j = " << std::endl;

	std::cout << "------------------------------ CHOSEN 2D Fermi Hubbard OPTIONS:" << std::endl;
    std::cout 
          << "Lx  = " << this->Lx << std::endl
		  << "Lxn = " << this->Lxn << std::endl
		  << "Lxs = " << this->Lxs << std::endl
          << "Ly  = " << this->Ly << std::endl
		  << "Lyn = " << this->Lyn << std::endl
		  << "Lys = " << this->Lys << std::endl
		  << "Nup = " << this->syms.Nup << std::endl
		  << "Ndown = " << this->syms.Ndown << std::endl
		  << "\u03C4  = " << this->tau << std::endl
		  << "\u03C4n = " << this->taun << std::endl
		  << "\u03C4s = " << this->taus << std::endl
		  << "m  = " << this->mass << std::endl
		  << "mn = " << this->massn << std::endl
		  << "ms = " << this->masss << std::endl;

    // std::cout << "k  = " << this->syms.k_sym << std::endl;
    // std::cout << "p  = " << this->syms.p_sym << std::endl;
    // std::cout << "Sz  = " << this->syms.Sz << std::endl;
    // std::cout << "n  = " << this->syms.n << std::endl;
		                                            
        std::cout << std::endl;
    printSeparated(std::cout, "\t", 16, true, "----------------------------------------------------------------------------------------------------");
}   

};



    // arma::vec gauge_violation2(times.size(), arma::fill::zeros);
    // arma::vec density_forbidden2(times.size(), arma::fill::zeros);
    // for(long t_idx = 0; t_idx < times.size(); t_idx++){
    //     double time = times(t_idx);
    //     arma::cx_vec init_state(dim, arma::fill::zeros);
    //     for(long alfa = 0; alfa < dim; alfa++)
    //     {
    //         auto state = V.col(alfa);
    //         init_state += std::exp(-1i * time * E(alfa)) * state * coeff(alfa);
    //     }
    //     init_state = arma::normalise(init_state);
    //     arma::mat density(this->Ly, this->Lx, arma::fill::zeros);
    //     arma::mat gauge(std::ceil(this->Ly/2.), std::ceil(this->Lx/2.), arma::fill::zeros);
    //     arma::mat density_forbidden_loc(std::ceil(this->Ly/2.), std::ceil(this->Lx/2.), arma::fill::zeros);
    //     for(long alfa = 0; alfa < dim; alfa++)
    //     {
    //         u64 base_state = _hilbert_space(alfa);
        
    //     // #pragma omp parallel for
    //         for(long int j = 0; j < this->L; j++)
    //         {
    //             int n_particle = check_spin(base_state, j);
    //             int xcoord = j % this->Lx;
    //             int ycoord = j / this->Lx;
    //             density(ycoord, xcoord) += _fermion_number(n_particle) * std::abs(init_state(alfa)) * std::abs(init_state(alfa));
    //         }
    //         for(long int i = 0; i < std::ceil(this->Lx/2.); i+=1)
    //             for(long int j = 0; j < std::ceil(this->Ly/2.); j+=1)
    //                 density_forbidden_loc(j, i) = density(2*j + 1, 2*i + 1);

    //     // #pragma omp parallel for
    //         for(long int i = 0; i < std::ceil(this->Lx/2.); i+=1)
    //         {
    //             for(long int j = 0; j < std::ceil(this->Ly/2.); j+=1)
    //             {
    //                 int site = 2*i + Lx * 2*j;
    //                 int n_particle = check_spin(base_state, site);
    //                 if(n_particle == 2) n_particle--; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
    //                 if(n_particle > 2) n_particle-=2; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
    //                 auto neis = _lattice.get_neighbours(site);
    //                 for(auto& nei : neis)
    //                 {
    //                     if(nei > 0){
    //                         int n_dimer = check_spin(base_state, nei);
    //                         if(n_dimer == 2) n_dimer--; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
    //                         if(n_dimer > 2) n_dimer-=2; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
    //                         n_particle += n_dimer;
    //                     }
    //                 }
    //                 gauge(j, i) += (double)n_particle * std::abs(init_state(alfa)) * std::abs(init_state(alfa));
    //             }
    //         }
    //     }
    //     gauge_violation2(t_idx) = std::sqrt( arma::accu( arma::square(gauge - gauge_init) ) / (this->Lx * this->Ly / 4) );
    //     density_forbidden2(t_idx) = arma::accu( density_forbidden_loc ) / (this->Lx * this->Ly / 4);
    //     printSeparated(std::cout, "\t", 20, true, std::abs(gauge_violation(t_idx) - gauge_violation2(t_idx)), std::abs(density_forbidden(t_idx) - density_forbidden2(t_idx)));
    // }