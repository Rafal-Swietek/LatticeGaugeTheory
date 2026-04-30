#include "includes/FermiHubbard_UI.hpp"

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

namespace FermiHubbard_UI{



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
		test_gauge_violation();
	case 2:
		test_gauge_violation_states();
		break;
	default:
		#define generate_scaling_array(name) arma::linspace(this->name, this->name + this->name##s * (this->name##n - 1), this->name##n)
        #define for_loop(param, var) for (auto& param : generate_scaling_array(var))

        for_loop(system_size, L){ 
            for_loop(Jx, J){    
            for_loop(Fx, F){    
        {
            this->L = system_size;
            this->J = Jx;
            this->F = Fx;
            this->site = this->L / 2.;
            
            this->reset_model_pointer();
            auto start_loop = std::chrono::system_clock::now();
            std::cout << " - - START NEW ITERATION:\t\t par = "; // simulation end
            printSeparated(std::cout, "\t", 16, true, this->L, this->J);
            test_gauge_violation();

            std::cout << "\t\t - - - - - - FINISHED ITERATION IN : " << tim_s(start_loop) << " seconds\n\t\t\t Total time : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
        }}}}
        std::cout << "Add default function" << std::endl;
	}
	std::cout << " - - - - - - FINISHED CALCULATIONS IN : " << tim_s(start) << " seconds - - - - - - " << std::endl; // simulation end
}



// -------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------- MODEL DEPENDENT FUNCTIONS
/// @brief Compare energie spactra for full model and all symmetry sectors combined
void ui::test_gauge_violation()
{          
    std::string dir = this->saving_dir + "GaugeViolation" + kPSep;
	createDirs(dir);
	
	size_t dim = this->ptr_to_model->get_hilbert_size();
	std::string info = this->set_info();

    arma::sp_cx_mat H = cast_cx_sparse(this->ptr_to_model->get_hamiltonian());
	const size_t size = dim > 1e5? this->l_steps : dim;

    elem_ty val;
    u64 initial_state = 0;
    std::tie(val, initial_state) = operators::fermions::spin_half::create_up<elem_ty>(initial_state, this->L, 0);
    std::tie(val, initial_state) = operators::fermions::spin_half::create_down<elem_ty>(initial_state, this->L, 2);
    std::tie(val, initial_state) = operators::fermions::spin_half::create_down<elem_ty>(initial_state, this->L, 8);
    std::tie(val, initial_state) = operators::fermions::spin_half::create_up<elem_ty>(initial_state, this->L, 10);
    
    auto check_spin = QOps::__builtins::get_digit(this->L);

    arma::mat state_in_lattice(this->Ly, this->Lx, arma::fill::zeros);
    for(long int j = 0; j < this->L; j++)
    {
        int n_particle = check_spin(initial_state, j);
        int xcoord = j % this->Lx;
        int ycoord = j / this->Lx;
        // printSeparated(std::cout, "\t", 20, true, state, j, xcoord, ycoord, n_particle);
        state_in_lattice(ycoord, xcoord) = _fermion_number(n_particle);
    }
    auto _hilbert_space = this->ptr_to_model->get_model_ref().get_hilbert_space();
    u64 index = _hilbert_space.find(initial_state);
    std::cout << initial_state << "\t" << index << "\t" << boost::dynamic_bitset<>(2*this->L, initial_state) << "\n" << state_in_lattice << std::endl << std::endl;

    arma::vec times = arma::linspace(0, 50, 5001);
    arma::vec gauge_violation(times.size(), arma::fill::zeros);
    arma::vec density_forbidden(times.size(), arma::fill::zeros);
    
    // this->ptr_to_model->diagonalization();
    // const arma::vec E = this->ptr_to_model->get_eigenvalues();
    // const auto& V = this->ptr_to_model->get_eigenvectors();
    // arma::Col<element_type> coeff = V.row(index).t();

    lattice::lattice2D _lattice(this->Lx, this->Ly, 1-this->boundary_conditions);
    // #pragma omp parallel for
    arma::mat gauge_init(std::ceil(this->Ly/2.), std::ceil(this->Lx/2.), arma::fill::zeros);
    for(long int i = 0; i < std::ceil(this->Lx/2.); i+=1)
    {
        for(long int j = 0; j < std::ceil(this->Ly/2.); j+=1)
        {
            int site = 2*i + Lx * 2*j;
            int n_particle = check_spin(initial_state, site);
            if(n_particle == 2) n_particle--; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
            if(n_particle > 2) n_particle-=2; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
            auto neis = _lattice.get_neighbours(site);
            for(auto& nei : neis)
            {
                if(nei > 0){
                    int n_dimer = check_spin(initial_state, nei);
                    // if(n_dimer == 2) n_dimer--; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
                    // if(n_dimer > 2) n_dimer-=2; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
                    if(n_dimer == 3) n_particle += 1;
                }
            }
            gauge_init(j, i) += (double)n_particle;
        }
    }

	gauge_init.save(arma::hdf5_name(dir + info + ".hdf5", "InitialGauge"));
    state_in_lattice.save(arma::hdf5_name(dir + info + ".hdf5", "InitialState", arma::hdf5_opts::append));
    
    times.save(arma::hdf5_name(dir + info + ".hdf5", "times", arma::hdf5_opts::append));
    
    arma::cx_vec init_state(dim, arma::fill::zeros);
    init_state(index) = 1.0;
    for(long t_idx = 0; t_idx < times.size(); t_idx++){
        double time = times(t_idx);
        // arma::cx_vec init_state(dim, arma::fill::zeros);
        // for(long alfa = 0; alfa < dim; alfa++)
        // {
        //     auto state = V.col(alfa);
        //     init_state += std::exp(-1i * time * E(alfa)) * state * coeff(alfa);
        // }
        // init_state = arma::normalise(init_state);



        // std::cout << arma::abs(init_state - init_state2).t();

        arma::mat density(this->Ly, this->Lx, arma::fill::zeros);
        arma::mat gauge(std::ceil(this->Ly/2.), std::ceil(this->Lx/2.), arma::fill::zeros);
        arma::mat density_forbidden_loc(std::ceil(this->Ly/2.), std::ceil(this->Lx/2.), arma::fill::zeros);
        for(long alfa = 0; alfa < dim; alfa++)
        {
            u64 base_state = _hilbert_space(alfa);
        
        // #pragma omp parallel for
            for(long int j = 0; j < this->L; j++)
            {
                int n_particle = check_spin(base_state, j);
                int xcoord = j % this->Lx;
                int ycoord = j / this->Lx;
                density(ycoord, xcoord) += _fermion_number(n_particle) * std::abs(init_state(alfa)) * std::abs(init_state(alfa));
            }
            for(long int i = 0; i < std::ceil(this->Lx/2.); i+=1)
                for(long int j = 0; j < std::ceil(this->Ly/2.); j+=1)
                    density_forbidden_loc(j, i) = density(2*j + 1, 2*i + 1);

        // #pragma omp parallel for
            for(long int i = 0; i < std::ceil(this->Lx/2.); i+=1)
            {
                for(long int j = 0; j < std::ceil(this->Ly/2.); j+=1)
                {
                    int site = 2*i + Lx * 2*j;
                    int n_particle = check_spin(base_state, site);
                    if(n_particle == 2) n_particle--; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
                    if(n_particle > 2) n_particle-=2; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
                    auto neis = _lattice.get_neighbours(site);
                    for(auto& nei : neis)
                    {
                        if(nei > 0){
                            int n_dimer = check_spin(base_state, nei);
                            if(n_dimer == 2) n_dimer--; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
                            if(n_dimer > 2) n_dimer-=2; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
                            n_particle += n_dimer;
                        }
                    }
                    gauge(j, i) += (double)n_particle * std::abs(init_state(alfa)) * std::abs(init_state(alfa));
                }
            }
        }
        gauge_violation(t_idx) = std::sqrt( arma::accu( arma::square(gauge - gauge_init) ) / (this->Lx * this->Ly / 4) );
        density_forbidden(t_idx) = arma::accu( density_forbidden_loc ) / (this->Lx * this->Ly / 4);

        lanczos::Lanczos<cpx> lancz(H, this->l_steps, -1, this->tol, this->seed, true, false, init_state);
        lancz.time_evolution_step(init_state, this->dt);
        // std::cout << "-------- Time: t = " << time << std::endl << std::endl;
        // std::cout << density << std::endl;
        // std::cout << "---" << std::endl;
        // std::cout << gauge << std::endl;
        // std::cout << "-------- ----------" << std::endl;
	    if( t_idx%50 == 0){
            gauge.save(arma::hdf5_name(dir + info + ".hdf5", "Gauge/t = " + to_string_prec(time), arma::hdf5_opts::append));
            density.save(arma::hdf5_name(dir + info + ".hdf5", "Density/t = " + to_string_prec(time), arma::hdf5_opts::append));
            density_forbidden_loc.save(arma::hdf5_name(dir + info + ".hdf5", "DensityForrbidden/t = " + to_string_prec(time), arma::hdf5_opts::append));
        }
        // arma::vec state_abs = arma::abs(init_state);
        // state_abs.save(arma::hdf5_name("GaugeViolation.hdf5", "state: t = " + to_string_prec(time), arma::hdf5_opts::append));
    }
    gauge_violation.save(arma::hdf5_name(dir + info + ".hdf5", "gauge_violation", arma::hdf5_opts::append));
    density_forbidden.save(arma::hdf5_name(dir + info + ".hdf5", "density_forbidden", arma::hdf5_opts::append));

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
}

/// @brief Compare energie spactra for full model and all symmetry sectors combined
void ui::test_gauge_violation_states()
{          
    std::string dir = this->saving_dir + "GaugeViolationState" + kPSep;
	createDirs(dir);
	
	size_t dim = this->ptr_to_model->get_hilbert_size();
	std::string info = this->set_info();

    arma::sp_cx_mat H = cast_cx_sparse(this->ptr_to_model->get_hamiltonian());
	const size_t size = dim > 1e5? this->l_steps : dim;

    auto check_spin = QOps::__builtins::get_digit(this->L);

    
    // this->ptr_to_model->diagonalization();
    // const arma::vec E = this->ptr_to_model->get_eigenvalues();
    // const auto& V = this->ptr_to_model->get_eigenvectors();
    // arma::Col<element_type> coeff = V.row(index).t();

    lattice::lattice2D _lattice(this->Lx, this->Ly, 1-this->boundary_conditions);
    this->ptr_to_model->diag_sparse(this->l_steps, this->l_bundle, this->tol, this->seed);	
    auto _hilbert_space = this->ptr_to_model->get_model_ref().get_hilbert_space();
    const arma::vec E = this->ptr_to_model->get_eigenvalues();
    E.save(arma::hdf5_name(dir + info + ".hdf5", "energies"));

    for(int k = 0; k < this->l_steps; k++)
    {
        arma::Col<ui::element_type> state = this->ptr_to_model->get_eigenState(k);
        arma::mat density(this->Ly, this->Lx, arma::fill::zeros);
        arma::mat gauge(std::ceil(this->Ly/2.), std::ceil(this->Lx/2.), arma::fill::zeros);
        arma::mat density_forbidden_loc(std::ceil(this->Ly/2.), std::ceil(this->Lx/2.), arma::fill::zeros);
        for(long alfa = 0; alfa < dim; alfa++)
        {
            u64 base_state = _hilbert_space(alfa);
        
        // #pragma omp parallel for
            for(long int j = 0; j < this->L; j++)
            {
                int n_particle = check_spin(base_state, j);
                int xcoord = j % this->Lx;
                int ycoord = j / this->Lx;
                density(ycoord, xcoord) += _fermion_number(n_particle) * std::abs(state(alfa)) * std::abs(state(alfa));
            }
            for(long int i = 0; i < std::ceil(this->Lx/2.); i+=1)
                for(long int j = 0; j < std::ceil(this->Ly/2.); j+=1)
                    density_forbidden_loc(j, i) = density(2*j + 1, 2*i + 1);

        // #pragma omp parallel for
            for(long int i = 0; i < std::ceil(this->Lx/2.); i+=1)
            {
                for(long int j = 0; j < std::ceil(this->Ly/2.); j+=1)
                {
                    int site = 2*i + Lx * 2*j;
                    int n_particle = check_spin(base_state, site);
                    if(n_particle == 2) n_particle--; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
                    if(n_particle > 2) n_particle-=2; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
                    auto neis = _lattice.get_neighbours(site);
                    for(auto& nei : neis)
                    {
                        if(nei > 0){
                            int n_dimer = check_spin(base_state, nei);
                            if(n_dimer == 2) n_dimer--; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
                            if(n_dimer > 2) n_dimer-=2; // 0 - empty, 1 - spin up, 2 - spin down, 3 - doublon
                            n_particle += n_dimer;
                        }
                    }
                    gauge(j, i) += (double)n_particle * std::abs(state(alfa)) * std::abs(state(alfa));
                }
            }
        }
	    {
            gauge.save(arma::hdf5_name(dir + info + ".hdf5", "Gauge/n = " + std::to_string(k), arma::hdf5_opts::append));
            density.save(arma::hdf5_name(dir + info + ".hdf5", "Density/n = " + std::to_string(k), arma::hdf5_opts::append));
            density_forbidden_loc.save(arma::hdf5_name(dir + info + ".hdf5", "DensityForrbidden/n = " + std::to_string(k), arma::hdf5_opts::append));
        }
        // arma::vec state_abs = arma::abs(init_state);
        // state_abs.save(arma::hdf5_name("GaugeViolation.hdf5", "state: t = " + to_string_prec(time), arma::hdf5_opts::append));
    }
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
            ",J=" + to_string_prec(this->J) + \
            ",U=" + to_string_prec(this->U) + \
            ",delta=" + to_string_prec(this->delta) + \
            ",F=" + to_string_prec(this->F);
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
    return std::make_unique<QHS::QHamSolver<FermiHubbard>>(this->boundary_conditions, this->Lx, this->Ly, this->J, this->U, this->delta, this->F, this->syms.Nup, this->syms.Ndown);
}

/// @brief Reset member unique pointer to model with current parameters in class
void ui::reset_model_pointer(){
    this->ptr_to_model.reset(new QHS::QHamSolver<FermiHubbard>(this->boundary_conditions, this->Lx, this->Ly, this->J, this->U, this->delta, this->F, this->syms.Nup, this->syms.Ndown) );
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
    set_param(J);
    set_param(U);
    set_param(delta);
    set_param(F);
    set_param(Lx);
    set_param(Ly);

    this->L = this->Lx * this->Ly;

    //<! SYMMETRIES
    choosen_option = "-Nup";
    this->set_option(this->syms.Nup, argv, choosen_option);
    choosen_option = "-Ndown";
    this->set_option(this->syms.Ndown, argv, choosen_option);

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
    this->J = 1.0;
	this->Js = 0.0;
	this->Jn = 1;
    
    this->delta = 0.0;
	this->deltas = 0.0;
	this->deltan = 1;

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
    std::cout << "H = \u03A3_j H_j" << std::endl << std::endl;
	std::cout << "H_j = P_j (X_j X_j+1 + Y_j Y_j+1 + \u0394 Z_j Z_j+1) P_j+1" << std::endl;

	std::cout << "------------------------------ CHOSEN ConstrainedXXZ OPTIONS:" << std::endl;
    std::cout 
          << "Lx  = " << this->Lx << std::endl
		  << "Lxn = " << this->Lxn << std::endl
		  << "Lxs = " << this->Lxs << std::endl
          << "Ly  = " << this->Ly << std::endl
		  << "Lyn = " << this->Lyn << std::endl
		  << "Lys = " << this->Lys << std::endl
		  << "Nup = " << this->syms.Nup << std::endl
		  << "Ndown = " << this->syms.Ndown << std::endl
		  << "J  = " << this->J << std::endl
		  << "Jn = " << this->Jn << std::endl
		  << "Js = " << this->Js << std::endl
		  << "\u0394  = " << this->delta << std::endl
		  << "\u0394n = " << this->deltan << std::endl
		  << "\u0394s = " << this->deltas << std::endl
		  << "U  = " << this->U << std::endl
		  << "Un = " << this->Un << std::endl
		  << "Us = " << this->Us << std::endl
		  << "F  = " << this->F << std::endl
		  << "Fn = " << this->Fn << std::endl
		  << "Fs = " << this->Fs << std::endl;

    // std::cout << "k  = " << this->syms.k_sym << std::endl;
    // std::cout << "p  = " << this->syms.p_sym << std::endl;
    // std::cout << "Sz  = " << this->syms.Sz << std::endl;
    // std::cout << "n  = " << this->syms.n << std::endl;
		                                            
        std::cout << std::endl;
    printSeparated(std::cout, "\t", 16, true, "----------------------------------------------------------------------------------------------------");
}   

};