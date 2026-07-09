#include "includes/FermionicQLM_UI.hpp"

int outer_threads = 10;
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
    compare_energies(); return;
    
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
		dimer_dynamics();
		break;
    case 3:
        dynamical_structure_factors();
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
            auto kernel = [this](int kx, int ky, int z, int px, int py, int pd)
            {
                this->syms.kx_sym = kx;
                this->syms.ky_sym = ky;
                this->syms.px_sym = px;
                this->syms.py_sym = py;
                this->syms.pd_sym = pd;
                this->reset_model_pointer();
                if(this->ptr_to_model->get_hilbert_size() > 0)
                    this->diagonalize();
            };
            loopSymmetrySectors(kernel); continue;
            // spin_impurity(); continue;
            // dimer_dynamics(); continue;


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

    arma::vec times = arma::linspace(0, 50, 1001);
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

/// @brief Compare energie spactra for full model and all symmetry sectors combined
void ui::dimer_dynamics()
{          
    std::string dir = this->saving_dir + "DimerDynamics_Test" + kPSep;
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
    
    int x0 = (this->Lx - 1.0)/2;
    int y0 = (this->Ly - 1.0)/2;
    initial_state.set_block(x0, y0, DIMER_RIGHT);
    initial_state.set_block(x0+1, y0, DIMER_LEFT);
    
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
    
    arma::vec _ndimer_density(size, arma::fill::zeros);
    arma::vec _ndimer_density_center(size, arma::fill::zeros);
    arma::vec _nmonomer_density(size, arma::fill::zeros);
    for(int k = 0; k < size; k++)
    {
        arma::Col<ui::element_type> state = this->ptr_to_model->get_eigenState(k);
        arma::mat _ndimer_density_state(this->Ly, this->Lx, arma::fill::zeros);
        arma::mat _nmonomer_density_state(this->Ly, this->Lx, arma::fill::zeros);
        for(long alfa = 0; alfa < dim; alfa++)
        {
            auto base_state = _hilbert_space(alfa);
        // #pragma omp parallel for
            for(std::size_t y = 0; y < this->Ly; y++){
                for(std::size_t x = 0; x < this->Lx; x++){
                    auto _conf = base_state.get_block(x, y);
                    _ndimer_density_state(y, x) += _dimer_number(_conf) * std::abs(state(alfa)) * std::abs(state(alfa));
                    _nmonomer_density_state(y, x) += _fermion_number(_conf) * std::abs(state(alfa)) * std::abs(state(alfa));
                }
            }
            auto _conf = base_state.get_block(x0, y0);
            if(_conf == DIMER_RIGHT)
                _ndimer_density_center(k) += std::abs(state(alfa)) * std::abs(state(alfa));
        }
	    {
            _ndimer_density_state.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDown_State/n = " + std::to_string(k), arma::hdf5_opts::append));
            _nmonomer_density_state.save(arma::hdf5_name(dir + info + ".hdf5", "Monomer_State/n = " + std::to_string(k), arma::hdf5_opts::append));
        }
        _ndimer_density(k) = arma::accu( _ndimer_density_state ) / (this->Lx * this->Ly);
        _nmonomer_density(k) = arma::accu( _nmonomer_density_state ) / (this->Lx * this->Ly);
    }
    _ndimer_density.save(arma::hdf5_name(dir + info + ".hdf5", "DimerDensity_states", arma::hdf5_opts::append));
    _nmonomer_density.save(arma::hdf5_name(dir + info + ".hdf5", "MonomerDensity_states", arma::hdf5_opts::append));
    _ndimer_density_center.save(arma::hdf5_name(dir + info + ".hdf5", "DimerCenterDensity_states", arma::hdf5_opts::append));
    std::cout << " - - - - - - finished eigenstate measures in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();

    arma::vec times = arma::linspace(0, 50, 1001);
    // arma::vec times = arma::logspace(1e-2, int(std::log10(dim)), 1001);
    this->dt = times(1) - times(0);
    _ndimer_density = arma::vec(times.size(), arma::fill::zeros);
    _nmonomer_density = arma::vec(times.size(), arma::fill::zeros);
    _ndimer_density_center = arma::vec(times.size(), arma::fill::zeros);
    
    arma::vec _ndimer_spreading(times.size(), arma::fill::zeros);
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

        arma::mat _nmonomer_density_t(this->Ly, this->Lx, arma::fill::zeros);
        for(long alfa = 0; alfa < dim; alfa++)
        {
            auto base_state = _hilbert_space(alfa);
        // #pragma omp parallel for
            for(std::size_t y = 0; y < this->Ly; y++){
                for(std::size_t x = 0; x < this->Lx; x++){
                    auto _conf = base_state.get_block(x, y);
                    _ndimer_density(t_idx) += _dimer_number(_conf) / (this->Lx * this->Ly) * std::abs(init_state(alfa)) * std::abs(init_state(alfa));
                    _nmonomer_density_t(y, x) += _fermion_number(_conf) * std::abs(init_state(alfa)) * std::abs(init_state(alfa));
                    
                    // if(x < this->Lx - 1 && y < this->Ly-1)
                    {
                        if(_conf == DIMER_UP){
                            _ndimer_spreading(t_idx) += ((x - (x0+0.5))*(x - (x0+0.5)) + (y+0.5 - y0)*(y+0.5 - y0)) * std::abs(init_state(alfa)) * std::abs(init_state(alfa));
                        } else if(_conf == DIMER_RIGHT){
                            _ndimer_spreading(t_idx) += ((x - x0)*(x - x0) + (y - y0)*(y - y0)) * std::abs(init_state(alfa)) * std::abs(init_state(alfa));
                        }

                    }
                }
            }

            auto _conf = base_state.get_block(x0, y0);
            if(_conf == DIMER_RIGHT)
                _ndimer_density_center(t_idx) += std::abs(init_state(alfa)) * std::abs(init_state(alfa));
        }
        // _ndimer_density(t_idx) = arma::accu( _ndimer_density_t ) / (this->Lx * this->Ly);
        _nmonomer_density(t_idx) = arma::accu( _nmonomer_density_t ) / (this->Lx * this->Ly);

        if(dim > dim_cut){
            lanczos::Lanczos<cpx> lancz(H, this->l_steps, -1, this->tol, this->seed, true, false, init_state);
            lancz.time_evolution_step(init_state, this->dt);
        }
        
	    if( t_idx%20 == 0){
            // _ndimer_density_t.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDimerDensity_snapshot/t = " + to_string_prec(time), arma::hdf5_opts::append));
            _nmonomer_density_t.save(arma::hdf5_name(dir + info + ".hdf5", "MonomerDensity_snapshot/t = " + to_string_prec(time), arma::hdf5_opts::append));
        }
        // arma::vec state_abs = arma::abs(init_state);
        // state_abs.save(arma::hdf5_name("GaugeViolation.hdf5", "state: t = " + to_string_prec(time), arma::hdf5_opts::append));
    }
    _ndimer_density.save(arma::hdf5_name(dir + info + ".hdf5", "DimerDensity", arma::hdf5_opts::append));
    _nmonomer_density.save(arma::hdf5_name(dir + info + ".hdf5", "MonomerDensity", arma::hdf5_opts::append));
    _ndimer_density_center.save(arma::hdf5_name(dir + info + ".hdf5", "DimerCenterDensity", arma::hdf5_opts::append));
    return_probability.save(arma::hdf5_name(dir + info + ".hdf5", "return_probability", arma::hdf5_opts::append));
    _ndimer_spreading.save(arma::hdf5_name(dir + info + ".hdf5", "spreading", arma::hdf5_opts::append));

    std::cout << " - - - - - - finished time evolution in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    std::cout << " - - - - - - finished evolution of dimer test in : " << tim_s(start_tot) << " s - - - - - - " << std::endl; // simulation end
}


void ui::dynamical_structure_factors()
{
	clk::time_point start = std::chrono::system_clock::now();

    std::string dir = this->saving_dir + "DynamicalStructureFactor" + kPSep;
	createDirs(dir);
    std::string dir_spin = dir + "Spin" + kPSep;
    std::string dir_scharge = dir + "Charge" + kPSep;
	createDirs(dir_spin, dir_scharge);

	size_t dim = this->ptr_to_model->get_hilbert_size();
	std::string info = this->set_info();
    std::size_t dim_cut = 5e3;

    auto _hilbert = this->ptr_to_model->get_model_ref().get_hilbert_space();
    auto _group = _hilbert.get_symmetry_group();


    arma::Col<elem_ty> GS;
    double dE, E0;
    arma::vec E;
    arma::SpMat<elem_ty> H;
    if(dim > dim_cut)
    {
        std::cout << " - - - - - - Large Hilbert-space, chosen Lanczos - - - - - - " << std::endl; // simulation end
        H = this->ptr_to_model->get_hamiltonian();
        auto lancz0 = lanczos::Lanczos<elem_ty, converge::energies>(H, 20, 1000, 1e-14, this->seed, 1);
        lancz0.diagonalization();
        E = lancz0.get_eigenvalues();
        GS = lancz0.get_eigenstates().col(0);
        GS = arma::normalise(GS);
    } else {
        this->ptr_to_model->diagonalization();
        E = this->ptr_to_model->get_eigenvalues();
        GS = this->ptr_to_model->get_eigenState(0);
    }
    E0 = E(0);
    dE = E(E.size() - 1) - E0;
    std::cout << " - - - - - - finished ground state in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();

    //<! create path in Brillouin Zone
    std::vector<std::pair<double, double>> path;
    for(int i = 0; i <= this->Lx / 2; i++)
        path.push_back( std::make_pair(i, 0) );
    for(int i = 1; i <= this->Ly / 2; i++)
        path.push_back( std::make_pair(this->Lx/2, i) );
    
    if(this->Lx == this->Ly){
        for(int i = this->Lx/2-1; i >= 0; i--)
            path.push_back( std::make_pair(i, i) );
    } else {
        for(int i = this->Lx/2-1; i >= 0; i--)
            path.push_back( std::make_pair(i, this->Ly / 2) );
        for(int i = this->Ly/2-1; i >= 0; i--)
            path.push_back( std::make_pair(0, i) );
    }

    std::cout << " - - - - - - finished path construction in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();
    //<! Get entire hilbert space
    QHS::Symmetric2DLattice<QOps::particle::fermion> _hilbert_full(this->Lx, this->Ly, this->syms.Nup, this->syms.Ndown, v_1d<QOps::generic_operator<cpx, PackedLattice2D>>(), 0, 0, 0, 1, 1, true, false);
    u64 dim_tot = _hilbert_full.get_hilbert_space_size();
    // Step 1: Expand GS to full Hilbert space
    arma::Col<elem_ty> fullGS(dim_tot, arma::fill::zeros); //<! full GS
    for(u64 k = 0; k < dim; k++){
        auto state = _hilbert(k);          // representative state in sym. basis
        cpx _norm = _hilbert.get_norm(k);  // normalization of sym. state
        // loop over all group elements to reconstruct full state
        for(auto& g : _group){
            auto [new_state, val] = g(state);
            u64 idx = _hilbert_full.find(new_state);
            fullGS(idx) += GS(k) * val / _norm;  // accumulate with correct weight
        }
    }
    std::cout << " - - - - - - finished mapping GS into full space in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();

    double domega = dE / (double)E.size() / 4.;
    double eta = 2 * domega;
    arma::vec omegas = arma::regspace(0, domega, 2*dE);

    arma::mat Sqw(path.size(), omegas.size(), arma::fill::zeros);
    start = std::chrono::system_clock::now();
#pragma omp parallel for num_threads(outer_threads)
    for(int j = 0; j < path.size(); j++)
    {
        int kx, ky;
        std::tie(kx, ky) = path[j];
        const double qx = constants<double>::two_pi * kx / double(this->Lx);
        const double qy = constants<double>::two_pi * ky / double(this->Ly);

        // Step 2: Apply S^z_q to fullGS
        arma::Col<elem_ty> SqGS = fullGS;
    // #pragma omp parallel for
        for(u64 k = 0; k < dim_tot; k++){
            if(std::abs(SqGS(k)) < 1e-14) continue;  // skip zeros for efficiency
            auto state = _hilbert_full(k);
            cpx sq_val = 0.0;
            for(int x = 0; x < this->Lx; x++){
                for(int y = 0; y < this->Ly; y++){
                    auto _conf = state.get_block(x, y);
                    if( _conf == UP_MONOMER )           sq_val += std::exp(-1i * (qx*x + qy*y)) / 2.0;
                    else if( _conf == DOWN_MONOMER )    sq_val -= std::exp(-1i * (qx*x + qy*y)) / 2.0;
                }
            }
            sq_val /= double(this->L);  // normalization
            SqGS(k) = SqGS(k) * sq_val;
        }
        
        double norm_SqGS = arma::norm(SqGS);
        // std::cout << norm_SqGS << std::endl;
        if(norm_SqGS < 1e-12){
            std::cout << "S^z_q|GS> = 0 for q=(" << kx << "," << ky 
                    << ") -- likely S^z_tot=0 at q=(0,0)" << std::endl;
            continue;
        }
        SqGS /= norm_SqGS;  // normalize before checking symmetry eigenvalues
        // auto new_GS = SqGS;
        std::cout << " - - - - - - finished applying Sq to GS in full space in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
        start = std::chrono::system_clock::now();

        SymmetrySectors sym2 = this->syms;
        sym2.z_sym = -this->syms.z_sym;
        sym2.kx_sym = (sym2.kx_sym + kx) % this->Lx;
        sym2.ky_sym = (sym2.ky_sym + ky) % this->Ly;
        //<! px and py change only in \pi sector
        if( this->Lx % 2? false : kx  == int(this->Lx / 2) ) sym2.px_sym = -this->syms.px_sym;
        if( this->Ly % 2? false : ky  == int(this->Ly / 2) ) sym2.py_sym = -this->syms.py_sym;
        
        
        // auto Tx = QOps::_translation_X_symmetry(this->Lx, this->Ly, 0, 1);
        // auto Ty = QOps::_translation_Y_symmetry(this->Lx, this->Ly, 0, 1);
        // auto Px = QOps::_mirror_X_symmetry(this->Lx, this->Ly, 1);
        // auto Py = QOps::_mirror_Y_symmetry(this->Lx, this->Ly, 1);
        // auto D = QOps::_mirror_diagonal_symmetry(this->Lx, this->Ly, 1);
        // auto Z = QOps::_spin_inversion_symmetry(this->Lx, this->Ly, 1);
        
        // cpx Tx_SqGS = 0;
        // cpx Ty_SqGS = 0;
        // cpx Px_SqGS = 0;
        // cpx Py_SqGS = 0;
        // cpx Pd_SqGS = 0;
        // cpx Z_SqGS = 0;
        // for(u64 k = 0; k < dim_tot; k++)
        // {
        //     auto state = _hilbert_full(k);
        //     PackedLattice2D new_state(this->Lx, this->Ly);
        //     cpx val;
        //     std::tie(new_state, val) = Tx(state);
        //     u64 idx = _hilbert_full.find(new_state);
        //     Tx_SqGS += std::conj(SqGS(k)) * val * SqGS(idx);

        //     std::tie(new_state, val) = Ty(state);
        //     idx = _hilbert_full.find(new_state);
        //     Ty_SqGS += std::conj(SqGS(k)) * val * SqGS(idx);

        //     std::tie(new_state, val) = Px(state);
        //     idx = _hilbert_full.find(new_state);
        //     Px_SqGS += std::conj(SqGS(k)) * val * SqGS(idx);

        //     std::tie(new_state, val) = Py(state);
        //     idx = _hilbert_full.find(new_state);
        //     Py_SqGS += std::conj(SqGS(k)) * val * SqGS(idx);

        //     std::tie(new_state, val) = D(state);
        //     idx = _hilbert_full.find(new_state);
        //     Pd_SqGS += std::conj(SqGS(k)) * val * SqGS(idx);

        //     std::tie(new_state, val) = Z(state);
        //     idx = _hilbert_full.find(new_state);
        //     Z_SqGS += std::conj(SqGS(k)) * val * SqGS(idx);
        // }
        // std::cout << "Symmetry eigenvalues of the new state are for momentum q - ( " << kx << "," << ky << " ):" << std::endl;
        // printSeparated(std::cout, "\t", 16, true, "symmetry", "", "Old", "\u2192", "New");
        // printSeparated(std::cout, "\t", 16, true, "k_x", "", this->syms.kx_sym, "\u2192", -std::arg( Tx_SqGS ) * double(this->Lx) / constants<double>::two_pi);
        // printSeparated(std::cout, "\t", 16, true, "k_y", "", this->syms.ky_sym, "\u2192", -std::arg( Ty_SqGS ) * double(this->Ly) / constants<double>::two_pi);
        
        // if(kx == 0 || kx == this->Lx/2.) printSeparated(std::cout, "\t", 16, true, "px", "", this->syms.px_sym, "\u2192", std::real(Px_SqGS));
        // if(ky == 0 || ky == this->Ly/2.) printSeparated(std::cout, "\t", 16, true, "py", "", this->syms.py_sym, "\u2192", std::real(Py_SqGS));
        // if(kx == ky)                     printSeparated(std::cout, "\t", 16, true, "pd", "", this->syms.pd_sym, "\u2192", std::real(Pd_SqGS));
        // if(this->syms.Nup == this->syms.Ndown) printSeparated(std::cout, "\t", 16, true, "z", "", this->syms.z_sym, "\u2192", std::real(Z_SqGS));
        
        
        auto Sq_model = std::make_unique<QHS::QHamSolver<FermionicQLM>>(this->boundary_conditions, this->Lx, this->Ly, this->tau, this->mass, sym2);
        auto _new_hilbert = Sq_model->get_model_ref().get_hilbert_space();
        auto _new_group = _new_hilbert.get_symmetry_group();
        u64 dim_new = _new_hilbert.get_hilbert_space_size();
        // H = Sq_model->get_hamiltonian();

        std::cout << " - - - - - - finished generating model in shifted sector in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
        start = std::chrono::system_clock::now();
        // //<! Step 3: Cast the state into the symmetric subspace
        arma::Col<elem_ty> new_GS(dim_new);
        for(u64 k = 0; k < dim_new; k++)
        {
            auto state = _new_hilbert(k);
            cpx _norm = _new_hilbert.get_norm(k);
            for(auto& g : _group){
                auto [new_state, val] = g(state);
                u64 idx = _hilbert_full.find(new_state);
                new_GS(k) += SqGS(idx) * val / _norm;
            }
        }
        
        double norm_newGS = arma::norm(new_GS);
        if(norm_newGS < 1e-14){
            std::cout << "S^z_q|GS> = 0 in symmetric subspace for q=(" << kx << "," << ky 
                    << ") -- likely S^z_tot=0 at q=(0,0)" << std::endl;
            continue;
        }
        new_GS /= norm_newGS;
        std::cout << " - - - - - - finished casting Sq|GS> into new reduced sector in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
        start = std::chrono::system_clock::now();

        cpx alfa = arma::cdot(new_GS, new_GS);

        start = std::chrono::system_clock::now();
        if(dim > dim_cut){
            //<! Step 4: Build lanczos matrix starting from Sq|GS>
            H = Sq_model->get_hamiltonian();
            auto lancz = lanczos::Lanczos<elem_ty, converge::energies>(H, 20, 1000, 1e-14, this->seed, 1);
            lancz.build(new_GS);
            int lsteps = lancz.get_lanczossteps();
            auto HL = lancz.get_lanczos_matrix();
            std::cout << " - - - - - - finished lanczos iteration in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end

            //<! Step 5: Continous fraction expansion
        // #pragma omp parallel for
            for(int w = 0; w < omegas.size(); w++){
                cpx z = omegas(w) + 1i*eta + E0;
                cpx continous_fraction = z - HL(lsteps - 1, lsteps - 1);
                continous_fraction = continous_fraction - std::sqrt(continous_fraction*continous_fraction - 4.0*HL(lsteps - 1, lsteps - 2)*HL(lsteps - 1, lsteps - 2));
                continous_fraction /= 2.0;
                for (int m = lsteps - 2; m >= 0; m--) {
                    continous_fraction = z - HL(m, m) - HL(m, m + 1) * HL(m, m + 1) / continous_fraction;
                }
                Sqw(j, w) = -std::imag( alfa / continous_fraction ) / constants<double>::pi;
            }
        } else {
            Sq_model->diagonalization();
            std::cout << " - - - - - - finished ED in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end

            arma::vec E_Sq = Sq_model->get_eigenvalues();
            // arma::vec E_Sq = E;
            for(int n = 0; n < E_Sq.size(); n++)
            {
                auto state = Sq_model->get_eigenState(n);
                // auto state = this->ptr_to_model->get_eigenState(n);
                cpx overlap = arma::cdot(state, new_GS);
            // #pragma omp parallel for
                for(int w = 0; w < omegas.size(); w++)
                {
                    cpx z = omegas(w) + 1i*eta + E0 - E_Sq(n);
                    Sqw(j, w) += -std::imag( std::abs( overlap * std::conj(overlap) ) / z ) / constants<double>::pi;
                }
            }
        }
        std::cout << " - - - - - - finished continous fraction for (kx,ky) = ( " << kx << "," << ky << " ) in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
        start = std::chrono::system_clock::now();
    }
    E.save(arma::hdf5_name(dir_spin + info + ".hdf5", "lanczos energies"));
    omegas.save(arma::hdf5_name(dir_spin + info + ".hdf5", "omegas", arma::hdf5_opts::append));
    Sqw.save(arma::hdf5_name(dir_spin + info + ".hdf5", "Sqw", arma::hdf5_opts::append));
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

            if(this->boundary_conditions == 0){
                name += ",kx=" + std::to_string(this->syms.kx_sym);
                name += ",ky=" + std::to_string(this->syms.ky_sym);
            }
            if(this->syms.Nup == this->syms.Ndown) 
                name += ",z=" + std::to_string(this->syms.z_sym);
            if(this->kx_real_sec(this->syms.kx_sym)) 
                name += ",px=" + std::to_string(this->syms.px_sym);
            if(this->ky_real_sec(this->syms.ky_sym)) 
                name += ",py=" + std::to_string(this->syms.py_sym);
            if( (this->Lx == this->Ly) && (this->syms.kx_sym == this->syms.ky_sym) && this->syms.px_sym == this->syms.py_sym )
                name += ",pd=" + std::to_string(this->syms.pd_sym);

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

    auto full_model = std::make_unique<QHS::QHamSolver<FermionicQLM>>(this->boundary_conditions, this->Lx, this->Ly, this->tau, this->mass, this->syms, false);
    full_model->diagonalization(false);
    arma::vec E_full = full_model->get_eigenvalues();
    
    v_1d<double> Esym;
	v_1d<std::string> symms;
    int count = 0;
    double E0 = 0;
    std::string sectorr ;
    auto kernel = [&Esym, &symms, this, &count, &E0, &sectorr](int kx, int ky, int z, int px, int py, int pd)
    {
        SymmetrySectors sym2 = this->syms;
        sym2.kx_sym = kx;
        sym2.ky_sym = ky;
        sym2.z_sym = z;
        sym2.px_sym = px;
        sym2.py_sym = py;
        sym2.pd_sym = pd;
        auto symmetric_model = std::make_unique<QHS::QHamSolver<FermionicQLM>>(this->boundary_conditions, this->Lx, this->Ly, this->tau, this->mass, sym2, true);
        u64 dim = symmetric_model->get_hilbert_size();
        if(dim > 0){
            symmetric_model->diagonalization(false);
            arma::vec E = symmetric_model->get_eigenvalues();
            std::cout << "DONE\t" << E.size() << std::endl;

            Esym.insert(Esym.end(), std::make_move_iterator(E.begin()), std::make_move_iterator(E.end()));
            std::string _sec = "kx=" + std::to_string(kx) + ",ky=" + std::to_string(ky) + ",z=" + to_string(z) + ",px=" + to_string(px) + ",py=" + to_string(py) + ",pd=" + to_string(pd);
            v_1d<std::string> temp_str = v_1d<std::string>(E.size(), _sec);
            symms.insert(symms.end(), std::make_move_iterator(temp_str.begin()), std::make_move_iterator(temp_str.end()));

            if(count == 0){
                E0 = E(0);
                sectorr = _sec;
                count ++;
            }
            if(E(0) < E0 && std::abs(E(0) - E0) > 1e-14)
                sectorr = _sec;
        }
    };
    loopSymmetrySectors(kernel);
    std::cout << "Found ground state in:\n" << sectorr << std::endl << std::endl;
    
    auto permut = sort_permutation(Esym, [](const double a, const double b)
								   { return a < b; });
	apply_permutation(Esym, permut);
	apply_permutation(symms, permut);
	std::cout << std::endl << Esym.size() << std::endl << E_full.size() << std::endl;
	printSeparated(std::cout, "\t", 20, true, "symmetry sector", "Energy sym", "Energy total", "difference");
	for (int k = 0; k < min((int)E_full.size(), (int)Esym.size()); k++){
        // if(std::abs(Esym[k] - E_dis(k)) > 1e-13)
		    printSeparated(std::cout, "\t", 20, true, symms[k], Esym[k], E_full(k), Esym[k] - E_full(k));
    }
}

// -------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------- IMPLEMENTATION OF UI

/// @brief Create unique pointer to model with current parameters in class
typename ui::model_pointer ui::create_new_model_pointer(){
    return std::make_unique<QHS::QHamSolver<FermionicQLM>>(this->boundary_conditions, this->Lx, this->Ly, this->tau, this->mass, this->syms);
}

/// @brief Reset member unique pointer to model with current parameters in class
void ui::reset_model_pointer(){
    this->ptr_to_model.reset(new QHS::QHamSolver<FermionicQLM>(this->boundary_conditions, this->Lx, this->Ly, this->tau, this->mass, this->syms) );
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