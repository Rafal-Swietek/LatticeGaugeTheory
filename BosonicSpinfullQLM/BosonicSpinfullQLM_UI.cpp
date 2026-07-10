#include "includes/BosonicSpinfullQLM_UI.hpp"

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

namespace BosonicSpinfullQLM_UI{



void ui::make_sim(){
    printAllOptions();

    bool stripe_or_shape_calculation = this->fun == 5 || this->fun == 6 || this->fun == 7;
    if( !stripe_or_shape_calculation )
        this->ptr_to_model = create_new_model_pointer();
    
    // arma::SpMat<ui::element_type> H = this->ptr_to_model->get_hamiltonian();
    // // // std::cout << arma::mat(H) << std::endl;
    // arma::vec E = arma::eig_sym(arma::mat(H));
    // arma::uvec x = arma::find(arma::abs(E) < 1e-10);
    // // std::cout << E << std::endl;
    // auto _hilbert_space = this->ptr_to_model->get_model_ref().get_hilbert_space();
    // u64 dim = _hilbert_space.get_hilbert_space_size();
    
    // struct PlaquetteTerm {
    //     int dx1, dy1;
    //     ui::element_type d1;

    //     int dx2, dy2;
    //     ui::element_type d2;

    //     double sign;
    // };

    // const std::array<PlaquetteTerm,4> plaquette = {{
    //     {0,0, DIMER_RIGHT, 1,0, DIMER_LEFT,  +1},
    //     {0,0, DIMER_UP,    0,1, DIMER_DOWN,  -1},
    //     {1,0, DIMER_UP,    1,1, DIMER_DOWN,  -1},
    //     {0,1, DIMER_RIGHT, 1,1, DIMER_LEFT,  +1}
    // }};

    // auto addPlaquette =
    // [&](PackedLattice2D &lat, int x, int y, const PlaquetteTerm &p)
    // {
    //     lat.set_block((x + p.dx1) % Lx, (y + p.dy1) % Ly, p.d1);
    //     lat.set_block((x + p.dx2) % Lx, (y + p.dy2) % Ly, p.d2);
    // };
    // PackedLattice2D initial_state(this->Lx, this->Ly);
    // int num_0_mode = 0;
    // // for (int x0 = 0; x0 < Lx - int(boundary_conditions); ++x0)
    // // for (int y0 = 0; y0 < Ly - int(boundary_conditions); ++y0)
    // // {
    // //     arma::Col<ui::element_type> state(dim, arma::fill::zeros);

    // //     for (const auto &p : plaquette)
    // //     {
    // //         initial_state.fill_all(UP_MONOMER);

    // //         addPlaquette(initial_state, x0, y0, p);

    // //         auto index = _hilbert_space.find(initial_state);
    // //         state(index) = p.sign;
    // //     }

    // //     state = arma::normalise(state);

    // //     printSeparated(std::cout, "\t", 20, true, x0, y0, "|| H|psi> || = ", arma::norm(H * state));

    // //     ++num_0_mode;
    // // }
    // int Nx = Lx - int(boundary_conditions);
    // int Ny = Ly - int(boundary_conditions);
    // auto overlap = [&](int x0, int y0, int x1, int y1)
    // {
    //     std::array<std::pair<int,int>,4> v0 = {{
    //         { x0 % Lx,             y0 % Ly             },
    //         { (x0+1) % Lx,         y0 % Ly             },
    //         { x0 % Lx,            (y0+1) % Ly          },
    //         { (x0+1) % Lx,        (y0+1) % Ly          }
    //     }};

    //     std::array<std::pair<int,int>,4> v1 = {{
    //         { x1 % Lx,             y1 % Ly             },
    //         { (x1+1) % Lx,         y1 % Ly             },
    //         { x1 % Lx,            (y1+1) % Ly          },
    //         { (x1+1) % Lx,        (y1+1) % Ly          }
    //     }};

    //     for (const auto &a : v0)
    //         for (const auto &b : v1)
    //             if (a == b)
    //                 return true;

    //     return false;
    // };
    // for (int x0 = 0; x0 < Nx; ++x0)
    // for (int y0 = 0; y0 < Ny; ++y0)
    // {
    //     int id0 = x0*Ny + y0;

    //     for (int x1 = 0; x1 < Nx; ++x1)
    //     for (int y1 = 0; y1 < Ny; ++y1)
    //     {
    //         int id1 = x1*Ny + y1;

    //         if (id1 <= id0)
    //             continue;

    //         // don't allow overlapping plaquettes
    //         if (overlap(x0, y0, x1, y1))
    //             continue;

    //         arma::Col<ui::element_type> state(dim, arma::fill::zeros);

    //         for (const auto &p0 : plaquette)
    //         for (const auto &p1 : plaquette)
    //         {
    //             initial_state.fill_all(UP_MONOMER);

    //             addPlaquette(initial_state, x0, y0, p0);
    //             addPlaquette(initial_state, x1, y1, p1);

    //             auto index = _hilbert_space.find(initial_state);
    //             if(index >= dim || index < 0)
    //                 initial_state.print_state();
    //             state(index) += p0.sign * p1.sign;
    //         }

    //         state = arma::normalise(state);

    //         printSeparated(std::cout, "\t", 20, true, x0, y0, x1, y1, "|| H|psi> || = ", arma::norm(H * state));
    //         ++num_0_mode;
    //     }
    // }
    // std::cout << "------------" << std::endl;
    // std::cout << "N0 = " << num_0_mode << "\t\t" << x.size() << std::endl;
    // return;

    // auto _hilbert_space = this->ptr_to_model->get_model_ref().get_hilbert_space();
    // size_t dim = this->ptr_to_model->get_hilbert_size();
    // arma::Mat<ui::element_type> H = this->ptr_to_model->get_dense_hamiltonian();
    // // std::cout << arma::Mat<ui::element_type>(H) << std::endl;
    // std::string info = this->set_info();
    // H.save(arma::hdf5_name(this->saving_dir + "Hamiltonian" + info + ".hdf5", "H"));
    // arma::vec DimerNumber(dim);
    // for(long alfa = 0; alfa < dim; alfa++)
    // {
    //     auto base_state = _hilbert_space(alfa);
    //     for(std::size_t y = 0; y < this->Ly; y++){
    //         for(std::size_t x = 0; x < this->Lx; x++){
    //             auto _conf = base_state.get_block(x, y);
    //             if(_conf == DIMER_RIGHT) 
    //                 DimerNumber(alfa) += 1.0;
    //             if(_conf == DIMER_UP) 
    //                 DimerNumber(alfa) += 1.0;
    //         }
    //     }
    // }
    // DimerNumber.save(arma::hdf5_name(this->saving_dir + "Hamiltonian" + info + ".hdf5", "DimerNumber", arma::hdf5_opts::append));
    // auto _hilbert_space = this->ptr_to_model->get_model_ref().get_hilbert_space();
    // u64 dim = _hilbert_space.get_hilbert_space_size();
    // auto _group = _hilbert_space.get_symmetry_group();
    // for(long alfa = 0; alfa < dim; alfa++)
    // {
    //     auto base_state = _hilbert_space(alfa);
    //     for(auto& g : _group){
    //         auto [Gstate, val] = g(base_state);
    //         Gstate.print_state();
    //         std::cout << std::endl << std::endl;
    //     }
    //     std::cout << " - - - - - -  - - - - - -  - - - - - - " << std::endl << std::endl;
    // }
    // arma::sp_cx_mat H = cast_cx_sparse(this->ptr_to_model->get_hamiltonian());
    // auto starta = std::chrono::system_clock::now();
    // auto a = arma::trace(H * H);
    // std::cout << " - - - - - - finished trace n = 2 in : " << tim_s(starta) << " s - - - - - - " << std::endl; // simulation end
    // starta = std::chrono::system_clock::now();
    // a = arma::trace(H * H * H * H * H * H * H * H * H * H);
    // std::cout << " - - - - - - finished trace n = 10 in : " << tim_s(starta) << " s - - - - - - " << std::endl; // simulation end
    // starta = std::chrono::system_clock::now();
    // a = arma::trace(H * H * H * H * H * H * H * H * H * H * H * H * H * H * H * H * H * H * H * H * H * H * H * H * H * H * H * H * H * H);
    // std::cout << " - - - - - - finished trace n = 30 in : " << tim_s(starta) << " s - - - - - - " << std::endl; // simulation end
    // return;

    // compare_energies(); return;
    // dynamical_structure_factors();
    // neel_order();
    // return;
    
    QHS::Symmetric2DLattice<> _hilbert_full;
    u64 dimtot;
    arma::SpMat<ui::element_type> DimerNumber_full(dimtot, dimtot);
    
    if( !stripe_or_shape_calculation )
    {
        _hilbert_full = QHS::Symmetric2DLattice<>(this->Lx, this->Ly, this->syms.Nup, this->syms.Ndown, v_1d<QOps::generic_operator<cpx, PackedLattice2D>>(), this->boundary_conditions, 0, 0, 1, 1, true, false);
        dimtot = _hilbert_full.get_hilbert_space_size();
        
        for(long alfa = 0; alfa < dimtot; alfa++)
        {
            auto base_state = _hilbert_full(alfa);
            for(std::size_t y = 0; y < this->Ly; y++){
                for(std::size_t x = 0; x < this->Lx; x++){
                    auto _conf = base_state.get_block(x, y);
                    if(_conf == DIMER_RIGHT) 
                        DimerNumber_full(alfa, alfa) += 1.0;
                    if(_conf == DIMER_UP) 
                        DimerNumber_full(alfa, alfa) += 1.0;
                }
            }
        }
    }

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
    case 4:
        neel_order();
        break;
    case 5:
        LiebLatticeStrip();
        break;
    case 6:
        LiebLatticeRandomShape();
        break;
    case 7:
        spectral_form_factor_shape();
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
            auto kernel = [this, &dimtot, &_hilbert_full, &DimerNumber_full](int kx, int ky, int z, int px, int py, int pd)
            {
                // this->syms.kx_sym = kx;
                // this->syms.ky_sym = ky;
                // this->syms.px_sym = px;
                // this->syms.py_sym = py;
                // this->syms.pd_sym = pd;
                // this->reset_model_pointer();
                printSeparated(std::cout, "\t", 16, true, kx, ky, z, px, py, pd);
                // arma::SpMat<ui::element_type> H = this->ptr_to_model->get_hamiltonian();
                // std::cout << arma::Mat<ui::element_type>(H) << std::endl;
                std::cout << " - - - - - -  - - - - - -  - - - - - - " << std::endl << std::endl;
                SymmetrySectors sym2 = this->syms;
                sym2.kx_sym = kx;
                sym2.ky_sym = ky;
                std::string info = "_Lx=" + std::to_string(this->Lx) + \
                    ",Ly=" + std::to_string(this->Ly) + \
                    ",Nu=" + std::to_string(this->syms.Nup) + \
                    ",Nd=" + std::to_string(this->syms.Ndown) + \
                    ",tau=" + to_string_prec(this->tau) + \
                    ",m=" + to_string_prec(this->mass);
                    #ifdef USE_SYMMETRIES
                        if(this->boundary_conditions == 0){
                            info += ",kx=" + std::to_string(kx);
                            info += ",ky=" + std::to_string(ky);
                        }
                        #ifndef USE_ONLY_TRANSLATION
                            if(this->syms.Nup == this->syms.Ndown) 
                                info += ",z=" + std::to_string(z);
                            if(this->kx_real_sec(kx)) 
                                info += ",px=" + std::to_string(px);
                            if(this->ky_real_sec(ky)) 
                                info += ",py=" + std::to_string(py);
                            if( (this->Lx == this->Ly) && (kx == ky) && px == py )
                                info += ",pd=" + std::to_string(pd);
                        #endif
                    #endif
                auto start = std::chrono::system_clock::now();
                auto Sq_model = std::make_unique<QHS::QHamSolver<BosonicSpinfullQLM>>(this->boundary_conditions, this->Lx, this->Ly, this->tau, this->mass, sym2);
                u64 dim = Sq_model->get_hilbert_size();
                if(dim > 0)
                {
                    auto _hilbert_space = Sq_model->get_model_ref().get_hilbert_space();
                    arma::SpMat<ui::element_type> U = _hilbert_space.symmetry_rotation(_hilbert_full);

                    arma::SpMat<ui::element_type> DimerNumber = U.t() * DimerNumber_full * U;
                    Sq_model->diagonalization(true);
                    std::cout << " - - - - - - finished diagonalization in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
                    start = std::chrono::system_clock::now();
                    arma::vec eigenvalues = Sq_model->get_eigenvalues();
                    std::string name = this->saving_dir + "DIAGONALIZATION" + kPSep + info + ".hdf5";
                    eigenvalues.save(arma::hdf5_name(name, "energies"));

                    arma::vec _dimer_num(dim, arma::fill::zeros);
                // #pragma omp parallel for
                    for(int k = 0; k < dim; k++)
                    {
                        arma::Col<ui::element_type> state = Sq_model->get_eigenState(k);
                        _dimer_num(k) = std::real( arma::cdot(state, DimerNumber * state) );
                    }
                    _dimer_num.save(arma::hdf5_name(name, "DimerNumber_n", arma::hdf5_opts::append));
                    std::cout << " - - - - - - finished dimer number in eigenstates in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
                }
            };
            loopSymmetrySectors(kernel); continue;
            // diagonalize(); continue;
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
void ui::LiebLatticeRandomShape(){
    
    const int NA = this->L;     // number of A sites

    std::mt19937 rng(this->seed);
    for(int real = 0; real < this->l_realis; real++){
        clk::time_point start = std::chrono::system_clock::now();
        clk::time_point start_re = std::chrono::system_clock::now();
        // ---------- random connected A cluster ----------
        std::vector<std::pair<int,int>> A{{0,0}};
        std::map<std::pair<int,int>,bool> occ;
        occ[{0,0}] = true;

        const int dx[4]={1,-1,0,0};
        const int dy[4]={0,0,1,-1};

        while((int)A.size()<NA)
        {
            auto s=A[rng()%A.size()];
            int d=rng()%4;

            std::pair<int,int> n={s.first+dx[d],s.second+dy[d]};

            if(!occ[n])
            {
                occ[n]=true;
                A.push_back(n);
            }
        }

        // ---------- create site list ----------
        std::map<std::pair<double,double>,int> id;
        std::vector<std::pair<double,double>> pos;
        std::vector<int> sub;          // 0=A,1=B,2=C

        for(auto a:A)
        {
            id[{a.first,a.second}] = pos.size();
            pos.push_back({a.first,a.second});
            sub.push_back(0);
        }

        for(auto a:A)
            if(occ[{a.first+1,a.second}])
            {
                id[{a.first+0.5,a.second}] = pos.size();
                pos.push_back({double(a.first)+0.5, double(a.second)});
                sub.push_back(1);
            }

        for(auto a:A)
            if(occ[{a.first,a.second+1}])
            {
                id[{a.first,a.second+0.5}] = pos.size();
                pos.push_back({double(a.first), double(a.second)+0.5});
                sub.push_back(2);
            }

        // ---------- coordinates ----------
        arma::vec x, y;
        x.set_size(pos.size());
        y.set_size(pos.size());

        for(size_t i=0;i<pos.size();i++)
        {
            x(i)=pos[i].first;
            y(i)=pos[i].second;
        }
        std::cout << " - - - - - - finished lattice construction in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
        start = std::chrono::system_clock::now();

        // ---------- Hamiltonian ----------
        arma::mat H(pos.size(),pos.size());
        std::vector<int> boundary;

        for(auto a:A)
        {
            int ia=id[{a.first,a.second}];
            int z=0;

            if(occ[{a.first-1,a.second}])
            {
                int j=id[{a.first-0.5,a.second}];
                H(ia,j)=H(j,ia)=-tau;
                z++;
            }

            if(occ[{a.first+1,a.second}])
            {
                int j=id[{a.first+0.5,a.second}];
                H(ia,j)=H(j,ia)=-tau;
                z++;
            }

            if(occ[{a.first,a.second-1}])
            {
                int j=id[{a.first,a.second-0.5}];
                H(ia,j)=H(j,ia)=-tau;
                z++;
            }

            if(occ[{a.first,a.second+1}])
            {
                int j=id[{a.first,a.second+0.5}];
                H(ia,j)=H(j,ia)=-tau;
                z++;
            }

            H(ia,ia)=z*mass;   // 4m bulk, 3m edge, 2m corner
            // printSeparated(std::cout, "\t", 20, true, ia, x(ia), y(ia), z);
            if (z < 4)
                boundary.push_back(ia);
        }
        std::cout << " - - - - - - finished matrix construction in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
        start = std::chrono::system_clock::now();

        arma::vec E; arma::mat V;
        arma::eig_sym(E, V, H);
        std::cout << " - - - - - - finished diagonalization in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
        start = std::chrono::system_clock::now();

        arma::vec PA(E.size());

        for (arma::uword n = 0; n < E.size(); n++)
            PA(n) = arma::accu(arma::square(arma::abs(V(arma::span(0,NA-1), n))));

        arma::vec Pboundary(E.size(), arma::fill::zeros);

        for (arma::uword n=0; n<E.size(); n++)
            for (int i : boundary)
                Pboundary(n) += std::norm(V(i,n));
        std::string dir = this->dir_prefix + "results" + kPSep + "LiebLatticeRandomShape" + kPSep + "realization=" + std::to_string(real + this->jobid) + kPSep;
        createDirs(dir);
        std::string info = "_L=" + std::to_string(this->L) + \
                        ",tau=" + to_string_prec(this->tau) + \
                        ",m=" + to_string_prec(this->mass);
        E.save(arma::hdf5_name(dir + info + ".hdf5", "E"));
        PA.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDownNumber", arma::hdf5_opts::append));
        Pboundary.save(arma::hdf5_name(dir + info + ".hdf5", "EdgeDownLoc", arma::hdf5_opts::append));

        arma::Col<int> sub_arma = arma::conv_to<arma::Col<int>>::from(sub);
        sub_arma.save(arma::hdf5_name(dir + info + ".hdf5", "sublattice", arma::hdf5_opts::append));

        // V.save(arma::hdf5_name(dir + info + ".hdf5", "V", arma::hdf5_opts::append));

        x.save(arma::hdf5_name(dir + info + ".hdf5", "x coords", arma::hdf5_opts::append));
        y.save(arma::hdf5_name(dir + info + ".hdf5", "y coords", arma::hdf5_opts::append));
        arma::vec energies = arma::linspace(E(0) - 2, E(E.size()-1) + 2, 3000);
        arma::vec DOS(energies.size(), arma::fill::zeros);
        arma::vec LDOS(energies.size(), arma::fill::zeros);

        // arma::vec initial_state(E.size(), arma::fill::zeros);
        // initial_state(0) = 1.0;
        double _eta = 0.02;
        double norm_inv = 1. / std::sqrt(constants<double>::two_pi * _eta*_eta);
        arma::rowvec coeff = V.row(0);
    #pragma omp for schedule(dynamic)
        for(long e_idx = 0; e_idx < energies.size(); e_idx++){
            // doubl
            // ldos(e_idx) += arma::sum(arma::square(coeff.t()) % _eta / ( arma::square(Esym))
            for(long alfa = 0; alfa < E.size(); alfa++){
                double om = E(alfa) - energies(e_idx);
                double gauss = std::exp( - om * om / (2.*_eta * _eta) ) * norm_inv;
                DOS(e_idx) += gauss;
                LDOS(e_idx) += gauss * std::norm(coeff(alfa));
            }
        }
        coeff.save(arma::hdf5_name(dir + info + ".hdf5", "coeff", arma::hdf5_opts::append));
        DOS.save(arma::hdf5_name(dir + info + ".hdf5", "dos", arma::hdf5_opts::append));
        LDOS.save(arma::hdf5_name(dir + info + ".hdf5", "LDOS", arma::hdf5_opts::append));
        energies.save(arma::hdf5_name(dir + info + ".hdf5", "dos_energies", arma::hdf5_opts::append));
        arma::Col<u64>({this->seed}).save(arma::hdf5_name(dir + info + ".hdf5", "seed", arma::hdf5_opts::append));
        std::cout << " - - - - - - finished eigenstate measures and DOS in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
        start = std::chrono::system_clock::now();
        // arma::vec initial_state(E.size(), arma::fill::zeros);
        // initial_state(0) = 1.0;
        arma::vec times = arma::logspace(-3, std::log10(50*this->L) + 1, 1000);
        arma::vec times_append({0, 0.1, 0.2, 0.5, 1, 1.5, 2., 3., 4., 6., 8., 10., 20., 30., 40., 50., 75., 100., 150., 200.});
        times = arma::sort(arma::join_cols(times, times_append));
        times.save(arma::hdf5_name(dir + info + ".hdf5", "times", arma::hdf5_opts::append));

        arma::cx_mat init_state(E.size(), times.size(), arma::fill::zeros);
    #pragma omp for schedule(dynamic)
        for(long t_idx = 0; t_idx < (long)times.size(); t_idx++)
        {
            const double time = times(t_idx);
            init_state.col(t_idx).zeros();
            for(long alfa = 0; alfa < E.size(); alfa++)
            {
                auto state = V.col(alfa);
                init_state.col(t_idx) += std::exp(-1i * time * E(alfa)) * state * coeff(alfa);
            }
        }
        arma::vec spread_coeff(E.size(), arma::fill::zeros);
        for(long alfa = 0; alfa < E.size(); alfa++)
        {
            spread_coeff(alfa) = ( x(alfa) - x(0) )*( x(alfa) - x(0) ) + ( y(alfa) - y(0) )*( y(alfa) - y(0) );
        }
        // std::cout << "Agadahbak i" << std::endl;
        arma::vec _ndown_spreading(times.size(), arma::fill::zeros);
        arma::vec return_probability(times.size(), arma::fill::zeros);
        arma::vec SpinDownNumber(times.size(), arma::fill::zeros);
    #pragma omp for schedule(dynamic)
        for(long t_idx = 0; t_idx < (long)times.size(); t_idx++)
        {
            const double time = times(t_idx);
            return_probability(t_idx) = std::norm(init_state(0, t_idx));
            
            for(long alfa = 0; alfa < E.size(); alfa++)
            {
                _ndown_spreading(t_idx) += std::norm(init_state(alfa, t_idx)) * spread_coeff(alfa);
                if(alfa < NA)
                    SpinDownNumber(t_idx) += std::norm(init_state(alfa, t_idx));
            }
        }

        return_probability.save(arma::hdf5_name(dir + info + ".hdf5", "return_probability", arma::hdf5_opts::append));
        _ndown_spreading.save(arma::hdf5_name(dir + info + ".hdf5", "spreading", arma::hdf5_opts::append));
        SpinDownNumber.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDownNumber(t)", arma::hdf5_opts::append));

        std::cout << " - - - - - - finished particle spreading in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
        std::cout << " - - - - - - finished realisation = " << real << " in : " << tim_s(start_re) << " s - - - - - - " << std::endl; // simulation end
    }
}

void ui::spectral_form_factor_shape(){
    std::string dir_out = this->dir_prefix + "results" + kPSep + "LiebLatticeRandomShape_SFF" + kPSep;
	createDirs(dir_out);
    std::string info = "_L=" + std::to_string(this->L) + \
                    ",tau=" + to_string_prec(this->tau) + \
                    ",m=" + to_string_prec(this->mass);
	
	//------- PREAMBLE
	// std::string info = this->set_info();

	u64 dim = 3 * this->L; // approx dimension

	double r1 = 0.0, r2 = 0.0;

	arma::vec times = arma::logspace(std::log10(1.0 / (two_pi * dim)), 1, this->num_of_points);
	arma::vec times_fold(this->num_of_points, arma::fill::zeros);
	
	arma::vec etas = arma::regspace(0.1, 0.1, 0.7);

	arma::vec energy_densities = arma::regspace(0.02, 0.02, 0.98);
	arma::vec gap_ratio(energy_densities.size(), arma::fill::zeros);
	arma::vec wH_density(energy_densities.size(), arma::fill::zeros);
	arma::vec wH_typ_density(energy_densities.size(), arma::fill::zeros);


	arma::mat sff_eps(energy_densities.size(), this->num_of_points, arma::fill::zeros);
	arma::mat sff_eps_folded(energy_densities.size(), this->num_of_points, arma::fill::zeros);
	arma::vec Z_eps(energy_densities.size(), arma::fill::zeros);
	arma::vec Z_eps_folded(energy_densities.size(), arma::fill::zeros);
	
	arma::vec sff_raw(this->num_of_points, arma::fill::zeros);
	arma::vec sff_raw_folded(this->num_of_points, arma::fill::zeros);
	double Z_raw = 0.0, Z_raw_folded = 0.0;

	arma::mat sff_eta(etas.size(), this->num_of_points, arma::fill::zeros);
	arma::vec Z_eta(etas.size(), arma::fill::zeros);
	
    arma::mat sff_eta_folded(etas.size(), this->num_of_points, arma::fill::zeros);
	arma::vec Z_eta_folded(etas.size(), arma::fill::zeros);

	double wH_mean = 0.0;
	double wH_typ  = 0.0;

    double wH_mean_dimer = 0., wH_typ_dimer = 0., wH_mean_mon = 0., wH_typ_mon = 0.;
    double r_dimer = 0., r_mon = 0.;

	arma::vec sff_dimer(this->num_of_points, arma::fill::zeros);
	arma::vec sff_dimer_folded(this->num_of_points, arma::fill::zeros);
	double Z_dimer = 0.0, Z_dimer_folded = 0.0;
	
    arma::vec sff_mon(this->num_of_points, arma::fill::zeros);
	arma::vec sff_mon_folded(this->num_of_points, arma::fill::zeros);
	double Z_mon = 0.0, Z_mon_folded = 0.0;

	u64 counter = 0;
	bool set_folded_times = 0;
// #pragma omp parallel for num_threads(outer_threads) // schedule(dynamic)
	for(int realis = 0; realis < this->l_realis; realis++)
	{
		clk::time_point start = std::chrono::system_clock::now();
		auto start_re = start;
	
		std::string prefix = "realisation=" + std::to_string(realis + this->jobid) + kPSep;
		// if(realis > 0)
		// 	this->ptr_to_model->generate_hamiltonian();
		arma::vec eigenvalues;// = this->get_eigenvalues(prefix, false);
		std::string dir = this->dir_prefix + "results" + kPSep + "LiebLatticeRandomShape" + kPSep + "realization=" + std::to_string(realis) + kPSep;
        // createDirs(dir);
        std::string name = dir + info;
        bool loaded = eigenvalues.load(arma::hdf5_name(name + ".hdf5", "E"));
		if( !loaded ){
            std::cout << "Not loaded:\t" << name << std::endl;
            continue;
        }
        arma::uvec Enon0_indices = arma::find(arma::abs(eigenvalues) > 1e-10);
        eigenvalues = eigenvalues(Enon0_indices);
        arma::vec E_save = eigenvalues;

		if(this->fun == 7) std::cout << "\t\t	--> finished loading eigenvalues for " << prefix + info << " - in time : " << tim_s(start) << "s" << std::endl;
		if(eigenvalues.empty()) continue;
		dim = eigenvalues.size();
		const double dE = eigenvalues(dim - 1) - eigenvalues(0);
		const double E0 = eigenvalues(0);
		const arma::vec gaps = arma::diff(eigenvalues);
		arma::vec _gap_ratios_ = arma::min( gaps.rows(0,dim-3), gaps.rows(1,dim-2)) / arma::max( gaps.rows(0,dim-3), gaps.rows(1,dim-2));

		start = std::chrono::system_clock::now();

		u64 E_av_idx = spectrals::get_mean_energy_index(eigenvalues);
		const u64 num = dim / 2;
		const u64 num2 = std::min( u64(500), dim/10);
		printSeparated(std::cout, "\t", 16, true, E_av_idx, arma::trace(eigenvalues) / double(eigenvalues.size()), num, num2);
		// ------------------------------------- calculate level statistics
			double r1_tmp = 0, r2_tmp = 0, wH_mean_r = 0, wH_typ_r = 0;
			int count = 0;
			for(int i = 1; i < dim - 1; i++){
				const double gap1 = eigenvalues(i) - eigenvalues(i - 1);
				const double gap2 = eigenvalues(i + 1) - eigenvalues(i);
				const double min = std::min(gap1, gap2);
				const double max = std::max(gap1, gap2);
        		if (abs(gap1) <= 1e-15 || abs(gap2) <= 1e-15){ 
                    printSeparated(std::cout, "\t", 20, true, "Found degeneracy", i, eigenvalues(i-1), eigenvalues(i), eigenvalues(i+1));
        		    // _assert_(false, "Found degeneracy, while doing r-statistics!\n");
                    continue;
        		}
				wH_mean_r += gap2;
				wH_typ_r += std::log(gap2);
				if(i >= (E_av_idx - num / 2) && i < (E_av_idx + num / 2))
					r1_tmp += min / max;
				if(i >= (E_av_idx - num2 / 2) && i < (E_av_idx + num2 / 2))
					r2_tmp += min / max;
				count++;
			}
			if(this->fun == 7) std::cout << "\t\t	--> finished unfolding for " << prefix + info << " - in time : " << tim_s(start) << "s" << std::endl;
			start = std::chrono::system_clock::now();

			wH_mean_r /= double(count);
			r1_tmp /= double(num);
			r2_tmp /= double(num2);

			#pragma omp critical
			{
				r1 += r1_tmp;
				r2 += r2_tmp;

				wH_mean += wH_mean_r;
				wH_typ  += wH_typ_r / double(count);
				counter++;
			}
			wH_typ_r = std::exp(wH_typ_r / double(count));

			if( !set_folded_times )
			{
				set_folded_times = 1;
				double tH_typ_temp = two_pi / std::exp(wH_typ_r / count);
				double time_end = std::log10(tH_typ_temp);
				time_end = (time_end / std::log10(tH_typ_temp) < 1.5) ? time_end + 0.8 : time_end;
				printSeparated(std::cout, "\t", 16, true, tH_typ_temp, time_end);
				times_fold = arma::logspace(-2.0, time_end, this->num_of_points);
			}


			arma::vec eigenvalues_unfolded = statistics::unfolding(eigenvalues).rows(5, dim - 5);

		// ------------------------------------- calculate sff
			statistics::SFF<statistics::filters::raw, statistics::ensemble::GC> SFF_raw(1.0, 0.0);
			arma::cx_vec x = SFF_raw.calculate(eigenvalues, times_fold);
            #pragma omp critical
			{
			sff_raw_folded += arma::abs(x) % arma::abs(x);
			Z_raw_folded += std::get<0>(SFF_raw.get_norms());
            }

			statistics::SFF<statistics::filters::raw, statistics::ensemble::GC> SFF_raw2(1.0, 0.0);
			x = SFF_raw2.calculate(eigenvalues_unfolded, times * two_pi);
            #pragma omp critical
			{
			sff_raw += arma::abs(x) % arma::abs(x);
			Z_raw += std::get<0>(SFF_raw2.get_norms());
            }

			if(this->fun == 7) std::cout << "\t\t	--> finished raw SFF for " << prefix + info << " - in time : " << tim_s(start) << "s" << std::endl;
			start = std::chrono::system_clock::now();

		// #pragma omp parallel for
			for(int e_idx = 0; e_idx < etas.size(); e_idx++)
			{
				statistics::SFF<statistics::filters::gauss, statistics::ensemble::GC> SFF_eta( etas(e_idx) );
				x = SFF_eta.calculate(eigenvalues_unfolded, times);
                #pragma omp critical
			    {
				sff_eta.row(e_idx) += (arma::abs(x) % arma::abs(x)).t();
				Z_eta(e_idx) += std::get<0>(SFF_eta.get_norms());
                }

				statistics::SFF<statistics::filters::gauss, statistics::ensemble::GC> SFF_eta2( etas(e_idx) );
				x = SFF_eta2.calculate(eigenvalues, times_fold);
                #pragma omp critical
			    {
				sff_eta_folded.row(e_idx) += (arma::abs(x) % arma::abs(x)).t();
				Z_eta_folded(e_idx) += std::get<0>(SFF_eta2.get_norms());
                }
			}

			if(this->fun == 7) std::cout << "\t\t	--> finished filtered SFF for " << prefix + info << " - in time : " << tim_s(start) << "s" << std::endl;
			start = std::chrono::system_clock::now();

		// #pragma omp parallel for
			for(int e_idx = 0; e_idx < energy_densities.size(); e_idx++)
			{
				double epsilon = energy_densities(e_idx);
				statistics::SFF<statistics::filters::gauss, statistics::ensemble::MC> SFF_MC(0.1, 0.0, epsilon);
				x = SFF_MC.calculate(eigenvalues_unfolded, times);
                #pragma omp critical
			    {
				sff_eps.row(e_idx) += (arma::abs(x) % arma::abs(x)).t();
				Z_eps(e_idx) += std::get<0>(SFF_MC.get_norms());
                }
                
				statistics::SFF<statistics::filters::gauss, statistics::ensemble::MC> SFF_MC2(0.1, 0.0, epsilon);
				x = SFF_MC2.calculate(eigenvalues, times_fold);
                #pragma omp critical
			    {
				sff_eps_folded.row(e_idx) += (arma::abs(x) % arma::abs(x)).t();
				Z_eps_folded(e_idx) += std::get<0>(SFF_MC2.get_norms());
                }

				u64 idx1 = arma::uvec( arma::sort_index( arma::vec(arma::abs( (eigenvalues - E0) / dE - (epsilon - 0.01))) ))(0);
				u64 idx2 = arma::uvec( arma::sort_index( arma::vec(arma::abs( (eigenvalues - E0) / dE - (epsilon + 0.01))) ))(0);
				if(idx2 - idx1 < 2){
					idx1 -= 5;
					idx2 += 5;
				}
				if(idx1 < 0 || idx1 > dim) idx1 = 0;
				if(idx2 > dim-3) idx2 = dim-3;
				// printSeparated(std::cout, "\t", 12, true, "IDX:=", idx1, idx2, dim);

                #pragma omp critical
			    {
                gap_ratio(e_idx) += arma::mean( _gap_ratios_.rows(idx1, idx2) );
				wH_density(e_idx) += arma::mean( gaps.rows(idx1, idx2) );
				wH_typ_density(e_idx) += std::exp( arma::mean( arma::log(gaps.rows(idx1, idx2) )) );
                }
			}

			if(this->fun == 7) std::cout << "\t\t	--> finished microcanonical SFF for " << prefix + info << " - in time : " << tim_s(start) << "s" << std::endl;
			start = std::chrono::system_clock::now();
            eigenvalues = E_save;
			arma::uvec Ediemr_ind = arma::find( eigenvalues < -1e-10 * this->mass);
            eigenvalues = eigenvalues(Ediemr_ind);
			r1_tmp = 0, wH_mean_r = 0, wH_typ_r = 0;
			count = 0;
			for(int i = 1; i < eigenvalues.size() - 1; i++){
				const double gap1 = eigenvalues(i) - eigenvalues(i - 1);
				const double gap2 = eigenvalues(i + 1) - eigenvalues(i);
				const double min = std::min(gap1, gap2);
				const double max = std::max(gap1, gap2);
        		if (abs(gap1) <= 1e-15 || abs(gap2) <= 1e-15){
                    printSeparated(std::cout, "\t", 20, true, "Found degeneracy", i, eigenvalues(i-1), eigenvalues(i), eigenvalues(i+1));
        		    // _assert_(false, "Found degeneracy, while doing r-statistics!\n");
                    continue;
        		}
				wH_mean_r += gap2;
				wH_typ_r += std::log(gap2);
				r1_tmp += min / max;
				count++;
			}

			#pragma omp critical
			{
            r_dimer += r1_tmp / double(count);
            wH_mean_dimer += wH_mean_r / double(count);
            wH_typ_dimer += wH_typ_r / double(count);
            }

            statistics::SFF<statistics::filters::raw, statistics::ensemble::GC> SFF_dimer(1.0, 0.0);
			x = SFF_dimer.calculate(eigenvalues, times_fold);
            #pragma omp critical
			{
			sff_dimer_folded += arma::abs(x) % arma::abs(x);
			Z_dimer_folded += std::get<0>(SFF_dimer.get_norms());
            }
			
            eigenvalues_unfolded = statistics::unfolding(eigenvalues).rows(5, eigenvalues.size() - 5);
			statistics::SFF<statistics::filters::raw, statistics::ensemble::GC> SFF_dimer2(1.0, 0.0);
			x = SFF_dimer2.calculate(eigenvalues_unfolded, times * two_pi);
            #pragma omp critical
			{
			sff_dimer += arma::abs(x) % arma::abs(x);
			Z_dimer += std::get<0>(SFF_dimer2.get_norms());
            }

			if(this->fun == 7) std::cout << "\t\t	--> finished dimer band SFF for " << prefix + info << " - in time : " << tim_s(start) << "s" << std::endl;
			start = std::chrono::system_clock::now();
            eigenvalues = E_save;
			arma::uvec Emon_ind = arma::find( arma::abs(eigenvalues) > 4 * this->mass);
            eigenvalues = eigenvalues(Emon_ind);
			r1_tmp = 0, wH_mean_r = 0, wH_typ_r = 0;
			count = 0;
			for(int i = 1; i < eigenvalues.size() - 1; i++){
				const double gap1 = eigenvalues(i) - eigenvalues(i - 1);
				const double gap2 = eigenvalues(i + 1) - eigenvalues(i);
				const double min = std::min(gap1, gap2);
				const double max = std::max(gap1, gap2);
        		if (abs(gap1) <= 1e-15 || abs(gap2) <= 1e-15){ 
                    printSeparated(std::cout, "\t", 20, true, "Found degeneracy", i, eigenvalues(i-1), eigenvalues(i), eigenvalues(i+1));
        		    // _assert_(false, "Found degeneracy, while doing r-statistics!\n");
                    continue;
        		}
				wH_mean_r += gap2;
				wH_typ_r += std::log(gap2);
				r1_tmp += min / max;
				count++;
			}

			#pragma omp critical
			{
            r_mon += r1_tmp / double(count);
            wH_mean_mon += wH_mean_r / double(count);
            wH_typ_mon += wH_typ_r / double(count);
            }
            statistics::SFF<statistics::filters::raw, statistics::ensemble::GC> SFF_mon(1.0, 0.0);
			x = SFF_mon.calculate(eigenvalues, times_fold);
            #pragma omp critical
			{
			sff_mon_folded += arma::abs(x) % arma::abs(x);
			Z_mon_folded += std::get<0>(SFF_mon.get_norms());
            }
            eigenvalues_unfolded = statistics::unfolding(eigenvalues).rows(5, eigenvalues.size() - 5);
			statistics::SFF<statistics::filters::raw, statistics::ensemble::GC> SFF_mon2(1.0, 0.0);
			x = SFF_mon2.calculate(eigenvalues_unfolded, times * two_pi);
            #pragma omp critical
			{
			sff_mon += arma::abs(x) % arma::abs(x);
			Z_mon += std::get<0>(SFF_mon2.get_norms());
            }
			if(this->fun == 7) std::cout << "\t\t	--> finished monomer band SFF for " << prefix + info << " - in time : " << tim_s(start) << "s" << std::endl;
			start = std::chrono::system_clock::now();

		if(this->fun == 7) std::cout << "--> finished realisation for " << prefix + info << " - in time : " << tim_s(start_re) << "s" << std::endl << std::endl;
	}

	// --------------------------------------------------------------- AVERAGE CURRENT REALISATIONS
	if(sff_eps.is_empty()) return;
	if(sff_eps.is_zero()) return;
	if(this->jobid > 0) return;
	if(counter == 0) return;

	double norm = counter;
	r1 /= norm;
	r2 /= norm;
    r_dimer /= norm;
    r_mon /= norm;

	gap_ratio /= norm;
	wH_density /= norm;
	wH_typ_density /= norm;

	for(int e_idx = 0; e_idx < energy_densities.size(); e_idx++){
		sff_eps.row(e_idx) = sff_eps.row(e_idx) / Z_eps(e_idx);
		sff_eps_folded.row(e_idx) = sff_eps_folded.row(e_idx) / Z_eps_folded(e_idx);
	}
	for(int e_idx = 0; e_idx < etas.size(); e_idx++){
		sff_eta.row(e_idx) = sff_eta.row(e_idx) / Z_eta(e_idx);
		sff_eta_folded.row(e_idx) = sff_eta_folded.row(e_idx) / Z_eta_folded(e_idx);
	}

	sff_raw = sff_raw / Z_raw;
	sff_raw_folded = sff_raw_folded / Z_raw_folded;

	sff_mon = sff_mon / Z_mon;
	sff_mon_folded = sff_mon_folded / Z_mon_folded;
	sff_dimer = sff_dimer / Z_dimer;
	sff_dimer_folded = sff_dimer_folded / Z_dimer_folded;

	wH_mean /= norm;
	wH_typ /= norm;

	wH_mean_dimer /= norm;
	wH_typ_dimer /= norm;
	wH_mean_mon /= norm;
	wH_typ_mon /= norm;

	// #ifdef MY_MAC
	times.save(arma::hdf5_name(dir_out + info + ".hdf5", "times"));
	times_fold.save(arma::hdf5_name(dir_out + info + ".hdf5", "times_fold", arma::hdf5_opts::append));

	energy_densities.save(arma::hdf5_name(dir_out + info + ".hdf5", "energy_densities", arma::hdf5_opts::append));
	sff_eps.save(arma::hdf5_name(dir_out + info + ".hdf5", "sff_eps", arma::hdf5_opts::append));
	sff_eps_folded.save(arma::hdf5_name(dir_out + info + ".hdf5", "sff_eps_folded", arma::hdf5_opts::append));
	gap_ratio.save(arma::hdf5_name(dir_out + info + ".hdf5", "gap_ratio density", arma::hdf5_opts::append));
	wH_density.save(arma::hdf5_name(dir_out + info + ".hdf5", "wH density", arma::hdf5_opts::append));
	wH_typ_density.save(arma::hdf5_name(dir_out + info + ".hdf5", "wH_typ density", arma::hdf5_opts::append));

	etas.save(arma::hdf5_name(dir_out + info + ".hdf5", "etas", arma::hdf5_opts::append));
	sff_eta.save(arma::hdf5_name(dir_out + info + ".hdf5", "sff_eta", arma::hdf5_opts::append));
	sff_eta_folded.save(arma::hdf5_name(dir_out + info + ".hdf5", "sff_eta_folded", arma::hdf5_opts::append));

	sff_raw.save(arma::hdf5_name(dir_out + info + ".hdf5", "sff_raw", arma::hdf5_opts::append));
	sff_raw_folded.save(arma::hdf5_name(dir_out + info + ".hdf5", "sff_raw_folded", arma::hdf5_opts::append));

    sff_mon.save(arma::hdf5_name(dir_out + info + ".hdf5", "sff_mon", arma::hdf5_opts::append));
	sff_mon_folded.save(arma::hdf5_name(dir_out + info + ".hdf5", "sff_mon_folded", arma::hdf5_opts::append));
	sff_dimer.save(arma::hdf5_name(dir_out + info + ".hdf5", "sff_dimer", arma::hdf5_opts::append));
	sff_dimer_folded.save(arma::hdf5_name(dir_out + info + ".hdf5", "sff_dimer_folded", arma::hdf5_opts::append));

	arma::vec({r2}).save(arma::hdf5_name(dir_out + info + ".hdf5", "r_500", arma::hdf5_opts::append));
	arma::vec({r1}).save(arma::hdf5_name(dir_out + info + ".hdf5", "r_D_2", arma::hdf5_opts::append));
    arma::vec({r_dimer}).save(arma::hdf5_name(dir_out + info + ".hdf5", "r_dimer", arma::hdf5_opts::append));
    arma::vec({r_mon}).save(arma::hdf5_name(dir_out + info + ".hdf5", "r_mon", arma::hdf5_opts::append));
	arma::uvec({dim}).save(arma::hdf5_name(dir_out + info + ".hdf5", "D", arma::hdf5_opts::append));
	arma::vec({two_pi / (wH_mean)}).save(arma::hdf5_name(dir_out + info + ".hdf5", "tH", arma::hdf5_opts::append));
	arma::vec({two_pi / std::exp(wH_typ)}).save(arma::hdf5_name(dir_out + info + ".hdf5", "tH_typ", arma::hdf5_opts::append));

	arma::vec({two_pi / (wH_mean_dimer)}).save(arma::hdf5_name(dir_out + info + ".hdf5", "tH_dimer", arma::hdf5_opts::append));
	arma::vec({two_pi / std::exp(wH_typ_dimer)}).save(arma::hdf5_name(dir_out + info + ".hdf5", "tH_typ_dimer", arma::hdf5_opts::append));
	arma::vec({two_pi / (wH_mean_mon)}).save(arma::hdf5_name(dir_out + info + ".hdf5", "tH_mon", arma::hdf5_opts::append));
	arma::vec({two_pi / std::exp(wH_typ_mon)}).save(arma::hdf5_name(dir_out + info + ".hdf5", "tH_typ_mon", arma::hdf5_opts::append));
	arma::uvec({counter}).save(arma::hdf5_name(dir_out + info + ".hdf5", "realisations", arma::hdf5_opts::append));
}

void ui::LiebLatticeStrip()
{
    // std::string dir = this->dir_prefix + "results" + kPSep + "LiebLatticeStrip" + kPSep + "Straight" + kPSep;
    // std::string dir = this->dir_prefix + "results" + kPSep + "LiebLatticeStrip" + kPSep + "Bearded" + kPSep;
    std::string dir = this->dir_prefix + "results" + kPSep + "LiebLatticeStrip" + kPSep + "Assymetric" + kPSep;
	createDirs(dir);
	std::string info = "_Lx=" + std::to_string(this->Lx) + \
                    ",Ly=" + std::to_string(this->Ly) + \
                    ",tau=" + to_string_prec(this->tau) + \
                    ",m=" + to_string_prec(this->mass);

    auto start = std::chrono::system_clock::now();

	arma::vec kx_points = arma::linspace(0, this->Lx, this->Lx + 1) * constants<double>::two_pi / double(this->Lx);
    
    arma::mat energies(3 * this->Ly, kx_points.size(), arma::fill::zeros); 
    arma::mat central_states(3 * this->Ly, kx_points.size(), arma::fill::zeros);
    for(int ik = 0; ik < kx_points.size(); ik++)
    {
        auto start_loop = std::chrono::system_clock::now();
        double kx = kx_points(ik);
        cpx fx = -this->tau * (1.0 + std::exp(1i * kx));
        // arma::cx_mat H(3 * this->Ly + 1, 3 * this->Ly + 1, arma::fill::zeros);
        // for(int i = 0; i < this->Ly; i++){
        //     H(3*i + 1, 3*i + 1) = 4 * this->mass;
            
        //     H(3*i, 3*i + 1) = -this->tau;
        //     H(3*i + 1, 3*i) = -this->tau;
        //     {
        //         H(3*i+1, 3*i+2) = fx;
        //         H(3*i+2, 3*i+1) = std::conj(fx);

        //         //<! coupling to next unit cell
        //         H(3*i+1, 3*i+3) = -this->tau;
        //         H(3*i+3, 3*i+1) = -this->tau;
        //     }
        // }
        // arma::cx_mat H(3 * this->Ly+2, 3 * this->Ly+2, arma::fill::zeros);
        // for(int i = 0; i < this->Ly; i++){
        //     H(3*i, 3*i) = 4 * this->mass;
            
        //     H(3*i, 3*i + 1) = fx;
        //     H(3*i + 1, 3*i) = std::conj(fx);

        //     H(3*i, 3*i + 2) = -this->tau;
        //     H(3*i + 2, 3*i) = -this->tau;

        //     //<! coupling to next unit cell
        //     H(3*i + 2, 3*i + 3) = -this->tau;
        //     H(3*i + 3, 3*i + 2) = -this->tau;
        // }
        // H(3*this->Ly, 3*this->Ly) = 4*this->mass;
        // H(3*this->Ly, 3*this->Ly+1) = fx;
        // H(3*this->Ly+1, 3*this->Ly) = std::conj(fx);


        arma::cx_mat H(3 * this->Ly, 3 * this->Ly, arma::fill::zeros);
        for(int i = 0; i < this->Ly; i++){
            if(i == 0)  H(3*i, 3*i) = 3 * this->mass;
            else        H(3*i, 3*i) = 4 * this->mass;
            
            H(3*i, 3*i + 1) = fx;
            H(3*i + 1, 3*i) = std::conj(fx);

            H(3*i, 3*i + 2) = -this->tau;
            H(3*i + 2, 3*i) = -this->tau;

            //<! coupling to next unit cell
            if(i < this->Ly - 1){
                H(3*i + 2, 3*i + 3) = -this->tau;
                H(3*i + 3, 3*i + 2) = -this->tau;
            }
        }
        // std::cout << H << std::endl;

        arma::vec E; arma::cx_mat V;
        arma::eig_sym(E, V, H);

        energies.col(ik) = E;
        central_states.col(ik) = arma::abs(V.col(3 * this->Ly / 2));

        std::cout << " - - - - - - finished diagonalization of matrix with dimension D = " << E.size() << " in : " << tim_s(start_loop) << " s - - - - - - " << std::endl; // simulation end
    }
    energies.save(arma::hdf5_name(dir + info + ".hdf5", "E"));
    kx_points.save(arma::hdf5_name(dir + info + ".hdf5", "kx_points", arma::hdf5_opts::append));
    central_states.save(arma::hdf5_name(dir + info + ".hdf5", "central_states", arma::hdf5_opts::append));
    std::cout << " - - - - - - finished diagonalization of strip in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
}

/// @brief Compare energie spactra for full model and all symmetry sectors combined
void ui::spin_impurity()
{          
    // #define DIMER_SPREADING
    
    #ifdef DIMER_SPREADING
        std::string dir = this->saving_dir + "DimerSpreading" + kPSep;
    #else
        std::string dir = this->saving_dir + "SpinImpurity" + kPSep;
    #endif
	createDirs(dir);
	
	size_t dim = this->ptr_to_model->get_hilbert_size();
	std::string info = this->set_info();

    size_t dim_cut = 1e5;
    arma::SpMat<ui::element_type> H = this->ptr_to_model->get_hamiltonian();

	const size_t size = dim > dim_cut? this->l_steps : dim;

    auto _hilbert_space = this->ptr_to_model->get_model_ref().get_hilbert_space();
    auto _group = _hilbert_space.get_symmetry_group();

	clk::time_point start = std::chrono::system_clock::now();
    clk::time_point start_tot = start;
    if(dim > dim_cut){
        this->ptr_to_model->diag_sparse(this->l_steps, this->l_bundle, this->tol, this->seed);	
    }
    else{
        this->ptr_to_model->diagonalization();
    }
    const double Emin = this->ptr_to_model->get_eigenValue(0);
    const double Emax = this->ptr_to_model->get_eigenValue(dim-1);
    std::cout << " - - - - - - finished diagonalization in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();
    
	QHS::Symmetric2DLattice<> _hilbert_full(this->Lx, this->Ly, this->syms.Nup, this->syms.Ndown, v_1d<QOps::generic_operator<cpx, PackedLattice2D>>(), this->boundary_conditions, 0, 0, 1, 1, true, false);
    u64 dimtot = _hilbert_full.get_hilbert_space_size();
    
    std::cout << " - - - - - - finished generating total Hilbert space of size D = " << dimtot << " in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();
    // return;
    arma::sp_mat n_dimer_center(dimtot, dimtot);
    arma::sp_mat n_spindown_center(dimtot, dimtot);
    arma::sp_mat DimerNumber(dimtot, dimtot);
    arma::sp_mat SpinDownNumber(dimtot, dimtot);
    
    // arma::sp_mat DimerEdgeNumber(dimtot, dimtot);
    arma::sp_mat SpinDownEdgeNumber(dimtot, dimtot);
    for(long alfa = 0; alfa < dimtot; alfa++)
    {
        auto base_state = _hilbert_full(alfa);
        auto _conf1 = base_state.get_block(this->Lx/2, this->Ly/2);
        n_spindown_center(alfa, alfa) += _fermion_down_number(_conf1);
        if(_conf1 == DIMER_RIGHT)
            n_dimer_center(alfa, alfa) += 1.0;
        for(std::size_t y = 0; y < this->Ly; y++){
            for(std::size_t x = 0; x < this->Lx; x++){
                auto _conf = base_state.get_block(x, y);
                SpinDownNumber(alfa, alfa) += _fermion_down_number(_conf);
                if(_conf == DIMER_RIGHT) 
                    DimerNumber(alfa, alfa) += 1.0;
                if(_conf == DIMER_UP) 
                    DimerNumber(alfa, alfa) += 1.0;
            }
        }
        for(std::size_t y = 0; y < this->Ly; y++){
            auto _conf = base_state.get_block(0, y);
            if(_conf == DOWN_MONOMER) 
                SpinDownEdgeNumber(alfa, alfa) += 1.0;
            
            _conf = base_state.get_block(this->Lx-1, y);
            if(_conf == DOWN_MONOMER) 
                SpinDownEdgeNumber(alfa, alfa) += 1.0;
        }
        for(std::size_t x = 0; x < this->Lx; x++){
            auto _conf = base_state.get_block(x, 0);
            if(_conf == DOWN_MONOMER) 
                SpinDownEdgeNumber(alfa, alfa) += 1.0;
            
            _conf = base_state.get_block(x, this->Ly-1);
            if(_conf == DOWN_MONOMER) 
                SpinDownEdgeNumber(alfa, alfa) += 1.0;
        }
    }
    #ifdef USE_SYMMETRIES
        auto U = _hilbert_space.symmetry_rotation(_hilbert_full);
    #endif
    std::cout << " - - - - - - finished setting up operators in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();

    PackedLattice2D initial_state(this->Lx, this->Ly);
    // initial_state.fill_all(EMPTY_SITE);
    // if(this->syms.Nup + this->syms.Ndown == this->L)
    initial_state.fill_all(UP_MONOMER);
    // int N = std::min(this->)
    #ifdef DIMER_SPREADING
        int x0 = this->Lx/2-1;
        int y0 = this->Ly/2;
        initial_state.set_block(x0, y0, DIMER_RIGHT);
        initial_state.set_block(x0+1, y0, DIMER_LEFT);
        if(this->syms.Ndown == 2){
            initial_state.set_block(x0, y0+1, DOWN_MONOMER);
        }
    #else
        int x0 = this->Lx/2 - int(Lx > Ly && this->syms.Ndown <= 4);
        int y0 = this->Ly/2;
        int nd = 1;

        // initial_state.set_block(0, 1, DOWN_MONOMER);
        // initial_state.set_block(1, 1, DOWN_MONOMER);
        initial_state.set_block(x0, y0, DOWN_MONOMER);
        if(this->syms.Ndown > 1){
            initial_state.set_block(x0+1, y0, DOWN_MONOMER);
        }
        if(this->syms.Ndown > 2){
            initial_state.set_block(x0, y0+1, DOWN_MONOMER);
        }
        if(this->syms.Ndown > 3){
            initial_state.set_block(x0, y0-1, DOWN_MONOMER);
        }
        if(this->syms.Ndown > 4){
            initial_state.set_block(x0+1, y0, DOWN_MONOMER);
        }
        // if(this->syms.Nup == 1)
        //     initial_state.set_block(x0-1, y0, UP_MONOMER);
        // for(int xp : {-1, 1}){
        //     if(nd >= this->syms.Ndown) continue;
        //     initial_state.set_block(x0 + xp, y0, DOWN_MONOMER);
        //     nd++;
        // }
        // for(int yp : {-1, 1}){
        //     if(nd >= this->syms.Ndown) continue;
        //     initial_state.set_block(x0, y0 + yp, DOWN_MONOMER);
        //     nd++;
        // }
    #endif

    auto initial_state_matrix = initial_state._lattice;
    initial_state_matrix.save(arma::hdf5_name(dir + info + ".hdf5", "InitialState"));

    u64 index = _hilbert_full.find(initial_state);
    initial_state.print_state(); std::cout << std::endl;
    // _hilbert_full(index).print_state(); std::cout << std::endl;
    
    std::cout << "\n" << index << std::endl << std::endl;
    arma::vec energies = arma::linspace(Emin - 2, Emax + 2, 3000);
    arma::vec ldos(energies.size(), arma::fill::zeros);
    arma::vec DOS(energies.size(), arma::fill::zeros);
    double _eta = 0.05;

    arma::vec times = this->syms.Ndown>1? arma::logspace(-3, 7, 1000) : arma::logspace(-3, std::log10(10*this->L) + 1, 1000);
    arma::vec times_append({0, 0.1, 0.2, 0.5, 1, 1.5, 2., 3., 4., 6., 8., 10., 20., 30., 40., 50., 75., 100., 150., 200.});
    times = arma::sort(arma::join_cols(times, times_append));
    times.save(arma::hdf5_name(dir + info + ".hdf5", "times", arma::hdf5_opts::append));
    arma::Col<u64>({dimtot}).save(arma::hdf5_name(dir + info + ".hdf5", "DIMTOT", arma::hdf5_opts::append));
    
    energies.save(arma::hdf5_name(dir + info + ".hdf5", "ldos_energies", arma::hdf5_opts::append));
    arma::cx_mat init_state(dimtot, times.size(), arma::fill::zeros);
    #ifdef USE_SYMMETRIES
        // arma::sp_cx_mat iden(dimtot, dimtot);

        //<! Collect coefficients and energies to calculate DOS and LDOS to get exact energy bounds for the calculations
        auto kernel = [&](int kx, int ky, int z, int px, int py, int pd)
        {
            SymmetrySectors sym2 = this->syms;
            sym2.kx_sym = kx;
            sym2.ky_sym = ky;
            sym2.z_sym = z;
            sym2.px_sym = px;
            sym2.py_sym = py;
            sym2.pd_sym = pd;
            auto symmetric_model = std::make_unique<QHS::QHamSolver<BosonicSpinfullQLM>>(this->boundary_conditions, this->Lx, this->Ly, this->tau, this->mass, sym2, true);
            auto _hilbert_space_sym = symmetric_model->get_model_ref().get_hilbert_space();
            u64 dim = symmetric_model->get_hilbert_size();
            if(dim > 0){
                std::string syminfo = "";
                if(this->boundary_conditions == 0){
                    syminfo += ",kx=" + std::to_string(sym2.kx_sym);
                    syminfo += ",ky=" + std::to_string(sym2.ky_sym);
                }
                #ifndef USE_ONLY_TRANSLATION
                    if(sym2.Nup == sym2.Ndown) 
                        syminfo += ",z=" + std::to_string(sym2.z_sym);
                    if(this->kx_real_sec(sym2.kx_sym)) 
                        syminfo += ",px=" + std::to_string(sym2.px_sym);
                    if(this->ky_real_sec(sym2.ky_sym)) 
                        syminfo += ",py=" + std::to_string(sym2.py_sym);
                    if( (this->Lx == this->Ly) && (sym2.kx_sym == sym2.ky_sym) && sym2.px_sym == sym2.py_sym )
                        syminfo += ",pd=" + std::to_string(sym2.pd_sym);
                #endif
                std::cout << syminfo << std::endl;
                symmetric_model->diagonalization(true);
                arma::vec Esym = symmetric_model->get_eigenvalues();
                Esym.save(arma::hdf5_name(dir + info + ".hdf5", "energies/" + syminfo , arma::hdf5_opts::append));

                auto [new_state, sym_eig] = _hilbert_space_sym.find_SEC_representative(initial_state);
                u64 indexsym = _hilbert_space_sym.find(new_state);

                if(indexsym >= 0 && indexsym < dim)
                {
                    auto Usym = _hilbert_space_sym.symmetry_rotation(_hilbert_full);
                    // arma::Row<ui::element_type> coeff = Usym.row(index) * symmetric_model->get_eigenvectors();
                    // std::cout << Usym.row(index) << std::endl;
                    arma::Mat<ui::element_type> Phi = Usym * symmetric_model->get_eigenvectors();
                    auto coeff = Phi.row(index);
                    double norm_inv = 1. / std::sqrt(constants<double>::two_pi * _eta*_eta);
                #pragma omp for schedule(dynamic)
                    for(long e_idx = 0; e_idx < energies.size(); e_idx++){
                        // doubl
                        // ldos(e_idx) += arma::sum(arma::square(coeff.t()) % _eta / ( arma::square(Esym))
                        for(long alfa = 0; alfa < Esym.size(); alfa++){
                            double om = Esym(alfa) - energies(e_idx);
                            double gauss = std::exp( - om * om / (2.*_eta * _eta) ) * norm_inv;
                            ldos(e_idx) += std::norm(coeff(alfa)) * gauss;
                            DOS(e_idx) += gauss;
                        }
                    }
                #pragma omp for schedule(dynamic)
                    for(long t_idx = 0; t_idx < (long)times.size(); t_idx++)
                    {
                        const double time = times(t_idx);
                        init_state.col(t_idx) += Phi * (coeff.t() % arma::exp(-1i*time*Esym));
                        // for(long alfa = 0; alfa < Esym.size(); alfa++)
                        // {
                        //     arma::Col<ui::element_type> state = Usym * symmetric_model->get_eigenState(alfa);
                        //     init_state.col(t_idx) += std::exp(-1i * time * Esym(alfa)) * state * ui::element_type( coeff(alfa) );
                        // }
                        // init_state.col(t_idx) = arma::normalise(init_state.col(t_idx));
                    }
                } else {
                    std::cout << "In sector: " + syminfo + " the initial state does not exist!" << std::endl;
                }
            }
        };
        loopSymmetrySectors(kernel);
    #else
        const arma::vec E = this->ptr_to_model->get_eigenvalues();
        // const auto& V = this->ptr_to_model->get_eigenvectors();
        auto coeff = this->ptr_to_model->get_eigenStateRow(index);
        E.save(arma::hdf5_name(dir + info + ".hdf5", "energies", arma::hdf5_opts::append));
        double norm_inv = 1. / std::sqrt(constants<double>::two_pi * _eta*_eta);
    #pragma omp for schedule(dynamic)
        for(long e_idx = 0; e_idx < energies.size(); e_idx++){
            // doubl
            // ldos(e_idx) += arma::sum(arma::square(coeff.t()) % _eta / ( arma::square(Esym))
            for(long alfa = 0; alfa < E.size(); alfa++){
                double om = E(alfa) - energies(e_idx);
                double gauss = std::exp( - om * om / (2.*_eta * _eta) ) * norm_inv;
                ldos(e_idx) += std::norm(coeff(alfa)) * gauss;
                DOS(e_idx) += gauss;
            }
        }
    #pragma omp for schedule(dynamic)
        for(long t_idx = 0; t_idx < (long)times.size(); t_idx++)
        {
            const double time = times(t_idx);
            init_state.col(t_idx).zeros();
            for(long alfa = 0; alfa < E.size(); alfa++)
            {
                auto state = this->ptr_to_model->get_eigenState(alfa);
                init_state.col(t_idx) += std::exp(-1i * time * E(alfa)) * state * coeff(alfa);
            }
            // init_state.col(t_idx) = arma::normalise(init_state.col(t_idx));
        }
    #endif
    for(int k = 0; k < dimtot; k++){
        if( std::norm(init_state.col(0)(k)) > 1e-12 ){
            std::cout << std::endl << init_state.col(0)(k) << std::endl;
            _hilbert_full(k).print_state();
            
        }
    }
    DOS.save(arma::hdf5_name(dir + info + ".hdf5", "DOS", arma::hdf5_opts::append));
    ldos.save(arma::hdf5_name(dir + info + ".hdf5", "LDOS", arma::hdf5_opts::append));
    std::cout << " - - - - - - finished setting up initial state in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();

    arma::vec _ndown_center(size, arma::fill::zeros);
    arma::vec _spindown_num(size, arma::fill::zeros);
    arma::vec _ndimer_center(size, arma::fill::zeros);
    arma::vec _dimer_num(size, arma::fill::zeros);
    arma::vec _ipr(size, arma::fill::zeros);
    arma::vec _edge_down_localization(size, arma::fill::zeros);
    
    #ifndef USE_SYMMETRIES
        arma::vec _overlap_init_state(E.size(), arma::fill::zeros);
        for(int k = 0; k < E.size(); k++)
        {
            double x = std::norm(coeff(k));
            _overlap_init_state(k) += x;
        }
        _overlap_init_state.save(arma::hdf5_name(dir + info + ".hdf5", "overlap_init_state", arma::hdf5_opts::append));
    #endif
// #pragma omp parallel for
#pragma omp parallel
    {
        // Per-thread preallocated buffer → no allocation inside loop
        arma::mat _nspindown_density_state(this->Ly, this->Lx, arma::fill::zeros);
        arma::mat _ndimer_right_density_state(this->Ly, this->Lx, arma::fill::zeros);
        arma::mat _ndimer_up_density_state(this->Ly, this->Lx, arma::fill::zeros);

    #pragma omp for schedule(dynamic)
        for(int k = 0; k < size; k++)
        {
            #ifdef USE_SYMMETRIES
                arma::Col<ui::element_type> state = U * this->ptr_to_model->get_eigenState(k);
            #else
                arma::Col<ui::element_type> state = this->ptr_to_model->get_eigenState(k);
            #endif
            // arma::mat _ndown_center_state(this->Ly, this->Lx, arma::fill::zeros);
            // arma::mat _nmonomer_density_state(this->Ly, this->Lx, arma::fill::zeros);
            double tmp = 0;
            for(long alfa = 0; alfa < dimtot; alfa++)
                tmp += std::norm(state(alfa)) * std::norm(state(alfa));
            
            _ipr(k) = tmp;
            _dimer_num(k) = std::real( arma::cdot(state, DimerNumber * state) );
            _ndimer_center(k) = std::real( arma::cdot(state, n_dimer_center * state) );
            _ndown_center(k) = std::real( arma::cdot(state, n_spindown_center * state) );
            _spindown_num(k) = std::real( arma::cdot(state, SpinDownNumber * state) );
            _edge_down_localization(k) = std::real( arma::cdot(state, SpinDownEdgeNumber * state) );
            
            _nspindown_density_state.zeros();
            _ndimer_right_density_state.zeros();
            _ndimer_up_density_state.zeros();
            for(long alfa = 0; alfa < dimtot; alfa++)
            {
                auto base_state = _hilbert_full(alfa);
                double weight = std::norm(state(alfa));
            // #pragma omp parallel for
                for(std::size_t y = 0; y < this->Ly; y++){
                    for(std::size_t x = 0; x < this->Lx; x++){
                        auto _conf = base_state.get_block(x, y);
                        if(_conf == DOWN_MONOMER) 
                            _nspindown_density_state(y, x) += weight;
                        if(_conf == DIMER_RIGHT) 
                            _ndimer_right_density_state(y, x) += weight;
                        if(_conf == DIMER_UP) 
                            _ndimer_up_density_state(y, x) += weight;
                    }
                }
            }
        #pragma omp critical
            {
            _nspindown_density_state.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDown_State/n = " + std::to_string(k), arma::hdf5_opts::append));
            _ndimer_right_density_state.save(arma::hdf5_name(dir + info + ".hdf5", "DimerRight_State/n = " + std::to_string(k), arma::hdf5_opts::append));
            _ndimer_up_density_state.save(arma::hdf5_name(dir + info + ".hdf5", "DimerUp_State/n = " + std::to_string(k), arma::hdf5_opts::append));
            }
        }
    }

    _edge_down_localization.save(arma::hdf5_name(dir + info + ".hdf5", "EdgeDownLoc", arma::hdf5_opts::append));
    _spindown_num.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDownNumber_n", arma::hdf5_opts::append));
    _ndimer_center.save(arma::hdf5_name(dir + info + ".hdf5", "DimerCenter_n", arma::hdf5_opts::append));
    _ndown_center.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDownCenter_n", arma::hdf5_opts::append));
    _dimer_num.save(arma::hdf5_name(dir + info + ".hdf5", "DimerNumber_n", arma::hdf5_opts::append));
    _ipr.save(arma::hdf5_name(dir + info + ".hdf5", "IPR", arma::hdf5_opts::append));
    std::cout << " - - - - - - finished eigenstate measures in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();

    _ndown_center = arma::vec(times.size(), arma::fill::zeros);
    _ndimer_center = arma::vec(times.size(), arma::fill::zeros);
    _dimer_num = arma::vec(times.size(), arma::fill::zeros);
    _spindown_num = arma::vec(times.size(), arma::fill::zeros);
    _edge_down_localization = arma::vec(times.size(), arma::fill::zeros);
    
    arma::vec _ndown_spreading(times.size(), arma::fill::zeros);

    arma::vec return_probability(times.size(), arma::fill::zeros);
    arma::vec spread_coeff(dimtot, arma::fill::zeros);
    for(long alfa = 0; alfa < dimtot; alfa++)
    {
        auto base_state = _hilbert_full(alfa);
        for(std::size_t y = 0; y < this->Ly; y++){
            for(std::size_t x = 0; x < this->Lx; x++){
                auto _conf = base_state.get_block(x, y);
                #ifdef DIMER_SPREADING
                    if(_conf == DIMER_RIGHT) 
                        spread_coeff(alfa) += ((x - x0)*(x - x0) + (y - y0)*(y - y0));
                    if(_conf == DIMER_UP) 
                        spread_coeff(alfa) += ((x - (x0+0.5))*(x - (x0+0.5)) + (y+0.5 - y0)*(y+0.5 - y0));
                #else
                    if(_conf == DOWN_MONOMER) 
                        spread_coeff(alfa) += ((x - x0)*(x - x0) + (y - y0)*(y - y0));
                    if(_conf == DIMER_RIGHT) 
                        spread_coeff(alfa) += ((x+0.5 - x0)*(x+0.5 - x0) + (y - y0)*(y - y0));
                    if(_conf == DIMER_UP) 
                        spread_coeff(alfa) += ((x - x0)*(x - x0) + (y+0.5 - y0)*(y+0.5 - y0));
                #endif
            }
        }
    }
#pragma omp parallel
    {
        // Per-thread preallocated buffer → no allocation inside loop
        arma::mat _ndown_density_t(this->Ly, this->Lx, arma::fill::zeros);
        arma::mat _nmonomer_density_t(this->Ly, this->Lx, arma::fill::zeros);
        arma::mat _dimer_right_density_t(this->Ly, this->Lx, arma::fill::zeros);

    #pragma omp for schedule(dynamic)
        for(long t_idx = 0; t_idx < (long)times.size(); t_idx++)
        {
            const double time = times(t_idx);
            return_probability(t_idx) = std::norm(init_state(index, t_idx));
            
            _edge_down_localization(t_idx) = std::real( arma::cdot(init_state.col(t_idx), SpinDownEdgeNumber * init_state.col(t_idx)) );

            _spindown_num(t_idx) = std::real( arma::cdot(init_state.col(t_idx), SpinDownNumber * init_state.col(t_idx)) );
            _dimer_num(t_idx) = std::real( arma::cdot(init_state.col(t_idx), DimerNumber * init_state.col(t_idx)) );
            _ndimer_center(t_idx) = std::real( arma::cdot(init_state.col(t_idx), n_dimer_center * init_state.col(t_idx)) );
            _ndown_center(t_idx) = std::real( arma::cdot(init_state.col(t_idx), n_spindown_center * init_state.col(t_idx)) );

            bool generate_real_space_slice = std::find(times_append.begin(), times_append.end(), time) != times_append.end();
            if(generate_real_space_slice)
                _ndown_density_t.zeros();
            // arma::mat _nmonomer_density_t(this->Ly, this->Lx, arma::fill::zeros);
            for(long alfa = 0; alfa < dimtot; alfa++)
            {
                auto base_state = _hilbert_full(alfa);
                _ndown_spreading(t_idx) += std::norm(init_state(alfa, t_idx)) * spread_coeff(alfa);
            }
            
            if (generate_real_space_slice) {
                _ndown_density_t.zeros();
                _nmonomer_density_t.zeros();
                _dimer_right_density_t.zeros();
                for(long alfa = 0; alfa < dimtot; alfa++)
                {
                    auto base_state = _hilbert_full(alfa);
                    double weight = std::norm(init_state(alfa, t_idx));// * std::abs(init_state(alfa));
                    for(std::size_t y = 0; y < this->Ly; y++){
                        for(std::size_t x = 0; x < this->Lx; x++){
                            auto _conf = base_state.get_block(x, y);
                            _ndown_density_t(y, x) += _fermion_down_number(_conf) * weight;
                            _nmonomer_density_t(y, x) += _fermion_number(_conf) * weight;
                            if(_conf == DIMER_RIGHT) 
                                _dimer_right_density_t(y, x) += weight;
                        }
                    }
                }
                #pragma omp critical
                {
                _ndown_density_t.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDownDensity_snapshot/t = " + to_string_prec(time), arma::hdf5_opts::append));
                _nmonomer_density_t.save(arma::hdf5_name(dir + info + ".hdf5", "MonomerDensity_snapshot/t = " + to_string_prec(time), arma::hdf5_opts::append));
                _dimer_right_density_t.save(arma::hdf5_name(dir + info + ".hdf5", "DimerRightDensity_snapshot/t = " + to_string_prec(time), arma::hdf5_opts::append));
                }
            }
        }
    }
    _ndimer_center.save(arma::hdf5_name(dir + info + ".hdf5", "DimerCenter(t)", arma::hdf5_opts::append));
    _ndown_center.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDownCenter(t)", arma::hdf5_opts::append));
    _dimer_num.save(arma::hdf5_name(dir + info + ".hdf5", "DimerNumber(t)", arma::hdf5_opts::append));
    _spindown_num.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDownNumber(t)", arma::hdf5_opts::append));
    _edge_down_localization.save(arma::hdf5_name(dir + info + ".hdf5", "EdgeDownLoc(t)", arma::hdf5_opts::append));
    
    return_probability.save(arma::hdf5_name(dir + info + ".hdf5", "return_probability", arma::hdf5_opts::append));
    _ndown_spreading.save(arma::hdf5_name(dir + info + ".hdf5", "spreading", arma::hdf5_opts::append));

//     if(dim <= dim_cut)
//     {
//         arma::Mat<ui::element_type> mat_elem = V.t() * n_spinup_center * V;
//         std::cout << " - - - - - - finished matrix elements SpinUp in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
//         start = std::chrono::system_clock::now();
//         // arma::Col<ui::element_type> diag_mat = arma::diagvec(mat_elem);
//         // diag_mat.save(arma::hdf5_name(dir + info + ".hdf5", "Diagonal_spinup", arma::hdf5_opts::append));
//         start = std::chrono::system_clock::now();
//         arma::vec autocorrelation(times.size());
//         double DE = 0;
//         for(long alfa = 0; alfa < dim; alfa++)
//             DE += std::abs(mat_elem(alfa, alfa) * std::conj(mat_elem(alfa, alfa)));
//         DE = DE / dim;
//     #pragma omp parallel
//         for(long t_idx = 0; t_idx < times.size(); t_idx++)
//         {
//             double time = times(t_idx);
//             double auto_evo = 0.0;
//             for(long alfa = 0; alfa < dim; alfa++){
//                 for(long beta = alfa+1; beta < dim; beta++){
//                     double wnm = E(beta) - E(alfa);
//                     auto_evo += std::abs(mat_elem(alfa, beta) * std::conj(mat_elem(alfa, beta))) * std::cos(wnm * time);
//                 }   
//             }
//             autocorrelation(t_idx) = auto_evo / double(dim) + DE;
//         }
//         autocorrelation.save(arma::hdf5_name(dir + info + ".hdf5", "Autocorrelation_spinup", arma::hdf5_opts::append));
//         std::cout << " - - - - - - finished autocorrelation SpinUp in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
//         start = std::chrono::system_clock::now();
//         std::cout << " - - - - - - finished spectral function SpinDown in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
//         start = std::chrono::system_clock::now();
//         mat_elem = V.t() * n_spindown_center * V;
//         std::cout << " - - - - - - finished matrix elements SpinDown in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
//         start = std::chrono::system_clock::now();
//         // diag_mat = arma::diagvec(mat_elem);
//         // diag_mat.save(arma::hdf5_name(dir + info + ".hdf5", "Diagonal_spindown", arma::hdf5_opts::append));
//         start = std::chrono::system_clock::now();
//         autocorrelation = arma::vec(times.size());
//         DE = 0;
//         for(long alfa = 0; alfa < dim; alfa++)
//             DE += std::abs(mat_elem(alfa, alfa) * std::conj(mat_elem(alfa, alfa)));
//         DE = DE / dim;
//     #pragma omp parallel
//         for(long t_idx = 0; t_idx < times.size(); t_idx++)
//         {
//             double time = times(t_idx);
//             double auto_evo = 0.0;
//             for(long alfa = 0; alfa < dim; alfa++){
//                 for(long beta = alfa+1; beta < dim; beta++){
//                     double wnm = E(beta) - E(alfa);
//                     auto_evo += std::abs(mat_elem(alfa, beta) * std::conj(mat_elem(alfa, beta))) * std::cos(wnm * time);
//                 }   
//             }
//             autocorrelation(t_idx) = auto_evo / double(dim) + DE;
//         }
//         autocorrelation.save(arma::hdf5_name(dir + info + ".hdf5", "Autocorrelation_spindown", arma::hdf5_opts::append));
//         std::cout << " - - - - - - finished autocorrelation SpinDown in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
//         start = std::chrono::system_clock::now();
//         std::cout << " - - - - - - finished spectral function SpinDown in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
// }
    std::cout << " - - - - - - finished time evolution in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    std::cout << " - - - - - - finished evolution of spin down and spreading in dilute limit in : " << tim_s(start_tot) << " s - - - - - - " << std::endl; // simulation end
}

/// @brief Compare energie spactra for full model and all symmetry sectors combined
void ui::dimer_dynamics()
{  
     std::string dir = this->saving_dir + "DimerDecay" + kPSep;
	createDirs(dir);
	
	size_t dim = this->ptr_to_model->get_hilbert_size();
	std::string info = this->set_info();

    size_t dim_cut = 1e5;
    // arma::sp_cx_mat H = cast_cx_sparse(this->ptr_to_model->get_hamiltonian());
	const size_t size = dim > dim_cut? this->l_steps : dim;
    
    auto _hilbert_space = this->ptr_to_model->get_model_ref().get_hilbert_space();
    auto _group = _hilbert_space.get_symmetry_group();


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

    QHS::Symmetric2DLattice<> _hilbert_full(this->Lx, this->Ly, this->syms.Nup, this->syms.Ndown, v_1d<QOps::generic_operator<cpx, PackedLattice2D>>(), this->boundary_conditions, 0, 0, 1, 1, true, false);
    u64 dimtot = _hilbert_full.get_hilbert_space_size();

    arma::sp_mat DimerNumber(dimtot, dimtot);
    arma::sp_mat SpinDownNumber(dimtot, dimtot);
    
    for(long alfa = 0; alfa < dimtot; alfa++)
    {
        auto base_state = _hilbert_full(alfa);
        for(std::size_t y = 0; y < this->Ly; y++){
            for(std::size_t x = 0; x < this->Lx; x++){
                auto _conf = base_state.get_block(x, y);
                if(_conf == DIMER_RIGHT) 
                    DimerNumber(alfa, alfa) += 1.0;
                if(_conf == DIMER_UP) 
                    DimerNumber(alfa, alfa) += 1.0;

                if(_conf == DOWN_MONOMER) 
                    SpinDownNumber(alfa, alfa) += 1.0;
            }
        }
    }
    
    #ifdef USE_SYMMETRIES
        auto U = _hilbert_space.symmetry_rotation(_hilbert_full);
    #endif

    arma::vec _spin_down_num(size, arma::fill::zeros);
    arma::vec _dimer_num(size, arma::fill::zeros);
    arma::vec _ipr(size, arma::fill::zeros);

#pragma omp parallel for
    for(int k = 0; k < size; k++)
    {
        #ifdef USE_SYMMETRIES
            arma::Col<ui::element_type> state = U * this->ptr_to_model->get_eigenState(k);
        #else
            arma::Col<ui::element_type> state = this->ptr_to_model->get_eigenState(k);
        #endif
        double tmp = 0;
        for(long alfa = 0; alfa < dimtot; alfa++)
            tmp += std::norm(state(alfa)) * std::norm(state(alfa));
        
        _ipr(k) = tmp;
        _dimer_num(k) = std::real( arma::cdot(state, DimerNumber * state) );
        _spin_down_num(k) = std::real( arma::cdot(state, SpinDownNumber * state) );
        
    }
    E.save(arma::hdf5_name(dir + info + ".hdf5", "energies"));
    _dimer_num.save(arma::hdf5_name(dir + info + ".hdf5", "DimerNumber_n", arma::hdf5_opts::append));
    _spin_down_num.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDownNumber_n", arma::hdf5_opts::append));
    _ipr.save(arma::hdf5_name(dir + info + ".hdf5", "IPR", arma::hdf5_opts::append));

    std::cout << " - - - - - - finished eigenstate measures in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    for(int Ndimer = 0; Ndimer <= std::min(this->syms.Nup, this->syms.Ndown); Ndimer++){
        start = std::chrono::system_clock::now();
        PackedLattice2D initial_state(this->Lx, this->Ly);
        // if(this->syms.Nup > this->syms.Ndown)   initial_state.fill_all(UP_MONOMER);
        // else                                    initial_state.fill_all(DOWN_MONOMER);
        u64 index;
        bool succeed = false;
        while(!succeed){
            initial_state.fill_all(UP_MONOMER);

            std::mt19937 rng(std::random_device{}());

            std::vector<std::vector<bool>> occupied(
                Lx, std::vector<bool>(Ly, false));

            if(Ndimer == 0){

                // all lattice sites
                std::vector<std::pair<int,int>> sites;
                for(int x = 0; x < Lx; ++x)
                    for(int y = 0; y < Ly; ++y)
                        sites.emplace_back(x,y);

                std::shuffle(sites.begin(), sites.end(), rng);

                for(int n = 0; n < this->syms.Ndown; ++n){
                    auto [x,y] = sites[n];
                    initial_state.set_block(x, y, DOWN_MONOMER);
                }

            } else {

                std::vector<std::tuple<int,int,int,int>> dimers;

                std::uniform_int_distribution<int> distX(0, Lx - 1);
                std::uniform_int_distribution<int> distY(0, Ly - 1);
                std::uniform_int_distribution<int> distDir(0, 1);

                while ((int)dimers.size() < Ndimer)
                {
                    int x = distX(rng);
                    int y = distY(rng);

                    int dir = distDir(rng);

                    int xp = x;
                    int yp = y;

                    if(dir == 0) xp = x + 1;   // horizontal
                    else         yp = y + 1;   // vertical

                    if (xp >= Lx || yp >= Ly)
                        continue;

                    if (occupied[x][y] || occupied[xp][yp])
                        continue;

                    occupied[x][y]   = true;
                    occupied[xp][yp] = true;

                    dimers.emplace_back(x,y,xp,yp);
                }

                // place dimers
                for(auto [x1,y1,x2,y2] : dimers){

                    if(x2 != x1){
                        initial_state.set_block(x1, y1, DIMER_RIGHT);
                        initial_state.set_block(x2, y2, DIMER_LEFT);
                    } else {
                        initial_state.set_block(x1, y1, DIMER_UP);
                        initial_state.set_block(x2, y2, DIMER_DOWN);
                    }
                }

                // ---- place remaining DOWN_MONOMERs ----

                int Nmono = this->syms.Ndown - Ndimer;

                std::vector<std::pair<int,int>> empty_sites;

                for(int x = 0; x < Lx; ++x)
                    for(int y = 0; y < Ly; ++y)
                        if(!occupied[x][y])
                            empty_sites.emplace_back(x,y);

                std::shuffle(empty_sites.begin(), empty_sites.end(), rng);

                for(int n = 0; n < Nmono; ++n){
                    auto [x,y] = empty_sites[n];
                    initial_state.set_block(x, y, DOWN_MONOMER);
                }
            }

            initial_state.print_state();
            std::cout << std::endl;
            // continue;

            #ifdef USE_SYMMETRIES
                auto [new_state, sym_eig] = _hilbert_space.find_SEC_representative(initial_state);
                index = _hilbert_space.find(new_state);
                if(index >= dim){
                    std::cout << "ERROR: Initial state does not exist in this symmetry sector.." << std::endl;
                    // return;
                    succeed = false;
                } else
                    succeed = true;
            #else 
                index = _hilbert_space.find(initial_state);
                succeed = true;
            #endif
        }
        // goto
        u64 index_full = _hilbert_full.find(initial_state);

        std::cout << "\n" << index << std::endl << std::endl;
        arma::Col<element_type> coeff = V.row(index).t();

        arma::cx_vec init_state(dim, arma::fill::zeros);
        init_state(index) = 1.0;
        std::cout << " - - - - - - finished setting up initial state in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
        start = std::chrono::system_clock::now();
        arma::Col<int> DimerNumber_dummy(dim);
        for(int k = 0; k < dim; k++)
            DimerNumber_dummy(k) = int( DimerNumber(k, k) );

        arma::uvec idx = arma::find(DimerNumber_dummy == Ndimer);
        u64 idx_start = idx.front();
        u64 idx_stop  = idx.back();
        std::cout << " - - - - - - finished findinf initial states with fixed dimer number with start = " << idx_start << " and stop = " << idx_stop << " in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
        start = std::chrono::system_clock::now();
        arma::vec _ipr_init_state(size, arma::fill::zeros);
        for(int k = 0; k < size; k++)
        {
            #ifdef USE_SYMMETRIES
                arma::Col<ui::element_type> state = U * this->ptr_to_model->get_eigenState(k);
            #else
                arma::Col<ui::element_type> state = this->ptr_to_model->get_eigenState(k);
            #endif
            double x = std::norm(state(index));
            _ipr_init_state(k) += x;
        }
        _ipr_init_state.save(arma::hdf5_name(dir + info + ".hdf5", "Ndimer=" + std::to_string(Ndimer) + "/IPR_initial_state", arma::hdf5_opts::append));
        std::cout << " - - - - - - finished setting IPR of initial state in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
        
        start = std::chrono::system_clock::now();
        
        arma::vec times = arma::logspace(-2, std::log10(dim) + 1, 1000);
        this->dt = times(1) - times(0);
        _spin_down_num = arma::vec(times.size(), arma::fill::zeros);
        _dimer_num = arma::vec(times.size(), arma::fill::zeros);

        arma::vec return_probability(times.size(), arma::fill::zeros);
        times.save(arma::hdf5_name(dir + info + ".hdf5", "times", arma::hdf5_opts::append));
        
        arma::vec weights(dim, arma::fill::zeros);
    
    #pragma omp for schedule(dynamic)
        for(long alfa = 0; alfa < dim; ++alfa)
        {
        #ifdef USE_SYMMETRIES
            arma::Col<ui::element_type> state =
                U * this->ptr_to_model->get_eigenState(alfa);
        #else
            arma::Col<ui::element_type> state =
                this->ptr_to_model->get_eigenState(alfa);
        #endif

            weights(alfa) = arma::accu(arma::square(arma::abs(state.rows(idx_start, idx_stop))));
        }
    #pragma omp for schedule(dynamic)
        for(long t_idx = 0; t_idx < (long)times.size(); t_idx++)
        {
            const double time = times(t_idx);
            cpx overlap = 0.0;
            for(long alfa = 0; alfa < dim; alfa++)
                overlap += std::exp(-1i * time * E(alfa)) * weights(alfa);
            return_probability(t_idx) = std::norm(overlap);
        }
        return_probability = return_probability / double(idx_stop - idx_start);
        arma::Col<unsigned long>({idx_start}).save(arma::hdf5_name(dir + info + ".hdf5", "Ndimer=" + std::to_string(Ndimer) + "/idx_start", arma::hdf5_opts::append));
        arma::Col<unsigned long>({idx_stop}).save(arma::hdf5_name(dir + info + ".hdf5", "Ndimer=" + std::to_string(Ndimer) + "/idx_stop", arma::hdf5_opts::append));
        return_probability.save(arma::hdf5_name(dir + info + ".hdf5", "Ndimer=" + std::to_string(Ndimer) + "/return_probability", arma::hdf5_opts::append));
        std::cout << " - - - - - - finished return probability of Ndimer = " << Ndimer << "in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    //     start = std::chrono::system_clock::now();
    // // #pragma omp parallel
    //     {
    //         // Per-thread preallocated buffer → no allocation inside loop
    //         // arma::cx_vec init_state(dim, arma::fill::zeros);

    //         // #pragma omp for schedule(dynamic)
    //         for(long t_idx = 0; t_idx < (long)times.size(); t_idx++)
    //         {
    //             const double time = times(t_idx);

    //             // ── Build init_state ─────────────────────────────────────────────────
    //             // init_state.zeros();  // reuse buffer, just zero it
    //             arma::cx_vec init_state(dim, arma::fill::zeros);
    //             if(dim <= dim_cut){
    //                 for(long alfa = 0; alfa < dim; alfa++)
    //                 {
    //                     // auto state = V.col(alfa);
    //                     auto state = this->ptr_to_model->get_eigenState(alfa);
    //                     init_state += std::exp(-1i * time * E(alfa)) * state * coeff(alfa);
    //                 }
    //                 init_state = arma::normalise(init_state);
    //             }
    //             #ifdef USE_SYMMETRIES
    //                 init_state = U * init_state;
    //             #endif

    //             _dimer_num(t_idx) = std::real( arma::cdot(init_state, DimerNumber * init_state) );
    //             _spin_down_num(t_idx) = std::real( arma::cdot(init_state, SpinDownNumber * init_state) );
    //         }
    //     }
    //     _spin_down_num.save(arma::hdf5_name(dir + info + ".hdf5", "Ndimer=" + std::to_string(Ndimer) + "/SpinDownNumber(t)", arma::hdf5_opts::append));
    //     _dimer_num.save(arma::hdf5_name(dir + info + ".hdf5", "Ndimer=" + std::to_string(Ndimer) + "/DimerNumber(t)", arma::hdf5_opts::append));
    //     std::cout << " - - - - - - finished time evolution of Ndimer = " << Ndimer << " in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    }
    std::cout << " - - - - - - finished evolution of dimer test in : " << tim_s(start_tot) << " s - - - - - - " << std::endl; // simulation end        
    // std::string dir = this->saving_dir + "DimerDynamics_Test" + kPSep;
	// createDirs(dir);
	
	// size_t dim = this->ptr_to_model->get_hilbert_size();
	// std::string info = this->set_info();

    // size_t dim_cut = 1e5;
    // arma::sp_cx_mat H = cast_cx_sparse(this->ptr_to_model->get_hamiltonian());
	// const size_t size = dim > dim_cut? this->l_steps : dim;

	// clk::time_point start = std::chrono::system_clock::now();
    // clk::time_point start_tot = start;
    // if(dim > dim_cut){
    //     this->ptr_to_model->diag_sparse(this->l_steps, this->l_bundle, this->tol, this->seed);	
    // }
    // else{
    //     this->ptr_to_model->diagonalization();
    // }
    // const arma::vec E = this->ptr_to_model->get_eigenvalues();
    // const auto& V = this->ptr_to_model->get_eigenvectors();
    // std::cout << " - - - - - - finished diagonalization in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    // start = std::chrono::system_clock::now();

    // PackedLattice2D initial_state(this->Lx, this->Ly);
    // initial_state.fill_all(UP_MONOMER);
    // // int N = std::min(this->)
    
    // int x0 = (this->Lx - 1.0)/2;
    // int y0 = (this->Ly - 1.0)/2;
    // initial_state.set_block(x0, y0, DIMER_RIGHT);
    // initial_state.set_block(x0+1, y0, DIMER_LEFT);
    
    // auto _hilbert_space = this->ptr_to_model->get_model_ref().get_hilbert_space();
    // u64 index = _hilbert_space.find(initial_state);
    // initial_state.print_state();
    // std::cout << "\n" << index << std::endl << std::endl;
    // arma::Col<element_type> coeff = V.row(index).t();

    // E.save(arma::hdf5_name(dir + info + ".hdf5", "energies"));
	// // auto H = this->ptr_to_model->get_dense_hamiltonian();
    // // H = arma::abs(H);
    // // H.save(arma::hdf5_name(dir + info + ".hdf5", "Hamiltonian", arma::hdf5_opts::append));


    // // Lochscmidt echo, total Nd on monomeres, Nd Nd correlation function, nd(k, t) ? energies as function of m

    // arma::cx_vec init_state(dim, arma::fill::zeros);
    // init_state(index) = 1.0;
    // std::cout << " - - - - - - finished setting up initial state in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    // start = std::chrono::system_clock::now();
    
    // arma::vec _ndimer_density(size, arma::fill::zeros);
    // arma::vec _ndimer_density_center(size, arma::fill::zeros);
    // arma::vec _nmonomer_density(size, arma::fill::zeros);
    // for(int k = 0; k < size; k++)
    // {
    //     arma::Col<ui::element_type> state = this->ptr_to_model->get_eigenState(k);
    //     arma::mat _ndimer_density_state(this->Ly, this->Lx, arma::fill::zeros);
    //     arma::mat _nmonomer_density_state(this->Ly, this->Lx, arma::fill::zeros);
    //     for(long alfa = 0; alfa < dim; alfa++)
    //     {
    //         auto base_state = _hilbert_space(alfa);
    //     // #pragma omp parallel for
    //         for(std::size_t y = 0; y < this->Ly; y++){
    //             for(std::size_t x = 0; x < this->Lx; x++){
    //                 auto _conf = base_state.get_block(x, y);
    //                 _ndimer_density_state(y, x) += _dimer_number(_conf) * std::abs(state(alfa)) * std::abs(state(alfa));
    //                 _nmonomer_density_state(y, x) += _fermion_number(_conf) * std::abs(state(alfa)) * std::abs(state(alfa));
    //             }
    //         }
    //         auto _conf = base_state.get_block(x0, y0);
    //         if(_conf == DIMER_RIGHT)
    //             _ndimer_density_center(k) += std::abs(state(alfa)) * std::abs(state(alfa));
    //     }
	//     {
    //         _ndimer_density_state.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDown_State/n = " + std::to_string(k), arma::hdf5_opts::append));
    //         _nmonomer_density_state.save(arma::hdf5_name(dir + info + ".hdf5", "Monomer_State/n = " + std::to_string(k), arma::hdf5_opts::append));
    //     }
    //     _ndimer_density(k) = arma::accu( _ndimer_density_state ) / (this->Lx * this->Ly);
    //     _nmonomer_density(k) = arma::accu( _nmonomer_density_state ) / (this->Lx * this->Ly);
    // }
    // _ndimer_density.save(arma::hdf5_name(dir + info + ".hdf5", "DimerDensity_states", arma::hdf5_opts::append));
    // _nmonomer_density.save(arma::hdf5_name(dir + info + ".hdf5", "MonomerDensity_states", arma::hdf5_opts::append));
    // _ndimer_density_center.save(arma::hdf5_name(dir + info + ".hdf5", "DimerCenterDensity_states", arma::hdf5_opts::append));
    // std::cout << " - - - - - - finished eigenstate measures in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    // start = std::chrono::system_clock::now();

    // arma::vec times = arma::linspace(0, 50, 1001);
    // // arma::vec times = arma::logspace(1e-2, int(std::log10(dim)), 1001);
    // this->dt = times(1) - times(0);
    // _ndimer_density = arma::vec(times.size(), arma::fill::zeros);
    // _nmonomer_density = arma::vec(times.size(), arma::fill::zeros);
    // _ndimer_density_center = arma::vec(times.size(), arma::fill::zeros);
    
    // arma::vec _ndimer_spreading(times.size(), arma::fill::zeros);
    // arma::vec return_probability(times.size(), arma::fill::zeros);

    // times.save(arma::hdf5_name(dir + info + ".hdf5", "times", arma::hdf5_opts::append));
    // for(long t_idx = 0; t_idx < times.size(); t_idx++)
    // {
    //     double time = times(t_idx);
    //     if(dim <= dim_cut){
    //         init_state = arma::cx_vec(dim, arma::fill::zeros);
    //         for(long alfa = 0; alfa < dim; alfa++)
    //         {
    //             // auto state = V.col(alfa);
    //             arma::Col<ui::element_type> state = this->ptr_to_model->get_eigenState(alfa);
    //             init_state += std::exp(-1i * time * E(alfa)) * state * coeff(alfa);
    //         }
    //         init_state = arma::normalise(init_state);
    //     }
    //     return_probability(t_idx) = std::abs(init_state(index) * init_state(index));

    //     arma::mat _nmonomer_density_t(this->Ly, this->Lx, arma::fill::zeros);
    //     for(long alfa = 0; alfa < dim; alfa++)
    //     {
    //         auto base_state = _hilbert_space(alfa);
    //     // #pragma omp parallel for
    //         for(std::size_t y = 0; y < this->Ly; y++){
    //             for(std::size_t x = 0; x < this->Lx; x++){
    //                 auto _conf = base_state.get_block(x, y);
    //                 _ndimer_density(t_idx) += _dimer_number(_conf) / (this->Lx * this->Ly) * std::abs(init_state(alfa)) * std::abs(init_state(alfa));
    //                 _nmonomer_density_t(y, x) += _fermion_number(_conf) * std::abs(init_state(alfa)) * std::abs(init_state(alfa));
                    
    //                 // if(x < this->Lx - 1 && y < this->Ly-1)
    //                 {
    //                     if(_conf == DIMER_UP){
    //                         _ndimer_spreading(t_idx) += ((x - (x0+0.5))*(x - (x0+0.5)) + (y+0.5 - y0)*(y+0.5 - y0)) * std::abs(init_state(alfa)) * std::abs(init_state(alfa));
    //                     } else if(_conf == DIMER_RIGHT){
    //                         _ndimer_spreading(t_idx) += ((x - x0)*(x - x0) + (y - y0)*(y - y0)) * std::abs(init_state(alfa)) * std::abs(init_state(alfa));
    //                     }

    //                 }
    //             }
    //         }

    //         auto _conf = base_state.get_block(x0, y0);
    //         if(_conf == DIMER_RIGHT)
    //             _ndimer_density_center(t_idx) += std::abs(init_state(alfa)) * std::abs(init_state(alfa));
    //     }
    //     // _ndimer_density(t_idx) = arma::accu( _ndimer_density_t ) / (this->Lx * this->Ly);
    //     _nmonomer_density(t_idx) = arma::accu( _nmonomer_density_t ) / (this->Lx * this->Ly);

    //     if(dim > dim_cut){
    //         lanczos::Lanczos<cpx> lancz(H, this->l_steps, -1, this->tol, this->seed, true, false, init_state);
    //         lancz.time_evolution_step(init_state, this->dt);
    //     }
        
	//     if( t_idx%20 == 0){
    //         // _ndimer_density_t.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDimerDensity_snapshot/t = " + to_string_prec(time), arma::hdf5_opts::append));
    //         _nmonomer_density_t.save(arma::hdf5_name(dir + info + ".hdf5", "MonomerDensity_snapshot/t = " + to_string_prec(time), arma::hdf5_opts::append));
    //     }
    //     // arma::vec state_abs = arma::abs(init_state);
    //     // state_abs.save(arma::hdf5_name("GaugeViolation.hdf5", "state: t = " + to_string_prec(time), arma::hdf5_opts::append));
    // }
    // _ndimer_density.save(arma::hdf5_name(dir + info + ".hdf5", "DimerDensity", arma::hdf5_opts::append));
    // _nmonomer_density.save(arma::hdf5_name(dir + info + ".hdf5", "MonomerDensity", arma::hdf5_opts::append));
    // _ndimer_density_center.save(arma::hdf5_name(dir + info + ".hdf5", "DimerCenterDensity", arma::hdf5_opts::append));
    // return_probability.save(arma::hdf5_name(dir + info + ".hdf5", "return_probability", arma::hdf5_opts::append));
    // _ndimer_spreading.save(arma::hdf5_name(dir + info + ".hdf5", "spreading", arma::hdf5_opts::append));

    // std::cout << " - - - - - - finished time evolution in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    // std::cout << " - - - - - - finished evolution of dimer test in : " << tim_s(start_tot) << " s - - - - - - " << std::endl; // simulation end
}

/// @brief Compare energie spactra for full model and all symmetry sectors combined
void ui::neel_order()
{          
    std::string dir = this->saving_dir + "NeelOrder" + kPSep;
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
    // const arma::SpMat<ui::element_type> H = this->ptr_to_model->get_hamiltonian();
	
    cpx meanH = arma::trace(H) / double(dim);
    H = H - meanH * arma::eye<arma::SpMat<ui::element_type>>(dim, dim);
    arma::sp_cx_mat H2 = H * H;

    arma::vec gec(dim, arma::fill::zeros);
#pragma omp parallel for
    for(u64 k = 0; k < dim; k++)
        gec(k) = std::real(2.0 * (cpx)H2(k,k) - (cpx)H(k,k) * (cpx)H(k,k));

    std::cout << " - - - - - - finished diagonalization in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();

    auto _hilbert_space = this->ptr_to_model->get_model_ref().get_hilbert_space();
    auto _group = _hilbert_space.get_symmetry_group();

    PackedLattice2D initial_state(this->Lx, this->Ly);
    initial_state.fill_all(DOWN_MONOMER);
    // int N = std::min(this->)
    int Nu = this->syms.Nup;
    for(int ell = 0; ell < this->L; ell++)
    {
        if(Nu <= 0) continue;
        std::size_t x = ell % this->Lx;
        std::size_t y = ell / this->Lx;
        initial_state.set_block(x, y, UP_MONOMER);
        Nu--;
    }
        
    initial_state.print_state();
    std::cout << std::endl;

    PackedLattice2D initial_state2(this->Lx, this->Ly);
    initial_state2.fill_all(DOWN_MONOMER);
    Nu = this->syms.Nup;
    for(int x = 0; x < this->Lx; x++){
        for(int y = 0; y < this->Ly; y++){
            if(Nu <= 0) continue;
            if(x%2 ^ y%2){
                initial_state2.set_block(x, y, UP_MONOMER);
                Nu--;
            }
        }
    }
        
    initial_state2.print_state();
    std::cout << std::endl;

    #ifdef USE_SYMMETRIES
        auto [new_state, sym_eig] = _hilbert_space.find_SEC_representative(initial_state);
        u64 index = _hilbert_space.find(new_state);
        if(index >= dim){
            std::cout << "ERROR: Initial state does not exist in this symmetry sector.." << std::endl;
            return;
        }
        auto [new_state2, sym_eig2] = _hilbert_space.find_SEC_representative(initial_state2);
        u64 index2 = _hilbert_space.find(new_state2);
        if(index2 >= dim){
            std::cout << "ERROR: Initial state does not exist in this symmetry sector.." << std::endl;
            return;
        }
    #else 
        u64 index = _hilbert_space.find(initial_state);
        u64 index2 = _hilbert_space.find(initial_state2);
    #endif
    std::cout << "\n" << index << "\n" << index2 << std::endl << std::endl;
    std::cout << " - - - - - - " << std::endl;

    arma::Col<element_type> coeff = V.row(index).t();
    arma::Col<element_type> coeff2 = V.row(index2).t();

    E.save(arma::hdf5_name(dir + info + ".hdf5", "energies"));
    gec.save(arma::hdf5_name(dir + info + ".hdf5", "gec", arma::hdf5_opts::append));
    // {
    //     arma::mat x = arma::square( arma::abs(V) );
    //     x.save(arma::hdf5_name(dir + info + ".hdf5", "|c^i_n|^2", arma::hdf5_opts::append));
    // }
    arma::vec x = arma::real(coeff); x.save(arma::hdf5_name(dir + info + ".hdf5", "DomainWall/overlaps_re", arma::hdf5_opts::append));
    x = arma::imag(coeff); x.save(arma::hdf5_name(dir + info + ".hdf5", "DomainWall/overlaps_im", arma::hdf5_opts::append));
    x = arma::real(coeff2); x.save(arma::hdf5_name(dir + info + ".hdf5", "Neel/overlaps_re", arma::hdf5_opts::append));
    x = arma::imag(coeff2); x.save(arma::hdf5_name(dir + info + ".hdf5", "Neel/overlaps_im", arma::hdf5_opts::append));
	// auto H = this->ptr_to_model->get_dense_hamiltonian();Neel/
    // H = arma::abs(H);
    // H.save(arma::hdf5_name(dir + info + ".hdf5", "Hamiltonian", arma::hdf5_opts::append));


    // Lochscmidt echo, total Nd on monomeres, Nd Nd correlation function, nd(k, t) ? energies as function of m
    std::cout << " - - - - - - finished setting up initial state in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();
    
    arma::sp_mat SpinSpinCorr(dim, dim);
    arma::sp_mat SpinUpDens(dim, dim);
    arma::sp_mat SpinDownDens(dim, dim);
    for(long alfa = 0; alfa < dim; alfa++)
    {
        auto base_state = _hilbert_space(alfa);
        #ifdef USE_SYMMETRIES
            cpx _norm = _hilbert_space.get_norm(alfa);  // normalization of sym. state
            // double _denom = std::abs(_norm * std::conj(_norm)) * double(_group.size());
            double _denom = double(_group.size());
            for(auto& g : _group){
                auto [Gstate, val] = g(base_state);
                auto _conf1 = Gstate.get_block(0, 0);
                SpinUpDens(alfa, alfa) += _fermion_up_number(_conf1) / _denom;
                SpinDownDens(alfa, alfa) += _fermion_down_number(_conf1) / _denom;
                auto _conf = Gstate.get_block(this->Lx/2, this->Ly/2);
                if(_fermion_number(_conf))
                    SpinSpinCorr(alfa, alfa) += (_fermion_number(_conf1) - 0.5) * (_fermion_number(_conf) - 0.5) / _denom;
                // for(int ell = 0; ell < this->L; ell++)
                // {
                //     std::size_t x = ell % this->Lx;
                //     std::size_t y = ell / this->Lx;
                //     // double sign = (x%2 ^ y%2)? -1. : 1.;
                //     auto _conf = Gstate.get_block(x, y);
                //     if(_fermion_number(_conf))
                //         SpinSpinCorr(alfa, alfa) += (_fermion_number(_conf1) - 0.5) * (_fermion_number(_conf) - 0.5) / _denom;
                // }
            }
        #else
            auto _conf1 = base_state.get_block(0, 0);
            // MonDens(alfa, alfa) += _fermion_number(_conf1);
            SpinUpDens(alfa, alfa) += _fermion_up_number(_conf1);
            SpinDownDens(alfa, alfa) += _fermion_down_number(_conf1);
            auto _conf = base_state.get_block(this->Lx/2, this->Ly/2);
            if(_fermion_number(_conf))
                SpinSpinCorr(alfa, alfa) += (_fermion_number(_conf1) - 0.5) * (_fermion_number(_conf) - 0.5);
            // for(int ell = 0; ell < this->L; ell++)
            // {
            //     std::size_t x = ell % this->Lx;
            //     std::size_t y = ell / this->Lx;
            //     // double sign = (x%2 ^ y%2)? -1. : 1.;
            //     auto _conf = base_state.get_block(x, y);
            //     if(_fermion_number(_conf))
            //         SpinSpinCorr(alfa, alfa) += (_fermion_number(_conf1) - 0.5) * (_fermion_number(_conf) - 0.5);
            // }
        #endif
        // std::cout << SpinSpinCorr(alfa, alfa) << std::endl;
        // base_state.print_state();
        // std::cout << std::endl;
    }
    // SpinSpinCorr /= std::sqrt(this->L);
    // SpinSpinCorr *= SpinSpinCorr;
    // std::cout << "Norm MonDens\t" << arma::trace(MonDens * MonDens.t()) / dim << std::endl;
    // std::cout << "Norm SpinSpinCorr\t" << arma::trace(SpinSpinCorr * SpinSpinCorr.t()) / dim << std::endl;
    // MonDens = MonDens / std::sqrt(arma::trace(MonDens * MonDens.t()) / dim);
    // SpinSpinCorr = SpinSpinCorr / std::sqrt(arma::trace(SpinSpinCorr * SpinSpinCorr.t()) / dim);
    // MonDens = MonDens / (arma::norm(MonDens, "fro") / std::sqrt(dim) );
    // std::cout << arma::mat(MonDens) << std::endl << std::endl;
    // std::cout << arma::mat(SpinSpinCorr) << std::endl << std::endl;
    
    std::cout << " - - - - - - finished setting up matrices in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();
    
    arma::vec _spinspin_corr_fun(size, arma::fill::zeros);
    arma::vec _nspinup_density(size, arma::fill::zeros);
    arma::vec _nspindown_density(size, arma::fill::zeros);
    arma::vec _ipr(size, arma::fill::zeros);
    arma::vec _ipr_red(size, arma::fill::zeros);
    // arma::vec _imbalance(size, arma::fill::zeros);
#pragma omp parallel for
    for(int k = 0; k < size; k++)
    {
        arma::Col<ui::element_type> state = this->ptr_to_model->get_eigenState(k);
        _nspinup_density(k) = std::real( arma::cdot(state, SpinUpDens * state) );
        _nspindown_density(k) = std::real( arma::cdot(state, SpinDownDens * state) );
        _spinspin_corr_fun(k) = std::real( arma::cdot(state, SpinSpinCorr * state) );
        for(long alfa = 0; alfa < dim; alfa++)
        {
            auto base_state = _hilbert_space(alfa);
            #ifdef USE_SYMMETRIES
                double tmp = 0;
                cpx _norm = _hilbert_space.get_norm(alfa);  // normalization of sym. state
                double _denom2 = std::abs(_norm * std::conj(_norm)) * double(_group.size());
                double coeff4 = std::norm(state(alfa)) * std::norm(state(alfa));
                // fullstate(idxtot) += chi_g / _norm_k / std::sqrt(_group.size());
                for(auto& g : _group){
                    auto [Gstate, chi_g] = g(base_state);
                    cpx x = chi_g * state(alfa);
                    tmp += std::abs(x * std::conj(x)) * std::abs(x * std::conj(x));
                }
                _ipr(k) += tmp / (_denom2 * _denom2);
                _ipr_red(k) += coeff4;
            #else
                ui::element_type x = state(alfa);
                _ipr(k) = std::abs(x * std::conj(x)) * std::abs(x * std::conj(x));
            #endif
        }
    }
    _ipr.save(arma::hdf5_name(dir + info + ".hdf5", "IPR", arma::hdf5_opts::append));
    _ipr_red.save(arma::hdf5_name(dir + info + ".hdf5", "IPR_red", arma::hdf5_opts::append));
    arma::vec _ipr_ratio = _ipr / _ipr_red;
    _ipr_ratio.save(arma::hdf5_name(dir + info + ".hdf5", "IPR_ratio", arma::hdf5_opts::append));
    
    _nspinup_density.save(arma::hdf5_name(dir + info + ".hdf5", "SpinUpDensity_states", arma::hdf5_opts::append));
    _nspindown_density.save(arma::hdf5_name(dir + info + ".hdf5", "SpinDownDensity_states", arma::hdf5_opts::append));
    _spinspin_corr_fun.save(arma::hdf5_name(dir + info + ".hdf5", "OrderParameter_states", arma::hdf5_opts::append));
    std::cout << " - - - - - - finished eigenstate measures in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();
    // dim_cut = 2e3;

    arma::cx_vec init_state(dim, arma::fill::zeros);
    init_state(index) = 1.0;
    arma::cx_vec init_state2(dim, arma::fill::zeros);
    init_state2(index2) = 1.0;

    arma::vec times = dim > dim_cut? arma::linspace(0, 50, 1001) : arma::logspace(-2, std::log10(dim)+1, 1001);
    // arma::vec times = arma::logspace(1e-2, int(std::log10(dim)), 1001);
    this->dt = times(1) - times(0);
    _nspinup_density = arma::vec(times.size(), arma::fill::zeros);
    _nspindown_density = arma::vec(times.size(), arma::fill::zeros);
    _spinspin_corr_fun = arma::vec(times.size(), arma::fill::zeros);
    
    arma::vec return_probability(times.size(), arma::fill::zeros);
    arma::vec return_probability2(times.size(), arma::fill::zeros);
    arma::vec _nspinup_density2(times.size(), arma::fill::zeros);
    arma::vec _nspindown_density2(times.size(), arma::fill::zeros);
    arma::vec _spinspin_corr_fun2(times.size(), arma::fill::zeros);

    times.save(arma::hdf5_name(dir + info + ".hdf5", "times", arma::hdf5_opts::append));
// #pragma omp parallel for
    for(long t_idx = 0; t_idx < times.size(); t_idx++)
    {
        double time = times(t_idx);
        if(dim <= dim_cut){
            init_state = arma::cx_vec(dim, arma::fill::zeros);
            init_state2 = arma::cx_vec(dim, arma::fill::zeros);
            for(long alfa = 0; alfa < dim; alfa++)
            {
                // auto state = V.col(alfa);
                arma::Col<ui::element_type> state = this->ptr_to_model->get_eigenState(alfa);
                init_state += std::exp(-1i * time * E(alfa)) * state * coeff(alfa);
                init_state2 += std::exp(-1i * time * E(alfa)) * state * coeff2(alfa);
            }
            // init_state = arma::normalise(init_state);
        }
        return_probability(t_idx) = std::abs(init_state(index) * init_state(index));
        return_probability2(t_idx) = std::abs(init_state2(index2) * init_state2(index2));

        _nspinup_density(t_idx) = std::real( arma::cdot(init_state, SpinUpDens * init_state) );
        _nspindown_density(t_idx) = std::real( arma::cdot(init_state, SpinDownDens * init_state) );
        _spinspin_corr_fun(t_idx) = std::real( arma::cdot(init_state, SpinSpinCorr * init_state) );

        _nspinup_density2(t_idx) = std::real( arma::cdot(init_state2, SpinUpDens * init_state2) );
        _nspindown_density2(t_idx) = std::real( arma::cdot(init_state2, SpinDownDens * init_state2) );
        _spinspin_corr_fun2(t_idx) = std::real( arma::cdot(init_state2, SpinSpinCorr * init_state2) );
        if(dim > dim_cut){
            lanczos::Lanczos<cpx> lancz(H, this->l_steps, -1, this->tol, this->seed, true, false, init_state);
            lancz.time_evolution_step(init_state, this->dt);

            lancz = lanczos::Lanczos<cpx>(H, this->l_steps, -1, this->tol, this->seed, true, false, init_state2);
            lancz.time_evolution_step(init_state2, this->dt);
        }
    }
    _spinspin_corr_fun.save(arma::hdf5_name(dir + info + ".hdf5", "DomainWall/OrderParameter", arma::hdf5_opts::append));
    _nspinup_density.save(arma::hdf5_name(dir + info + ".hdf5", "DomainWall/SpinUpDensity", arma::hdf5_opts::append));
    _nspindown_density.save(arma::hdf5_name(dir + info + ".hdf5", "DomainWall/SpinDownDensity", arma::hdf5_opts::append));
    return_probability.save(arma::hdf5_name(dir + info + ".hdf5", "DomainWall/return_probability", arma::hdf5_opts::append));
    
    _spinspin_corr_fun2.save(arma::hdf5_name(dir + info + ".hdf5", "Neel/OrderParameter", arma::hdf5_opts::append));
    _nspinup_density2.save(arma::hdf5_name(dir + info + ".hdf5", "Neel/SpinUpDensity", arma::hdf5_opts::append));
    _nspindown_density2.save(arma::hdf5_name(dir + info + ".hdf5", "Neel/SpinDownDensity", arma::hdf5_opts::append));
    return_probability2.save(arma::hdf5_name(dir + info + ".hdf5", "Neel/return_probability", arma::hdf5_opts::append));
    std::cout << " - - - - - - finished time evolution in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    
    if(dim <= dim_cut)
    {
        arma::Mat<ui::element_type> mat_elem = V.t() * SpinUpDens * V;
        std::cout << " - - - - - - finished matrix elements SpinUp in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
        start = std::chrono::system_clock::now();
        arma::Col<ui::element_type> diag_mat = arma::diagvec(mat_elem);
        diag_mat.save(arma::hdf5_name(dir + info + ".hdf5", "Diagonal_spinup", arma::hdf5_opts::append));
        start = std::chrono::system_clock::now();
        arma::vec autocorrelation(times.size());
        double DE = 0;
        for(long alfa = 0; alfa < dim; alfa++)
            DE += std::abs(mat_elem(alfa, alfa) * std::conj(mat_elem(alfa, alfa)));
        DE = DE / dim;
        for(long t_idx = 0; t_idx < times.size(); t_idx++)
        {
            double time = times(t_idx);
            double auto_evo = 0.0;
            for(long alfa = 0; alfa < dim; alfa++){
                for(long beta = alfa+1; beta < dim; beta++){
                    double wnm = E(beta) - E(alfa);
                    auto_evo += std::abs(mat_elem(alfa, beta) * std::conj(mat_elem(alfa, beta))) * std::cos(wnm * time);
                }   
            }
            autocorrelation(t_idx) = auto_evo / double(dim) + DE;
        }
        autocorrelation.save(arma::hdf5_name(dir + info + ".hdf5", "Autocorrelation_spinup", arma::hdf5_opts::append));
        std::cout << " - - - - - - finished autocorrelation SpinUp in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
        start = std::chrono::system_clock::now();


        std::cout << " - - - - - - finished spectral function SpinDown in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
        start = std::chrono::system_clock::now();
        mat_elem = V.t() * SpinDownDens * V;

        std::cout << " - - - - - - finished matrix elements SpinDown in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
        start = std::chrono::system_clock::now();
        diag_mat = arma::diagvec(mat_elem);
        diag_mat.save(arma::hdf5_name(dir + info + ".hdf5", "Diagonal_spindown", arma::hdf5_opts::append));
        start = std::chrono::system_clock::now();
        autocorrelation = arma::vec(times.size());
        DE = 0;
        for(long alfa = 0; alfa < dim; alfa++)
            DE += std::abs(mat_elem(alfa, alfa) * std::conj(mat_elem(alfa, alfa)));
        DE = DE / dim;
        for(long t_idx = 0; t_idx < times.size(); t_idx++)
        {
            double time = times(t_idx);
            double auto_evo = 0.0;
            for(long alfa = 0; alfa < dim; alfa++){
                for(long beta = alfa+1; beta < dim; beta++){
                    double wnm = E(beta) - E(alfa);
                    auto_evo += std::abs(mat_elem(alfa, beta) * std::conj(mat_elem(alfa, beta))) * std::cos(wnm * time);
                }   
            }
            autocorrelation(t_idx) = auto_evo / double(dim) + DE;
        }
        autocorrelation.save(arma::hdf5_name(dir + info + ".hdf5", "Autocorrelation_spindown", arma::hdf5_opts::append));
        std::cout << " - - - - - - finished autocorrelation SpinDown in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
        start = std::chrono::system_clock::now();


        std::cout << " - - - - - - finished spectral function SpinDown in : " << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
    }
    std::cout << " - - - - - - finished Neel/DomainWall decay in : " << tim_s(start_tot) << " s - - - - - - " << std::endl; // simulation end
}


/// @brief Decay of Neel order
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
    std::size_t dim_cut = 1e4;

    auto _hilbert = this->ptr_to_model->get_model_ref().get_hilbert_space();
    auto _group = _hilbert.get_symmetry_group();


    arma::Col<elem_ty> GS;
    double dE, E0;
    arma::vec E;
    arma::Mat<ui::element_type> V;
    // arma::SpMat<elem_ty> H;
    if(dim > dim_cut)
    {
        std::cout << " - - - - - - Large Hilbert-space, chosen Lanczos - - - - - - " << std::endl; // simulation end
        arma::SpMat<elem_ty> H = this->ptr_to_model->get_hamiltonian();
        auto lancz0 = lanczos::Lanczos<elem_ty, converge::energies>(H, 20, 1000, 1e-14, this->seed, 1);
        lancz0.diagonalization();
        E = lancz0.get_eigenvalues();
        GS = lancz0.get_eigenstates().col(0);
        GS = arma::normalise(GS);
        E.save(arma::hdf5_name(dir_spin + info + ".hdf5", "lanczos energies"));
    } else {
        this->ptr_to_model->diagonalization();
        E = this->ptr_to_model->get_eigenvalues();
        V = this->ptr_to_model->get_eigenvectors();
        GS = arma::normalise(this->ptr_to_model->get_eigenState(0));
        E.save(arma::hdf5_name(dir_spin + info + ".hdf5", "energies"));
    }
    E0 = E(0);
    dE = E(E.size() - 1) - E0;
    std::cout << " - - - - - - finished ground state in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
    start = std::chrono::system_clock::now();

    //<! create path in Brillouin Zone
    std::vector<std::pair<double, double>> path;
    if(this->Lx > 1 & this->Ly > 1){
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
    } else {
        if(this->Ly == 1){
        for(int i = 0; i <= this->Lx / 2; i++)
            path.push_back( std::make_pair(i, 0) );
        } else if (this->Lx == 1){
        for(int i = 0; i <= this->Ly / 2; i++)
            path.push_back( std::make_pair(0, i) );
        }
    }

    std::cout << " - - - - - - finished path construction in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
    auto orbit_size = [&](int kx, int ky) -> int
    {
        bool square   = (this->Lx == this->Ly);
        
        // Check stabilizer conditions:
        // Px fixes q iff ky = 0 or ky = Ly/2
        // Py fixes q iff kx = 0 or kx = Lx/2
        // Pd fixes q iff kx = ky (only for square)
        
        bool Px_fixes = (ky == 0 || (this->Ly%2==0 && ky == this->Ly/2));
        bool Py_fixes = (kx == 0 || (this->Lx%2==0 && kx == this->Lx/2));
        bool Pd_fixes = square && (kx == ky);

        // Point group size
        int G_size = square ? 8 : 4;   // D4 or Z2xZ2

        // Stabilizer size: count how many generators fix q
        // Each independent fixing condition halves the orbit
        int stab_size = 1;
        if(Px_fixes) stab_size *= 2;
        if(Py_fixes) stab_size *= 2;
        if(Pd_fixes) stab_size *= 2;
        // Cap at group size (over-counting for high-symmetry points)
        stab_size = std::min(stab_size, G_size);

        return G_size / stab_size;
    };
    start = std::chrono::system_clock::now();

    double domega = dE / (double)E.size() / 4.;
    double eta = 3 * domega;
    arma::vec omegas = arma::regspace(-1.5*dE, domega, 1.5*dE);
    const arma::vec omegax = arma::logspace(std::log10(0.5/dim) - 1, std::log10( 3*dE ) + 1, 7 * this->L);

    using Func = std::function<std::pair<PackedLattice2D, double>(const PackedLattice2D&, int, int)>;
    Func Sz = [](const PackedLattice2D& state, int x, int y){
        auto _conf = state.get_block(x, y);
        if( _conf == UP_MONOMER )           return std::make_pair(state, 0.5);
        else if( _conf == DOWN_MONOMER )    return std::make_pair(state, -0.5);
        else //if( _conf == UP_MONOMER ) and dimers
            return std::make_pair(state, 0.0);
    };
    Func Splus = [](const PackedLattice2D& state, int x, int y){
        auto _conf = state.get_block(x, y);
        if( _conf == DOWN_MONOMER ){
            PackedLattice2D new_state = state;
            new_state.set_block(x, y, UP_MONOMER);
            return std::make_pair(new_state, 1.0);
        } else //if( _conf == UP_MONOMER ) and dimers
            return std::make_pair(state, 0.0);
    };
    Func Sminus = [](const PackedLattice2D& state, int x, int y){
        auto _conf = state.get_block(x, y);
        if( _conf == UP_MONOMER ){
            PackedLattice2D new_state = state;
            new_state.set_block(x, y, DOWN_MONOMER);
            return std::make_pair(new_state, 1.0);
        }
        else // if( _conf == DOWN_MONOMER ) and dimers
            return std::make_pair(state, 0.0);
    };

    std::vector<Func> functions = {Sz, Splus, Sminus};
    std::vector<std::string> names = {"zz", "-+", "+-"};
    for(int ifun = 0; ifun < names.size(); ifun++){
        start = std::chrono::system_clock::now();
        //<! Get entire hilbert space
        int Nu = this->syms.Nup;
        int Nd = this->syms.Ndown;
        if(ifun == 1) { Nu++; Nd--; }
        if(ifun == 2) { Nu--; Nd++; }
        QHS::Symmetric2DLattice<QOps::particle::fermion> _hilbert_full(this->Lx, this->Ly, Nu, Nd, v_1d<QOps::generic_operator<cpx, PackedLattice2D>>(), 0, 0, 0, 1, 1, true, false);
        u64 dim_tot = _hilbert_full.get_hilbert_space_size();
        
        std::cout << " - - - - - - finished hilnert space for full space for S" + names[ifun] + "_q in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
        arma::mat Sqw(path.size(), omegas.size(), arma::fill::zeros);
        arma::vec static_structure_factor(path.size(), arma::fill::zeros);
        start = std::chrono::system_clock::now();

    #pragma omp parallel for num_threads(outer_threads)
        for(int j = 0; j < path.size(); j++)
        {
            int kx, ky;
            std::tie(kx, ky) = path[j];
            const double qx = constants<double>::two_pi * kx / double(this->Lx);
            const double qy = constants<double>::two_pi * ky / double(this->Ly);
            if(kx == 0 && ky == 0 && ifun == 0){
                double magn_density = ( (Nu - Nd) / 2. ) / double(this->L);
                static_structure_factor(j) = magn_density * magn_density;
                continue;
            }
            // AFTER: precompute once per q-point
            // Outside the k-loop:
            arma::Mat<cpx> phases(this->Lx, this->Ly);
            for(int x = 0; x < this->Lx; x++)
                for(int y = 0; y < this->Ly; y++)
                    phases(x,y) = std::exp(-1i*(qx*x + qy*y)) / double(this->L);
            
            start = std::chrono::system_clock::now();

            SymmetrySectors sym2 = this->syms;
            sym2.z_sym = -this->syms.z_sym;
            sym2.kx_sym = (sym2.kx_sym + kx) % this->Lx;
            sym2.ky_sym = (sym2.ky_sym + ky) % this->Ly;
            //<! px and py change only in \pi sector
            // if( this->Lx % 2? false : kx  == int(this->Lx / 2) ) sym2.px_sym = -this->syms.px_sym;
            // if( this->Ly % 2? false : ky  == int(this->Ly / 2) ) sym2.py_sym = -this->syms.py_sym;
            if( (this->Lx % 2 == 0) && (kx == this->Lx / 2) ) sym2.px_sym = -this->syms.px_sym;
            if( (this->Ly % 2 == 0) && (ky == this->Ly / 2) ) sym2.py_sym = -this->syms.py_sym;
            sym2.Nup = Nu; 
            sym2.Ndown = Nd;
            
            #ifdef USE_SYMMETRIES
                auto Sq_model = std::make_unique<QHS::QHamSolver<BosonicSpinfullQLM>>(this->boundary_conditions, this->Lx, this->Ly, this->tau, this->mass, sym2);
                auto _new_hilbert = Sq_model->get_model_ref().get_hilbert_space();
                auto _new_group = _new_hilbert.get_symmetry_group();
                u64 dim_new = _new_hilbert.get_hilbert_space_size();
                if(dim_new == 0) continue;
                // sym2.print();
                std::cout << " - - - - - - finished generating model in shifted sector in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
            #else
                u64 dim_new = dim;
            #endif
            start = std::chrono::system_clock::now();
            arma::sp_cx_mat Sq(dim_new, dim);
        // #pragma omp parallel for
            #ifdef USE_SYMMETRIES
                for(u64 k = 0; k < dim; k++){
                    auto state = _hilbert(k);          // representative state in sym. basis
                    ui::element_type _norm_k = _hilbert.get_norm(k);  // normalization of sym. state
                    // arma::cx_vec fullstate(dim_tot, arma::fill::zeros);
                    // for(auto& g : _group)
                    {
                        // auto [Gstate, chi_g] = g(state);
                        // u64 idxtot = _hilbert_full.find(Gstate);
                        // fullstate(idxtot) += chi_g / _norm_k / std::sqrt(_group.size());
                        for(int x = 0; x < this->Lx; x++){
                            for(int y = 0; y < this->Ly; y++){
                                auto [new_state, op_val] = functions[ifun](state, x, y);
                                if(std::abs(op_val) < 1e-14) continue;
                                {
                                    // auto [rep, phase_to_rep] = _new_hilbert.find_SEC_representative(new_state);
                                    //<! For operators that obey the same symmetries as the Hamiltonian
                                    auto [idx, phase_to_rep_with_norms] = _new_hilbert.find_matrix_element(new_state, _norm_k);
                                    Sq(idx, k) += op_val * phase_to_rep_with_norms * phases(x,y);
                                }
                            }
                        }
                    }
                    // printSeparated(std::cout, "\t", 20, true, k, some_val, _state_norm  / double(_group.size()), arma::norm(fullstate));
                }
                // Sq = Sq / double(_group.size());
            #else
                for(u64 k = 0; k < dim; k++){
                    auto state = _hilbert(k);          // representative state in sym. basis
                    for(int x = 0; x < this->Lx; x++){
                        for(int y = 0; y < this->Ly; y++){
                            auto [new_state, op_val] = functions[ifun](state, x, y);
                            if(std::abs(op_val) < 1e-14) continue;
                            {
                                u64 idx = _hilbert.find(new_state);
                                if(idx < dim_new)
                                    Sq(idx, k) += op_val * phases(x,y); 
                            }
                        }
                    }
                }
            #endif
            std::cout << std::real((cpx)Sq(0,0)) << std::endl;
            // std::cout << arma::mat(arma::real(Sq)) << std::endl;
            arma::cx_vec new_GS = Sq * GS;
            std::cout << " - - - - - - finished constructing operator between sectors in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
            start = std::chrono::system_clock::now();

            int orbit_s = orbit_size(kx, ky);
            double alfa = std::real( arma::cdot(new_GS, new_GS) * (double)orbit_s );
            static_structure_factor(j) = alfa;

            start = std::chrono::system_clock::now();
            if(dim_new > dim_cut){
                //<! Step 3: Build lanczos matrix starting from Sq|GS>
                #ifdef USE_SYMMETRIES
                arma::sp_cx_mat H = cast_cx_sparse(Sq_model->get_hamiltonian());
                #else
                arma::sp_cx_mat H = cast_cx_sparse(this->ptr_to_model->get_hamiltonian());
                #endif
                auto lancz = lanczos::Lanczos<cpx, converge::energies>(H, 10, 1000, 1e-14, this->seed, 1);
                lancz.build(new_GS);
                int lsteps = lancz.get_lanczossteps();
                auto HL = lancz.get_lanczos_matrix();
                std::cout << " - - - - - - finished lanczos iteration in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end

                //<! Step 4: Continous fraction expansion
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
                //<! Step 3: Exact diagonalization
                #ifdef USE_SYMMETRIES
                    Sq_model->diagonalization();
                    std::cout << " - - - - - - finished ED in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end

                    arma::vec E_Sq = Sq_model->get_eigenvalues();
                    auto V_Sq = Sq_model->get_eigenvectors();
                #else
                    arma::vec E_Sq = E;
                    auto& V_Sq = V;
                #endif
                #pragma omp critical
                {
                    E_Sq.save(arma::hdf5_name(dir_spin + info + ".hdf5", "E/" + names[ifun] + "/j=" + std::to_string(j), arma::hdf5_opts::append));
                }
                // arma::vec E_Sq = E;

                //<! Step 4: Continous fraction expansion
                for(int n = 0; n < E_Sq.size(); n++)
                {
                    #ifdef USE_SYMMETRIES
                    auto state = Sq_model->get_eigenState(n);
                    #else
                    auto state = V.col(n);
                    #endif
                    // auto state = this->ptr_to_model->get_eigenState(n);
                    cpx overlap = dot_prod(state, new_GS);
                // #pragma omp parallel for
                    for(int w = 0; w < omegas.size(); w++)
                    {
                        cpx z = omegas(w) + 1i*eta + E0 - E_Sq(n);
                        Sqw(j, w) += -std::abs( overlap * std::conj(overlap) ) * std::imag( 1. / z ) / constants<double>::pi;
                    }
                    // const auto idx = int( ((E0 - E_Sq(n)) - omegas(0)) / (omegas(1) - omegas(0)) );
                    // Sqw(j, idx) += std::abs( overlap * std::conj(overlap) );
                }

                arma::vec _spectral_fun(omegax.size()-1, arma::fill::zeros);
                arma::vec _element_count(omegax.size()-1, arma::fill::zeros);
                arma::vec _spectral_fun_neg(omegax.size()-1, arma::fill::zeros);
                arma::vec _element_count_neg(omegax.size()-1, arma::fill::zeros);
                const double dw_log = std::log10(omegax[1]) - std::log10(omegax[0]);
                const double w0_log = std::log10(omegax[0]);

                {
                    // std::cout << "V.shape = (" << V.n_rows << ", " << V.n_cols << ")" << std::endl;
                    // std::cout << "Sq.shape = (" << Sq.n_rows << ", " << Sq.n_cols << ")" << std::endl;
                    // std::cout << "V_Sq.shape = (" << V_Sq.n_rows << ", " << V_Sq.n_cols << ")" << std::endl;
                    arma::cx_mat mat_elem = V_Sq.t() * Sq * V;
                    for(int n = 0; n < E.size(); n++){
                        for(int m = 0; m < E_Sq.size(); m++)
                        {
                            double wnm = E_Sq(m) - E(n);
                            if(wnm >= 0){
                                const auto idx = int( (std::log10(wnm) - w0_log) / dw_log);
                                if(idx < omegax.size()-1 && idx >= 0){
                                    _spectral_fun(idx) += std::abs(mat_elem(m, n) * std::conj(mat_elem(m, n))) * (double)orbit_s;
                                    _element_count(idx) += 1;
                                }
                            } else {
                                wnm = -wnm;
                                const auto idx = int( (std::log10(wnm) - w0_log) / dw_log);
                                if(idx < omegax.size()-1 && idx >= 0){
                                    _spectral_fun_neg(idx) += std::abs(mat_elem(m, n) * std::conj(mat_elem(m, n))) * (double)orbit_s;
                                    _element_count_neg(idx) += 1;
                                }
                            }
                        }	
                    }
                    // _spectral_fun = _spectral_fun / _element_count;
                    // _spectral_fun_neg = _spectral_fun_neg / _element_count_neg;
                    std::cout << " - - - - - - finished spectral function of momentum mode in time:" << tim_s(start) << " s - - - - - - " << std::endl; // simulation end
                    #pragma omp critical
                    {
                        _spectral_fun.save(arma::hdf5_name(dir_spin + info + ".hdf5", "Sqw_inf/" + names[ifun] + "/j=" + std::to_string(j) + "/spec_fun", arma::hdf5_opts::append));
                        _spectral_fun_neg.save(arma::hdf5_name(dir_spin + info + ".hdf5", "Sqw_inf/" + names[ifun] + "/j=" + std::to_string(j) + "/spec_fun_neg", arma::hdf5_opts::append));
                        _element_count.save(arma::hdf5_name(dir_spin + info + ".hdf5", "Sqw_inf/" + names[ifun] + "/j=" + std::to_string(j) + "/element_count", arma::hdf5_opts::append));
                        _element_count_neg.save(arma::hdf5_name(dir_spin + info + ".hdf5", "Sqw_inf/" + names[ifun] + "/j=" + std::to_string(j) + "/element_count_neg", arma::hdf5_opts::append));
                    }
                }
            }
            std::cout << " - - - - - - finished continous fraction for (kx,ky) = ( " << kx << "," << ky << " ) in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
            start = std::chrono::system_clock::now();
        }
        #pragma omp critical
        {
            omegas.save(arma::hdf5_name(dir_spin + info + ".hdf5", "omegas", arma::hdf5_opts::append));
            Sqw.save(arma::hdf5_name(dir_spin + info + ".hdf5", "Sqw/" + names[ifun], arma::hdf5_opts::append));
            static_structure_factor.save(arma::hdf5_name(dir_spin + info + ".hdf5", "static_SF/" + names[ifun], arma::hdf5_opts::append));
        }
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
            ",tau=" + to_string_prec(this->tau) + \
            ",m=" + to_string_prec(this->mass);
            #ifdef USE_SYMMETRIES
                if(this->boundary_conditions == 0){
                    name += ",kx=" + std::to_string(this->syms.kx_sym);
                    name += ",ky=" + std::to_string(this->syms.ky_sym);
                }
                #ifndef USE_ONLY_TRANSLATION
                    if(this->syms.Nup == this->syms.Ndown) 
                        name += ",z=" + std::to_string(this->syms.z_sym);
                    if(this->kx_real_sec(this->syms.kx_sym)) 
                        name += ",px=" + std::to_string(this->syms.px_sym);
                    if(this->ky_real_sec(this->syms.ky_sym)) 
                        name += ",py=" + std::to_string(this->syms.py_sym);
                    if( (this->Lx == this->Ly) && (this->syms.kx_sym == this->syms.ky_sym) && this->syms.px_sym == this->syms.py_sym )
                        name += ",pd=" + std::to_string(this->syms.pd_sym);
                #endif
            #endif
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

    auto full_model = std::make_unique<QHS::QHamSolver<BosonicSpinfullQLM>>(this->boundary_conditions, this->Lx, this->Ly, this->tau, this->mass, this->syms, false);
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
        auto symmetric_model = std::make_unique<QHS::QHamSolver<BosonicSpinfullQLM>>(this->boundary_conditions, this->Lx, this->Ly, this->tau, this->mass, sym2, true);
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
    return std::make_unique<QHS::QHamSolver<BosonicSpinfullQLM>>(this->boundary_conditions, this->Lx, this->Ly, this->tau, this->mass, this->syms);
}

/// @brief Reset member unique pointer to model with current parameters in class
void ui::reset_model_pointer(){
    this->ptr_to_model.reset(new QHS::QHamSolver<BosonicSpinfullQLM>(this->boundary_conditions, this->Lx, this->Ly, this->tau, this->mass, this->syms) );
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
    BoseUI::parse_cmd_options(argc, argv);

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

    choosen_option = "-jobid";
    this->set_option(this->jobid, argv, choosen_option);

    //<! SYMMETRIES
    choosen_option = "-Nup";
    this->set_option(this->syms.Nup, argv, choosen_option);
    choosen_option = "-Ndown";
    this->set_option(this->syms.Ndown, argv, choosen_option);
    
    if(this->syms.Ndown < 0)
        this->syms.Ndown = this->L - this->syms.Nup;
    if(this->syms.Nup < 0)
        this->syms.Nup = this->L - this->syms.Ndown;

    choosen_option = "-kx";
    this->set_option(this->syms.kx_sym, argv, choosen_option);
    choosen_option = "-ky";
    this->set_option(this->syms.ky_sym, argv, choosen_option);

    choosen_option = "-px";
    this->set_option(this->syms.px_sym, argv, choosen_option);
    choosen_option = "-py";
    this->set_option(this->syms.py_sym, argv, choosen_option);
    choosen_option = "-pd";
    this->set_option(this->syms.pd_sym, argv, choosen_option);
    choosen_option = "-z";
    this->set_option(this->syms.z_sym, argv, choosen_option);

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
    BoseUI::set_default();
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
    BoseUI::print_help();

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
    BoseUI::printAllOptions();
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

    if(this->boundary_conditions == 0){
        std::cout << "kx  = " << this->syms.kx_sym << std::endl;
        std::cout << "ky  = " << this->syms.ky_sym << std::endl;
    }
    if(this->syms.Nup == this->syms.Ndown) 
        std::cout << "z  = " << this->syms.z_sym << std::endl;
    
    bool is_piX_sector = this->Lx % 2? false : this->syms.kx_sym  == int(this->Lx / 2);
    bool is_piY_sector = this->Ly % 2? false : this->syms.ky_sym  == int(this->Ly / 2);
    bool is_kx_real_sector = this->syms.kx_sym==0 || is_piX_sector;
    bool is_ky_real_sector = this->syms.ky_sym==0 || is_piY_sector;
    if(this->boundary_conditions || (is_kx_real_sector && this->Lx > 1))
        std::cout << "px  = " << this->syms.px_sym << std::endl;
    if(this->boundary_conditions || (is_ky_real_sector && this->Ly > 1))
        std::cout << "py  = " << this->syms.py_sym << std::endl;
    
    if( (this->Lx == this->Ly) && (this->syms.kx_sym == this->syms.ky_sym) && this->syms.px_sym == this->syms.py_sym )
        std::cout << "pd  = " << this->syms.pd_sym << std::endl;
                               
    std::cout << std::endl;
    printSeparated(std::cout, "\t", 16, true, "----------------------------------------------------------------------------------------------------");
}   

};


        // Step 1: Apply S^z_q to fullGS
        //     arma::Col<elem_ty> SqGS(dim_tot, arma::fill::zeros);
        // // #pragma omp parallel for
        //     for(u64 k = 0; k < dim; k++){
        //         if(std::abs(GS(k)) < 1e-14) continue;  // skip zeros for efficiency
        //         auto state = _hilbert(k);          // representative state in sym. basis
        //         cpx _norm = _hilbert.get_norm(k);  // normalization of sym. state
        //         for(auto& g : _group){
        //             auto [Gstate, val] = g(state);
        //             for(int x = 0; x < this->Lx; x++){
        //                 for(int y = 0; y < this->Ly; y++){
        //                     auto [new_state, val2] = functions[ifun](Gstate, x, y);
        //                     if(std::abs(val2) > 0){
        //                         u64 idx = _hilbert_full.find(new_state);
        //                         if(idx >= dim_tot){
        //                             new_state.print_state(); continue;
        //                         }
        //                         SqGS(idx) += GS(k) * val * val2 * phases(x,y) / _norm;  // accumulate with correct weight
        //                     }
        //                 }
        //             }
        //         }
        //     }
        //     SqGS /= std::sqrt(_group.size());

        //     double norm_SqGS = arma::norm(SqGS);
        //     // std::cout << norm_SqGS << std::endl;
        //     if(norm_SqGS < 1e-12){
        //         std::cout << "S^z_q|GS> = 0 for q=(" << kx << "," << ky << ")" << std::endl;
        //         continue;
        //     }
        //     // auto new_GS = SqGS;
        //     std::cout << " - - - - - - finished mapping GS into full space for S" + names[ifun] + "_q in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
        //     start = std::chrono::system_clock::now();

        //     SymmetrySectors sym2 = this->syms;
        //     sym2.z_sym = -this->syms.z_sym;
        //     sym2.kx_sym = (sym2.kx_sym + kx) % this->Lx;
        //     sym2.ky_sym = (sym2.ky_sym + ky) % this->Ly;
        //     //<! px and py change only in \pi sector
        //     // if( this->Lx % 2? false : kx  == int(this->Lx / 2) ) sym2.px_sym = -this->syms.px_sym;
        //     // if( this->Ly % 2? false : ky  == int(this->Ly / 2) ) sym2.py_sym = -this->syms.py_sym;
        //     if( (this->Lx % 2 == 0) && (kx == this->Lx / 2) ) sym2.px_sym = -this->syms.px_sym;
        //     if( (this->Ly % 2 == 0) && (ky == this->Ly / 2) ) sym2.py_sym = -this->syms.py_sym;
        //     sym2.Nup = Nu; 
        //     sym2.Ndown = Nd;
            
            
        //     auto Sq_model = std::make_unique<QHS::QHamSolver<BosonicSpinfullQLM>>(this->boundary_conditions, this->Lx, this->Ly, this->tau, this->mass, sym2);
        //     auto _new_hilbert = Sq_model->get_model_ref().get_hilbert_space();
        //     auto _new_group = _new_hilbert.get_symmetry_group();
        //     u64 dim_new = _new_hilbert.get_hilbert_space_size();
        //     // sym2.print();
        //     std::cout << " - - - - - - finished generating model in shifted sector in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
        //     start = std::chrono::system_clock::now();
        //     arma::sp_cx_mat Sq(dim_new, dim);
        // #pragma omp parallel for
        //     for(u64 k = 0; k < dim; k++){
        //         auto state = _hilbert(k);          // representative state in sym. basis
        //         cpx _norm_k = _hilbert.get_norm(k);  // normalization of sym. state
        //         for(auto& g : _group){
        //             auto [Gstate, chi_g] = g(state);
        //             for(int x = 0; x < this->Lx; x++){
        //                 for(int y = 0; y < this->Ly; y++){
        //                     auto [new_state, op_val] = functions[ifun](Gstate, x, y);
        //                     if(std::abs(op_val) < 1e-14) continue;
        //                     {
        //                         auto [rep, phase_to_rep] = _new_hilbert.find_SEC_representative(new_state);
        //                         u64 idx = _new_hilbert.find(rep);
        //                         if(idx < dim_new){
        //                             cpx _norm_idx = _new_hilbert.get_norm(idx);
        //                         #pragma omp critical
        //                             Sq(idx, k) += chi_g * op_val * std::conj(phase_to_rep) * phases(x,y) / (_norm_k * _norm_idx);  // accumulate with correct weight
        //                         } 
        //                     }
        //                 }
        //             }
        //         }
        //     }
        //     Sq = Sq / double(_group.size());
        //     std::cout << " - - - - - - finished constructing operator between sectors in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
        //     start = std::chrono::system_clock::now();

        //     // //<! Step 2: Cast the state into the symmetric subspace
        //     arma::Col<elem_ty> new_GS(dim_new);
        //     for(u64 k = 0; k < dim_tot; k++)
        //     {
        //         if(std::abs(SqGS(k)) < 1e-14) continue;  // skip zeros for efficiency
        //         auto state = _hilbert_full(k);
        //         auto [new_state, sym_eig] = _new_hilbert.find_SEC_representative(state);
        //         u64 idx = _new_hilbert.find(new_state);
        //         if(idx < dim_new){
        //             new_GS(idx) += SqGS(k) * std::conj(sym_eig) / _new_hilbert.get_norm(idx);
        //         } 
        //         // else 
        //         //     { new_state.print_state(); }
        //     }
        //     new_GS /= std::sqrt(double(_new_group.size()));  // same normalization fix
            
        //     // Sanity check
        //     double norm_newGS = arma::norm(new_GS);
        //     if(norm_newGS < 1e-14){
        //         std::cout << "S^z_q|GS> = 0 in symmetric subspace for q=(" << kx << "," << ky << ")" << std::endl;
        //         continue;
        //     }
        //     // new_GS /= norm_newGS;
        //     std::cout << " - - - - - - finished casting Sq|GS> into new reduced sector in : " << tim_s(start) << " - - - - - - " << std::endl; // simulation end
            
        
            // SqGS /= norm_SqGS;  // normalize before checking symmetry eigenvalues
            // auto Tx = QOps::_translation_X_symmetry(this->Lx, this->Ly, 0, 1);
            // auto Ty = QOps::_translation_Y_symmetry(this->Lx, this->Ly, 0, 1);
            // auto Px = QOps::_mirror_X_symmetry(this->Lx, this->Ly, 1);
            // auto Py = QOps::_mirror_Y_symmetry(this->Lx, this->Ly, 1);
            // // auto D = QOps::_mirror_diagonal_symmetry(this->Lx, this->Ly, 1);
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

            //     // std::tie(new_state, val) = D(state);
            //     // idx = _hilbert_full.find(new_state);
            //     // Pd_SqGS += std::conj(SqGS(k)) * val * SqGS(idx);

            //     if(Nu == Nd){
            //         std::tie(new_state, val) = Z(state);
            //         idx = _hilbert_full.find(new_state);
            //         Z_SqGS += std::conj(SqGS(k)) * val * SqGS(idx);
            //     }
            // }
            // std::cout << "Symmetry eigenvalues of the new state are for momentum q - ( " << kx << "," << ky << " ):" << std::endl;
            // printSeparated(std::cout, "\t", 16, true, "symmetry", "", "Old", "\u2192", "New");
            // printSeparated(std::cout, "\t", 16, true, "k_x", "", this->syms.kx_sym, "\u2192", -std::arg( Tx_SqGS ) * double(this->Lx) / constants<double>::two_pi);
            // printSeparated(std::cout, "\t", 16, true, "k_y", "", this->syms.ky_sym, "\u2192", -std::arg( Ty_SqGS ) * double(this->Ly) / constants<double>::two_pi);
            
            // if(kx == 0 || kx == this->Lx/2.) printSeparated(std::cout, "\t", 16, true, "px", "", this->syms.px_sym, "\u2192", std::real(Px_SqGS));
            // if(ky == 0 || ky == this->Ly/2.) printSeparated(std::cout, "\t", 16, true, "py", "", this->syms.py_sym, "\u2192", std::real(Py_SqGS));
            // // if(kx == ky)                     printSeparated(std::cout, "\t", 16, true, "pd", "", this->syms.pd_sym, "\u2192", std::real(Pd_SqGS));
            // if(Nu == Nd)                     printSeparated(std::cout, "\t", 16, true, "z", "", this->syms.z_sym, "\u2192", std::real(Z_SqGS));
            