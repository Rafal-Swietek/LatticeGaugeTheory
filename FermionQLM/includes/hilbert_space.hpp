#pragma once

#ifndef _HILBERT_BASE
    #include "_base.hpp"
#endif

#ifndef _HILBERT_SYM
#define _HILBERT_SYM
#include "symmetry_generators2D.hpp"
// #include "../../../QHamSolver/include/lattices/_base.hpp"
// #include <functional>

static constexpr std::uint8_t DOWN_MONOMER = 0;
static constexpr std::uint8_t UP_MONOMER   = 1;

// dimer orientation stored on EACH occupied site
static constexpr std::uint8_t DIMER_UP     = 2; // connected to y+1
static constexpr std::uint8_t DIMER_RIGHT  = 3; // connected to x+1
static constexpr std::uint8_t DIMER_DOWN   = 4; // connected to y-1
static constexpr std::uint8_t DIMER_LEFT   = 5; // connected to x-1
static constexpr std::uint8_t EMPTY_SITE   = 6; // empty side for states away from half-filling
static constexpr std::uint8_t SKIP_SITE    = 7; // disallowed config, denote skiped config

namespace QHS{

	#ifdef USE_REAL_SECTORS
		using elem_ty = double;
		#pragma message ("--> Using real symmetry transformation")
	#else
		using elem_ty = cpx;
	#endif

	enum class Direction {
        Right,
        Left,
        Up,
        Down
    };
	/// @brief 
	template <QOps::particle _particle_type = QOps::particle::boson>
	class Symmetric2DLattice : public hilbert_space_base<PackedLattice2D>
	{
		typedef std::pair<PackedLattice2D, elem_ty> return_type;						//<! return type of operator, resulting state and value
		typedef QOps::generic_operator<cpx, PackedLattice2D> symmetry_ty;	//<! type of symemtry generator
		typedef v_1d<symmetry_ty> symmetry_group_ty;						//<! type of symemtry generator array/group

		std::unordered_map<PackedLattice2D, u64, PackedLattice2DHash> index_of;	//<! hash map of states mapping to indices
		symmetry_group_ty _symmetry_group;										//<! vector of symmetry generators
		v_1d<elem_ty> _normalisation;											//<! vector containing state normalizations
		v_1d<int> _sectors;														//<! vector of symmetry sectors
		
		unsigned int _Lx = 1;			//<! number of spin-up particles
		unsigned int _Ly = 1;			//<! number of spin-down particles
		unsigned int _Nup = 1;			//<! number of spin-up particles
		unsigned int _Ndown = 1;		//<! number of spin-down particles
		int _trans_shift_x = 1;			//<! translation shift, i.e. one-, two-,.. sided translation in x direction
		int _trans_shift_y = 1;			//<! translation shift, i.e. one-, two-,.. sided translation in y direction
		int _kx_sector = 0;				//<! quasimomentum symmetry sector
		int _ky_sector = 0;				//<! quasimomentum symmetry sector
		int _boundary_cond = 1;			//<! 1-OBC, 0-PBC
		// bool _include_parity = true;	//<! inclusion of parity generator in input symmetries?
		bool _include_U1U1 = true;		//<! include U(1) x U(1) conservation (spin x particle number)
		bool _real_kx_sector = 1;		//<! is the k_sector real or complex?, i.e. k = 0, pi
		bool _real_ky_sector = 1;		//<! is the k_sector real or complex?, i.e. k = 0, pi
		
		void generate_symmetry_group(const symmetry_group_ty& sym_gen);
		auto get_symmetry_normalization(const PackedLattice2D& state) const -> elem_ty;

		/// @brief 
		virtual void init() override 
			{ this->create_basis(); }
		
    public:
		Symmetric2DLattice() = default;
		Symmetric2DLattice(unsigned int Lx, unsigned int Ly, unsigned int Nu, unsigned int Nd, const symmetry_group_ty& sym_gen, 
			int _BC = 1, int kx_sector = 0, int ky_sector = 0, int trans_shift_x = 1, int trans_shift_y = 1, bool include_U1U1 = true);
		
		virtual void create_basis() override;

		virtual PackedLattice2D operator()(u64 idx) 		const override;
		virtual u64 find(const PackedLattice2D& element) 	const override;

		return_type find_SEC_representative(const PackedLattice2D& state) const;
		std::pair<u64, elem_ty> find_matrix_element(const PackedLattice2D& new_state, elem_ty norm) const;

		auto get_symmetry_group() const { return this-> _symmetry_group; }
		auto get_normalisation()  const { return this-> _normalisation; }
		auto get_norm(u64 idx)	  const { return this-> _normalisation[idx]; }
		auto get_sectors() 		  const { return this-> _sectors; }
		
		u64 count_dimension() const;
		int move(int site, Direction dir) const;
		
		auto set_dimensions(unsigned int Lx, unsigned int Ly, unsigned int Nu, unsigned int Nd, int _BC = 0)
		{
			this->_Lx = Lx;		this->_Ly = Ly;
			this->system_size = Lx * Ly;
			this->_Nup = Nu;	this->_Ndown = Nd;
			this->_boundary_cond = _BC;
		}
		
		arma::SpMat<elem_ty> symmetry_rotation() const { return arma::SpMat<elem_ty>(this->dim, this->dim); }
		// arma::SpMat<elem_ty> symmetry_rotation(const hilbert_space_base& global_hilbert_space) const;
	};


	//<! ---------------------------------------------------------------------------------------------------------------------------------------
	//<! ------------------------------------------------------------------------------------------------------------------------ IMPLEMENTATION

	/// @brief Constructor for point symmetric hilbert space with input generators (combinations are set within)
	/// @param L system size
	/// @param sym_gen list of symmetry generators (i.e. T, P, Z, ...)
	/// @param _BC boundary condition (0-PBC, 1-OBC, ...) (default = 1)
	/// @param k_sec quasimomentum symmetyr sector
	/// @param pos_of_parity position of parity in sym_gen (if not present -> -1, by default -1)
	template <QOps::particle _particle_type>
	inline
	Symmetric2DLattice<_particle_type>::Symmetric2DLattice(unsigned int Lx, unsigned int Ly, unsigned int Nu, unsigned int Nd, const symmetry_group_ty& sym_gen, 
			int _BC, int kx_sector, int ky_sector, int trans_shift_x, int trans_shift_y, bool include_U1U1)
	{
		this->_Lx = Lx;
		this->_Ly = Ly;
		this->system_size = Lx * Ly;
		this->_Nup = Nu;
		this->_Ndown = Nd;
		this->_boundary_cond = _BC;
		this->_trans_shift_x = trans_shift_x;
		this->_trans_shift_y = trans_shift_y;

		if(this->_boundary_cond == 0){ // if PBC
			this->_kx_sector = kx_sector;
			this->_ky_sector = ky_sector;
			bool is_pi_sector = Lx % 2? false : this->_kx_sector  == int(this->_Lx / 2);
			this->_real_kx_sector = (this->_kx_sector == 0) || is_pi_sector;

			is_pi_sector = Ly % 2? false : this->_ky_sector  == int(this->_Ly / 2);
			this->_real_ky_sector = (this->_ky_sector == 0) || is_pi_sector;
		}
		
		// include U(1) x U(1) conservation
		this->_include_U1U1 = include_U1U1;

		this->generate_symmetry_group(sym_gen);
		this->init();
	}

	//<! ------------------------------------------------------------------------------------------ SYMMETRY GROUP
	/// @brief Generate symmetry group with all combinations of symmetry generators
	/// @param sym_gen list of symmetry generators (shall not include translation! )
	template <QOps::particle _particle_type>
	inline
	void Symmetric2DLattice<_particle_type>::generate_symmetry_group(const symmetry_group_ty& sym_gen_in)
	{
		this->_symmetry_group = symmetry_group_ty();
		symmetry_group_ty sym_gen = sym_gen_in;

		// // remove parity for complex quasimomentum sectors
		// if (this->_real_kx_sector){
		// 	#ifndef NODEBUG
		// 		std::cout << "Working in real kx sector. Adding parity" << std::endl;
		// 	#endif
		// 	// sym_gen.erase(sym_gen.begin() + this->_pos_of_parity);
		// }

		// //<! collect all sectors (Z2 sectors are +-1)
		// for(auto& G : sym_gen)
		// 	this->_sectors.emplace_back((int)std::real(chi(G)));

		// // add neutral element
		// this->_symmetry_group.push_back(symmetry_ty(this->system_size));
		
		// // set combinations of available symmetries
		// const int NUM_OF_GENERATORS = sym_gen.size();
		// for (int k = 1; k <= NUM_OF_GENERATORS; k++) {  // loop over all product length
		// 	std::vector<int> bitmask(k, 1);             // K leading 1's
		// 	bitmask.resize(NUM_OF_GENERATORS, 0);   	// N - K trailing 0's
		// 	// std:: cout << k << "\t\t" << bitmask << std::endl;
		// 	// std:: cout << "aa\t" << k << std::endl;
		// 	do {
		// 		symmetry_ty sym_temp(this->system_size);
		// 		for (int i = 0; i < NUM_OF_GENERATORS; ++i) // [0..N-1] integers
		// 			if (bitmask[i] == 1)
		// 				sym_temp %= sym_gen[i];
					
				
		// 		this->_symmetry_group.push_back(sym_temp);
		// 	} while (std::prev_permutation(bitmask.begin(), bitmask.end())); // loop over all combinations with bitmask
		// }
		
		// // set combination of all syms with all translations
		// if (this->_boundary_cond == 0) {
		// 	symmetry_group_ty sym_group_copy = this->_symmetry_group;
		// 	symmetry_ty Tx = QOps::_translation_X_symmetry(this->_Lx, this->_Ly, this->_kx_sector);
		// 	symmetry_ty Ty = QOps::_translation_Y_symmetry(this->_Lx, this->_Ly, this->_ky_sector);
		// 	symmetry_ty Tx_pow = symmetry_ty(this->system_size);

		// 	for (int lx = 0; lx < this->_Lx; ++lx) {
		// 		symmetry_ty T = Tx_pow;   // T = Tx^lx
		// 		for (int ly = 0; ly < this->_Ly; ++ly) {
		// 			// skip lx=ly=0 if already present
		// 			if (!(lx == 0 && ly == 0)) {
		// 				for (const auto& G : sym_group_copy)
		// 					this->_symmetry_group.push_back(T % G);
		// 			}
		// 			T %= Ty;   // T = Tx^lx Ty^(ly+1)
		// 		}
		// 		Tx_pow %= Tx; // Tx^(lx+1)
		// 	}

		// 	//<! append quasimomentum sector
		// 	this->_sectors.emplace_back(this->_kx_sector);
		// 	this->_sectors.emplace_back(this->_ky_sector);
		// }
		// #ifndef NODEBUG
		// 	std::cout << this->_sectors.size() << "\t\t" << this->_symmetry_group.size() << std::endl;
		// #endif
	}

	/// @brief Find super-equivalent class (SEC) representative for given set of states related by symmetry transformations
	/// @param input_state find SEC for given input state
	/// @return SEC
	template <QOps::particle _particle_type>
	inline
	typename Symmetric2DLattice<_particle_type>::return_type
	Symmetric2DLattice<_particle_type>::find_SEC_representative(const PackedLattice2D& input_state) const 
	{
		PackedLattice2D SEC(this->_Lx, this->_Ly);
		SEC.fill_all();
		cpx return_val = 1.0;
		for( auto &G : this->_symmetry_group){
			auto [state, return_value] = G(input_state);
			if (state < SEC) {
				SEC = state;
				return_val = return_value;
			}
		}
		#ifdef USE_REAL_SECTORS
			return std::make_pair(SEC, std::real(return_val));
		#else
			return std::make_pair(SEC, return_val);
		#endif
	}

	/// @brief Calculate normalisation for input state (sum off all symmetry eigenvalues for generators not changing input state)
	/// @param input_state input state
	/// @return normalisation
	template <QOps::particle _particle_type>
	inline
	elem_ty 
	Symmetric2DLattice<_particle_type>::get_symmetry_normalization(const PackedLattice2D& input_state) const 
	{
		elem_ty normalisation = 0.0;
		//for (unsigned int L = 0; l < this->_symmetry_group.size(); l++) {
		for( auto &G : this->_symmetry_group){
			auto [state, return_value] = G(input_state);
			#ifdef USE_REAL_SECTORS
				if (state == input_state)
					normalisation += std::real(return_value);
			#else
				if (state == input_state)
					normalisation += return_value;
			#endif
		}
		if( std::abs(normalisation) < 1e-13 )
			normalisation = 0;
		
		return std::sqrt(normalisation);
	}

	/// @brief Find symmetry generator returning to SEC state
	/// @param new_state input state
	/// @return tuple with SEC state and symmetry return value
	template <QOps::particle _particle_type>
	inline
	typename  std::pair<u64, elem_ty> 
	Symmetric2DLattice<_particle_type>::find_matrix_element(const PackedLattice2D& new_state, elem_ty norm) const
	{
		if( std::abs( get_symmetry_normalization(new_state) ) < 1e-12 )
			return std::make_pair(0, 0.0);

		//<! Look for index in reduced basis (maybe its the SEC already)
		u64 idx = this->find(new_state);
		if (idx < this->dim)
			return std::make_pair(idx, this->_normalisation[idx] / norm);
		
		//<! find SEC for input state
		auto [min, sym_eig] = this->find_SEC_representative(new_state);
		idx = this->find(min);

		#ifndef USE_REAL_SECTORS
			sym_eig = std::conj(sym_eig);
		#endif
		// std::cout << this->_normalisation << std::endl;
		// printSeparated(std::cout, "\t", 16, true, sym_eig, norm, idx, this->_normalisation[idx]);
		// input norm, cause can be used between sectors
		//	return std::make_pair(idx, this->_normalisation[idx] / this->_normalisation[base] * sym_eig);
		
		if (idx < dim)	return std::make_pair(idx, this->_normalisation[idx] / norm * sym_eig);
		else			return std::make_pair(0, 0.0);
			
	}

	//<! ------------------------------------------------------------------------------------------ BASIS CONSTRUCTION
	template <QOps::particle _particle_type>
	inline
	void Symmetric2DLattice<_particle_type>::create_basis()
	{
		this->mapping.clear();
		this->_normalisation.clear();
		this->index_of.clear();
		this->dim = 0;

		const int Lx = static_cast<int>(this->_Lx);
		const int Ly = static_cast<int>(this->_Ly);
		const int V  = Lx * Ly;

		if (Lx <= 0 || Ly <= 0) return;

		// arma::Cube<std::uint8_t> local_configs_reset(this->_Lx, this->_Ly, config);
		// for(int i = 0; i < this->_Lx; i++)
		// 	for(int j = 0; j < this->_Ly; j++)
		// 		for(int k = 0; k < config; k++)
		// 			local_configs(i, j, k) = allowed_configs[k];

		auto idx_to_xy = [Lx](int s)
		{
			return std::make_pair( static_cast<std::uint8_t>(s % Lx), static_cast<std::uint8_t>(s / Lx) );
		};


		//----------------------------------------------------------------------
		// Define possible configurations on next lattice site
		//----------------------------------------------------------------------
		auto get_next_configs = [this, &idx_to_xy](const PackedLattice2D& state, int site, int _nup_left, int _ndown_left)
		{
			auto [x, y] = idx_to_xy(site);

			std::vector<std::uint8_t> next_configs;

			// if (_nempty_left > 0)	next_configs.push_back(EMPTY_SITE);
			if (_ndown_left > 0)	next_configs.push_back(DOWN_MONOMER);
			if (_nup_left > 0)		next_configs.push_back(UP_MONOMER);

			if(_nup_left >= 1 && _ndown_left >= 1){
				// up
				int nei = int(y) + 1;
				if(this->_boundary_cond == 0 && nei >= this->_Ly) nei = 0;
				if (this->_Ly > 1 && nei < this->_Ly && state.get_block(x, nei) == SKIP_SITE)
					next_configs.push_back(DIMER_UP);

				// right
				nei = int(x) + 1;
				if(this->_boundary_cond == 0 && nei >= this->_Lx) nei = 0;
				if (this->_Lx > 1 && nei < this->_Lx && state.get_block(nei, y) == SKIP_SITE)
					next_configs.push_back(DIMER_RIGHT);

				// down
				nei = int(y) - 1;
				if(this->_boundary_cond == 0 && nei < 0) nei = this->_Ly-1;
				if (this->_Ly > 1 && nei >= 0 && state.get_block(x, nei) == SKIP_SITE)
					next_configs.push_back(DIMER_DOWN);

				// left
				nei = int(x) - 1;
				if(this->_boundary_cond == 0 && nei < 0) nei = this->_Lx-1;
				if (this->_Lx > 1 && nei >= 0  && state.get_block(nei, y) == SKIP_SITE)
					next_configs.push_back(DIMER_LEFT);
			}
			return next_configs;
		};
		//----------------------------------------------------------------------
		// Go over all possibilities
		//----------------------------------------------------------------------

		//<! loop over all on-site configurations at sites[j]
		std::function<void(PackedLattice2D&, std::vector<int>&, int, int, int)> loop_configs;
		loop_configs = [this, &idx_to_xy, &loop_configs, &get_next_configs](
				PackedLattice2D& state, std::vector<int>& sites, int _conf, int Nu, int Nd) -> void
		{
			if (Nu < 0 || Nd < 0) return;

			//<! Fill remaining sites away from half-filling
			if (Nu == 0 && Nd == 0) 
			{
				//<! sites now contains remaining sites, if empty then it goes immeadeately to include state
				for (int s : sites) {
					auto [xx, yy] = idx_to_xy(s);
					state.set_block(xx, yy, EMPTY_SITE);
				}
				this->mapping.push_back(state);
				return;
			}

			if (sites.empty()) {
				return;
			}

			auto [x, y] = idx_to_xy(sites[0]);
			int site = sites[0];
			state.set_block(x, y, static_cast<std::uint8_t>(_conf));
			sites.erase(std::remove(sites.begin(), sites.end(), site), sites.end());

			if (_conf > 1 && _conf < 6) {
				Direction dir;
				std::uint8_t partner_conf;

				switch (_conf) {
					case 2:
						dir = Direction::Up;
						partner_conf = DIMER_DOWN;
						break;
					case 3:
						dir = Direction::Right;
						partner_conf = DIMER_LEFT;
						break;
					case 4:
						dir = Direction::Down;
						partner_conf = DIMER_UP;
						break;
					case 5:
						dir = Direction::Left;
						partner_conf = DIMER_RIGHT;
						break;
					default:
						return; // should never happen
				}

				int partner = this->move(site, dir);
				if (partner >= 0) {
					auto [xp, yp] = idx_to_xy(partner);

					state.set_block(xp, yp, partner_conf);
					sites.erase(std::remove(sites.begin(), sites.end(), partner), sites.end());

					if (this->_include_U1U1) { Nu--; Nd--; }

					if (sites.empty()) {
						if (Nu == 0 && Nd == 0)
							this->mapping.push_back(state);
						return;
					}

					std::vector<std::uint8_t> _next_configs = get_next_configs(state, sites[0], Nu, Nd);
					for (auto& _conf2 : _next_configs) {
						PackedLattice2D new_state = state;
						std::vector<int> new_sites = sites;
						loop_configs(new_state, new_sites, _conf2, Nu, Nd);
					}
				}
			}
			else {
				if(this->_include_U1U1){
					if(_conf == DOWN_MONOMER) Nd--;
					else if(_conf == UP_MONOMER) Nu--;
				}

				if (sites.empty()) {
					if (Nu == 0 && Nd == 0)
						this->mapping.push_back(state);
					return;
				}

				std::vector<std::uint8_t> _next_configs = get_next_configs(state, sites[0], Nu, Nd);
				for(auto& _conf2 : _next_configs){
					PackedLattice2D new_state = state;
					std::vector<int> new_sites = sites;
					loop_configs(new_state, new_sites, _conf2, Nu, Nd);
				}
			}
		};

		std::vector<std::uint8_t> allowed_configs;
		if (this->_Ndown > 0) allowed_configs.push_back(DOWN_MONOMER);
		if (this->_Nup > 0)   allowed_configs.push_back(UP_MONOMER);
		if (!this->_include_U1U1 || (this->_Nup > 0 && this->_Ndown > 0)) {
			allowed_configs.push_back(DIMER_UP);
			allowed_configs.push_back(DIMER_RIGHT);
			allowed_configs.push_back(DIMER_DOWN);
			allowed_configs.push_back(DIMER_LEFT);
		}
		for(auto& _conf : allowed_configs)
		{
			//<! create a set of indices to go over
			std::vector<int> sites(this->system_size);
			for (int i = 0; i < static_cast<int>(this->system_size); i++)
				sites[i] = i;
			// for (int i = 0; i < static_cast<int>(this->system_size); i++)
			// 	std::cout << sites[i] << std::endl;
			
			//<! Skip perpendicular dimers for 1D lattices in either X or Y direction
			if(this->_Lx == 1 && (_conf ==  DIMER_LEFT || _conf == DIMER_RIGHT)) continue;
			if(this->_Ly == 1 && (_conf ==  DIMER_UP   || _conf == DIMER_DOWN )) continue;

			//<! create an empty state to operate on
			PackedLattice2D state(this->_Lx, this->_Ly);
			state.fill_all();
			// state.set_block(0, 0, static_cast<std::uint8_t>(_conf) );
			loop_configs(state, sites, _conf, this->_Nup, this->_Ndown);
		}

		//----------------------------------------------------------------------
		// hash map
		//----------------------------------------------------------------------

		this->dim = this->mapping.size();
		std::sort(this->mapping.begin(), this->mapping.end());
		this->index_of.reserve(this->dim);

		for (u64 i = 0; i < this->dim; ++i)
			this->index_of.emplace(this->mapping[i], i);

		_extra_debug(
			for (const PackedLattice2D& elem : this->mapping) {
				elem.print_state();
				std::cout << std::endl;
			}
			std::cout << "Hilbert-space size = " << this->dim << std::endl;
		);
	}

	//<! ------------------------------------------------------------------------------------------ ACCESS TOOLS
	/// @brief Overloaded operator to access elements in hilbert space
	/// @param idx Index of element in hilbert space
	/// @return Element of hilbert space at position 'index'
	template <QOps::particle _particle_type>
	inline
	PackedLattice2D 
	Symmetric2DLattice<_particle_type>::operator()(u64 idx) const
	{ 
		_assert_((idx < this->dim), OUT_OF_MAP);
		return this->mapping[idx]; 
	}


	/// @brief Find index of element in hilbert space
	/// @param element element to find its index
	/// @return index of element 'element'
	template <QOps::particle _particle_type>
	inline
	u64 
	Symmetric2DLattice<_particle_type>::find(const PackedLattice2D& element) const
	{ 
		auto it = this->index_of.find(element);
		return (it != this->index_of.end()) ? it->second : this->dim;
	}


	/// @brief Count Hilbert-space dimension for fixed particle numbers.
	/// @return dimension of Hilbert space
	template <QOps::particle _particle_type>
	inline
	u64 Symmetric2DLattice<_particle_type>::count_dimension() const
	{
		if (this->_Lx <= 0 || this->_Ly <= 0) return 0;

		const int N = this->system_size;

		if (this->_Nup < 0 || this->_Ndown < 0) return 0;
		if (this->_Nup > N || this->_Ndown > N) return 0;
		if (this->_Nup + this->_Ndown > N) return 0;

		std::vector<char> occupied(N, 0);

		auto next_free = [&occupied, N]() -> int {
			for (int i = 0; i < N; ++i)
				if (!occupied[i]) return i;
			return -1;
		};

		auto count_free = [&occupied, N]() -> int {
			int nfree = 0;
			for (int i = 0; i < N; ++i)
				if (!occupied[i]) ++nfree;
			return nfree;
		};

		std::function<u64(int,int)> dfs;
		dfs = [&](int up_left, int down_left) -> u64
		{
			if (up_left < 0 || down_left < 0) return 0;

			int s = next_free();
			if (s == -1) {
				return (up_left == 0 && down_left == 0) ? 1ULL : 0ULL;
			}

			int nfree = count_free();
			if (up_left > nfree || down_left > nfree) return 0;

			u64 count = 0;
			occupied[s] = 1;

			// 1) empty site
			count += dfs(up_left, down_left);

			// 2) up monomer
			count += dfs(up_left - 1, down_left);

			// 3) down monomer
			count += dfs(up_left, down_left - 1);

			// 4) dimers: {right, left, up, down}
			for (Direction dir : {Direction::Right, Direction::Left, Direction::Up, Direction::Down}) {
				int n = this->move(s, dir);
				if (n >= 0 && n != s && !occupied[n]) {
					occupied[n] = 1;
					count += dfs(up_left - 1, down_left - 1);
					occupied[n] = 0;
				}
			}

			occupied[s] = 0;
			return count;
		};

		return dfs(static_cast<int>(this->_Nup), static_cast<int>(this->_Ndown));
	}

	/// @brief Move from a site in a given lattice direction.
	/// @param site linear site index
	/// @param dir direction of motion
	/// @return new site index, or -1 if move is not allowed under OBC
	template <QOps::particle _particle_type>
	inline
	int Symmetric2DLattice<_particle_type>::move(int site, Direction dir) const
	{
		switch (dir) {
			case Direction::Right:
			{
				if (site % this->_Lx == this->_Lx - 1) {
					if (this->_boundary_cond == 0)
						return site + 1 - this->_Lx;   // PBC
					else
						return -1;                     // OBC
				}
				return site + 1;
			}

			case Direction::Left:
			{
				if (site % this->_Lx == 0) {
					if (this->_boundary_cond == 0)
						return site - 1 + this->_Lx;   // PBC
					else
						return -1;                     // OBC
				}
				return site - 1;
			}

			case Direction::Up:
			{
				if (site >= this->system_size - this->_Lx) {
					if (this->_boundary_cond == 0)
						return site % this->_Lx;       // PBC
					else
						return -1;                     // OBC
				}
				return site + this->_Lx;
			}

			case Direction::Down:
			{
				if (site < this->_Lx) {
					if (this->_boundary_cond == 0)
						return site + this->system_size - this->_Lx; // PBC
					else
						return -1;                                   // OBC
				}
				return site - this->_Lx;
			}
		}

		return -1;
	}
}

#endif


		// auto mapping_kernel = [this](const PackedLattice2D& state)
		// {
		// 	this->mapping.push_back(state);
		// 	// auto SEC = std::get<0>(find_SEC_representative(state));
		// 	// if (SEC == state) {
		// 	// 	elem_ty N = get_symmetry_normalization(state);
		// 	// 	if (std::abs(N) > 1e-10) {
		// 	// 		this->mapping.push_back(state);
		// 	// 		this->_normalisation.push_back(N);
		// 	// 	}
		// 	// }
		// };