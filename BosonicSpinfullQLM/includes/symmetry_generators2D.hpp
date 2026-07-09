#pragma once
#include "symmetry_helpers.hpp"
//#define NOP
namespace QOps{
    /// @brief Creates translation in X generator for given sector and hilbert space
    /// @param _Lx linear system size in x direction
    /// @param _Ly linear system size in y direction
    /// @param sector symmetry sector (int)
    /// @param shift translation shift in x direction
    /// @return translation generator
    template <particle _particle_type = particle::boson>
    inline
    auto _translation_X_symmetry(unsigned int _Lx, unsigned int _Ly, int sector, int shift = 1)
    {
        _assert_((sector >= 0 && sector < static_cast<int>(_Lx)), NOT_ALLOWED_SYM_SECTOR);
        _assert_((shift > 0 && shift < static_cast<int>(_Lx)), NOT_ALLOWED_SYM_SECTORX);

        const double ksym =
            constants<long double>::two_pi * static_cast<double>(sector) / static_cast<double>(_Lx)
            * static_cast<double>(shift);

        const int sx = ((shift % static_cast<int>(_Lx)) + static_cast<int>(_Lx)) % static_cast<int>(_Lx);
        auto _kernel = [_Lx, _Ly, sx](const PackedLattice2D& state)
        {
            _assert_((state._Lx == _Lx && state._Ly == _Ly), INPUT_STATE_MISMATCH);

            PackedLattice2D new_state(state._Lx, state._Ly);

            if (sx == 0) {
                new_state._lattice = state._lattice;
            } else {
                new_state._lattice.cols(sx, _Lx - 1) = state._lattice.cols(0, _Lx - sx - 1);
                new_state._lattice.cols(0, sx - 1)  = state._lattice.cols(_Lx - sx, _Lx - 1);
            }

            if constexpr(_particle_type == particle::fermion){
                double sign = fermionic_sign_translation_x(state, sx);
                return std::make_pair(new_state, sign);
            } else
                return std::pair(new_state, 1.0);
            // return std::pair(new_state, 1.0);
        };

        return generic_operator<cpx, PackedLattice2D>(_Lx * _Ly, _kernel, std::exp(im * ksym));
    }


    /// @brief Creates translation in Y generator for given sector and hilbert space
    /// @param _Lx linear system size in X direction
    /// @param _Ly linear system size in Y direction
    /// @param sector symmetry sector (int)
    /// @param shift translation shift in y direction
    /// @return translation generator
    template <particle _particle_type = particle::boson>
    inline
    auto _translation_Y_symmetry(unsigned int _Lx, unsigned int _Ly, int sector, int shift = 1)
    {
        _assert_((sector >= 0 && sector < static_cast<int>(_Ly)), NOT_ALLOWED_SYM_SECTOR);
        _assert_((shift > 0 && shift < static_cast<int>(_Ly)), NOT_ALLOWED_SYM_SECTORX);

        const double ksym =
            constants<long double>::two_pi * static_cast<double>(sector) / static_cast<double>(_Ly)
            * static_cast<double>(shift);
        
        const int sy = ((shift % static_cast<int>(_Ly)) + static_cast<int>(_Ly)) % static_cast<int>(_Ly);

        auto _kernel = [_Lx, _Ly, sy](const PackedLattice2D& state)
        {
            _assert_((state._Lx == _Lx && state._Ly == _Ly), INPUT_STATE_MISMATCH);

            PackedLattice2D new_state(state._Lx, state._Ly);

            if (sy == 0) {
                new_state._lattice = state._lattice;
            } else {
                // move all rows by +sy
                new_state._lattice.rows(sy, _Ly - 1) = state._lattice.rows(0, _Ly - sy - 1);
                new_state._lattice.rows(0, sy - 1)  = state._lattice.rows(_Ly - sy, _Ly - 1);
            }

            if constexpr(_particle_type == particle::fermion){
                double sign = fermionic_sign_translation_y(state, sy);
                return std::make_pair(new_state, sign);
            } else
                return std::pair(new_state, 1.0);
            // return std::pair(new_state, 1.0);
        };

        return generic_operator<cpx, PackedLattice2D>(_Lx * _Ly, _kernel, std::exp(im * ksym));
    }

    /// @brief Creates combined translation symmetry T_x^shiftX T_y^shiftY.
    /// @param _Lx lattice size in x
    /// @param _Ly lattice size in y
    /// @param kx_sector momentum sector in x
    /// @param ky_sector momentum sector in y
    /// @param shiftX translation shift in x
    /// @param shiftY translation shift in y
    template <particle _particle_type = particle::boson>
    inline
    auto _translation_XY_symmetry(unsigned int _Lx, unsigned int _Ly,
        int kx_sector, int ky_sector,
        int sx, int sy)
    {
        _assert_((kx_sector >= 0 && kx_sector < static_cast<int>(_Lx)), NOT_ALLOWED_SYM_SECTOR);
        _assert_((ky_sector >= 0 && ky_sector < static_cast<int>(_Ly)), NOT_ALLOWED_SYM_SECTOR);
        _assert_((sx >= 0 && sx < static_cast<int>(_Lx)), NOT_ALLOWED_SYM_SECTORX);
        _assert_((sy >= 0 && sy < static_cast<int>(_Ly)), NOT_ALLOWED_SYM_SECTORY);

        const double phase =
            constants<long double>::two_pi *
            (
                static_cast<double>(kx_sector) * static_cast<double>(sx) / static_cast<double>(_Lx)
                + static_cast<double>(ky_sector) * static_cast<double>(sy) / static_cast<double>(_Ly)
            );

        auto _kernel = [_Lx, _Ly, sx, sy](const PackedLattice2D& state)
        {
            _assert_((state._Lx == _Lx && state._Ly == _Ly), INPUT_STATE_MISMATCH);

            if (sx == 0 && sy == 0) {
                return std::make_pair(state, 1.0);
            }
            else {
                PackedLattice2D new_state(state._Lx, state._Ly);

                if constexpr(_particle_type == particle::fermion)
                {
                    std::vector<int> monomer_images;
                    monomer_images.reserve(_Lx * _Ly);
                    
                    for (std::size_t y = 0; y < _Ly; ++y) {
                        const std::size_t ny = (y + sy) % _Ly;
                        for (std::size_t x = 0; x < _Lx; ++x) {
                            const std::size_t nx = (x + sx) % _Lx;
                            const auto v = state.get_block(x, y);
                            new_state.set_block(nx, ny, v);
                            
                            // If v is monomer, v=0 spin down, v=1 is spin up
                            if (v < 2) {
                                const int new_index = static_cast<int>(ny * state._Lx + nx);
                                monomer_images.push_back(new_index);
                            }
                        }
                    }
                    const int sign = fermionic_sign_from_images(monomer_images);
                    return std::pair(new_state, static_cast<double>(sign));
                } else {
                    // for (std::size_t y = 0; y < state._Ly; ++y) {
                    //     const std::size_t ny = (y + sy) % state._Ly;
                    //     for (std::size_t x = 0; x < state._Lx; ++x) {
                    //         const std::size_t nx = (x + sx) % state._Lx;
                    //         new_state.set_block(nx, ny, state.get_block(x, y));
                    //     }
                    // }
                    arma::Mat<std::uint8_t> tmp(_Ly, _Lx);

                    // Shift columns by sx
                    if (sx == 0) {
                        tmp = state._lattice;
                    } else {
                        tmp.cols(sx, _Lx - 1) = state._lattice.cols(0, _Lx - sx - 1);
                        tmp.cols(0, sx - 1)  = state._lattice.cols(_Lx - sx, _Lx - 1);
                    }

                    // Shift rows by sy
                    if (sy == 0) {
                        new_state._lattice = std::move(tmp);
                    } else {
                        new_state._lattice.rows(sy, _Ly - 1) = tmp.rows(0, _Ly - sy - 1);
                        new_state._lattice.rows(0, sy - 1)  = tmp.rows(_Ly - sy, _Ly - 1);
                    }
                    return std::pair(new_state, 1.0);
                }
            }
            // return std::pair(new_state, 1.0);
        };

        return generic_operator<cpx, PackedLattice2D>(
            _Lx * _Ly, _kernel,
            std::exp(im * phase)
        );
    }
    
    /// @brief Creates spin inversion generator
    /// @param _Lx linear system size in X direction
    /// @param _Ly linear system size in Y direction
    /// @param sector symmetry sector z = +-1
    /// @return transformed state with \up <-> \down, while dimers unchanged
    template <particle _particle_type = particle::boson>
    inline
    auto _spin_inversion_symmetry(unsigned int _Lx, unsigned int _Ly, int sector)
    {
        _assert_((sector == -1 || sector == 1), NOT_ALLOWED_SYM_SECTOR);

        auto _kernel = [_Lx, _Ly](const PackedLattice2D& state)
        {
            _assert_((state._Lx == _Lx && state._Ly == _Ly), INPUT_STATE_MISMATCH);

            PackedLattice2D new_state(state._Lx, state._Ly);

            // Copy whole lattice first
            new_state._lattice = state._lattice;

            // Armadillo matrices are contiguous in memory.
            // This is faster than get_block/set_block loops.
            auto* ptr = new_state._lattice.memptr();
            const arma::uword n_elem = new_state._lattice.n_elem;

            //<! v = 0: 0 ^ 1 = 1
            //<! v = 1: 1 ^ 1 = 0
            //<! v >= 2: v ^ 0 = v
            for (arma::uword i = 0; i < n_elem; ++i)
                ptr[i] ^= static_cast<std::uint8_t>(ptr[i] < 2);

            return std::pair(new_state, 1.0);
        };

        return generic_operator<cpx, PackedLattice2D>(_Lx * _Ly, _kernel, double(sector));
    }

    /// @brief Creates Mrror in X generator for given sector and hilbert space
    /// @param _Lx linear system size in X direction
    /// @param _Ly linear system size in Y direction
    /// @param sector symmetry sector Px = +-1
    /// @return transformed state with (x,y) -> (Lx - x - 1, y)
    template <particle _particle_type = particle::boson>
    inline
    auto _mirror_X_symmetry(unsigned int _Lx, unsigned int _Ly, int sector)
    {
        _assert_((sector == -1 || sector == 1), NOT_ALLOWED_SYM_SECTOR);
        static constexpr std::array<std::uint8_t, 8> label_map_mirrorX = {0, 1, 2, 5, 4, 3, 6, 7};

        auto _kernel = [_Lx, _Ly](const PackedLattice2D& state)
        {
            _assert_((state._Lx == _Lx && state._Ly == _Ly), INPUT_STATE_MISMATCH);

            PackedLattice2D new_state(state._Lx, state._Ly);

            if constexpr(_particle_type == particle::fermion){
                std::vector<int> monomer_images;
                monomer_images.reserve(_Lx * _Ly);
                for (std::size_t y = 0; y < _Ly; ++y) {
                    for (std::size_t x = 0; x < _Lx; ++x) {
                        const std::size_t nx = _Lx - 1 - x;
                        const std::size_t ny = y;
                        const std::uint8_t v = state._lattice(y, x);
                        new_state._lattice(ny, nx) = label_map_mirrorX[v];

                        if (v < 2) {
                            const int new_index = static_cast<int>(ny * _Lx + nx);
                            monomer_images.push_back(new_index);
                        }
                    }
                }
                const int sign = fermionic_sign_from_images(monomer_images);
                return std::pair(new_state, static_cast<double>(sign));
            } else { 
                new_state._lattice = arma::fliplr(state._lattice);
                for (auto& v : new_state._lattice)
                    v = label_map_mirrorX[v];

                return std::pair(new_state, 1.0);
            }
        };

        return generic_operator<cpx, PackedLattice2D>(_Lx * _Ly, _kernel, double(sector));
    }

    /// @brief Creates Mrror in Y generator for given sector and hilbert space
    /// @param _Lx linear system size in X direction
    /// @param _Ly linear system size in Y direction
    /// @param sector symmetry sector Py = +-1
    /// @return transformed state with (x,y) -> (x, Ly - y - 1)
    template <particle _particle_type = particle::boson>
    inline
    auto _mirror_Y_symmetry(unsigned int _Lx, unsigned int _Ly, int sector)
    {
        _assert_((sector == -1 || sector == 1), NOT_ALLOWED_SYM_SECTOR);

        static constexpr std::array<std::uint8_t, 8> label_map_mirrorY =
            {0, 1, 4, 3, 2, 5, 6, 7};

        auto _kernel = [_Lx, _Ly](const PackedLattice2D& state)
        {
            _assert_((state._Lx == _Lx && state._Ly == _Ly), INPUT_STATE_MISMATCH);

            PackedLattice2D new_state(state._Lx, state._Ly);
            if constexpr (_particle_type == particle::fermion) {
                
                std::vector<int> monomer_images;
                monomer_images.reserve(state._Lx * state._Ly);
                for (std::size_t y = 0; y < state._Ly; ++y) {
                    for (std::size_t x = 0; x < state._Lx; ++x) {
                        const std::size_t nx = x;
                        const std::size_t ny = state._Ly - 1 - y;
                        const std::uint8_t v = state._lattice(y, x);
                        new_state._lattice(ny, nx) = label_map_mirrorY[v];

                        // Only monomers contribute to fermionic sign.
                        if (v < 2) {
                            const int new_index = static_cast<int>(ny * state._Lx + nx);
                            monomer_images.push_back(new_index);
                        }
                    }
                }
                const int sign = fermionic_sign_from_images(monomer_images);
                return std::pair(new_state, static_cast<double>(sign));
            }
            else {
                new_state._lattice = arma::flipud(state._lattice);
                for (auto& v : new_state._lattice)
                    v = label_map_mirrorY[v];
                return std::pair(new_state, 1.0);
            }
        };

        return generic_operator<cpx, PackedLattice2D>(_Lx * _Ly, _kernel, double(sector) );
    }

    /// @brief Creates Mirror along the diagonal for given sector and hilbert space
    /// @param _Lx linear system size in X direction
    /// @param _Ly linear system size in Y direction (has to be equal for diagonal mirror!)
    /// @param sector symmetry sector D = +-1
    /// @return transformed state with (x,y) -> (y, x)
    template <particle _particle_type = particle::boson>
    inline
    auto _mirror_diagonal_symmetry(unsigned int _Lx, unsigned int _Ly, int sector)
    {
        _assert_((sector == -1 || sector == 1), NOT_ALLOWED_SYM_SECTOR);
        _assert_((_Lx == _Ly), INPUT_STATE_MISMATCH);

        static constexpr std::array<std::uint8_t, 8> label_map_mirrorDiag =
            {0, 1, 3, 2, 5, 4, 6, 7};

        auto _kernel = [_Lx, _Ly](const PackedLattice2D& state)
        {
            _assert_((state._Lx == _Lx && state._Ly == _Ly), INPUT_STATE_MISMATCH);
            _assert_((state._Lx == state._Ly), INPUT_STATE_MISMATCH);
            PackedLattice2D new_state(state._Lx, state._Ly);

            if constexpr (_particle_type == particle::fermion) {
                std::vector<int> monomer_images;
                monomer_images.reserve(state._Lx * state._Ly);
                const std::size_t L = state._Lx;
                for (std::size_t y = 0; y < L; ++y) {
                    for (std::size_t x = 0; x < L; ++x) {
                        const std::size_t nx = y;
                        const std::size_t ny = x;
                        const std::uint8_t v = state._lattice(y, x);
                        new_state._lattice(ny, nx) = label_map_mirrorDiag[v];

                        // Only monomers contribute to fermionic sign.
                        if (v < 2) {
                            const int new_index =static_cast<int>(ny * state._Lx + nx);
                            monomer_images.push_back(new_index);
                        }
                    }
                }
                const int sign = fermionic_sign_from_images(monomer_images);
                return std::pair(new_state, static_cast<double>(sign));
            }
            else {
                new_state._lattice = state._lattice.t();
                for (auto& v : new_state._lattice)
                    v = label_map_mirrorDiag[v];
                return std::pair(new_state, 1.0);
            }
        };

        return generic_operator<cpx, PackedLattice2D>(_Lx * _Ly, _kernel, double(sector));
    }

    /// @brief Creates Mirror along the anti-diagonal for given sector and hilbert space
    /// @param _Lx linear system size in X direction
    /// @param _Ly linear system size in Y direction (has to be equal for (anti)diagonal mirror!)
    /// @param sector symmetry sector AD = +-1
    /// @return transformed state with (x,y) -> (Ly - 1 - y, Lx - 1 - x)
    template <particle _particle_type = particle::boson>
    inline
    auto _mirror_anti_diagonal_symmetry(unsigned int _Lx, unsigned int _Ly, int sector)
    {
        _assert_((sector == -1 || sector == 1), NOT_ALLOWED_SYM_SECTOR);
        _assert_((_Lx == _Ly), INPUT_STATE_MISMATCH);

        auto _kernel = [_Lx, _Ly](const PackedLattice2D& state)
        {
            _assert_((state._Lx == _Lx && state._Ly == _Ly), INPUT_STATE_MISMATCH);
            _assert_((state._Lx == state._Ly), INPUT_STATE_MISMATCH);

            PackedLattice2D new_state(state._Lx, state._Ly);

            const std::size_t L = state._Lx;

            for (std::size_t y = 0; y < L; ++y) {
                for (std::size_t x = 0; x < L; ++x) {
                    std::size_t nx = L - 1 - y;
                    std::size_t ny = L - 1 - x;

                    new_state._lattice(ny, nx) =
                        mirror_anti_diagonal_label(state._lattice(y, x));
                }
            }

            return std::pair(new_state, 1.0);
        };

        return generic_operator<cpx, PackedLattice2D>(_Lx * _Ly, _kernel, double(sector));
    }
    
    /// @brief Creates rotation symmetry in \pi/2 angle for given sector and hilbert space
    /// @param _Lx linear system size in X direction
    /// @param _Ly linear system size in Y direction (has to be equal to _Lx for rotation symmetry!)
    /// @param sector symmetry sector C^4 = 1
    /// @return transformed state with (x,y) -> (y, x)
    template <particle _particle_type = particle::boson>
    inline
    auto _rotation_90_symmetry(unsigned int _Lx, unsigned int _Ly, int sector)
    {
        _assert_((sector >= 0 && sector < 4), NOT_ALLOWED_SYM_SECTOR);
        _assert_((_Lx == _Ly), NOT_ALLOWED_OPERATION);

         auto _kernel = [_Lx, _Ly](const PackedLattice2D& state)
        {
            _assert_((state._Lx == _Lx && state._Ly == _Ly), INPUT_STATE_MISMATCH);
            _assert_((state._Lx == state._Ly), INPUT_STATE_MISMATCH);

            PackedLattice2D new_state(state._Lx, state._Ly);

            const std::size_t L = state._Lx;

            for (std::size_t y = 0; y < L; ++y) {
                for (std::size_t x = 0; x < L; ++x) {
                    std::size_t nx = y;
                    std::size_t ny = L - 1 - x;

                    new_state._lattice(ny, nx) =
                        rotate90_label(state._lattice(y, x));
                }
            }
            return std::pair(new_state, 1.0);
        };
        cpx C4sym;
        switch(sector){
            case 0: C4sym = 1.0;
            case 1: C4sym = +im;
            case 2: C4sym = -1.0;
            case 3: C4sym = -im;
        }
        return generic_operator<cpx, PackedLattice2D>(_Lx * _Ly, _kernel, C4sym);
    }

    /// @brief Creates rotation symmetry by pi angle for given sector and Hilbert space.
    /// @param _Lx linear system size in X direction
    /// @param _Ly linear system size in Y direction
    /// @param sector symmetry sector for C2^2 = 1; sector = 0 gives +1, sector = 1 gives -1
    /// @return C2 rotation symmetry generator
    template <particle _particle_type = particle::boson>
    inline
    auto _rotation_180_symmetry(unsigned int _Lx, unsigned int _Ly, int sector)
    {
        _assert_((sector == -1 || sector == 1), NOT_ALLOWED_SYM_SECTOR);
        _assert_((_Lx == _Ly), NOT_ALLOWED_OPERATION);

        auto _kernel = [_Lx, _Ly](const PackedLattice2D& state)
        {
            _assert_((state._Lx == _Lx && state._Ly == _Ly), INPUT_STATE_MISMATCH);

            PackedLattice2D new_state(state._Lx, state._Ly);

            const std::size_t Lx = state._Lx;
            const std::size_t Ly = state._Ly;

            for (std::size_t y = 0; y < Ly; ++y) {
                for (std::size_t x = 0; x < Lx; ++x) {
                    std::size_t nx = Lx - 1 - x;
                    std::size_t ny = Ly - 1 - y;

                    new_state._lattice(ny, nx) =
                        rotate180_label(state._lattice(y, x));
                }
            }

            return std::pair(new_state, 1.0);
        };

        return generic_operator<cpx, PackedLattice2D>(_Lx * _Ly, _kernel, double(sector));
    }
    
};


    // /// @brief Creates translation in X generator for given sector and hilbert space
    // /// @param Lx linear system size in x direction
    // /// @param sector symmetry sector (int) --> quasimomentum sector is calculated within
    // /// @return translation geerator
	// template <particle _particle_type = particle::boson>
    // inline
    // auto _translation_X_symmetry(unsigned int _Lx, unsigned int _Ly, int sector, int shift = 1)
    // { 
    //     _assert_((sector >= 0 && sector < _Lx), NOT_ALLOWED_SYM_SECTOR);
    //     _assert_((shift > 0 && shift < _Lx), NOT_ALLOWED_SYM_SECTORX);
    //     const double ksym = constants<long double>::two_pi * (double)sector / _Lx * double(shift);    // eigenvalue of T^l is e^(-ikl) with l = shift
    //     auto _kernel = [_Lx, _Ly, shift](const PackedLattice2D& state)
    //             {
    //                 _assert_( (state._Lx == _Lx && state._Ly == _Ly), INPUT_STATE_MISMATCH);

    //                 PackedLattice2D new_state(state._Lx, state._Ly);
    //                 const std::size_t row_bits   = state._Lx * block_size;
    //                 const std::size_t shift_bits = shift * block_size;

    //                 const u64 row_mask =
    //                     (row_bits == 64) ? ~u64(0) : ((u64(1) << row_bits) - 1ULL);

    //                 const u64 low_mask =
    //                     (shift_bits == 64) ? ~u64(0) : ((u64(1) << shift_bits) - 1ULL);

    //                 const std::size_t rotate = row_bits - shift_bits;
    //                 for (std::size_t y = 0; y < state._Ly; ++y) 
    //                 {
    //                     const u64 row = state._rows[y] & row_mask;
    //                     const u64 wrapped = row & low_mask;
    //                     const u64 rest = row >> shift_bits;
    //                     new_state._rows[y] = (rest | (wrapped << rotate)) & row_mask;
    //                 }

    //                 return std::pair(new_state, 1.0);
    //             };
    //     return generic_operator<cpx, PackedLattice2D>(_Lx*_Ly, _kernel, std::exp(-im * ksym));
    // }

    // /// @brief Creates translation in Y generator for given sector and hilbert space
    // /// @param Ly linear system size in Y direction
    // /// @param sector symmetry sector (int) --> quasimomentum sector is calculated within
    // /// @return translation geerator
	// template <particle _particle_type = particle::boson>
    // inline
    // auto _translation_Y_symmetry(unsigned int _Lx, unsigned int _Ly, int sector, int shift = 1)
    // { 
    //     _assert_((sector >= 0 && sector < _Ly), NOT_ALLOWED_SYM_SECTOR);
    //     _assert_((shift > 0 && shift < _Ly), NOT_ALLOWED_SYM_SECTORX);
    //     const double ksym = constants<long double>::two_pi * (double)sector / _Ly * double(shift);    // eigenvalue of T^l is e^(-ikl) with l = shift
    //     auto _kernel = [_Lx, _Ly, shift](const PackedLattice2D& state)
    //             {
    //                 _assert_( (state._Lx == _Lx && state._Ly == _Ly), INPUT_STATE_MISMATCH);

    //                 PackedLattice2D new_state(state._Lx, state._Ly);
    //                 for (std::size_t y = 0; y < state._Ly; ++y)
    //                     new_state._rows[(y + shift) % state._Ly] = state._rows[y];

    //                 return std::pair(new_state, 1.0);
    //             };
    //     return generic_operator<cpx, PackedLattice2D>(_Lx*_Ly, _kernel, std::exp(-im * ksym));
    // }