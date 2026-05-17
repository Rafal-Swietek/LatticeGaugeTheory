#pragma once

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

        auto _kernel = [_Lx, _Ly, shift](const PackedLattice2D& state)
        {
            _assert_((state._Lx == _Lx && state._Ly == _Ly), INPUT_STATE_MISMATCH);

            PackedLattice2D new_state(state._Lx, state._Ly);

            const std::size_t Lx = state._Lx;
            const int sx = ((shift % static_cast<int>(Lx)) + static_cast<int>(Lx)) % static_cast<int>(Lx);

            if (sx == 0) {
                new_state._lattice = state._lattice;
            } else {
                new_state._lattice.cols(sx, Lx - 1) = state._lattice.cols(0, Lx - sx - 1);
                new_state._lattice.cols(0, sx - 1)  = state._lattice.cols(Lx - sx, Lx - 1);
            }

            return std::pair(new_state, 1.0);
        };

        return generic_operator<cpx, PackedLattice2D>(_Lx * _Ly, _kernel, std::exp(-im * ksym));
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

        auto _kernel = [_Lx, _Ly, shift](const PackedLattice2D& state)
        {
            _assert_((state._Lx == _Lx && state._Ly == _Ly), INPUT_STATE_MISMATCH);

            PackedLattice2D new_state(state._Lx, state._Ly);

            const std::size_t Ly = state._Ly;
            const int sy = ((shift % static_cast<int>(Ly)) + static_cast<int>(Ly)) % static_cast<int>(Ly);

            if (sy == 0) {
                new_state._lattice = state._lattice;
            } else {
                // move all rows by +sy
                new_state._lattice.rows(sy, Ly - 1) = state._lattice.rows(0, Ly - sy - 1);
                new_state._lattice.rows(0, sy - 1)  = state._lattice.rows(Ly - sy, Ly - 1);
            }

            return std::pair(new_state, 1.0);
        };

        return generic_operator<cpx, PackedLattice2D>(_Lx * _Ly, _kernel, std::exp(-im * ksym));
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