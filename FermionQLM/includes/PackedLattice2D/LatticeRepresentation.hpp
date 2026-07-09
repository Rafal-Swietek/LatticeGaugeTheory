#pragma once

class PackedLattice2D {
public:
    using value_type = std::uint8_t;

    /// @brief Constructor of the 2D lattice type
    /// @param Lx size of lattice in the x direction
    /// @param Ly size of lattice in the y direction
    PackedLattice2D(std::size_t Lx, std::size_t Ly)
        : _Lx(Lx), _Ly(Ly), _lattice(Ly, Lx, arma::fill::zeros)
    {
        if (block_size == 0)
            throw std::invalid_argument("block_size must be > 0");
    }

    /// @brief Return the local value at coordinates (x,y)
    value_type get_block(std::size_t x, std::size_t y) const {
        check_xy(x, y);
        return _lattice(y, x);
    }

    /// @brief Return the local value at linear site index
    value_type get_block(std::size_t site) const {
        std::size_t x = site % _Lx;
        std::size_t y = site / _Lx;
        check_xy(x, y);
        return _lattice(y, x);
    }

    /// @brief Set the local value at coordinates (x,y)
    void set_block(std::size_t x, std::size_t y, value_type value) {
        check_xy(x, y);
        _lattice(y, x) = value;
    }

    /// @brief Set the local value at linear site index
    void set_block(std::size_t site, value_type value) {
        std::size_t x = site % _Lx;
        std::size_t y = site / _Lx;
        check_xy(x, y);
        _lattice(y, x) = value;
    }

    /// @brief Fill every lattice site with the same local value given by global config
    void fill_all() {
        _lattice.fill(static_cast<value_type>(config));
    }

    /// @brief Fill every lattice site with a specified value
    void fill_all(value_type value) {
        _lattice.fill(value);
    }

    /// @brief Print the current state in symbolic form
    void print_state() const {
        for (std::size_t y = _Ly; y-- > 0; ) {
            for (std::size_t x = 0; x < _Lx; ++x) {
                value_type v = get_block(x, y);
                switch (v) {
                    case 0: std::cout << "\u2193"; break; // ↓
                    case 1: std::cout << "\u2191"; break; // ↑
                    case 2: std::cout << "\u25BC"; break; // ▼
                    case 3: std::cout << "\u25C0"; break; // ◀
                    case 4: std::cout << "\u25B2"; break; // ▲
                    case 5: std::cout << "\u25B6"; break; // ▶
                    case 6: std::cout << "\u2022"; break; // •
                    default: std::cout << int(v); break;
                }
                std::cout << ' ';
            }
            std::cout << '\n';
        }
    }

    std::size_t _Lx, _Ly;
    arma::Mat<value_type> _lattice;

private:
    /// @brief Check whether coordinates are inside lattice bounds
    void check_xy(std::size_t x, std::size_t y) const {
        std::string dims = "(" + std::to_string(_Lx) + "," + std::to_string(_Ly) + ")";
        std::string pos  = "(" + std::to_string(x)   + "," + std::to_string(y)   + ")";
        if (x >= _Lx || y >= _Ly) {
            throw std::out_of_range(
                "PackedLattice2D: index out of range:\n"
                "\t\tLattice dimensions: (Lx,Ly) = " + dims + "\n"
                "\t\tTried position: (x,y) = " + pos
            );
        }
    }
};

#include "operator_overload.hpp"
#include "hash.hpp"

// class PackedLattice2D {
// public:
//     /// @brief Constructor of the 2D lattice type
//     /// @param Lx size of lattice in the x direction
//     /// @param Ly size of lattice in the x direction 
//     PackedLattice2D(std::size_t Lx, std::size_t Ly)
//         : _Lx(Lx), _Ly(Ly), _rows(Ly, 0)
//     {
//         if (block_size == 0)
//             throw std::invalid_argument("block_size must be > 0");
//         if (_Lx * block_size >= 64)
//             throw std::invalid_argument("required system size exceed current bit implementation, where we need Lx * block_size < 64");
//     }

//     // std::size_t Lx() const { return _Lx; }
//     // std::size_t Ly() const { return _Ly; }

//     /// @brief Return the local bits at coordinates x and y
//     /// @param x coordinate on the X-axis
//     /// @param y coordinate on the Y-axis 
//     /// @return local bits at given coordinates
//     std::uint8_t get_block(std::size_t x, std::size_t y) const {
//         check_xy(x, y);
//         return static_cast<std::uint8_t>( (_rows[y] >> (x * block_size)) & block_mask() );
//     }

//     /// @brief Set the local bits at coordinates x and y to specific value
//     /// @param x coordinate on the X-axis
//     /// @param y coordinate on the Y-axis 
//     /// @param value value to be set (0, ..., 2^block_size - 1)
//     void set_block(std::size_t x, std::size_t y, std::uint8_t value) {
//         check_xy(x, y);
//         u64 shift = x * block_size;
//         u64 mask = block_mask() << shift;
//         _rows[y] = (_rows[y] & ~mask) | (( u64(value) & block_mask()) << shift);
//     }

//     /// @brief Set the local bits at given site to specific value
//     /// @param site site index, which determins coordinates x, y
//     /// @param value value to be set (0, ..., 2^block_size - 1)
//     void set_block(std::size_t site, std::uint8_t value) {
//         std::size_t x = site % this->_Lx;
//         std::size_t y = site / this->_Lx;
//         check_xy(x, y);
//         u64 shift = x * block_size;
//         u64 mask = block_mask() << shift;
//         _rows[y] = (_rows[y] & ~mask) | (( u64(value) & block_mask()) << shift);
//     }

//     // const std::vector<u64>& rows() const { return _rows; }
//     // std::vector<u64>& rows() { return _rows; }

//     /// @brief Write out current state to console in bit format
//     void print_state() const {
//         // for (std::size_t i = _Ly; i-- > 0; ) {
//         //     std::cout << boost::dynamic_bitset<>(this->_Lx * block_size, _rows[i]) << std::endl;
//         // }
//         for (std::size_t y = _Ly; y-- > 0; ) {
//             for (std::size_t x = 0; x < _Lx; ++x) {
//                 std::uint8_t v = get_block(x, y);
//                 // std::cout << int(v);
//                 if(v == 0) {
//                     std::cout << "\u2193";
//                 } else if(v == 1) {
//                     std::cout << "\u2191";
//                 } else if(v == 2) {
//                     std::cout << "\u25BC";
//                 } else if(v == 3) {
//                     std::cout << "\u25C0";
//                 } else if(v == 4) {
//                     std::cout << "\u25B2";
//                 } else if(v == 5) {
//                     std::cout << "\u25B6";
//                 } else if(v == 6) {
//                     std::cout << "\u2022";
//                 }
//                 else { std::cout << int(v);}
//                 // boost::dynamic_bitset<> bs(block_size, static_cast<u64>(v));
//                 // for (std::size_t b = block_size; b-- > 0; ) {
//                 //     if(bs[b] == 0) {
//                 //         std::cout << "\u2191";
//                 //     } else if(bs[b] == 1) {
//                 //         std::cout << "\u2193";
//                 //     } else if(bs[b] == 2) {
//                 //         std::cout << "\u25BC";
//                 //     } else if(bs[b] == 3) {
//                 //         std::cout << "\u25C0";
//                 //     } else if(bs[b] == 4) {
//                 //         std::cout << "\u25B2";
//                 //     } else if(bs[b] == 5) {
//                 //         std::cout << "\u25B6";
//                 //     }
//                 //     else { std::cout << bs[b];}
//                 // }
//                 std::cout << ' ';
//             }
//             std::cout << '\n';
//         }
//     }
    

//     /// @brief Fill every lattice site with the same local value `config`.
//     /// 
//     /// More efficient version:
//     /// - constructs one packed row once,
//     /// - copies it to all rows.
//     /// 
//     /// Requires:
//     /// - _Lx * block_size < 64
//     void fill_all()
//     {
//         u64 row_pattern = (1ULL << this->_Lx * block_size) - 1ULL;
//         std::fill(_rows.begin(), _rows.end(), row_pattern);
//     }

//     std::vector<u64> _rows;
//     std::size_t _Lx, _Ly;
// private:
//     /// @brief Check if the input coordinates are within lattice bounds
//     /// @param x coordinate on the X-axis
//     /// @param y  coordinate on the Y-axis
//     void check_xy(std::size_t x, std::size_t y) const {
//         std::string dims = "(" + std::to_string(this->_Lx) + "," + std::to_string(this->_Ly) + ")";
//         std::string pos = "(" + std::to_string(this->_Lx) + "," + std::to_string(this->_Ly) + ")";
//         if (x >= _Lx || y >= _Ly)
//             throw std::out_of_range("PackedLattice2D: index out of range:\n\t\t Lattice dimensions: (Lx,Ly) = " + dims + "\n\t\t Tried position: (x,y) = " + pos);
//     }

//     /// @brief Generate mask on block_size bits at first site
//     /// @return mask
//     static constexpr u64 block_mask()
//         { return (u64(1) << block_size) - 1ULL; }
    
//     /// @brief Generate mask on block_size bits at first site
//     // /// @return mask
//     // u64 row_mask() const {
//     //     return (u64(1) << (_Lx * block_size)) - 1ULL;
//     // }
// };

// #include "operator_overload.hpp"
// #include "hash.hpp"