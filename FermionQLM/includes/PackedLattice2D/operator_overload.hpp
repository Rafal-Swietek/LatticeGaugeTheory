#include <cstdint>
#include <cstddef>
#include <stdexcept>

/// @brief Equality comparison for PackedLattice2D.
///
/// Two lattices are equal iff:
/// - they have the same dimensions
/// - and all local entries are equal
inline bool operator==(const PackedLattice2D& a, const PackedLattice2D& b)
{
    if (a._Lx != b._Lx || a._Ly != b._Ly)
        throw std::invalid_argument("PackedLattice2D::operator== dimension mismatch");

    for (std::size_t y = 0; y < a._Ly; ++y)
        for (std::size_t x = 0; x < a._Lx; ++x)
            if (a._lattice(y, x) != b._lattice(y, x))
                return false;

    return true;
}

/// @brief Inequality comparison for PackedLattice2D.
inline bool operator!=(const PackedLattice2D& a, const PackedLattice2D& b)
{
    return !(a == b);
}

/// @brief Lexicographic ordering for PackedLattice2D.
///
/// Significance convention:
/// - least significant site:  (x=0, y=0)
/// - most significant site:   (x=_Lx-1, y=_Ly-1)
///
/// Therefore we compare:
/// 1. rows from top to bottom: y = _Ly-1 ... 0
/// 2. inside each row from right to left: x = _Lx-1 ... 0
///
/// At the first site where the two lattices differ, the comparison of the
/// local entries decides the result.
inline bool operator<(const PackedLattice2D& a, const PackedLattice2D& b)
{
    if (a._Lx != b._Lx || a._Ly != b._Ly)
        throw std::invalid_argument("PackedLattice2D::operator< dimension mismatch");

    for (std::size_t y = a._Ly; y-- > 0; ) {
        for (std::size_t x = a._Lx; x-- > 0; ) {
            auto av = a._lattice(y, x);
            auto bv = b._lattice(y, x);

            // first differing site determines the ordering
            if (av < bv) return true;
            if (av > bv) return false;
        }
    }

    return false;
}

/// @brief Greater-than ordering for PackedLattice2D.
inline bool operator>(const PackedLattice2D& a, const PackedLattice2D& b)
{
    return b < a;
}

/// @brief Less-than-or-equal comparison for PackedLattice2D.
inline bool operator<=(const PackedLattice2D& a, const PackedLattice2D& b)
{
    return !(b < a);
}

/// @brief Greater-than-or-equal comparison for PackedLattice2D.
inline bool operator>=(const PackedLattice2D& a, const PackedLattice2D& b)
{
    return !(a < b);
}

// #include <cstdint>
// #include <cstddef>

// /// @brief Equality comparison for PackedLattice2D.
// /// 
// /// Two lattices are equal iff:
// /// - they have the same dimensions _Lx, _Ly
// /// - and all packed rows are equal
// ///
// /// This is the natural structural equality.
// inline bool operator==(const PackedLattice2D& a, const PackedLattice2D& b)
// {
//     return a._Lx == b._Lx
//         && a._Ly == b._Ly
//         && a._rows == b._rows;
// }

// /// @brief Inequality comparison for PackedLattice2D.
// inline bool operator!=(const PackedLattice2D& a, const PackedLattice2D& b)
// {
//     return !(a == b);
// }

// /// @brief Lexicographic ordering for PackedLattice2D.
// ///
// /// Ordering convention:
// /// - (x=0, y=0) is the least significant site
// /// - (x=Lx-1, y=Ly-1) is the most significant site
// ///
// /// Since each row is packed into one u64 and higher y corresponds to
// /// more significant rows, we compare:
// ///   1. dimensions first
// ///   2. then packed rows from top to bottom, i.e. y = _Ly-1 down to 0
// ///
// /// At the first row where they differ, the smaller packed u64 determines
// /// the ordering.
// ///
// /// This is efficient: O(_Ly) comparisons of u64 values.
// inline bool operator<(const PackedLattice2D& a, const PackedLattice2D& b)
// {
//     if (a._Lx != b._Lx) return a._Lx < b._Lx;
//     if (a._Ly != b._Ly) return a._Ly < b._Ly;

//     for (std::size_t y = a._Ly; y-- > 0; ) {
//         if (a._rows[y] != b._rows[y])
//             return a._rows[y] < b._rows[y];
//     }
//     return false;
// }

// /// @brief Greater-than ordering for PackedLattice2D.
// ///
// /// Uses the same lexicographic convention as operator<.
// inline bool operator>(const PackedLattice2D& a, const PackedLattice2D& b)
// {
//     if (a._Lx != b._Lx) return a._Lx > b._Lx;
//     if (a._Ly != b._Ly) return a._Ly > b._Ly;

//     for (std::size_t y = a._Ly; y-- > 0; ) {
//         if (a._rows[y] != b._rows[y])
//             return a._rows[y] > b._rows[y];
//     }
//     return false;
// }

// /// @brief Less-than-or-equal comparison for PackedLattice2D.
// inline bool operator<=(const PackedLattice2D& a, const PackedLattice2D& b)
// {
//     return !(a > b);
// }

// /// @brief Greater-than-or-equal comparison for PackedLattice2D.
// inline bool operator>=(const PackedLattice2D& a, const PackedLattice2D& b)
// {
//     return !(a < b);
// }
