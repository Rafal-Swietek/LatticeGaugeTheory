struct PackedLattice2DHash
{
    std::size_t operator()(const PackedLattice2D& lat) const noexcept
    {
        std::size_t seed = 0;

        auto hash_combine = [&seed](std::size_t value) {
            seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        };

        hash_combine(std::hash<std::size_t>{}(lat._Lx));
        hash_combine(std::hash<std::size_t>{}(lat._Ly));

        // Hash lattice entries in the same order as comparison:
        // top to bottom, right to left
        for (std::size_t y = lat._Ly; y-- > 0; ) {
            for (std::size_t x = lat._Lx; x-- > 0; ) {
                hash_combine(std::hash<std::uint8_t>{}(lat._lattice(y, x)));
            }
        }

        return seed;
    }
};

// /// @brief Hash functor for PackedLattice2D.
// ///
// /// This allows PackedLattice2D to be used as key in std::unordered_map
// /// and std::unordered_set.
// ///
// /// The hash combines:
// /// - lattice dimensions _Lx, _Ly
// /// - all packed row values
// ///
// /// The combination formula is a standard hash-combine pattern similar to
// /// boost::hash_combine.
// ///
// /// Example usage:
// /// @code
// /// std::unordered_map<PackedLattice2D, int, PackedLattice2DHash> table;
// /// @endcode
// struct PackedLattice2DHash
// {
//     std::size_t operator()(const PackedLattice2D& lat) const noexcept
//     {
//         std::size_t seed = 0;

//         auto hash_combine = [&seed](std::size_t value) {
//             seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
//         };

//         hash_combine(std::hash<std::size_t>{}(lat._Lx));
//         hash_combine(std::hash<std::size_t>{}(lat._Ly));

//         for (u64 row : lat._rows) {
//             hash_combine(std::hash<u64>{}(row));
//         }

//         return seed;
//     }
// };