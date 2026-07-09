#pragma once

namespace QOps{
    
    /// @brief Mirror the appropriate local configuration in X (mirroring dimer orientations)
    /// @param v input local configuration
    /// @return mirrored configuration
    inline std::uint8_t mirror_x_label(std::uint8_t v)
    {
        switch (v) {
            case 3: return 5; // right -> left
            case 5: return 3; // left  -> right
            default: return v;
        }
    }
    
    /// @brief Mirror the appropriate local configuration in Y (mirroring dimer orientations)
    /// @param v input local configuration
    /// @return mirrored configuration
    inline std::uint8_t mirror_y_label(std::uint8_t v)
    {
        switch (v) {
            case 2: return 4; // up   -> down
            case 4: return 2; // down -> up
            default: return v;
        }
    }

    /// @brief Mirror the appropriate local configuration along diagonal (mirroring dimer orientations)
    /// @param v input local configuration
    /// @return mirrored configuration
    inline std::uint8_t mirror_diagonal_label(std::uint8_t v)
    {
        switch (v) {
            case 2: return 3; // up    -> right
            case 3: return 2; // right -> up
            case 4: return 5; // down  -> left
            case 5: return 4; // left  -> down
            default: return v;
        }
    }

    /// @brief Mirror the appropriate local configuration along anti-diagonal (mirroring dimer orientations)
    /// @param v input local configuration
    /// @return mirrored configuration
    inline std::uint8_t mirror_anti_diagonal_label(std::uint8_t v)
    {
        // Reflection across anti-diagonal x + y = L - 1:
        // up <-> left, right <-> down
        switch (v) {
            case 2: return 5; // up    -> left
            case 5: return 2; // left  -> up
            case 3: return 4; // right -> down
            case 4: return 3; // down  -> right
            default: return v;
        }
    }

    /// @brief Rotates the appropriate local configuration clockwise in 90` (rotate dimer orientation)
    /// @param v input local configuration
    /// @return rotated configuration
    inline std::uint8_t rotate90_label(std::uint8_t v)
    {
        switch (v) {
            case 2: return 3; // up    -> right
            case 3: return 4; // right -> down
            case 4: return 5; // down  -> left
            case 5: return 2; // left  -> up
            default: return v;
        }
    }
    
    /// @brief Rotates the appropriate local configuration clockwise in 180` (rotate dimer orientation)
    /// @param v input local configuration
    /// @return rotated configuration
    inline std::uint8_t rotate180_label(std::uint8_t v)
    {
        switch (v) {
            case 2: return 4; // up    -> down
            case 4: return 2; // down  -> up
            case 3: return 5; // right -> left
            case 5: return 3; // left  -> right
            default: return v;
        }
    }


    /// @brief Calculate the sign for lattice translation in X direction
    /// @param state input state
    /// @param sx shift considered
    /// @return fermionic sign
    inline int fermionic_sign_translation_x(const PackedLattice2D& state, int sx)
    {
        int parity = 0;
        for (int y = 0; y < state._Ly; ++y) {
            int n_wrap = 0;
            int n_rest = 0;
            // Translation right by sx:
            // old columns x = Lx-sx,...,Lx-1 wrap to the beginning.
            for (int x = 0; x < state._Lx; ++x) {
                if (x >= state._Lx - sx)   n_wrap += int(state._lattice(y, x) < 2);
                else                       n_rest += int(state._lattice(y, x) < 2);
            }
            // Row sign = (-1)^(n_wrap * n_rest)
            parity ^= ((n_wrap & 1) & (n_rest & 1));
        }

        return parity ? -1 : +1;
    }

    inline int monomer_parity_rows(
        const PackedLattice2D& state,
        int y_begin,
        int y_end)
    {
        // Computes monomer parity in rows y_begin, ..., y_end-1
        int parity = 0;

        const int Lx = static_cast<int>(state._Lx);

        for (int y = y_begin; y < y_end; ++y) {
            for (int x = 0; x < Lx; ++x) {
                parity ^= static_cast<int>(state._lattice(y, x) < 2);
            }
        }

        return parity;
    }

    /// @brief Calculate the sign for lattice translation in Y direction
    /// @param state input state
    /// @param sx shift considered
    /// @return fermionic sign
    inline int fermionic_sign_translation_y(
        const PackedLattice2D& state,
        int shiftY)
    {
        const int Ly = static_cast<int>(state._Ly);

        if (Ly <= 1)
            return +1;

        int sy = shiftY % Ly;
        if (sy < 0) sy += Ly;

        if (sy == 0)
            return +1;

        // For y -> y + sy:
        // rows [Ly-sy, Ly) wrap around to the beginning.
        const int wrap_begin = Ly - sy;
        const int wrap_end   = Ly;

        const int rest_begin = 0;
        const int rest_end   = Ly - sy;

        const int wrap_rows = sy;
        const int rest_rows = Ly - sy;

        int first_parity;
        int second_parity;

        // Scan the smaller block first.
        if (wrap_rows <= rest_rows) {
            first_parity = monomer_parity_rows(state, wrap_begin, wrap_end);

            // If first block has even monomer parity, exponent is even.
            if (first_parity == 0)
                return +1;

            second_parity = monomer_parity_rows(state, rest_begin, rest_end);
        } else {
            first_parity = monomer_parity_rows(state, rest_begin, rest_end);

            if (first_parity == 0)
                return +1;

            second_parity = monomer_parity_rows(state, wrap_begin, wrap_end);
        }

        return second_parity ? -1 : +1;
    }

    /// @brief Calculate the sign for lattice translation in both X & Y direction from collected images
    /// @param images input images
    /// @return fermionic sign
    inline int fermionic_sign_from_images(std::vector<int>& images)
    {
        // images are new positions of monomers listed in old canonical order.
        // The fermionic sign is parity of this permutation.
        int parity = 0;

        const int n = static_cast<int>(images.size());
        // For dilute sectors, O(n^2) is usually fastest/simplest.
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                parity ^= (images[i] > images[j]);
            }
        }

        return parity ? -1 : +1;
    }
}