#pragma once

#include <numbers>

namespace em::math::constants {

/* --- Simulation Constants --- */

// Divisions per wavelength in memory grid
inline constexpr std::size_t divs{10};

/* --- Mathematical Constants --- */

// Speed of light in a vacuum
inline constexpr float c0{299'792'458.0};

// permeability of free space
inline constexpr float mu0{4e-7 / std::numbers::pi};

// permittivity of free space
inline constexpr float eps0{1.0 / (mu0 * c0 * c0)};

} // namespace em::math::constants
