#pragma once

#include <numbers>

namespace em::math::constants {

/* --- Simulation Constants --- */

// Divisions per wavelength in memory grid
inline constexpr std::size_t divs{20};

/* --- Mathematical Constants --- */

// Speed of light in a vacuum
inline constexpr float c0{299'792'458.0f}; // NOLINT

// Permeability of free space
inline constexpr float mu0{4e-7f * std::numbers::pi_v<float>};

// Permittivity of free space
inline constexpr float eps0{1.0f / (mu0 * c0 * c0)};

// Pi
inline constexpr float pi{std::numbers::pi_v<float>}; // NOLINT

} // namespace em::math::constants
