#pragma once

#include <numbers>

namespace em::math::constants {

/* --- Simulation Constants --- */

// Divisions per wavelength in memory grid
inline constexpr std::size_t divs{20};

inline constexpr std::size_t pml_cells{12};

inline constexpr float pml_m{3.5f};
inline constexpr float pml_r0{1e-6f};

inline constexpr float pml_kappa_max{5.0f};
inline constexpr float pml_alpha_max{0.05f};
inline constexpr float pm_alpha_max{pml_alpha_max}; // Backward compatibility alias

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
