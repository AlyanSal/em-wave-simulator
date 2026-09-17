#pragma once

#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#include "Constants.hpp"
#include "StructureOfArrays.hpp"

namespace em::mem {

enum struct MainField : std::uint8_t {
  EzField,
  HxField,
  HyField,
  DzField,
  IzField,
  SzField,
  Permittivity,
  Conductivity,
  inv_E_Coeffs,
  I_Coeffs,
  SM_Coeffs,
  SD_Coeffs,
  Real_E,
  Imag_E,
  Chi_1,
  T_0,
  Psi_Ez_x,
  Psi_Ez_y,
  Psi_Hx_y,
  Psi_Hy_x,
  SIZE,
};

inline constexpr auto MainFieldDefaults =
    [] -> std::array<float, std::to_underlying(MainField::SIZE)> {
  std::array<float, std::to_underlying(MainField::SIZE)> defs{};
  defs[std::to_underlying(MainField::Permittivity)] = 1.0f;
  defs[std::to_underlying(MainField::T_0)] = 1.0f;
  return defs;
};

enum struct PMLX : std::uint8_t { BE, AE, Inv_KE, BH, AH, Inv_KH, SIZE };
enum struct PMLY : std::uint8_t { BE, AE, Inv_KE, BH, AH, Inv_KH, SIZE };

inline constexpr auto PMLDefaults =
    [] -> std::array<float, std::to_underlying(PMLX::SIZE)> {
  std::array<float, std::to_underlying(PMLX::SIZE)> defs{};
  defs[std::to_underlying(PMLX::Inv_KE)] = 1.0f;
  defs[std::to_underlying(PMLX::Inv_KH)] = 1.0f;
  return defs;
};

/**
 * Memory handler for the simulation's memory
 */
class grid {
public:
  explicit grid(std::size_t rows, std::size_t cols);

  [[nodiscard]] auto constexpr Dims() const noexcept
      -> std::array<std::size_t, 3> {
    return {cells_, rows_, cols_};
  }

  [[nodiscard]]
  auto constexpr operator[](MainField array) const noexcept -> float* {
    return mainField_[array];
  }

  [[nodiscard]]
  auto constexpr operator[](PMLX array) const noexcept -> float* {
    return pmlx_[array];
  }

  [[nodiscard]]
  auto constexpr operator[](PMLY array) const noexcept -> float* {
    return pmly_[array];
  }

private:
  const std::size_t rows_;
  const std::size_t cols_;
  const std::size_t cells_;

  const SoA<float, std::to_underlying(MainField::SIZE), alignof(float)>
      mainField_;
  const SoA<float, std::to_underlying(PMLX::SIZE), alignof(float)> pmlx_;
  const SoA<float, std::to_underlying(PMLY::SIZE), alignof(float)> pmly_;
};

} // namespace em::mem