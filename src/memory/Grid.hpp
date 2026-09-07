#pragma once

#include <cstdint>
#include <iostream>
#include <vector>

#include "Constants.hpp"

namespace em::memory {

/**
 * Memory handler for the simulation's memory
 */
class grid {
public:
  grid(std::size_t cells);
  grid(std::size_t rows, std::size_t cols);

  [[nodiscard]] auto Rows() const noexcept -> std::size_t { return rows_; }
  [[nodiscard]] auto Cols() const noexcept -> std::size_t { return cols_; }
  [[nodiscard]] auto Cells() const noexcept -> std::size_t { return cells_; }

  auto Ezfield() -> std::vector<float>& { return Ezfield_; }
  auto Hxfield() -> std::vector<float>& { return Hxfield_; }
  auto Hyfield() -> std::vector<float>& { return Hyfield_; }
  auto Dfield() -> std::vector<float>& { return Dfield_; }
  auto Ifield() -> std::vector<float>& { return Ifield_; }
  auto Sfield() -> std::vector<float>& { return Sfield_; }
  auto Permittivity() -> std::vector<float>& { return permittivity_; }
  auto Conductance() -> std::vector<float>& { return conductance_; }
  auto ECoeff() -> std::vector<float>& { return e_den_coeff_; }
  auto ICoeff() -> std::vector<float>& { return i_mult_coeff_; }
  auto SMCoeff() -> std::vector<float>& { return s_mult_coeff_; }
  auto SDCoeff() -> std::vector<float>& { return s_decay_coeff_; }
  auto RealE() -> std::vector<float>& { return real_E_; }
  auto ImagE() -> std::vector<float>& { return imag_E_; }
  auto Chi_1() -> std::vector<float>& { return chi1_; }
  auto T_0() -> std::vector<float>& { return t0_; }

  auto PsiEzx() -> std::vector<float>& { return psi_Ezx_; }
  auto PsiEzy() -> std::vector<float>& { return psi_Ezy_; }
  auto PsiHxy() -> std::vector<float>& { return psi_Hxy_; }
  auto PsiHyx() -> std::vector<float>& { return psi_Hyx_; }

  auto BeX() -> std::vector<float>& { return be_x_; }
  auto AeX() -> std::vector<float>& { return ae_x_; }
  auto InvKeX() -> std::vector<float>& { return inv_ke_x_; }
  auto BhX() -> std::vector<float>& { return bh_x_; }
  auto AhX() -> std::vector<float>& { return ah_x_; }
  auto InvKhX() -> std::vector<float>& { return inv_kh_x_; }

  auto BeY() -> std::vector<float>& { return be_y_; }
  auto AeY() -> std::vector<float>& { return ae_y_; }
  auto InvKeY() -> std::vector<float>& { return inv_ke_y_; }
  auto BhY() -> std::vector<float>& { return bh_y_; }
  auto AhY() -> std::vector<float>& { return ah_y_; }
  auto InvKhY() -> std::vector<float>& { return inv_kh_y_; }

private:
  const std::size_t rows_;
  const std::size_t cols_;
  const std::size_t cells_;

  std::vector<float> Ezfield_;
  std::vector<float> Hxfield_;
  std::vector<float> Hyfield_;
  std::vector<float> Dfield_;
  std::vector<float> Ifield_;
  std::vector<float> Sfield_;
  std::vector<float> permittivity_;
  std::vector<float> conductance_;
  std::vector<float> e_den_coeff_;
  std::vector<float> i_mult_coeff_;
  std::vector<float> s_mult_coeff_;
  std::vector<float> s_decay_coeff_;
  std::vector<float> real_E_;
  std::vector<float> imag_E_;
  std::vector<float> chi1_;
  std::vector<float> t0_;

  std::vector<float> psi_Ezx_;
  std::vector<float> psi_Ezy_;
  std::vector<float> psi_Hxy_;
  std::vector<float> psi_Hyx_;

  std::vector<float> be_x_;
  std::vector<float> ae_x_;
  std::vector<float> inv_ke_x_;
  std::vector<float> bh_x_;
  std::vector<float> ah_x_;
  std::vector<float> inv_kh_x_;

  std::vector<float> be_y_;
  std::vector<float> ae_y_;
  std::vector<float> inv_ke_y_;
  std::vector<float> bh_y_;
  std::vector<float> ah_y_;
  std::vector<float> inv_kh_y_;
};

} // namespace em::memory