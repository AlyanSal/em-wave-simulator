#include "Grid.hpp"

namespace em::memory {

/**
 * Constructor for the allocator
 */
grid::grid(const std::size_t cells)
    : grid(1uz, cells) {}

grid::grid(const std::size_t rows, const std::size_t cols)
    : rows_{rows},
      cols_{cols},
      cells_{rows * cols},
      Ezfield_(cells_, 0.0f),
      Hxfield_(cells_, 0.0f),
      Hyfield_(cells_, 0.0f),
      Dfield_(cells_, 0.0f),
      Ifield_(cells_, 0.0f),
      Sfield_(cells_, 0.0f),
      permittivity_(cells_, 1.0f),
      conductance_(cells_, 0.0f),
      e_den_coeff_(cells_, 0.0f),
      i_mult_coeff_(cells_, 0.0f),
      s_mult_coeff_(cells_, 0.0f),
      s_decay_coeff_(cells_, 0.0f),
      real_E_(cells_, 0.0f),
      imag_E_(cells_, 0.0f),
      chi1_(cells_, 0.0f),
      t0_(cells_, 1.0f),
      psi_Ezx_(cells_, 0.0f),
      psi_Ezy_(cells_, 0.0f),
      psi_Hxy_(cells_, 0.0f),
      psi_Hyx_(cells_, 0.0f),
      be_x_(cols_, 0.0f),
      ae_x_(cols_, 0.0f),
      inv_ke_x_(cols_, 1.0f),
      bh_x_(cols_, 0.0f),
      ah_x_(cols_, 0.0f),
      inv_kh_x_(cols_, 1.0f),
      be_y_(rows_, 0.0f),
      ae_y_(rows_, 0.0f),
      inv_ke_y_(rows_, 1.0f),
      bh_y_(rows_, 0.0f),
      ah_y_(rows_, 0.0f),
      inv_kh_y_(rows_, 1.0f) {}

} // namespace em::memory