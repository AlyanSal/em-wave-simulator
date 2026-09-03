#include "Grid.hpp"

namespace em::memory {

/**
 * Constructor for the allocator
 */
grid::grid(const std::size_t cells)
    : cells_{cells},
      Efield_(cells, 0.0f),
      Hfield_(cells, 0.0f),
      Dfield_(cells, 0.0f),
      Ifield_(cells, 0.0f),
      Sfield_(cells, 0.0f),
      permittivity_(cells, 1.0f),
      conductance_(cells, 0.0f),
      e_den_coeff_(cells, 0.0f),
      i_mult_coeff_(cells, 0.0f),
      s_mult_coeff_(cells, 0.0f),
      s_decay_coeff_(cells, 0.0f),
      real_E_(cells, 0.0f),
      imag_E_(cells, 0.0f),
      chi1_(cells, 0.0f),
      t0_(cells, 1.0f) {}

} // namespace em::memory