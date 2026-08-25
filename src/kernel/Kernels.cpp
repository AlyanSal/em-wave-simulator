#include "Kernels.hpp"

namespace em::kernel {

void calculateFutureEField(std::vector<memory::vector>& EField,
                           std::vector<memory::vector>& HField,
                           const float del_t, const float del_x) {
  for (auto i{0uz}; i <= EField.size(); ++i) {
    EField[i].x_ -= (del_t / (del_x * math::constants::eps0)) *
                    (HField[i + 1].x_ - HField[i - 1].x_);
  }
}

void calculateFutureHField(std::vector<memory::vector>& HField,
                           std::vector<memory::vector>& EField,
                           const float del_t, const float del_x) {
  for (auto i{0uz}; i <= HField.size(); ++i) {
    HField[i].x_ -= (del_t / (del_x * math::constants::mu0)) *
                    (EField[i + 1].x_ - EField[i - 1].x_);
  }
}

} // namespace em::kernel
