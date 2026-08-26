#include "Kernels.hpp"

namespace em::kernel {

void calculateFutureEField(std::vector<float>& EField,
                           std::vector<float>& HField, const float del_t,
                           const float del_x) {
  const float factor = del_t / (del_x * math::constants::eps0);
  for (auto i{0uz}; i < EField.size() - 1; ++i) {
    EField[i] -= factor * (HField[i] - HField[i - 1]);
  }
}

void calculateFutureHField(std::vector<float>& HField,
                           std::vector<float>& EField, const float del_t,
                           const float del_x) {
  const float factor = del_t / (del_x * math::constants::mu0);
  for (auto i{0uz}; i < HField.size() - 1; ++i) {
    HField[i] -= factor * (EField[i + 1] - EField[i]);
  }
}

} // namespace em::kernel
