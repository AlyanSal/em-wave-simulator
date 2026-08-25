#pragma once

#include <vector>

#include "math/Constants.hpp"
#include "memory/Vector.hpp"

namespace em::kernel {

/**
 * Calculates the H-Field vectors for a time step into the future from its
 * current state
 *
 * TODO: Implement the function
 */
void calculateFutureEField(std::vector<memory::vector>& EField,
                           std::vector<memory::vector>& HField, float del_t,
                           float del_x);

/**
 * Calculates the E-Field vectors for a time step into the future from its
 * current state
 */
void calculateFutrueHField(std::vector<memory::vector>& HField,
                           std::vector<memory::vector>& EField, float del_t,
                           float del_x);

} // namespace em::kernel
