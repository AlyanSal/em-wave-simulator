#pragma once

#include <vector>

#include "Constants.hpp"

namespace em::kernel {

/**
 * Calculates the H-Field vectors for a time step into the future from its
 * current state
 */
void calculateFutureEField(std::vector<float>& EField,
                           std::vector<float>& HField, float del_t,
                           float del_x);

/**
 * Calculates the E-Field vectors for a time step into the future from its
 * current state
 */
void calculateFutureHField(std::vector<float>& HField,
                           std::vector<float>& EField, float del_t,
                           float del_x);

} // namespace em::kernel
