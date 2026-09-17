#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <vector>

#include "Constants.hpp"
#include "Simulation.hpp"
#include "Sources.hpp"

namespace {

[[nodiscard]] auto
find_peak(const std::vector<float>& waveform,
          const float dt) noexcept -> std::pair<float, float> {
  std::size_t max_idx{0uz};
  float max_val{-1e9f};

  for (auto i{0uz}; i < waveform.size(); ++i) {
    if (waveform[i] > max_val) {
      max_val = waveform[i];
      max_idx = i;
    }
  }

  if (max_idx == 0uz || max_idx >= waveform.size() - 1uz) {
    return {static_cast<float>(max_idx) * dt, max_val};
  }

  // 3-point parabolic sub-timestep refinement:
  const float ym1{waveform[max_idx - 1uz]};
  const float y0{waveform[max_idx]};
  const float yp1{waveform[max_idx + 1uz]};
  const float denom{ym1 - (2.0f * y0) + yp1};
  const float delta{(denom != 0.0f) ? 0.5f * (ym1 - yp1) / denom : 0.0f};

  const float refined_time{(static_cast<float>(max_idx) + delta) * dt};
  const float refined_val{y0 - (0.25f * (ym1 - yp1) * delta)};

  return {refined_time, refined_val};
}

[[nodiscard]] auto measure_phase_velocity(const std::size_t divs) -> float {
  const auto FRQ{1e9f};
  const auto WIDTH{2.0f};
  const auto LENGTH{2.0f};
  const auto LARGEST_EPS{1.0f};
  const auto PML_FACTOR{0.20f};
  const auto SIM_TIME{divs * 40uz};

  std::function<float(float, float)> perm{
      []([[maybe_unused]] float x, [[maybe_unused]] float y) -> float {
        return 1.0f;
      }};
  std::function<float(float, float)> cond{
      []([[maybe_unused]] float x, [[maybe_unused]] float y) -> float {
        return 0.0f;
      }};
  std::function<float(float, float)> chi1{
      []([[maybe_unused]] float x, [[maybe_unused]] float y) -> float {
        return 0.0f;
      }};
  std::function<float(float, float)> t0{
      []([[maybe_unused]] float x, [[maybe_unused]] float y) -> float {
        return 1.0f;
      }};

  em::sim::simulation sim(FRQ, WIDTH, LENGTH, LARGEST_EPS, PML_FACTOR, divs);
  sim.setup_simulation(perm, cond, chi1, t0);

  const float dx{(em::math::constants::c0 / FRQ) / static_cast<float>(divs)};
  const float dt{dx / (2.0f * em::math::constants::c0)};

  const auto [ROWS, COLS, INT_ROWS, INT_COLS, PML_X, PML_Y] = sim.Dims();
  const float x0{LENGTH / 2.0f};
  const float y0{WIDTH / 2.0f};

  sim.add_source(em::source::gaussian, x0, y0);

  const std::size_t j0{PML_Y + static_cast<std::size_t>(y0 / dx)};
  const std::size_t i1{PML_X + static_cast<std::size_t>((x0 + 0.20f) / dx)};
  const std::size_t i2{PML_X + static_cast<std::size_t>((x0 + 0.70f) / dx)};

  const std::size_t idx1{(j0 * COLS) + i1};
  const std::size_t idx2{(j0 * COLS) + i2};

  const auto [EZ, HX, HY] = sim.getEmData();
  std::vector<float> w1(SIM_TIME);
  std::vector<float> w2(SIM_TIME);

  for (auto s{0uz}; s < SIM_TIME; ++s) {
    sim.step_simulation();
    w1[s] = EZ[idx1];
    w2[s] = EZ[idx2];
  }

  const auto [t1, a1] = find_peak(w1, dt);
  const auto [t2, a2] = find_peak(w2, dt);

  const float delta_r{static_cast<float>(i2 - i1) * dx};
  return delta_r / (t2 - t1);
}

} // namespace

TEST_CASE("FDTD Numerical Order: Grid Convergence Rate O(dx^2)",
          "[physics][convergence]") {
  // Discretization levels: coarse (divs=10), medium (divs=20), fine (divs=40)
  const float c10{measure_phase_velocity(10uz)};
  const float c20{measure_phase_velocity(20uz)};
  const float c40{measure_phase_velocity(40uz)};
  constexpr float c0{em::math::constants::c0};

  const float err10{std::abs(c10 - c0)};
  const float err20{std::abs(c20 - c0)};
  const float err40{std::abs(c40 - c0)};

  // Cauchy / Richardson error difference ratio:
  // (c20 - c10) / (c40 - c20) = 2^p
  const float diff12{c20 - c10};
  const float diff23{c40 - c20};
  const float p{std::log2(diff12 / diff23)};

  SECTION(
      "Numerical dispersion error strictly decreases under grid refinement") {
    REQUIRE(err20 < err10);
    REQUIRE(err40 < err20);
  }

  SECTION("Spatial convergence order p is second-order (1.8 <= p <= 2.2)") {
    REQUIRE(p >= 1.80f);
    REQUIRE(p <= 2.20f);
  }
}
