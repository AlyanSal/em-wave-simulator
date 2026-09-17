#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <vector>

#include "Constants.hpp"
#include "Simulation.hpp"
#include "Sources.hpp"

TEST_CASE("FDTD Absorbing Boundary: CPML Reflection Suppression",
          "[physics][cpml]") {
  const auto FRQ{1e9f};
  const auto LARGEST_EPS{1.0f};
  const auto PML_FACTOR{0.30f};
  const auto SIM_TIME{600uz};

  // Reference simulation: large domain 4.0m x 4.0m
  // The PML boundaries are 2.0m away from the center source, so any boundary
  // reflections cannot return to the probe before time step 600.
  const auto L_REF{4.0f};
  const auto W_REF{4.0f};

  // Test simulation: small domain 2.0m x 2.0m with CPML boundaries
  // The boundary is 1.0m away from the center source.
  const auto L_TEST{2.0f};
  const auto W_TEST{2.0f};

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

  em::sim::simulation sim_ref(FRQ, W_REF, L_REF, LARGEST_EPS, PML_FACTOR);
  sim_ref.setup_simulation(perm, cond, chi1, t0);

  em::sim::simulation sim_test(FRQ, W_TEST, L_TEST, LARGEST_EPS, PML_FACTOR);
  sim_test.setup_simulation(perm, cond, chi1, t0);

  const float dx{(em::math::constants::c0 / FRQ) /
                 static_cast<float>(em::math::constants::divs)};

  const float x0_ref{L_REF / 2.0f};
  const float y0_ref{W_REF / 2.0f};
  const float x0_test{L_TEST / 2.0f};
  const float y0_test{W_TEST / 2.0f};

  sim_ref.add_source(em::source::gaussian, x0_ref, y0_ref);
  sim_test.add_source(em::source::gaussian, x0_test, y0_test);

  const auto [ROWS_ref, COLS_ref, INT_ROWS_ref, INT_COLS_ref, PML_X_ref,
              PML_Y_ref] = sim_ref.Dims();
  const auto [ROWS_test, COLS_test, INT_ROWS_test, INT_COLS_test, PML_X_test,
              PML_Y_test] = sim_test.Dims();

  const std::size_t idx_ref{
      (PML_Y_ref + static_cast<std::size_t>(y0_ref / dx)) * COLS_ref +
      (PML_X_ref + static_cast<std::size_t>(x0_ref / dx))};
  const std::size_t idx_test{
      (PML_Y_test + static_cast<std::size_t>(y0_test / dx)) * COLS_test +
      (PML_X_test + static_cast<std::size_t>(x0_test / dx))};

  const auto [EZ_ref, HX_ref, HY_ref] = sim_ref.getEmData();
  const auto [EZ_test, HX_test, HY_test] = sim_test.getEmData();

  float peak_signal{0.0f};
  float max_diff_after_bounce{0.0f};

  // The excitation source injects during the first ~100 time steps.
  // After step 120, any difference between the small domain and the large
  // domain is due to reflection from the CPML boundary returning to the center
  // probe.
  constexpr std::size_t bounce_eval_start{120uz};

  for (auto step{0uz}; step < SIM_TIME; ++step) {
    sim_ref.step_simulation();
    sim_test.step_simulation();

    const float ref_val{EZ_ref[idx_ref]};
    const float test_val{EZ_test[idx_test]};
    const float diff{std::abs(test_val - ref_val)};

    if (std::abs(ref_val) > peak_signal) {
      peak_signal = std::abs(ref_val);
    }

    if (step >= bounce_eval_start && diff > max_diff_after_bounce) {
      max_diff_after_bounce = diff;
    }
  }

  const float reflection_ratio{max_diff_after_bounce / peak_signal};
  const float reflection_dB{20.0f * std::log10(reflection_ratio)};

  SECTION("CPML boundary reflection is suppressed below -40 dB") {
    REQUIRE(reflection_dB < -40.0f);
  }

  SECTION("Peak reflected amplitude is less than 1.0% of incident signal") {
    REQUIRE(reflection_ratio < 0.01f);
  }
}
