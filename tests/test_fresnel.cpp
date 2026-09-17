#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <vector>

#include "Constants.hpp"
#include "Simulation.hpp"
#include "Sources.hpp"

TEST_CASE("FDTD Material Interface: Fresnel Reflection and Transmission",
          "[physics][fresnel]") {
  const auto FRQ{1e9f};
  const auto WIDTH{2.0f};
  const auto LENGTH{4.0f};
  const auto LARGEST_EPS{4.0f}; // Region 2 has eps=4.0 -> n2=2.0
  const auto PML_FACTOR{0.15f};
  const auto SIM_TIME{1000uz};

  // Vacuum baseline simulation
  std::function<float(float, float)> perm_vac{
      []([[maybe_unused]] float x, [[maybe_unused]] float y) -> float {
        return 1.0f;
      }};
  std::function<float(float, float)> cond_zero{
      []([[maybe_unused]] float x, [[maybe_unused]] float y) -> float {
        return 0.0f;
      }};
  std::function<float(float, float)> chi1_zero{
      []([[maybe_unused]] float x, [[maybe_unused]] float y) -> float {
        return 0.0f;
      }};
  std::function<float(float, float)> t0_one{
      []([[maybe_unused]] float x, [[maybe_unused]] float y) -> float {
        return 1.0f;
      }};

  em::sim::simulation sim_vac(FRQ, WIDTH, LENGTH, LARGEST_EPS, PML_FACTOR);
  sim_vac.setup_simulation(perm_vac, cond_zero, chi1_zero, t0_one);

  // Dielectric interface: eps = 4.0 for x >= 2.0m, eps = 1.0 for x < 2.0m
  std::function<float(float, float)> perm_diel{
      []([[maybe_unused]] float x, [[maybe_unused]] float y) -> float {
        return (x >= 2.0f) ? 4.0f : 1.0f;
      }};

  em::sim::simulation sim_diel(FRQ, WIDTH, LENGTH, LARGEST_EPS, PML_FACTOR);
  sim_diel.setup_simulation(perm_diel, cond_zero, chi1_zero, t0_one);

  const auto [ROWS, COLS, INT_ROWS, INT_COLS, PML_X, PML_Y] = sim_vac.Dims();
  const float dx{(em::math::constants::c0 / (std::sqrt(LARGEST_EPS) * FRQ)) /
                 static_cast<float>(em::math::constants::divs)};

  // Inject a plane wave along vertical line at x = 0.8m
  const float xs{0.8f};
  for (float y{0.0f}; y < WIDTH; y += dx) {
    sim_vac.add_source(em::source::gaussian, xs, y);
    sim_diel.add_source(em::source::gaussian, xs, y);
  }

  // Reflection probe at x = 1.4m (in vacuum region, before interface)
  // Transmission probe at x = 2.6m (in dielectric region, after interface)
  const float x_refl{1.4f};
  const float x_trans{2.6f};
  const float y_mid{WIDTH / 2.0f};

  const std::size_t j_mid{PML_Y + static_cast<std::size_t>(y_mid / dx)};
  const std::size_t i_refl{PML_X + static_cast<std::size_t>(x_refl / dx)};
  const std::size_t i_trans{PML_X + static_cast<std::size_t>(x_trans / dx)};

  const std::size_t idx_refl{(j_mid * COLS) + i_refl};
  const std::size_t idx_trans{(j_mid * COLS) + i_trans};

  const auto [EZ_vac, HX_vac, HY_vac] = sim_vac.getEmData();
  const auto [EZ_diel, HX_diel, HY_diel] = sim_diel.getEmData();

  float max_inc{0.0f};
  float max_refl{0.0f};
  float max_trans{0.0f};

  for (auto step{0uz}; step < SIM_TIME; ++step) {
    sim_vac.step_simulation();
    sim_diel.step_simulation();

    const float ez_inc{EZ_vac[idx_refl]};
    const float ez_diel_probe{EZ_diel[idx_refl]};
    const float ez_refl{ez_diel_probe - ez_inc};
    const float ez_trans{EZ_diel[idx_trans]};

    if (std::abs(ez_inc) > std::abs(max_inc)) {
      max_inc = ez_inc;
    }
    if (std::abs(ez_refl) > std::abs(max_refl)) {
      max_refl = ez_refl;
    }
    if (std::abs(ez_trans) > std::abs(max_trans)) {
      max_trans = ez_trans;
    }
  }

  const float measured_R{max_refl / max_inc};
  const float measured_T{max_trans / max_inc};

  // Exact Fresnel analytical values for n1=1.0, n2=2.0
  // R = (1 - 2) / (1 + 2) = -1/3
  // T = 2(1) / (1 + 2) = 2/3
  constexpr float expected_R{-1.0f / 3.0f};
  constexpr float expected_T{2.0f / 3.0f};

  SECTION(
      "Measured reflection coefficient R matches Fresnel formula within 1.5%") {
    REQUIRE_THAT(measured_R, Catch::Matchers::WithinRel(expected_R, 0.015f));
  }

  SECTION("Measured transmission coefficient T matches Fresnel formula within "
          "1.0%") {
    REQUIRE_THAT(measured_T, Catch::Matchers::WithinRel(expected_T, 0.01f));
  }

  SECTION("Power conservation across interface is preserved (|R|^2 + "
          "(n2/n1)*|T|^2 = 1.0)") {
    constexpr float n1{1.0f};
    constexpr float n2{2.0f};
    const float power_sum{(measured_R * measured_R) +
                          ((n2 / n1) * measured_T * measured_T)};
    REQUIRE_THAT(power_sum, Catch::Matchers::WithinRel(1.0f, 0.015f));
  }
}
