#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <numbers>
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

  // 3-point parabolic / quadratic sub-timestep refinement:
  // delta = 0.5 * (y[-1] - y[+1]) / (y[-1] - 2*y[0] + y[+1])
  const float ym1{waveform[max_idx - 1uz]};
  const float y0{waveform[max_idx]};
  const float yp1{waveform[max_idx + 1uz]};
  const float denom{ym1 - (2.0f * y0) + yp1};
  const float delta{(denom != 0.0f) ? 0.5f * (ym1 - yp1) / denom : 0.0f};

  const float refined_time{(static_cast<float>(max_idx) + delta) * dt};
  const float refined_val{y0 - (0.25f * (ym1 - yp1) * delta)};

  return {refined_time, refined_val};
}

} // namespace

TEST_CASE("FDTD Wave Accuracy: Speed of Light and Decay",
          "[physics][speed_of_light]") {
  const auto FRQ{1e9f};
  const auto WIDTH{3.0f};
  const auto LENGTH{3.0f};
  const auto LARGEST_EPS{1.0f};
  const auto PML_FACTOR{0.2f};
  const auto SIM_TIME{800uz};

  std::function<float(float, float)> perm_func{
      []([[maybe_unused]] float x, [[maybe_unused]] float y) -> float {
        return 1.0f;
      }};
  std::function<float(float, float)> cond_func{
      []([[maybe_unused]] float x, [[maybe_unused]] float y) -> float {
        return 0.0f;
      }};
  std::function<float(float, float)> chi1_func{
      []([[maybe_unused]] float x, [[maybe_unused]] float y) -> float {
        return 0.0f;
      }};
  std::function<float(float, float)> t0_func{
      []([[maybe_unused]] float x, [[maybe_unused]] float y) -> float {
        return 1.0f;
      }};

  em::sim::simulation sim(FRQ, WIDTH, LENGTH, LARGEST_EPS, PML_FACTOR);
  sim.setup_simulation(perm_func, cond_func, chi1_func, t0_func);

  const float x0{LENGTH / 2.0f};
  const float y0{WIDTH / 2.0f};

  sim.add_source(em::source::gaussian, x0, y0);

  const auto [ROWS, COLS, INT_ROWS, INT_COLS, PML_X, PML_Y] = sim.Dims();
  const float dx{(em::math::constants::c0 / FRQ) /
                 static_cast<float>(em::math::constants::divs)};
  const float dt{dx / (2.0f * em::math::constants::c0)};

  const std::size_t j_center{PML_Y + static_cast<std::size_t>(y0 / dx)};

  // Probes along +x axis
  const std::vector<float> axial_radii{0.30f, 0.60f, 0.90f, 1.20f};
  std::vector<std::size_t> axial_indices;
  for (const auto r : axial_radii) {
    const auto i_probe{PML_X + static_cast<std::size_t>((x0 + r) / dx)};
    axial_indices.push_back((j_center * COLS) + i_probe);
  }

  // Probes along diagonal (+x, +y, 45 degrees)
  const std::vector<float> diag_offsets{0.20f, 0.40f, 0.60f, 0.80f};
  std::vector<float> diag_radii;
  std::vector<std::size_t> diag_indices;
  for (const auto d : diag_offsets) {
    const auto i_probe{PML_X + static_cast<std::size_t>((x0 + d) / dx)};
    const auto j_probe{PML_Y + static_cast<std::size_t>((y0 + d) / dx)};
    diag_indices.push_back((j_probe * COLS) + i_probe);
    diag_radii.push_back(d * std::numbers::sqrt2_v<float>);
  }

  const auto [EZ, HX, HY] = sim.getEmData();

  std::vector<std::vector<float>> axial_waveforms(axial_radii.size(),
                                                  std::vector<float>(SIM_TIME));
  std::vector<std::vector<float>> diag_waveforms(diag_radii.size(),
                                                 std::vector<float>(SIM_TIME));

  for (auto step{0uz}; step < SIM_TIME; ++step) {
    sim.step_simulation();
    for (auto k{0uz}; k < axial_radii.size(); ++k) {
      axial_waveforms[k][step] = EZ[axial_indices[k]];
    }
    for (auto k{0uz}; k < diag_radii.size(); ++k) {
      diag_waveforms[k][step] = EZ[diag_indices[k]];
    }
  }

  SECTION("Axial propagation velocity matches physical c0 within 1.0%") {
    const auto [t1, a1] = find_peak(axial_waveforms.front(), dt);
    const auto [t4, a4] = find_peak(axial_waveforms.back(), dt);

    const float delta_r{axial_radii.back() - axial_radii.front()};
    const float delta_t{t4 - t1};
    const float measured_c{delta_r / delta_t};

    REQUIRE_THAT(measured_c,
                 Catch::Matchers::WithinRel(em::math::constants::c0, 0.01f));
  }

  SECTION("Axial propagation velocity matches theoretical FDTD dispersion "
          "within 0.5%") {
    const auto [t1, a1] = find_peak(axial_waveforms.front(), dt);
    const auto [t4, a4] = find_peak(axial_waveforms.back(), dt);

    const float delta_r{axial_radii.back() - axial_radii.front()};
    const float delta_t{t4 - t1};
    const float measured_c{delta_r / delta_t};

    const float s_factor{0.5f};
    const float lambda_div{1.0f /
                           static_cast<float>(em::math::constants::divs)};
    const float theoretical_ratio{
        1.0f - ((em::math::constants::pi * em::math::constants::pi / 24.0f) *
                (1.0f - (s_factor * s_factor)) * (lambda_div * lambda_div))};
    const float expected_c{em::math::constants::c0 * theoretical_ratio};

    REQUIRE_THAT(measured_c, Catch::Matchers::WithinRel(expected_c, 0.005f));
  }

  SECTION("Diagonal propagation velocity matches physical c0 within 1.0%") {
    const auto [t1, a1] = find_peak(diag_waveforms.front(), dt);
    const auto [t4, a4] = find_peak(diag_waveforms.back(), dt);

    const float delta_r{diag_radii.back() - diag_radii.front()};
    const float delta_t{t4 - t1};
    const float measured_c{delta_r / delta_t};

    REQUIRE_THAT(measured_c,
                 Catch::Matchers::WithinRel(em::math::constants::c0, 0.01f));
  }

  SECTION("2D cylindrical wave amplitude decay follows 1/sqrt(r) to < 1.0%") {
    std::vector<float> decay_invariants;
    for (auto k{0uz}; k < axial_radii.size(); ++k) {
      const auto [t_peak, a_peak] = find_peak(axial_waveforms[k], dt);
      decay_invariants.push_back(a_peak * std::sqrt(axial_radii[k]));
    }

    const float mean_far_field{
        (decay_invariants[1] + decay_invariants[2] + decay_invariants[3]) /
        3.0f};
    for (auto i{1uz}; i < 4uz; ++i) {
      const float dev{std::abs(decay_invariants[i] - mean_far_field) /
                      mean_far_field};
      REQUIRE(dev < 0.01f);
    }
  }
}
