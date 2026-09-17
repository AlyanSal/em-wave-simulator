cmake_minimum_required(VERSION 3.20)

project(
  em_wave_simulator
  VERSION 0.1.0
  DESCRIPTION "Electromagnetic Wave Simulation Engine"
  LANGUAGES CXX
)

# Enforce C++26 standard
set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile commands for clangd / clang-tidy / LSP support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------------------------------------------------------------------
# Compiler Warnings & Diagnostics
# ---------------------------------------------------------------------------
add_library(em_compiler_flags INTERFACE)
target_compile_options(em_compiler_flags INTERFACE
  $<$<CXX_COMPILER_ID:MSVC>:
    /W4
    /permissive-
  >
  $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
    -march=native
    -flto
    -fno-math-errno
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
    -Wnon-virtual-dtor
    -Wold-style-cast
    -Wcast-align
    -Wunused
    -Woverloaded-virtual
    -Wconversion
    -Wsign-conversion
  >
)

target_link_options(em_compiler_flags INTERFACE
  $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
    -flto
  >
)

# Optional Clang-Tidy Integration
option(ENABLE_CLANG_TIDY "Run clang-tidy with the compiler" OFF)
if(ENABLE_CLANG_TIDY)
  find_program(CLANG_TIDY_EXE NAMES clang-tidy)
  if(CLANG_TIDY_EXE)
    set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_EXE})
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
  else()
    message(WARNING "ENABLE_CLANG_TIDY was set to ON, but clang-tidy executable was not found.")
  endif()
endif()

# ---------------------------------------------------------------------------
# Core Library
# ---------------------------------------------------------------------------
add_library(em_core
  src/kernel/Kernels.hpp
  src/math/Constants.hpp
  src/memory/Grid.cpp
  src/memory/Grid.hpp
  src/output/BinaryOutput.hpp
  src/output/Output.hpp
  src/simulation/Simulation.cpp
  src/simulation/Simulation.hpp
  src/utilities/Timer.hpp
)

target_include_directories(em_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/kernel>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/math>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/memory>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/output>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/simulation>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/utilities>
)

target_link_libraries(em_core PUBLIC em_compiler_flags)

# ---------------------------------------------------------------------------
# Executable
# ---------------------------------------------------------------------------
add_executable(em_wave_simulator src/main.cpp)

target_link_libraries(em_wave_simulator PRIVATE em_core)

# Set binary output name to match project convention
set_target_properties(em_wave_simulator PROPERTIES
  OUTPUT_NAME "em-wave-simulator"
)

# ---------------------------------------------------------------------------
# Verification & Testing Framework (Catch2 + CTest)
# ---------------------------------------------------------------------------
option(BUILD_TESTING "Build verification and unit tests" ON)
if(BUILD_TESTING)
  enable_testing()
  include(FetchContent)

  set(CATCH_DEVELOPMENT_BUILD OFF CACHE BOOL "" FORCE)
  set(CATCH_INSTALL_DOCS OFF CACHE BOOL "" FORCE)

  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.7.0
    GIT_SHALLOW TRUE
  )
  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
  include(Catch)

  add_executable(em_wave_tests
    tests/test_speed_of_light.cpp
    tests/test_fresnel.cpp
    tests/test_cpml.cpp
    tests/test_convergence.cpp
  )
  target_link_libraries(em_wave_tests PRIVATE em_core Catch2::Catch2WithMain)

  catch_discover_tests(em_wave_tests)
endif()

