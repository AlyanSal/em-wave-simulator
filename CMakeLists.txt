cmake_minimum_required(VERSION 3.20)

project(
  em_wave_simulator
  VERSION 0.1.0
  DESCRIPTION "Electromagnetic Wave Simulation Engine"
  LANGUAGES CXX
)

# Enforce C++23 standard
set(CMAKE_CXX_STANDARD 23)
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
  src/kernel/Kernels.cpp
  src/kernel/Kernels.hpp
  src/math/Constants.hpp
  src/memory/Grid.cpp
  src/memory/Grid.hpp
  src/output/Output.hpp
  src/simulation/Simulation.cpp
  src/simulation/Simulation.hpp
)

target_include_directories(em_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/math>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/memory>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/output>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/kernel>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/simulation>
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
