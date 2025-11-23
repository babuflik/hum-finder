cmake_minimum_required(VERSION 3.14)
project(humming-localization-cpp LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -------------------------------------------------------------------
# CMake modules path (for FindFFTW3.cmake)
# -------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# -------------------------------------------------------------------
# Dependencies
# -------------------------------------------------------------------
find_package(Eigen3 REQUIRED)
find_package(FFTW3 REQUIRED)        # uses your custom module
find_package(GTest REQUIRED)

# -------------------------------------------------------------------
# Library: humming_core
# -------------------------------------------------------------------
add_library(humming_core STATIC
    src/microphone.cpp
    src/utils.cpp
    src/tdoa_calculator.cpp
    src/localizer.cpp
    src/gcc_phat.cpp
    src/estimators.cpp
)

target_include_directories(humming_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(humming_core PUBLIC
    Eigen3::Eigen
    FFTW3::fftw3
)

# -------------------------------------------------------------------
# Main executable
# -------------------------------------------------------------------
add_executable(humming_localization
    src/main.cpp
)

target_link_libraries(humming_localization PRIVATE
    humming_core
)

# -------------------------------------------------------------------
# Tests
# -------------------------------------------------------------------
# Enable testing
enable_testing()

# Hitta Google Test
find_package(GTest REQUIRED)

# Skapa test-executable
add_executable(run_tests
    tests/test_localizer.cpp
)

# Länka bibliotek
target_link_libraries(run_tests
    PRIVATE
        humming_core
        GTest::GTest
        GTest::Main
        pthread
)

# Lägg till CTest
add_test(NAME LocalizerTest COMMAND run_tests)
