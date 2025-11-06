cmake_minimum_required(VERSION 3.10)
project(humming-localization-cpp)

set(CMAKE_CXX_STANDARD 14)

# Add cmake modules path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# Find required packages
find_package(GTest REQUIRED)
find_package(FFTW3 REQUIRED)

# Include directories
include_directories(${GTEST_INCLUDE_DIRS})
include_directories(${FFTW3_INCLUDE_DIRS})
include_directories(include)

# Core library sources
add_library(humming_core STATIC
    src/microphone.cpp
    src/dft.cpp
    src/localizer.cpp
)

# Link core library with FFTW
target_link_libraries(humming_core PUBLIC ${FFTW3_LIBRARIES})

# Main executable
add_executable(humming_localization src/main.cpp)
target_link_libraries(humming_localization PRIVATE humming_core)

# Testing
enable_testing()

# Create test executable
add_executable(run_tests 
    tests/test_dft.cpp
    tests/test_localizer.cpp
)

# Link test executable
target_link_libraries(run_tests
    PRIVATE
        humming_core
        GTest::GTest
        GTest::Main
        pthread
)