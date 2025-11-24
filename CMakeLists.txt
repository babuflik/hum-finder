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
        src/nl.cpp
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
    tests/test_TDOA_and_ML.cpp
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

# Toolbox test executable
add_executable(test_toolbox tests/test_toolbox.cpp)
target_link_libraries(test_toolbox
    PRIVATE
        humming_core
        GTest::GTest
        GTest::Main
        pthread
)
add_test(NAME ToolboxTest COMMAND test_toolbox)

# NL class test executable
add_executable(test_nl tests/test_nl.cpp)
target_link_libraries(test_nl
    PRIVATE
        humming_core
        GTest::GTest
        GTest::Main
        pthread
)
add_test(NAME NLTest COMMAND test_nl)

# Sig class test executable
add_executable(test_sig tests/test_sig.cpp)
target_link_libraries(test_sig
    PRIVATE
        humming_core
        GTest::GTest
        GTest::Main
        pthread
)
add_test(NAME SigTest COMMAND test_sig)
