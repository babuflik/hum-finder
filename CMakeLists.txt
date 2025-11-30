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
# Sensor Fusion Toolbox (subdirectory)
# -------------------------------------------------------------------
set(SENSOR_FUSION_BUILD_EXAMPLES OFF CACHE BOOL "Don't build sensor fusion examples in main project")
add_subdirectory(sensor-fusion)

# -------------------------------------------------------------------
# Library: humming_core
# -------------------------------------------------------------------
add_library(humming_core STATIC
    src/microphone.cpp
    src/utils.cpp
    src/tdoa_calculator.cpp
    src/localizer.cpp
    src/gcc_phat.cpp
    src/audio_file_loader.cpp
)

target_include_directories(humming_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(humming_core PUBLIC
    Eigen3::Eigen
    FFTW3::fftw3
    SensorFusion::sensor_fusion
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
# Real-time localization executable
# -------------------------------------------------------------------
add_executable(realtime_localizer
    src/realtime_localizer.cpp
)

target_link_libraries(realtime_localizer PRIVATE
    humming_core
    pthread
)

# -------------------------------------------------------------------
# Streaming example
# -------------------------------------------------------------------
add_executable(streaming_example
    examples/streaming_example.cpp
)

target_link_libraries(streaming_example PRIVATE
    humming_core
    pthread
)

# -------------------------------------------------------------------
# File-based localization
# -------------------------------------------------------------------
add_executable(localize_from_files
    src/localize_from_files.cpp
)

target_link_libraries(localize_from_files PRIVATE
    humming_core
)

# -------------------------------------------------------------------
# Multi-recording localization (3 mics, multiple recordings)
# -------------------------------------------------------------------
add_executable(localize_multi_recording
    src/localize_multi_recording.cpp
)

target_link_libraries(localize_multi_recording PRIVATE
    humming_core
)

# -------------------------------------------------------------------
# Real-time microphone localization
# -------------------------------------------------------------------
add_executable(localize_realtime
    src/localize_realtime.cpp
)

target_link_libraries(localize_realtime PRIVATE
    humming_core
    pthread
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

# All estimators test executable (comprehensive test suite)
add_executable(test_all_estimators tests/test_all_estimators.cpp)
target_link_libraries(test_all_estimators
    PRIVATE
        humming_core
        GTest::GTest
        GTest::Main
        pthread
)
add_test(NAME AllEstimatorsTest COMMAND test_all_estimators)

# Real-world drone localization test
add_executable(test_real_world_drone tests/test_real_world_drone.cpp)
target_link_libraries(test_real_world_drone
    PRIVATE
        humming_core
        GTest::GTest
        GTest::Main
        pthread
)
add_test(NAME RealWorldDroneTest COMMAND test_real_world_drone)

# Broadband/multi-frequency drone test
add_executable(test_broadband_drone tests/test_broadband_drone.cpp)
target_link_libraries(test_broadband_drone
    PRIVATE
        humming_core
        GTest::GTest
        GTest::Main
        pthread
)
add_test(NAME BroadbandDroneTest COMMAND test_broadband_drone)
