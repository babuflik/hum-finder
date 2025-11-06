project(humming-localization-cpp)

cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

file(GLOB SOURCES
    src/*.cpp
)

add_executable(humming_localization ${SOURCES})

enable_testing()

file(GLOB TEST_SOURCES
    tests/*.cpp
)

add_executable(run_tests ${TEST_SOURCES})
target_link_libraries(run_tests humming_localization)

add_test(NAME DFTTests COMMAND run_tests)
add_test(NAME LocalizerTests COMMAND run_tests)