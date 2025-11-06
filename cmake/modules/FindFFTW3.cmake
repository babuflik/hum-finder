# Find the FFTW library
#
# Usage:
#   find_package(FFTW3 [REQUIRED] [QUIET])
#
# Variables defined by this module:
#   FFTW3_FOUND          - True if FFTW3 was found
#   FFTW3_INCLUDE_DIRS   - Include directories for FFTW3
#   FFTW3_LIBRARIES      - Libraries for FFTW3

find_path(FFTW3_INCLUDE_DIR
    NAMES fftw3.h
    PATHS /usr/include
          /usr/local/include
)

find_library(FFTW3_LIBRARY
    NAMES fftw3
    PATHS /usr/lib
          /usr/local/lib
          /usr/lib/x86_64-linux-gnu
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3
    REQUIRED_VARS
        FFTW3_LIBRARY
        FFTW3_INCLUDE_DIR
)

if(FFTW3_FOUND)
    set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
    set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
endif()

mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARY)