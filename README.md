# Hum-Finder

Acoustic source localization for finding annoying low-frequency drone sounds in apartments.

## Quick Start

```bash
# Build
mkdir build && cd build
cmake .. && make

# Test
ctest --output-on-failure

# Find that annoying sound
./localize_from_files mic1.wav mic2.wav mic3.wav mic4.wav
```

## What It Does

Pinpoints the 3D location of continuous sounds (hum, buzz, drone) using 4+ microphones:
- ✅ HVAC systems
- ✅ Refrigerators  
- ✅ Electrical transformers
- ✅ Water leaks
- ✅ Any continuous noise

**Accuracy**: 0.1-1.0m depending on source type and microphone placement.

## Documentation

- **`DOCS.md`** - Complete documentation (user guide, API, architecture)
- **`PRACTICAL_DRONE_FINDING_GUIDE.md`** - Step-by-step instructions for end users
- **`FILE_INPUT_GUIDE.md`** - Using audio files
- **`REALTIME_GUIDE.md`** - Live microphone streaming
- **`CONTINUATION_GUIDE.md`** - Developer notes

## Features

- 6 estimation algorithms (LS, WLS, ML, EKF, UKF, CRLB)
- File-based and realtime modes
- Works with any continuous sound (doesn't need to be pure tone)
- Comprehensive test suite (60 tests)
- Validated with real-world scenarios

## Requirements

- CMake 3.14+
- C++17
- Eigen3
- FFTW3
- Google Test (optional, for tests)

## License

See LICENSE file.
