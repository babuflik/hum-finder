# Developer Continuation Guide

Quick reference for developers picking up this project.

## Latest Work (Nov 27, 2025)

**Added broadband source support**: System now validated for multi-frequency sounds (HVAC, water flow, etc.), not just pure tones. See `tests/test_broadband_drone.cpp`.

## Current State

**Test Status**: 7 test suites, 60 tests, 100% passing ✅

```
1. LocalizerTest (2)          - End-to-end localization
2. ToolboxTest (21)           - Distribution & utility functions  
3. NLTest (6)                 - Kalman filters (EKF, UKF, PF)
4. SigTest (5)                - Signal class operations
5. AllEstimatorsTest (12)     - LS, WLS, ML, CRLB comparison
6. RealWorldDroneTest (7)     - Apartment scenarios
7. BroadbandDroneTest (7)     - Multi-frequency validation
```

**Core Features Working**:
- ✅ 6 estimators (LS, WLS, ML, EKF, UKF, CRLB)
- ✅ File input mode (audio files → localization)
- ✅ Realtime mode (streaming microphones)
- ✅ 3D positioning
- ✅ Broadband source handling

## Architecture Overview

```
Audio → GCC-PHAT → TDOA → Estimator → 3D Position
```

**Key Classes**:
- `Localizer`: Orchestrates localization process
- `TDOACalculator`: Extracts time-difference-of-arrival from audio
- `GCC_PHAT`: Cross-correlation algorithm
- `NL`: Kalman filter implementations
- `SensorMod`: Sensor modeling and CRLB calculations

**Estimator Performance** (from `test_all_estimators.cpp`):
- EKF: 0.007m RMSE (best) ⭐
- LS/WLS: 0.012m RMSE
- UKF: 0.023m RMSE
- ML: 0.28m RMSE (grid search)

## Building

```bash
mkdir build && cd build
cmake ..
make

# Run tests
ctest --output-on-failure

# Run specific test
./test_broadband_drone
```

## File Organization

```
src/               - Core implementation
  estimators.cpp   - LS, WLS, ML algorithms
  localizer.cpp    - Main localization class
  nl.cpp           - EKF, UKF, PF implementations
  tdoa_calculator.cpp
  gcc_phat.cpp

include/           - Headers
  estimators.h
  localizer.h
  nl.h
  sensormod.h

tests/             - Test suites
  test_all_estimators.cpp
  test_broadband_drone.cpp
  test_real_world_drone.cpp
  test_localizer.cpp
  test_toolbox.cpp

matlab_resources/  - Original MATLAB toolbox (reference)
```

## Next Steps / TODO

**High Priority**:
- [ ] Connect audio file I/O to complete pipeline
- [ ] Command-line interface for `localize_from_files`
- [ ] Add position visualization/plotting

**Nice to Have**:
- [ ] Particle filter (PF) integration
- [ ] Real audio file examples
- [ ] Performance profiling

**Documentation**:
- [x] Consolidated into DOCS.md
- [x] Real-world validation tests
- [x] Broadband source support documented

## Common Tasks

**Add new estimator**:
1. Declare in `include/estimators.h`
2. Implement in `src/estimators.cpp`
3. Add test in `tests/test_all_estimators.cpp`
4. Update comparison test

**Add new test scenario**:
1. Add test case in appropriate `tests/test_*.cpp`
2. Run: `cd build && cmake .. && make && ctest`

**Debug localization issue**:
1. Check `test_real_world_drone.cpp` for similar scenario
2. Verify mic positions and speed of sound
3. Try different estimators (EKF usually best)
4. Check CRLB - if high, geometry is poor

## Key Insights

1. **EKF is best overall** - 0.007m RMSE, handles time-series well
2. **Broadband works fine** - TOA-based method is frequency-agnostic
3. **Microphone geometry matters** - Avoid linear/planar arrangements
4. **Initial guess affects convergence** - Room center usually good start
5. **Sources outside mic array** - Give direction but lower position accuracy

## References

- MATLAB SigSys Toolbox (in `matlab_resources/`)
- `DOCS.md` - Complete API and architecture docs
- Test suites - Examples of every feature
