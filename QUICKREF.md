# Quick Reference Card

## Verify Everything Works
```bash
cd /home/william/projects/hum-finder/build
make -j4 && ctest --output-on-failure
# Expected: 28/28 tests passing
```

## Key Files to Know

### Documentation
- `CONTINUATION_GUIDE.md` ← **START HERE** for next session
- `TOOLBOX_STATUS.md` - Detailed implementation tracking
- `IMPLEMENTATION_REPORT.md` - What's been completed
- `README.md` - Project overview

### Core Implementation
- `include/pdfclass.h` - Base class for all distributions
- `include/nl.h` + `src/nl.cpp` - Nonlinear filtering (EKF/UKF/PF)
- `include/sig.h` - Signal processing class
- `include/utils_sigsys.h` - Numerical utilities

### Tests
- `tests/test_toolbox.cpp` - 16 tests (distributions + utilities)
- `tests/test_nl.cpp` - 5 tests (filtering)
- `tests/test_sig.cpp` - 5 tests (signal processing)

## What Works Right Now

✅ 10 probability distributions (all header-only)
✅ Extended Kalman Filter
✅ Unscented Kalman Filter  
✅ Particle Filter
✅ Signal statistics and processing
✅ Numerical differentiation (grad/hess/jac)
✅ 28/28 tests passing

## Next Easy Wins

1. **NCChi2Dist** - Copy `chi2dist.h`, add non-centrality parameter (~30 min)
2. **Sig::fft()** - Add FFT using FFTW3 (already linked) (~1 hour)
3. **SensorMod::ls()** - Least squares estimation (~30 min)

## MATLAB Reference
All MATLAB source in: `matlab_resources/classes/*.m`

## Common Commands
```bash
# Build
cd build && make -j4

# Test
ctest --output-on-failure
./test_toolbox  # Run specific test suite

# Clean rebuild
cd build && rm -rf * && cmake .. && make -j4

# Plot results
make plot
make plot3d ELEV=30 AZIM=45
```

## Progress: ~45% Complete
- Distributions: 10/11 ✅
- NL filtering: 4/6 methods ✅
- Utilities: 8 functions ✅
- Sig processing: Basic methods ✅
