# Implementation Continuation Guide

**Last Updated**: November 24, 2025  
**Current Progress**: ~45% of MATLAB SigSys Toolbox

## Quick Status Check

Run these commands to verify current state:
```bash
cd /home/william/projects/hum-finder/build
make -j4                      # Should build cleanly
ctest --output-on-failure     # Should show 28/28 tests passing
```

## What's Currently Working

### ✅ Fully Implemented (28 tests passing)

1. **10 Probability Distributions** (all header-only in `include/`)
   - NDist, UDist, ExpDist, Chi2Dist, GammaDist, GMDist
   - TDist, BetaDist, EmpDist, LogNDist
   - All inherit from `PDFClass` base class
   - Tests: `test_toolbox.cpp` (16 tests)

2. **NL Class** - Nonlinear filtering framework
   - Files: `include/nl.h`, `src/nl.cpp`
   - Methods: `simulate()`, `ekf()`, `ukf()`, `pf()`
   - Tests: `tests/test_nl.cpp` (5 tests)

3. **Sig Class** - Signal processing
   - File: `include/sig.h` (header-only)
   - Methods: `mean_y()`, `var_y()`, `cov_y()`, `extract()`, `downsample()`
   - Tests: `tests/test_sig.cpp` (5 tests)

4. **Utility Functions** in `include/utils_sigsys.h`
   - `numgrad()`, `numhess()`, `numjac()` - numerical differentiation
   - `sqrtcov()`, `iscov()` - covariance utilities
   - `confellipse()`, `getwindow()`, `condition()`

## Next Implementation Priorities

### Priority 1: Complete Remaining Distribution (Easy)

**Task**: Implement `NCChi2Dist` (Non-central chi-squared)
- **Reference**: `matlab_resources/classes/ncchi2dist.m`
- **Template**: Copy `chi2dist.h` and modify
- **Effort**: ~30 minutes
- **Key difference**: Additional non-centrality parameter λ

### Priority 2: LTI System Classes (Medium)

**Task**: Implement Linear Time-Invariant system framework
- **Files to create**:
  - `include/ltimod.h` - Base class for LTI systems
  - `include/lss.h` - State-space models: x' = Ax + Bu, y = Cx + Du
  - `include/ltf.h` - Transfer functions
- **Reference**: `matlab_resources/classes/ltimod.m`, `lss.m`, `ltf.m`
- **Methods needed**:
  - `simulate()` - Time/frequency response
  - `bode()`, `nyquist()` - Frequency analysis
  - `step()`, `impulse()` - Standard responses
- **Effort**: ~2-3 hours
- **Note**: Can reuse NL class patterns

### Priority 3: Signal Processing Extensions (Medium)

**Task**: Add FFT/spectral methods to Sig class
- **Methods to add**:
  - `fft()`, `ifft()` - Using FFTW3 (already linked)
  - `psd()` - Power spectral density
  - `filter(b, a)` - Digital filtering
- **Effort**: ~1-2 hours
- **Benefit**: Makes Sig class much more useful

### Priority 4: Enhance SensorMod (Easy-Medium)

**Task**: Add estimation methods to SensorMod
- **Methods to implement**:
  - `ls()` - Least squares
  - `wls()` - Weighted least squares  
  - `ml()` - Maximum likelihood
  - `crlb()` - Cramér-Rao lower bound
- **Reference**: `matlab_resources/classes/sensormod.m`
- **Effort**: ~1 hour
- **Note**: Can use `numjac()` for Jacobians

## Implementation Patterns

### Pattern 1: Adding a Distribution Class

1. Create `include/newdist.h`
2. Inherit from `PDFClass`
3. Implement pure virtual methods: `rand()`, `pdf()`, `cdf()`, `mean()`, `cov()`, `desc()`, `length()`
4. Override optional methods if needed: `var()`, `skew()`, `kurt()`
5. Add test to `tests/test_toolbox.cpp`
6. Rebuild and test

**Template**:
```cpp
#ifndef NEWDIST_H
#define NEWDIST_H

#include "pdfclass.h"
#include <Eigen/Dense>
#include <random>

class NewDist : public PDFClass {
private:
    double param_;
    
public:
    NewDist(double param) : param_(param) {
        state = generate_random_state();
    }
    
    Matrix rand(int n = 1) override { /* ... */ }
    Vector pdf(const Matrix& x) const override { /* ... */ }
    Vector cdf(const Vector& x) const override { /* ... */ }
    Vector mean() const override { /* ... */ }
    Matrix cov() const override { /* ... */ }
    std::string desc() const override { /* ... */ }
    int length() const override { return 1; }
};

#endif
```

### Pattern 2: Adding Methods to Existing Classes

1. Add declaration to header file
2. If header-only: implement inline
3. If needs .cpp: add to existing .cpp file
4. Add test case
5. Rebuild and verify

### Pattern 3: Adding Utility Functions

1. Add template function to `include/utils_sigsys.h`
2. Add test to `tests/test_toolbox.cpp` in UtilsTest suite
3. Rebuild and test

## Common Issues and Solutions

### Issue 1: Template Compilation Errors
**Symptom**: "undefined reference" for template functions
**Solution**: Templates must be in header files, not .cpp files

### Issue 2: Eigen Matrix Dimension Mismatches
**Symptom**: Runtime assertion failures
**Solution**: Always check dimensions with `.rows()`, `.cols()` before operations

### Issue 3: PDFClass Interface Mismatches
**Symptom**: "cannot instantiate abstract class"
**Solution**: Ensure all pure virtual methods are implemented with correct signatures:
- `Matrix rand(int n = 1)` not `Vector rand()`
- `Vector pdf(const Matrix& x) const` not `double pdf(const Vector& x)`
- `double var()` not `Vector var()` for scalar returns

### Issue 4: CMake Not Finding New Test
**Symptom**: Test not listed in `ctest`
**Solution**: Add to `CMakeLists.txt`:
```cmake
add_executable(test_name tests/test_name.cpp)
target_link_libraries(test_name PRIVATE humming_core GTest::GTest GTest::Main pthread)
add_test(NAME TestName COMMAND test_name)
```

## File Organization

### Header Files (`include/`)
- `pdfclass.h` - Base class for distributions
- `*dist.h` - Individual distribution classes (header-only)
- `nl.h` - Nonlinear systems
- `sig.h` - Signal data class
- `sensormod.h` - Sensor models
- `utils_sigsys.h` - Utility functions (header-only templates)
- Existing: `microphone.h`, `localizer.h`, `tdoa_calculator.h`, `gcc_phat.h`, `estimators.h`, `utils.h`, `ndist.h`

### Source Files (`src/`)
- `nl.cpp` - NL class implementations
- Existing: `microphone.cpp`, `localizer.cpp`, `tdoa_calculator.cpp`, `gcc_phat.cpp`, `estimators.cpp`, `utils.cpp`

### Test Files (`tests/`)
- `test_toolbox.cpp` - Distribution and utility tests (16 tests)
- `test_nl.cpp` - NL filtering tests (5 tests)
- `test_sig.cpp` - Signal processing tests (5 tests)
- Existing: `test_localizer.cpp`, `test_TDOA_and_ML.cpp` (2 tests)

## MATLAB Reference Location

All MATLAB source code is in `matlab_resources/`:
- `classes/*.m` - Class definitions
- `mfiles/*.m` - Utility functions
- Use these as reference for implementation details

## Building and Testing

### Full Rebuild
```bash
cd /home/william/projects/hum-finder/build
rm -rf *
cmake ..
make -j4
ctest --output-on-failure
```

### Incremental Build
```bash
cd /home/william/projects/hum-finder/build
make -j4
```

### Run Specific Test
```bash
cd /home/william/projects/hum-finder/build
./test_toolbox    # Distribution and utility tests
./test_nl         # Filtering tests
./test_sig        # Signal processing tests
./run_tests       # Original localizer tests
```

### Run All Tests
```bash
cd /home/william/projects/hum-finder/build
ctest --output-on-failure
```

## Documentation Files

- `TOOLBOX_STATUS.md` - Comprehensive status tracking
- `IMPLEMENTATION_REPORT.md` - Summary of completed work
- `CONTINUATION_GUIDE.md` - This file
- `README.md` - Project overview with plotting commands

## Key Design Decisions

1. **Header-only distributions**: Easier to use, no .cpp needed
2. **PDFClass base**: All distributions inherit common interface
3. **Eigen for linear algebra**: Matrix/Vector operations
4. **Smart pointers for distributions**: `std::shared_ptr<PDFClass>`
5. **MATLAB compatibility**: Function names, behavior match MATLAB
6. **Template utilities**: Generic numerical algorithms in headers

## Current Test Coverage

- LocalizerTest: 2/2 tests ✅
- ToolboxTest: 16/16 tests ✅
- NLTest: 5/5 tests ✅
- SigTest: 5/5 tests ✅
- **Total: 28/28 passing (100%)**

## Suggested Work Session

**Session 1 (1-2 hours)**: Complete distributions + Sig FFT
1. Implement `NCChi2Dist` (~30 min)
2. Add FFT methods to Sig class (~1 hour)
3. Test and document

**Session 2 (2-3 hours)**: LTI framework
1. Create `ltimod.h` base class
2. Implement `lss.h` state-space
3. Add `simulate()` and `step()` methods
4. Create test suite

**Session 3 (1-2 hours)**: Enhance SensorMod
1. Add estimation methods (ls, wls, ml)
2. Implement CRLB computation
3. Test with existing sensor models

## Notes for Tomorrow

- All current code compiles with no warnings
- Test suite is comprehensive and well-structured
- Architecture is solid and extensible
- MATLAB sources available for reference
- Build system (CMake) is properly configured
- Consider adding more signal processing utilities next
- LTI classes would significantly expand capabilities
