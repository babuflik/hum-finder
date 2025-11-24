# MATLAB SigSys Toolbox C++ Implementation

This document tracks the implementation status of the MATLAB Signals and Systems toolbox in C++.

## Implementation Status

### ✅ Core Infrastructure (DONE)

1. **PDFClass** (`include/pdfclass.h`)
   - Base class for all probability distributions
   - Pure virtual methods: `rand()`, `pdf()`, `cdf()`, `mean()`, `cov()`, `desc()`, `length()`
   - Optional methods: `median()`, `mode()`, `var()`, `std()`, `skew()`, `kurt()`, `symbolic()`

2. **NDist** (`include/ndist.h`)
   - Gaussian/Normal distribution N(mu, P)
   - ✅ Implemented: `rand()`, `pdf()`, `cdf()`, `mean()`, `cov()`, `symbolic()`
   - ✅ Setters: `set_mu()`, `set_P()`
   - ✅ Indexing: `subsref()` for MATLAB-style X(i) indexing
   - ✅ Validation: `is_valid_covariance()` helper

3. **Sig** (`include/sig.h`)
   - Signal data class with time series y, input u, state x
   - ✅ Properties: `y`, `t`, `u`, `x`, `fs`, `Py`, `Px`, `yMC`, `xMC`, `MC`
   - ✅ Labels: `ylabel`, `ulabel`, `xlabel`, `tlabel`
   - ✅ Metadata: `nn`, `name`, `desc`, `marker`, `markerlabel`
   - ✅ Constructors: empty, (y,fs), (y,t), (y,t,u), (y,t,u,x)
   - ✅ Methods: `length()`, `size(dim)`, automatic fs detection
   - ⏳ TODO: plot(), resample(), filter(), FFT, arithmetic operators

4. **SensorMod** (`include/sensormod.h`)
   - Sensor model class y = h(t,x,u,th) + e
   - ✅ Properties: `h`, `nn`, `x0`, `th`, `pv`, `pe`, `fs`
   - ✅ Methods: `simulate()`, `likelihood_function()`
   - ⏳ TODO: `ls()`, `wls()`, `ml()`, `estimate()`, `crlb()` as methods

### ✅ Distribution Classes (10/11 DONE)

5. **NDist** - Gaussian/Normal distribution ✅ COMPLETE
6. **UDist** (`include/udist.h`) ✅ COMPLETE
   - Uniform distribution U(a, b)
   - Supports scalar and multidimensional
   
7. **ExpDist** (`include/expdist.h`) ✅ COMPLETE
   - Exponential distribution Exp(λ)
   - All moments implemented
   
8. **Chi2Dist** (`include/chi2dist.h`) ✅ COMPLETE
   - Chi-squared distribution χ²(n)
   - Custom incomplete gamma for CDF
   
9. **GammaDist** (`include/gammadist.h`) ✅ COMPLETE
   - Gamma distribution Γ(k, θ)
   - Log-space computations for stability
   
10. **GMDist** (`include/gmdist.h`) ✅ COMPLETE
    - Gaussian mixture distribution
    - K components with weights
    
11. **TDist** (`include/tdist.h`) ✅ COMPLETE
    - Student's t-distribution T(nu)
    - Multivariate support
    
12. **BetaDist** (`include/betadist.h`) ✅ COMPLETE
    - Beta distribution Beta(α, β)
    - Incomplete beta function for CDF
    
13. **EmpDist** (`include/empdist.h`) ✅ COMPLETE
    - Empirical distribution from samples
    - Kernel density estimation
    
14. **LogNDist** (`include/logndist.h`) ✅ COMPLETE
    - Log-normal distribution LogN(μ, σ)
    - All moments closed-form
    
15. **NCChi2Dist** ⏳ TODO
    - Non-central chi-squared distribution

### ✅ Utility Functions (7/50+ DONE)

16. **utils_sigsys.h** ✅ PARTIAL
    - ✅ `numgrad()` - Numerical gradient (central differences)
    - ✅ `numhess()` - Numerical Hessian
    - ✅ `sqrtcov()` - Cholesky decomposition
    - ✅ `iscov()` - Enhanced covariance validation with error codes
    - ✅ `confellipse()` - Confidence ellipse for 2D plots
    - ✅ `getwindow()` - Hamming, Hann, Blackman windows
    - ✅ `condition()` - Matrix condition number
    - ⏳ `filtfilt()` - Zero-phase filtering
    - ⏳ `ncfilter()` - NC filter
    - ⏳ `interp()` - Interpolation
    - ⏳ `resample()` - Resampling

### ✅ NL Class - Nonlinear Systems (COMPLETE)

6. **NL - Nonlinear System Class** (`nl.h`, `nl.cpp`) ✅ COMPLETE
   - Parent class for nonlinear dynamic systems
   - State equation: x(t+1) = f(t, x(t), u(t), th, v(t))
   - Measurement: y(t) = h(t, x(t), u(t), th, e(t))
   - ✅ `simulate()` - Monte Carlo simulation
   - ✅ `ekf()` - Extended Kalman Filter
   - ✅ `ukf()` - Unscented Kalman Filter
   - ✅ `pf()` - Particle Filter with systematic resampling
   - ⏳ `pmf()` - Point Mass Filter (TODO)
   - ⏳ `crlb()` - Cramér-Rao Lower Bound (TODO)

### ❌ Not Yet Implemented

7. **LTI System Classes**
   - ❌ `ltimod.h` - Base class for Linear Time-Invariant systems
   - ❌ `lss.h` - State Space models (continuous/discrete, SISO/MIMO)
   - ❌ `ltf.h` - Transfer Functions
   - Methods: `simulate()`, `estimate()`, `bode()`, `nyquist()`, `step()`, `impulse()`

8. **Signal Processing Classes**
   - ❌ `ft.h` - Fourier Transform object
   - ❌ `spec.h` - Spectrum object (PSD, cross-spectrum)
   - ❌ `covfun.h` - Covariance function object
   - ❌ `tfd.h` - Time-Frequency Description (STFT)

9. **Model Estimation Classes**
   - ❌ `arx.h` - ARX models (Auto-Regressive with eXogenous input)
   - ❌ `ar.h` - AR models (time-series)
   - ❌ `fir.h` - Finite Impulse Response models
   - ❌ `rarx.h` - Recursive ARX (time-varying)
   - ❌ `jarx.h` - Jump ARX (switching)

10. **Utility Functions** (create `utils_sigsys.h/cpp`)
    - ✅ `numgrad()` - Numerical gradient ✅ DONE
    - ✅ `numhess()` - Numerical Hessian ✅ DONE
    - ✅ `confellipse()` - Confidence ellipse computation ✅ DONE
    - ✅ `iscov()` - Enhanced covariance validation ✅ DONE
    - ✅ `sqrtcov()` - Square root of covariance (Cholesky) ✅ DONE
    - ✅ `getwindow()` - Window functions ✅ DONE
    - ✅ `condition()` - Matrix condition number ✅ DONE
    - ❌ `getfilter()` - Standard filter designs (basic done, needs enhancement)
    - ❌ `filtfilt()` - Zero-phase filtering
    - ❌ `ncfilter()` - NC filter
    - ❌ `interp()` - Interpolation
    - ❌ `resample()` - Resampling

## Current Test Status

- ✅ All existing tests pass (2/2 test suites):
  - ✅ LocalizerTest (1 test)
  - ✅ ToolboxTest (16 tests):
    - 10 distribution tests (NDist, UDist, ExpDist, Chi2Dist, GammaDist, GMDist, TDist, BetaDist, EmpDist, LogNDist)
    - 6 utility tests (NumGrad, NumHess, SqrtCov, IsCov, ConfEllipse, GetWindow)
- ✅ Backward compatibility maintained
- ✅ Build system updated for new headers
- ✅ All tests passing (17/17)

## Next Steps

### Priority 1: NL Class for Nonlinear Filtering ⏭️ NEXT
1. Create `nl.h` as parent of `sensormod`
2. Implement EKF, UKF basic versions
3. Add particle filter support

### Priority 2: Complete Distribution Classes
1. Implement `ncchi2dist` (non-central chi-squared)

### Priority 3: Essential Utility Functions
1. Enhanced `getfilter()` for Butterworth, Chebyshev
2. `filtfilt()` for zero-phase filtering
3. `resample()`, `interp()` for signal processing

### Priority 4: Signal Processing Extensions
1. Add FFT/IFFT to Sig class
2. Implement `ft` class for frequency domain
3. Add `spec` for spectral analysis

## Usage Examples

### NDist (Normal Distribution)
```cpp
#include "ndist.h"

// Create 2D Gaussian
Eigen::Vector2d mu(1.0, 2.0);
Eigen::Matrix2d P;
P << 1.0, 0.3,
     0.3, 2.0;
NDist X(mu, P);

// Generate samples
Eigen::MatrixXd samples = X.rand(100);  // 100 samples

// Evaluate PDF
Eigen::MatrixXd test_points(10, 2);
// ... fill test_points ...
Eigen::VectorXd pdf_values = X.pdf(test_points);

// Get statistics
Eigen::Vector2d mean = X.mean();
Eigen::Matrix2d cov = X.cov();
```

### Sig (Signal Class)
```cpp
#include "sig.h"

// Uniformly sampled signal
int N = 100;
double fs = 0.01;
Eigen::MatrixXd y = Eigen::MatrixXd::Random(N, 2);  // 2 outputs
Sig s1(y, fs);

// Non-uniformly sampled
Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(N, 0, 1);
Sig s2(y, t);

// With input and state
Eigen::MatrixXd u = Eigen::MatrixXd::Random(N, 1);
Eigen::MatrixXd x = Eigen::MatrixXd::Random(N, 3);
Sig s3(y, t, u, x);

// Access properties
int length = s3.length();
double sampling_freq = s3.fs;
std::cout << s3.ylabel[0] << std::endl;  // "y1"
```

### SensorMod (Sensor Model)
```cpp
#include "sensormod.h"

// Define measurement function: range measurements
auto h_func = [](double t, const Eigen::VectorXd& x, 
                  const Eigen::VectorXd& u, const Eigen::VectorXd& th) {
    // Compute ranges from sensors to target at x
    Eigen::VectorXd y(3);
    // ... range calculations ...
    return y;
};

Eigen::Vector4i nn(2, 0, 3, 0);  // nx=2, nu=0, ny=3, nth=0
SensorMod sensor(h_func, nn);

// Set parameters
sensor.x0 << 1.0, 0.5;  // Initial target position
sensor.pe = 0.1 * Eigen::Matrix3d::Identity();  // Measurement noise

// Simulate
Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(10, 0, 1);
Sig y = sensor.simulate(t);
```

## Build Instructions

The toolbox is integrated into the existing CMake build:

```bash
cd build
cmake ..
make -j4
ctest
```

All new header files are header-only where possible, minimizing build complexity.

## Architecture Notes

1. **Header-Only Design**: Most distribution classes are header-only for ease of use
2. **Eigen Integration**: All linear algebra uses Eigen3 for performance and compatibility
3. **MATLAB Compatibility**: Method names and signatures match MATLAB where possible
4. **Modern C++17**: Uses C++17 features (std::optional, std::variant where appropriate)
5. **Exception Safety**: Proper error handling with descriptive error messages

## Contributing

When adding new classes:
1. Inherit from appropriate base class (`PDFClass`, `NL`, etc.)
2. Implement all pure virtual methods
3. Add constructor validation
4. Include usage examples in comments
5. Add unit tests in `tests/`
6. Update this README

## References

- MATLAB SigSys Toolbox v2025.2 (Fredrik Gustafsson, Sigmoid AB)
- Original MATLAB code in `matlab_resources/`
