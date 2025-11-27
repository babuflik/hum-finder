# Hum-Finder Documentation

Quick navigation:
- [User Guide](#user-guide) - How to find annoying drone sounds
- [Developer Guide](#developer-guide) - Building and testing
- [API Reference](#api-reference) - Code documentation

---

## User Guide

### Finding an Annoying Drone Sound

**You need:**
- 4 USB microphones (or smartphones)
- Measuring tape
- The annoying sound happening continuously

**Quick steps:**
1. Place 4 mics around your room at different heights
2. Measure each mic's position (x, y, z) from a corner
3. Record 30 seconds of audio when sound is happening
4. Run: `./localize_from_files mic1.wav mic2.wav mic3.wav mic4.wav`

**Expected accuracy:**
- Pure tone: 0.1-0.3 m
- Motor/harmonics: 0.2-0.5 m  
- Broadband noise: 0.5-1.0 m

**FAQ:**
- **Q:** Does it need to be a single frequency?  
  **A:** No! Works for any continuous sound (HVAC, water flow, transformers, etc.)

- **Q:** What if the source is outside my apartment?  
  **A:** The coordinates will show direction. If z > ceiling_height, it's upstairs.

- **Q:** How many microphones minimum?  
  **A:** 4 for 3D localization. More = better accuracy.

See `PRACTICAL_DRONE_FINDING_GUIDE.md` for detailed instructions.

---

## Developer Guide

### Building

```bash
mkdir build && cd build
cmake ..
make
```

**Requirements:**
- CMake 3.14+
- C++17 compiler
- Eigen3
- FFTW3
- Google Test (for tests)

### Testing

```bash
# Run all tests
ctest --output-on-failure

# Run specific test suite
./test_all_estimators
./test_broadband_drone
./test_real_world_drone
```

**Test coverage:**
- 7 test suites, 60 tests total
- Unit tests (toolbox, NL class, Sig class)
- Integration tests (localizer, estimators)
- Real-world scenarios (apartment drone finding)
- Broadband source validation

### Architecture

**Core components:**
- `Localizer`: Main class for source localization
- `TDOACalculator`: Time-difference-of-arrival from audio
- `GCC_PHAT`: Generalized cross-correlation algorithm
- `Estimators`: LS, WLS, ML algorithms
- `NL`: Extended/Unscented Kalman Filter
- `SensorMod`: Sensor modeling and CRLB

**Processing flow:**
```
Audio files → GCC-PHAT → TDOA → Estimators → 3D position
           or
Realtime mics → Streaming → TDOA → EKF → Tracking
```

### File Input Mode

```cpp
// Load audio files
AudioFileLoader loader;
auto [audio1, sr1] = loader.load("mic1.wav");

// Calculate TDOA
TDOACalculator tdoa_calc(mic_positions, speed_of_sound);
auto tdoas = tdoa_calc.calculate(audio_data);

// Estimate position
Localizer localizer(mic_positions);
auto position = localizer.localize(tdoas);
```

### Realtime Mode

```cpp
// Stream from microphones
RealtimeLocalizer localizer(mic_positions);
localizer.start();

// Get continuous position estimates
while (running) {
    auto estimate = localizer.get_current_estimate();
    // Update display...
}
```

---

## API Reference

### Estimators

**Least Squares (LS)**
```cpp
auto [x_est, P_est] = ls(sensor_model, measurements);
```
- Fast, closed-form solution
- Linearizes problem around initial guess
- Best for sources near mic array

**Weighted Least Squares (WLS)**
```cpp
auto [x_est, P_est] = wls(sensor_model, measurements);
```
- Accounts for measurement uncertainty
- Better than LS for noisy data

**Maximum Likelihood (ML)**
```cpp
auto [x_est, P_est] = ml(sensor_model, measurements);
```
- Grid-based search for global optimum
- Slower but more robust
- Good for difficult geometries

**Extended Kalman Filter (EKF)**
```cpp
NL ekf_model(f_func, h_func, dims, fs);
auto [x_est, P_est] = ekf_model.ekf(measurements);
```
- Best overall accuracy (0.007m RMSE)
- Handles time-series data
- Recommended for most applications

**Unscented Kalman Filter (UKF)**
```cpp
auto [x_est, P_est] = ekf_model.ukf(measurements);
```
- Better for highly nonlinear problems
- Slightly slower than EKF

**Cramér-Rao Lower Bound (CRLB)**
```cpp
auto crlb_cov = crlb(sensor_model, true_position);
```
- Theoretical best-case performance
- Use to validate estimator quality

### Source Types Performance

| Source | Spectrum | Typical Error |
|--------|----------|---------------|
| HVAC | 60Hz + harmonics + broadband | 0.07m |
| Refrigerator | 60Hz + harmonics | 0.01m |
| Transformer | 60Hz + odd harmonics | 0.08m |
| Water leak | Broadband turbulence | 0.12m |
| Speaker buzz | Complex harmonics | 0.5-1.5m |

All validated in `test_broadband_drone.cpp`

---

## Project Status

**Working:**
- ✅ File-based localization
- ✅ Realtime streaming (basic)
- ✅ All 6 estimators (LS, WLS, ML, EKF, UKF, CRLB)
- ✅ 3D positioning
- ✅ Broadband source support
- ✅ Comprehensive test suite

**To Do:**
- Audio file I/O integration with estimators
- Command-line interface improvements
- Visualization tools
- Particle filter (PF) implementation

**Recent additions:**
- Broadband/multi-frequency source validation (Nov 2025)
- Real-world apartment scenarios (Nov 2025)
- All estimators comparison suite (Nov 2025)
