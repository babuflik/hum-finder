# Real-Time Estimation Guide

This document explains how to implement real-time localization for the humming sound source tracker.

## Overview

Your system already supports **real-time estimation** through the `Localizer` class, which implements:
- **Extended Kalman Filter (EKF)** - Linearizes the nonlinear TDOA measurement model
- **Unscented Kalman Filter (UKF)** - Uses sigma points for better nonlinear handling
- **State persistence** - Maintains position/velocity state (`x_`) and uncertainty (`P_`) between updates

**Note:** Real-time mode requires **4 microphones** for continuous tracking. The 3-microphone multi-recording approach (`localize_multi_recording`) is designed for offline/batch processing only, as it requires repositioning mics between recordings.

## Microphone Requirements

### Real-Time Mode: 4 Microphones (Standard)

**Required for continuous tracking:**
- 4 synchronized microphones
- Fixed positions (no movement during tracking)
- Can track moving sources
- Updates at ~43 Hz (every 23ms)

### Batch Mode: 3 Microphones (Alternative)

**For stationary sources only:**
- Use `localize_multi_recording` executable instead
- Take multiple recordings at different positions
- Combine offline for improved accuracy
- See `MULTI_RECORDING_GUIDE.md` for details

**Not suitable for real-time** because:
- Requires repositioning mics between measurements
- No temporal tracking (batch processing only)
- Cannot track moving sources

## Current Architecture

```
Audio Capture → TDOA Calculation → Kalman Filter → Position Estimate
    (4 mics)    (GCC-PHAT)          (EKF/UKF)        (x, y, z)
```

### Key Components

1. **TdoaCalculator** (`src/tdoa_calculator.cpp`)
   - Processes 4 audio buffers (1024 samples each)
   - Uses GCC-PHAT algorithm for robust TDOA estimation
   - Outputs 6 TDOA measurements (one per microphone pair)

2. **Localizer** (`src/localizer.cpp`)
   - Takes TDOA measurements
   - Runs Kalman filter update
   - Maintains 6D state: position (x,y,z) + velocity (vx,vy,vz)
   - Returns updated position estimate

## Real-Time Implementation

### Option 1: Continuous Loop (Simplest)

```cpp
#include "localizer.h"
#include "microphone.h"

int main() {
    // Setup
    std::array<Microphone, 4> mics = { /* ... */ };
    Localizer localizer(mics, 0.01, 343.0, FilterType::EKF1);
    
    // Real-time loop
    while (true) {
        // 1. Capture audio from all mics
        std::array<std::vector<double>, 4> audioBuffers;
        captureAudio(audioBuffers);  // Your audio capture function
        
        // 2. Update estimate (Kalman filter update)
        auto position = localizer.locateSource(audioBuffers);
        
        // 3. Use the estimate
        std::cout << "Position: " << position[0] << ", " 
                  << position[1] << ", " << position[2] << std::endl;
        
        // 4. Sleep until next buffer is ready
        std::this_thread::sleep_for(std::chrono::milliseconds(23));
    }
}
```

### Option 2: Threaded with PortAudio (Production)

The provided `src/realtime_localizer.cpp` demonstrates:
- **Thread-safe** position updates
- **Callback-based** audio capture (placeholder for PortAudio)
- **Continuous estimation** at audio buffer rate (~43 Hz for 1024 samples @ 44.1kHz)

**To use with real audio:**

1. Install PortAudio:
   ```bash
   sudo apt-get install portaudio19-dev  # Ubuntu/Debian
   brew install portaudio                 # macOS
   ```

2. Replace `captureAudioBuffers()` with PortAudio callback:
   ```cpp
   #include <portaudio.h>
   
   static int audioCallback(
       const void* inputBuffer,
       void* outputBuffer,
       unsigned long framesPerBuffer,
       const PaStreamCallbackTimeInfo* timeInfo,
       PaStreamCallbackFlags statusFlags,
       void* userData)
   {
       auto* localizer = (RealtimeLocalizer*)userData;
       float* input = (float*)inputBuffer;
       
       // De-interleave 4-channel input into separate buffers
       std::array<std::vector<double>, 4> buffers;
       for (int i = 0; i < 4; i++) {
           buffers[i].resize(framesPerBuffer);
           for (size_t j = 0; j < framesPerBuffer; j++) {
               buffers[i][j] = input[j * 4 + i];
           }
       }
       
       // Process
       localizer->processBuffers(buffers);
       return paContinue;
   }
   ```

### Option 3: Embedded/Microcontroller

For ESP32 or similar:
```cpp
// In FreeRTOS task
void localizationTask(void* params) {
    Localizer localizer(mics, 0.01, 343.0, FilterType::EKF1);
    
    while (1) {
        // Read from I2S microphone array
        std::array<std::vector<double>, 4> buffers;
        i2s_read(buffers);  // Your I2S capture
        
        auto pos = localizer.locateSource(buffers);
        
        // Send via WiFi, Serial, etc.
        sendPosition(pos);
        
        vTaskDelay(pdMS_TO_TICKS(23));
    }
}
```

## Performance Considerations

### Update Rate
- **Audio buffer size**: 1024 samples @ 44.1 kHz = ~23ms per update
- **Maximum update rate**: ~43 Hz
- **Latency**: Single buffer delay (~23ms) + processing time

### Computational Cost per Update
1. **TDOA calculation** (GCC-PHAT):
   - 6 FFT operations (1024-point)
   - 6 IFFT operations
   - Peak finding
   - **Total**: ~5-10ms on modern CPU

2. **Kalman filter**:
   - Matrix operations (6×6 state)
   - Jacobian computation (numerical or analytical)
   - **Total**: <1ms

**Total processing**: ~6-11ms per update (well within 23ms budget)

### Optimization Tips

1. **Use FFTW wisdomodel**:
   ```cpp
   // One-time setup
   fftw_import_wisdom_from_filename("fft_wisdom.dat");
   ```

2. **Pre-allocate buffers**:
   ```cpp
   // Member variables in class
   std::array<std::vector<double>, 4> reusable_buffers_;
   ```

3. **Analytical Jacobian** (instead of numerical):
   ```cpp
   // Current: numerical differentiation in calculateJacobianH()
   // Better: Derive analytical ∂h/∂x for TDOA measurement model
   ```

## State Estimation Details

### State Vector (6D)
```
x = [x, y, z, vx, vy, vz]ᵀ
```

### Process Model (Constant Velocity)
```
x(k+1) = F·x(k) + w,  w ~ N(0, Q)

F = [I₃  dt·I₃]
    [0₃   I₃  ]
```

### Measurement Model (TDOA)
```
z = h(x) + v,  v ~ N(0, R)

h(x) = [
  (‖x - m₁‖ - ‖x - m₂‖) / c,
  (‖x - m₁‖ - ‖x - m₃‖) / c,
  (‖x - m₁‖ - ‖x - m₄‖) / c,
  (‖x - m₂‖ - ‖x - m₃‖) / c,
  (‖x - m₂‖ - ‖x - m₄‖) / c,
  (‖x - m₃‖ - ‖x - m₄‖) / c
]
```

Where:
- `mᵢ` = microphone i position
- `c` = speed of sound (343 m/s)
- `‖·‖` = Euclidean distance

## Tuning Parameters

### Process Noise (Q)
Controls how much the state can change between updates:
```cpp
Q_.setIdentity(); 
Q_.block<3,3>(0,0) *= 0.01;  // Position process noise
Q_.block<3,3>(3,3) *= 0.05;  // Velocity process noise
```
- **Increase** if target moves erratically
- **Decrease** if target moves smoothly

### Measurement Noise (R)
Reflects TDOA measurement uncertainty:
```cpp
R_.setIdentity(); 
R_ *= 1e-6;  // Very confident in TDOA measurements
```
- **Increase** if TDOA is noisy (reverberant environment)
- **Decrease** if TDOA is accurate (anechoic conditions)

### Initial Covariance (P)
Starting uncertainty:
```cpp
P_.setIdentity(); 
P_.block<3,3>(0,0) *= 1.0;    // 1m² position uncertainty
P_.block<3,3>(3,3) *= 0.01;   // 0.01 m²/s² velocity uncertainty
```

## Testing Real-Time Performance

### Build and Run
```bash
cd build
cmake ..
make realtime_localizer
./realtime_localizer
```

### Measure Latency
```cpp
auto start = std::chrono::high_resolution_clock::now();
auto pos = localizer.locateSource(buffers);
auto end = std::chrono::high_resolution_clock::now();
auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
std::cout << "Processing time: " << duration.count() << " μs" << std::endl;
```

### Expected Output
```
Starting real-time localization...
[REALTIME] Position: (0.142, 0.089, 0.023)
[REALTIME] Position: (0.145, 0.091, 0.024)
[REALTIME] Position: (0.147, 0.092, 0.025)
...
```

## Integration with Hardware

### USB Microphone Array
```cpp
// Use PortAudio to open 4-channel USB device
PaStreamParameters inputParams;
inputParams.device = Pa_GetDefaultInputDevice();
inputParams.channelCount = 4;
inputParams.sampleFormat = paFloat32;

Pa_OpenStream(&stream, &inputParams, NULL, 
              44100, 1024, paClipOff, audioCallback, this);
```

### I2S Microphones (ESP32)
```cpp
i2s_config_t i2s_config = {
    .mode = I2S_MODE_MASTER | I2S_MODE_RX,
    .sample_rate = 44100,
    .bits_per_sample = I2S_BITS_PER_SAMPLE_32BIT,
    .channel_format = I2S_CHANNEL_FMT_MULTIPLE,
    .communication_format = I2S_COMM_FORMAT_I2S,
    // ... configure 4 TDM channels
};
```

## Visualization

### Real-Time Plot (Python + ZMQ)
Send position estimates to Python for live plotting:

```cpp
// In C++ (sender)
#include <zmq.hpp>
zmq::context_t ctx;
zmq::socket_t socket(ctx, ZMQ_PUB);
socket.bind("tcp://*:5555");

// Send position
std::stringstream ss;
ss << pos[0] << "," << pos[1] << "," << pos[2];
socket.send(zmq::buffer(ss.str()));
```

```python
# In Python (receiver)
import zmq
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

context = zmq.Context()
socket = context.socket(zmq.SUB)
socket.connect("tcp://localhost:5555")
socket.subscribe("")

fig, ax = plt.subplots()
positions = []

def update(frame):
    msg = socket.recv_string(zmq.NOBLOCK)
    x, y, z = map(float, msg.split(','))
    positions.append([x, y])
    ax.clear()
    ax.plot(*zip(*positions))
    
ani = FuncAnimation(fig, update, interval=50)
plt.show()
```

## Troubleshooting

### Filter Diverges
- **Symptom**: Position estimates grow unbounded
- **Fix**: Increase measurement noise `R` or decrease process noise `Q`

### Estimates Too Smooth
- **Symptom**: Filter doesn't track fast movements
- **Fix**: Increase process noise `Q`

### High Latency
- **Symptom**: Position lags behind actual source
- **Fix**: 
  - Reduce buffer size (trade-off: less frequency resolution)
  - Use EKF instead of UKF (faster)
  - Optimize FFT (use FFTW with wisdom)

### Noisy Estimates
- **Symptom**: Position jumps around
- **Fix**: Decrease measurement noise `R` or increase state covariance `P`

## Next Steps

1. **Integrate real audio capture** (PortAudio/ALSA)
2. **Tune filter parameters** for your environment
3. **Add outlier rejection** (Mahalanobis distance gating)
4. **Implement analytical Jacobian** for speed
5. **Add multiple target tracking** (Multi-Hypothesis Tracking)

## References

- GCC-PHAT: Knapp & Carter (1976) "The generalized correlation method for estimation of time delay"
- EKF: Julier & Uhlmann (2004) "Unscented Filtering and Nonlinear Estimation"
- TDOA Localization: Huang et al. (2001) "Passive Acoustic Source Localization"
