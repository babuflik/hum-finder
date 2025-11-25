# Real-Time Estimation - Quick Reference

## What You Have

Your system **already supports real-time estimation**! The `Localizer` class maintains state between calls.

## Basic Usage

```cpp
#include "localizer.h"
#include "microphone.h"

// 1. Setup (once)
std::array<Microphone, 4> mics = { /* positions */ };
Localizer localizer(mics, 0.01, 343.0, FilterType::EKF1);

// 2. Real-time loop
while (running) {
    // Capture 1024 samples from each mic
    std::array<std::vector<double>, 4> audioBuffers;
    captureAudio(audioBuffers);
    
    // Update estimate (Kalman filter)
    auto position = localizer.locateSource(audioBuffers);
    
    // Use position
    std::cout << position[0] << ", " << position[1] << std::endl;
}
```

## Key Points

1. **State is maintained** - Each call to `locateSource()` updates the filter state
2. **No initialization needed** - Filter starts with initial guess and converges
3. **Update rate** - ~43 Hz (limited by audio buffer size: 1024 samples @ 44.1 kHz)

## Performance

- **Processing time**: ~6-11 ms per update
- **Buffer time**: ~23 ms
- **Total latency**: ~23-34 ms (one buffer + processing)

## Examples

### Run the streaming demo
```bash
cd build
cmake ..
make streaming_example
./streaming_example
```

### Run real-time localizer
```bash
make realtime_localizer
./realtime_localizer
```

## Tuning

Adjust in `localizer.cpp` constructor:

```cpp
// Position uncertainty (higher = less smooth, more responsive)
Q_.block<3,3>(0,0) *= 0.01;

// Velocity uncertainty
Q_.block<3,3>(3,3) *= 0.05;

// Measurement noise (higher = trust TDOA less)
R_ *= 1e-6;
```

## Integration with Hardware

### PortAudio (USB Microphones)
```cpp
#include <portaudio.h>

int callback(const void* in, void* out, unsigned long frames, 
             const PaStreamCallbackTimeInfo* info,
             PaStreamCallbackFlags flags, void* userData) {
    auto* loc = (Localizer*)userData;
    // Process 4-channel audio...
    return paContinue;
}

Pa_OpenStream(&stream, &inputParams, NULL, 44100, 1024, 
              paClipOff, callback, &localizer);
Pa_StartStream(stream);
```

### ESP32 (I2S Microphones)
```cpp
void task(void* params) {
    while (1) {
        std::array<std::vector<double>, 4> buffers;
        i2s_read(I2S_NUM_0, buffers, 1024, portMAX_DELAY);
        auto pos = localizer.locateSource(buffers);
        vTaskDelay(pdMS_TO_TICKS(23));
    }
}
```

## See Also

- `REALTIME_GUIDE.md` - Comprehensive guide
- `examples/streaming_example.cpp` - Working example with simulated moving source
- `src/realtime_localizer.cpp` - Thread-based implementation
