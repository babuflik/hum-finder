# File Input and Real-Time Microphone Usage Guide

This guide shows how to use the localization system with both file inputs and real-time microphones.

## Quick Start

### Build Everything
```bash
cd build
cmake ..
make -j4
```

This creates three executables:
- `localize_from_files` - Process pre-recorded audio files
- `localize_realtime` - Real-time processing with live microphones
- `streaming_example` - Demo with simulated moving source

---

## Option 1: File Input (WAV or CSV)

### Generate Test Files

**CSV files:**
```bash
python3 tools/generate_test_audio.py --source-x 0.15 --source-y 0.10 --output-dir test_data
```

**WAV files:**
```bash
python3 tools/generate_test_wav.py --source-x 0.15 --source-y 0.10 --output-dir test_data
```

This creates:
- `test_data/mic1.csv` (or `.wav`)
- `test_data/mic2.csv` (or `.wav`)
- `test_data/mic3.csv` (or `.wav`)
- `test_data/mic4.csv` (or `.wav`)

### Run Localization from Files

```bash
./localize_from_files test_data/mic1.csv test_data/mic2.csv test_data/mic3.csv test_data/mic4.csv
```

Or with WAV files:
```bash
./localize_from_files test_data/mic1.wav test_data/mic2.wav test_data/mic3.wav test_data/mic4.wav
```

### Expected Output
```
=== Sound Source Localization from Files ===

Microphone positions:
  M1: (0, 0, 0)
  M2: (0.2, 0, 0)
  M3: (0.2, 0.2, 0)
  M4: (0, 0.2, 0)

Loading audio files...
Loading CSV: test_data/mic1.csv
  Loaded 1024 samples
Loading CSV: test_data/mic2.csv
  Loaded 1024 samples
...

Running localization...

=== Results ===
Estimated position: (0.148, 0.102, 0.051) meters
Position uncertainty: 0.0234 m (RMS)

Localization complete!
```

### File Formats Supported

#### CSV Format
One sample per line, plain text:
```
0.12345678
-0.23456789
0.34567890
...
```

#### WAV Format
Standard PCM WAV files:
- **Supported bit depths**: 8, 16, or 32 bit
- **Channels**: Mono or stereo (converted to mono)
- **Sample rate**: Any (recommended: 44.1 kHz)
- **Format**: PCM uncompressed

---

## Option 2: Real-Time Microphones

### Using Simulated Audio (Current)

```bash
./localize_realtime
```

This runs with simulated microphone data. Output:
```
=== Real-Time Sound Source Localization ===

NOTE: Using simulated audio. See REALTIME_GUIDE.md for real microphone setup.

Microphone positions:
  M1: (0, 0, 0)
  M2: (0.2, 0, 0)
  M3: (0.2, 0.2, 0)
  M4: (0, 0.2, 0)

Starting real-time localization...
Update rate: ~43 Hz
Buffer size: 1024 samples

Press Ctrl+C to stop.

[  10] Position: ( 0.142,  0.089,  0.023) m  |  Uncertainty: 0.0456 m  |  Process:  8234 μs
[  20] Position: ( 0.145,  0.091,  0.024) m  |  Uncertainty: 0.0312 m  |  Process:  7891 μs
...
```

### Using Real Microphones with PortAudio

#### 1. Install PortAudio
```bash
# Ubuntu/Debian
sudo apt-get install portaudio19-dev

# macOS
brew install portaudio

# Fedora
sudo dnf install portaudio-devel
```

#### 2. Add PortAudio to CMake

Edit `CMakeLists.txt`:
```cmake
find_package(PkgConfig REQUIRED)
pkg_check_modules(PORTAUDIO REQUIRED portaudio-2.0)

target_link_libraries(localize_realtime PRIVATE
    humming_core
    pthread
    ${PORTAUDIO_LIBRARIES}
)

target_include_directories(localize_realtime PRIVATE
    ${PORTAUDIO_INCLUDE_DIRS}
)
```

#### 3. Replace Audio Capture in `src/localize_realtime.cpp`

Replace the `captureAudioFromMicrophones()` function with:

```cpp
#include <portaudio.h>

// Global buffer for callback
std::array<std::vector<double>, MIC_COUNT> g_audio_buffers;
std::mutex g_buffer_mutex;
bool g_buffer_ready = false;

int audioCallback(const void* inputBuffer, void* outputBuffer,
                  unsigned long framesPerBuffer,
                  const PaStreamCallbackTimeInfo* timeInfo,
                  PaStreamCallbackFlags statusFlags,
                  void* userData)
{
    const float* input = (const float*)inputBuffer;
    
    std::lock_guard<std::mutex> lock(g_buffer_mutex);
    
    // De-interleave 4-channel input
    for (size_t i = 0; i < framesPerBuffer; ++i) {
        for (int ch = 0; ch < MIC_COUNT; ++ch) {
            g_audio_buffers[ch][i] = input[i * MIC_COUNT + ch];
        }
    }
    
    g_buffer_ready = true;
    return paContinue;
}

void initPortAudio() {
    Pa_Initialize();
    
    PaStreamParameters inputParams;
    inputParams.device = Pa_GetDefaultInputDevice();
    if (inputParams.device == paNoDevice) {
        throw std::runtime_error("No default input device found");
    }
    
    inputParams.channelCount = MIC_COUNT;
    inputParams.sampleFormat = paFloat32;
    inputParams.suggestedLatency = Pa_GetDeviceInfo(inputParams.device)->defaultLowInputLatency;
    inputParams.hostApiSpecificStreamInfo = NULL;
    
    PaStream* stream;
    PaError err = Pa_OpenStream(
        &stream,
        &inputParams,
        NULL,  // No output
        44100,
        1024,
        paClipOff,
        audioCallback,
        NULL
    );
    
    if (err != paNoError) {
        throw std::runtime_error(std::string("PortAudio error: ") + Pa_GetErrorText(err));
    }
    
    Pa_StartStream(stream);
    std::cout << "PortAudio stream started successfully" << std::endl;
}

void captureAudioFromMicrophones(std::array<std::vector<double>, MIC_COUNT>& buffers) {
    std::lock_guard<std::mutex> lock(g_buffer_mutex);
    while (!g_buffer_ready) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    buffers = g_audio_buffers;
    g_buffer_ready = false;
}
```

Then in `main()`, add before the loop:
```cpp
// Initialize buffers
for (int i = 0; i < MIC_COUNT; ++i) {
    g_audio_buffers[i].resize(Microphone::BUFFER_SIZE);
}

// Start audio capture
try {
    initPortAudio();
} catch (const std::exception& e) {
    std::cerr << "Failed to initialize audio: " << e.what() << std::endl;
    return 1;
}
```

#### 4. Rebuild and Run
```bash
cd build
cmake ..
make localize_realtime
./localize_realtime
```

---

## Hardware Setup

### Recommended Microphone Array

**USB Audio Interface** (4+ channels):
- Focusrite Scarlett 4i4
- Behringer U-PHORIA UMC404HD
- MOTU M4
- Any multi-channel USB audio interface

**Microphones**:
- 4x omnidirectional microphones
- Recommendation: Small diaphragm condensers or measurement mics
- Examples: Behringer ECM8000, Dayton Audio EMM-6

### Physical Setup

```
    M3 (0.2, 0.2, 0)  ●────────● M2 (0.2, 0, 0)
                      │        │
                      │  20cm  │
                      │        │
    M4 (0, 0.2, 0)    ●────────● M1 (0, 0, 0)
```

**Tips:**
- Place microphones in a square formation (20cm spacing)
- Keep array on a flat surface (all z=0)
- Ensure microphones are securely mounted
- Minimize reflective surfaces nearby
- Place sound source within ~1-2 meters

---

## Command-Line Options

### localize_from_files
```bash
./localize_from_files <mic1_file> <mic2_file> <mic3_file> <mic4_file>

Arguments:
  mic1_file  Audio file for microphone 1 (.wav or .csv)
  mic2_file  Audio file for microphone 2 (.wav or .csv)
  mic3_file  Audio file for microphone 3 (.wav or .csv)
  mic4_file  Audio file for microphone 4 (.wav or .csv)
```

### localize_realtime
```bash
./localize_realtime [options]

Options:
  --duration <seconds>  Run for specified duration (default: infinite)
  --help               Show help message
```

---

## Recording Your Own Audio Files

### Using Audacity
1. File → Import → Audio...
2. Select all 4 channels
3. File → Export → Export Multiple...
4. Format: WAV (16-bit PCM)
5. Split by: Tracks

### Using arecord (Linux)
```bash
# Record 4-channel audio
arecord -D hw:1,0 -f S16_LE -c 4 -r 44100 -d 5 recording.wav

# Split into separate files
ffmpeg -i recording.wav -map_channel 0.0.0 mic1.wav \
                        -map_channel 0.0.1 mic2.wav \
                        -map_channel 0.0.2 mic3.wav \
                        -map_channel 0.0.3 mic4.wav
```

### Using SoX
```bash
# Record 4 channels for 5 seconds
sox -t alsa hw:1,0 -c 4 recording.wav trim 0 5

# Split channels
sox recording.wav mic1.wav remix 1
sox recording.wav mic2.wav remix 2
sox recording.wav mic3.wav remix 3
sox recording.wav mic4.wav remix 4
```

---

## Troubleshooting

### File Loading Issues

**Error: "Cannot open file"**
- Check file path is correct
- Verify file permissions

**Error: "Invalid WAV format"**
- Ensure file is uncompressed PCM WAV
- Try re-exporting from audio editor

**Error: "No samples loaded"**
- Check file isn't empty
- Verify CSV has one number per line

### Real-Time Issues

**Error: "No default input device"**
- Check microphone is connected
- Run `arecord -l` (Linux) or check System Preferences (macOS)

**Noisy/Unstable Estimates**
- Check microphone connections
- Increase `R` (measurement noise) in `localizer.cpp`
- Verify microphone positions are accurate

**Filter Diverges**
- Decrease `Q` (process noise)
- Check for outliers in TDOA measurements

---

## Next Steps

1. **Test with sample files**: Generate CSV/WAV and verify localization works
2. **Connect real microphones**: Set up 4-channel USB audio interface
3. **Integrate PortAudio**: Follow instructions above
4. **Calibrate**: Measure actual microphone positions
5. **Tune parameters**: Adjust Q and R for your environment

For more details, see:
- `REALTIME_GUIDE.md` - Comprehensive real-time guide
- `ARCHITECTURE_REALTIME.md` - System architecture diagrams
- `QUICKREF_REALTIME.md` - Quick reference
