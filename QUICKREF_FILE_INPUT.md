# Quick Reference: File Input & Real-Time Usage

## File Input (WAV or CSV)

### One-Liner Test
```bash
python3 tools/generate_test_audio.py --output-dir . && \
./localize_from_files mic1.csv mic2.csv mic3.csv mic4.csv
```

### Usage
```bash
./localize_from_files <mic1> <mic2> <mic3> <mic4>
```

### Generate Test Files
```bash
# CSV files
python3 tools/generate_test_audio.py --source-x 0.15 --source-y 0.10 --output-dir data

# WAV files  
python3 tools/generate_test_wav.py --source-x 0.15 --source-y 0.10 --output-dir data
```

## Real-Time Microphones

### Usage
```bash
./localize_realtime [--duration <seconds>]
```

### Examples
```bash
./localize_realtime                 # Run until Ctrl+C
./localize_realtime --duration 10   # Run for 10 seconds
```

## File Formats

### CSV Format
```
0.12345678
-0.23456789
0.34567890
...
```

### WAV Format
- PCM uncompressed (8/16/32-bit)
- Mono or stereo
- Any sample rate (44.1 kHz recommended)

## Build Commands

```bash
cd build
cmake ..
make localize_from_files localize_realtime -j4
```

## Microphone Array Setup

```
M3 (0.2, 0.2, 0) ●────● M2 (0.2, 0, 0)
                 │    │
                 │ 20 │ cm
                 │    │
M4 (0, 0.2, 0)   ●────● M1 (0, 0, 0)
```

## PortAudio Integration

See `FILE_INPUT_GUIDE.md` section "Using Real Microphones with PortAudio"

**Quick steps:**
1. Install: `sudo apt-get install portaudio19-dev`
2. Update CMakeLists.txt
3. Replace `captureAudioFromMicrophones()` function
4. Rebuild

## Tools

| Tool | Purpose |
|------|---------|
| `generate_test_audio.py` | Generate CSV test files |
| `generate_test_wav.py` | Generate WAV test files |
| `localize_from_files` | Process audio files |
| `localize_realtime` | Real-time microphones |

## Documentation

- `FILE_INPUT_GUIDE.md` - **Complete guide**
- `REALTIME_GUIDE.md` - Real-time details
- `FILE_INPUT_IMPLEMENTATION.md` - What was added

## Help

```bash
./localize_from_files --help
./localize_realtime --help
```
