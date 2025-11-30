# Real-World Testing Protocol

Complete step-by-step guide for conducting acoustic source localization in a real apartment/building.

---

## Table of Contents

1. [Pre-Test Planning](#pre-test-planning)
2. [Equipment Setup](#equipment-setup)
3. [Environmental Preparation](#environmental-preparation)
4. [Microphone Positioning](#microphone-positioning)
5. [Recording Setup](#recording-setup)
6. [Data Collection](#data-collection)
7. [Running the Localization](#running-the-localization)
8. [Interpreting Results](#interpreting-results)
9. [Troubleshooting](#troubleshooting)
10. [Advanced: Tuning Parameters](#advanced-tuning-parameters)

---

## Pre-Test Planning

### 1. Identify the Sound Characteristics

**Before setting up equipment:**

- [ ] **Verify the sound is continuous** (not intermittent clicks/pops)
- [ ] **Note when it occurs** (24/7, only at night, specific times)
- [ ] **Estimate frequency range** (low hum ~60Hz, mid buzz ~120-500Hz, etc.)
- [ ] **Check if it's directional** (walk around - does it seem louder in one area?)
- [ ] **Record on phone** (quick test to verify it's audible on recording)

**Sound types that work well:**
- ✅ HVAC hum (continuous motor noise)
- ✅ Refrigerator/freezer (60Hz hum + harmonics)
- ✅ Transformer buzz (electrical equipment)
- ✅ Water pipe leak/flow (steady broadband)
- ✅ Mechanical equipment (pumps, motors)

**Challenging cases:**
- ⚠ Very intermittent sounds (< 10 seconds duration)
- ⚠ Multiple simultaneous sources at same level
- ⚠ Pure high-frequency sounds (> 2000 Hz)
- ⚠ Sources more than 20m away

### 2. Room Reconnaissance

**Measure and sketch your room:**
```
     ← 5.2m →
  ┌──────────┐  ↑
  │          │  │
  │          │  │ 4.1m
  │          │  │
  └──────────┘  ↓
  
  Height: 2.8m
```

**Note:**
- Wall materials (concrete, drywall, etc.)
- Large furniture positions
- Windows and doors
- Known noise sources (refrigerator, AC vents, etc.)

---

## Equipment Setup

### Minimum Equipment

#### Standard 4-Mic Setup (Budget: ~$100-150)

**4 USB Microphones:**
- Recommended: [Blue Snowball](https://www.bluemic.com/) (~$50 each) or similar
- Alternative: 4 identical USB mics (consistency matters more than quality)
- Budget option: 4 smartphones with recording apps

**Computer:**
- Laptop with 4+ USB ports (or use USB hub)
- Linux/macOS/Windows (Linux preferred for this software)

**Accessories:**
- USB hub (powered, if using 4 USB mics)
- Extension cables (to position mics far apart)
- Measuring tape (at least 5m)
- Notebook and pen (to record positions)
- Optional: Mic stands/clips (for precise positioning)
- Optional: Small weights (to stabilize mics)

#### NEW: Budget 3-Mic Setup (Budget: ~$75-120)

**For stationary sources only!**

**3 USB Microphones or Smartphones:**
- Same options as above, but only need 3 devices
- **25% cost savings** compared to 4-mic setup
- Requires taking multiple recordings (2-3 sessions)

**Tradeoffs:**
- ✅ Lower hardware cost
- ✅ Can achieve similar accuracy with multiple recordings
- ⚠️ Source MUST be stationary between recordings
- ⚠️ Takes longer (multiple recording sessions)
- ⚠️ Must reposition microphones between recordings

**See "3-Microphone Multi-Recording" section below for details**

### Recommended Equipment (Budget: ~$300-500)

**4-Channel Audio Interface:**
- Focusrite Scarlett 4i4 (~$250)
- PreSonus AudioBox 1818VSL (~$400)
- Benefit: Better synchronization than USB mics

**4 Matched Microphones:**
- Same model for all 4 mics
- Omnidirectional pattern preferred
- Examples: Behringer ECM8000, Dayton Audio EMM-6

---

## Environmental Preparation

### 1. Minimize Background Noise

**Critical for good results:**

- [ ] **Turn off all unnecessary noise sources**
  - Air conditioning / heating (if possible)
  - Fans, air purifiers
  - TV, radio, computer fans
  - Refrigerator (temporarily unplug if NOT the source)

- [ ] **Close windows and doors** (reduce traffic noise)

- [ ] **Choose quiet time** (late night, early morning)

- [ ] **Warn household members** (no talking, walking during recording)

**Note:** The drone sound you're tracking should be the LOUDEST sound in the recording. Ideal SNR (Signal-to-Noise Ratio) > 10 dB.

### 2. Measure Room Dimensions

**Establish coordinate system:**

```
Origin (0,0,0) = Front-left corner at floor level

     Y-axis (width) →
  ┌──────────────────┐
X │ (0,0,0)          │ ↑ Z-axis
│ │                  │ │ (height)
a │                  │ │
x │                  │ ↓
i │                  │
s └──────────────────┘
↓
```

**Measure:**
- Room length (X-axis): _______ m
- Room width (Y-axis): _______ m
- Room height (Z-axis): _______ m

---

## Microphone Positioning

### Critical Rules for Mic Placement

**Golden Rules:**
1. **Non-planar arrangement** - All 4 mics must NOT be in same plane
2. **Spread out** - Maximum distance between mics (use whole room)
3. **Different heights** - Vary Z-coordinates (floor, table, shelf, wall)
4. **Surround suspected area** - Form a 3D "cage" around where you think the sound is

### Recommended Configuration

**For a 5m × 4m × 2.8m room:**

```
Mic 1: [1.5, 1.0, 0.4]   - Coffee table height
Mic 2: [4.0, 3.5, 1.2]   - Shelf on opposite wall
Mic 3: [4.5, 0.8, 0.5]   - Low stand, different corner
Mic 4: [0.5, 3.0, 2.0]   - High on wall (mounted or tall stand)
```

**Visualization (Top View):**
```
     ← 4m →
  ┌─────────┐
5 │ 4     2 │  Heights:
m │    ?    │  Mic1: 0.4m
  │ 1     3 │  Mic2: 1.2m
  └─────────┘  Mic3: 0.5m
               Mic4: 2.0m
  ? = suspected source location
```

### Step-by-Step Positioning

**For each microphone:**

1. **Place microphone** at desired location
2. **Measure position from origin corner:**
   - X (length): _______ m
   - Y (width): _______ m
   - Z (height from floor): _______ m
3. **Record in notebook:**
   ```
   Mic 1: X=_____, Y=_____, Z=_____
   Mic 2: X=_____, Y=_____, Z=_____
   Mic 3: X=_____, Y=_____, Z=_____
   Mic 4: X=_____, Y=_____, Z=_____
   ```
4. **Mark position** with tape (in case mic moves)

**Pro tip:** Take photos of mic positions for reference.

---

## 3-Microphone Multi-Recording (Alternative Method)

### When to Use This Method

**✅ Use 3-mic multi-recording when:**
- You only have 3 recording devices (phones, USB mics, etc.)
- The sound source is definitely stationary (HVAC, transformer, refrigerator)
- You want to minimize hardware cost
- You have time for 2-3 recording sessions

**❌ Don't use this method when:**
- Source might be moving or intermittent
- You need results immediately
- You already have 4 microphones

### Concept

Since the source doesn't move, you can take multiple recordings with 3 mics at different positions and combine them for better accuracy.

**Standard approach:** 4 mics × 1 recording = 1 estimate  
**Multi-recording approach:** 3 mics × 2 recordings = 2 estimates → averaged result

### Step-by-Step Procedure

#### Recording 1: Ground Level Triangle

1. **Position 3 mics in equilateral triangle:**
   ```
   Mic 1: [1.0, 1.0, 0.0]   - Floor level
   Mic 2: [1.5, 1.0, 0.0]   - 50cm away
   Mic 3: [1.25, 1.43, 0.0] - Forms triangle
   ```

2. **Record 30 seconds** simultaneously

3. **Save files:** `rec1_mic1.wav`, `rec1_mic2.wav`, `rec1_mic3.wav`

4. **Write down exact positions** in notebook

#### Recording 2: Elevated Triangle

1. **Move all 3 mics up** (0.5-1.0m higher):
   ```
   Mic 1: [1.0, 1.0, 0.8]   - Same X,Y, +0.8m Z
   Mic 2: [1.5, 1.0, 0.8]   
   Mic 3: [1.25, 1.43, 0.8]
   ```

2. **Record another 30 seconds**

3. **Save files:** `rec2_mic1.wav`, `rec2_mic2.wav`, `rec2_mic3.wav`

4. **Write down new positions**

#### Optional Recording 3: Different Location

**For maximum accuracy**, take a third recording at a completely different location:
```
Mic 1: [3.0, 2.0, 0.5]   - Different area
Mic 2: [3.5, 2.0, 0.5]
Mic 3: [3.25, 2.43, 0.5]
```

### Running 3-Mic Localization

```bash
cd /path/to/recordings
./localize_multi_recording \
  rec1_mic1.wav 1.0,1.0,0.0 rec1_mic2.wav 1.5,1.0,0.0 rec1_mic3.wav 1.25,1.43,0.0 \
  rec2_mic1.wav 1.0,1.0,0.8 rec2_mic2.wav 1.5,1.0,0.8 rec2_mic3.wav 1.25,1.43,0.8
```

**Format:** Each recording needs 6 arguments:  
`file1 x1,y1,z1 file2 x2,y2,z2 file3 x3,y3,z3`

### Interpreting Multi-Recording Results

The program will show:
```
Individual estimates:
  Recording 1: (2.3, 1.8, 1.2) m, weight: 0.55
  Recording 2: (2.4, 1.7, 1.3) m, weight: 0.45

Final estimate (weighted average):
  Position: (2.35, 1.75, 1.25) meters
  Consistency (std dev): 0.07 m
```

**Good result:** Consistency < 0.2m  
**Poor result:** Consistency > 0.5m (source may have moved, or bad positions)

### Tips for Best Results

1. **Use equilateral triangle spacing** (all sides equal, ~30-50cm)
2. **Vary height between recordings** (at least 50cm difference)
3. **Ensure source stays constant** (check sound level on phone before each recording)
4. **Measure positions carefully** (accuracy of results depends on position accuracy)
5. **Keep same mic order** (mic1 is always the same device across recordings)

See `MULTI_RECORDING_GUIDE.md` for comprehensive details.

### Common Positioning Mistakes

❌ **Bad:** All mics on floor (planar)
❌ **Bad:** All mics in straight line (linear)
❌ **Bad:** All mics clustered together (< 1m apart)
❌ **Bad:** Mics only in one half of room

✅ **Good:** Mics spread across room at different heights
✅ **Good:** Tetrahedral or rectangular 3D arrangement
✅ **Good:** Mics 2-5m apart

---

## Recording Setup

### 1. Connect Microphones

**USB Microphone Setup:**
```bash
# Check mics are recognized (Linux)
arecord -l

# Should show:
# card 0: Microphone1
# card 1: Microphone2
# card 2: Microphone3
# card 3: Microphone4
```

**Audio Interface Setup:**
- Connect all 4 mics to interface
- Connect interface to computer via USB
- Install drivers if needed
- Set phantom power ON if using condenser mics

### 2. Test Recording

**Quick test to verify setup:**

```bash
# Record 5 seconds from mic 1 (USB mic)
arecord -D hw:0,0 -d 5 -f S16_LE -r 48000 test_mic1.wav

# Play back
aplay test_mic1.wav
```

**Check:**
- [ ] Can hear the drone sound clearly
- [ ] Minimal clipping (volume not too high)
- [ ] All 4 mics recording simultaneously
- [ ] Sample rate 48000 Hz (recommended)

### 3. Set Recording Levels

**Critical for good results:**

- **Aim for:** -12 dB to -6 dB peak levels
- **Avoid:** Clipping (red lights on interface)
- **Test:** Record 10 seconds, check waveform

**Audacity/SoX visualization:**
```bash
# Install sox if needed
sudo apt-get install sox

# View waveform stats
sox test_mic1.wav -n stats
```

Look for: `Maximum amplitude: 0.5-0.8` (avoid 1.0 = clipping)

### 4. Synchronization Check

**CRITICAL:** All 4 mics must start recording at EXACT same time.

**Methods:**

**Option A: Multi-track recording software**
```bash
# Using Audacity (GUI):
# 1. Set all 4 mics as separate tracks
# 2. Hit record once
# 3. All tracks start simultaneously

# Or use sox for command-line:
rec -c 4 -r 48000 all_mics.wav
```

**Option B: Synchronized USB recording**
```bash
# Start all recordings simultaneously with &
arecord -D hw:0,0 -d 60 -r 48000 -f S16_LE mic1.wav &
arecord -D hw:1,0 -d 60 -r 48000 -f S16_LE mic2.wav &
arecord -D hw:2,0 -d 60 -r 48000 -f S16_LE mic3.wav &
arecord -D hw:3,0 -d 60 -r 48000 -f S16_LE mic4.wav &
wait
```

**Option C: Clap test (if slight async)**
- Clap loudly before/after recording
- Use clap peak to align files in post-processing

---

## Data Collection

### Recording Session

**Environment:**
- [ ] All background noise minimized
- [ ] Drone sound is clearly audible
- [ ] No one moving/talking
- [ ] Doors/windows closed

**Recording parameters:**
- **Duration:** 30-60 seconds (minimum 10 seconds)
- **Sample rate:** 48000 Hz (recommended)
- **Bit depth:** 16-bit (minimum)
- **Format:** WAV (uncompressed)
- **Channels:** Mono per mic (4 separate files or 1 multi-channel)

**Procedure:**

1. **Pre-recording:**
   - Verify all mics are on and positioned correctly
   - Check recording levels (-12 to -6 dB)
   - Note exact time of recording (for reference)

2. **Recording:**
   - Start recording on all 4 mics simultaneously
   - Stay silent and still for entire duration
   - Monitor levels (no clipping)
   - Record for 30-60 seconds

3. **Post-recording:**
   - Save files with clear names: `mic1.wav`, `mic2.wav`, `mic3.wav`, `mic4.wav`
   - Note any anomalies (unexpected noises, mic movement, etc.)
   - Do NOT move mics yet (in case you need to re-record)

### Recording Checklist

- [ ] Duration: _____ seconds
- [ ] Sample rate: 48000 Hz
- [ ] Format: WAV
- [ ] All 4 files same length
- [ ] Drone sound clearly audible
- [ ] No clipping occurred
- [ ] Files saved with clear names

---

## Running the Localization

### 1. Prepare Configuration File

**Create `mic_positions.txt`:**
```
# Microphone positions in meters (X, Y, Z)
# Format: one line per microphone
1.5 1.0 0.4
4.0 3.5 1.2
4.5 0.8 0.5
0.5 3.0 2.0
```

**Create `room_config.txt`:**
```
# Room dimensions (length, width, height in meters)
5.0 4.0 2.8

# Speed of sound (m/s) - adjust for temperature
343.0

# Temperature (Celsius) - optional, for automatic speed adjustment
20.0
```

**Speed of sound adjustment for temperature:**
```
c = 331.3 + 0.606 * T  (where T is temperature in Celsius)

Examples:
15°C → 340.4 m/s
20°C → 343.4 m/s
25°C → 346.5 m/s
```

### 2. Run the Localization

**Using the file-based localizer:**

```bash
cd /home/will/projects/hum-finder/build

# Run with default parameters (uses EKF - recommended)
./localize_from_files ../mic_positions.txt mic1.wav mic2.wav mic3.wav mic4.wav

# Or specify estimator
./localize_from_files --estimator ekf ../mic_positions.txt mic1.wav mic2.wav mic3.wav mic4.wav
```

**Expected output:**
```
Loading microphone positions...
  Mic 1: [1.5, 1.0, 0.4] m
  Mic 2: [4.0, 3.5, 1.2] m
  Mic 3: [4.5, 0.8, 0.5] m
  Mic 4: [0.5, 3.0, 2.0] m

Loading audio files...
  mic1.wav: 48000 Hz, 60.0 seconds
  mic2.wav: 48000 Hz, 60.0 seconds
  mic3.wav: 48000 Hz, 60.0 seconds
  mic4.wav: 48000 Hz, 60.0 seconds

Computing TDOA using GCC-PHAT...
  TDOA(1,2): 0.00234 s
  TDOA(1,3): 0.00567 s
  TDOA(1,4): -0.00123 s
  TDOA(2,3): 0.00333 s
  TDOA(2,4): -0.00357 s
  TDOA(3,4): -0.00690 s

Running EKF estimator...
  Initial guess: [2.5, 2.0, 1.4] m (room center)
  
  Convergence over time:
    0.0s: [2.45, 1.87, 1.56] m - error: 0.52 m
    1.0s: [2.38, 1.94, 1.48] m - error: 0.31 m
    2.0s: [2.34, 1.98, 1.43] m - error: 0.18 m
    5.0s: [2.31, 2.01, 1.41] m - error: 0.09 m

FINAL ESTIMATE:
  Position: [2.31, 2.01, 1.41] m
  Uncertainty: ±[0.08, 0.12, 0.15] m
  
  Interpretation:
    X: 2.31 m from front wall
    Y: 2.01 m from left wall
    Z: 1.41 m from floor (approximately waist height)
    
  Location: INSIDE the room (all coordinates within bounds)
  Confidence: HIGH (uncertainty < 0.2m)
```

### 3. Choosing the Right Estimator

**Estimator comparison:**

| Estimator | Speed | Accuracy | Best For |
|-----------|-------|----------|----------|
| **EKF** ⭐ | Fast | Best (0.01m) | Most situations - RECOMMENDED |
| LS | Fast | Good (0.02m) | Quick estimates, sources near mics |
| WLS | Fast | Good (0.02m) | When you know measurement noise |
| ML | Slow | Variable | Difficult geometry, multiple local minima |
| UKF | Medium | Good (0.02m) | Highly nonlinear problems |
| CRLB | N/A | N/A | Theoretical bound only (not an estimator) |

**Recommendation: Start with EKF**

```bash
./localize_from_files --estimator ekf mic_positions.txt mic1.wav mic2.wav mic3.wav mic4.wav
```

If EKF gives poor results, try:
1. ML (slower but more robust)
2. UKF (better for nonlinear cases)

---

## Interpreting Results

### 1. Understanding the Output

**Position estimate: [X, Y, Z]**

```
Example: [2.5, 6.1, 1.2] m
Room: 5.0 × 4.0 × 2.8 m

Analysis:
  X = 2.5 m → Within room (0 to 5.0)
  Y = 6.1 m → OUTSIDE room (room width only 4.0)
  Z = 1.2 m → Normal height (counter-level)

Conclusion: Source is 2.1m beyond the right wall (6.1 - 4.0 = 2.1)
            Likely: Neighbor's apartment appliance
```

**Uncertainty values:**

```
Uncertainty: ±[0.05, 0.08, 0.12] m

Interpretation:
  Small (< 0.2 m): HIGH confidence
  Medium (0.2-0.5 m): MODERATE confidence  
  Large (> 0.5 m): LOW confidence - check setup
```

### 2. Locating the Source

**In-room source:**
```
Estimate: [1.2, 0.8, 0.3] m
Room: 5.0 × 4.0 × 2.8 m

All coordinates within room bounds
→ Look near front-left corner, floor level
→ Check behind furniture, under appliances
```

**Out-of-room source:**
```
Estimate: [3.0, -0.5, 1.5] m
Y = -0.5 (negative!) → Beyond left wall

Direction: Left wall, mid-height
→ Check neighbor's apartment on that side
→ Look for shared walls, pipes, vents
```

**Above/below source:**
```
Estimate: [2.5, 2.0, 4.5] m
Z = 4.5 m (room height only 2.8)

→ Source is UPSTAIRS (4.5 - 2.8 = 1.7m into apartment above)
→ Check for HVAC, pipes, appliances in upstairs unit
```

### 3. Validation

**Sanity checks:**

- [ ] **Position is physically plausible** (not in outer space)
- [ ] **Direction makes sense** (matches where you hear it loudest)
- [ ] **Uncertainty is reasonable** (< 0.5m for in-room, < 2m for distant)
- [ ] **Multiple runs give similar results** (within uncertainty bounds)

**If results seem wrong:**
1. Verify mic positions are correct
2. Check speed of sound value (temperature-dependent)
3. Try different estimator (EKF → ML)
4. Record again with better SNR

---

## Troubleshooting

### Poor Accuracy / High Uncertainty

**Possible causes:**

1. **Bad microphone geometry**
   - ❌ All mics in a plane (linear/planar arrangement)
   - ✅ Fix: Spread mics in 3D, vary heights significantly

2. **Low SNR (signal-to-noise ratio)**
   - ❌ Background noise too high
   - ✅ Fix: Record at quieter time, turn off other noise sources

3. **Incorrect mic positions**
   - ❌ Measurement errors, wrong coordinate system
   - ✅ Fix: Re-measure carefully, verify origin corner

4. **Wrong speed of sound**
   - ❌ Using default 343 m/s in very hot/cold room
   - ✅ Fix: Adjust for temperature: `c = 331.3 + 0.606 * T`

5. **Source too far away**
   - ❌ Source > 20m from mics
   - ✅ Fix: Add more microphones, accept lower accuracy

6. **Multiple sources**
   - ❌ Several drones at similar levels
   - ✅ Fix: Try filtering by frequency, record when only one active

### Estimator Diverges / Crashes

**EKF/UKF specific issues:**

1. **Bad initial guess**
   - ✅ Fix: Use room center as initial guess
   - ✅ Or: Run ML first, use result as EKF initial guess

2. **Process/measurement noise too small**
   - ✅ Fix: Increase covariance matrices (see Advanced section)

3. **Singular covariance matrix**
   - ✅ Fix: Check mic geometry (likely planar arrangement)

### Inconsistent Results

**Different runs give very different positions:**

1. **Source is intermittent** (turning on/off)
   - ✅ Fix: Verify sound is truly continuous, record longer

2. **Microphones moved between recordings**
   - ✅ Fix: Mark positions, don't move mics until testing complete

3. **Multiple sources** (algorithm finding different ones)
   - ✅ Fix: Identify and isolate by frequency

4. **Time synchronization issues**
   - ✅ Fix: Use proper multi-channel recording setup

---

## Advanced: Tuning Parameters

### Measurement Noise Covariance (R Matrix)

**Default values (in code):**
```cpp
double sigma_toa = 2.5 / fs;  // 2.5 samples at sample rate fs
Eigen::Matrix4d R = Eigen::Matrix4d::Identity() * (sigma_toa * sigma_toa);
```

**When to adjust:**

**Increase R (higher noise) if:**
- Broadband source (water flow, air turbulence)
- Poor SNR (< 10 dB)
- Multipath-rich environment (lots of reflections)

```cpp
// For broadband sources
double sigma_toa = 3.5 / fs;  // Increase from 2.5 to 3.5

// For very noisy environment  
double sigma_toa = 5.0 / fs;
```

**Decrease R (lower noise) if:**
- Pure tone source (single frequency)
- Excellent SNR (> 20 dB)
- Anechoic-like environment

```cpp
// For pure tones in quiet environment
double sigma_toa = 1.5 / fs;
```

### Process Noise Covariance (Q Matrix)

**For stationary (non-moving) sources:**

```cpp
// Very small - source doesn't move
Eigen::Matrix3d Q = Eigen::Matrix3d::Identity() * 1e-8;
```

**For moving sources** (if tracking):
```cpp
// Higher - source position varies
Eigen::Matrix3d Q = Eigen::Matrix3d::Identity() * 1e-4;
```

### Initial Covariance (P0 Matrix)

**Reflects uncertainty in initial guess:**

```cpp
// Default: moderate uncertainty (2m standard deviation)
Eigen::Matrix3d P0 = Eigen::Matrix3d::Identity() * 4.0;

// High uncertainty (bad initial guess)
Eigen::Matrix3d P0 = Eigen::Matrix3d::Identity() * 10.0;

// Low uncertainty (good initial guess from ML)
Eigen::Matrix3d P0 = Eigen::Matrix3d::Identity() * 1.0;
```

### Initial Guess Selection

**Strategies:**

1. **Room center** (default, usually good):
   ```cpp
   ekf_model->x0 << room_length/2, room_width/2, room_height/2;
   ```

2. **Suspected location** (if you have a clue):
   ```cpp
   ekf_model->x0 << 3.0, 1.5, 1.0;  // Near suspected appliance
   ```

3. **ML result** (most robust):
   ```cpp
   auto [x_ml, P_ml] = ml(sensor_model, measurements);
   ekf_model->x0 = x_ml;
   ```

### Frequency Filtering

**If multiple sources at different frequencies:**

```bash
# Pre-filter audio to isolate specific frequency range
sox mic1.wav mic1_filtered.wav sinc 50-200  # Bandpass 50-200 Hz

# Then run localization on filtered files
./localize_from_files mic_positions.txt mic1_filtered.wav mic2_filtered.wav ...
```

**Useful filters:**
- HVAC motor: 50-150 Hz
- Electrical transformer: 50-180 Hz (fundamental + harmonics)
- Water flow: 100-2000 Hz (broadband)

---

## Real-World Example Walkthrough

### Scenario: Mystery Hum in Apartment

**Problem:** Constant low-frequency hum, source unknown, worse at night.

### Step 1: Planning (30 minutes)

- Confirmed sound is 24/7 continuous
- Seems to come from left side of room
- Recorded on phone - definitely audible
- Room: 5.2m × 4.1m × 2.8m
- Scheduled test for 2 AM (quietest time)

### Step 2: Equipment (Borrowed/Purchased)

- 4× Blue Snowball USB mics (borrowed from friends)
- Laptop with 4 USB ports
- USB extension cables
- Measuring tape
- Notebook

### Step 3: Mic Positioning (20 minutes)

**Placed mics to form 3D array:**
```
Mic 1: [1.8, 1.2, 0.5] m - Coffee table
Mic 2: [4.5, 3.8, 1.1] m - Bookshelf  
Mic 3: [4.8, 0.9, 0.4] m - Floor near corner
Mic 4: [0.6, 3.5, 2.2] m - Wall-mounted (tape + thumbtack)
```

### Step 4: Environment Prep (10 minutes)

- Unplugged refrigerator (not the source)
- Turned off all fans, AC
- Closed windows
- Put phone in airplane mode
- Hum still clearly audible ✓

### Step 5: Recording (5 minutes)

```bash
# Test levels first
arecord -D hw:0,0 -d 5 -r 48000 -f S16_LE test.wav
aplay test.wav  # Hum is clear, no clipping ✓

# Record 60 seconds from all 4 mics
arecord -D hw:0,0 -d 60 -r 48000 -f S16_LE mic1.wav &
arecord -D hw:1,0 -d 60 -r 48000 -f S16_LE mic2.wav &
arecord -D hw:2,0 -d 60 -r 48000 -f S16_LE mic3.wav &
arecord -D hw:3,0 -d 60 -r 48000 -f S16_LE mic4.wav &
wait
```

### Step 6: Running Localization (2 minutes)

```bash
# Created mic_positions.txt
echo "1.8 1.2 0.5" > mic_positions.txt
echo "4.5 3.8 1.1" >> mic_positions.txt
echo "4.8 0.9 0.4" >> mic_positions.txt
echo "0.6 3.5 2.2" >> mic_positions.txt

# Ran EKF
cd ~/hum-finder/build
./localize_from_files --estimator ekf ../mic_positions.txt \
    mic1.wav mic2.wav mic3.wav mic4.wav
```

### Step 7: Results

```
FINAL ESTIMATE:
  Position: [3.2, -0.7, 1.3] m
  Uncertainty: ±[0.15, 0.22, 0.18] m
  
  Interpretation:
    X: 3.2 m from front wall (within room)
    Y: -0.7 m from left wall (OUTSIDE - 0.7m beyond wall!)
    Z: 1.3 m from floor (counter height)
    
  Location: Beyond LEFT WALL
  Direction: Neighbor's apartment
  Height: Counter/appliance level
  Confidence: HIGH
```

### Step 8: Investigation

- Knocked on neighbor's door
- Sound was their refrigerator against shared wall!
- They moved it 1m away from wall
- Hum reduced by ~80% ✓

**Total time:** ~1 hour setup + testing, solved long-standing mystery!

---

## Checklist Summary

### Pre-Test
- [ ] Sound is continuous (not intermittent)
- [ ] Room measured and sketched
- [ ] Equipment obtained (4 mics + computer)
- [ ] Quiet time scheduled

### Setup
- [ ] Mics positioned in 3D array (not planar!)
- [ ] Positions measured and recorded
- [ ] All background noise minimized
- [ ] Recording levels tested (-12 to -6 dB)

### Recording
- [ ] All 4 mics recording simultaneously
- [ ] Duration 30-60 seconds minimum
- [ ] Sample rate 48000 Hz
- [ ] No clipping occurred
- [ ] Files saved with clear names

### Processing
- [ ] Mic positions entered correctly
- [ ] Speed of sound adjusted for temperature
- [ ] Estimator chosen (EKF recommended)
- [ ] Results seem physically plausible

### Validation
- [ ] Multiple recordings give consistent results
- [ ] Direction matches subjective perception
- [ ] Uncertainty is reasonable (< 0.5m)
- [ ] Source located and verified!

---

## Quick Reference: Command Cheatsheet

```bash
# Test mic is working
arecord -D hw:0,0 -d 5 -r 48000 -f S16_LE test.wav && aplay test.wav

# Record 60s from 4 USB mics simultaneously
arecord -D hw:0,0 -d 60 -r 48000 -f S16_LE mic1.wav &
arecord -D hw:1,0 -d 60 -r 48000 -f S16_LE mic2.wav &
arecord -D hw:2,0 -d 60 -r 48000 -f S16_LE mic3.wav &
arecord -D hw:3,0 -d 60 -r 48000 -f S16_LE mic4.wav &
wait

# Check audio file info
soxi mic1.wav

# View waveform stats
sox mic1.wav -n stats

# Filter to specific frequency range (if needed)
sox mic1.wav mic1_filtered.wav sinc 50-200

# Run localization (EKF - recommended)
./localize_from_files --estimator ekf mic_positions.txt \
    mic1.wav mic2.wav mic3.wav mic4.wav

# Try different estimator if EKF fails
./localize_from_files --estimator ml mic_positions.txt \
    mic1.wav mic2.wav mic3.wav mic4.wav
```

---

## Further Reading

- `DOCS.md` - Complete API and architecture documentation
- `tests/test_real_world_drone.cpp` - Example scenarios with code
- `tests/test_broadband_drone.cpp` - Multi-frequency source validation
- `PRACTICAL_DRONE_FINDING_GUIDE.md` - User-friendly overview

---

**Good luck tracking down that annoying sound! 🎯**
