# Multi-Recording Localization Guide

## Overview

When your sound source is **stationary** (not moving), you can use **fewer microphones** by taking **multiple recordings**. This approach:

- ✅ Works with only **3 microphones** instead of 4
- ✅ Improves accuracy through averaging
- ✅ Allows you to change microphone positions between recordings
- ✅ Reduces hardware costs
- ⚠️ Requires the source to remain stationary between recordings

## How It Works

### Standard Approach (4 mics, 1 recording)
```
4 microphones → 6 TDOA measurements → 1 position estimate
```

### Multi-Recording Approach (3 mics, N recordings)
```
3 microphones × N recordings → 2N TDOA measurements → N estimates → averaged result
```

Each recording with 3 mics gives 2 independent TDOA measurements. Multiple recordings provide:
1. **Redundancy**: Reduce random measurement errors
2. **Geometric diversity**: Different mic positions improve constraint geometry
3. **Statistical improvement**: Uncertainty decreases with √N recordings

## Usage

### Basic Syntax
```bash
./localize_multi_recording \
  file1 x1,y1,z1 file2 x2,y2,z2 file3 x3,y3,z3 \
  [more recordings...]
```

### Example 1: Single Recording (3 mics)
```bash
./localize_multi_recording \
  mic1.wav 0,0,0 \
  mic2.wav 0.2,0,0 \
  mic3.wav 0.2,0.2,0
```

This is the minimum setup - works but less accurate than 4 mics.

### Example 2: Two Recordings (better accuracy)
```bash
./localize_multi_recording \
  rec1_mic1.wav 0,0,0 rec1_mic2.wav 0.2,0,0 rec1_mic3.wav 0.2,0.2,0 \
  rec2_mic1.wav 0,0,0.5 rec2_mic2.wav 0.2,0,0.5 rec2_mic3.wav 0.2,0.2,0.5
```

Here we moved the mic array up by 0.5m for the second recording, giving better vertical accuracy.

### Example 3: Three Recordings (even better)
```bash
./localize_multi_recording \
  rec1_mic1.wav 0,0,0 rec1_mic2.wav 0.2,0,0 rec1_mic3.wav 0.2,0.2,0 \
  rec2_mic1.wav 0,0,0.5 rec2_mic2.wav 0.2,0,0.5 rec2_mic3.wav 0.2,0.2,0.5 \
  rec3_mic1.wav 0.3,0,0.2 rec3_mic2.wav 0.5,0,0.2 rec3_mic3.wav 0.5,0.2,0.2
```

Three recordings with different positions provide excellent geometric diversity.

## Practical Workflow

### Scenario: Find a Drone in Your Apartment with Only 3 Phones

**Equipment needed:**
- 3 smartphones with voice recorder apps
- Measuring tape
- Notebook for positions

**Steps:**

1. **Recording 1 - Ground level array**
   ```
   Place phones on floor in triangle pattern:
   - Phone 1: (0, 0, 0) - corner of room
   - Phone 2: (0.5, 0, 0) - 50cm to the right
   - Phone 3: (0.25, 0.43, 0) - forms equilateral triangle
   
   Record 30 seconds simultaneously
   Transfer files: rec1_phone1.wav, rec1_phone2.wav, rec1_phone3.wav
   ```

2. **Recording 2 - Elevated array**
   ```
   Move phones up 1 meter (on table/chairs):
   - Phone 1: (0, 0, 1.0)
   - Phone 2: (0.5, 0, 1.0)
   - Phone 3: (0.25, 0.43, 1.0)
   
   Record 30 seconds simultaneously
   Transfer files: rec2_phone1.wav, rec2_phone2.wav, rec2_phone3.wav
   ```

3. **Recording 3 - Different position (optional, for best results)**
   ```
   Move to different location in room:
   - Phone 1: (1.0, 1.0, 0.5)
   - Phone 2: (1.5, 1.0, 0.5)
   - Phone 3: (1.25, 1.43, 0.5)
   
   Record 30 seconds simultaneously
   Transfer files: rec3_phone1.wav, rec3_phone2.wav, rec3_phone3.wav
   ```

4. **Run localization**
   ```bash
   ./localize_multi_recording \
     rec1_phone1.wav 0,0,0 rec1_phone2.wav 0.5,0,0 rec1_phone3.wav 0.25,0.43,0 \
     rec2_phone1.wav 0,0,1.0 rec2_phone2.wav 0.5,0,1.0 rec2_phone3.wav 0.25,0.43,1.0 \
     rec3_phone1.wav 1.0,1.0,0.5 rec3_phone2.wav 1.5,1.0,0.5 rec3_phone3.wav 1.25,1.43,0.5
   ```

## Microphone Placement Tips

### Same Positions (Easiest)
Keep microphones in the same positions, just take multiple recordings:
- **Pros**: Simple, repeatable
- **Cons**: Only averages measurement noise, doesn't improve geometry
- **Use when**: You want to reduce random errors

### Different Heights
Move array up/down between recordings:
- **Pros**: Improves vertical (z-axis) accuracy significantly
- **Cons**: Need stable mounting at different heights
- **Use when**: Source might be above/below array plane

### Different Locations
Move entire array to different room positions:
- **Pros**: Best geometric diversity, most accurate results
- **Cons**: More complex to measure and record
- **Use when**: You want maximum accuracy with 3 mics

### Optimal Geometry
For best results with 3 mics:
- Use **equilateral triangle** spacing (all sides equal)
- Separate recordings should have **different orientations** relative to source
- Include at least one recording with **significant height difference**

## Expected Accuracy

### Single Recording (3 mics)
- Pure tone: 0.3-0.5 m
- Harmonic sound: 0.5-0.8 m
- Broadband: 0.8-1.5 m

### Two Recordings (3 mics each)
- Pure tone: 0.2-0.3 m
- Harmonic sound: 0.3-0.5 m
- Broadband: 0.5-1.0 m

### Three+ Recordings (3 mics each)
- Pure tone: 0.1-0.2 m
- Harmonic sound: 0.2-0.4 m
- Broadband: 0.4-0.8 m

**Note**: Accuracy approaches 4-mic performance with 3+ diverse recordings.

## Understanding the Output

```
=== Results ===
Individual estimates:
  Recording 1: (1.234, 0.567, 1.890) m, weight: 0.35
  Recording 2: (1.245, 0.571, 1.885) m, weight: 0.40
  Recording 3: (1.238, 0.569, 1.892) m, weight: 0.25

Final estimate (weighted average):
  Position: (1.240, 0.569, 1.888) meters
  Consistency (std dev): 0.005 m
  Per-axis std dev: (0.004, 0.002, 0.003) m
```

**Interpretation**:
- **Individual estimates**: Each recording's result and its reliability weight
- **Final estimate**: Weighted average of all recordings
- **Consistency**: How well the recordings agree (lower is better)
  - < 0.1 m: Excellent agreement
  - 0.1-0.3 m: Good agreement
  - > 0.3 m: Poor agreement - check mic positions or source movement
- **Per-axis std dev**: Uncertainty in each dimension

## Troubleshooting

### High Inconsistency (std dev > 0.5m)
**Possible causes**:
- Source moved between recordings
- Incorrect microphone position measurements
- Different recording volumes/clipping
- Poor TDOA correlation (low SNR)

**Solutions**:
- Verify source is truly stationary
- Double-check all position measurements
- Normalize audio levels before recording
- Move closer to source or reduce background noise

### Poor Z-axis (vertical) Accuracy
**Cause**: All recordings at same height

**Solution**: Include recordings at different heights (0.5-1m difference)

### Estimates Are Weighted Unevenly
**Cause**: Some recordings have better geometry or SNR

**This is normal**: The algorithm automatically weights better measurements more heavily

## Comparison: 3 Mics vs 4 Mics

| Aspect | 4 Mics (1 recording) | 3 Mics (3 recordings) |
|--------|---------------------|----------------------|
| Microphones needed | 4 | 3 |
| Recording sessions | 1 | 3 |
| Total recording time | 30 sec | 90 sec |
| Setup complexity | Medium | Higher |
| Accuracy (stationary) | Very Good | Very Good |
| Accuracy (moving) | Very Good | ❌ Won't work |
| Cost | Higher | Lower |
| Best for | Real-time tracking | Stationary sources |

## Advanced: Optimal Recording Strategy

For **maximum accuracy** with **minimum recordings**:

1. **Recording 1**: Ground level, triangular array
   - Spacing: 0.3-0.5 m between mics
   
2. **Recording 2**: Same positions, elevated by 1 m
   - Improves vertical resolution
   
3. **Recording 3**: Different horizontal position, mid-height
   - Provides geometric diversity

This gives you comparable accuracy to a 4-mic array with just 3 microphones!

## When to Use This vs Standard 4-Mic Approach

**Use Multi-Recording (3 mics):**
- ✓ Source is stationary (HVAC, electrical hum, static drone)
- ✓ You only have 3 recording devices
- ✓ You want to minimize hardware cost
- ✓ You have time for multiple recordings
- ✓ You can move microphones between recordings

**Use Standard (4 mics):**
- ✓ Source might move
- ✓ You need real-time tracking
- ✓ You want instant results
- ✓ You have 4+ microphones available
- ✓ Cannot take multiple recordings

## Example Script for Automation

```bash
#!/bin/bash
# auto_locate.sh - Automate 3-mic multi-recording localization

echo "Recording 1 - Ground level"
echo "Press Enter when mics are at: (0,0,0), (0.5,0,0), (0.25,0.43,0)"
read
arecord -d 30 -f S16_LE -r 44100 -c 1 rec1_mic1.wav &
arecord -d 30 -f S16_LE -r 44100 -c 1 rec1_mic2.wav &
arecord -d 30 -f S16_LE -r 44100 -c 1 rec1_mic3.wav &
wait

echo "Recording 2 - Elevated 1m"
echo "Press Enter when mics are at: (0,0,1), (0.5,0,1), (0.25,0.43,1)"
read
arecord -d 30 -f S16_LE -r 44100 -c 1 rec2_mic1.wav &
arecord -d 30 -f S16_LE -r 44100 -c 1 rec2_mic2.wav &
arecord -d 30 -f S16_LE -r 44100 -c 1 rec2_mic3.wav &
wait

echo "Running localization..."
./localize_multi_recording \
  rec1_mic1.wav 0,0,0 rec1_mic2.wav 0.5,0,0 rec1_mic3.wav 0.25,0.43,0 \
  rec2_mic1.wav 0,0,1 rec2_mic2.wav 0.5,0,1 rec2_mic3.wav 0.25,0.43,1
```

## Summary

Multi-recording localization is a **cost-effective solution** for finding **stationary sound sources** with minimal hardware. By taking advantage of the source's immobility, you can achieve comparable accuracy to a 4-microphone array using only 3 microphones and a few minutes of recording time.

**Key takeaway**: For stationary sources, time can substitute for hardware!
