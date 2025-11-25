# Real-Time Estimation Architecture

## System Diagram

```
┌─────────────────────────────────────────────────────────────────────┐
│                        HUMMING LOCALIZATION                         │
│                         Real-Time System                            │
└─────────────────────────────────────────────────────────────────────┘

┌──────────────┐     ┌──────────────┐     ┌──────────────┐     ┌──────────────┐
│  Microphone  │     │  Microphone  │     │  Microphone  │     │  Microphone  │
│      #1      │     │      #2      │     │      #3      │     │      #4      │
│  (0,0,0)     │     │  (0.2,0,0)   │     │  (0.2,0.2,0) │     │  (0,0.2,0)   │
└──────┬───────┘     └──────┬───────┘     └──────┬───────┘     └──────┬───────┘
       │                    │                    │                    │
       │                    │                    │                    │
       │    Capture 1024 samples @ 44.1kHz (23ms buffer)             │
       │                    │                    │                    │
       ▼                    ▼                    ▼                    ▼
┌────────────────────────────────────────────────────────────────────────────┐
│                          AUDIO BUFFERS (4 channels)                        │
│  Buffer[0]: [s₀, s₁, s₂, ..., s₁₀₂₃]   (1024 samples)                    │
│  Buffer[1]: [s₀, s₁, s₂, ..., s₁₀₂₃]   (1024 samples)                    │
│  Buffer[2]: [s₀, s₁, s₂, ..., s₁₀₂₃]   (1024 samples)                    │
│  Buffer[3]: [s₀, s₁, s₂, ..., s₁₀₂₃]   (1024 samples)                    │
└────────────────────────────┬───────────────────────────────────────────────┘
                             │
                             ▼
┌────────────────────────────────────────────────────────────────────────────┐
│                        TDOA CALCULATOR (GCC-PHAT)                          │
│  ┌──────────────────────────────────────────────────────────────────┐     │
│  │ For each mic pair (i,j):                                         │     │
│  │   1. Bandpass filter 100-300 Hz (FFT domain)                     │     │
│  │   2. FFT both signals                                            │     │
│  │   3. Cross-correlation: R[k] = Y_i[k] · Y_j*[k] / |Y_i·Y_j*|    │     │
│  │   4. IFFT to get correlation                                     │     │
│  │   5. Find peak → TDOA in samples                                 │     │
│  └──────────────────────────────────────────────────────────────────┘     │
│                                                                            │
│  Output: z = [τ₀₁, τ₀₂, τ₀₃, τ₁₂, τ₁₃, τ₂₃]ᵀ  (6 TDOA measurements)     │
└────────────────────────────┬───────────────────────────────────────────────┘
                             │
                             ▼
┌────────────────────────────────────────────────────────────────────────────┐
│                    KALMAN FILTER (EKF or UKF)                              │
│                                                                            │
│  State: x = [x, y, z, vx, vy, vz]ᵀ  (position + velocity)                │
│                                                                            │
│  ┌──────────────────────────────────────────────────────────────────┐     │
│  │ PREDICTION STEP                                                  │     │
│  │   x̂⁻(k) = F · x̂(k-1)                                            │     │
│  │   P⁻(k)  = F · P(k-1) · Fᵀ + Q                                  │     │
│  │                                                                  │     │
│  │   where F = [I₃  Δt·I₃]  (constant velocity model)              │     │
│  │             [0₃   I₃  ]                                          │     │
│  └──────────────────────────────────────────────────────────────────┘     │
│                                                                            │
│  ┌──────────────────────────────────────────────────────────────────┐     │
│  │ UPDATE STEP                                                      │     │
│  │   Innovation: y = z - h(x̂⁻)                                     │     │
│  │                                                                  │     │
│  │   where h(x) = [(‖x-m₁‖ - ‖x-m₂‖)/c, ...]  (TDOA model)        │     │
│  │                                                                  │     │
│  │   Kalman Gain: K = P⁻ · Hᵀ · (H·P⁻·Hᵀ + R)⁻¹                   │     │
│  │                                                                  │     │
│  │   State update: x̂(k) = x̂⁻(k) + K · y                           │     │
│  │   Cov update:   P(k)  = P⁻(k) - K · H · P⁻(k)                  │     │
│  └──────────────────────────────────────────────────────────────────┘     │
└────────────────────────────┬───────────────────────────────────────────────┘
                             │
                             ▼
                    ┌────────────────┐
                    │ Position (x,y,z)│
                    │ + Uncertainty P │
                    └────────────────┘
```

## Data Flow Timeline

```
Time (ms)  │ Activity
───────────┼─────────────────────────────────────────────────────────────
    0      │ ┌─ Capture audio buffer (1024 samples)
           │ │
   ~23     │ └─ Buffer full → TDOA calculation begins
           │    ├─ Bandpass filter (FFT)        [~3 ms]
           │    ├─ Cross-correlation (6 pairs)  [~4 ms]
           │    └─ Peak detection                [~1 ms]
           │
   ~31     │ TDOA complete → Kalman filter update begins
           │    ├─ Compute Jacobian              [~1 ms]
           │    ├─ Matrix operations (6×6)       [<1 ms]
           │    └─ State/Covariance update       [<1 ms]
           │
   ~34     │ Position estimate ready ✓
           │    └─ Return (x, y, z)
           │
           │ [Next buffer starts capturing...]
   ~46     │ ┌─ Next buffer ready
           │ ...
```

## State Persistence (Key to Real-Time)

```
Call #1:                    Call #2:                    Call #3:
┌──────────┐               ┌──────────┐               ┌──────────┐
│ Buffer 1 │               │ Buffer 2 │               │ Buffer 3 │
└────┬─────┘               └────┬─────┘               └────┬─────┘
     │                          │                          │
     ▼                          ▼                          ▼
┌─────────────┐            ┌─────────────┐            ┌─────────────┐
│ locateSource│            │ locateSource│            │ locateSource│
│             │            │             │            │             │
│ x̂ = [0.5,  │            │ x̂ = [0.48, │            │ x̂ = [0.49, │
│      0.5,  │  ───────>  │      0.51, │  ───────>  │      0.50, │
│      0.0]  │  (update)  │      0.02] │  (update)  │      0.01] │
│             │            │             │            │             │
│ P = 1.0·I  │            │ P = 0.8·I  │            │ P = 0.6·I  │
└─────────────┘            └─────────────┘            └─────────────┘
     │                          │                          │
     └─────────────┬────────────┴───────────────┬──────────┘
                   │                            │
              State persists between calls      │
              (Kalman filter memory)            │
                                    Uncertainty decreases
                                    as more data arrives
```

## Performance Metrics

```
┌────────────────────────────────────────────┐
│ Processing Budget (23 ms per buffer)      │
├────────────────────────────────────────────┤
│ Audio capture:          0 ms (async)       │
│ TDOA (GCC-PHAT):       ~8 ms               │
│   ├─ FFT/IFFT:         ~5 ms               │
│   └─ Peak finding:     ~3 ms               │
│ Kalman filter:         ~2 ms               │
│   ├─ Jacobian:         ~1 ms               │
│   └─ Matrix ops:       ~1 ms               │
│ Total:                ~10 ms               │
│                                            │
│ Margin:               ~13 ms (56%)         │
└────────────────────────────────────────────┘

Update Rate: 43.1 Hz (1000 / 23.2 ms)
Latency:     23-34 ms (buffer + processing)
```

## Thread Model (realtime_localizer.cpp)

```
┌─────────────────────────────────────────────────────────────┐
│                       Main Thread                           │
│                                                             │
│  ┌──────────────────────────────────────────────────────┐  │
│  │ RealtimeLocalizer::start()                           │  │
│  │   └─ Spawn processing thread                         │  │
│  └──────────────────────────────────────────────────────┘  │
│                                                             │
│  ┌──────────────────────────────────────────────────────┐  │
│  │ RealtimeLocalizer::getLatestPosition()               │  │
│  │   └─ Read shared state (mutex protected)             │  │
│  └──────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────┘
                              │
                              │ spawns
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                    Processing Thread                        │
│                                                             │
│  while (running) {                                          │
│    ┌──────────────────────────────────────────────────┐    │
│    │ 1. captureAudioBuffers()                         │    │
│    └──────────────────────────────────────────────────┘    │
│                        │                                    │
│                        ▼                                    │
│    ┌──────────────────────────────────────────────────┐    │
│    │ 2. localizer.locateSource(buffers)               │    │
│    │    [TDOA + Kalman filter]                        │    │
│    └──────────────────────────────────────────────────┘    │
│                        │                                    │
│                        ▼                                    │
│    ┌──────────────────────────────────────────────────┐    │
│    │ 3. Update shared state                           │    │
│    │    { mutex_lock                                  │    │
│    │      latest_position = position;                 │    │
│    │      latest_covariance = P;                      │    │
│    │    }                                             │    │
│    └──────────────────────────────────────────────────┘    │
│                        │                                    │
│                        ▼                                    │
│    ┌──────────────────────────────────────────────────┐    │
│    │ 4. sleep(23 ms)                                  │    │
│    └──────────────────────────────────────────────────┘    │
│  }                                                          │
└─────────────────────────────────────────────────────────────┘
```

## Filter Convergence

```
Uncertainty (trace of P) over time:

P (m²)
  │
6 │ ╭─╮
  │ │ │                Initial high uncertainty
5 │ │ │╮
  │ │ ││╮
4 │ │ │││╮             Filter converges as
  │ │ │││╰╮            measurements arrive
3 │ │ │││ ╰╮
  │ │ │││  ╰╮
2 │ │ │││   ╰╮
  │ │ │││    ╰─╮
1 │ │ │││      ╰───────────────────  Steady state
  │ │ │││
0 └─┴─┴┴┴──────────────────────────────────────► Time
  0   1  2  3  4  5  6  7  8  9  10 (seconds)
      Updates (43 Hz)
```

This architecture enables **continuous, low-latency** position tracking!
