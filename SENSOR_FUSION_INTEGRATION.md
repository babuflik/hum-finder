# Sensor Fusion Integration Guide

> **Note**: For complete documentation, see **[sensor-fusion/README.md](sensor-fusion/README.md)**

The sensor fusion toolbox has been extracted as a standalone library that can be used in other projects.

## Quick Integration

```cmake
# Add to your CMakeLists.txt
add_subdirectory(path/to/sensor-fusion)
target_link_libraries(your_target SensorFusion::sensor_fusion)
```

```cpp
// Use in your code
#include <sensor_fusion/nl.h>          // EKF, UKF, Particle Filter
#include <sensor_fusion/estimators.h>  // LS, WLS, ML, CRLB
#include <sensor_fusion/sensormod.h>   // Sensor modeling
```

For complete documentation including:
- Installation methods (subdirectory, system install, FetchContent)
- API reference
- Working examples
- Integration patterns

See: **[sensor-fusion/README.md](sensor-fusion/README.md)**
