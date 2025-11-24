# Humming Localization Project

This project implements a sensor network model utilizing three microphone sensors to localize the source of a low-frequency humming sound. The system captures audio data, analyzes it using a Discrete Fourier Transform (DFT), and employs triangulation techniques to determine the sound source's location.

**Extended Implementation**: This project also includes a comprehensive C++ port of the MATLAB SigSys Toolbox, providing probability distributions, nonlinear filtering (EKF/UKF/PF), and signal processing capabilities.

## Quick Start

### Build
```bash
cd build
cmake ..
make -j4
```

### Run Tests
```bash
ctest --output-on-failure  # All tests: 28/28 passing
```

### Run Localization
```bash
./humming_localization
```

### Visualize Results
```bash
make plot      # 2D plot
make plot3d ELEV=30 AZIM=45  # 3D plot with custom viewing angles
```

## Current Implementation Status

- ✅ **10 Probability Distributions**: NDist, UDist, ExpDist, Chi2Dist, GammaDist, GMDist, TDist, BetaDist, EmpDist, LogNDist
- ✅ **NL Filtering Framework**: EKF, UKF, Particle Filter with simulation
- ✅ **Signal Processing**: Sig class with statistics, extraction, downsampling
- ✅ **Numerical Utilities**: Gradient, Hessian, Jacobian, covariance tools
- ✅ **28 Tests Passing**: Full test coverage for all new functionality

See `CONTINUATION_GUIDE.md` for implementation details and next steps.

## Project Structure

  - **main.cpp**: Entry point of the application.
  - **microphone.cpp**: Implements the Microphone class for audio capture.
  - **dft.cpp**: Implements the DFT class for frequency analysis.
  - **localizer.cpp**: Implements the Localizer class for sound source triangulation.
  - **utils.cpp**: Contains utility functions for data processing and calculations.

  - **microphone.h**: Declaration of the Microphone class.
  - **dft.h**: Declaration of the DFT class.
  - **localizer.h**: Declaration of the Localizer class.
  - **utils.h**: Declaration of utility functions.

  - **test_dft.cpp**: Unit tests for the DFT class.
  - **test_localizer.cpp**: Unit tests for the Localizer class.




## Setup Instructions

1. Clone the repository:
   ```
   git clone <repository-url>
   cd humming-localization-cpp
   ```

2. Create a build directory and navigate into it:
   ```
   mkdir build
   cd build
   ```

3. Run CMake to configure the project:
   ```
   cmake ..
   ```

4. Build the project:
   ```
   make
   ```

## Dependencies

This project requires several system and development libraries to build and run. Key dependencies are discovered in `CMakeLists.txt` and `cmake/modules/FindFFTW3.cmake`.


Quick install examples (Debian/Ubuntu):

```bash
sudo apt update
sudo apt install -y build-essential cmake libeigen3-dev libfftw3-dev libgtest-dev
# On some Ubuntu versions, libgtest-dev installs sources only; build and install the library:
cd /usr/src/gtest || true
sudo cmake .
sudo make
sudo cp lib/*.a /usr/lib || true
```

Fedora (example):

```bash
sudo dnf install -y cmake gcc-c++ eigen3-devel fftw-devel gtest-devel
```

macOS (Homebrew):

```bash
brew update
brew install cmake eigen fftw googletest
```

If `FFTW3` or `Eigen3` are installed in non-standard locations, set `CMAKE_PREFIX_PATH` or provide `-DFFTW3_INCLUDE_DIR`/`-DFFTW3_LIBRARY` when invoking `cmake`.

Important files referencing these dependencies:


## Usage

After building the project, you can run the application using the following command:
```
./humming-localization-cpp
```

The application will initialize the microphone sensors, capture audio data, and process it to localize the source of the humming sound.

## Contributing

Contributions are welcome! Please feel free to submit a pull request or open an issue for any suggestions or improvements.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.

## Plotting estimation results

I added a plotting utility in `tools/plot_estimates.py` and sample-data generators to visualize sensor layouts, target positions, estimates and confidence ellipses/ellipsoids.

### Quick Start (using Makefile)

The easiest way to generate plots after running tests:

```bash
# 2D plotting (using test artifacts)
make plot

# 3D plotting with custom viewing angle
make plot3d ELEV=20 AZIM=45

# Default viewing angle: ELEV=30 AZIM=-60
make plot3d
```

### 2D Plotting

Manual command:

```bash
# from repo root
python3 tools/generate_sample_data.py
python3 tools/plot_estimates.py --sensors artifacts/sensors.csv --targets artifacts/targets.csv \
   --estimates artifacts/estimates.csv --cov artifacts/covariances.npy --outfile artifacts/plot.png --show
```

Or simply:
```bash
make plot
```

### 3D Plotting

For 3D visualization with confidence ellipsoids:

```bash
# Generate 3D sample data
python3 tools/generate_sample_data_3d.py

# Plot in 3D with custom viewing angle
python3 tools/plot_estimates.py --sensors artifacts/sensors_3d.csv --targets artifacts/targets_3d.csv \
   --estimates artifacts/estimates_3d.csv --cov artifacts/covariances_3d.npy --3d \
   --elev 20 --azim 45 --outfile artifacts/plot_3d.png
```

Or use the Makefile for easier commands:
```bash
make plot3d              # Default view (elevation=30°, azimuth=-60°)
make plot3d ELEV=20 AZIM=45   # Custom viewing angle
```

The `--elev` parameter controls the elevation angle (vertical tilt), and `--azim` controls the azimuth (horizontal rotation) of the 3D plot.

Python dependencies:

```bash
pip3 install -r requirements.txt
```

How to export data from C++ tests:

- The project does not yet write CSV/numpy outputs by default. To export data from tests or code, add a small serializer where you have access to sensor positions, estimated positions and covariance matrices. Example C++ snippet (pseudo-code) you can drop into a test or helper:

```cpp
// example: write estimate and 2x2 covariance to CSV / binary
#include <fstream>
#include <Eigen/Dense>

void write_estimate(const std::string &path, const Eigen::Vector2d &x, const Eigen::Matrix2d &P) {
      std::ofstream f(path);
      f << x(0) << "," << x(1) << "\n";
      f << P(0,0) << "," << P(0,1) << "\n";
      f << P(1,0) << "," << P(1,1) << "\n";
}
```

Then load the values into `plot_estimates.py` (you can convert the CSV into the `estimates.csv` / `covariances.npy` formats the script expects). If you'd like, I can add automatic exporting to the test runner or a small C++ helper that writes `sensors.csv`, `estimates.csv` and `covariances.npy` directly from `humming_core` types.

Virtual environment (recommended)

Create and use a Python virtual environment to isolate plotting dependencies:

```bash
# create venv in project root
python3 -m venv .venv
# activate (Linux / macOS bash)
source .venv/bin/activate
# install requirements
pip install -U pip
pip install -r requirements.txt

# run the demo
python3 tools/generate_sample_data.py
python3 tools/plot_estimates.py --sensors sensors.csv --targets targets.csv \
   --estimates estimates.csv --cov covariances.npy --show
```

To deactivate the venv:

```bash
deactivate
```

Makefile convenience

There is a top-level `Makefile` with convenience targets to manage the Python virtual environment:

- `make venv` — creates `.venv` (or the directory set via `VENV_DIR`) and installs `requirements.txt` using the Python interpreter in `PYTHON` (defaults: `.venv` and `python3`).
- `make clean-venv` — removes the venv directory (defaults to `.venv`).

Examples:

```bash
# create venv (default .venv)
make venv

# create venv with custom python or dir
make venv PYTHON=python3.11 VENV_DIR=.venv3.11

# remove venv
make clean-venv
```