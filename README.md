# Humming Localization Project

This project implements a sensor network model utilizing three microphone sensors to localize the source of a low-frequency humming sound. The system captures audio data, analyzes it using a Discrete Fourier Transform (DFT), and employs triangulation techniques to determine the sound source's location.

## Project Structure

- **src/**: Contains the source code files.
  - **main.cpp**: Entry point of the application.
  - **microphone.cpp**: Implements the Microphone class for audio capture.
  - **dft.cpp**: Implements the DFT class for frequency analysis.
  - **localizer.cpp**: Implements the Localizer class for sound source triangulation.
  - **utils.cpp**: Contains utility functions for data processing and calculations.

- **include/**: Contains header files for the classes and utility functions.
  - **microphone.h**: Declaration of the Microphone class.
  - **dft.h**: Declaration of the DFT class.
  - **localizer.h**: Declaration of the Localizer class.
  - **utils.h**: Declaration of utility functions.

- **tests/**: Contains unit tests for the project.
  - **test_dft.cpp**: Unit tests for the DFT class.
  - **test_localizer.cpp**: Unit tests for the Localizer class.

- **CMakeLists.txt**: Configuration file for building the project with CMake.

- **.gitignore**: Specifies files and directories to be ignored by Git.

- **LICENSE**: Licensing information for the project.

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