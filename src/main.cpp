#include <iostream>
#include <vector>
#include "microphone.h"
#include "dft.h"
#include "localizer.h"

int main() {
    
    // Initialize microphones at specific positions
    std::array<Microphone, 3> microphones = {
        Microphone(0.0, 0.0),    // First microphone at origin
        Microphone(1.0, 0.0),    // Second microphone 1m to the right
        Microphone(0.5, 0.866)   // Third microphone forming equilateral triangle
    };

    /*
    // Start capturing audio
    for (int i = 0; i < numMicrophones; ++i) {
        microphones[i].startCapture();
    }

    // Capture audio data
    std::vector<std::vector<double>> audioData(numMicrophones);
    for (int i = 0; i < numMicrophones; ++i) {
        audioData[i] = microphones[i].captureAudio();
    }

    // Stop capturing audio
    for (int i = 0; i < numMicrophones; ++i) {
        microphones[i].stopCapture();
    }
*/
    /*
    // Perform Discrete Fourier Transform
    DFT dft;
    std::vector<std::vector<std::complex<double>>> frequencyData(numMicrophones);
    for (int i = 0; i < numMicrophones; ++i) {
        frequencyData[i] = dft.performDFT(audioData[i]);
    }
    */

    // Create localizer with microphone array
    Localizer localizer(microphones);

    // Sample data processing (replace with actual data collection)
    std::vector<double> frequencyData = {0.5, 0.7, 0.3}; // Example magnitudes

    auto sourceLocation = localizer.locateSource(frequencyData);

    // Output the result
    std::cout << "Source located at: (" 
              << sourceLocation[0] << ", " 
              << sourceLocation[1] << ")" << std::endl;

    return 0;
}