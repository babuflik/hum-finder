#include <iostream>
#include <vector>
#include "microphone.h"
#include "dft.h"
#include "localizer.h"

int main() {
    const int numMicrophones = 3;
    Microphone microphones[numMicrophones];
    
    // Initialize microphones
    for (int i = 0; i < numMicrophones; ++i) {
        if (!microphones[i].initialize(i)) {
            std::cerr << "Failed to initialize microphone " << i << std::endl;
            return -1;
        }
    }

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

    // Perform Discrete Fourier Transform
    DFT dft;
    std::vector<std::vector<std::complex<double>>> frequencyData(numMicrophones);
    for (int i = 0; i < numMicrophones; ++i) {
        frequencyData[i] = dft.performDFT(audioData[i]);
    }

    // Localize the sound source
    Localizer localizer;
    auto sourceLocation = localizer.localize(frequencyData);

    // Output the result
    std::cout << "Estimated source location: (" 
              << sourceLocation.first << ", " 
              << sourceLocation.second << ")" << std::endl;

    return 0;
}