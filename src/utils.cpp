#include "utils.h"
#include <iostream>
#include <vector>
#include <cmath>

// Function to log messages to the console
void logMessage(const std::string& message) {
    std::cout << "[LOG]: " << message << std::endl;
}

// Function to calculate the distance between two points
double calculateDistance(double x1, double y1, double x2, double y2) {
    return std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
}

// Function to normalize a vector
std::vector<double> normalize(const std::vector<double>& vec) {
    double magnitude = 0.0;
    for (double val : vec) {
        magnitude += val * val;
    }
    magnitude = std::sqrt(magnitude);
    
    std::vector<double> normalizedVec;
    for (double val : vec) {
        normalizedVec.push_back(val / magnitude);
    }
    return normalizedVec;
}

// Function to convert frequency to wavelength
double frequencyToWavelength(double frequency) {
    const double speedOfSound = 343.0; // Speed of sound in air in m/s
    return speedOfSound / frequency;
}